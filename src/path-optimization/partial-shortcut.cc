// Copyright (c) 2015, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr), Mylene Campana
//
// This file is part of hpp-core.
// hpp-core is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-core is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-core. If not, see <http://www.gnu.org/licenses/>.

#include <hpp/core/path-optimization/partial-shortcut.hh>

// #include <limits>
// #include <deque>
// #include <cstdlib>
// #include <hpp/util/assertion.hh>
#include <hpp/util/debug.hh>
#include <hpp/model/joint.hh>
#include <hpp/model/joint-configuration.hh>
#include <hpp/model/device.hh>
#include <hpp/core/distance.hh>
#include <hpp/core/path-validation.hh>
#include <hpp/core/path-vector.hh>
#include <hpp/core/problem.hh>
#include <hpp/core/config-projector.hh>
#include <hpp/core/locked-joint.hh>

#define BILLION 1E9

namespace hpp {
  namespace core {
    namespace pathOptimization {
      namespace {
        void unpack (PathPtr_t path, PathVectorPtr_t out) {
          PathVectorPtr_t pv = HPP_DYNAMIC_PTR_CAST(PathVector, path);
          if (!pv) {
            out->appendPath (path);
          } else {
            for (std::size_t i = 0; i < pv->numberPaths (); ++i)
              unpack (pv->pathAtRank (i), out);
          }
        }

        // Compute the length of a vector of paths assuming that each element
        // is optimal for the given distance.
        static value_type pathLength (const PathVectorPtr_t& path,
				      const DistancePtr_t& distance)
        {
          value_type result = 0;
          for (std::size_t i=0; i<path->numberPaths (); ++i) {
            const PathPtr_t& element (path->pathAtRank (i));
            result += (*distance) (element->initial (), element->end ());
          }
          return result;
        }
      }

      PartialShortcut::Parameters::Parameters () :
        removeLockedJoints (true), onlyFullShortcut (true),
        numberOfConsecutiveFailurePerJoints (5), progressionMargin (1e-3)
      {}

      PartialShortcutPtr_t PartialShortcut::create (const Problem& problem)
      {
        return createWithTraits <PartialShortcutTraits> (problem);
      }

      PartialShortcut::PartialShortcut (const Problem& problem) :
        PathOptimizer (problem)
      {
      }

      PathVectorPtr_t PartialShortcut::optimize (const PathVectorPtr_t& path)
      {
	problem ().timeValues_.clear ();
	problem ().gainValues_.clear ();
	//clock_gettime(CLOCK_MONOTONIC, &start_); // building jv very long...
	initialLength_ = pathLength (path, problem ().distance ());
	
        PathVectorPtr_t unpacked =
	  PathVector::create (path->outputSize(),path->outputDerivativeSize ());
        unpack (path, unpacked);

        /// Step 1: Generate a suitable vector of joints
        JointVector_t straight_jv = generateJointVector (unpacked);
        JointVector_t jv;

        /// Step 2: First try to optimize each joint from beginning to end
        //PathVectorPtr_t result = optimizeFullPath (unpacked, straight_jv, jv);
        //if (parameters.onlyFullShortcut) return result;

        /// Step 3: Optimize randomly each joint
        //return optimizeRandom (result, jv);
	return optimizeRandom (path, straight_jv);
      }

      PathVectorPtr_t PartialShortcut::generatePath (PathVectorPtr_t path,
						     const JointPtr_t joint,
						     const value_type t1,
						     ConfigurationIn_t q1,
						     const value_type t2,
						     ConfigurationIn_t q2) const
      {
        value_type lt1, lt2;
        std::size_t rkAtP1 = path->rankAtParam (t1, lt1);
        std::size_t rkAtP2 = path->rankAtParam (t2, lt2);
        if (rkAtP2 == rkAtP1) return PathVectorPtr_t ();

        PathVectorPtr_t pv = PathVector::create (path->outputSize (),
						 path->outputDerivativeSize ());
        PathPtr_t last;

        std::size_t rkCfg = joint->rankInConfiguration ();
        Configuration_t qi = q1;
        Configuration_t q_inter (path->outputSize ());
        value_type t = - lt1;
        for (std::size_t i = rkAtP1; i < rkAtP2; ++i) {
          t += path->pathAtRank (i)->timeRange().second;
          q_inter = path->pathAtRank (i)->end (),
	    joint->configuration()->interpolate ( q1, q2,
						  t / (t2-t1), rkCfg, q_inter);
          if (path->pathAtRank (i)->constraints ()) {
            if (!path->pathAtRank (i)->constraints ()->apply (q_inter)) {
              hppDout (warning, "PartialShortcut could not apply "
		       "the constraints");
              return PathVectorPtr_t ();
            }
          }
          last = (steer) (qi, q_inter);
          if (!last) return PathVectorPtr_t ();
          pv->appendPath (last);
          qi = q_inter;
        }
        last = steer (qi, q2);
        if (!last) return PathVectorPtr_t ();
        pv->appendPath (last);
        PathVectorPtr_t out = PathVector::create
	  (path->outputSize (), path->outputDerivativeSize ());
        pv->flatten (out);
        return out;
      }

      JointVector_t PartialShortcut::generateJointVector
      (const PathVectorPtr_t& pv) const
      {
        const JointVector_t& rjv = problem().robot()->getJointVector ();
        JointVector_t jv;
	ConfigProjectorPtr_t proj;
	LockedJoints_t lj;
	if (pv->pathAtRank (0)->constraints ()){
	  proj = pv->pathAtRank (0)->constraints ()->configProjector ();
	  if (proj) lj = proj->lockedJoints ();
	}

        for (JointVector_t::const_iterator it = rjv.begin ();
	     it != rjv.end (); ++it) {
          if ((*it)->numberDof () > 0) {
            bool lock = false;
            if (parameters.removeLockedJoints && proj) {
              const std::size_t rkCfg = (*it)->rankInConfiguration ();
              for (LockedJoints_t::const_iterator itLJ = lj.begin ();
		   itLJ != lj.end (); ++itLJ) {
                if ((*itLJ)->rankInConfiguration () == rkCfg) {
                  lock = true;
                  break;
                }
              }
            }
            if (!lock) jv.push_back (*it);
          }
        }
        return jv;
      }

      PathVectorPtr_t PartialShortcut::optimizeFullPath
      (const PathVectorPtr_t& pv, const JointVector_t& jvIn,
       JointVector_t& jvOut) const
      {
        Configuration_t q0 = pv->initial ();
        Configuration_t q3 = pv->end ();
        const value_type t0 = 0;
        value_type t3;
        PathVectorPtr_t opted = pv;

        /// First try to optimize each joint from beginning to end
        for (std::size_t iJ = 0; iJ < jvIn.size(); ++iJ) {
          t3 = opted->timeRange ().second;
          JointPtr_t joint = jvIn[iJ];

          // Validate sub parts
          bool valid;
          PathVectorPtr_t straight;
          straight = generatePath (opted, joint, t0, q0, t3, q3);
          {
            PathPtr_t validPart;
            PathValidationReportPtr_t report;
            if (!straight) valid = false;
            else {
              valid = problem ().pathValidation ()->validate
                (straight, false, validPart, report);
            }
          }
          if (!valid) {
            jvOut.push_back (joint);
            continue;
          }
          opted = straight;

          hppDout (info, "length = "
		   << pathLength (opted, problem ().distance ())
		   << ", joint " << joint->name());
        }
        return opted;
      }

      PathVectorPtr_t PartialShortcut::optimizeRandom
      (const PathVectorPtr_t& pv, const JointVector_t& jv) const
      {
	//srand (time (NULL));
	struct timespec start, now;
	clock_gettime(CLOCK_MONOTONIC, &start);
	std::size_t iVec = 0; // index for timeValues and gainValues vectors
	value_type currentTime;
        PathVectorPtr_t current = pv, result = pv;
        const value_type t0 = 0;
        value_type t3;
        Configuration_t q0 = pv->initial ();
        Configuration_t q3 = pv->end ();
	Configuration_t q1 (pv->outputSize ()), q2 (pv->outputSize ());
        value_type length = pathLength (pv, problem ().distance ()),
	  newLength = std::numeric_limits <value_type>::infinity ();
	hppDout (info, "length (after fullOptimize) = " << length);
        hppDout (info, "random partial shorcut on " << jv.size ()
		 << " joints.");

	// Fill first value after fullOptimize
	/*clock_gettime(CLOCK_MONOTONIC, &now);
	currentTime = now.tv_sec - start_.tv_sec +
	  (now.tv_nsec - start_.tv_nsec) / BILLION;
	problem ().timeValues_.resize (1);
	problem ().gainValues_.resize (1);
	problem ().timeValues_ [0] = currentTime;
	problem ().gainValues_ [0] = (initialLength_ - length)/initialLength_;
	iVec++;*/

        std::size_t iJ;
	bool finished = false;
        while (!finished) {
	  // select random joint
	  iJ = floor (((value_type) jv.size()) * rand ()/RAND_MAX); 
          JointPtr_t joint = jv[iJ];
	  hppDout (info, "shorcut on joint: " << joint->name ()
		   << ", index: " << iJ);

          t3 = current->timeRange ().second;
          value_type u2 = t3 * rand ()/RAND_MAX;
          value_type u1 = t3 * rand ()/RAND_MAX;

          value_type t1, t2;
          if (u1 < u2) {t1 = u1; t2 = u2;} else {t1 = u2; t2 = u1;}
	  bool success = (*current) (q1, t1) &&(*current) (q2, t2);
	  // Warning: Constraints satisfaction (success) NOT checked
          // Generate sub parts
          bool valid [3];
          PathVectorPtr_t straight [3];
          straight [0] = generatePath (current, joint, t0, q0, t1, q1);
          straight [1] = generatePath (current, joint, t1, q1, t2, q2);
          straight [2] = generatePath (current, joint, t2, q2, t3, q3);
	  // Validate sub parts
          for (unsigned i=0; i<3; ++i) {
            PathPtr_t validPart;
            PathValidationReportPtr_t report;
            if (!straight [i]) valid[i] = false;
            else {
              valid [i] = problem ().pathValidation ()->validate
                (straight [i], false, validPart, report);
            }
          }
          // Replace valid parts
          result = PathVector::create (pv->outputSize (),
				       pv->outputDerivativeSize ());
          if (valid [0])
            result->concatenate (*straight [0]);
          else
            result->concatenate (*(current->extract
				   (std::make_pair (t0, t1))->
				   as <PathVector> ()));
          if (valid [1])
            result->concatenate (*straight [1]);
          else
            result->concatenate (*(current->extract
				   (std::make_pair (t1, t2))->
				   as <PathVector> ()));
          if (valid [2])
            result->concatenate (*straight [2]);
          else
            result->concatenate (*(current->extract
				   (std::make_pair (t2, t3))->
				   as <PathVector> ()));
	  // Check length improvement
          newLength = pathLength (result, problem ().distance ());
          
          hppDout (info, "previous length = " << length
		   << ", new length = " << newLength);
			  
	  clock_gettime(CLOCK_MONOTONIC, &now);
	  currentTime = now.tv_sec - start.tv_sec +
	    (now.tv_nsec - start.tv_nsec) / BILLION;
	  hppDout (info, "currentTime = " << currentTime);
	  if (newLength < length) {
	    problem ().timeValues_.resize (iVec+1);
	    problem ().gainValues_.resize (iVec+1);
	    problem ().timeValues_ [iVec] = currentTime;
	    problem ().gainValues_ [iVec] = (initialLength_ -
					     newLength)/initialLength_;
	    iVec++;
	  }
	  length = newLength; // previousLength = currentLength;
	  current = result; // previous = current;
	  finished = currentTime > problem ().tGB_;
        }
        return result;
      }
    } // namespace pathOptimization
  } // namespace core
} // namespace hpp
