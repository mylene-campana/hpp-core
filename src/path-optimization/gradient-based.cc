//
// Copyright (c) 2014 CNRS
// Authors: Florent Lamiraux, Mylene Campana
//
// This file is part of hpp-core
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
// hpp-core  If not, see
// <http://www.gnu.org/licenses/>.

#include <hpp/fcl/collision_object.h>
#include <hpp/fcl/distance.h>
#include <hpp/fcl/collision.h>

#include <hpp/util/timer.hh>

#include <Eigen/SVD>
#include <Eigen/Dense>
#include <hpp/util/debug.hh>
#include <hpp/model/body.hh>
#include <hpp/model/joint.hh>
#include <hpp/model/configuration.hh>
#include <hpp/model/collision-object.hh>
#include <hpp/core/collision-path-validation-report.hh>
#include <hpp/core/config-projector.hh>
#include <hpp/core/config-validations.hh>
#include <hpp/core/constraint-set.hh>
#include <hpp/core/numerical-constraint.hh>
#include <hpp/core/path-optimization/gradient-based.hh>
#include <hpp/core/path-optimization/path-length.hh>
#include <hpp/core/path-validation.hh>
#include <hpp/core/problem.hh>
#include <hpp/constraints/generic-transformation.hh>
#include "path-optimization/collision-constraints-result.hh"

// For debug
#include <hpp/core/problem-solver.hh>

#define BILLION 1E9

namespace hpp {
  namespace core {
    namespace pathOptimization {
      namespace {
        HPP_DEFINE_TIMECOUNTER(GBO_computeIterate);
        HPP_DEFINE_TIMECOUNTER(GBO_oneStep);
      }

      // Compute the length of a vector of paths assuming that each element
      // is optimal for the given distance.
      static value_type pathLength (const PathVectorPtr_t& path,
				    const DistancePtr_t& distance)
      {
	value_type result = 0;
	for (std::size_t i=0; i<path->numberPaths (); ++i) {
	  const PathPtr_t& element (path->pathAtRank (i));
	  const value_type& tmin (element->timeRange ().first);
	  const value_type& tmax (element->timeRange ().second);
	  Configuration_t q1 ((*element) (tmin));
	  Configuration_t q2 ((*element) (tmax));
	  result += (*distance) (q1, q2);
	}
	return result;
      }

      using model::displayConfig;
      GradientBasedPtr_t GradientBased::create
      (const Problem& problem)
      {
	GradientBased* ptr = new GradientBased (problem);
	GradientBasedPtr_t shPtr (ptr);
	return shPtr;
      }

      GradientBased::GradientBased (const Problem& problem) :
	PathOptimizer (problem), cost_ (), robot_ (problem.robot ()),
	configSize_ (robot_->configSize ()), robotNumberDofs_
	(robot_->numberDof ()),	robotNbNonLockedDofs_ (robot_->numberDof ()),
	fSize_ (1),
	initial_ (), end_ (), epsilon_ (1e-3), iterMax_ (30),
	alphaInit_ (problem.alphaInit_), alphaMax_ (1.)
      {
	distance_ = HPP_DYNAMIC_PTR_CAST (WeighedDistance, problem.distance ());
	if (!distance_) {
	  throw std::runtime_error
	    ("Distance is not of type WeighedDistance.");
	}
	HPP_STATIC_CAST_REF_CHECK (SteeringMethodStraight,
				   *(problem.steeringMethod ()));
	steeringMethod_ = HPP_DYNAMIC_PTR_CAST
	  (SteeringMethodStraight, problem.steeringMethod ()->copy ());
        // Unset the constraints.
        steeringMethod_->constraints (ConstraintSetPtr_t());
      }

      void GradientBased::compressHessian (matrixIn_t normal, matrixOut_t small)
	const
      {
	ConstraintSetPtr_t constraints (problem ().constraints ());
	if (!constraints) return;
	size_type normalSize = robotNumberDofs_;
	size_type smallSize = constraints->numberNonLockedDof ();
	for (size_type i=0; i<nbWaypoints_; ++i) {
	  for (size_type j=0; j<nbWaypoints_; ++j) {
	    constraints->compressMatrix
	      (normal.block (i*normalSize, j*normalSize, normalSize,normalSize),
	       small.block (i*smallSize, j*smallSize, smallSize, smallSize));
	  }
	}
      }

      void GradientBased::compressVector (vectorIn_t normal,
					  vectorOut_t small) const
      {
	ConstraintSetPtr_t constraints (problem ().constraints ());
	if (!constraints) {
	  small = normal;
	  return;
	}
	size_type normalSize = robotNumberDofs_;
	size_type smallSize = constraints->numberNonLockedDof ();
	for (size_type i=0; i<nbWaypoints_; ++i) {
	  constraints->compressVector
	    (normal.segment (i*normalSize, normalSize),
	     small.segment (i*smallSize, smallSize));
	}
      }

      void GradientBased::uncompressVector (vectorIn_t small,
					    vectorOut_t normal) const
      {
	ConstraintSetPtr_t constraints (problem ().constraints ());
	if (!constraints) {
	  normal = small;
	  return;
	}
	size_type normalSize = robotNumberDofs_;
	size_type smallSize = constraints->numberNonLockedDof ();
	for (size_type i=0; i<nbWaypoints_; ++i) {
	  constraints->uncompressVector
	    (small.segment (i*smallSize, smallSize),
	     normal.segment (i*normalSize, normalSize));
	}
      }

      void GradientBased::initialize (const PathVectorPtr_t& path)
      {
	nbWaypoints_ = path->numberPaths () -1;
	hppDout (info, "nbWaypoints_ = " << nbWaypoints_);
	if (nbWaypoints_ == 0) return;
	vector_t onesVector (path->numberPaths ());
	for (int i = 0; i < onesVector.size (); i++) onesVector [i] = 1;
	/* Create cost */
	if (!cost_ || cost_->inputSize () !=
	    nbWaypoints_ * path->outputSize () ||
	    cost_->inputDerivativeSize () !=
	    nbWaypoints_ * path->outputDerivativeSize ())
	  {
	    hppDout (info, "creating cost");
	    cost_ = PathLength::create (distance_, path);
	    //cost_->setLambda (onesVector);
	  }
	hppDout (info, "size of x = " << cost_->inputSize ());
	numberDofs_ = cost_->inputDerivativeSize ();
	stepNormal_.resize (numberDofs_);
	stepNormal_.setZero ();
	p_.resize (numberDofs_);
	// Get problem constraints and locked degrees of freedom
	initializeProblemConstraints ();

	/* Create Hessian and inverse Hessian of cost */
	matrix_t H;
	H.resize (cost_->inputDerivativeSize (),
		  cost_->inputDerivativeSize ());
	cost_->hessian (H);
	H_.resize (numberDofs_, numberDofs_);
	rgrad_.resize (1, numberDofs_);
	compressHessian (H, H_);
	/* Store first and last way points */
	initial_ = path->initial ();
	end_ = path->end ();
	alphaInit_ = problem ().alphaInit_;
	alpha_ = alphaInit_;
      }

      vector_t GradientBased::computeIterate (vectorIn_t) const
      {
	if (J_.rows () == 0) {
	  // no constraints
	  //   grad (x)^T  = grad (x1)^T + H * (x-x1)
	  // - grad (x1)^T = H * (x-x1)
	  //      x-x1   = - Hinverse * grad (x1)^T
	  vector_t result (-alpha_ * Hinverse_ * rgrad_.transpose ());
	  return result;
	} else {
	  // cost
	  // p = x - x1
	  // C (x) = C (x1) + grad (x1) p + 1/2 p^T H p
	  // constraints
	  // f (x) = f0
	  // linearized
	  // f (x) = f (x1) + J p = rhs
	  // J p = rhs - f (x1)
	  // svd (J) -> V0
	  // p = p0 + V0 z
	  // with p0 = J^{+} (rhs - f(x1))
	  // C (x) = C (x1) + grad (x1) * (p0 + V0 z)
	  //         + 1/2 (p0 + V0 z)^T H (p0 + V0 z)
	  // C (x) = (grad (x1) + p0^T H) V0 z  + 1/2 z^T V0^T H V0 z + C (x1)
	  // z = - (V0^T H V0)^{-1} V0^T(grad (x1)^T + H p0)
          hppDout (info, "J_\n" << J_);
	  Jacobi_t svd (J_, Eigen::ComputeThinU | Eigen::ComputeFullV);
	  size_type rank = svd.rank();
	  hppDout (info, "J_ singular values = " <<
		   svd.singularValues ().transpose ());
	  hppDout (info, "Jrows = " << J_.rows ());
	  hppDout (info, "rank(J) = " << rank);
	  if (rank < J_.rows ()) {
	    p_.setZero ();
	    hppDout (error, "rank loss");
	    return p_;
	  }
	  hppDout (info, "rhs_ - value_ = " << rhs_ - value_);
	  p0_ = svd.solve (rhs_ - value_);
          const size_type nullity = numberDofs_ - rank; 
	  hppDout (info, "nullity = " << nullity);
	  if (nullity != 0) {
            const Jacobi_t::MatrixVType& V = svd.matrixV();
            // V0_ orthogonal base of Jf_ null space
            // V0_ = svd.matrixV().rightCols (numberDofs_-rank);
            // Hz_ = V0_.transpose () * H_ * V0_;
            Hz_ = V.rightCols(nullity).transpose() * H_ * V.rightCols(nullity);
            // gz_ = - V0_.transpose () * (rgrad_.transpose () + H_ * p0_);
            gz_ = - V.rightCols(nullity).transpose() * (rgrad_.transpose () + H_ * p0_);
	    Jacobi_t svd2 (Hz_, Eigen::ComputeThinU | Eigen::ComputeFullV);
	    /*rank = svd2.rank ();
	    hppDout (info, "Hz_ singular values = " <<
		     svd2.singularValues ().transpose ());
		     hppDout (info, "rank(Hz_) = " << rank);*/
	    vector_t z (svd2.solve (gz_));
	    hppDout (info, "Hz_ * z - gz_=" << (Hz_ * z - gz_).transpose ());
	    // p_ = p0_ + V0_ * z;
	    p_ = p0_ + V.rightCols(nullity) * z;
	    hppDout (info, "constraint satisfaction: " <<
		     (J_*p_ - (rhs_ - value_)).squaredNorm ());
	    hppDout (info, "norm(grad*V0)? = " <<
		     ((rgrad_ + (H_*p_).transpose ())*V.rightCols(nullity)).squaredNorm ());
	  } else {
	    p_ = p0_;
	  }
	  hppDout (info, "p_: " << p_.transpose ());
	  return alpha_*p_;
	}
      }

      typedef std::vector <std::pair <CollisionPathValidationReportPtr_t,
				      std::size_t> > Reports_t;

      // avoid adding the same collision twice in reports
      // because the subpath 2 'begins' in the obstacle
      static bool collisionRedundancy (const std::size_t& i1,
				       const std::size_t& i2,
				       const value_type& param1,
				       const value_type& param2) {
	if (i1 == i2 - 1 && param1 == 1. && param2 == 0.)
	  return true;
	return false;
      }

      bool validatePath (const PathValidationPtr_t& pathValidation,
			 const PathVectorPtr_t& path, Reports_t& reports)
      {
	bool valid = true;
	PathPtr_t validPart;
	PathValidationReportPtr_t report;
	reports.clear ();
	bool redundantColl = false;
	for (std::size_t i=0; i<path->numberPaths (); ++i) {
	  if (!pathValidation->validate
	      (path->pathAtRank (i), false, validPart, report)) {
	    if (reports.size () != 0) {
	      value_type tmpL = path->pathAtRank
		(reports.back ().second)->length ();
	      redundantColl =
		collisionRedundancy (reports.back ().second, i,
				     reports.back ().first->parameter/tmpL,
				     report->parameter);
	    }
	    if (!redundantColl) {
	      HPP_STATIC_CAST_REF_CHECK (CollisionPathValidationReport,
					 *report);
	      reports.push_back
		(std::make_pair (HPP_STATIC_PTR_CAST
				 (CollisionPathValidationReport, report), i));
	    } else
	      hppDout (info, "collision non added because already reported");
	    redundantColl = false;
	    valid = false;
	  }
	}
	return valid;
      }

      PathVectorPtr_t GradientBased::optimize (const PathVectorPtr_t& path)
      {
	problem ().timeValues_.clear ();
	problem ().gainValues_.clear ();
	struct timespec start, now, end;
	clock_gettime(CLOCK_MONOTONIC, &start);
	std::size_t iVec = 0; // index for timeValues and gainValues vectors
	const value_type initialLength = pathLength (path,
						     problem ().distance ());
	value_type currentTime, currentLength, previousLength = initialLength;
	interrupt_ = false;
	initialize (path);
	if (nbWaypoints_ == 0) // path is direct and optimal
	  return path;
	PathValidationPtr_t pathValidation (problem ().pathValidation ());
	PathVectorPtr_t result = PathVector::create (configSize_,
						     robotNumberDofs_);
	ConstraintSetPtr_t constraints (problem ().constraints ());
	Configuration_t qCollConstr, qCollConstr_old;
	ConfigValidationsPtr_t configValtions (problem ().configValidations());
	Reports_t reports, newReports;
	bool noCollision = true; // Assume there is no collision
	bool compute_iterate = true;
	/* Create initial path */
	vector_t x1, s; x1.resize (cost_->inputSize ());
	rowvector_t grad; grad.resize (cost_->inputDerivativeSize ());
	pathToVector (path, x1);
	vector_t x0 = x1;
	vector_t x = x1;
	Hinverse_ = H_.inverse ();
	//hppDout (info, "Hessian = " << H_);
	//hppDout (info, "Hinverse_ = " << Hinverse_);
	Jacobi_t svdHinv (Hinverse_, Eigen::ComputeThinU | Eigen::ComputeThinV);
	hppDout (info, "Hinv rows: " << svdHinv.rows ());
	hppDout (info, "Hinv rank: " << svdHinv.rank ());

	/* Fill jacobian J_ and Jf_ with constraints FROM Problem */
	cost_->jacobian (grad, x1);
	hppDout (info, "init grad: " << grad);
	compressVector (grad.transpose (), rgrad_.transpose ());
	bool minimumReached = (rgrad_.squaredNorm () <= epsilon_);
	bool isConstraintRedundant;

	PathVectorPtr_t path0 (path);
	CollisionConstraintsResults_t collisionConstraints;
	updateProblemConstraints (x0);
        while (!(noCollision && minimumReached) && (!interrupt_)) {
          HPP_START_TIMECOUNTER(GBO_oneStep);
	  if (compute_iterate) {
	    HPP_START_TIMECOUNTER(GBO_computeIterate);
	    s = computeIterate (x0);
	    hppDout (info, "s: " << s.transpose ());
	    HPP_STOP_TIMECOUNTER(GBO_computeIterate);
	    HPP_DISPLAY_TIMECOUNTER(GBO_computeIterate);
	  }
	  else
	    s = alpha_ * p_; // re-use optimal p_
          hppDout (info, "alpha_ = " << alpha_);
          value_type norm_s (s.norm ());
          minimumReached = (norm_s < epsilon_ || alpha_ == 1.);
          hppDout (info, "norm of s: " << norm_s);
          hppDout (info, "s: " << s.transpose ());
          if (norm_s == 0) {
            HPP_STOP_TIMECOUNTER(GBO_oneStep);
            return path0;
          }
          integrate (x0, s, x1);
          hppDout (info, "x1=" << x1.transpose ());
          PathVectorPtr_t path1 = PathVector::create (configSize_,
              robotNumberDofs_);
          vectorToPath (x1, path1);
          if (!solveConstraints (x1, path1, collisionConstraints)) {
            hppDout (info, "solveConstraints failed");
            alpha_ /= 2;
            HPP_STOP_TIMECOUNTER(GBO_oneStep);
            continue;
          }
          bool isPathValid = validatePath (pathValidation, path1, reports);
          // if new path is in collision, we add some constraints
          if (!isPathValid) {
            if (alpha_ != 1.) {
              for (Reports_t::const_iterator it = reports.begin ();
                  it != reports.end (); ++it) {
                CollisionConstraintsResult ccr (robot_, path0, path1,
                    *it, J_.rows (),
                    robotNbNonLockedDofs_);
                ConstraintSetPtr_t constraints (problem ().constraints ());
                if (constraints) {
                  ConfigProjectorPtr_t configProjector =
                    constraints->configProjector ();
                  if (configProjector)
                    ccr.add (configProjector->lockedJoints());
                }
		// Linearize ccr around path0 and add it to J_
                addCollisionConstraint (ccr, path0);

		// check redundancy (rank loss in J_)
		Jacobi_t svd (J_, Eigen::ComputeThinU | Eigen::ComputeFullV);
		size_type rank = svd.rank();
		hppDout (info, "Jrows = " << J_.rows ());
		hppDout (info, "rank(J) = " << rank);
		isConstraintRedundant = rank < J_.rows ();
		
		// Avoid redundancy solving
		if (isConstraintRedundant) {
		  hppDout(error, "redundant constraint");
		  clock_gettime(CLOCK_MONOTONIC, &end);
		  problem ().tGB_ = end.tv_sec - start.tv_sec +
		    (end.tv_nsec - start.tv_nsec) / BILLION;
		  hppDout(info, "tGB_: " << problem ().tGB_);
		  return path0;
		}

		// Solve redundancy if necessary (finding a new constraint)
		/*while (isConstraintRedundant) {
		  hppDout (info, "redundancy found");
		  J_.conservativeResize(J_.rows ()-1, J_.cols ());
		  // execute step while halving alpha
		  alpha_ *= 0.5;
		  hppDout (info, "alpha_ = " << alpha_);
		  integrate (x0, alpha_*p_, x);
		  PathVectorPtr_t pathTmp =
		    PathVector::create (configSize_, robotNumberDofs_);
		  vectorToPath (x, pathTmp);
		  // findNewConstraint with temporary path
		  if (validatePath (pathValidation, pathTmp, newReports)) {
		    // compute and linearize new constraint around temp. path
		    hppDout (info, "temporary path collision-free");
		    path0 = pathTmp;
		    x0 = x;
		    CollisionConstraintsResult ccr (robot_, path0, path1,
						    *(reports.begin ()),
						    J_.rows (),
						    robotNbNonLockedDofs_);
		  } else {
		    // compute new constraint with temp. path collision
		    hppDout (info, "intermediaite step in collision");
		    path1 = pathTmp;
		    CollisionConstraintsResult ccr (robot_, path0, path1,
						    *(newReports.begin ()),
						    J_.rows (),
						    robotNbNonLockedDofs_);
		  }
		  addCollisionConstraint (ccr, path0);
		  // check if new constraint solves redundancy
		  Jacobi_t svd (J_, Eigen::ComputeThinU | Eigen::ComputeFullV);
		  rank = svd.rank();
		  isConstraintRedundant = rank < J_.rows ();
		  hppDout (info, "Jrows = " << J_.rows ());
		  hppDout (info, "rank(J) = " << rank);
		}*/


                collisionConstraints.push_back (ccr);
                hppDout (info, "Number of collision constraints: "
                    << collisionConstraints.size ());
                // When adding a new constraint, try first minimum under this
                // constraint. If this latter minimum is in collision,
                // re-initialize alpha_ to alphaInit_.
                alpha_ = 1.;
		compute_iterate = true;
              }
            } else {
              alpha_ = alphaInit_;
	      compute_iterate = false;
            }
            noCollision = false;
          } else { // path valid
            rowvector_t rgrad0 (rgrad_);
            x0 = x1;
            path0 = path1;
            cost_->jacobian (grad, x0);
            compressVector (grad.transpose (), rgrad_.transpose ());
            noCollision = true;
	    compute_iterate = true;
          }
          HPP_STOP_TIMECOUNTER(GBO_oneStep);
          HPP_DISPLAY_TIMECOUNTER(GBO_oneStep);
	  // gather result data
	  clock_gettime(CLOCK_MONOTONIC, &now);
	  currentTime = now.tv_sec - start.tv_sec +
	    (now.tv_nsec - start.tv_nsec) / BILLION;
	  hppDout (info, "currentTime= " << currentTime);
	  currentLength = pathLength (path0, problem ().distance ());
	  hppDout (info, "currentLength= " << currentLength);
	  if (currentLength < previousLength) {
	    problem ().timeValues_.resize (iVec+1);
	    problem ().gainValues_.resize (iVec+1);
	    problem ().timeValues_ [iVec] = currentTime;
	    problem ().gainValues_ [iVec] = (initialLength -
					     currentLength)/initialLength;
	    iVec++;
	  }
	  previousLength = currentLength;
        } // while (!(noCollision && minimumReached) && (!interrupt_))
	clock_gettime(CLOCK_MONOTONIC, &end);
	const value_type endTime = end.tv_sec - start.tv_sec +
	  (end.tv_nsec - start.tv_nsec) / BILLION;
	problem ().tGB_ = endTime;
	hppDout(info, "tGB_: " << endTime);
	return path0;
      }

      void GradientBased::integrate (vectorIn_t x0, vectorIn_t step,
				     vectorOut_t x1) const
      {
	assert (x0.size () == x1.size ());
	size_type indexConfig = 0;
	size_type indexVelocity = 0;
	// uncompress step
	uncompressVector (step, stepNormal_);
	while (indexConfig < x0.size ()) {
	  model::integrate (robot_, x0.segment (indexConfig, configSize_),
			    stepNormal_.segment
			    (indexVelocity, robotNumberDofs_),
			    x1.segment (indexConfig, configSize_));
	  indexConfig += configSize_;
	  indexVelocity += robotNumberDofs_;
	}
	assert (indexConfig == x0.size ());
	assert (indexVelocity == stepNormal_.size ());
      }

      void GradientBased::initializeProblemConstraints ()
      {
	ConstraintSetPtr_t constraints (problem ().constraints ());
	if (!constraints){
	  J_.resize (0, numberDofs_);
	  return;
	}
	robotNbNonLockedDofs_ = constraints->numberNonLockedDof ();
	numberDofs_ = constraints->numberNonLockedDof () * nbWaypoints_;
	ConfigProjectorPtr_t configProjector (constraints->configProjector ());
	size_type rows = 0;
	if (configProjector) {
	  // Apply problem constraints to each waypoint.
	  size_type sizeConstraint = configProjector->dimension ();
	  rows = sizeConstraint * nbWaypoints_;
	  rhs_.resize (rows); rhs_.setZero ();
	}
	value_.resize (rows);
	J_.resize (rows, numberDofs_); J_.setZero ();
      }

      /// Compute value and Jacobian of the problem constraints
      ///
      /// \param x input path as a vector
      ///
      /// Constraints of the problem are applied to each waypoint. Therefore,
      /// the first 6n lines of the Jacobian and of the value are computed,
      /// where n is the number of waypoints.
      void GradientBased::updateProblemConstraints (vectorIn_t x)
      {
	ConstraintSetPtr_t constraints (problem ().constraints ());
	if (!constraints){
	  return;
	}
	ConfigProjectorPtr_t configProjector (constraints->configProjector ());
	if (configProjector) {
	  // Apply problem constraints to each waypoint.
          size_type sizeConstraint = configProjector->dimension ();
	  for (size_type i=0; i<nbWaypoints_; ++i) {
	    configProjector->computeValueAndJacobian
	      (x.segment (i*configSize_, configSize_),
	       value_.segment (i*sizeConstraint, sizeConstraint),
	       J_.block (i*sizeConstraint, i*constraints->numberNonLockedDof (),
			 sizeConstraint, constraints->numberNonLockedDof ()));
	  }
	}
      }

      /// Check whether constraints are satisfied
      ///
      /// \param x input path as a vector,
      /// \param path input path (same path as x)
      /// \param collisionConstraints vector of collision constraints
      ///
      /// Check that
      ///   \li each waypoint satisfies the problem constraints, and
      ///   \li each collision constraint drift is within bounds that ensures
      ///       the user that there is no collision
      ///
      bool GradientBased::constraintsSatisfied
      (vectorOut_t& x, PathVectorPtr_t&,
       CollisionConstraintsResults_t&)
      {
	bool satisfied = true;
	if (problem ().constraints ()) {
	  for (size_type i=0; i<nbWaypoints_; ++i) {
	    if (!problem ().constraints ()->isSatisfied
		(x.segment (i*configSize_, configSize_))) {
	      satisfied = false;
	    }
	  }
	}
	return satisfied;
      }

      /// Solve problem constraints and collision constraints
      ///
      /// \retval x path that satisfies the constraints as a vector,
      /// \retval path path that satisfies the constraints (same as x)
      /// \param collisionConstraints vector of collision constraints
      ///
      /// While constraints are not satisfied
      /// (method GradientBased::constraintsSatisfied), linearize constraints
      /// (problem and collision) around current path and solve linearized
      /// constraints regardless of the cost
      bool GradientBased::solveConstraints
      (vectorOut_t x, PathVectorPtr_t& path,
       CollisionConstraintsResults_t& collisionConstraints)
      {
	size_type iter = 0;
	while (!constraintsSatisfied (x, path, collisionConstraints)) {
	  vector_t prevX_ = x;
	  updateProblemConstraints (x);
	  for (CollisionConstraintsResults_t::iterator it =
		 collisionConstraints.begin (); it !=
		 collisionConstraints.end (); ++it) {
	    if (!it->linearize (path, J_, value_)) return false;
	  }
	  Eigen::JacobiSVD <matrix_t> svd (J_,
					   Eigen::ComputeThinU |
					   Eigen::ComputeThinV);
	  svd.setThreshold (1e-3);
	  vector_t dx = svd.solve(rhs_ - value_);
	  //hppDout (info, "rhs_ - value = " << (rhs_ - value_).transpose ());
	  vector_t x1 (x);
	  integrate (x1, dx, x);
	  path = PathVector::create (configSize_, robotNumberDofs_);
	  vectorToPath (x, path);
	  ++ iter;
	  if (iter > 100) {
	    return false;
	  }
	}
	updateProblemConstraints (x);
	for (CollisionConstraintsResults_t::iterator it =
	       collisionConstraints.begin (); it !=
	       collisionConstraints.end (); ++it) {
	  if (!it->linearize (path, J_, value_)) return false;
	}
	return true;
      }

      void GradientBased::addCollisionConstraint
	(const CollisionConstraintsResult& ccr, const PathVectorPtr_t& path)
	const
      {
	size_type Jrows = J_.rows ();
	J_.conservativeResize (Jrows + fSize_, J_.cols ());
	J_.middleRows (ccr.rowInJacobian (), fSize_).setZero ();
	value_.conservativeResize (Jrows + fSize_);
	value_.segment (ccr.rowInJacobian (), fSize_).setZero ();
	rhs_.conservativeResize (Jrows + fSize_);
	rhs_ [ccr.rowInJacobian ()] = 0;
	ccr.linearize (path, J_, value_);
      }

      void GradientBased::updateRightHandSide
      (const CollisionConstraintsResults_t& collisionConstraints,
       const PathVectorPtr_t& path) const
      {
	for (CollisionConstraintsResults_t::const_iterator it =
	       collisionConstraints.begin (); it != collisionConstraints.end ();
	     ++it) {
	  it->updateRightHandSide (path, rhs_);
	}

      }

    } // namespace pathOptimization
  }  // namespace core
} // namespace hpp
