//
// Copyright (c) 2015 CNRS
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

//#include <limits>
//#include <deque>
//#include <cstdlib>
#include <hpp/util/assertion.hh>
#include <hpp/util/debug.hh>
#include <hpp/core/distance.hh>
#include <hpp/core/config-validations.hh>
#include <hpp/core/path-validation.hh>
#include <hpp/core/path-vector.hh>
#include <hpp/core/problem.hh>
#include <hpp/core/prune.hh>
#include <hpp/core/steering-method.hh>

namespace hpp {
  namespace core {
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

    // Compute a vector of the added path sub-paths lengths
    static vector_t addedPathLengths (const PathVectorPtr_t& path)
    {
      vector_t result (path->numberPaths () + 1);
      value_type length = 0;
      result (0) = 0;
      for (std::size_t i=0; i<path->numberPaths (); ++i) {
	const PathPtr_t& element (path->pathAtRank (i));
	length += element->length ();
	result (i+1) = length;
      }
      return result;
    }

    // Re-create the path-vector replacing the index part with the local path
    static PathVectorPtr_t replaceInPathVector (const PathVectorPtr_t& path,
						const PathPtr_t& lp,
						const std::size_t index)
    {
      PathVectorPtr_t result =
	PathVector::create (path->outputSize (), path->outputDerivativeSize ());
      for (std::size_t i=0; i<path->numberPaths (); ++i) {
	const PathPtr_t& element (path->pathAtRank (i));
	if (element->length () > 1e-6) {
	  if (i != index && i != index + 1) {
	    result->appendPath (element);
	  } else if (i == index) {
	    result->appendPath (lp);
	  }
	}
      }
      return result;
    }

    PrunePtr_t Prune::create (const Problem& problem)
    {
      Prune* ptr = new Prune (problem);
      return PrunePtr_t (ptr);
    }

    Prune::Prune (const Problem& problem) :
      PathOptimizer (problem)
    {
    }

    PathVectorPtr_t Prune::optimize (const PathVectorPtr_t& path)
    {
      using std::make_pair;
      const DistancePtr_t& distance = problem ().distance ();
      PathValidationPtr_t pathValidation (problem ().pathValidation ());
      const SteeringMethodPtr_t& sm (problem ().steeringMethod ());
      PathValidationReportPtr_t report;
      PathPtr_t validPart;
      size_type nbSubPaths = path->numberPaths ();
      std::size_t i = 0;
      vector_t pathLengths;
      PathVectorPtr_t result = PathVector::createCopy (*path);
      
      while (i < nbSubPaths -1) {
	pathLengths = addedPathLengths (result);
	const value_type t1 = pathLengths (i);
	const value_type t2 = pathLengths (i + 2);
	Configuration_t q1 = (*result) (t1);
	Configuration_t q2 = (*result) (t2);
	PathPtr_t lp = (*sm) (q1,q2);
	
	PathVectorPtr_t subPath = result->extract
	  (make_pair <value_type, value_type> (t1, t2))->
	  as <PathVector> ();
	value_type subPathLength = pathLength (subPath, distance);

	if (pathValidation->validate (lp, false, validPart, report) &&
	    lp->length () < subPathLength){
	  // Replace path portion
	  result = PathVector::createCopy (*replaceInPathVector (result, 
								 lp, i));
	  nbSubPaths -= 1;
	} else {
	  // Path portion in collision or not improving global path length
	  i += 1;
	}
      }
      return result;
    }

  } // namespace core
} // namespace hpp
