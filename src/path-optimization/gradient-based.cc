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
#include <hpp/constraints/differentiable-function.hh>
#include <hpp/constraints/position.hh>
#include <hpp/constraints/relative-position.hh>
#include <hpp/constraints/orientation.hh>
#include <hpp/constraints/relative-orientation.hh>
#include <hpp/constraints/transformation.hh>
#include <hpp/constraints/relative-transformation.hh>
#include "path-optimization/collision-constraints-result.hh"

// For debug
#include <hpp/core/problem-solver.hh>

namespace hpp {
  namespace core {
    namespace pathOptimization {
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
	fSize_ (6),
	initial_ (), end_ (), epsilon_ (1e-3), iterMax_ (30), alphaInit_ (0.02),
	alphaMax_ (1.)
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
	/* Create cost */
	if (!cost_ || cost_->inputSize () !=
	    nbWaypoints_ * path->outputSize () ||
	    cost_->inputDerivativeSize () !=
	    nbWaypoints_ * path->outputDerivativeSize ())
	  {
	    hppDout (info, "creating cost");
	    cost_ = PathLength::create (distance_, path);
	  }
	hppDout (info, "size of x = " << cost_->inputSize ());
	numberDofs_ = cost_->inputDerivativeSize ();
	stepNormal_.resize (numberDofs_);
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
	alpha_ = alphaInit_;
      }

      vector_t GradientBased::computeIterate () const
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
	  Jacobi_t svd (J_, Eigen::ComputeThinU | Eigen::ComputeFullV);
	  svd.setThreshold (1e-6);
	  size_type rank = svd.rank();
	  hppDout (info, "J_ singular values = " <<svd.singularValues ().transpose ());
	  hppDout (info, "rank(J) = " << rank);
	  V0_ = svd.matrixV ().rightCols (numberDofs_-rank);
	  p0_ = svd.solve (rhs_ - value_);
	  if (V0_.cols () != 0) {
	    Hz_ = V0_.transpose () * H_ * V0_;
	    gz_ = - V0_.transpose () * (rgrad_.transpose () + H_ * p0_);
	    Jacobi_t svd2 (Hz_, Eigen::ComputeThinU | Eigen::ComputeFullV);
	    rank = svd2.rank ();
	    hppDout (info, "Hz_ singular values = " <<svd2.singularValues ().transpose ());
	    hppDout (info, "rank(Hz_) = " << rank);
	    p_ = p0_ + V0_ * svd2.solve (gz_);
	    hppDout (info, "constraint satisfaction: " <<(J_*p_ - (rhs_ - value_)).squaredNorm ());
	    hppDout (info, "norm(grad*V0)? = " <<((rgrad_ + (H_*p_).transpose ())*V0_).squaredNorm ());
	  } else {
	    p_ = p0_;
	  }
	  return alpha_*p_;
	}
      }

      typedef std::vector <std::pair <CollisionPathValidationReportPtr_t,
				      std::size_t> > Reports_t;

      bool validatePath (const PathValidationPtr_t& pathValidation,
			 const PathVectorPtr_t& path, Reports_t& reports)
      {
	bool valid = true;
	PathPtr_t validPart;
	PathValidationReportPtr_t report;
	reports.clear ();
	for (std::size_t i=0; i<path->numberPaths (); ++i) {
	  if (!pathValidation->validate
	      (path->pathAtRank (i), false, validPart, report)) {
	    HPP_STATIC_CAST_REF_CHECK (CollisionPathValidationReport, *report);
	    reports.push_back
	      (std::make_pair (HPP_STATIC_PTR_CAST
			       (CollisionPathValidationReport, report), i));
	    valid = false;
	  }
	}
	return valid;
      }

      PathVectorPtr_t GradientBased::optimize (const PathVectorPtr_t& path)
      {
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
	Reports_t reports;
	bool noCollision;
	bool validConstraints = false;
	/* Create initial path */
	vector_t x1; x1.resize (cost_->inputSize ());
	rowvector_t grad; grad.resize (cost_->inputDerivativeSize ());
	pathToVector (path, x1);
	vector_t x0 = x1;
	Hinverse_ = H_.inverse ();
	hppDout (info, "inverse Hessian = " << Hinverse_);
	Jacobi_t svdH (Hinverse_, Eigen::ComputeThinU | Eigen::ComputeFullV);
	hppDout (info, "inverse Hessian rank = " << svdH.rank ());
	hppDout (info, "inv Hessian singular values = "
		 << svdH.singularValues ().transpose ());
	vector_t cost (1);
	(*cost_) (cost, x0);
	hppDout (info, "initial cost: " << cost);

	/* Fill jacobian J_ and Jf_ with constraints FROM Problem */
	cost_->jacobian (grad, x1);
	compressVector (grad.transpose (), rgrad_.transpose ());
	bool minimumReached = (rgrad_.squaredNorm () <= epsilon_);

	PathVectorPtr_t path0 (path);
	CollisionConstraintsResults_t collisionConstraints;
	size_type pathId = 0;

	alpha_ = 1; // ROTATION PROBLEM ONLY !
	if (ProblemSolver::latest ()) {
	  pathId = ProblemSolver::latest ()->paths ().size ();
	}
	if (!minimumReached) {
	  do {
	    // Compute value and Jacobian in x0
	    updateProblemConstraints (x0);
	    vector_t s = computeIterate ();
	    hppDout (info, "alpha_ = " << alpha_);
	    minimumReached = (s.norm () < epsilon_ || alpha_ == 1.);
	    //minimumReached = s.norm () < epsilon_;
	    hppDout (info, "norm of s: " << s.norm ());
	    hppDout (info, "s: " << s.transpose ());
	    integrate (x0, s, x1);
	    hppDout (info, "x1=" << x1.transpose ());
	    PathVectorPtr_t path1 = PathVector::create (configSize_,
							robotNumberDofs_);
	    vectorToPath (x1, path1);
	    if (ProblemSolver::latest ()) {
	      hppDout (info, "applied iteration to lastest valid path:" " path id = " << pathId);
	      ProblemSolver::latest ()->addPath (path1);
	      ++pathId;
	    }

	    // Solve problem constraints and linearize collision constraints
	    // around new x1, if alpha != 1
	    // do not increase latest valid path either if alpha = 1 ?
	    if (alpha_ != 1.) {
	      if (!solveConstraints (x1, path1, collisionConstraints)) {
		if (ProblemSolver::latest ()) {
		  hppDout (info, "failed to apply constraints: path id = " << pathId);
		}
		validConstraints = false;
		alpha_ *= .5;
		hppDout (info, "Halve alpha: " << alpha_);
	      } else { // (!solveConstraints (x1, path1, collisionConstraints))
		hppDout (info, "x1=" << x1.transpose ());
		validConstraints = true;
		if (ProblemSolver::latest ()) {
		  hppDout (info, "successfully applied constraints: path id = " << pathId);
		}
	      }
	    }
	    if (ProblemSolver::latest ()) {
	      ProblemSolver::latest ()->addPath (path1);
	      ++pathId;
	    }
	    bool isPathValid = validatePath (pathValidation, path1, reports);
	    // if new path is in collision and alpha != 1,
	    // we add some constraints
	    if (!isPathValid) {
	      if (alpha_ != 1.) {
		for (Reports_t::const_iterator it = reports.begin ();
		     it != reports.end (); ++it) {
		  CollisionConstraintsResult ccr (robot_, path0, path1,
						  *it, J_.rows (),
						  robotNbNonLockedDofs_);
		  // Compute J_
		  addCollisionConstraint (ccr, path0);
		  collisionConstraints.push_back (ccr);
		  hppDout (info, "Number of collision constraints: "<< collisionConstraints.size ());
		}
	      }
	      // When adding a new constraint, try first minimum under this
	      // constraint. If this latter minimum is in collision,
	      // re-initialize alpha_ to alphaInit_.
	      if (alpha_ == 1.)
		alpha_ = alphaInit_;
	      else
		if (validConstraints) // else do not modify alpha
		  alpha_ = 1.;
	      //alpha_ = alphaInit_;
	      noCollision = false;
	    } else { // path valid
	      x0 = x1;
	      path0 = path1;
	      hppDout (info, "latest valid path: " << pathId - 1);
	      cost_->jacobian (grad, x0);
	      compressVector (grad.transpose (), rgrad_.transpose ());
	      updateReference (collisionConstraints, path0);
	      //assert (solveConstraints (x0, path0, collisionConstraints));
	      if (!solveConstraints (x0, path0, collisionConstraints)) {
		throw std::runtime_error ("violated constraints");
		hppDout (error, "violated constraints");
	      }
	      noCollision = true;
	      alpha_ = .5*(1 + alpha_);
	    }
	  } while (!(noCollision && minimumReached) && (!interrupt_));
	} // while (!minimumReached)
	(*cost_) (cost, x0);
	hppDout (info, "final cost: " << cost);
	return path0;
      }

      void GradientBased::integrate (vectorIn_t x0, vectorIn_t step,
				     vectorOut_t x1) const
      {
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
	  size_type sizeConstraint = configProjector->rightHandSide ().size ();
	  rows = sizeConstraint * nbWaypoints_;
	  rhs_.resize (rows);
	  // Fill right hand side
	  for (size_type i=0; i<nbWaypoints_; ++i) {
	    rhs_.segment (i*sizeConstraint, sizeConstraint) =
	      configProjector->rightHandSide ();
	  }
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
	  size_type sizeConstraint = configProjector->rightHandSide ().size ();
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
      (vectorOut_t& x, PathVectorPtr_t& path,
       CollisionConstraintsResults_t& collisionConstraints)
      {
	bool satisfied = true;
	if (problem ().constraints ()) {
	  for (size_type i=0; i<nbWaypoints_; ++i) {
	    if (!problem ().constraints ()->isSatisfied
		(x.segment (i*robotNumberDofs_, robotNumberDofs_))) {
	      satisfied = false;
	    }
	  }
	}
	path = PathVector::create (configSize_, robotNumberDofs_);
	vectorToPath (x, path);
	// Check that collision constraints are within bounds, if not
	// linearize around new value.
	for (CollisionConstraintsResults_t::iterator it =
	       collisionConstraints.begin (); it != collisionConstraints.end ();
	     ++it) {
	  if (!it->satisfied (path)) satisfied = false;
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
	  hppDout (info, "rhs_ - value = " << (rhs_ - value_).transpose ());
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
	hppDout (info, "Jrows = " << Jrows);
	J_.conservativeResize (Jrows + fSize_, J_.cols ());
	J_.middleRows (ccr.rowInJacobian (), fSize_).setZero ();
	value_.conservativeResize (Jrows + fSize_);
	value_.segment (ccr.rowInJacobian (), fSize_).setZero ();
	rhs_.conservativeResize (Jrows + fSize_);
	rhs_.segment (ccr.rowInJacobian (), fSize_).setZero ();
	ccr.linearize (path, J_, value_);
      }

      void GradientBased::updateReference
      (CollisionConstraintsResults_t& collisionConstraints,
       const PathVectorPtr_t& path)
      {
	for (CollisionConstraintsResults_t::iterator it =
	       collisionConstraints.begin (); it != collisionConstraints.end ();
	     ++it) {
	  it->updateReference (path, rhs_);
	}

      }

    } // namespace pathOptimization
  }  // namespace core
} // namespace hpp
