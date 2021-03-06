//crabe8pince
// Copyright (c) 2015 CNRS
// Authors: Mylene Campana
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

# include <sstream>
# include <hpp/util/debug.hh>
# include <hpp/model/collision-object.hh>
# include <hpp/model/configuration.hh>
# include <hpp/model/device.hh>
# include <hpp/model/joint.hh>
# include <hpp/model/joint-configuration.hh>
# include <hpp/core/configuration-shooter.hh>
# include <hpp/core/configuration-projection-shooter.hh>
# include <hpp/core/config-validations.hh>
# include <hpp/core/distance-between-objects.hh>
#include <hpp/core/parabola/parabola-library.hh>
# include <fcl/distance.h>

namespace hpp {
  namespace core {
    using model::displayConfig;
    
    ConfigurationPtr_t ConfigurationProjectionShooter::shoot () const
      {
	ConfigValidationsPtr_t configValidations (problem_.configValidations());
	ValidationReportPtr_t validationReport;
	ConfigurationPtr_t config (new Configuration_t (robot_->configSize ()));
	// Sample a collision-free configuration
	do {
	  *config = uniformlySample ();
	}
	while (!configValidations->validate (*config, validationReport));

	// Project on nearest obstacle and shift away
	*config = project (*config);
	return config;
      }

      Configuration_t ConfigurationProjectionShooter::uniformlySample () const {
	Configuration_t q (robot_->configSize ());
	JointVector_t jv = robot_->getJointVector ();
	for (JointVector_t::const_iterator itJoint = jv.begin ();
	     itJoint != jv.end (); itJoint++) {
	  std::size_t rank = (*itJoint)->rankInConfiguration ();
	  (*itJoint)->configuration ()->uniformlySample (rank, q);
	}
	// Shoot extra configuration variables
	size_type extraDim = robot_->extraConfigSpace ().dimension ();
	size_type offset = robot_->configSize () - extraDim;
	for (size_type i=0; i<extraDim; ++i) {
	  value_type lower = robot_->extraConfigSpace ().lower (i);
	  value_type upper = robot_->extraConfigSpace ().upper (i);
	  value_type range = upper - lower;
	  if ((range < 0) ||
	      (range == std::numeric_limits<double>::infinity())) {
	    std::ostringstream oss
	      ("Cannot uniformy sample extra config variable ");
	    oss << i << ". min = " <<lower<< ", max = " << upper << std::endl;
	    throw std::runtime_error (oss.str ());
	  }
	  q [offset + i] = lower + (upper - lower) * rand ()/RAND_MAX;
	}
	return q;
      }

    Configuration_t ConfigurationProjectionShooter::project
    (const Configuration_t q) const {
      Configuration_t qout = q;
      const size_type ecsDim = robot_->extraConfigSpace ().dimension ();
      const size_type index = robot_->configSize() - ecsDim; // ecs index
      fcl::Vec3f pi, pj, dir, n; // fcl nearest points of collision pairs
      value_type minDistance = std::numeric_limits <value_type>::infinity();
      value_type distance = minDistance;
      CollisionObjectPtr_t nearestRob, nearestObst;
      DistanceBetweenObjectsPtr_t distanceBetweenObjects
	(problem_.distanceBetweenObjects ());
      ConfigValidationsPtr_t configValidations (problem_.configValidations());
      ValidationReportPtr_t validationReport;

      // Step 1: get nearest obstacle and surface-normale at config q
      robot_->currentConfiguration (q);
      robot_->computeForwardKinematics ();
      distanceBetweenObjects->computeDistances (); // only outers !
      const model::DistanceResults_t& dr =
	distanceBetweenObjects->distanceResults ();

      for (model::DistanceResults_t::const_iterator itDistance = 
	     dr.begin (); itDistance != dr.end (); itDistance++) {
	distance = itDistance->distance ();
	if (distance < minDistance){
	  minDistance = distance;
	  nearestRob = itDistance->innerObject;
	  nearestObst = itDistance->outerObject;
	  pi = itDistance->closestPointInner (); // point Body
	  pj = itDistance->closestPointOuter (); // point Obst
	  dir = pi - pj; // obstacle normale direction
	}
      }
      const value_type dir_norm = sqrt (dir [0]*dir [0] + dir [1]*dir [1]
					+ dir [2]*dir [2]);
      n = dir/dir_norm;

      // Step 2: set orientation with n
      qout (index) = n [0];
      qout (index + 1) = n [1];
      if (ecsDim != 2)
	qout (index + 2) = n [2];
      qout = setOrientation (robot_, qout);
      if (!configValidations->validate (qout, validationReport))
	return qout; // thrown by planner

      // Step 3: re-compute new distance to nearestObst
      fcl::DistanceRequest distanceRequest (true, 0, 0, fcl::GST_INDEP);
      model::DistanceResult dr1;
      robot_->currentConfiguration (qout);
      robot_->computeForwardKinematics (); // may not be necessary
      fcl::distance (nearestRob->fcl ().get (), nearestObst->fcl ().get (),
		     distanceRequest, dr1.fcl);
      const value_type dist = dr1.distance ();

      // Step 4: project and shift
      if (ecsDim == 2) { /* 2D*/
	const value_type gamma = atan2(q (3), q (2)) - M_PI/2;
	qout (0) -= dir [0] + shiftDistance_*sin(gamma); // x
	qout (1) -= dir [1] - shiftDistance_*cos(gamma); // y
      }
      else { /* 3D */
	qout (0) += - dist * n [0] + shiftDistance_ * n [0]; // x
	qout (1) += - dist * n [1] + shiftDistance_ * n [1]; // y
	qout (2) += - dist * n [2] + shiftDistance_ * n [2]; // z
      }
      return qout;
    }
    /// \}
  } //   namespace core
} // namespace hpp

