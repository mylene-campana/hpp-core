//
// Copyright (c) 2014 CNRS
// Authors: Florent Lamiraux
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

#define BOOST_TEST_MODULE gradient_based

#include <cmath>
#include <boost/test/included/unit_test.hpp>

#include <hpp/fcl/collision_object.h>
#include <hpp/fcl/math/transform.h>
#include <hpp/fcl/shape/geometric_shapes.h>

#include <hpp/model/joint.hh>
#include <hpp/model/collision-object.hh>
#include <hpp/model/device.hh>
#include <hpp/model/object-factory.hh>

#include <hpp/core/steering-method-straight.hh>
#include <hpp/core/path-optimization/gradient-based.hh>
#include <hpp/core/path-vector.hh>
#include <hpp/core/problem.hh>

using hpp::model::BodyPtr_t;
using hpp::model::Body;
using hpp::model::CollisionObject;
using hpp::model::CollisionObjectPtr_t;
using hpp::model::Configuration_t;
using hpp::model::Device;
using hpp::model::DevicePtr_t;
using hpp::model::Transform3f;
using hpp::model::JointPtr_t;
using hpp::model::JointTranslation;
using hpp::model::ObjectFactory;
using hpp::model::value_type;
using fcl::Quaternion3f;
using fcl::Box;
using hpp::core::ConfigurationPtr_t;
using hpp::core::PathVector;
using hpp::core::PathVectorPtr_t;
using hpp::core::Problem;
using hpp::core::SteeringMethodPtr_t;
using hpp::core::SteeringMethodStraight;
using hpp::core::PathOptimizerPtr_t;
using hpp::core::pathOptimization::GradientBased;

BOOST_AUTO_TEST_SUITE( test_hpp_core )

DevicePtr_t createRobot ()
{
  DevicePtr_t robot = Device::create ("planar-robot");
  Transform3f position; position.setIdentity ();
  ObjectFactory factory;

  // planar root joint
  JointPtr_t root = factory.createJointTranslation2 (position);
  robot->rootJoint (root);
  // Rotation around z
  position.setQuatRotation (Quaternion3f (sqrt (2)/2, 0, -sqrt (2)/2, 0));
  JointPtr_t joint = factory.createUnBoundedJointRotation (position);
  root->addChildJoint (joint);
  position.setIdentity ();
  boost::shared_ptr <Box> box (new Box (1,2,1));
  CollisionObjectPtr_t object = CollisionObject::create (box, position, "box");
  BodyPtr_t body = new Body ();
  joint->setLinkedBody (body);
  body->addInnerObject (object, true, true);

  return robot;
}

BOOST_AUTO_TEST_CASE (BFGS)
{
  DevicePtr_t robot = createRobot ();
  Configuration_t q0 (robot->configSize ());
  Configuration_t q1 (robot->configSize ());
  Configuration_t q2 (robot->configSize ());
  Configuration_t q3 (robot->configSize ());
  Configuration_t q4 (robot->configSize ());
  value_type s = sqrt (2)/2;
  q0 (0) = -1; q0 (1) = 0; q0 (2) = 1; q0 (3) = 0;
  q1 (0) = -s; q1 (1) = s; q1 (2) = s; q1 (3) = s;
  q2 (0) = 0; q2 (1) = 1; q2 (2) = 0; q2 (3) = 1;
  q3 (0) = s; q3 (1) = s; q3 (2) = s; q3 (3) = s;
  q4 (0) = 1; q4 (1) = 0; q4 (2) = 1; q4 (3) = 0;

  Problem problem (robot);
  SteeringMethodPtr_t sm = problem.steeringMethod ();
  PathVectorPtr_t path = PathVector::create (robot->configSize (),
					     robot->numberDof ());
  path->appendPath ((*sm) (q0, q1));
  path->appendPath ((*sm) (q1, q2));
  path->appendPath ((*sm) (q2, q3));
  path->appendPath ((*sm) (q3, q4));

  PathOptimizerPtr_t pathOptimizer (GradientBased::create (problem));
  pathOptimizer->optimize (path);
}
BOOST_AUTO_TEST_SUITE_END()
