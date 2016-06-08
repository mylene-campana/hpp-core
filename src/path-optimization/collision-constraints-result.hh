//
// Copyright (c) 2015 CNRS
// Author: Mylene Campana
//
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

#ifndef HPP_CORE_PATH_OPTIMIZATION_COLLISION_CONSTRAINTS_RESULT_HH
# define HPP_CORE_PATH_OPTIMIZATION_COLLISION_CONSTRAINTS_RESULT_HH

# include <hpp/fcl/distance.h>
# include <hpp/model/fcl-to-eigen.hh>
# include <hpp/constraints/generic-transformation.hh>

namespace hpp {
  namespace core {
    using model::displayConfig;
    namespace pathOptimization {
      HPP_PREDEF_CLASS (CollisionConstraintsResult);
      HPP_PREDEF_CLASS (CollisionConstraint);
      typedef boost::shared_ptr <CollisionConstraint> CollisionConstraintPtr_t;
      typedef model::Transform3f Transform3f;
      namespace eigen {
	typedef Eigen::Matrix <value_type, 3, 1> vector3_t;
      } // namespace eigen


      class CollisionConstraint : public DifferentiableFunction
      {
      public:
       virtual ~CollisionConstraint () {}
       static CollisionConstraintPtr_t create
       (const DevicePtr_t& robot, const Configuration_t& qFree,
	const Configuration_t& qColl, const CollisionObjectPtr_t& object1,
	const CollisionObjectPtr_t& object2)
       {
         CollisionConstraint* ptr = new CollisionConstraint
           (robot, qFree, qColl, object1, object2);
         CollisionConstraintPtr_t shPtr (ptr);
         return shPtr;
       }
      protected:
       CollisionConstraint (const DevicePtr_t& robot,
                            const Configuration_t& qFree,
                            const Configuration_t& qColl,
			    const CollisionObjectPtr_t& object1,
			    const CollisionObjectPtr_t& object2)
         : DifferentiableFunction (robot->configSize (), robot->numberDof (),
                                   1, ""), robot_ (robot), qFree_ (qFree),
           J_ (), difference_ ()
	   
       {
	 difference_.resize (robot->numberDof ());
	 // Compute contact point in configuration qColl
	 robot_->currentConfiguration (qColl);
	 robot_->computeForwardKinematics ();
	 fcl::CollisionResult result;
	 fcl::CollisionRequest collisionRequest (1, true, false, 1, false, true,
						 fcl::GST_INDEP);
	 fcl::collide (object1->fcl ().get (), object2->fcl ().get (),
		       collisionRequest, result);
	 assert (result.numContacts () == 1);
	 const vector3_t& contactPoint (result.getContact (0).pos);
	 hppDout (info, "contact point = " << contactPoint);

	 JointPtr_t joint1 = object1->joint ();
	 Transform3f M1 (joint1->currentTransformation ());
	 if (object2->joint ()) { // object2 = body part
	   JointPtr_t joint2 = object2->joint ();
	   Transform3f M2 (joint2->currentTransformation ());
	   // Position of contact point in each object local frame
	   vector3_t x1_J1 = inverse (M1).transform (contactPoint);
	   vector3_t x2_J2 = inverse (M2).transform (contactPoint);
	   // Compute contact points in configuration qFree
	   robot_->currentConfiguration (qFree);
	   robot_->computeForwardKinematics ();
	   M2 = joint2->currentTransformation ();
	   M1 = joint1->currentTransformation ();
	   // Position of x2 in local frame of joint1
	   vector3_t x2_J1 (inverse (M1).transform (M2.transform (x2_J2)));
	   hppDout (info, "x1 in J1 = " << x1_J1);
	   hppDout (info, "x2 in J1 = " << x2_J1);
	   eigen::vector3_t u; model::toEigen (x2_J1 - x1_J1, u);
	   DifferentiableFunctionPtr_t f = constraints::RelativePosition::create
	     ("", robot_, joint1, joint2, x1_J1, x2_J2);
	   matrix_t Jpos (f->outputSize (), f->inputDerivativeSize ());
	   f->jacobian (Jpos, qFree);
	   J_ = u.transpose () * Jpos;
	   assert (J_.rows () == 1);
	 } else{ // object2 = fixed obstacle and has no joint
	   vector3_t x1_J1 (inverse (M1).transform (contactPoint));
	   vector3_t x2_J2 (contactPoint);
	   // Compute contact points in configuration qFree
	   robot_->currentConfiguration (qFree);
	   robot_->computeForwardKinematics ();
	   Transform3f M1 (joint1->currentTransformation ());
	   // position of x1 in global frame
	   vector3_t x1_J2 (M1.transform (x1_J1));
	   hppDout (info, "x1 in J2 = " << x1_J2);
	   eigen::vector3_t u; model::toEigen (x1_J2 - x2_J2, u);

	   DifferentiableFunctionPtr_t f = constraints::Position::create
	     ("", robot_, joint1, x1_J1, x2_J2);
	   matrix_t Jpos (f->outputSize (), f->inputDerivativeSize ());
	   f->jacobian (Jpos, qFree);
	   J_ = u.transpose () * Jpos;
	   assert (J_.rows () == 1);
	 }
       }

       virtual void impl_compute (vectorOut_t result, vectorIn_t argument)
         const
       {
         model::difference (robot_, argument, qFree_, difference_);
         result = J_ * difference_;
       }
       virtual void impl_jacobian (matrixOut_t jacobian, vectorIn_t) const
       {
         jacobian = J_;
       }
      private:
	DevicePtr_t robot_;
	Configuration_t qFree_;
	matrix_t J_;
	mutable vector_t difference_;
      }; // class CollisionConstraint

      /// Storing data about one collision, to facilitate building
      /// collision-constraints in Gradient-Based path-optimizer

      struct HPP_CORE_DLLAPI CollisionConstraintsResult
      {
	static size_type fSize_;
	CollisionConstraintsResult
	  (const DevicePtr_t& robot, const PathVectorPtr_t& previousPath,
	   const PathVectorPtr_t& currentPath,
	   const std::pair <CollisionPathValidationReportPtr_t,
	   std::size_t>& report,
	   size_type rowInJacobian, size_type nbNonLockedDofs) :
	  robot_ (robot), robotNumberDofs_ (robot_->numberDof ()),
	  nbNonLockedDofs_ (nbNonLockedDofs), rowInJacobian_ (rowInJacobian),
	  nbWaypoints_ ((size_type) previousPath->numberPaths () - 1)
	  {
	    Jcompressed_.resize (fSize_, nbNonLockedDofs_);
	    Jcompressed_.setZero ();

	    value_type t_local = report.first->parameter;
	    localPathId_ = report.second;
	    hppDout (info, "localPathId_ = " << localPathId_);
	    PathPtr_t localPath = currentPath->pathAtRank(localPathId_);
	    bool success;
	    Configuration_t qColl = (*localPath) (t_local, success);
	    posAlongLocalPath_ = t_local/localPath->length ();
	    HPP_STATIC_CAST_REF_CHECK (CollisionValidationReport,
				       *(report.first->configurationReport));
	    object1_ =
	      HPP_STATIC_PTR_CAST (CollisionValidationReport,
				   report.first->configurationReport)->object1;
	    object2_ =
	      HPP_STATIC_PTR_CAST (CollisionValidationReport,
				   report.first->configurationReport)->object2;
	    hppDout (info, "posAlongLocalPath_ (in [0,1]) = "
		     << posAlongLocalPath_);
	    hppDout (info, "obj1 = " << object1_->name()
		     << " and obj2 = " << object2_->name());
	    hppDout (info, "qColl = " << displayConfig (qColl));

	    // Backtrack collision in previous path (x0) to create constraint
	    PathPtr_t prevLocalPath = previousPath->pathAtRank (localPathId_);
	    value_type t_local_new =
	      prevLocalPath->length () * posAlongLocalPath_;
	    qFree_ = (*prevLocalPath) (t_local_new, success);
	    hppDout (info, "qFree_ = "
		     << displayConfig (qFree_));

	    f_ = CollisionConstraint::create (robot_, qFree_, qColl,
					      object1_, object2_);
	    configProjector_ = ConfigProjector::create (robot_,
							"collision constraint",
							1e-6, 30);
	    configProjector_->add (NumericalConstraint::create (f_));
	  }

	value_type distance () const { return distance_;}

	/// Get id of row in Jacobian
	size_type rowInJacobian () const { return rowInJacobian_;}

	bool linearize (const PathVectorPtr_t& path, matrixOut_t jacobian,
			vectorOut_t value) const
	{
	  PathPtr_t localPath (path->pathAtRank (localPathId_));
	  value_type t_local = localPath->length () * posAlongLocalPath_;
	  bool success;
	  Configuration_t q = (*localPath) (t_local, success);
	  size_type rank = localPathId_ - 1;
	  value_type position = 0; // 0=not extremity, 1=begin, 2=end
	  if (rank + 1 == 0)
	    position = 1;
	  if (rank + 1 == (size_type) nbWaypoints_)
	    position = 2;

	  configProjector_->computeValueAndJacobian
	    (q, value.segment (rowInJacobian_, fSize_), Jcompressed_);
	  hppDout (info, "value = " << value [rowInJacobian_]);
	  // Complete 2 blocks in J_ except for begin and end segments
	  if (position != 1) {
	    jacobian.block (rowInJacobian_, rank*nbNonLockedDofs_,
			    fSize_, nbNonLockedDofs_) =
	      (1-posAlongLocalPath_)*Jcompressed_;
	  }
	  if (position != 2) {
	    jacobian.block (rowInJacobian_, (rank+1)*nbNonLockedDofs_,
			    fSize_, nbNonLockedDofs_) =
	      posAlongLocalPath_ * Jcompressed_;
	  }
	  return true;
	}

	void updateRightHandSide (const PathVectorPtr_t& path, vectorOut_t rhs)
	  const
	{
	  PathPtr_t localPath (path->pathAtRank (localPathId_));
	  value_type t_local = localPath->length () * posAlongLocalPath_;
	  bool success;
	  Configuration_t q = (*localPath) (t_local, success);
	  (*f_) (rhs.segment (rowInJacobian_, fSize_), q);
	}

        void add (const LockedJoints_t& lj) {
            for (LockedJoints_t::const_iterator _lj = lj.begin ();
                _lj != lj.end (); ++_lj)
              configProjector_->add (*_lj);
        }

      private:
	void computeConstraint ()
	{
	}

      private:
	DevicePtr_t robot_;
	size_type robotNumberDofs_;
	size_type nbNonLockedDofs_;
	CollisionObjectPtr_t object1_;
	CollisionObjectPtr_t object2_;
	value_type posAlongLocalPath_; // posAlongLocalPath_ \in [0,1]
	size_type localPathId_; // index in global path
	ConfigProjectorPtr_t configProjector_;
	DifferentiableFunctionPtr_t f_;
	Configuration_t qFree_;
	value_type distance_;
	value_type radius_;
	size_type rowInJacobian_;
	size_type nbWaypoints_;
	mutable matrix_t Jconstraint_;
	mutable matrix_t Jcompressed_;
      }; // struct CollisionConstraintsResult
      size_type CollisionConstraintsResult::fSize_ = 1;
      typedef std::vector <CollisionConstraintsResult>
      CollisionConstraintsResults_t;
    } // namespace pathOptimization
  } // namespace core
} // namespace hpp
#endif // HPP_CORE_PATH_OPTIMIZATION_COLLISION_CONSTRAINTS_RESULT_HH
