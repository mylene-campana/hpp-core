/*
 *  Copyright (c) 2010 LAAS-CNRS
 *
 *  Author: Florent Lamiraux
 */

#ifndef HPP_CORE_ROTATION_JOINT_HH
#define HPP_CORE_ROTATION_JOINT_HH

#include <KineoModel/kppRotationJointComponent.h>

#include "joint.hh"

namespace hpp {
  namespace core {
    namespace io {
      KIT_PREDEF_CLASS(RotationJoint);
      ///
      /// \brief Intermediate rotation joint.
      ///
      /// This class implements a rotation joint deriving from
      /// CkppRotationJointComponent only for kxml input output purposes.
      ///
      /// Joints of this class contain inertia data as CkppDoubleProperty
      /// attributes.
      ///
      /// Once read from a file the device containing this joint will be
      /// transformed into a ChppHumanoidRobot.

      class RotationJoint : public CkppRotationJointComponent,
			    public Joint
      {
      public:
	virtual bool isComponentClonable () const
	{
	  return false;
	}
	static RotationJointShPtr create(const std::string& inName)
	{
	  RotationJoint *ptr = new RotationJoint();
	  RotationJointShPtr shPtr = RotationJointShPtr(ptr);
	  RotationJointWkPtr wkPtr = RotationJointWkPtr(shPtr);

	  if (ptr->init(wkPtr, inName) != KD_OK) {
	    shPtr.reset();
	    return shPtr;
	  }
	  return shPtr;
	}
	void
	fillPropertyVector(std::vector<CkppPropertyShPtr>& inOutPropertyVector)
	  const
	{
	  CkppRotationJointComponent::fillPropertyVector(inOutPropertyVector);
	  inOutPropertyVector.push_back(mass);
	  inOutPropertyVector.push_back(comX);
	  inOutPropertyVector.push_back(comY);
	  inOutPropertyVector.push_back(comZ);
	  inOutPropertyVector.push_back(inertiaMatrixXX);
	  inOutPropertyVector.push_back(inertiaMatrixYY);
	  inOutPropertyVector.push_back(inertiaMatrixZZ);
	  inOutPropertyVector.push_back(inertiaMatrixXY);
	  inOutPropertyVector.push_back(inertiaMatrixXZ);
	  inOutPropertyVector.push_back(inertiaMatrixYZ);
	}

	virtual CkwsJointShPtr kwsJoint() const
	{
	  return CkppJointComponent::kwsJoint();
	}

      protected:
	RotationJoint() : CkppRotationJointComponent() {}
	ktStatus init (const RotationJointWkPtr &inWeakPtr,
		       const std::string &inName)
	{
	  ktStatus status = KD_OK;
	  status = CkppRotationJointComponent::init(inWeakPtr, inName);
	  if (status == KD_ERROR) return KD_ERROR;
	  return Joint::init(inWeakPtr);
	}
      };
    } // namespace io
  } // namespace core
} // namespace hpp
#endif //HPP_CORE_ROTATION_JOINT_HH
