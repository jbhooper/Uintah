#ifndef __STABILITY_CHECK_H__
#define __STABILITY_CHECK_H__

#include <Packages/Uintah/Core/Math/Matrix3.h>
#include <Packages/Uintah/Core/Math/TangentModulusTensor.h>

namespace Uintah {

  /*! \class StabilityCheck
   *  \brief  A generic wrapper for various methods of checking the stability of
   *  \brief  the motion, e.g. loss of hyperbolicity/ellipticity, Drucker 
   *  \brief  stability criterion, Hill condition etc.
   *  \brief  Provides an abstract base class for various methods of checking 
   *  \brief  the stability of motion/bifurcation points
   *  \author  Biswajit Banerjee, 
   *  \author  C-SAFE and Department of Mechanical Engineering,
   *  \author  University of Utah.
   *  \author  Copyright (C) 2003 Container Dynamics Group
  */
  class StabilityCheck {

  public:
	 
    //! Construct an object that can be used to check stability
    StabilityCheck();

    //! Destructor of stability check
    virtual ~StabilityCheck();
	 
    /*! Check the stability and return the direction of instability
        if any */
    virtual bool checkStability(const Matrix3& cauchyStress,
                               const Matrix3& strainMeasure,
                               TangentModulusTensor& tangentModulus,
                               Vector& direction) = 0;
  };
} // End namespace Uintah
      
#endif  // __STABILITY_CHECK_H__

