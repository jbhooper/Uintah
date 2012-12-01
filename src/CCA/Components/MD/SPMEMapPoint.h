/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef UINTAH_SPME_MAP_POINT_H
#define UINTAH_SPME_MAP_POINT_H

#include <Core/Geometry/Vector.h>
#include <Core/Geometry/BBox.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Grid/Variables/Array3.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>

#include <vector>
#include <list>

namespace Uintah {

using SCIRun::Vector;
using SCIRun::IntVector;

template<class T> class SPMEMapPoint {

  public:
    SPMEMapPoint();

    ~SPMEMapPoint();

    inline const double GetParticleCharge() const
    {
      // TODO Figure out what charge is
//      return Particle->Charge();
      return this->globalParticleSubset;
    }

    inline const double ChargeWeight() const
    {
      return this->d_chargeWeight;
    }

    inline void setParticleCharge(double charge)
    {

    }

  private:
    // Some manner of mapping the provided particle pointer/reference to a storable way to index the particle
    ParticleSubset& globalParticleSubset;
    double d_chargeWeight;
    Vector d_forceGRadient;

};

}  // End namespace Uintah

#endif
