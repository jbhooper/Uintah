
#ifndef Uintah_Component_Arches_PhysicalConstants_h
#define Uintah_Component_Arches_PhysicalConstants_h

/*
 * Placeholder - nothing here yet
 */

#include <Uintah/Parallel/UintahParallelComponent.h>
#include <Uintah/Interface/CFDInterface.h>
#include <Uintah/Grid/Region.h>
#include <Uintah/Parallel/ProcessorContext.h>
#include <SCICore/Geometry/Vector.h>

namespace Uintah {
namespace Components {
  using Uintah::Interface::CFDInterface;
  using Uintah::Parallel::UintahParallelComponent;
  using Uintah::Interface::ProblemSpecP;
  using Uintah::Interface::DataWarehouseP;
  using SCICore::Geometry::Vector;

class PhysicalConstants {
public:
    PhysicalConstants();
    ~PhysicalConstants();

    void problemSetup(const ProblemSpecP& params);
    const Vector& getGravity(){
      return d_gravity;
    }
  double getMolecularViscosity() {
      return d_viscosity;
    }
  double getabsPressure() {
      return d_absPressure;
    }
    
private:
    PhysicalConstants(const PhysicalConstants&);
    PhysicalConstants& operator=(const PhysicalConstants&);
    Vector d_gravity;
    double d_viscosity;
    double d_absPressure;
    
};

} // end namespace Components
} // end namespace Uintah


#endif

