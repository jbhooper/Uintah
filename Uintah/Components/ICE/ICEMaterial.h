#ifndef __ICE_MATERIAL_H__
#define __ICE_MATERIAL_H__

#include <Uintah/Interface/DataWarehouseP.h>
#include <Uintah/Grid/Material.h>
#include <Uintah/Interface/ProblemSpecP.h>
#include <Uintah/Grid/ParticleVariable.h>
#include <Uintah/Grid/PerPatch.h>
#include <vector>
#include <Uintah/Components/ICE/ICELabel.h>
#include <Uintah/Components/ICE/EquationOfState.h>

namespace SCICore {
   namespace Geometry {
      class Point;
      class Vector;
   }
}

namespace Uintah {
   class Patch;
   class VarLabel;
   namespace ICESpace {
      using SCICore::Geometry::Point;
      using SCICore::Geometry::Vector;
      
/**************************************
     
CLASS
   ICEMaterial

   Short description...

GENERAL INFORMATION

   ICEMaterial.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)

   Copyright (C) 2000 SCI Group

KEYWORDS
   ICE

DESCRIPTION
   Long description...

WARNING

****************************************/

      class ICEMaterial : public Material {
      public:
	 ICEMaterial(ProblemSpecP&);
	 
	 ~ICEMaterial();
	 
	 //////////
	 // Return correct EOS model pointer for this material
	 EquationOfState * getEOS() const;

         //for HeatConductionModel
         double getThermalConductivity() const;
         double getSpecificHeat() const;
         double getHeatTransferCoefficient() const;

      private:

	 // Specific constitutive model associated with this material
	 EquationOfState *d_eos;

	 double d_density;
         double d_thermalConductivity;
         double d_specificHeat;
	 double d_heatTransferCoefficient;
         
	 ICELabel* lb;

	 // Prevent copying of this class
	 // copy constructor
	 ICEMaterial(const ICEMaterial &icem);
	 ICEMaterial& operator=(const ICEMaterial &icem);
      };

} // end namespace ICE
} // end namespace Uintah

#endif // __ICE_MATERIAL_H__

// $Log$
// Revision 1.2  2000/10/04 20:17:52  jas
// Change namespace ICE to ICESpace.
//
// Revision 1.1  2000/10/04 19:26:14  guilkey
// Initial commit of some classes to help mainline ICE.
//
