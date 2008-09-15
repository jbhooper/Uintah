//  AngioMaterial.cc

#include <Packages/Uintah/CCA/Components/Angio/AngioParticleCreator.h>
#include <Packages/Uintah/CCA/Components/Angio/AngioMaterial.h>
#include <Packages/Uintah/CCA/Components/Angio/AngioFlags.h>
#include <Packages/Uintah/Core/GeometryPiece/GeometryObject.h>
#include <Core/Geometry/IntVector.h>
#include <Packages/Uintah/Core/Grid/Box.h>
#include <Packages/Uintah/Core/Grid/Patch.h>
#include <Packages/Uintah/Core/Grid/Variables/VarLabel.h>
#include <Packages/Uintah/Core/Grid/Variables/PerPatch.h>
#include <Packages/Uintah/Core/Labels/AngioLabel.h>
#include <Packages/Uintah/Core/GeometryPiece/GeometryPieceFactory.h>
#include <Packages/Uintah/Core/GeometryPiece/UnionGeometryPiece.h>
#include <Packages/Uintah/Core/GeometryPiece/NullGeometryPiece.h>
#include <Packages/Uintah/Core/Exceptions/ParameterNotFound.h>
#include <Packages/Uintah/CCA/Ports/DataWarehouse.h>
#include <Packages/Uintah/Core/ProblemSpec/ProblemSpecP.h>
#include <sgi_stl_warnings_off.h>
#include <iostream>
#include <string>
#include <list>
#include <sgi_stl_warnings_on.h>

#define d_TINY_RHO 1.0e-12 // also defined  ICE.cc and ICEMaterial.cc 

using namespace std;
using namespace Uintah;
using namespace SCIRun;

// Standard Constructor
AngioMaterial::AngioMaterial(ProblemSpecP& ps, SimulationStateP& ss,AngioFlags* flags) : Material(ps)
{
  d_lb = scinew AngioLabel();
  // The standard set of initializations needed
  standardInitialization(ps,flags);
  
  // Create a ParticleCreator object
  d_particle_creator = scinew AngioParticleCreator(this,flags);
}

void
AngioMaterial::standardInitialization(ProblemSpecP& ps, AngioFlags* flags)
{
  ps->require("density",d_density);
  ps->require("fragment_file", d_init_frag_file);
}

// Default constructor
AngioMaterial::AngioMaterial() : d_particle_creator(0)
{
  d_lb = scinew AngioLabel();
}

AngioMaterial::~AngioMaterial()
{
  delete d_lb;
  delete d_particle_creator;

  for (int i = 0; i<(int)d_geom_objs.size(); i++) {
    delete d_geom_objs[i];
  }
}

void AngioMaterial::registerParticleState(SimulationState* sharedState)
{
  sharedState->d_particleState.push_back(d_particle_creator->returnParticleState());
  sharedState->d_particleState_preReloc.push_back(d_particle_creator->returnParticleStatePreReloc());
}

ProblemSpecP AngioMaterial::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP angio_ps = Material::outputProblemSpec(ps);
  angio_ps->appendElement("density",d_density);

  for (vector<GeometryObject*>::const_iterator it = d_geom_objs.begin();
       it != d_geom_objs.end(); it++) {
    (*it)->outputProblemSpec(angio_ps);
  }

  return angio_ps;
}

void AngioMaterial::copyWithoutGeom(ProblemSpecP& ps,const AngioMaterial* mat, 
                                    AngioFlags* flags)
{
  d_density = mat->d_density;

  // Check to see which ParticleCreator object we need
  d_particle_creator = scinew AngioParticleCreator(this,flags);
}

particleIndex AngioMaterial::countParticles(const Patch* patch)
{
  return d_particle_creator->countParticles(patch, d_init_frag_file);
}

void AngioMaterial::createParticles(particleIndex numParticles,
                                    CCVariable<short int>& cellNAPID,
                                    const Patch* patch,
                                    DataWarehouse* new_dw)
{
  d_particle_creator->createParticles(this,d_init_frag_file,
                                      numParticles,cellNAPID,
                                      patch,new_dw);
}

AngioParticleCreator* AngioMaterial::getParticleCreator()
{
  return  d_particle_creator;
}

double AngioMaterial::getInitialDensity() const
{
  return d_density;
}

int AngioMaterial::nullGeomObject() const
{
  for (int obj = 0; obj <(int)d_geom_objs.size(); obj++) {
    GeometryPieceP piece = d_geom_objs[obj]->getPiece();
    NullGeometryPiece* null_piece = dynamic_cast<NullGeometryPiece*>(piece.get_rep());
    if (null_piece)
      return obj;
  }
  return -1;
}
