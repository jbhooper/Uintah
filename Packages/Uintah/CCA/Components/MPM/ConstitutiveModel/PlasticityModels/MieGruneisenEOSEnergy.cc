
#include "MieGruneisenEOSEnergy.h"
#include <cmath>

using namespace Uintah;
using namespace SCIRun;

MieGruneisenEOSEnergy::MieGruneisenEOSEnergy(ProblemSpecP& ps)
{
  ps->require("C_0",d_const.C_0);
  ps->require("Gamma_0",d_const.Gamma_0);
  ps->require("S_alpha",d_const.S_1);
  ps->getWithDefault("S_2",d_const.S_2,0.0);
  ps->getWithDefault("S_3",d_const.S_3,0.0);
} 
	 
MieGruneisenEOSEnergy::MieGruneisenEOSEnergy(const MieGruneisenEOSEnergy* cm)
{
  d_const.C_0 = cm->d_const.C_0;
  d_const.Gamma_0 = cm->d_const.Gamma_0;
  d_const.S_1 = cm->d_const.S_1;
  d_const.S_2 = cm->d_const.S_2;
  d_const.S_3 = cm->d_const.S_3;
} 
	 
MieGruneisenEOSEnergy::~MieGruneisenEOSEnergy()
{
}
	 
void MieGruneisenEOSEnergy::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP eos_ps = ps->appendChild("equation_of_state");
  eos_ps->setAttribute("type","mie_gruneisen");

  eos_ps->appendElement("C_0",d_const.C_0);
  eos_ps->appendElement("Gamma_0",d_const.Gamma_0);
  eos_ps->appendElement("S_alpha",d_const.S_1);
  eos_ps->appendElement("S_2",d_const.S_2);
  eos_ps->appendElement("S_3",d_const.S_3);
}

//////////
// Calculate the pressure using the Mie-Gruneisen equation of state
double 
MieGruneisenEOSEnergy::computePressure(const MPMMaterial* matl,
                                 const PlasticityState* state,
                                 const Matrix3& ,
                                 const Matrix3& ,
                                 const double& )
{
  // Get the current density
  double rho = state->density;

  // Get original density
  double rho_0 = matl->getInitialDensity();
   
  // Calc. eta
  double eta = 1. - rho_0/rho;

  // Retrieve specific internal energy e
  double e = state->energy;

  // Calculate the pressure
  double denom = 
               (1.-d_const.S_1*eta-d_const.S_2*eta*eta-d_const.S_3*eta*eta*eta)
              *(1.-d_const.S_1*eta-d_const.S_2*eta*eta-d_const.S_3*eta*eta*eta);
  double p;
  p = rho_0*d_const.Gamma_0*e 
    + rho_0*(d_const.C_0*d_const.C_0)*eta*(1. - .5*d_const.Gamma_0*eta)/denom;

  return -p;
}


double 
MieGruneisenEOSEnergy::computeIsentropicTemperatureRate(const double T,
                                                        const double rho_0,
                                                        const double rho_cur,
                                                        const double Dtrace)
{
  double dTdt = -T*d_const.Gamma_0*rho_0*Dtrace/rho_cur;

  return dTdt;
}

double 
MieGruneisenEOSEnergy::eval_dp_dJ(const MPMMaterial* matl,
                            const double& detF, 
                            const PlasticityState* state)
{
  double rho_0 = matl->getInitialDensity();
  double C_0 = d_const.C_0;
  double S_1 = d_const.S_1;
  double S_2 = d_const.S_2;
  double S_3 = d_const.S_3;
  double Gamma_0 = d_const.Gamma_0;

  double J = detF;
  double numer = rho_0*C_0*C_0*(1.0 + (S_1 - Gamma_0)*(1.0-J));
  double denom = (1.0 - S_1*(1.0-J));
  double denom3 = (denom*denom*denom);
  if (denom3 == 0.0) {
    cout << "rh0_0 = " << rho_0 << " J = " << J 
           << " numer = " << numer << endl;
    denom3 = 1.0e-5;
  }

  return (numer/denom);
}
