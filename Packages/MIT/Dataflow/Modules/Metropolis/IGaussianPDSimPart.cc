/*
 *  IGaussianPDSimPart.cc:
 *
 *  Written by:
 *   Yarden Livnat
 *   Sep-2001
 *
 */

#include <Packages/MIT/Dataflow/Modules/Metropolis/Sampler.h>
#include <Packages/MIT/Dataflow/Modules/Metropolis/IGaussianPDSimPart.h>

namespace MIT {

using namespace SCIRun;

static double sqr( double x ) { return x*x; }

IGaussianPDSimPart::IGaussianPDSimPart( Sampler *sampler, 
							PartInterface *parent,
							const string &name )
  :PDSimPart( sampler, parent, name )
{
  UNUR_DISTR *distr = unur_distr_normal(0,0);
  UNUR_PAR *par = unur_arou_new(distr);
  gen_ = unur_init(par);
  if ( ! gen_ ) 
    cerr << "Error, cannot create generator object\n";
  unur_distr_free(distr);
}

IGaussianPDSimPart::~IGaussianPDSimPart()
{
}
  
void
IGaussianPDSimPart::compute( vector<double> &theta, 
				     vector<double> &star )
{
  int n = theta.size();

  vector<double> r(n);

  for (int i=0; i<n; i++)
    r[i] = unur_sample_cont(gen_);

  Array2<double> &lkappa = sampler_->get_lkappa();

  for (int i=0; i<n; i++) {
    star[i] = theta[i];
    for (int j=0; j<n; j++ )
      star[i] += lkappa(i,j) * r[j];
  }
}

double
IGaussianPDSimPart::lpr( vector<double> &)
{
  return 0;
}

} // End namespace MIT


