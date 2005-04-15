#include "EquationOfStateFactory.h"
#include "IdealGas.h"
#include <Uintah/Interface/ProblemSpec.h>
#include <SCICore/Malloc/Allocator.h>
#include <fstream>
#include <iostream>
#include <string>
using std::cerr;
using std::ifstream;
using std::ofstream;

using namespace Uintah::ICESpace;

EquationOfState* EquationOfStateFactory::create(ProblemSpecP& ps)
{
   for (ProblemSpecP child = ps->findBlock(); child != 0;
	child = child->findNextBlock()) {
      std::string mat_type = child->getNodeName();
      if (mat_type == "ideal_gas") {
		 return(scinew IdealGas(child));
      } else {
	 cerr << "Unknown Material Type R (" << mat_type << ")" << std::endl;;
	 //      exit(1);
      }
   }
   return 0;
}
//$Log$
//Revision 1.1.2.1  2000/10/19 05:17:41  sparker
//Merge changes from main branch into csafe_risky1
//
//Revision 1.1  2000/10/06 04:02:16  jas
//Move into a separate EOS directory.
//
//Revision 1.3  2000/10/04 23:42:29  jas
//Add IdealGas to the EOS factory.
//
//Revision 1.2  2000/10/04 20:17:52  jas
//Change namespace ICE to ICESpace.
//
//Revision 1.1  2000/10/04 19:26:14  guilkey
//Initial commit of some classes to help mainline ICE.
//




