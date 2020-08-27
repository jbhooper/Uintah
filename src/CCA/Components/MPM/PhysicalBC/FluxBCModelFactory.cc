/*
 * The MIT License
 *
 * Copyright (c) 1997-2020 The University of Utah
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
#include <CCA/Components/MPM/PhysicalBC/FluxBCModelFactory.h>
#include <CCA/Components/MPM/PhysicalBC/FluxBCModel.h>
#include <CCA/Components/MPM/PhysicalBC/AutoCycleFluxBC.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Malloc/Allocator.h>
#include <string>

using namespace std;
using namespace Uintah;

FluxBCModel* FluxBCModelFactory::create(MaterialManagerP& ss, MPMFlags* flags)
{
  if (flags->d_doAutoCycleBC) {
	  std::cerr << "Instantiating autocycle flux BC." << std::endl;
	  return(scinew AutoCycleFluxBC(ss, flags));
  } else {
	  std::cerr << "Instantiating regular flux BC." << std::endl;
    return(scinew FluxBCModel(ss, flags));
  }

  return 0;
}
