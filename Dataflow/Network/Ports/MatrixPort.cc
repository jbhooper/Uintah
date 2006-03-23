/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2004 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/



/*
 *  MatrixPort.cc
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   July 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#include <Dataflow/Network/Ports/MatrixPort.h>
#include <Core/Malloc/Allocator.h>

#undef SCISHARE
#ifdef _WIN32
#define SCISHARE __declspec(dllexport)
#else
#define SCISHARE
#endif

namespace SCIRun {

extern "C" {
SCISHARE IPort* make_MatrixIPort(Module* module, const string& name) {
  return scinew SimpleIPort<MatrixHandle>(module,name);
}
SCISHARE OPort* make_MatrixOPort(Module* module, const string& name) {
  return scinew SimpleOPort<MatrixHandle>(module,name);
}
}

template<> string SimpleIPort<MatrixHandle>::port_type_("Matrix");
template<> string SimpleIPort<MatrixHandle>::port_color_("dodgerblue");

template class SimpleIPort<MatrixHandle>;
template class SimpleOPort<MatrixHandle>;

} // End namespace SCIRun

