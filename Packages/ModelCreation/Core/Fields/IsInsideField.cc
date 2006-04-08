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

#include <Packages/ModelCreation/Core/Fields/IsInsideField.h>

namespace ModelCreation {

using namespace SCIRun;

bool IsInsideFieldAlgo::IsInsideField(ProgressReporter *pr, FieldHandle input, FieldHandle& output, FieldHandle objfield)
{
  if (input.get_rep() == 0)
  {
    pr->error("IsInsideField: No input field");
    return (false);
  }

  if (objfield.get_rep() == 0)
  {
    pr->error("IsInsideField: No object field");
    return (false);
  }

  // no precompiled version available, so compile one

  FieldInformation fi(input);
  FieldInformation fo(input);
  FieldInformation fobj(objfield);
  
  if (fi.is_nonlinear())
  {
    pr->error("IsInsideField: This function has not yet been defined for non-linear elements");
    return (false);
  }
    
  fo.set_data_type("double");
  
  SCIRun::CompileInfoHandle ci = scinew CompileInfo(
    "IsInsideField."+fi.get_field_filename()+"."+fobj.get_field_filename()+".",
    "IsInsideFieldAlgo","IsInsideFieldAlgoT",
    fi.get_field_name() + "," + fo.get_field_name() + "," + fobj.get_field_name() );

  ci->add_include(TypeDescription::cc_to_h(__FILE__));
  ci->add_namespace("ModelCreation");
  ci->add_namespace("SCIRun");
  
  fi.fill_compile_info(ci);
  fo.fill_compile_info(ci);
  fobj.fill_compile_info(ci);
  
  // Handle dynamic compilation
  SCIRun::Handle<IsInsideFieldAlgo> algo;
  if(!(SCIRun::DynamicCompilation::compile(ci,algo,pr)))
  {
    pr->compile_error(ci->filename_);
//    SCIRun::DynamicLoader::scirun_loader().cleanup_failed_compile(ci);  
    return(false);
  }

  return(algo->IsInsideField(pr,input,output,objfield));
}

} // End namespace ModelCreation
