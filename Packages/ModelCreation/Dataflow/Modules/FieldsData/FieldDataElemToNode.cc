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
 *  FieldDataElemToNode.cc:
 *
 *  Written by:
 *   jeroen
 *   TODAY'S DATE HERE
 *
 */

#include <Dataflow/Network/Module.h>
#include <Core/Malloc/Allocator.h>

#include <Dataflow/Network/Ports/FieldPort.h>
#include <Dataflow/Network/Ports/MatrixPort.h>

#include <Core/Datatypes/Field.h>
#include <Core/Datatypes/Matrix.h>

#include <Core/Algorithms/Fields/FieldsAlgo.h>

namespace ModelCreation {

using namespace SCIRun;

class FieldDataElemToNode : public Module {
  public:
    FieldDataElemToNode(GuiContext*);
    virtual void execute();

  private:
    GuiString method_;
};


DECLARE_MAKER(FieldDataElemToNode)
FieldDataElemToNode::FieldDataElemToNode(GuiContext* ctx)
  : Module("FieldDataElemToNode", ctx, Source, "FieldsData", "ModelCreation"),
    method_(get_ctx()->subVar("method"))  
{
}

void
 FieldDataElemToNode::execute()
{   
  FieldHandle input;
  FieldHandle output;

  if (!(get_input_handle("Field",input,true))) return;
  
  std::string method = method_.get();
  SCIRunAlgo::FieldsAlgo algo(this);
  

  if(!(algo.FieldDataElemToNode(input,output,method))) return;
 
  send_output_handle("Field",output,true); 
}

} // End namespace ModelCreation


