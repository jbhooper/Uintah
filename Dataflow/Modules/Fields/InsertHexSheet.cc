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
 *  InsertHexSheet.cc:  Insert a layer of hexes.
 *
 *  Written by:
 *   Jason Shepherd
 *   Department of Computer Science
 *   University of Utah
 *   March 2006
 *
 *  Copyright (C) 2006 SCI Group
 */

#include <Dataflow/Network/Module.h>
#include <Core/Datatypes/Field.h>
#include <Dataflow/Network/Ports/FieldPort.h>
#include <Core/Datatypes/FieldInterface.h>
#include <Dataflow/Modules/Fields/InsertHexSheet.h>
#include <Core/Containers/StringUtil.h>
#include <iostream>

namespace SCIRun {

class InsertHexSheet : public Module
{
private:
  int       last_field_generation_;

public:
  InsertHexSheet(GuiContext* ctx);
  virtual ~InsertHexSheet();

  virtual void execute();
};


DECLARE_MAKER(InsertHexSheet)


InsertHexSheet::InsertHexSheet(GuiContext* ctx)
        : Module("InsertHexSheet", ctx, Filter, "FieldsCreate", "SCIRun"),
          last_field_generation_(0)
{
}

InsertHexSheet::~InsertHexSheet()
{
}

void InsertHexSheet::execute()
{
    // Get input fields.
  FieldIPort *hexfp = (FieldIPort *)get_iport("HexField");
  FieldHandle hexfieldhandle;
  if (!(hexfp->get(hexfieldhandle) && hexfieldhandle.get_rep()))
  {
    return;
  }

  FieldIPort *trifp = (FieldIPort *)get_iport("TriField");
  FieldHandle trifieldhandle;
  if (!(trifp->get(trifieldhandle) && trifieldhandle.get_rep()))
  {
    return;
  }

  if (last_field_generation_ == hexfieldhandle->generation &&
      last_field_generation_ == trifieldhandle->generation &&
      oport_cached( "IntersectField" )&&
      oport_cached( "Side1Field" )&&
      oport_cached( "Side2Field" ) )
  {
    // We're up to date, return.
    return;
  }
  last_field_generation_ = hexfieldhandle->generation;

  string ext = "";
  const TypeDescription *mtd = hexfieldhandle->mesh()->get_type_description();
  if (mtd->get_name().find("HexVolMesh") != string::npos)
  {
    ext = "Hex";
  }
  else
  {
    error( "Only HexVolFields are currently supported in the InsertHexSheet module.");
    return;
  }

  const TypeDescription *tri_mtd = trifieldhandle->mesh()->get_type_description();
  if (tri_mtd->get_name().find("TriSurfMesh") != string::npos)
  {
      //just checking... do nothing...
  }
  else
  {
    error( "Only TriSurfFields can be input to the InsertHexSheet module.");
    return;
  }

  const TypeDescription *ftd = hexfieldhandle->get_type_description();
  CompileInfoHandle ci = InsertHexSheetAlgo::get_compile_info(ftd, ext);
  Handle<InsertHexSheetAlgo> algo;
  if (!DynamicCompilation::compile(ci, algo, false, this))
  {
    error("Unable to compile InsertHexSheet algorithm.");
    return;
  }

  FieldHandle side1field, side2field, intersectfield;
  algo->execute( this, hexfieldhandle, trifieldhandle, 
                 side1field, side2field, intersectfield );
  
  FieldOPort *ofield_port1 = (FieldOPort *)get_oport("IntersectField");
  ofield_port1->send_and_dereference(intersectfield);  
  FieldOPort *ofield_port2 = (FieldOPort *)get_oport("Side1Field");
  ofield_port2->send_and_dereference(side1field);
  FieldOPort *ofield_port3 = (FieldOPort *)get_oport("Side2Field");
  ofield_port3->send_and_dereference(side2field);
}

CompileInfoHandle
InsertHexSheetAlgo::get_compile_info(const TypeDescription *fsrc,
			      string ext)
{
  // Use cc_to_h if this is in the .cc file, otherwise just __FILE__
  static const string include_path(TypeDescription::cc_to_h(__FILE__));
  const string template_class_name("InsertHexSheetAlgo" + ext);
  static const string base_class_name("InsertHexSheetAlgo");

  CompileInfo *rval = 
    scinew CompileInfo(template_class_name + "." +
		       fsrc->get_filename() + ".",
                       base_class_name, 
                       template_class_name, 
                       fsrc->get_name());

  // Add in the include path to compile this obj
  rval->add_include(include_path);
  fsrc->fill_compile_info(rval);

  return rval;
}

} // End namespace SCIRun

