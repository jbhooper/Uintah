/*
  The contents of this file are subject to the University of Utah Public
  License (the "License"); you may not use this file except in compliance
  with the License.
  
  Software distributed under the License is distributed on an "AS IS"
  basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
  License for the specific language governing rights and limitations under
  the License.
  
  The Original Source Code is SCIRun, released March 12, 2001.
  
  The Original Source Code was developed by the University of Utah.
  Portions created by UNIVERSITY are Copyright (C) 2001, 1994 
  University of Utah. All Rights Reserved.
*/

//    File   : SelectField.h
//    Author : Michael Callahan
//    Date   : August 2001

#include <Core/Persistent/Pstreams.h>
#include <Dataflow/Network/Module.h>
#include <Dataflow/Ports/FieldPort.h>
#include <Dataflow/Ports/GeometryPort.h>
#include <Dataflow/Modules/Fields/SelectField.h>
#include <Dataflow/Widgets/BoxWidget.h>
#include <Core/Datatypes/Clipper.h>
#include <Core/GuiInterface/GuiVar.h>

#include <iostream>
#include <stdio.h>

namespace SCIRun {

class SelectField : public Module
{
private:
  FieldHandle output_field_;
  CrowdMonitor widget_lock_;
  BoxWidget *box_;

  GuiInt value_;
  GuiInt mode_;  // 0 nothing 1 accumulate 2 replace

  int  last_generation_;
  BBox last_bounds_;

public:
  SelectField(const string& id);
  virtual ~SelectField();
  virtual void execute();
};


extern "C" Module* make_SelectField(const string& id)
{
  return new SelectField(id);
}

SelectField::SelectField(const string& id)
  : Module("SelectField", id, Filter, "Fields", "SCIRun"),
    widget_lock_("SelectField widget lock"),
    value_("stampvalue", id, this),
    mode_("runmode", id, this),
    last_generation_(0)
{
  box_ = scinew BoxWidget(this, &widget_lock_, 1.0, false, false);
}



SelectField::~SelectField()
{
}



void
SelectField::execute()
{
  // Get input field.
  FieldIPort *ifp = (FieldIPort *)get_iport("Input Field");
  FieldHandle ifieldhandle;
  if (!ifp) {
    postMessage("Unable to initialize "+name+"'s iport\n");
    return;
  }
  if (!(ifp->get(ifieldhandle) && ifieldhandle.get_rep()))
  {
    return;
  }

  bool forward_p = false;

  if (output_field_.get_rep() == NULL ||
      last_generation_ != ifieldhandle->generation)
  {
    const TypeDescription *mtd = ifieldhandle->mesh()->get_type_description();
    const TypeDescription *ftd = ifieldhandle->get_type_description();
    CompileInfo *ci = SelectFieldCreateAlgo::get_compile_info(mtd, ftd);
    DynamicAlgoHandle algo_handle;
    if (! DynamicLoader::scirun_loader().get(*ci, algo_handle))
    {
      cout << "Could not compile algorithm." << std::endl;
      return;
    }
    SelectFieldCreateAlgo *algo =
      dynamic_cast<SelectFieldCreateAlgo *>(algo_handle.get_rep());
    if (algo == 0)
    {
      cout << "Could not get algorithm." << std::endl;
      return;
    }
    output_field_ =
      algo->execute(ifieldhandle->mesh(), ifieldhandle->data_at());

    last_generation_ = ifieldhandle->generation;

    BBox obox = output_field_->mesh()->get_bounding_box();
    if (!(last_bounds_.valid() && obox.valid() &&
	  obox.min() == last_bounds_.min() &&
	  obox.max() == last_bounds_.max()))
    {
      // update the widget
      const BBox bbox = output_field_->mesh()->get_bounding_box();
      const Point &bmin = bbox.min();
      const Point &bmax = bbox.max();
#if 0
      const Point center = bmin + Vector(bmax - bmin) * 0.5;
      const Point right = center + Vector((bmax.x()-bmin.x())/2.,0,0);
      const Point down = center + Vector(0,(bmax.y()-bmin.y())/2.,0);
      const Point in = center + Vector(0,0,(bmax.z()-bmin.z())/2.);
#else
      const Point center = bmin + Vector(bmax - bmin) * 0.25;
      const Point right = center + Vector((bmax.x()-bmin.x())/4.0, 0, 0);
      const Point down = center + Vector(0, (bmax.y()-bmin.y())/4.0, 0);
      const Point in = center + Vector(0, 0, (bmax.z()-bmin.z())/4.0);
#endif
      const double l2norm = (bmax - bmin).length();

      box_->SetScale(l2norm * 0.015);
      box_->SetPosition(center, right, down, in);

      GeomGroup *widget_group = scinew GeomGroup;
      widget_group->add(box_->GetWidget());

      GeometryOPort *ogport=0;
      ogport = (GeometryOPort*)get_oport("Selection Widget");
      if (!ogport) {
	postMessage("Unable to initialize "+name+"'s oport\n");
	return;
      }
      ogport->addObj(widget_group, "SelectField Selection Widget",
		     &widget_lock_);
      ogport->flushViews();

      last_bounds_ = obox;
    }
    forward_p = true;
  }

  if (mode_.get() == 1 || mode_.get() == 2)
  {
    output_field_.detach();
    const TypeDescription *oftd = output_field_->get_type_description();
    const TypeDescription *oltd = output_field_->data_at_type_description();
    CompileInfo *ci = SelectFieldFillAlgo::get_compile_info(oftd, oltd);
    DynamicAlgoHandle algo_handle;
    if (! DynamicLoader::scirun_loader().get(*ci, algo_handle))
    {
      cout << "Could not compile algorithm." << std::endl;
      return;
    }
    SelectFieldFillAlgo *algo =
      dynamic_cast<SelectFieldFillAlgo *>(algo_handle.get_rep());
    if (algo == 0)
    {
      cout << "Could not get algorithm." << std::endl;
      return;
    }
    bool replace_p = false;
    if (mode_.get() == 2) { replace_p = true; }

    ClipperHandle clipper = box_->get_clipper();
    algo->execute(output_field_, clipper, value_.get(), replace_p, 0);

    forward_p = true;
  }

  if (forward_p)
  {
    FieldOPort *ofield_port = (FieldOPort *)get_oport("Output Field");
    if (!ofield_port) {
      postMessage("Unable to initialize "+name+"'s oport\n");
      return;
    }
    
    ofield_port->send(output_field_);
  }
}



CompileInfo *
SelectFieldCreateAlgo::get_compile_info(const TypeDescription *msrc,
					const TypeDescription *fsrc)
{
  // use cc_to_h if this is in the .cc file, otherwise just __FILE__
  static const string include_path(TypeDescription::cc_to_h(__FILE__));
  static const string template_class_name("SelectFieldCreateAlgoT");
  static const string base_class_name("SelectFieldCreateAlgo");

  const string::size_type loc = fsrc->get_name().find_first_of('<');
  const string fout = fsrc->get_name().substr(0, loc) + "<int> ";

  CompileInfo *rval = 
    scinew CompileInfo(template_class_name + "." +
		       msrc->get_filename() + "." +
		       to_filename(fout) + ".",
                       base_class_name, 
                       template_class_name, 
                       msrc->get_name() + ", " + fout);

  // Add in the include path to compile this obj
  rval->add_include(include_path);
  fsrc->fill_compile_info(rval);
  return rval;
}


CompileInfo *
SelectFieldFillAlgo::get_compile_info(const TypeDescription *fsrc,
				      const TypeDescription *lsrc)
{
  // use cc_to_h if this is in the .cc file, otherwise just __FILE__
  static const string include_path(TypeDescription::cc_to_h(__FILE__));
  static const string template_class_name("SelectFieldFillAlgoT");
  static const string base_class_name("SelectFieldFillAlgo");

  CompileInfo *rval = 
    scinew CompileInfo(template_class_name + "." +
		       fsrc->get_filename() + "." +
		       lsrc->get_filename() + ".",
                       base_class_name, 
                       template_class_name, 
                       fsrc->get_name() + ", " + lsrc->get_name());

  // Add in the include path to compile this obj
  rval->add_include(include_path);
  fsrc->fill_compile_info(rval);
  return rval;
}



} // End namespace SCIRun
