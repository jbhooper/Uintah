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

#ifndef TEXCUTTINGPLANES_H
#define TEXCUTTINGPLANES_H
/*
 * TexCuttingPlanes.cc
 *
 * Simple interface to volume rendering stuff
 */

#include <Dataflow/Network/Module.h>
#include <Core/Datatypes/ColorMap.h>
#include <Dataflow/Ports/ColorMapPort.h>
#include <Dataflow/Ports/GeometryPort.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Parts/GuiVar.h>
#include <Core/Thread/CrowdMonitor.h>
#include <Dataflow/Widgets/PointWidget.h>

#include <Core/Datatypes/GLVolumeRenderer.h>
#include <Dataflow/Ports/GLTexture3DPort.h>
#include <Core/Datatypes/GLTexture3D.h>

namespace SCIRun {

class GeomObj;



class TexCuttingPlanes : public Module {

public:
  TexCuttingPlanes( const string& id);

  virtual ~TexCuttingPlanes();
  virtual void widget_moved(int last);    
  virtual void execute();
  void tcl_command( TCLArgs&, void* );

private:

   
  GLTexture3DHandle tex;

  ColorMapIPort* incolormap;
  GLTexture3DIPort* intexture;
  GeometryOPort* ogeom;
   
  CrowdMonitor control_lock; 
  PointWidget *control_widget;
  GeomID control_id;


  int cmap_id;  // id associated with color map...
  

  GuiInt drawX;
  GuiInt drawY;
  GuiInt drawZ;
  GuiInt drawView;
  GuiInt interp_mode;

  GLVolumeRenderer* volren;
  Vector ddv;
  double ddview;
};


} // End namespace SCIRun

#endif
