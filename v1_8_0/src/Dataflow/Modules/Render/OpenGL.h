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


/*
 *  OpenGL.cc: Render geometry using opengl
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   December 1994
 *
 *  Copyright (C) 1994 SCI Group
 */
//milan was here

#ifndef SCI_project_module_OpenGL_h
#define SCI_project_module_OpenGL_h

#include <tcl.h>
#include <tk.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <map>
#include <vector>

#include <sci_defs.h>

#include <Core/Geom/GeomObj.h>
#include <Core/Util/Timer.h>
#include <Core/Geom/GeomObj.h>
#include <Core/Geom/GeomOpenGL.h>
#include <Core/Geom/Light.h>
#include <Core/Geom/Lighting.h>
#include <Core/Geom/RenderMode.h>
#include <Core/Geom/View.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Trig.h>
#include <Core/Datatypes/Image.h>
#include <Dataflow/Ports/GeometryPort.h>
#include <Dataflow/Modules/Render/Ball.h>
#include <Dataflow/Modules/Render/ViewWindow.h>
#include <Dataflow/Modules/Render/Viewer.h>
#include <Core/Thread/FutureValue.h>
#include <Core/Thread/Runnable.h>
#include <Core/Thread/Thread.h>

#include <Core/Geom/GeomCone.h>
#include <Core/Geom/GeomCylinder.h>
#include <Core/Geom/GeomSphere.h>
#include <Core/Geom/GeomTri.h>
#include <Core/Geom/GeomText.h>

#ifdef HAVE_MPEG
#include <mpege.h>
#endif // HAVE_MPEG

#ifdef __sgi
#include <X11/extensions/SGIStereo.h>
#endif

namespace SCIRun {

using std::ostringstream;
using std::ofstream;
using std::vector;




class OpenGLHelper;
  class GuiArgs;

/* #define DO_REDRAW 0 */
/* #define DO_PICK 1 */
/* #define DO_GETDATA 2 */
/* #define REDRAW_DONE 4 */
/* #define PICK_DONE 5 */
/* #define DO_IMAGE 6 */
/* #define DO_HIRES_IMAGE 7 */
/* #define IMAGE_DONE 8 */

struct GetReq {
    int datamask;
    FutureValue<GeometryData*>* result;
    GetReq(int, FutureValue<GeometryData*>* result);
    GetReq();
};

struct ImgReq {
  string name;
  string type;
  int resx;
  int resy;
  ImgReq( const string& n, const string& t, int rx, int ry);
  ImgReq();
};

struct Frustum {
  double znear;
  double zfar;
  double left;
  double right;
  double bottom;
  double top;
  double width;
  double height;
};

struct HiRes {
  double nrows;
  double ncols;
  int row;
  int col;
  int resx;
  int resy;
};

class OpenGL {
  
protected:
  GuiInterface* gui;
  bool compute_depth(ViewWindow* viewwindow, const View& view, double& near, double& far);
  int xres, yres;
  Runnable *helper;
    Tk_Window tkwin;
    Window win;
    Display* dpy;
    GLXContext cx;
    int maxlights;
    DrawInfoOpenGL* drawinfo;
    WallClockTimer fpstimer;
    double current_time;

    int old_stereo;
    GLuint imglist;
    void make_image();

    void redraw_obj(Viewer* viewer, ViewWindow* viewwindow, GeomObj* obj);
    void pick_draw_obj(Viewer* viewer, ViewWindow* viewwindow, GeomObj* obj);
    //OpenGLHelper* helper;
    string my_openglname;
    vector<XVisualInfo*> visuals;

   /* Call this each time an mpeg frame is generated. */
  void StartMpeg(const string& fname);
  void AddMpegFrame();
  void EndMpeg();
#ifdef __sgi
#if (_MIPS_SZPTR == 64)
#else
  //CLhandle compressorHdl;
  //int compressedBufferSize;
  //CLbufferHdl compressedBufferHdl;
  //ofstream os;
#endif
#endif
  bool do_hi_res;
  bool encoding_mpeg;
#ifdef HAVE_MPEG
  FILE *output;
  MPEGe_options options;
#endif // HAVE_MPEG

public:
  static bool query(GuiInterface* gui);
    OpenGL(GuiInterface* gui);
     ~OpenGL();
     void redraw(Viewer*, ViewWindow*, double tbeg, double tend,
      int ntimesteps, double frametime);
    void real_get_pick(Viewer*, ViewWindow*, int, int, GeomObj*&, GeomPick*&, int&);
     void get_pick(Viewer*, ViewWindow*, int, int,
      GeomObj*&, GeomPick*&, int& );
     void dump_image(const string& fname,
			    const string& type = "raw");
     void put_scanline(int y, int width, Color* scanline, int repeat=1);

   void saveImage(const string& fname,
			 const string& type = "ppm", int x=640, int y=512);
  void render_and_save_image(int x, int y,
			     const string& fname,
			     const string& type = "ppm");
  void setFrustumToWindowPortion();
  void deriveFrustum();

// HACK -- support data for get_pixel_depth
float  pixel_depth_data[1310720];  // assume no screen is > 1280x1024
GLdouble get_depth_model[16];
GLdouble get_depth_proj[16];
GLint    get_depth_view[4];


    // compute world space point under cursor (x,y).  If successful,
    // set 'p' to that value & return true.  Otherwise, return false.
     int    pick_scene(int x, int y, Point *p);
     void kill_helper();

    string myname;
     void redraw_loop();
    Mailbox<int> send_mb;
    Mailbox<int> recv_mb;
    Mailbox<GetReq> get_mb;
    Mailbox<ImgReq> img_mb;
    
    Frustum frustum;
    HiRes hi_res;

    Viewer* viewer;
    ViewWindow* viewwindow;
    double tbeg;
    double tend;
    int nframes;
    double framerate;
    void redraw_frame();

    int send_pick_x;
    int send_pick_y;

    GeomObj* ret_pick_obj;
    GeomPick* ret_pick_pick;
    int ret_pick_index;

     void listvisuals(GuiArgs&);
     void setvisual(const string&, int i, int width, int height);

    View lastview;
    double znear, zfar;

    GeomCappedCylinder* stylusCylinder[2];
    GeomTri* stylusTriangle[4];
    GeomSphere* pinchSphere;
    Material* stylusMaterial[16], *pinchMaterial;
    
    GeomText* pinchText[2];
    GeomCappedCylinder* pinchCylinder[4];

    Thread *helper_thread;
    bool dead;

    // these functions were added to clean things up a bit...

//protected:
     void getData(int datamask,
      FutureValue<GeometryData*>* result);
  
protected:
    
     void real_getData(int datamask,
      FutureValue<GeometryData*>* result);
  
};

} // End namespace SCIRun

#endif // SCI_project_module_OpenGL_h
