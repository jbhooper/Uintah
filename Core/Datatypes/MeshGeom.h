//  MeshGeom.h - A base class for regular geometries with alligned axes
//  Written by:
//   Eric Kuehne
//   Department of Computer Science
//   University of Utah
//   April 2000
//  Copyright (C) 2000 SCI Institute


#ifndef SCI_project_MeshGeom_h
#define SCI_project_MeshGeom_h 1

#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Point.h>
#include <Core/Geom/GeomTriangles.h>
#include <Core/Geom/GeomPolyline.h>
#include <Core/Datatypes/SurfaceGeom.h>
#include <Core/Containers/LockingHandle.h>
#include <Core/Math/MiscMath.h>
#include <Core/Util/DebugStream.h>

#include <list>
#include <vector>
#include <string>


namespace SCIRun {

using std::list;
using std::vector;
using std::string;

class MeshGeom:public SurfaceGeom
{
public:

  MeshGeom();
  virtual ~MeshGeom() {}
  
  virtual string getInfo();
  virtual string getTypeName(int=0);
  
  ///////////
  // Persistent representation...
  virtual void io(Piostream&);
  static PersistentTypeID type_id;
  static string typeName(int);
protected:

  vector<list<int> > cell_;

private:
  static DebugStream dbg;
};

} // End namespace SCIRun


#endif
