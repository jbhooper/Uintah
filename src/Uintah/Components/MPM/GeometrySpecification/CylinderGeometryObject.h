#ifndef __CYLINDER_GEOMETRY_OBJECT_H__
#define __CYLINDER_GEOMETRY_OBJECT_H__

#include "GeometryObject.h"
#include <SCICore/Geometry/Point.h>
#include <Uintah/Grid/Box.h>


using SCICore::Geometry::Point;
using Uintah::Grid::Box;

namespace Uintah {
namespace Components {

class CylinderGeometryObject : public GeometryObject {

 public:

  enum AXIS {X = 1, Y = 2, Z = 3};

  CylinderGeometryObject(ProblemSpecP &);
  virtual ~CylinderGeometryObject();

  virtual bool inside(const Point &p) const;
  virtual Box getBoundingBox() const;
 
 private:
  AXIS  d_axis;
  Point d_origin;
  double d_length;
  double d_radius;
 
  

};

} // end namespace Uintah
} // end namespace Components

#endif // __CYLINDER_GEOMTRY_OBJECT_H__

// $Log$
// Revision 1.3  2000/04/20 22:37:13  jas
// Fixed up the GeometryObjectFactory.  Added findBlock() and findNextBlock()
// to ProblemSpec stuff.  This will iterate through all of the nodes (hopefully).
//
// Revision 1.2  2000/04/20 15:09:25  jas
// Added factory methods for GeometryObjects.
//
// Revision 1.1  2000/04/19 21:31:07  jas
// Revamping of the way objects are defined.  The different geometry object
// subtypes only do a few simple things such as testing whether a point
// falls inside the object and also gets the bounding box for the object.
// The constructive solid geometry objects:union,difference, and intersection
// have the same simple operations.
//
