
#ifndef SUBMATERIAL_H
#define SUBMATERIAL_H

#include <Packages/rtrt/Core/Material.h>
#include <vector>

namespace rtrt {

class SubMaterial : public Material
{

 protected:

  vector<Material*> materials_;

 public:

  SubMaterial() {}
  virtual ~SubMaterial() {}

  Material *operator[](unsigned index) { 
    if (index>=materials_.size())
      return 0;
    else
      return materials_[index]; 
  }

  add_material(Material *mat) { materials_.push_back(mat); }

  virtual void shade(Color& result, const Ray& ray,
                     const HitInfo& hit, int depth, 
                     double atten, const Color& accumcolor,
                     Context* cx) { 
    double r = drand48();
    result = Color(r,r,r); 
  }
};

} // end namespace

#endif
