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


#ifndef Datatypes_LevelField_h
#define Datatypes_LevelField_h


#include "LevelMesh.h"
#include <Core/Datatypes/GenericField.h>
#include <Core/Containers/LockingHandle.h>
#include <Core/Math/MiscMath.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Util/NotFinished.h>
#include <Core/Util/Assert.h>
#include <Packages/Uintah/Core/Grid/Array3.h>
#include <Packages/Uintah/CCA/Components/MPM/Util/Matrix3.h>

#include <string>
#include <vector>
using std::string;
using std::vector;


/* //Specialization needed for InterpFunctor<LevelField<Matrix3> > */
/* #include <Core/Datatypes/FieldAlgo.h> */
/* namespace SCIRun { */
/*   using Uintah::Matrix3; */
/*   using Uintah::LevelField; */
/*   template<> InterpFunctor<LevelField<Matrix3> >::InterpFunctor() : */
/*     result_(double(0)) {} */
/* } // namespace SCIRun */

namespace Uintah {

using SCIRun::GenericField;
using SCIRun::LockingHandle;
using SCIRun::Interpolate;

// class Data must inherit from Packages/Uintah/Core/Grid/Array3 or
// this will not compile

template <class Data>
class LevelData : public vector<Array3<Data> > 
{
public:
  typedef Data value_type;
  //  typedef  iterator;
  typedef vector<Array3<Data> > parent;



  LevelData():vector<Array3<Data> >(){}
  LevelData(const LevelData& data) :
vector<Array3<Data> >(data) {} 
  virtual ~LevelData(){ }
  
  const value_type &operator[](typename LevelMesh::cell_index idx) const 
{ return parent::operator[](idx.patch_->getLevelIndex())
    [IntVector(idx.i_,idx.j_,idx.k_)]; } 
  const value_type &operator[](typename LevelMesh::face_index idx) const
{ return parent::operator[](idx.patch_->getLevelIndex())
    [IntVector(idx.i_,0,0)];}
const value_type &operator[](typename LevelMesh::edge_index idx) const 
{ return parent::operator[](idx.patch_->getLevelIndex())
    [IntVector(idx.i_, 0, 0)]; }
const value_type &operator[](typename LevelMesh::node_index idx) const
{ return parent::operator[](idx.patch_->getLevelIndex())
    [IntVector(idx.i_,idx.j_,idx.k_)]; }

value_type &operator[](typename LevelMesh::cell_index idx)
{ return parent::operator[](idx.patch_->getLevelIndex())
    [IntVector(idx.i_,idx.j_,idx.k_)]; } 
value_type &operator[](typename LevelMesh::face_index idx)
{ return parent::operator[](idx.patch_->getLevelIndex())
    [IntVector(idx.i_, 0, 0)]; }
value_type &operator[](typename LevelMesh::edge_index idx)
{ return parent::operator[](idx.patch_->getLevelIndex())
    [IntVector(idx.i_, 0, 0)]; }
value_type &operator[](typename LevelMesh::node_index idx)
{ return parent::operator[](idx.patch_->getLevelIndex())
    [IntVector(idx.i_,idx.j_,idx.k_)]; }

static const string type_name(int n = -1);
virtual const string get_type_name(int n = -1) const { return type_name(n); }

void resize(const LevelMesh::node_size_type &) {}
void resize(LevelMesh::edge_size_type) {}
void resize(LevelMesh::face_size_type) {}
void resize(const LevelMesh::cell_size_type &) {}

  class iterator
  {
  public:
    iterator(const vector<Array3<Data> >* data, IntVector index) 
      : it_( (*data)[0].begin() ), vit_(data->begin()), vitend_(data->end())
      {
	for(; vit_ != vitend_; vit_++){
	  IntVector low = (*vit_).getLowIndex();
	  IntVector high = (*vit_).getHighIndex();
	  if( index.x() >= low.x() && index.y() >= low.y() &&
	      index.z() >= low.z() && index.x() < high.x() &&
	      index.y() < high.y() && index.z() < high.z())
	  {
	    it_ = Array3<Data>::iterator( &(*vit_), index);
	    break;
	  }
	}
      }
    
    iterator(const iterator& iter) 
      : it_(iter.it_), vit_(iter.vit_), vitend_(iter.vitend_){}
    iterator& operator=(const iterator& it)
    { it_ = it.it_; vit_ = it.vit_; vitend_ == it.vitend_; }
    inline bool operator==(const iterator& it)
    { return it_ == it.it_ && vit_ == it.vit_ && vitend_ == it.vitend_;}
    inline bool operator!=(const iterator& it)
    { return !(operator==(it)); }
    inline Data& operator*() const {return *it_;}
    inline iterator& operator++(){
      if( ++it_ ==  (*vit_).end() ){
	if(++vit_ != vitend_){
	  it_ = (*vit_).begin();
	}
      }
      return *this;
    }
    inline iterator operator++(int){
      iterator result( *this );
      ++(*this);
      return result;
    }
  private:
    typename Array3<Data>::iterator it_;
    typename vector<Array3<Data> >::const_iterator vit_;
    typename vector<Array3<Data> >::const_iterator vitend_;
  };
  
iterator begin() { 
  const Array3<Data>& a = parent::operator[](0);
  return iterator(this, a.getLowIndex());
}
iterator end() {
  const Array3<Data>& a = parent::operator[](size()-1);
  iterator it(this, a.getHighIndex() - IntVector(1,1,1));
  return ++it;
}
  
};

    
template <class Data>
const string
LevelData<Data>::type_name(int n)
{
  ASSERT((n >= -1) && n <= 1);
  if (n == -1)
  {
    static const string name = type_name(0) +
      FTNS + "Array3" + FTNS + type_name(1) + FTNE + " " + FTNE;
    return name;
  }
  else if (n == 0)
  {
    return "LevelData";
  }
  else
  {
    return find_type_name((Data *)0);
  }
}
#define LEVELDATA_VERSION 1

template<class T>
void Pio(Piostream& stream, LevelData<T>& data)
{
#ifdef __GNUG__
#else
#endif

  stream.begin_class("Level", LEVELDATA_VERSION);
  NOT_FINISHED("LevelData::io");

  stream.end_class();
}


template <class Data>
class LevelField : public GenericField< LevelMesh, LevelData<Data>  > 
{ 

public:

  LevelField() :
    GenericField<LevelMesh, LevelData<Data> >() {}
  LevelField(Field::data_location data_at) :
    GenericField<LevelMesh, LevelData<Data> >(data_at) {}
  LevelField(LevelMeshHandle mesh, Field::data_location data_at) : 
    GenericField<LevelMesh, LevelData<Data> >(mesh, data_at) {}
  
  virtual ~LevelField(){}

  virtual LevelField<Data> *clone() const 
    { return new LevelField<Data>(*this); }
 
  static const string type_name(int n = -1);
  virtual const string get_type_name(int n = -1) const { return type_name(n); }
  static PersistentTypeID type_id;
  virtual void io(Piostream &stream);
  bool get_gradient(Vector &, Point &);
private:
  static Persistent* maker();
};


#define LEVELFIELD_VERSION 1

template <class Data>
Persistent* 
LevelField<Data>::maker()
{
  return scinew LevelField<Data>;
}

template <class Data>
PersistentTypeID
LevelField<Data>::type_id(type_name(),
		GenericField<LevelMesh, LevelData<Data> >::type_name(),
                maker); 

template <class Data>
void
LevelField<Data>::io(Piostream &stream)
{
  stream.begin_class(type_name(), LEVELFIELD_VERSION);
  GenericField<LevelMesh, LevelData<Data> >::io(stream);
  stream.end_class();                                                         
}


template <class Data>
const string
LevelField<Data>::type_name(int n)
{
  ASSERT((n >= -1) && n <= 1);
  if (n == -1)
  {
    static const string name = type_name(0) + FTNS + type_name(1) + FTNE;
    return name;

  }
  else if (n == 0)
  {
    return "LevelField";
  }
  else
  {
    return find_type_name((Data *)0);
  }
} 




//! compute the gradient g, at point p
/* template <>  */
/* bool LevelField<Matrix3>::get_gradient(Vector &, Point &p); */

/* template <>  */
/* bool LevelField<Vector>::get_gradient(Vector &, Point &p); */


template <class Data>
bool LevelField<Data>::get_gradient(Vector &g, Point &p) {
  // for now we only know how to do this for fields with scalars at the nodes
  ASSERT(data_at() == Field::NODE)
  ASSERT(type_name(1) == "double")
  LevelField<double> *lvd = dynamic_cast<LevelField<double> *>(this);

  mesh_handle_type mesh = get_typed_mesh();
  Vector pn=p-mesh->get_min();
  Vector diagonal = mesh->diagonal();
  int nx=mesh->get_nx();
  int ny=mesh->get_ny();
  int nz=mesh->get_nz();
  double diagx=diagonal.x();
  double diagy=diagonal.y();
  double diagz=diagonal.z();
  double x=pn.x()*(nx-1)/diagx;
  double y=pn.y()*(ny-1)/diagy;
  double z=pn.z()*(nz-1)/diagz;
  int ix0=(int)x;
  int iy0=(int)y;
  int iz0=(int)z;
  int ix1=ix0+1;
  int iy1=iy0+1;
  int iz1=iz0+1;
  if(ix0<0 || ix1>=nx)return false;
  if(iy0<0 || iy1>=ny)return false;
  if(iz0<0 || iz1>=nz)return false;
  double fx=x-ix0;
  double fy=y-iy0;
  double fz=z-iz0;
  double d000=lvd->value(LevelMesh::node_index(ix0,iy0,iz0));
  double d100=lvd->value(LevelMesh::node_index(ix1,iy0,iz0));
  double d010=lvd->value(LevelMesh::node_index(ix0,iy1,iz0));
  double d110=lvd->value(LevelMesh::node_index(ix1,iy1,iz0));
  double d001=lvd->value(LevelMesh::node_index(ix0,iy0,iz1));
  double d101=lvd->value(LevelMesh::node_index(ix1,iy0,iz1));
  double d011=lvd->value(LevelMesh::node_index(ix0,iy1,iz1));
  double d111=lvd->value(LevelMesh::node_index(ix1,iy1,iz1));
  double z00=Interpolate(d000, d001, fz);
  double z01=Interpolate(d010, d011, fz);
  double z10=Interpolate(d100, d101, fz);
  double z11=Interpolate(d110, d111, fz);
  double yy0=Interpolate(z00, z01, fy);
  double yy1=Interpolate(z10, z11, fy);
  double dx=(yy1-yy0)*(nx-1)/diagx;
  double x00=Interpolate(d000, d100, fx);
  double x01=Interpolate(d001, d101, fx);
  double x10=Interpolate(d010, d110, fx);
  double x11=Interpolate(d011, d111, fx);
  double y0=Interpolate(x00, x10, fy);
  double y1=Interpolate(x01, x11, fy);
  double dz=(y1-y0)*(nz-1)/diagz;
  double z0=Interpolate(x00, x01, fz);
  double z1=Interpolate(x10, x11, fz);
  double dy=(z1-z0)*(ny-1)/diagy;
  g = Vector(dx, dy, dz);
  return true;
}

}
#endif // Datatypes_LevelField_h
