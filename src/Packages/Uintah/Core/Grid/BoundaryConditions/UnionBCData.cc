#include <Packages/Uintah/Core/Grid/BoundaryConditions/UnionBCData.h>
#include <Core/Geometry/Point.h>
#include <Packages/Uintah/Core/Grid/Box.h>
#include <Packages/Uintah/Core/Grid/BoundaryConditions/BoundCondFactory.h>
#include <Packages/Uintah/Core/Grid/Variables/UnionIterator.h>
#include <Packages/Uintah/Core/Grid/Level.h>
#include <Core/Malloc/Allocator.h>
#include <iostream>
#include <algorithm>

using namespace SCIRun;
using namespace Uintah;

UnionBCData::UnionBCData()
{
}

UnionBCData::~UnionBCData()
{
  for (vector<BCGeomBase*>::const_iterator bc = child.begin();
       bc != child.end(); ++bc)
    delete (*bc);
  
  child.clear();
}

UnionBCData::UnionBCData(const UnionBCData& mybc)
{

  vector<BCGeomBase*>::const_iterator itr;
  for (itr=mybc.child.begin(); itr != mybc.child.end(); ++itr)
    child.push_back((*itr)->clone());
  
}

UnionBCData& UnionBCData::operator=(const UnionBCData& rhs)
{
  BCGeomBase::operator=(rhs);

  if (this == &rhs)
    return *this;

  // Delete the lhs
  vector<BCGeomBase*>::const_iterator itr;
  for(itr=child.begin(); itr != child.end();++itr)
    delete *itr;

  child.clear();
  
  // copy the rhs to the lhs
  for (itr=rhs.child.begin(); itr != rhs.child.end();++itr)
    child.push_back((*itr)->clone());
  
  return *this;
}


bool UnionBCData::operator==(const BCGeomBase& rhs) const
{
  const UnionBCData* p_rhs = 
    dynamic_cast<const UnionBCData*>(&rhs);

  if (p_rhs == NULL)
    return false;
  else {
    if (this->child.size() != p_rhs->child.size())
      return false;

    return equal(this->child.begin(),this->child.end(),p_rhs->child.begin());
  }

}

UnionBCData* UnionBCData::clone()
{
  return scinew UnionBCData(*this);

}

void UnionBCData::addBCData(BCData& bc)
{
}

void UnionBCData::addBC(BoundCondBase* bc)
{
  
}

void UnionBCData::addBCData(BCGeomBase* bc)
{
  child.push_back(bc);
}

void UnionBCData::getBCData(BCData& bc) const
{
  child[0]->getBCData(bc);
}

bool UnionBCData::inside(const Point &p) const 
{
  for (vector<BCGeomBase*>::const_iterator i = child.begin(); i != child.end();
       ++i){
    if ((*i)->inside(p))
      return true;
  }
  return false;
}

void UnionBCData::print()
{
  cout << "Geometry type = " << typeid(this).name() << endl;
  for (vector<BCGeomBase*>::const_iterator i = child.begin(); i != child.end();
       ++i)
    (*i)->print();

}


void UnionBCData::determineIteratorLimits(Patch::FaceType face, 
					  const Patch* patch, 
					  vector<Point>& test_pts)
{
#if 0
  cout << "UnionBC determineIteratorLimits()" << endl;
#endif

  for (vector<BCGeomBase*>::const_iterator bc = child.begin();
       bc != child.end(); ++bc) {
    (*bc)->determineIteratorLimits(face,patch,test_pts);
  }
  
  UnionIterator cells,nodes;

  for (vector<BCGeomBase*>::const_iterator bc = child.begin();
       bc != child.end(); ++bc) {
    Iterator cell_itr,node_itr;
    (*bc)->getCellFaceIterator(cell_itr);
    (*bc)->getNodeFaceIterator(node_itr);
    Iterator base_ci(cells),base_ni(nodes);
    cells = UnionIterator(base_ci,cell_itr);
    nodes = UnionIterator(base_ni,node_itr);
  }

  d_cells = UnionIterator(cells);   
  d_nodes = UnionIterator(nodes); 

}

