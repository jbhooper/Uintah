#ifndef UINTAH_GRID_LEVEL_H
#define UINTAH_GRID_LEVEL_H

#include <Uintah/Grid/RefCounted.h>
#include <Uintah/Grid/Box.h>
#include <Uintah/Interface/DataWarehouseP.h>
#include <Uintah/Grid/GridP.h>
#include <Uintah/Grid/LevelP.h>
#include <Uintah/Grid/Handle.h>
#include <SCICore/Geometry/Point.h>
#include <SCICore/Geometry/BBox.h>
#include <SCICore/Geometry/IntVector.h>
#include <Uintah/Interface/ProblemSpec.h>
#include <Uintah/Interface/ProblemSpecP.h>
#include <Uintah/Grid/fixedvector.h>

#include <string>
#include <vector>

#define SELECT_GRID
//#define SELECT_RANGETREE

#ifdef SELECT_RANGETREE
class PatchRangeTree;
#endif

namespace Uintah {

   using SCICore::Geometry::Point;
   using SCICore::Geometry::Vector;
   using SCICore::Geometry::IntVector;
   class Patch;
   class Task;

   
/**************************************

CLASS
   Level
   
   Just a container class that manages a set of Patches that
   make up this level.

GENERAL INFORMATION

   Level.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
  
   Copyright (C) 2000 SCI Group

KEYWORDS
   Level

DESCRIPTION
   Long description...
  
WARNING
  
****************************************/

   class Level : public RefCounted {
   public:
      Level(Grid* grid, const Point& anchor, const Vector& dcell);
      virtual ~Level();
      
      void setPatchDistributionHint(const IntVector& patchDistribution);
      typedef std::vector<Patch*>::iterator patchIterator;
      typedef std::vector<Patch*>::const_iterator const_patchIterator;
      const_patchIterator patchesBegin() const;
      const_patchIterator patchesEnd() const;
      patchIterator patchesBegin();
      patchIterator patchesEnd();
      
      Patch* addPatch(const SCICore::Geometry::IntVector& lowIndex,
		      const SCICore::Geometry::IntVector& highIndex);
      
      Patch* addPatch(const SCICore::Geometry::IntVector& lowIndex,
		      const SCICore::Geometry::IntVector& highIndex,
		      int ID);

     //////////
     // Find a patch containing the point, return 0 if non exists
      Patch* getPatchFromPoint( const Point& );

      void finalizeLevel();
      void assignBCS(const ProblemSpecP& ps);
      
      int numPatches() const;
      long totalCells() const;

      void getIndexRange(SCICore::Geometry::BBox& b) const;
      void getSpatialRange(SCICore::Geometry::BBox& b) const;

      // it is better that IndexRange expressed by IntVector
      // instead of BBox, and the function is const. --tan
      void getIndexRange(IntVector& lowIndex,IntVector& highIndex) const;
      
      void performConsistencyCheck() const;
      GridP getGrid() const;

      Vector dCell() const {
	 return d_dcell;
      }
      Point getAnchor() const {
	 return d_anchor;
      }
      Point getNodePosition(const IntVector&) const;
      Point getCellPosition(const IntVector&) const;
      IntVector getCellIndex(const Point&) const;
      Point positionToIndex(const Point&) const;

      Box getBox(const IntVector&, const IntVector&) const;

      static const int MAX_PATCH_SELECT = 32;
      typedef fixedvector<const Patch*, MAX_PATCH_SELECT> selectType;

      void selectPatches(const IntVector&, const IntVector&,
			 selectType&) const;

      bool containsPoint(const Point&) const;
   private:
      Level(const Level&);
      Level& operator=(const Level&);
      
      std::vector<Patch*> d_patches;
     
      Grid* grid;
      Point d_anchor;
      Vector d_dcell;
      bool d_finalized;
      IntVector d_patchDistribution;

#ifdef SELECT_GRID
      IntVector d_idxLow;
      IntVector d_idxHigh;
      IntVector d_idxSize;
      IntVector d_gridSize;
      vector<int> d_gridStarts;
      vector<Patch*> d_gridPatches;
#else
#ifdef SELECT_RANGETREE
      PatchRangeTree* d_rangeTree;
#endif
#endif
   };
   
} // end namespace Uintah

//
// $Log$
// Revision 1.21.4.3  2000/10/25 20:36:31  witzel
// Added RangeTree option for selectPatches implementation.
//
// Revision 1.21.4.2  2000/10/10 05:28:08  sparker
// Added support for NullScheduler (used for profiling taskgraph overhead)
//
// Revision 1.21.4.1  2000/10/07 06:10:36  sparker
// Optimized implementation of Level::selectPatches
// Cured g++ warnings
//
// Revision 1.21  2000/09/25 20:37:42  sparker
// Quiet g++ compiler warnings
// Work around g++ compiler bug instantiating vector<NCVariable<Vector> >
// Added computeVariableExtents to (eventually) simplify data warehouses
//
// Revision 1.20  2000/08/22 18:36:40  bigler
// Added functionality to get a cell's position with the index.
//
// Revision 1.19  2000/07/11 19:32:16  kuzimmer
// Added getPatchFromPoint(const Point& p)
//
// Revision 1.18  2000/07/11 15:53:45  tan
// Changed the following to be const member function.
//       void getIndexRange(SCICore::Geometry::BBox& b) const;
//       void getSpatialRange(SCICore::Geometry::BBox& b) const;
//
// Revision 1.17  2000/07/07 03:08:14  tan
// It is better that IndexRange expressed by IntVector instead of BBox,
// and the function is const.
//
// Revision 1.16  2000/06/27 22:49:03  jas
// Added grid boundary condition support.
//
// Revision 1.15  2000/06/23 19:20:19  jas
// Added in the early makings of Grid bcs.
//
// Revision 1.14  2000/06/15 21:57:17  sparker
// Added multi-patch support (bugzilla #107)
// Changed interface to datawarehouse for particle data
// Particles now move from patch to patch
//
// Revision 1.13  2000/05/30 20:19:29  sparker
// Changed new to scinew to help track down memory leaks
// Changed region to patch
//
// Revision 1.12  2000/05/20 08:09:22  sparker
// Improved TypeDescription
// Finished I/O
// Use new XML utility libraries
//
// Revision 1.11  2000/05/20 02:36:06  kuzimmer
// Multiple changes for new vis tools and DataArchive
//
// Revision 1.10  2000/05/15 19:39:47  sparker
// Implemented initial version of DataArchive (output only so far)
// Other misc. cleanups
//
// Revision 1.9  2000/05/10 20:02:59  sparker
// Added support for ghost cells on node variables and particle variables
//  (work for 1 patch but not debugged for multiple)
// Do not schedule fracture tasks if fracture not enabled
// Added fracture directory to MPM sub.mk
// Be more uniform about using IntVector
// Made patches have a single uniform index space - still needs work
//
// Revision 1.8  2000/04/26 06:48:49  sparker
// Streamlined namespaces
//
// Revision 1.7  2000/04/12 23:00:47  sparker
// Starting problem setup code
// Other compilation fixes
//
// Revision 1.6  2000/03/22 00:32:12  sparker
// Added Face-centered variable class
// Added Per-patch data class
// Added new task constructor for procedures with arguments
// Use Array3Index more often
//
// Revision 1.5  2000/03/21 01:29:42  dav
// working to make MPM stuff compile successfully
//
// Revision 1.4  2000/03/17 18:45:42  dav
// fixed a few more namespace problems
//
// Revision 1.3  2000/03/16 22:07:59  dav
// Added the beginnings of cocoon docs.  Added namespaces.  Did a few other coding standards updates too
//
//

#endif




