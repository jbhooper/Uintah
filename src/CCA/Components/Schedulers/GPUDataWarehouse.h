/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/* GPU DataWarehouse device & host access*/

#ifndef GPU_DW_H
#define GPU_DW_H

#include <sci_defs/cuda_defs.h>
#include <Core/Grid/Variables/GPUVariable.h>
#include <Core/Grid/Variables/GPUGridVariable.h>
#include <Core/Grid/Variables/GPUReductionVariable.h>
#include <Core/Grid/Variables/GridVariableBase.h>
#include <Core/Grid/Variables/GPUPerPatch.h>

#include <map> //for host code only.
#include <string>
#include <vector>
#include <Core/Thread/CrowdMonitor.h>

#define MAX_ITEM 10000000
#define MAX_GHOST_CELLS 30000
#define MAX_MATERIALS 20
#define MAX_LVITEM 20

#define MAX_LABEL   20 //How big a particular label can be.

namespace Uintah {

enum materialType {
  IDEAL_GAS = 0/*,
  HARD_SPHERE_GAS_EOS = 1,
  TST_EOS = 2,
  JWL_EOS = 3,
  JWLC_EOS = 4,
  MURNAGHAN_EOS = 5,
  BIRCH_MURNAGHAN_EOS = 6,
  GRUNEISEN_EOS = 7,
  TILLOTSON_EOS = 8,
  THOMSEN_HARTKA_WATER_EOS = 9,
  KNAUSS_SEA_WATER_EOS = 10,
  KUMARI_DASS_EOS = 11 */
};

class GPUDataWarehouse {

public:

  GPUDataWarehouse(int id, void * placementNewBuffer): allocateLock("allocate lock"), varLock("var lock"){
    d_device_id = id;
    objectSizeInBytes = 0;
    this->placementNewBuffer = placementNewBuffer;
    init();

  }
  virtual ~GPUDataWarehouse() {};

  enum GhostType {
    None,
    AroundNodes,
    AroundCells,
    AroundFacesX,
    AroundFacesY,
    AroundFacesZ,
    AroundFaces,
    numGhostTypes // 7
  };

  struct materialItem {
    materialType material;
    char       simulationType[MAX_LABEL]; //ICE, MPM, etc.  Currently unused
  };


  //The dataItem can hold two kinds of data.
  //The first is information related to a regular data variable.
  //The second is information indicating how a ghost cell should be copied from one var to another
  //The two types of information are stored in one struct to allow us to make the size of the GPU
  //data warehouse dynamic.  For small problems the GPUDW can be small, and large problems it can
  //be large.
  //The biggest problem is that multiple CPU threads will be adding to the size of the GPUDW
  //IF we had two separate structs, then we would have to let both grow independently and
  //it would require two copies to the GPU DW instead of one.
  //So the solution is to allocate a large buffer of possible GPUDW data in init_device(),
  //one on the host RAM and device RAM
  //Then on the CPU side a thread running a task will collect various dataItems, and just before
  //it sends it to the GPU DW, it will dump its results into the host side buffer (using locks).
  //Then it can either copy in only as much of the GPUDW as needed, instead of the entire buffer.

  struct VarItem {
    GhostType       gtype;
    unsigned int    numGhostCells;
    bool            validOnGPU; //true if the GPU copy is the "live" copy and not an old version of the data. //TODO: Remove me
    bool            validOnCPU; //true if the CPU copy is the current "live" copy. (It's possible to be both.)  //TODO: Remove me
    bool            queueingOnGPU; //true if we've created the variable but we haven't yet set validOnGPU to true. //TODO: Remove me

  };

  struct GhostItem {
    //This assumes the ghost cell is already in the GPU's memory
    //We only need to know the source patch, the destination patch
    //the material number, and the shared coordinates.
    int3            sharedLowCoordinates;
    int3            sharedHighCoordinates;

    //Wasatch has virtual patches, which come as a result of periodic boundary cells,
    //which wrap around on each other.  (Like a label wrapping around a soup can, but
    //for all boundary faces of the grid).  While they're technically sharing the same
    //coordinates (just wrapped around once), from our perspective we need their actual indexes
    //this offset helps us get that.
    int3 virtualOffset;

    //So we can look up the size and offset information in the d_varDB
    int dest_varDB_index;  //Will be set to -1 if this d_varDB item isn't a ghost cell item.

    //Ghost cells can come from two locations.  The first way is that they are already in the
    //d_varDB (such as same GPU->same GPU ghost cell, or CPU ghost cell->GPU).
    //The other way is if the data comes from another GPU on the same node,
    //and if so, then the data will be found in d_ghostCellDB instead of d_varDB.
    //Fortunately it turns out that this isn't needed.
    //int source_varDB_index;

  };

  struct dataItem {   // flat array  //57 BYTES  combined 110 bytes

    char            label[MAX_LABEL];
    int             domainID;
    int             matlIndex;
    int3            var_offset;
    int3            var_size;
    size_t          num_elems;
    void*           var_ptr;
    unsigned int    xstride;
    VarItem         varItem;
    GhostItem       ghostItem;

  };


  struct contiguousArrayInfo {
    void * allocatedDeviceMemory;
    void * allocatedHostMemory;
    size_t sizeOfAllocatedMemory;
    size_t assignedOffset;
    size_t copiedOffset;
    //The default constructor
    contiguousArrayInfo() {
      allocatedDeviceMemory = NULL;
      allocatedHostMemory = NULL;
      sizeOfAllocatedMemory = 0;
      assignedOffset = 0; //To keep up to the point where data has been "put".  Computes data will be assigned
      copiedOffset = 0; //To keep up to the point where data will need to be copied.  Required data will be copied
    }
    //Overloaded constructor
    contiguousArrayInfo(double * allocatedDeviceMemory,
                        double * allocatedHostMemory,
                        size_t sizeOfAllocatedMemory) {
      this->allocatedDeviceMemory = allocatedDeviceMemory;
      this->allocatedHostMemory = allocatedHostMemory;
      this->sizeOfAllocatedMemory = sizeOfAllocatedMemory;
      assignedOffset = 0; //To keep up to the point where data has been "put"
      copiedOffset = 0; //To keep up to the point where data has been copied.
    }
  };


  struct allVarPointersInfo {
    void*      device_ptr;   //Where it is on the device
    int3       device_offset;
    int3       device_size;
    void*      host_contiguousArrayPtr;  //Use this address only if partOfContiguousArray is set to true.
    GridVariableBase* gridVar;  //TODO: Rename, not everything is a gridvar, e.g. perpatch.
                               //The host variable, which also holds the host pointer to the data, and

    unsigned int    xstride;
    GhostType       gtype;
    unsigned int    numGhostCells;

    int        varDB_index;     //Where this also shows up in the varDB.  We can use this to
                                //get the rest of the information we need.

    bool            validOnGPU; //true if the GPU copy is the "live" copy and not an old version of the data.
    bool            validOnCPU; //true if the CPU copy is the current "live" copy. (It's possible to be both.)
    bool            queueingOnGPU; //true if we've created the variable but we haven't yet set validOnGPU to true.

  };


  struct charlabelPatchMatl {
    std::string     label;
    int        patchID;
    int        matlIndex;
    charlabelPatchMatl(const char * label, int patchID, int matlIndex) {
      this->label = label;
      this->patchID = patchID;
      this->matlIndex = matlIndex;
    }
    //This so it can be used in an STL map
    bool operator<(const charlabelPatchMatl& right) const {
      if (this->label < right.label) {
        return true;
      } else if (this->label == right.label && (this->patchID < right.patchID)) {
        return true;
      } else if (this->label == right.label && (this->patchID == right.patchID) && (this->matlIndex < right.matlIndex)) {
        return true;
      } else {
        return false;
      }

    }

  };

  struct tempGhostCellInfo {  //We only need enough information to copy a linear chunk of data.
    std::string     label;
    int        patchID;
    int        matlIndex;

    void* cpuDetailedTaskOwner;   //Only the task that placed the item in should be the one that copies it.

    int        xstride;

    void*      device_ptr;
    int        memSize;
    int3 ghostCellLow;
    int3 ghostCellHigh;
    int toPatchID;
    int fromDeviceIndex;
    int toDeviceIndex;
    //int fromresource;
    //int toresource;
    //bool copied;
  };


  //______________________________________________________________________
  // GPU GridVariable methods
  HOST_DEVICE void get(const GPUGridVariableBase& var, char const* label, int patchID, int matlIndex);
  HOST_DEVICE void get(const GPUReductionVariableBase& var, char const* label, int patchID, int matlIndex);
  HOST_DEVICE void get(const GPUPerPatchBase& var, char const* label, int patchID, int matlIndex);

  HOST_DEVICE void getModifiable(GPUGridVariableBase& var, char const* label, int patchID, int matlIndex);
  HOST_DEVICE void getModifiable(GPUReductionVariableBase& var, char const* label, int patchID, int matlIndex);
  HOST_DEVICE void getModifiable(GPUPerPatchBase& var, char const* label, int patchID, int matlIndex);

  //HOST_DEVICE void put(GPUGridVariableBase& var, char const* label, int patchID, int matlIndex, bool overWrite=false);
  HOST_DEVICE void put(GPUReductionVariableBase& var, char const* label, int patchID, int matlIndex, bool overWrite=false);
  HOST_DEVICE void put(GPUGridVariableBase& var, char const* label, int patchID, int matlIndex, size_t xstride, GhostType gtype = None, int numGhostCells = 0, GridVariableBase* gridVar = NULL, void* hostPtr = NULL);
  HOST_DEVICE void put(GPUPerPatchBase& var, char const* label, int patchID, int matlIndex, size_t xstride, GPUPerPatchBase* gridVar = NULL, void* hostPtr = NULL);

  HOST_DEVICE void allocateAndPut(GPUGridVariableBase& var, char const* label, int patchID, int matlID, int3 low, int3 high, size_t sizeOfDataType, GhostType gtype = None, int numGhostCells = 0);
  HOST_DEVICE void allocateAndPut(GPUReductionVariableBase& var, char const* label, int patchID, int matlIndex);
  HOST_DEVICE void allocateAndPut(GPUPerPatchBase& var, char const* label, int patchID, int matlIndex, size_t sizeOfDataType);

  //HOST_DEVICE void* getPointer(char const* label, int patchID, int matlIndex);
  HOST_DEVICE void putContiguous(GPUGridVariableBase &var, char const* indexID, char const* label, int patchID, int matlID, int3 low, int3 high, size_t sizeOfDataType, GridVariableBase* gridVar, bool stageOnHost);
  HOST_DEVICE void allocate(const char* indexID, size_t size);
private:


public:

  //______________________________________________________________________
  // GPU DataWarehouse support methods
  
  HOST_DEVICE bool exist(char const* label, int patchID, int matlID);
  HOST_DEVICE bool exist(char const* label, int patchID, int matlID, int3 host_size, int3 host_offset, bool skipContiguous = false, bool onlyContiguous = false);
  HOST_DEVICE bool existContiguously(char const* label, int patchID, int matlIndex, int3 host_size, int3 host_offset);
  HOST_DEVICE bool remove(char const* label, int patchID, int matlID);
  HOST_DEVICE void init_device(size_t objectSizeInBytes );
  HOST_DEVICE void syncto_device(void *cuda_stream);
  HOST_DEVICE void clear();
  HOST_DEVICE GPUDataWarehouse* getdevice_ptr(){return d_device_copy;};
  HOST_DEVICE void setDebug(bool s){d_debug=s;}
  HOST_DEVICE cudaError_t copyDataHostToDevice(char const* indexID, void *cuda_stream);
  HOST_DEVICE cudaError_t copyDataDeviceToHost(char const* indexID, void *cuda_stream);
  HOST_DEVICE void copyHostContiguousToHost(GPUGridVariableBase& device_var, GridVariableBase* host_var, char const* label, int patchID, int matlIndex);

  //______________________________________________________________________
  //Additional support methods
  HOST_DEVICE void putMaterials(std::vector< std::string > materials);
  HOST_DEVICE materialType getMaterial(int i) const;
  HOST_DEVICE int getNumMaterials() const;
  HOST_DEVICE void putGhostCell(char const* label, int sourcePatchID, int destPatchID, int matlID,
                                int3 sharedLowCoordinates, int3 sharedHighCoordinates, int3 virtualOffset,
                                bool sourceIsInTempGhostCells, void * data_ptr, int3 var_offset, int3 var_size, int xstride);
  HOST_DEVICE bool getValidOnGPU(char const* label, int patchID, int matlID);
  HOST_DEVICE void setValidOnGPU(char const* label, int patchID, int matlID);
  HOST_DEVICE bool getValidOnCPU(char const* label, int patchID, int matlID);
  HOST_DEVICE void setValidOnCPU(char const* label, int patchID, int matlID);

  //Regarding ghost cells...
  //Scenario 1: ghost cell on CPU -> var on same CPU:  Managed by getGridVar().
  //Scenario 2: ghost cell on CPU -> var on different CPU:  Managed by sendMPI()
  //Scenario 3: ghost cell on different CPU -> var on CPU:  Managed by recvMPI() and then getGridVar().

  //Scenario 4: ghost cell on GPU -> var on different host's CPU:  Immediately after a task runs (send_old_data only?),
  // prepareGPUDependencies is invoked, and it discovers all external dependencies.  It then calls
  // the correct GPU DW's prepareGpuGhostCellIntoGpuArray which stages the ghost cell into an array on
  // that GPU.  Then sendMPI() is called, which calls copyTempGhostCellsToHostVar(), which finds the
  // correct array previously made and copies it to a host grid var.  From there the MPI engine
  // can manage it normally.

  //Scenario 5: ghost cell on GPU -> var on different host's GPU:  Immediately after a task runs (send_old_data only?),
  // prepareGPUDependencies is invoked, and it discovers all external dependencies.  It then calls
  // the correct GPU DW's prepareGpuGhostCellIntoGpuArray which stages the ghost cell into an array on
  // that GPU.  Then sendMPI() is called, which calls copyTempGhostCellsToHostVar(), which finds the
  // correct array previously made and copies it to a host grid var.  From there the receiving host
  // moves onto scenario #6.

  //Scenario 6: ghost cell different GPU -> var on current host's GPU: recvMPI() is called and processes
  // data from other nodes.  These are added to the CPU DW. Then it follows the steps listed in the prior
  // scenario. The UnifiedScheduler sees that the destination is
  // valid on the GPU but the ghost cell is not, it is instead on the CPU.  So it puts the CPU ghost cell
  // in the GPU DW.  It also adds the copying information and task ID information to correct GPU DW's d_ghostCellDB.
  // Once d_ghostCellDB is on the GPU then a GPU kernel is called via copyGpuGhostCellsToGpuVarsInvoker
  // which process all entries owned by the same task ID listed earlier.  All data is then copied into the
  // correct destination GPU var.

  //Scenario 7: ghost cell different CPU -> var on current host's GPU: Same as scenario 5.

  //Scenario 8: ghost cell on CPU -> var on current host's GPU: The UnifiedScheduler sees that the destination is
  // valid on the GPU but the ghost cell is not, it is instead on the CPU.  So it puts the CPU ghost cell
  // in the GPU DW.  It also adds the copying information and task ID information to correct GPU DW's d_ghostCellDB.
  // Once d_ghostCellDB is on the GPU then a GPU kernel is called via copyGpuGhostCellsToGpuVarsInvoker
  // which process all entries owned by the same task ID listed earlier.  All data is then copied into the
  // correct destination GPU var.

  //Scenario 9: ghost cell on GPU -> var on same GPU:  initiateH2D recognizes the destination var is valid
  // and so it adds the copying information and task ID information to correct GPU DW's d_ghostCellDB.
  // Once d_ghostCellDB is on the GPU then a GPU kernel is called via copyGpuGhostCellsToGpuVarsInvoker
  // which process all entries owned by the same task ID listed earlier.  All data is then copied into the
  // correct destination GPU var.

  //Scenario 10: var and ghost cell on CPU -> var on GPU: initiateH2D recognizes the data is on the CPU but not on the GPU.
  // It calls the CPU's getGridVar indicating how many ghost cells, resulting in a CPU grid var with ghost cells in it.
  // It is added to the GPU DW by adding it to the host's varPointers collection and d_varDB collection.  The
  // d_varDB is copied to the GPU.

  //Scenario 11: var and ghost cell on GPU -> var on CPU: The schedule recognizes the task is on the CPU, so initiateD2H
  // checks if any data is valid on the GPU (see TODO).  If the GPU var has larger dimensions, then resize the CPU var to match
  // (any task not needing ghost cells will just ignore it anyway).  Copy the data to the CPU.  Then mark CPU data
  // as valid.

  //Scenario 12: ghost cells and vars are on both CPU and a GPU -> var on different GPU: This happens after an data output
  //task runs.  Instead of doing a GPU->different GPU copy, instead treat it like Scenario 9.

  //Scenario 13: Ghost cells on one GPU -> var on different GPU but same host:  Immediately after a task runs,
  // prepareGPUDependencies is invoked, and it discovers all internal dependencies.  It then calls
  // the correct GPU DW's prepareGpuGhostCellIntoGpuArray which stages the ghost cell into an array on
  // that GPU.  postMPISends() is invoked.  Then copyGPUInternalDependencies() is invoked.  It gets a collection
  // of all tempGhostCells items belonging to that task.  For each ghost cell it then calls
  // prepareGpuToGpuGhostCellDestination(), which creates a same sized array on the destination GPU.  This info
  // is added into that GPU DW tempGhostCells.  Then a cudaMemcpyPeer() is called copying data between devices.
  //
  //

  //Useful for sendMPI.  Processes the tempGhostCells collection for the exact var and size that sendMPI requested.
  //It pulls data out of the GPU into a host array variable.  From there it's managed like any other CPU variable.
  HOST_DEVICE void  copyTempGhostCellsToHostVar(void* hostVarPointer, int3 ghostCellLow, int3 ghostCellHigh, char const* label, int patchID, int matlID);

  //Called by prepareGPUDependencies, handles the GPU->anywhere else but the same GPU
  //It does this by staging the ghost cell data into an array for something else to pick up and copy.
  //These arrays are tracked by the tempGhostCells collection.
  HOST_DEVICE void prepareGpuGhostCellIntoGpuArray(void* cpuDetailedTaskOwner,
                                                    int3 ghostCellLow, int3 ghostCellHigh,
                                                    int xstride,
                                                    char const* label, int matlID,
                                                    int fromPatchID, int toPatchID,
                                                    int fromDeviceIndex, int toDeviceIndex,
                                                    int fromresource, int toresource);
  //Called by prepareGpuGhostCellIntoGpuArray, this is a method for a kernel which
  //goes through everything listed in the d_ghostCellDB and copies them to a specified array.
  HOST_DEVICE void copyGhostCellsToArray(void* d_ghostCellData, int index, int3 ghostCellLow, int3 ghostCellHigh);


  //This and the function below go through the d_ghostCellData array and copies data into
  //the correct destination GPU var.  This would be the final step of a GPU ghost cell transfer.
  HOST_DEVICE void copyGpuGhostCellsToGpuVarsInvoker(cudaStream_t* stream);
  HOST_DEVICE void copyGpuGhostCellsToGpuVars();

  //Called by copyGPUGhostCellsBetweenDevices().  If the destination is another GPU on the same physical node
  //then this creates room on the destination GPU and stores that information in d_ghostCellDB to
  //later process.
  HOST_DEVICE void prepareGpuToGpuGhostCellDestination(void* cpuDetailedTaskOwner,
                                                      int3 ghostCellLow, int3 ghostCellHigh,
                                                      int xstride,
                                                      char const* label, int matlID,
                                                      int fromPatchID, int toPatchID,
                                                      int fromDeviceIndex, int toDeviceIndex,
                                                      void * &data_ptr);

  HOST_DEVICE bool ghostCellCopiesNeeded();
  HOST_DEVICE void getSizes(int3& low, int3& high, int3& siz, GhostType& gtype, int& numGhostCells, char const* label, int patchID, int matlIndex);
  HOST_DEVICE void getTempGhostCells(void * dtask, std::vector<tempGhostCellInfo>& temp);

  HOST_DEVICE void* getPlacementNewBuffer();

private:
  HOST_DEVICE void init();
  HOST_DEVICE dataItem* getItem(char const* label, int patchID, int matlID);
  HOST_DEVICE void resetdVarDB();
 // HOST_DEVICE void copyGpuGhostCellsToGpuVars();

private:

  HOST_DEVICE void printGetError(const char* msg, char const* label, int patchID, int matlIndex);

  void * placementNewBuffer;                        //For task DWs, we want to seraliaze and size this object as small as possible.
                                                    //So we create a buffer, and keep track of the start of that buffer here.
  mutable SCIRun::CrowdMonitor allocateLock;
  mutable SCIRun::CrowdMonitor varLock;


  materialItem d_materialDB[MAX_MATERIALS];
  dataItem d_levelDB[MAX_LVITEM];
  //TODO:, remove all references
  //ghostCellItem d_ghostCellDB[MAX_GHOST_CELLS];  //For the device side. This contains information about what ghost cell copies need to complete.
                                                     //This only covers GPU->same GPU scenarios.  It is meant to copy the data directly into the
                                                     //correct GPU var, thus completing the ghost cell merging process.
  int d_numItems;

  int d_numMaterials;
  int numGhostCellCopiesNeeded;
  //int d_numGhostCells;

  GPUDataWarehouse*  d_device_copy;             //The pointer to the copy of this object in the GPU.
  bool d_dirty;                                 //if this changes, we have to recopy the GPUDW.
  int d_device_id;
  bool d_debug;
  size_t objectSizeInBytes;

  //These being here do not pose a problem for the CUDA compiler

  std::map<charlabelPatchMatl, allVarPointersInfo> varPointers; //For the host side.  The variable database.
                                                   //Holds the host ptr, the device ptr, and a staging
                                                   //contiguous host ptr var location.  This is not part
                                                   //of d_varDB because the device doesn't need to know about
                                                   //the host pointers.  Being a map makes this much faster.


  std::map<std::string, contiguousArrayInfo> contiguousArrays;

  std::vector<tempGhostCellInfo> tempGhostCells; //For the host side.  Holds both ghost cell data and information about prepared ghosts cell.
                                                 //Data is loaded in here through prepareGPUDependencies, which checks ahead and
                                                 //loads all ghost cells into contiguous chunks ready for copying to the destination.
                                                 //See prepareGPUDependencies for more info.

  //dataItem *d_varDB;

  //This variable MUST be last.  C and C++ emphasize that
  //dataItem d_varDB[MAX_ITEM];                      //For the device side.  The variable database.
  dataItem d_varDB[MAX_ITEM];                      //A very large buffer.  Important note,
                                                  //we should never transfer a full GPUDataWarehouse object as is.  Instead
                                                  //we should use malloc and only use a section of it.  See here for more
                                                  //information:  http://www.open-std.org/Jtc1/sc22/wg14/www/docs/dr_051.html
                                                  //Doing it this way allows to only need one malloc instead of two.  And the
                                                  //object is serialized allowing for a single copy if needed instead of
                                                  //worrying about two copies (one for the GPU DW and one for the array).
                                                  //As another thread writes to this array, it first accumulates all items that
                                                  //will go into it, then count them up, write that section in one chunk,
                                                  //then that thread can copy only that chunk to the GPUDW.
                                                  //***This must be the last data member of the class***
                                                  //This follows C++ 98 standards "Nonstatic data members of a (non-union) class
                                                  //declared without an intervening access-specifier are allocated so that later
                                                  //members have higher addresses within a class object. The order of allocation of
                                                  //nonstatic data members separated by an access-specifier is unspecified (11.1).
                                                  //Implementation alignment requirements might cause two adjacent members not to
                                                  //be allocated immediately after each other; so might requirements for space for
                                                  //managing virtual functions (10.3) and virtual base classes (10.1)."

};


} //end namespace Uintah

#endif
