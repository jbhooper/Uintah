/****************************************
CLASS
    VectorFieldExtractor

    

OVERVIEW TEXT
    This module receives a DataArchive object.  The user
    interface is dynamically created based information provided by the
    DataArchive.  The user can then select which variables he/she
    wishes to view in a visualization.



KEYWORDS
    ParticleGridReader, Material/Particle Method

AUTHOR
    Kurt Zimmerman
    Department of Computer Science
    University of Utah
    June, 2000

    Copyright (C) 2000 SCI Group

LOG
    Created June 27, 2000
****************************************/
#ifndef VECTORFIELDEXTRACTOR_H
#define VECTORFIELDEXTRACTOR_H 1


#include <Uintah/Datatypes/Archive.h>
#include <Uintah/Datatypes/ArchivePort.h>
#include <PSECore/Datatypes/VectorFieldPort.h>
#include <PSECore/Dataflow/Module.h> 
#include <SCICore/TclInterface/TCLvar.h> 
#include <string>
#include <vector>


namespace Uintah {
namespace Modules {
using namespace Uintah::Datatypes;
using namespace PSECore::Dataflow;
using namespace PSECore::Datatypes;
using namespace SCICore::TclInterface;

class VectorFieldExtractor : public Module { 
  
public: 

  // GROUP: Constructors
  //////////
  VectorFieldExtractor(const clString& id); 

  // GROUP: Destructors
  //////////
  virtual ~VectorFieldExtractor(); 

  // GROUP: cloning and execution 
  ////////// 
  virtual void execute(); 

  //////////
  // overides tcl_command in base class Module
  //  void tcl_command( TCLArgs&, void* );

protected:
  
private:

  TCLstring tcl_status;

  TCLstring sVar;
  TCLint sMatNum;

  const TypeDescription *type;

  ArchiveIPort *in;
  VectorFieldOPort *sfout;
  
  std::string positionName;

  ArchiveHandle  archiveH;
  void setVars();

}; //class 

} // end namespace Modules
} // end namespace Uintah


#endif
