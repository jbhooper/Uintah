// PackageDB.cc - Interface to module-finding and loading mechanisms
// $Id$


#include <SCICore/Util/soloader.h>
#ifdef ASSERT
#undef ASSERT
#endif
#include <SCICore/Containers/AVLTree.h>
#include <SCICore/Containers/String.h>
#include <PSECore/Dataflow/PackageDB.h>
#include <iostream.h>

namespace PSECore {
namespace Dataflow {

using namespace SCICore::Containers;
using SCICore::Containers::AVLTree;
using SCICore::Containers::AVLTreeIter;

typedef struct ModuleInfo_tag {
  ModuleMaker maker;
  clString uiFile;
} ModuleInfo;

typedef AVLTree<clString,ModuleInfo*> Category;
typedef AVLTree<clString,Category*> Package;
typedef AVLTree<clString,Package*> Packages;

typedef AVLTreeIter<clString,ModuleInfo*> CategoryIter;
typedef AVLTreeIter<clString,Category*> PackageIter;
typedef AVLTreeIter<clString,Package*> PackagesIter;

PackageDB packageDB;

PackageDB::PackageDB(void) :
  d_db((void*)new Packages), d_packageList(0)
{
}

PackageDB::~PackageDB(void)
  { 
    delete (Packages*)d_db; 
  }

typedef void (*pkgInitter)(const clString& tclPath);

void PackageDB::loadPackage(const clString& packPath) {

  // The format of a package path element is either "soLib" or "soLib(tclPath)"
  // where soLib is the name of a package .so file, and tclPath is the path
  // to use to locate that package's tcl files.  If the first form is used,
  // the tclPath is constructed by substituting "TCL" for everything in the
  // soLib string after the final "/".
  //
  // A package path is a colon-separated list of package path elements.

  clString packagePath(packPath);     // Copy to deal with non-const methods

  while(packagePath!="") {

    // Strip off the first element, leave the rest in the path for the next
    // iteration.

    clString packageElt;
    int firstColon=packagePath.index(':');
    if(firstColon!=-1) {
      packageElt=packagePath.substr(0,firstColon);
cerr << "Before '" << packagePath << "'\n";
      packagePath=packagePath.substr(firstColon+1,-1);
cerr << "After '" << packagePath << "'\n";
    } else {
      packageElt=packagePath;
      packagePath="";
    }

    // Parse the element apart into soName and tclPath, using the default
    // tclpath (soName's directory + "/GUI") if there isn't one specified.

    clString soName;
    clString tclPath;

    int openParen=packageElt.index('(');
    if(openParen!=-1) {
      int closeParen=packageElt.index(')');
      soName=packageElt.substr(0,openParen);
      tclPath=packageElt.substr(openParen+1,closeParen-openParen-1);
    } else {
      soName=packageElt;
      if(pathname(packageElt)!="")
        tclPath=pathname(packageElt)+"/GUI";
      else
        tclPath=pathname(packageElt)+"GUI";
    }

    // Load the package
    {
      clString result;
      TCL::eval(clString(".top.errorFrame.text insert end \"Loading package '")
                +soName+"' with TCLPath '"+tclPath+"'\\n\"",result);

      // Tell Tcl where to find the .tcl files...
      TCL::eval(clString("lappend auto_path ") + tclPath, result);
    }

    //void* so=dlopen(soName(),RTLD_NOW);
	LIBRARY_HANDLE so = GetLibraryHandle(soName());
    if(!so) {
	  //cerr << dlerror() << '\n';
      cerr << "ERROR: Can't open package '" << soName << "'\n";
      continue;
    }
    //pkgInitter initFn=(pkgInitter)dlsym(so,"initPackage");
	pkgInitter initFn=(pkgInitter)GetHandleSymbolAddress(so,"initPackage");
    if(!initFn) {
      cerr << "ERROR: Package '" << soName << "' has no initPackage(...)\n";
      continue;
    }
    initFn(tclPath);

    // You can't close the sofile; it loads more stuff in when you instantiate
    // a module -- all the linking does not occur here.

    // XXX: We need to keep the handle around to avoid opening the same one
    //      a zillion times, and to close them off when you're in a development
    //      cycle so the old inodes can get freed.

  }

  clString result;
  TCL::eval("createCategoryMenu",result);
}

void PackageDB::registerModule(const clString& packageName,
                               const clString& categoryName,
                               const clString& moduleName,
                               ModuleMaker moduleMaker,
                               const clString& tclUIFile) {
  Packages* db=(Packages*)d_db;

  Package* package;
  if(!db->lookup(packageName,package))
    {
      db->insert(packageName,package=new Package);
      d_packageList.add( packageName );
    }

  Category* category;
  if(!package->lookup(categoryName,category))
    package->insert(categoryName,category=new Category);

  ModuleInfo* moduleInfo;
  if(!category->lookup(moduleName,moduleInfo)) {
    moduleInfo=new ModuleInfo;
    category->insert(moduleName,moduleInfo);
  } else cerr << "WARNING: Overriding multiply registered module "
              << packageName << "." << categoryName << "."
              << moduleName << "\n";

  moduleInfo->maker=moduleMaker;
  moduleInfo->uiFile=tclUIFile;
}

Module* PackageDB::instantiateModule(const clString& packageName,
                                     const clString& categoryName,
                                     const clString& moduleName,
                                     const clString& instanceName) const {
  Packages* db=(Packages*)d_db;

  Package* package;
  if(!db->lookup(packageName,package)) {
    cerr << "ERROR: Instantiating from nonexistant package " << packageName 
         << "\n";
    return 0;
  }

  Category* category;
  if(!package->lookup(categoryName,category)) {
    cerr << "ERROR: Instantiating from nonexistant category " << packageName
         << "." << categoryName << "\n";
    return 0;
  }

  ModuleInfo* moduleInfo;
  if(!category->lookup(moduleName,moduleInfo)) {
    cerr << "ERROR: Instantiating nonexistant module " << packageName 
         << "." << categoryName << "." << moduleName << "\n";
    return 0;
  }

#if 0
  // This was McQ's somewhat silly replacement for TCL's tclIndex/auto_path
  // mechanism.  The idea was that there would be a path in the index.cc
  // that pointed to a TCL file to source before instantiating a module
  // of some particular class for the frist time -- sortof a TCL-end class
  // constructor for the module's class.
  //
  // Steve understandably doesn't like new, fragile mechanisms where
  // perfectly good old, traditional ones already exist, so he if0'd this
  // away and added the "lappend auto_path" at package-load-time, above.
  //
  // This code is still here 'cause Some Day it might be nice to allow the
  // source of the TCL files to be stored in the .so (as strings) and eval'd
  // here.  This was the "faraway vision" that drove me to do things this way
  // in the first place, but since that vision seems to have stalled
  // indefinately in lieu of Useful Work, there's no reason not to use
  // auto_path (except that it produces yet one more file to maintain).  And
  // auto_path is useful if you write global f'ns and want to use them in lots
  // of your modules -- auto_path nicely handles this whereas the code below
  // doesn't handle it at all.
  //
  // Some day it might be nice to actually achieve the "package is one .so
  // and that's all" vision, but not today.  :)
  //
  //                                                      -mcq 99/10/6

  if(moduleInfo->uiFile!="") {
    clString result;
    if(!TCL::eval(clString("source ")+moduleInfo->uiFile,result)) {
      cerr << "Can't source UI file " << moduleInfo->uiFile << "...\n";
      cerr << "  TCL Error: " << result << "\n";
    }
    moduleInfo->uiFile="";                       // Don't do it again
  }
#endif

  Module *module = (moduleInfo->maker)(instanceName);
  module->packageName = packageName;
  module->moduleName = moduleName;
  module->categoryName = categoryName;

  return module;
}

Array1<clString> PackageDB::packageNames(void) const {

  // d_packageList is used to keep a list of the packages 
  // that are in this PSE IN THE ORDER THAT THEY ARE SPECIFIED
  // by the user in the Makefile (for main.cc) or in their
  // environment.

  return d_packageList;
}

Array1<clString>
PackageDB::categoryNames(const clString& packageName) const {
  Packages* db=(Packages*)d_db;

  {
    PackagesIter iter(db);
    for(iter.first();iter.ok();++iter) if(iter.get_key()==packageName) {
      Package* package=iter.get_data();
      Array1<clString> result(package->size());
      {
        PackageIter iter(package);
        int i=0;
        for(iter.first();iter.ok();++iter) result[i++]=iter.get_key();
      }
      return result;
    }
  }
  cerr << "WARNING: Unknown package " << packageName << "\n";

  Array1<clString> result(0);
  return result;
}

Array1<clString>
PackageDB::moduleNames(const clString& packageName,
		       const clString& categoryName) const {
  Packages* db=(Packages*)d_db;
  {
    PackagesIter iter(db);
    for(iter.first();iter.ok();++iter) if(iter.get_key()==packageName) {
      Package* package=iter.get_data();
      {
        PackageIter iter(package);
        for(iter.first();iter.ok();++iter) if(iter.get_key()==categoryName) {
          Category* category=iter.get_data();
          Array1<clString> result(category->size());
          {
            CategoryIter iter(category);
            int i=0;
            for(iter.first();iter.ok();++iter) result[i++]=iter.get_key();
          }
          return result;
        }
        cerr << "WARNING: Unknown category " << packageName << "."
             << categoryName << "\n";
      }
    }
  }
  cerr << "WARNING: Unknown package " << packageName << "\n";

  Array1<clString> result(0);
  return result;
}

} // Dataflow namespace
} // PSECore namespace

//
// $Log$
// Revision 1.13  1999/10/06 20:37:37  mcq
// Added memoirs.
//
// Revision 1.12  1999/09/22 23:51:46  dav
// removed debug print
//
// Revision 1.11  1999/09/22 22:39:50  dav
// updated to use tclIndex files
//
// Revision 1.10  1999/09/08 02:26:41  sparker
// Various #include cleanups
//
// Revision 1.9  1999/09/04 06:01:41  sparker
// Updates to .h files, to minimize #includes
// removed .icc files (yeah!)
//
// Revision 1.8  1999/08/31 23:27:53  sparker
// Added Log and Id entries
//
//
