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
 *  TCL.cc: Interface to TCL
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   March 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#include <Core/Containers/String.h>
//#include <Dataflow/Network/Module.h>
#include <Core/Malloc/Allocator.h>
#include <Core/GuiInterface/GuiManager.h>
#include <Core/GuiInterface/Remote.h>
#include <Core/GuiInterface/TCL.h>
#include <Core/GuiInterface/TCLTask.h>
#include <Core/GuiInterface/GuiVar.h>
#include <iostream>
using std::cerr;
using std::endl;
using std::ostream;
#include <tcl.h>
#include <tk.h>
#include <stdlib.h>
#include <string.h>
#ifndef _WIN32
#include <unistd.h>                             // defines read() and write()
#endif

extern "C" Tcl_Interp* the_interp;

namespace SCIRun {

struct TCLCommandData {
    TCL* object;
    void* userdata;
};

GuiManager *gm = NULL;

void set_guiManager (GuiManager *mgr)
{
    gm = mgr;
}
 
#if 0
static clString prefix()
{
    static int haveit=0;
    static clString pf;
    if(!haveit){
	char* p = getenv("PSE_WORK");
	if(p){
	    pf=clString(p);
	} else {
		printf("Error: PSE_WORK variable not set!\n");
	    Task::exit_all(-1);
	}
    }
    return pf;
}

static clString application()
{
    static int haveit=0;
    static clString app;
    if(!haveit){
	char* a = getenv("APPLICATION");
	if(a){
	    app=clString(a);
	} else {
		printf("Error; APPLICATION environment variable not set!\n");
	    Task::exit_all(-1);
	}
    }
    return app;
}

#endif

void TCL::execute(const clString& str)
{
#ifndef _WIN32
    if (gm != NULL) {
    	int skt = gm->getConnection();

	printf ("TCL::execute(%s): Got skt from gm->getConnection() = %d", 
	str(), skt);

	// format request - no TCL variable name, just a string to execute
	TCLMessage msg;
	msg.f = exec;
	strcpy (msg.un.tstring, str());

	// send request to server - no need for reply, error goes to Tk
	if (sendRequest (&msg, skt) == -1) {
	    // error case ???
	}
        gm->putConnection (skt);
    } else
#endif
    {
        TCLTask::lock();
        int code = Tcl_Eval(the_interp, const_cast<char *>(str()));
        if(code != TCL_OK)
	    Tk_BackgroundError(the_interp);
        TCLTask::unlock();
    }
}

void TCL::execute(char* str)
{
#ifndef _WIN32
    if (gm != NULL) {
    	int skt = gm->getConnection();
printf ("TCL::execute(%s): Got skt from gm->getConnection() = %d", str, skt);

	// format request - no TCL variable name, just a string to execute
	TCLMessage msg;
	msg.f = exec;
	strcpy (msg.un.tstring, str);

	// send request to server - no need for reply, error goes to Tk
	if (sendRequest (&msg, skt) == -1) {
	    // error case ???
	}
        gm->putConnection (skt);
    }
    else 
#endif
	{
	    //printf("TCL::execute() 1\n");
	    TCLTask::lock();
        int code = Tcl_Eval(the_interp, str);
        if(code != TCL_OK)
		{
			Tk_BackgroundError(the_interp);
			printf("Tcl_Eval(the_inter,%s) failed\n",str);
		}
        TCLTask::unlock();
    }
}

int TCL::eval(const clString& str, clString& result)
{
    TCLTask::lock();
    int code = Tcl_Eval(the_interp, const_cast<char *>(str()));
    if(code != TCL_OK){
	Tk_BackgroundError(the_interp);
	result="";
    } else {
	result=clString(the_interp->result);
    }
    TCLTask::unlock();
    return code == TCL_OK;
}

int TCL::eval(char* str, clString& result)
{
    TCLTask::lock();
    int code = Tcl_Eval(the_interp, str);
    if(code != TCL_OK){
	Tk_BackgroundError(the_interp);
	result="";
    } else {
	result=clString(the_interp->result);
    }
    TCLTask::unlock();
    return code != TCL_OK;
}

void TCL::source_once(const clString& filename)
{
  clString result;
  if(!eval((clString("source ")+filename)(),result)) {
    char* msg=const_cast<char *>
             ((clString("Couldn't source file '")+filename+"'")());
    Tcl_AddErrorInfo(the_interp,msg);
    Tk_BackgroundError(the_interp);
  }

#if 0
    int code;

    TCLTask::lock();

    clString pse_filename(filename);

    char* fn=const_cast<char *>(pse_filename());
    code = Tcl_EvalFile(the_interp, fn);

    if(code != TCL_OK) {
      char msg[256];
      sprintf(msg,"Failed on loading %s err %d\n", fn, code);
      Tcl_AddErrorInfo(the_interp,msg);
      Tk_BackgroundError(the_interp);
    }

    TCLTask::unlock();
#endif
}

static int do_command(ClientData cd, Tcl_Interp*, int argc, char* argv[])
{
    TCLCommandData* td=(TCLCommandData*)cd;
    TCLArgs args(argc, argv);
    td->object->tcl_command(args, td->userdata);
    if(args.have_result_)
	Tcl_SetResult(the_interp, strdup(args.string_()), (Tcl_FreeProc*)free);
    return args.have_error_?TCL_ERROR:TCL_OK;
}

void TCL::add_command(const clString& command, TCL* callback, void* userdata)
{
    TCLTask::lock();
    TCLCommandData* command_data=scinew TCLCommandData;
    command_data->object=callback;
    command_data->userdata=userdata;
    Tcl_CreateCommand(the_interp, const_cast<char *>(command()), do_command, command_data, 0);
    TCLTask::unlock();
}

void
TCL::delete_command( const clString& command )
{
    TCLTask::lock();
    Tcl_DeleteCommand(the_interp, const_cast<char *>(command()) );
    TCLTask::unlock();
}

TCL::TCL()
{
}

TCL::~TCL()
{
}

void TCL::emit_vars(ostream& out, clString& midx)
{
    for(int i=0;i<vars.size();i++)
      {
//	cerr << "emit: " << vars[i]->str() << endl;
        vars[i]->emit(out, midx);
      }
}


TCLArgs::TCLArgs(int argc, char* argv[])
: args_(argc)
{
    for(int i=0;i<argc;i++)
	args_[i] = clString(argv[i]);
    have_error_ = false;
    have_result_ = false;
}

TCLArgs::~TCLArgs()
{
}

int TCLArgs::count()
{
    return args_.size();
}

clString TCLArgs::operator[](int i)
{
    return args_[i];
}

void TCLArgs::error(const clString& e)
{
    string_ = e;
    have_error_ = true;
    have_result_ = true;
}

void TCLArgs::result(const clString& r)
{
    if(!have_error_){
	string_ = r;
	have_result_ = true;
    }
}

void TCLArgs::append_result(const clString& r)
{
    if(!have_error_){
	string_ += r;
	have_result_ = true;
    }
}

void TCLArgs::append_element(const clString& e)
{
    if(!have_error_){
	if(have_result_)
	    string_ += ' ';
	string_ += e;
	have_result_ = true;
    }
}

clString TCLArgs::make_list(const clString& item1, const clString& item2)
{
    char* argv[2];
    argv[0]=const_cast<char *>(item1());
    argv[1]=const_cast<char *>(item2());
    char* ilist=Tcl_Merge(2, argv);
    clString res(ilist);
    free(ilist);
    return res;
}

clString TCLArgs::make_list(const clString& item1, const clString& item2,
			const clString& item3)
{
    char* argv[3];
    argv[0]=const_cast<char *>(item1());
    argv[1]=const_cast<char *>(item2());
    argv[2]=const_cast<char *>(item3());
    char* ilist=Tcl_Merge(3, argv);
    clString res(ilist);
    free(ilist);
    return res;
}

clString TCLArgs::make_list(const Array1<clString>& items)
{
    char** argv=scinew char*[items.size()];
    for(int i=0;i<items.size();i++)
	argv[i]=const_cast<char *>(items[i]());
    char* ilist=Tcl_Merge(items.size(), argv);
    clString res(ilist);
    free(ilist);
    delete[] argv;
    return res;
}

void TCL::reset_vars()
{
    for(int i=0;i<vars.size();i++)
	vars[i]->reset();
}

void TCL::register_var(GuiVar* v)
{
    vars.add(v);
}

void TCL::unregister_var(GuiVar* v)
{
    for(int i=0;i<vars.size();i++){
	if(vars[i]==v){
	    vars.remove(i);
	    return;
	}
    }
}

int TCL::get_gui_stringvar(const clString& base, const clString& name,
			   clString& value)
{
    clString n(base+"-"+name);
    TCLTask::lock();
    char* l=Tcl_GetVar(the_interp, const_cast<char *>(n()), TCL_GLOBAL_ONLY);
    if(!l){
	TCLTask::unlock();
	return 0;
    }
    value=l;
    TCLTask::unlock();
    return 1;
}

int TCL::get_gui_boolvar(const clString& base, const clString& name,
			 int& value)
{
    clString n(base+"-"+name);
    TCLTask::lock();
    char* l=Tcl_GetVar(the_interp, const_cast<char *>(n()), TCL_GLOBAL_ONLY);
    if(!l){
	TCLTask::unlock();
	return 0;
    }
    value=0;
    if(!strcmp(l, "yes")){
	value=1;
    } else if(!strcmp(l, "true")){
	value=1;
    } else if(!strcmp(l, "1")){
	value=1;
    }
    TCLTask::unlock();
    return 1;
}

int TCL::get_gui_doublevar(const clString& base, const clString& name,
			   double& value)
{
    clString n(base+"-"+name);
    TCLTask::lock();
    char* l=Tcl_GetVar(the_interp, const_cast<char *>(n()), TCL_GLOBAL_ONLY);
    if(!l){
	TCLTask::unlock();
	return 0;
    }
    char* end;
    value=strtod(l, &end);
    if(*end != 0){
	// Error reading the double....
	TCLTask::unlock();
	return 0;
    } else {
	TCLTask::unlock();
	return 1;
    }
}

int TCL::get_gui_intvar(const clString& base, const clString& name,
			 int& value)
{
    clString n(base+"-"+name);
    TCLTask::lock();
    char* l=Tcl_GetVar(the_interp, const_cast<char *>(n()), TCL_GLOBAL_ONLY);
    if(!l){
	Tk_BackgroundError(the_interp);
	TCLTask::unlock();
	return 0;
    }
    char* end;
    value=(int)strtol(l, &end, 10);
    if(*end != 0){
	// Error reading the double....
	TCLTask::unlock();
	return 0;
    } else {
	TCLTask::unlock();
	return 1;
    }
}

void TCL::set_guivar(const clString& base, const clString& name,
		     const clString& value)
{
    clString n(base+"-"+name);
    TCLTask::lock();
    Tcl_SetVar(the_interp, const_cast<char *>(n()), const_cast<char *>(value()), TCL_GLOBAL_ONLY);
    TCLTask::unlock();
}

} // End namespace SCIRun

