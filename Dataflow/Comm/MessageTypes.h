
/*
 *  MessageTypes.h: enum definitions for Message Types...
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   March 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#ifndef SCI_project_MessageTypes_h
#define SCI_project_MessageTypes_h 1

#include <Dataflow/share/share.h>

namespace SCIRun {

class PSECORESHARE MessageTypes {
    MessageTypes();
public:
    void dummy(); // Keeps g++ happy...
    enum MessageType {
	DoCallback,
	DoDBCallback,
	ExecuteModule,
	TriggerPort,
	ReSchedule,
	MultiSend,
	Demand,

	MUIDispatch,

	GeometryInit,
	GeometryAddObj,
	GeometryDelObj,
	GeometryDelAll,
	GeometryFlush,
	GeometryFlushViews,
	GeometryGetData,
	GeometryGetNViewWindows,
	GeometrySetView,

	GeometryPick,
	GeometryRelease,
	GeometryMoved,

	ViewWindowRedraw,
	ViewWindowMouse,
	ViewWindowDumpImage,
	ViewWindowDumpObjects,

	AttachDialbox,
	DialMoved,

	TrackerMoved,

	ModuleGeneric1,
	ModuleGeneric2,
	ModuleGeneric3,
	ModuleGeneric4
	
    };
};

} // End namespace SCIRun


#endif
