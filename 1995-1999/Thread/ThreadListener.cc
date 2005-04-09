
#include "ThreadListener.h"
#include "Thread.h"
#include "ThreadEvent.h"
#include <iostream.h>

/*
 * An object may implement this interface and register with
 * <b>Thread::</b><i>addListener</i> in order to receive events
 * regarding thread creation, deletion, and other thread activity.
 * These events are received from the <tt>eventQueue</tt> mailbox.
 * This interface is still under design and not yet completely
 * specified.
 */

void ThreadListener::send_event(Thread*, ThreadEvent) {
    cerr << "ThreadListener::send_event not finished\n";
}

ThreadListener::ThreadListener() {
}

ThreadListener::~ThreadListener() {
}

