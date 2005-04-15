/*
 * This file was automatically generated by SCC - do NOT edit!
 * You should edit RealtimeThread.scc instead 
 */

#ifndef SCI_THREAD_REALTIMETHREAD_H
#define SCI_THREAD_REALTIMETHREAD_H 1

/*
 * Extended threads with realtime capabilities.  These realtime
 * capabilities may not be available on all machines, but the
 * functionality will be approximated on those that do not.
 *
 * <p>On multiprocessor machines, the realtime threads typically
 * take control of a single CPU.  On uniprocessor machines, they
 * can consume a large percentage of the processing resources
 * if the interval is too small.
 *
 * <p>Until this thread calls <i>frameSchedule</i>, this thread
 * behaves exactly like any other thread.
 */
class Thread;

/**************************************
 
CLASS
   RealtimeThread
   
KEYWORDS
   RealtimeThread
   
DESCRIPTION
   Extended threads with realtime capabilities.  These realtime
   capabilities may not be available on all machines, but the
   functionality will be approximated on those that do not.
  
   <p>On multiprocessor machines, the realtime threads typically
   take control of a single CPU.  On uniprocessor machines, they
   can consume a large percentage of the processing resources
   if the interval is too small.
  
   <p>Until this thread calls <i>frameSchedule</i>, this thread
   behaves exactly like any other thread.
PATTERNS


WARNING
   
****************************************/
class RealtimeThread {
    int d_interval;
public:
    //////////
    //Arrange to have the thread frame scheduled, so that it can execute every <i>interval</i>
    //microseconds.  An interval ranging from 1 to 100 will schedule the thread at that multiple
    //of the vertical retrace interval.  An interval of zero will cause frame scheduling to be
    //disabled.  Using this may cause the thread to seize control of one CPU, and should not 
    //be used for low frequency events.  <i>idleThread</i> will be executed for the raemainder of
    //the interval once the thread has yielded.  Both threads should  call <i>frameReady</i> when
    //they are ready to begin realtime execution, and the main thread should call frameYield when
    //the frame is completed.
    //
    //<p><b>This interface is sketchy and is likely to evolve, or it may be removed and we will
    //leave the complicated OS stuff in the graphical viewer unless somebody comes up with another
    //application that needs it.</b>
    void frameSchedule(int interval, RealtimeThread* idleThread);
    
    //////////
    //Called when Realtime threads are ready to begin frame scheduling.  The first frame begins
    //immediately after this returns.
    void frameReady();
    
    //////////
    //Called by the main thread (not the idle thread) when it has completed the work for this
    //frame.
    void frameYield();

    //////////
    //Flag which will get set if a frame scheduled thread exceeds the allocated time quantum.  
    //The thread should use the high resolution clock interface to determine the amount of overflow.
    bool quantumExceeded;

    //////////
    //Return the frame interval in microseconds.  For frames which are keyed to the vertical retrace
    //interval, this may be different than the interval supplied in the constructor.
    int frameInterval();
};

#endif


