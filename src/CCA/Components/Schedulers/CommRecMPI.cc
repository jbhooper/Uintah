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

#include <CCA/Components/Schedulers/CommRecMPI.h>

#include <Core/Thread/Mutex.h>
#include <Core/Util/DebugStream.h>
#include <Core/Util/FancyAssert.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/Thread/Time.h>

using namespace std;
using namespace Uintah;
using namespace SCIRun;

// sync cerr so it's readable when output by multiple threads
extern SCIRun::Mutex       cerrLock;
extern SCIRun::DebugStream mixedDebug;

static DebugStream dbg( "RecvTiming", false );

double CommRecMPI::WaitTimePerMessage = 0;


//______________________________________________________________________
//
//______________________________________________________________________

void
CommRecMPI::add(       MPI_Request                 id, 
                       int                         bytes, 
                       AfterCommunicationHandler * handler, 
                 const string                    & var, 
                       int                         message )
{
  ids_.push_back( id );
  handlers_.push_back( handler );
  byteCounts_.push_back( bytes );
  vars_.push_back( var );
  messageNums_.push_back( message );
  totalBytes_ += bytes;
}

//______________________________________________________________________
//
//______________________________________________________________________

void
CommRecMPI::print( const ProcessorGroup * pg )
{
  for( unsigned int i = 0; i < ids_.size(); i++ ) {
    cout << pg->myrank() << " Message: " << byteCounts_[i] << " vars: " << " num " << messageNums_[i] << " Vars: " << vars_[i] << "\n";
  }
}

//______________________________________________________________________
//
//______________________________________________________________________

bool
CommRecMPI::waitsome( const ProcessorGroup * pg )
{
  if( ids_.size() == 0 ) {
    return false; // no more to test
  }

  statii.resize( ids_.size() );
  indices.resize( ids_.size() );

  // debugging
  if (mixedDebug.active()) {
    int me = pg->myrank();
    mixedDebug << me << " Waitsome: " << ids_.size() << " waiters:\n";

    for (unsigned i = 0; i < messageNums_.size(); i++) {
      mixedDebug << me << "  Num: " << messageNums_[i] << " size: " << byteCounts_[i] << "\n";
    }
  }

  int     donecount;
  clock_t start = clock();
  
  MPI_Waitsome( (int)ids_.size(), &ids_[0], &donecount, &indices[0], &statii[0] );
  
  WaitTimePerMessage = (clock() - start) / (double)CLOCKS_PER_SEC / donecount;

  return donesome( pg, donecount, statii );
}

//______________________________________________________________________
//
//______________________________________________________________________

bool
CommRecMPI::testsome( const ProcessorGroup * pg )
{
  if (ids_.size() == 0) {
    return false;  // no more to test
  }
  statii.resize(  ids_.size() );
  indices.resize( ids_.size() );
  int me = pg->myrank();

  // debugging
  if (mixedDebug.active()) {
    cerrLock.lock();
    mixedDebug << me << " Calling testsome with " << ids_.size() << " waiters\n";
    cerrLock.unlock();
  }

  int     donecount;
  clock_t start = clock();
  
  MPI_Testsome( (int)ids_.size(), &ids_[0], &donecount, &indices[0], &statii[0] );
  
  if( donecount>0 ){
    WaitTimePerMessage = (clock() - start) / (double)CLOCKS_PER_SEC / donecount;
  }
  return donesome( pg, donecount,statii );
}

//______________________________________________________________________
//
//______________________________________________________________________

bool
CommRecMPI::donesome( const ProcessorGroup     * pg, 
                            int                  donecount, 
                            vector<MPI_Status> & statii )
{
  bool anyFinished = false;
  int  numReceived = 0;
  int  volReceived = 0;

  ASSERT(donecount != MPI_UNDEFINED);

  for( int i = 0; i < donecount; i++ ) {
    int idx = indices[i];

    if( handlers_[idx] ) {

      if (mixedDebug.active()) {
        cerrLock.lock();
        mixedDebug << "Actually received " << idx << "\n";
        cerrLock.unlock();
      }

      handlers_[idx]->finishedCommunication( pg, statii[i] );
      ASSERT( handlers_[idx] != 0 );

      delete handlers_[idx];
      handlers_[idx] = 0;
    }

    numReceived++;
    volReceived += byteCounts_[idx];
    ids_[idx] = MPI_REQUEST_NULL;
    totalBytes_ -= byteCounts_[idx];
    byteCounts_[idx] = 0;
  }

  if (dbg.active() && numReceived > 0) {
    if (pg->myrank() == pg->size() / 2) {
      cerrLock.lock();
      dbg << pg->myrank() << " Time: " << Time::currentSeconds() << " , NumReceived= " << numReceived << " , VolReceived: "
          << volReceived << endl;
      cerrLock.unlock();
    }
  }

  if (donecount == (int)ids_.size()) {
    ASSERT(totalBytes_ == 0);
    ids_.clear();
    handlers_.clear();
    byteCounts_.clear();
    vars_.clear();
    messageNums_.clear();
    return false;  // no more to test
  }

  // remove finished requests
  int j = 0;
  for (int i = 0; i < (int)ids_.size(); i++) {
    if (ids_[i] != MPI_REQUEST_NULL) {
      ids_[j] = ids_[i];
      handlers_[j] = handlers_[i];
      byteCounts_[j] = byteCounts_[i];
      messageNums_[j] = messageNums_[i];
      vars_[j] = vars_[i];
      ++j;
    }
  }
  ASSERT( (int )ids_.size() - donecount == j );

  ids_.resize(j);
  handlers_.resize(j);
  vars_.resize(j);
  messageNums_.resize(j);
  byteCounts_.resize(j);

  return !anyFinished;  // keep waiting until something finished
}

//______________________________________________________________________
//
//______________________________________________________________________

void
CommRecMPI::waitall( const ProcessorGroup * pg )
{
  if( ids_.size() == 0 ) {
    return;
  }

  statii.resize( ids_.size() );
//    mixedDebug << me << " Calling waitall with " << ids.size() << " waiters\n";
  clock_t start = clock();

  MPI_Waitall((int)ids_.size(), &ids_[0], &statii[0]);

  WaitTimePerMessage = (clock() - start) / (double)CLOCKS_PER_SEC / ids_.size();
  //  mixedDebug << me << " Done calling waitall with " << ids_.size() << " waiters\n";

  for (int i = 0; i < (int)ids_.size(); i++) {
    if( handlers_[i] ) {
      handlers_[i]->finishedCommunication( pg, statii[i] );
      ASSERT( handlers_[i] != 0 );
      delete handlers_[i];
      handlers_[i] = 0;
    }
  }

  ids_.clear();
  handlers_.clear();
  byteCounts_.clear();
  messageNums_.clear();
  vars_.clear();
  totalBytes_ = 0;
}
