/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2004 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/



/*
 *  CleanupManager.cc: Manage exitAll callbacks.
 *
 *  Written by:
 *   Michael Callahan
 *   Department of Computer Science
 *   University of Utah
 *   November 2004
 *
 *  Copyright (C) 2004 SCI Group
 */

#include <Core/Thread/CleanupManager.h>

#include <sgi_stl_warnings_off.h>
#include <algorithm>
#include <sgi_stl_warnings_on.h>


namespace SCIRun {

std::vector<std::pair<CleanupManagerCallback, void *> > CleanupManager::callbacks_;
Mutex CleanupManager::lock_("CleanupManager lock");

typedef std::pair<CleanupManagerCallback, void *> CMCPair;

void
CleanupManager::add_callback(CleanupManagerCallback cb, void *data)
{
  lock_.lock();
  if (std::find(callbacks_.begin(), callbacks_.end(), CMCPair(cb, data))
      == callbacks_.end())
  {
    callbacks_.push_back(CMCPair(cb, data));
  }
  lock_.unlock();
}

void
CleanupManager::invoke_remove_callback(CleanupManagerCallback callback,
				       void *data)
{
  lock_.lock();
  callback(data);
  callbacks_.erase(std::remove(callbacks_.begin(), callbacks_.end(),
			       CMCPair(callback, data)),
		   callbacks_.end());
  lock_.unlock();
}

void
CleanupManager::remove_callback(CleanupManagerCallback callback, void *data)
{
  lock_.lock();
  callbacks_.erase(std::remove(callbacks_.begin(), callbacks_.end(),
			       CMCPair(callback, data)),
		   callbacks_.end());
  lock_.unlock();
}

void
CleanupManager::call_callbacks()
{
  lock_.lock();
  for (unsigned int i = 0; i < callbacks_.size(); i++)
  {
    callbacks_[i].first(callbacks_[i].second);
  }
  callbacks_.clear();
  lock_.unlock();
}


} // End namespace SCIRun

