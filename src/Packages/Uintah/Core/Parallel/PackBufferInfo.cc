
#include <Packages/Uintah/Core/Parallel/PackBufferInfo.h>
#include <Packages/Uintah/Core/Parallel/Parallel.h>
#include <Packages/Uintah/Core/Util/RefCounted.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Util/Assert.h>

using namespace Uintah;

#include <iostream>

PackBufferInfo::PackBufferInfo()
  : BufferInfo()
{
  packedBuffer=0;
}

PackBufferInfo::~PackBufferInfo()
{
  if (packedBuffer && packedBuffer->removeReference())
    delete packedBuffer;
}

void
PackBufferInfo::get_type(void*& out_buf, int& out_count,
			 MPI_Datatype& out_datatype, MPI_Comm comm)
{
  MALLOC_TRACE_TAG_SCOPE("PackBufferInfo::get_type");
  ASSERT(count() > 0);
  if(!have_datatype){
    int packed_size;
    int total_packed_size=0;
    for (int i = 0; i < (int)startbufs.size(); i++) {
      MPI_Pack_size(counts[i], datatypes[i], comm, &packed_size);
      total_packed_size += packed_size;
    }
    packedBuffer = scinew PackedBuffer(total_packed_size);
    packedBuffer->addReference();

    datatype = MPI_PACKED;
    cnt=total_packed_size;
    buf = packedBuffer->getBuffer();
    have_datatype=true;
  }
  out_buf=buf;
  out_count=cnt;
  out_datatype=datatype;
}

void PackBufferInfo::get_type(void*&, int&, MPI_Datatype&)
{
  // Should use other overload for a PackBufferInfo
  SCI_THROW(SCIRun::InternalError("get_type(void*&, int&, MPI_Datatype&) should not be called on PackBufferInfo objects", __FILE__, __LINE__));
}


void
PackBufferInfo::pack(MPI_Comm comm, int& out_count)
{
  MALLOC_TRACE_TAG_SCOPE("PackBufferInfo::pack");
  ASSERT(have_datatype);

  int position = 0;
  int bufsize = packedBuffer->getBufSize();
  for (int i = 0; i < (int)startbufs.size(); i++) {
    if(counts[i]>0)
      MPI_Pack(startbufs[i], counts[i], datatypes[i], buf, bufsize,
	       &position, comm);
  }
  out_count = position;

  // When it is all packed, only the buffer necessarily needs to be kept
  // around until after it is sent.
  delete sendlist;
  sendlist = 0;
  addSendlist(packedBuffer);
}

void
PackBufferInfo::unpack(MPI_Comm comm)
{
  MALLOC_TRACE_TAG_SCOPE("PackBufferInfo::unpack");
  ASSERT(have_datatype);

  int position = 0;
  int bufsize = packedBuffer->getBufSize();
  for (int i = 0; i < (int)startbufs.size(); i++) {
    if(counts[i]>0)
      MPI_Unpack(buf, bufsize, &position, startbufs[i], counts[i], datatypes[i], comm);
  }
}

