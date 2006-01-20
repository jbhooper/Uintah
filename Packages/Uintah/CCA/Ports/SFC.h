#ifndef _SFC
#define _SFC


#include <vector>
#include <queue>
#include <iostream>
#include <cmath>
using namespace std;

#include<Packages/Uintah/Core/Parallel/ProcessorGroup.h>
#include<Core/Thread/Time.h>

namespace Uintah{

#ifdef _TIMESFC_
SCIRun::Time *timer;
double start, finish;
double ptime, cleantime, sertime, gtime;
#endif

enum Curve {HILBERT, MORTON, GREY};
enum CleanupType{BATCHERS,LINEAR};

template<class BITS>
struct History
{
	unsigned int i;
	BITS bits;
};

extern int dir3[8][3];
extern int dir2[4][2];

extern int hinv3[][8];
extern int ginv3[][8];
extern int hinv2[][4];
extern int ginv2[][4];

extern int horder3[][8];
extern int gorder3[][8];
extern int morder3[][8];
extern int horder2[][4];
extern int gorder2[][4];
extern int morder2[][4];

extern int horient3[][8];
extern int gorient3[][8];
extern int morient3[][8];
extern int horient2[][4];
extern int gorient2[][4];
extern int morient2[][4];



#define BINS (1<<DIM)
template<int DIM, class LOCS>
class SFC
{
public:
	SFC(int dir[][DIM], ProcessorGroup *d_myworld) : dir(dir),set(0), locsv(0), locs(0), orders(0), sendbuf(0), recievebuf(0), mergebuf(0), d_myworld(d_myworld), block_size(3000), blocks_in_transit(3), sample_percent(.1), cleanup(BATCHERS) {};
	void GenerateCurve(bool force_serial=false);
	void SetRefinements(unsigned int refinements);
	void SetLocalSize(unsigned int n);
	void SetLocations(vector<LOCS> *locs);
	void SetOutputVector(vector<unsigned int> *orders);
	void SetMergeParameters(unsigned int block_size, unsigned int blocks_in_transit, float sample_percent);
	void SetBlockSize(unsigned int b) {blocksize=b;};
	void SetBlocksInTransit(unsigned int b) {blocks_in_transit=b;};
	void SetSamplePercent(float p) {sample_percent=p;};
	void SetCleanup(CleanupType cleanup) {this->cleanup=cleanup;};
	void Balance();
	void Profile();

protected:
	
	//order and orientation arrays
	int (*order)[BINS];
	int (*orientation)[BINS];
	int (*inverse)[BINS];
	
	//direction array
	int (*dir)[DIM];

	//curve parameters
	LOCS dimensions[DIM];
	LOCS center[DIM];
	unsigned int refinements;
	unsigned int n;
	
	//byte variable used to determine what curve parameters are set
	unsigned char set;

	//XY(Z) locations of points
	vector<LOCS> *locsv;
	LOCS *locs;

	//Output vector
	vector<unsigned int> *orders;
	
	//Buffers
	void* sendbuf;
	void* recievebuf;
	void* mergebuf;
	
	ProcessorGroup *d_myworld;
	
	//Merge-Exchange Parameters
	unsigned int block_size;
	unsigned int blocks_in_transit;
	float sample_percent;
	
	unsigned int *n_per_proc;
	CleanupType cleanup;
	
	int rank, P;
	MPI_Comm Comm;	
	void Serial();
	void SerialR(unsigned int* orders,vector<unsigned int> *bin, unsigned int n, LOCS *center, LOCS *dimension, unsigned int o=0);

	template<class BITS> void SerialH();
	template<class BITS> void SerialHR(History<BITS>* orders,vector<unsigned int > *bin, unsigned int n, LOCS *center, LOCS *dimension, unsigned int o=0, unsigned int r=0, BITS history=0);
	template<class BITS> void Parallel();
	template<class BITS> int MergeExchange(int to);	
	template<class BITS> void PrimaryMerge();
	template<class BITS> void Cleanup();
	template<class BITS> void Batchers();
	template<class BITS> void Linear();

	inline virtual unsigned char Bin(LOCS *point, LOCS *center)=0;
};

template<class LOCS>
class SFC2D : public SFC<2,LOCS>
{

	public:
		SFC2D(Curve curve,ProcessorGroup *d_myworld) : SFC<2,LOCS>(dir2,d_myworld) {SetCurve(curve);};
		void SetCurve(Curve curve);
		void SetDimensions(LOCS wx, LOCS wy);
		void SetCenter(LOCS x, LOCS y);
		void SetRefinementsByDelta(LOCS deltax, LOCS deltay);

	private:
		 inline unsigned char Bin(LOCS *point, LOCS *center);

};

template<class LOCS>
class SFC3D : public SFC<3,LOCS>
{
	public:
		SFC3D(Curve curve,ProcessorGroup *d_myworld) : SFC<3,LOCS>(dir3, d_myworld) {SetCurve(curve);};
		void SetCurve(Curve curve);
		void SetDimensions(LOCS wx, LOCS wy, LOCS wz);
		void SetCenter(LOCS x, LOCS y, LOCS z);
		void SetRefinementsByDelta(LOCS deltax, LOCS deltay, LOCS deltaz);

private:
		inline unsigned char Bin(LOCS *point, LOCS *center);
};

class SFC2f : public SFC2D<float>
{
	public:
		SFC2f(Curve curve,ProcessorGroup *d_myworld) : SFC2D<float>(curve,d_myworld) {};
};

class SFC2d : public SFC2D<double>
{
	public:
		SFC2d(Curve curve,ProcessorGroup *d_myworld) : SFC2D<double>(curve,d_myworld) {} ;
};

class SFC3f : public SFC3D<float>
{
	public:
		SFC3f(Curve curve,ProcessorGroup *d_myworld) : SFC3D<float>(curve,d_myworld) {};
};

class SFC3d : public SFC3D<double>
{
	public:
		SFC3d(Curve curve,ProcessorGroup *d_myworld) : SFC3D<double>(curve,d_myworld) {} ;
};


/***********SFC**************************/
const char errormsg[6][30]={
	"Dimensions not set\n",
	"Center not set\n", 
	"Locations vector not set\n",
	"Output vector not set\n", 
	"Local size not set\n",
	"Refinements not set\n",
					 };

template<int DIM, class LOCS> 
void SFC<DIM,LOCS>::Balance()
{
	//for each edge of the hypercube
	
		//exchange size
		//if balance needed
			//exchange locs
			//set n
}


template<int DIM, class LOCS> 
void SFC<DIM,LOCS>::Profile()
{
	int  starti;
	rank=d_myworld->myrank(); P=d_myworld->size(), Comm=d_myworld->getComm();
	vector<unsigned int> n_per_proc(P);
	this->n_per_proc=&n_per_proc[0];	
	n=2000;
	

	//initialize list
	vector<History<unsigned int> > sbuf, rbuf, mbuf;
	for(int i=0;i<P;i++)
		n_per_proc[i]=n;

	sbuf.resize(n);
	rbuf.resize(n);
	mbuf.resize(n);

	sendbuf=(void*)&(sbuf[0]);
	recievebuf=(void*)&(rbuf[0]);
	mergebuf=(void*)&(mbuf[0]);
	
	starti=rank;
	for(unsigned i=0;i<n;i++)
	{
		mbuf[i].bits=sbuf[i].bits=starti;
		starti+=P;
		//starti+=rand()%(2*P)+1;
	}

	for(int b=25;b<2000;b+=5)
	{
		block_size=b;
		double ptime=0, cleantime=0;
		int r;
		MPI_Barrier(Comm);
		start=timer->currentSeconds();
		for(r=0;r<250;r++)
		{
			
			PrimaryMerge<unsigned int>();
			Cleanup<unsigned int>();
			swap(sendbuf,mergebuf);
		
		}
		finish=timer->currentSeconds();
		ptime=finish-start;	
		double sum,sum2;
		MPI_Reduce(&ptime,&sum,1,MPI_DOUBLE,MPI_SUM,0,Comm);
		MPI_Reduce(&cleantime,&sum2,1,MPI_DOUBLE,MPI_SUM,0,Comm);
				
		if(rank==0)
			cout << b << " " << ptime/P/r << " " << cleantime/P/r << endl;

	}

}

template<int DIM, class LOCS>
void SFC<DIM,LOCS>::GenerateCurve(bool force_serial)
{
#ifdef _TIMESFC_
	cleantime=sertime=ptime=gtime=0;
#endif
	int errors=0;
	unsigned char mask;
	P=d_myworld->size();
	if(P==1 || force_serial)
	{
		errors=5;
		mask=0x0f;
	}
	else 
	{
		errors=6;
		mask=0x1f;
	}	
		
	char res=mask&set;
	if(res!=mask)
	{	
		cout << "Error(s) forming SFC:\n************************\n";
		mask=1;
		for(int i=0;i<errors;i++)
		{
			res=mask&set;
			if(res!=mask)
			{
				cout << "  " << errormsg[i];
			}
			mask=mask<<1;
		}
		cout << "************************\n";	
		return;
	}
	if(P==1 || force_serial)
	{
		Serial();
	}
	else
	{
		//make new sub group if needed?
		rank=d_myworld->myrank();
		Comm=d_myworld->getComm();
		//Pick which generate to use
		if(refinements*DIM<=32)
		{
			Parallel<unsigned int>();
		}
		else if(refinements*DIM<=64)
		{
			Parallel<unsigned long long>();
		}
		else
		{
			cout << "SFC Error not enough bits for histories.... need to design a larger history class or lower refinements\n";
		}
	}
}


template<int DIM, class LOCS>
void SFC<DIM,LOCS>::Serial()
{
	orders->resize(n);

	unsigned int *o=&(*orders)[0];
	
	for(unsigned int i=0;i<n;i++)
	{
		o[i]=i;
	}

	vector<unsigned int> bin[BINS];
	for(int b=0;b<BINS;b++)
	{
		bin[b].reserve(n/BINS);
	}

	//Recursive call
	SerialR(o,bin,n,center,dimensions);

}


template<int DIM, class LOCS> template<class BITS>
void SFC<DIM,LOCS>::SerialH()
{
	orders->resize(n);

	History<BITS> *sbuf=(History<BITS>*)sendbuf;
	
	for(unsigned int i=0;i<n;i++)
	{
		sbuf[i].i=i;
	}

	vector<unsigned int> bin[BINS];
	for(int b=0;b<BINS;b++)
	{
		bin[b].reserve(n/BINS);
	}
	
	//Recursive call
	SerialHR<BITS>(sbuf,bin,n,center,dimensions);

}

template<int DIM, class LOCS> 
void SFC<DIM,LOCS>::SerialR(unsigned int* orders,vector<unsigned int> *bin, unsigned int n, LOCS *center, LOCS *dimension, unsigned int o)
{
	LOCS newcenter[BINS][DIM], newdimension[DIM];
	unsigned int size[BINS];

	unsigned char b;
	unsigned int i;
	unsigned int index=0;

	//Empty bins
	for(b=0;b<BINS;b++)
	{
		bin[b].clear();
	}

	//Bin points
	for(i=0;i<n;i++)
	{
		b=Bin(&locs[orders[i]*DIM],center);
		bin[inverse[o][b]].push_back(orders[i]);
	}

	//Reorder points by placing bins together in order
	for(b=0;b<BINS;b++)
	{
		size[b]=(unsigned int)bin[b].size();
		memcpy(&orders[index],&bin[b][0],sizeof(unsigned int)*size[b]);
		index+=size[b];
	}

	//Halve all dimensions
	for(int d=0;d<DIM;d++)
	{
		newdimension[d]=dimension[d]/2;
	}

	//recursivly call
	for(b=0;b<BINS;b++)
	{
		for(int d=0;d<DIM;d++)
		{
			newcenter[b][d]=center[d]+dimension[d]/4*dir[order[o][b]][d];		
		}

		if(size[b]==n)
		{
			bool same=true;
			//Check if locations are the same
			LOCS l[DIM];
			memcpy(l,&locs[orders[0]*DIM],sizeof(LOCS)*DIM);
			i=1;
			while(same && i<n)
			{
				for(int d=0;d<DIM;d++)
				{
					if(l[d]!=locs[orders[i]])
					{
						same=false;
						break;
					}
				}

			}

			if(!same)
			{
				SerialR(orders,bin,size[b],newcenter[b],newdimension,orientation[o][b]);
			}

		}
		else if(size[b]>1 )
		{
				SerialR(orders,bin,size[b],newcenter[b],newdimension,orientation[o][b]);
		}
		orders+=size[b];
	}
}

template<int DIM, class LOCS> template<class BITS> 
void SFC<DIM,LOCS>::SerialHR(History<BITS>* orders,vector<unsigned int> *bin, unsigned int n, LOCS *center, LOCS *dimension, unsigned int o, unsigned int r, BITS history)
{
	LOCS newcenter[BINS][DIM], newdimension[DIM];
	unsigned int size[BINS];

	unsigned char b;
	unsigned int i;
	unsigned int index=0;

	//Empty bins
	for(b=0;b<BINS;b++)
	{
		bin[b].clear();
	}

	//Bin points
	for(i=0;i<n;i++)
	{
		b=inverse[o][Bin(&locs[orders[i].i*DIM],center)];
		bin[b].push_back(orders[i].i);
	}

	//Reorder points by placing bins together in order
	for(b=0;b<BINS;b++)
	{
		size[b]=(unsigned int)bin[b].size();
		
		for(unsigned int i=0;i<size[b];i++)
		{
			orders[index++].i=bin[b][i];
		}
	}

	//Halve all dimensions
	for(int d=0;d<DIM;d++)
	{
		newdimension[d]=dimension[d]/2;
	}

	BITS NextHistory;
	//recursivly call
	for(b=0;b<BINS;b++)
	{
		for(int d=0;d<DIM;d++)
		{
			newcenter[b][d]=center[d]+dimension[d]/4*dir[order[o][b]][d];		
		}

		if(r==refinements)
		{
			NextHistory= ((history<<DIM)|b);
			//save history for each point in bucket
			for(unsigned int j=0;j<size[b];j++)
			{
				orders[j].bits=NextHistory;
			}
		}
		else if( size[b]>1)
		{
			NextHistory= ((history<<DIM)|b);
			SerialHR<BITS>(orders,bin,size[b],newcenter[b],newdimension,orientation[o][b],r+1,NextHistory);
		}
		else if (size[b]==1)
		{
			NextHistory= ((history<<DIM)|b);

			LOCS *loc=&locs[orders[0].i*DIM];
			LOCS Clocs[DIM];
			LOCS dims[DIM];
			int Co=orientation[o][b];

			for(int d=0;d<DIM;d++)
			{
				Clocs[d]=center[d]+dir[order[o][b]][d]*dimension[d]/4;

				dims[d]=newdimension[d];
			}
			
			unsigned int ref=r;
			int b;
			for(;ref<refinements;ref++)
			{
				b=inverse[Co][Bin(loc,Clocs)];

				NextHistory<<=DIM;
				NextHistory|=b;

				for(int d=0;d<DIM;d++)
				{
					dims[d]/=2;
					Clocs[d]=Clocs[d]+dir[order[o][b]][d]*dims[d]/2;
				}
				Co=orientation[Co][b];
			}
			orders[0].bits=NextHistory;
		}
		orders+=size[b];
	}
}

template<int DIM, class LOCS> template<class BITS>
void SFC<DIM,LOCS>::Parallel()
{
	vector<History<BITS> > sbuf, rbuf, mbuf;
	unsigned int i;
	vector<unsigned int> n_per_proc(P);
	this->n_per_proc=&n_per_proc[0];	
	sbuf.resize(n);
	mbuf.resize(n);
	sendbuf=(void*)&(sbuf[0]);
	recievebuf=(void*)&(rbuf[0]);
	mergebuf=(void*)&(mbuf[0]);
	
	//start sending n to every other processor
	//start recieving n from every other processor
	//or
	//preform allgather (threaded preferable)
#ifdef _TIMESFC_
	start=timer->currentSeconds();
#endif
	MPI_Allgather(&n,1,MPI_INT,&n_per_proc[0],1,MPI_INT,Comm);
#ifdef _TIMESFC_
	finish=timer->currentSeconds();
	gtime+=finish-start;
#endif

#ifdef _TIMESFC_
	start=timer->currentSeconds();
#endif
	SerialH<BITS>();	//Saves results in sendbuf
#ifdef _TIMESFC_
	finish=timer->currentSeconds();
	sertime+=finish-start;
#endif
	//find index start & max
	unsigned int max_n=n;
	unsigned int istart=0;

	for(int p=0;p<P;p++)
	{
		if(n_per_proc[p]>max_n)
			max_n=n_per_proc[p];

		if(p<rank)
			istart+=n_per_proc[p];
	}
	
	//make two pointers to internal buffers for a fast copy	
	History<BITS> *c=(History<BITS>*)sendbuf;
	
	//increment indexies 
	for(i=0;i<n;i++)
	{
		c[i].i+=istart;	
	}

	//resize buffers
//	sbuf.resize(max_n);
	rbuf.resize(max_n);
//	mbuf.resize(max_n);

	sendbuf=(void*)&(sbuf[0]);
	recievebuf=(void*)&(rbuf[0]);
	mergebuf=(void*)&(mbuf[0]);

#ifdef _TIMESFC_
	start=timer->currentSeconds();
#endif
	PrimaryMerge<BITS>();
#ifdef _TIMESFC_
	finish=timer->currentSeconds();
	ptime+=finish-start;
#endif
	
#ifdef _TIMESFC_
	start=timer->currentSeconds();
#endif
	Cleanup<BITS>();
#ifdef _TIMESFC_
	finish=timer->currentSeconds();
	cleantime+=finish-start;
#endif

	orders->resize(n);

	//make pointer to internal buffers for a fast copy	
	unsigned int* o=&(*orders)[0];
	c=(History<BITS>*)sendbuf;
		
	//copy permutation to orders
	for(i=0;i<n;i++)
	{
		//COPY to orders
		o[i]=c[i].i;	
	}
	
}

#define ASCENDING 0
#define DESCENDING 1
template<int DIM, class LOCS> template<class BITS>
int SFC<DIM,LOCS>::MergeExchange(int to)
{

//	cout << rank <<  ": Merge Exchange started with " << to << endl;
	int direction= (int) (rank>to);
	BITS emax, emin;
	queue<MPI_Request> squeue, rqueue;
	unsigned int tag=0;
	unsigned int n2=n_per_proc[to];
	MPI_Request srequest, rrequest;
	MPI_Status status;
	
	History<BITS> *sbuf=(History<BITS>*)sendbuf, *rbuf=(History<BITS>*)recievebuf, *mbuf=(History<BITS>*)mergebuf;
	History<BITS> *msbuf=sbuf, *mrbuf=rbuf;
	
	//min_max exchange
	if(direction==ASCENDING)
	{
		emax=sbuf[n-1].bits;
		MPI_Isend(&emax,sizeof(BITS),MPI_BYTE,to,0,Comm,&srequest);
		MPI_Irecv(&emin,sizeof(BITS),MPI_BYTE,to,0,Comm,&rrequest);
	}
	else
	{
		emin=sbuf[0].bits;
		MPI_Isend(&emin,sizeof(BITS),MPI_BYTE,to,0,Comm,&srequest);
		MPI_Irecv(&emax,sizeof(BITS),MPI_BYTE,to,0,Comm,&rrequest);
	}
	MPI_Wait(&rrequest,&status);
	MPI_Wait(&srequest,&status);
	
	if(emax<emin)	//if exchange not needed 
	{
		return 0;
	}
//	cout << rank << ": Max-min done\n";
	
	unsigned int nsend=n_per_proc[rank];
	unsigned int nrecv=n_per_proc[to];
	//sample exchange
	unsigned int minn=min(n,n2);
	unsigned int sample_size=(int)(minn*sample_percent);


	if(sample_size>=5)
	{
//		cout << rank << " creating samples\n";
		BITS *highsample=(BITS*)mbuf, *lowsample=(BITS*)rbuf, *mysample, *theirsample;
		float stridelow,stridehigh,mystride;
		unsigned int index=0, ihigh=0,ilow=0,count=0;
		if(direction==ASCENDING)
		{
			mysample=lowsample;
			theirsample=highsample;
			mystride=stridelow=n/(float)sample_size;
			stridehigh=n2/(float)sample_size;
		}
		else
		{
			mysample=highsample;
			theirsample=lowsample;
			stridelow=n2/(float)sample_size;
			mystride=stridehigh=n/(float)sample_size;
		}
	
		//create sample
		for(unsigned int i=0;i<sample_size;i++)
		{
			index=int(mystride*i);
			mysample[i]=sbuf[index].bits;
		}
//		cout << "exchanging samples\n";
		//exchange samples
		MPI_Isend(mysample,sample_size*sizeof(BITS),MPI_BYTE,to,1,Comm,&srequest);
		MPI_Irecv(theirsample,sample_size*sizeof(BITS),MPI_BYTE,to,1,Comm,&rrequest);
	
		MPI_Wait(&rrequest,&status);
		MPI_Wait(&srequest,&status);
		
//		cout << "done exchanging samples\n";
		//merge samples
	
		while(count<minn)
		{
			if(lowsample[ilow]<=highsample[ihigh])
			{
				ilow++;
			}
			else
			{
				ihigh++;
			}
			count=int(ilow*stridelow)+int(ihigh*stridehigh);
		}
		
		if(ilow>sample_size) //handle case where ilow goes to far
		{
			ihigh+=(ilow-sample_size);
		}
		nrecv=nsend=int((ihigh+1)*stridehigh);

		if(nsend>n)
		{
			nsend=n;
		}
		if(nrecv>n2)
		{
			nrecv=n2;
		}
	}	
	
	//final exchange
//	cout << rank << ": sample done\n";
	int b;
	unsigned int block_count=0;
	int sremaining=nsend;
	int rremaining=nrecv;
	
//	cout << sremaining << " " << rremaining << endl;
	unsigned int sent=0, recvd=0, merged=0;
//	cout << rank << " Block size: " << block_size << endl;	
	if(direction==ASCENDING)
	{
		//Merge Ascending
		//Send Descending
		//Recieve Ascending
		
		//position buffers
		sbuf+=n;
		
		while(block_count<blocks_in_transit)
		{
			//send
			if(sremaining>=(int)block_size)
			{
				sbuf-=block_size;
				MPI_Isend(sbuf,block_size*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&srequest);
				squeue.push(srequest);
				sent+=block_size;
//				cout << rank << ": Sending block of size: " << block_size << endl;
				sremaining-=block_size;
			}
			else if(sremaining>0)
			{
				sbuf-=sremaining;
				MPI_Isend(sbuf,sremaining*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&srequest);
				squeue.push(srequest);
				sent+=sremaining;
//				cout << rank << ": Sending block of size: " << sremaining << endl;
				sremaining=0;
			}
			
			//recieve
			if(rremaining>=(int)block_size)
			{
				MPI_Irecv(rbuf+recvd,block_size*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&rrequest);
				rqueue.push(rrequest);
				recvd+=block_size;
//				cout << rank << ": Recieving block of size: " << block_size << endl;
				rremaining-=block_size;
			}
			else if(rremaining>0)
			{
				MPI_Irecv(rbuf+recvd,rremaining*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&rrequest);
				rqueue.push(rrequest);
				recvd+=rremaining;

//				cout << rank << ": Recieving block of size: " << rremaining << endl;
				rremaining=0;
			}
		
			block_count++;
			tag++;
		}
		while(!rqueue.empty())
		{
			MPI_Wait(&(rqueue.front()),&status);
			rqueue.pop();
			
			MPI_Get_count(&status,MPI_BYTE,&b);
			b/=sizeof(History<BITS>);
//			cout << rank << " recieved block of size\n";
			while(b>0 && merged<n)
			{
				if(mrbuf[0].bits<msbuf[0].bits)
				{
					//pull from recieve buffer
					mbuf[merged]=mrbuf[0];
					mrbuf++;
					b--;
				}
				else
				{
					//pull from send buffer
					mbuf[merged]=msbuf[0];
					msbuf++;
				}
				merged++;
			}

			//send next block
			if(sremaining>=(int)block_size)
			{
				sbuf-=block_size;
				MPI_Isend(sbuf,block_size*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&srequest);
				squeue.push(srequest);
				sent+=block_size;
				sremaining-=block_size;
			}
			else if(sremaining>0)
			{
				sbuf-=sremaining;
				MPI_Isend(sbuf,sremaining*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&srequest);
				squeue.push(srequest);
				sent+=sremaining;
				sremaining=0;
			}
			
			//recieve
			if(rremaining>=(int)block_size)
			{
				MPI_Irecv(rbuf+recvd,block_size*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&rrequest);
				rqueue.push(rrequest);
				recvd+=block_size;
				rremaining-=block_size;
			}
			else if(rremaining>0)
			{
				MPI_Irecv(rbuf+recvd,rremaining*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&rrequest);
				rqueue.push(rrequest);
				recvd+=rremaining;
				rremaining=0;
			}
		
			block_count++;
			tag++;
		
			//wait for a send, it should be done now
			MPI_Wait(&(squeue.front()),&status);
			squeue.pop();
		}

		if(merged<n)	//merge additional elements from send buff
		{
			memcpy(mbuf+merged,msbuf,(n-merged)*sizeof(History<BITS>));
		}
	}
	else
	{
		//Merge Descending
		//Send Ascending
		//Recieve Descending

		//position buffers
		mbuf+=n;
		rbuf+=n2;

		msbuf+=n-1;
		mrbuf+=n2-1;

			
		while(block_count<blocks_in_transit)
		{
			//send
			if(sremaining>=(int)block_size)
			{
//				cout << rank << " sending block of size " << block_size << endl;
				MPI_Isend(sbuf+sent,block_size*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&srequest);
				squeue.push(srequest);
				sent+=block_size;
				sremaining-=block_size;
			}
			else if(sremaining>0)
			{
//				cout << rank << " sending block of size " << sremaining << endl;
				MPI_Isend(sbuf+sent,sremaining*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&srequest);
				squeue.push(srequest);
				sent+=sremaining;
				sremaining=0;
			}

			if(rremaining>=(int)block_size)
			{
//				cout << rank << " recieving block of size " << block_size << endl;
				rbuf-=block_size;
				MPI_Irecv(rbuf,block_size*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&rrequest);
				rqueue.push(rrequest);
				recvd+=block_size;
				rremaining-=block_size;
			}
			else if(rremaining>0)
			{
//				cout << rank << " recieving block of size " << rremaining << endl;
				rbuf-=rremaining;
				MPI_Irecv(rbuf,rremaining*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&rrequest);
				rqueue.push(rrequest);
				recvd+=rremaining;
				rremaining=0;
			}
			block_count++;
			tag++;
		}
		while(!rqueue.empty())
		{
			MPI_Wait(&(rqueue.front()),&status);
			rqueue.pop();

			MPI_Get_count(&status,MPI_BYTE,&b);
			b/=sizeof(History<BITS>);
//			cout << rank << " recieved block of size\n";

			while(b>0 && merged<n)
			{
				if(mrbuf[0].bits>msbuf[0].bits) //merge from recieve buff
				{
					mbuf--;
					mbuf[0]=mrbuf[0];
					mrbuf--;
					b--;
				}
				else	//merge from send buff
				{
					mbuf--;
					mbuf[0]=msbuf[0];
					msbuf--;
				}
				merged++;



			}
			
			//send more if needed
			if(sremaining>=(int)block_size)
			{
//				cout << rank << " sending block of size " << block_size << endl;
				MPI_Isend(sbuf+sent,block_size*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&srequest);
				squeue.push(srequest);
				sent+=block_size;
				sremaining-=block_size;
			}
			else if(sremaining>0)
			{
//				cout << rank << " sending block of size " << sremaining << endl;
				MPI_Isend(sbuf+sent,sremaining*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&srequest);
				squeue.push(srequest);
				sent+=sremaining;
				sremaining=0;
			}
			if(rremaining>=(int)block_size)
			{
//				cout << rank << " recieving block of size " << block_size << endl;
				rbuf-=block_size;
				MPI_Irecv(rbuf,block_size*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&rrequest);
				rqueue.push(rrequest);
				recvd+=block_size;
				rremaining-=block_size;
			}
			else if(rremaining>0)
			{
//				cout << rank << " recieving block of size " << rremaining << endl;
				rbuf-=rremaining;
				MPI_Irecv(rbuf,rremaining*sizeof(History<BITS>),MPI_BYTE,to,tag,Comm,&rrequest);
				rqueue.push(rrequest);
				recvd+=rremaining;
				rremaining=0;
			}
			tag++;

			MPI_Wait(&(squeue.front()),&status);
			squeue.pop();
		}
		if(merged<n) //merge additional elements off of msbuf
		{
			memcpy(mbuf-n+merged,msbuf-n+merged,(n-merged)*sizeof(History<BITS>));
		}
	}
	while(!rqueue.empty())
	{
		MPI_Wait(&(rqueue.front()),&status);
		rqueue.pop();
//		cout << rank << " recieved left over block\n";
	}
	while(!squeue.empty())
	{
		MPI_Wait(&(squeue.front()),&status);
		squeue.pop();
//		cout << rank << " sent left over block\n";
	}
	swap(mergebuf,sendbuf);
	return 1;
}

struct HC_MERGE
{
	unsigned int base;
	unsigned int P;
};

template<int DIM, class LOCS> template<class BITS>
void SFC<DIM,LOCS>::PrimaryMerge()
{
	queue<HC_MERGE> q;
	HC_MERGE cur;
	bool send;
	int to;

	cur.base=0;
	cur.P=P;

	q.push(cur);
	while(!q.empty())
	{
		int base, P;
		cur=q.front();
		q.pop();

		base=cur.base;
		P=cur.P;
		send=false;
		if(rank>=base && rank<base+P/2)
		{
			send=true;
			to=rank+(P+1)/2;
		}
		else if(rank-(P+1)/2>=base && rank-(P+1)/2<base+P/2)
		{
			send=true;
			to=rank-(P+1)/2;
		}

		if(send)
		{
			MergeExchange<BITS>(to);
		}

		//make next stages

		cur.P=(P+1)/2;
		if(cur.P>1)
		{
			cur.base=base+P/2;
			q.push(cur);
		}

		cur.P=P-(P+1)/2;
		if(cur.P>1)
		{
			cur.base=base;
			q.push(cur);
		}
	}

}

template<int DIM, class LOCS> template<class BITS>
void SFC<DIM,LOCS>::Cleanup()
{
	switch(cleanup)
	{
		case BATCHERS:
			Batchers<BITS>();
			break;
		case LINEAR:
			Linear<BITS>();
			break;
	};
}
template<int DIM, class LOCS> template <class BITS>
void SFC<DIM,LOCS>::Batchers()
{
	int p, r, t, q, d;

	t=1;

	while(t<P)
		t<<=1;

	p=t>>1;
	for(;p>=1;p>>=1)
	{
		q=t>>1;
		r=0;
		d=p;

		bool more;
		do
		{
			more=false;

			if(rank<P-d && (rank&p)==r)
			{
				MergeExchange<BITS>(rank+d);
			}
			else if(rank-d>=0 && ((rank-d)&p)==r)
			{
				MergeExchange<BITS>(rank-d);
			}
			if(q!=p)
			{
				more=true;
				d=q-p;
				q>>=1;
				r=p;
			}
		}while(more);
	}
}

template<int DIM, class LOCS> template <class BITS>
void SFC<DIM,LOCS>::Linear()
{
	unsigned int i=1, c=1, val=0;
	int mod=(int)ceil(log((float)P)/log(3.0f));

	while(c!=0)
	{
		val=0;
		if(rank%2==0)	//exchange right then left
		{
			if(rank!=P-1)
			{
				val+=MergeExchange<BITS>(rank+1);	
			}

			if(rank!=0)
			{
				val+=MergeExchange<BITS>(rank-1);
			}
			
		}
		else	//exchange left then right
		{
			if(rank!=0)
			{
				val+=MergeExchange<BITS>(rank-1);
			}
			
			if(rank!=P-1)
			{
				val+=MergeExchange<BITS>(rank+1);	
			}
		}
		i++;

		if(i%mod==0)
		{
			MPI_Allreduce(&val,&c,1,MPI_INT,MPI_MAX,Comm);
		}
	}
	

}
template<int DIM, class LOCS>
void SFC<DIM,LOCS>::SetMergeParameters(unsigned int block_size, unsigned int blocks_in_transit, float sample_percent)
{
	this->block_size=block_size;
	this->blocks_in_transit=blocks_in_transit;
	this->sample_percent=sample_percent;
}
template<int DIM, class LOCS>
void SFC<DIM,LOCS>::SetRefinements(unsigned int refinements)
{
	this->refinements=refinements;
	set=set|32;
}

template<int DIM, class LOCS>
void SFC<DIM,LOCS>::SetLocalSize(unsigned int n)
{
	this->n=n;
	set=set|16;
}


template<int DIM, class LOCS>
void SFC<DIM,LOCS>::SetLocations(vector<LOCS> *locsv)
{
	if(locsv!=0)
	{
		this->locsv=locsv;
		this->locs=&(*locsv)[0];
		set=set|4;
	}
}

template<int DIM, class LOCS>
void SFC<DIM,LOCS>::SetOutputVector(vector<unsigned int> *orders)
{
	if(orders!=0)
	{
		this->orders=orders;
		set=set|8;
	}
}

template<class LOCS>
void SFC2D<LOCS>::SetCurve(Curve curve)
{
	switch(curve)
	{
		case HILBERT:
			order=horder2;
		        orientation=horient2;
			inverse=hinv2;
		      	break;
		case MORTON:
			order=morder2;
		        orientation=morient2;
			inverse=morder2;
			break;
		case GREY:
			order=gorder2;
		        orientation=gorient2;
			inverse=ginv2;
		      	break;
	}
}
		 
template<class LOCS>
void SFC2D<LOCS>::SetDimensions(LOCS wx, LOCS wy)
{
	dimensions[0]=wx;
	dimensions[1]=wy;
	set=set|1;
}

template<class LOCS>
void SFC2D<LOCS>::SetCenter(LOCS x, LOCS y)
{
	center[0]=x;
	center[1]=y;
	set=set|2;
}

template<class LOCS>
void SFC2D<LOCS>::SetRefinementsByDelta(LOCS deltax, LOCS deltay)
{
	char mask=1;
	if( (mask&set) != mask)
	{
		cout << "SFC Error: Cannot set refinements by delta until dimensions have been set\n";
	}

	refinements=(unsigned int)ceil(log(dimensions[0]/deltax)/log(2.0));
	refinements=max(refinements,(unsigned int)ceil(log(dimensions[1]/deltay)/log(2.0)));
	set=set|32;
}
	
template<class LOCS>
unsigned char  SFC2D<LOCS>::Bin(LOCS *point, LOCS *center)
{

	unsigned char bin=0;
	
	if(point[0]<center[0])
		bin|=1;

	if(point[1]<center[1])
		bin|=2;

	return bin;

}
/*****************SFC3D*********************/
template<class LOCS>
void SFC3D<LOCS>::SetCurve(Curve curve)
{
	switch(curve)
	{
		case HILBERT:
			order=horder3;
		        orientation=horient3;
			inverse=hinv3;
		      	break;
		case MORTON:
			order=morder3;
		        orientation=morient3;
			inverse=morder3;
			break;
		case GREY:
			order=gorder3;
		        orientation=gorient3;
			inverse=ginv3;
		        break;
	}
}
template<class LOCS>
void SFC3D<LOCS>::SetDimensions(LOCS wx, LOCS wy, LOCS wz)
{
	dimensions[0]=wx;
	dimensions[1]=wy;
	dimensions[2]=wz;
	set=set|1;
}
template<class LOCS>
void SFC3D<LOCS>::SetCenter(LOCS x, LOCS y, LOCS z)
{
	center[0]=x;
	center[1]=y;
	center[2]=z;

	set=set|2;
}
template<class LOCS>
void SFC3D<LOCS>::SetRefinementsByDelta(LOCS deltax, LOCS deltay, LOCS deltaz)
{
	char mask=1;
	if( (mask&set) != mask)
	{
		cout << "SFC Error: Cannot set refinements by delta until dimensions have been set\n";
	}

	refinements=(unsigned int)ceil(log(dimensions[0]/deltax)/log(2.0));
	refinements=max(refinements,(unsigned int)ceil(log(dimensions[1]/deltay)/log(2.0)));
	refinements=max(refinements,(unsigned int)ceil(log(dimensions[2]/deltaz)/log(2.0)));
	set=set|32;
}
template<class LOCS>
unsigned char  SFC3D<LOCS>::Bin(LOCS *point, LOCS *center)
{
	unsigned char bin=0;
	
	if(point[0]>=center[0])
		bin|=4;

	if(point[1]<center[1])
		bin|=2;

	if(point[2]<center[2])
		bin|=1;

	return bin;
}

} //End Namespace Uintah
#endif
