
/*
 *  Array1.cc: Implementation of dynamic 1D array class
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   Feb. 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#include "Array1.h"

template<class T>
Array1<T>::Array1(const Array1<T>& a)
{
    _size=a._size;
    nalloc=_size;
    objs=new T[_size];
    for(int i=0;i<_size;i++)objs[i]=a.objs[i];
    nalloc=_size;
    default_grow_size=a.default_grow_size;
}

template<class T>
Array1<T>& Array1<T>::operator=(const Array1<T>& copy)
{
    if (objs)delete [] objs;
    _size=copy._size;
    nalloc=_size;
    objs=new T[_size];
    for(int i=0;i<_size;i++)objs[i]=copy.objs[i];
    nalloc=_size;
    default_grow_size=copy.default_grow_size;
    return(*this);
}

template<class T>
Array1<T>::Array1(int size, int gs, int asize)
{
    ASSERT(size >= 0);
    default_grow_size=gs;
    if(size){
	if(asize==-1){
	    objs=new T[size];
	    _size=size;
	    nalloc=_size;
	} else {
	    objs=new T[asize];
	    _size=size;
	    nalloc=asize;
	}
    } else {
	if(asize==-1){
	    objs=0;
	    _size=0;
	    nalloc=0;
	} else {
	    objs=new T[asize];
	    _size=0;
	    nalloc=asize;
	}
    }
    nalloc=_size;
}	

template<class T>
Array1<T>::~Array1()
{
    if(objs)delete [] objs;
}

template<class T>
void Array1<T>::grow(int count, int grow_size)
{
    int newsize=_size+count;
    if(newsize>nalloc){
	// Reallocate...
	int gs1=newsize>>2;
	int gs=gs1>grow_size?gs1:grow_size;
	int newalloc=newsize+gs;
	T* newobjs=new T[newalloc];
	if(objs){
	    for(int i=0;i<_size;i++){
		newobjs[i]=objs[i];
	    }
	    delete[] objs;
	}
	objs=newobjs;
	nalloc=newalloc;
    }
    _size=newsize;
}

template<class T>
void Array1<T>::add(const T& obj)
{
    grow(1, default_grow_size);
    objs[_size-1]=obj;
}

template<class T>
void Array1<T>::insert(int idx, const T& obj)
{
    grow(1, default_grow_size);
    for(int i=_size-1;i>idx;i--)objs[i]=objs[i-1];
    objs[idx]=obj;
}

template<class T>
void Array1<T>::remove(int idx)
{
    _size--;
    for(int i=idx;i<_size;i++)objs[i]=objs[i+1];
}

template<class T>
void Array1<T>::remove_all()
{
    _size=0;
}

template<class T>
void Array1<T>::resize(int newsize)
{
    if(newsize > _size)
	grow(newsize-_size);
    else
	_size=newsize;
}
