
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                   MV++ Numerical Matrix/Vector C++ Library                */
/*                             MV++ Version 1.5                              */
/*                                                                           */
/*                                  R. Pozo                                  */
/*               National Institute of Standards and Technology              */
/*                                                                           */
/*                                  NOTICE                                   */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that this permission notice appear in all copies and             */
/* supporting documentation.                                                 */
/*                                                                           */
/* Neither the Institution (National Institute of Standards and Technology)  */
/* nor the author makes any representations about the suitability of this    */
/* software for any purpose.  This software is provided ``as is''without     */
/* expressed or implied warranty.                                            */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//
//          Basic vector class (int precision)
//


                                 
#include "mvvi.h"


MV_Vector_int::MV_Vector_int()  : p_(0), dim_(0) , ref_(0){};

MV_Vector_int::MV_Vector_int(unsigned int n) : p_(new int[n]), dim_(n), 
            ref_(0)
{
    if (n > 0 && p_ == NULL)
    {
        cerr << "Error: NULL pointer in MV_Vector_int(int) constructor " << endl;
        cerr << "       Most likely out of memory... " << endl;
        exit(1);
    }
}


MV_Vector_int::MV_Vector_int(unsigned int n, const int& v) : 
        p_(new int[n]), dim_(n), ref_(0)
{
    if (n > 0 && p_ == NULL)
    {
        cerr << "Error: NULL pointer in MV_Vector_int(int) constructor " << endl;
        cerr << "       Most likely out of memory... " << endl;
        exit(1);
    }
    for (int i=0; i<n; i++)
        p_[i] = v;
}

// operators and member functions
//


MV_Vector_int& MV_Vector_int::operator=(const int & m) 
{
#ifdef TRACE_VEC
    cout << "> MV_Vector_int::operator=(const int & m)  " << endl;
#endif

    // unroll loops to depth of length 4

    int N = size();

    int Nminus4 = N-4;
    int i;

    for (i=0; i<Nminus4; )
    {
        p_[i++] = m;
        p_[i++] = m;
        p_[i++] = m;
        p_[i++] = m;
    }

    for (; i<N; p_[i++] = m);   // finish off last piece...

#ifdef TRACE_VEC
    cout << "< MV_Vector_int::operator=(const int & m)  " << endl;
#endif
    return *this;
}


MV_Vector_int& MV_Vector_int::newsize(unsigned int n)
{
#ifdef TRACE_VEC
    cout << "> MV_Vector_int::newsize(unsigned int n) " << endl;
#endif
    if (ref_ )                  // is this structure just a pointer?
    {
        {
            cerr << "MV_Vector::newsize can't operator on references.\n";
            exit(1);
        }
    }
    else
    if (dim_ != n )                     // only delete and new if
    {                                   // the size of memory is really
        if (p_) delete [] p_;           // changing, otherwise just
        p_ = new int[n];              // copy in place.
        if (n > 0 && p_ == NULL)
        {
            cerr << "Error : NULL pointer in operator= " << endl;
            exit(1);
        }
        dim_ = n;
    }

#ifdef TRACE_VEC
    cout << "< MV_Vector_int::newsize(unsigned int n) " << endl;
#endif

    return *this;
}


    


MV_Vector_int& MV_Vector_int::operator=(const MV_Vector_int & m) 
{

    int N = m.dim_;
    int i;

    if (ref_ )                  // is this structure just a pointer?
    {
        if (dim_ != m.dim_)     // check conformance,
        {
            cerr << "MV_VectorRef::operator=  non-conformant assignment.\n";
            exit(1);
        }

        // handle overlapping matrix references
        if ((m.p_ + m.dim_) >= p_)
        {
            // overlap case, copy backwards to avoid overwriting results
            for (i= N-1; i>=0; i--)
                p_[i] = m.p_[i];
        }
        else
        {
            for (i=0; i<N; i++)
                p_[i] = m.p_[i];
        }
                
    }
    else
    {
        newsize(N);

        // no need to test for overlap, since this region is new
        for (i =0; i< N; i++)       // careful not to use bcopy()
            p_[i] = m.p_[i];        // here, but int::operator= int.
    }
    return *this;   
}


MV_Vector_int::MV_Vector_int(const MV_Vector_int & m) : p_(new int[m.dim_]), 
    dim_(m.dim_) , ref_(0)
{
    if (m.dim_ > 0 && p_ == NULL)
    {
        cerr << "Error:  Null pointer in MV_Vector_int(const MV_Vector&); " << endl;
        exit(1);
    }

    int N = m.dim_;

    for (int i=0; i<N; i++)
        p_[i] = m.p_[i];
}

MV_Vector_int::MV_Vector_int(int* d, unsigned int n) : p_(new int[n]), 
      dim_(n) , ref_(0)
{
    if (n > 0 && p_ == NULL)
    {
        cerr << "Error: Null pointer in MV_Vector_int(int*, int) " << endl;
        exit(1);
    }
    for (int i=0; i<n; i++)
        p_[i] = d[i];

}



MV_Vector_int::MV_Vector_int(const int* d, unsigned int n) : p_(new int[n]), 
      dim_(n) , ref_(0)
{
    if (n > 0 && p_ == NULL)
    {
        cerr << "Error: Null pointer in MV_Vector_int(int*, int) " << endl;
        exit(1);
    }
    for (int i=0; i<n; i++)
        p_[i] = d[i];

}


MV_Vector_int MV_Vector_int::operator()(void)
{
    return MV_Vector_int(p_, dim_, MV_Vector_::ref);
}


const MV_Vector_int MV_Vector_int::operator()(void) const
{
    return MV_Vector_int(p_, dim_, MV_Vector_::ref);
}


MV_Vector_int MV_Vector_int::operator()(const MV_VecIndex &I) 
{
    // default parameters
    if (I.all())
        return MV_Vector_int(p_, dim_, MV_Vector_::ref);
    else
    {
    // check that index is not out of bounds
    //
        if ( I.end() >= dim_)
        {
            cerr << "MV_VecIndex: (" << I.start() << ":" << I.end() << 
                ") too big for matrix (0:" << dim_ - 1 << ") " << endl;
            exit(1);
        }
        return MV_Vector_int(p_+ I.start(), I.end() - I.start() + 1,
            MV_Vector_::ref);
    }
}


const MV_Vector_int MV_Vector_int::operator()(const MV_VecIndex &I) const
{
    // check that index is not out of bounds
    //
    if (I.all())
        return MV_Vector_int(p_, dim_, MV_Vector_::ref);
    else
    {
      if ( I.end() >= dim_)
      {
        cerr << "MV_VecIndex: (" << I.start() << ":" << I.end() << 
                ") too big for matrix (0:" << dim_ - 1 << ") " << endl;
        exit(1);
      }
      return MV_Vector_int(p_+ I.start(), I.end() - I.start() + 1,
            MV_Vector_::ref);
    }
}


MV_Vector_int::~MV_Vector_int()
{
        if (p_ && !ref_ ) delete [] p_;
}


ostream&   operator<<(ostream& s, const MV_Vector_int& V)
{
    int N = V.size();

    for (int i=0; i< N; i++)
        s << V(i) << endl;
    
    return s;
}


