
/*
 *  ColumnMatrix.h: for RHS and LHS
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   July 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#ifndef SCI_project_ColumnMatrix_h
#define SCI_project_ColumnMatrix_h 1

#include <Datatypes/Datatype.h>
#include <Classlib/LockingHandle.h>

class ColumnMatrix;
class ostream;
typedef LockingHandle<ColumnMatrix> ColumnMatrixHandle;

class ColumnMatrix : public Datatype {
    int rows;
    double* data;
public:
    ColumnMatrix(int);
    ~ColumnMatrix();
    ColumnMatrix(const ColumnMatrix&);
    virtual ColumnMatrix* clone() const;
    ColumnMatrix& operator=(const ColumnMatrix&);
    int nrows() const;
    inline double& operator[](int) const;

    double vector_norm();
    double vector_norm(int& flops, int& memrefs);

    friend void Mult(ColumnMatrix&, const ColumnMatrix&, const ColumnMatrix&);
    friend void Mult(ColumnMatrix&, const ColumnMatrix&, const ColumnMatrix&,
		     int& flops, int& memrefs);
    friend void Sub(ColumnMatrix&, const ColumnMatrix&, const ColumnMatrix&);
    friend void Sub(ColumnMatrix&, const ColumnMatrix&, const ColumnMatrix&,
		    int& flops, int& memrefs);
    friend double Dot(const ColumnMatrix&, const ColumnMatrix&);
    friend double Dot(const ColumnMatrix&, const ColumnMatrix&,
		      int& flops, int& memrefs);
    friend void ScMult_Add(ColumnMatrix&, double s, const ColumnMatrix&,
			   const ColumnMatrix&);
    friend void ScMult_Add(ColumnMatrix&, double s, const ColumnMatrix&,
			   const ColumnMatrix&, int& flops, int& memrefs);

    void zero();
    void print(ostream&);

    // Persistent representation...
    virtual void io(Piostream&);
    static PersistentTypeID type_id;
};

#include <Classlib/Assert.h>

inline double& ColumnMatrix::operator[](int i) const
{
    return data[i];
}

#endif
