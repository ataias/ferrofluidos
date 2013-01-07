//  Copyright Joel de Guzman 2002-2004. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)
//  Hello World Example from the tutorial
//  [Joel de Guzman 10/9/2002]

#include <stdheader.hpp>

#define WRAP_PYTHON 1
#if WRAP_PYTHON
#include <boost/python/detail/wrap_python.hpp>
#include <boost/python.hpp>
#include <numpy/arrayobject.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#endif

using namespace std;
using namespace Eigen;
using namespace boost::python;

class FooClass
{
public:
        FooClass( int new_m , int new_n);
        ~FooClass();

        template<typename Derived>
        int foo(const MatrixBase<Derived>& barIn, MatrixBase<Derived>& barOut);
#if WRAP_PYTHON
        int foo_python(PyObject* barIn, PyObject* barOut);
#endif
private:
        int m, n;
};

FooClass::FooClass( int new_m , int new_n){
        m = new_m;
        n = new_n;
}
FooClass::~FooClass(){
}

template<typename Derived>
int FooClass::foo(const MatrixBase<Derived>& barIn, MatrixBase<Derived>& barOut){
		MatrixXd bar(m,n);
		Matrix2d hey;
		hey << 1,2,
				3,4;
		bar = barIn;
		cout << endl << "check" << bar << endl;
        barOut = barIn*3.0;  // Some trivial placeholder computation.
}

#if WRAP_PYTHON
int FooClass::foo_python(PyObject* barIn, PyObject* barOut){
        Map<MatrixXd> _barIn((double *) PyArray_DATA(barIn),m,n);
        Map<MatrixXd> _barOut((double *) PyArray_DATA(barOut),m,n);
        return foo(_barIn, _barOut);
}

BOOST_PYTHON_MODULE(libfooclass)
{
    class_<FooClass>("FooClass", init<int,int>(args("m","n")))
        .def("foo", &FooClass::foo_python)
    ;
}
#endif
