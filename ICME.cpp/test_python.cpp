#include <iostream>
#include <Python.h>
#include <numpy/arrayobject.h>
//#include "C:\Python26\Lib\site-packages\numpy\core\include\numpy\arrayobject.h"

void output_result( PyObject* rslt )
{
   if( !rslt ) {
      printf("output_result: argument must be a valid pointer\n");
      return;
   }

   // output scalar result
   if( PyFloat_Check(rslt) ) {
      printf( "result: %f\n", PyFloat_AsDouble(rslt) );
   }

   if( PyArray_Check(rslt) ) {
      PyArrayObject* obj = PyArray_GETCONTIGUOUS( (PyArrayObject*) rslt );
      int ndims = obj->nd;
      npy_intp* dims = obj->dimensions; // not copying data
      double* data = (double*) obj->data; // not copying data
      int i, j, k = 0;

      // output vector result
      if( ndims == 1 ) {
         for( i=0; i<dims[0]; i++ )
            printf( "element: %i  value: %f\n", i, data[k++] );
         printf("\n");
      }

      // output matrix reult
      else if( ndims == 2 ) {
         for( i=0; i<dims[0]; i++ ) {
            for( j=0; j<dims[1]; j++ )
               printf( "%f  ", data[k++] );
            printf( "\n" );
         }
      }

      // output N-D result
      else {
         for( i=0; i<ndims; i++ )
            for( j=0; j<dims[i]; j++ )
               printf( "dimension: %i  element: %i  value: %f\n", i, j, data[k++] );
      }

      // clean
      Py_XDECREF(obj);
   }
}

void handle_error( PyObject* fe )
{
   // get exception info
   PyObject *type, *value, *traceback;
   PyErr_Fetch( &type, &value, &traceback );
   PyErr_NormalizeException( &type, &value, &traceback );

   // create a argument for "format exception"
   PyObject* args = PyTuple_New(3);
   PyTuple_SetItem( args, 0, type );
   PyTuple_SetItem( args, 1, value );
   PyTuple_SetItem( args, 2, traceback );

   // get a list of string describing what went wrong
   PyObject* output = PyObject_CallObject( fe, args );

   // print error message
   int i, n = PyList_Size( output );
   for( i=0; i<n; i++ ) printf( "%s", PyString_AsString( PyList_GetItem( output, i ) ) );

   // clean up
   Py_XDECREF( args );
   Py_XDECREF( output );
}


int main( int argc, char* argv[] )
{
   PyObject *mod1, *mod2, *dict1, *dict2, *func1, *func2, *fexcp,
            *expr, *script,
            *scalar1, *scalar2, *vec, *mat,
            *vars, *args1, *args2, *rslt1, *rslt2, *rslt3;

   // launch the python interpreter
   Py_Initialize();

   // this macro is defined be NumPy and must be included
   import_array1(-1);

   // load module "traceback" (for error handling) and module from file "Test.py"

   PyObject *sys_module = PyImport_ImportModule("sys");  /* New Reference */
   PyObject *sys_dict = PyModule_GetDict(sys_module);               /* Borrowed Reference, never fails */
   PyObject *sys_path = PyMapping_GetItemString(sys_dict, "path");  /* New reference */
   PyObject *add_value = PyString_FromString("./");  /* New reference */
   PyList_Append(sys_path, add_value);
   Py_DECREF(add_value);
   Py_DECREF(sys_path);
   Py_DECREF(sys_module);


   mod1 = PyImport_ImportModule( "traceback" );
   mod2 = PyImport_ImportModule( "Test" );
   if( mod1 && mod2 )
   {
      // get dictionary of available items in the modules
      dict1 = PyModule_GetDict(mod1);
      dict2 = PyModule_GetDict(mod2);

      // grab the functions we are interested in
      fexcp = PyDict_GetItemString(dict1, "format_exception");
      func1 = PyDict_GetItemString(dict2, "run_expression");
      func2 = PyDict_GetItemString(dict2, "run_script");
      if( func1 && func2 && fexcp )
      {
         if( PyCallable_Check(func1) && PyCallable_Check(func2) && PyCallable_Check(fexcp) )
         {
            // define a python expression and a python script
            expr = PyString_FromString("m[0:3,0]");
            script = PyString_FromString( "import Test\n"
                                          "import numpy\n"
                                          "def myfunction( a ):\n"
                                          "   '''do something'''\n"
                                          "   return a[:,0]\n"
                                          "\n"
                                          "M = numpy.array( [[ 1, 2, 3 ], [ 4, 5, 6 ]], float )\n"
                                          "value = myfunction( M )\n"
                                          "print M\n"
                                          "print value\n" );

            // define some scalars
            scalar1 = PyFloat_FromDouble(5.0);
            scalar2 = PyFloat_FromDouble(4.0);

            // define a vector
            double* v = new double[3];
            v[0] = 1.0;
            v[1] = 2.0;
            v[2] = 3.0;
            npy_intp vdim[] = { 3 };
            vec = PyArray_SimpleNewFromData( 1, vdim, PyArray_DOUBLE, v );

            // define a matrix using the classical C-format:
            // pointer-to-pointer-to-numeric type, last dimension running fastest
            //
            // setup pointers
            const int nrow = 3, ncol = 4, nelem = nrow*ncol;
            double** m = new double*[nrow];
            m[0] = new double[nelem];
            m[1] = m[0] + ncol;
            m[2] = m[1] + ncol;
            //
            // fill in values
            m[0][0] = 1.0; m[0][1] = 2.0; m[0][2] = 5.0; m[0][3] = 34;
            m[1][0] = 5.0; m[1][1] = 1.0; m[1][2] = 8.0; m[1][3] = 64;
            m[2][0] = 8.0; m[2][1] = 0.0; m[2][2] = 3.0; m[2][3] = 12;
            npy_intp mdim[] = { nrow, ncol };
            mat = PyArray_SimpleNewFromData( 2, mdim, PyArray_DOUBLE, m[0] );

            // put variables into dictionary of local variables for python
            vars = PyDict_New();
            PyDict_SetItemString( vars, "x", scalar1 );
            PyDict_SetItemString( vars, "y", scalar2 );
            PyDict_SetItemString( vars, "v", vec );
            PyDict_SetItemString( vars, "m", mat );

            // create arguments for "run_expression"
            args1 = PyTuple_New(2);
            PyTuple_SetItem(args1, 0, expr);
            PyTuple_SetItem(args1, 1, vars);

            // create argument for "run_script"
            Py_INCREF(vars);
            args2 = PyTuple_New(2);
            PyTuple_SetItem(args2, 0, script);
            PyTuple_SetItem(args2, 1, vars);

            // execute single expression
            rslt1 = PyObject_CallObject(func1, args1);
            if( !rslt1 )
               handle_error( fexcp );
            else
               output_result( rslt1 );
            Py_XDECREF(rslt1);

            // execute script
            rslt2 = PyObject_CallObject(func2, args2);
            if( !rslt2 )
               handle_error( fexcp );
            else {
               rslt3 = PyDict_GetItemString( rslt2, "value" );
               if( !rslt3 )
                  printf( "Requested variable not defined during evaluation of script\n" );
               else
                  output_result( rslt3 );
               Py_XDECREF(rslt3);
            }
            Py_XDECREF(rslt2);

            //
            Py_XDECREF(args1); // also decrements expr and vars
            Py_XDECREF(args2); // also decrements script and vars
            Py_XDECREF(mat);
            delete[] m[0];
            delete[] m;
            Py_XDECREF(vec);
            delete[] v;
            Py_XDECREF(scalar2);
            Py_XDECREF(scalar1);
            Py_XDECREF(script);
            Py_XDECREF(expr);
         }
      }
      Py_XDECREF(func1);
      Py_XDECREF(func2);
      Py_XDECREF(fexcp);
   }
   Py_XDECREF(mod1);
   Py_XDECREF(mod2);

   // (to keep debug console open on windows)
//   getchar();

//   Py_Finalize();

   return 0;
}

