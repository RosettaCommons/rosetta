// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/py/PyHelper.hh
/// @author Sergey Lyskov
///
/// @note Defining various helper function for PyRosetta.


#ifndef INCLUDED_utility_py_PyHelper_hh
#define INCLUDED_utility_py_PyHelper_hh

#include <boost/python.hpp>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

#include <boost/python/ptr.hpp>
#include <Python.h>

namespace utility { namespace py {

/// @brief Function to create argument for Python subclassing
// template< class T > utility::pointer::access_ptr<T> PyAP( T & o) { return utility::pointer::access_ptr<T> ( & o ); }


/// @brief Function to create argument for Python subclassing
//template< class T > PyObject * to_py_object( T & o) { return boost::python::object(o).ptr(); }
//template< class T > boost::python::object &to_py_object( T & o) { return boost::python::object(o); }

//template< class T > utility::pointer::weak_ptr<T> to_py_object( T & o ) { return utility::pointer::weak_ptr<T> (o); }

//template< class T > T * to_py_object( T & o ) { return &o; }


//template< class T > boost::python::object & to_py_object( T & o) { return boost::python::object( boost::python::handle<>(boost::python::borrowed(o) ) ); }


// template< class T >  T * wrap_access_pointer_get_function( pointer::access_ptr<T> rs ) {  return rs.get(); }
template< class T > utility::pointer::shared_ptr<T> wrap_access_pointer_get_function( utility::pointer::weak_ptr<T> rs ) {  return rs.lock(); }
template< class T > utility::pointer::shared_ptr<T const> wrap_access_const_pointer_get_function( utility::pointer::weak_ptr<T const> rs ) {  return rs.lock(); }


template< class T >
void wrap_access_pointer(std::string class_name)
{
    // //boost::python::implicitly_convertible< utility::pointer::access_ptr< T >
    // //                                     , utility::pointer::access_ptr< T const > >();

    boost::python::implicitly_convertible< utility::pointer::weak_ptr< T >, utility::pointer::weak_ptr< T const > >();


    // boost::python::class_< utility::pointer::access_ptr< T > >( std::string(class_name+"AP").c_str() )
    //     .def("get", (  T * (*)( utility::pointer::access_ptr<T> )  )( & wrap_access_pointer_get_function<T> )
    //          , boost::python::return_value_policy< boost::python::reference_existing_object >() );

    boost::python::class_< utility::pointer::weak_ptr< T > >( std::string(class_name+"AP").c_str() )
        .def("get", (  utility::pointer::shared_ptr<T> (*)( utility::pointer::weak_ptr<T> )  )( & wrap_access_pointer_get_function<T> ) );
	//, boost::python::return_value_policy< boost::python::reference_existing_object >() );


    // boost::python::class_< utility::pointer::access_ptr< T const > >( std::string(class_name+"CAP").c_str() )
    //     .def("get", (  T const * (*)( utility::pointer::access_ptr<T const > )  )( & wrap_access_pointer_get_function<T const> )
    //          , boost::python::return_value_policy< boost::python::reference_existing_object >() );

    boost::python::class_< utility::pointer::weak_ptr< T const > >( std::string(class_name+"CAP").c_str() )
        .def("get", (  utility::pointer::shared_ptr<T const> (*)( utility::pointer::weak_ptr<T const> ) )( & wrap_access_const_pointer_get_function<T> ) );
}


template< class T>
struct COP_to_Python_converter
{
	// static PyObject *convert( utility::pointer::owning_ptr< T const > const & o ) {
	// 	return  boost::python::incref( boost::python::object( new utility::pointer::owning_ptr< T > ( (T*)o.get() ) ).ptr() );
	// }

	static PyObject *convert( utility::pointer::shared_ptr< T const > const & o ) {
		return boost::python::incref( boost::python::object( boost::const_pointer_cast< T >(o) ).ptr() );
	}
};


/*
// COP_OP_convertor
template< class T >
operator utility::pointer::access_ptr< T > (utility::pointer::access_ptr< T const > const & cop )
{
	return utility::pointer::access_ptr< T >( (T*)cop.get() );
} */


} } // namespace py namespace utility

#endif // INCLUDED_utility_py_PyHelper_hh
