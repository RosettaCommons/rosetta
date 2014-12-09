// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include "boost/python.hpp"

#include "utility/vector1.hh"
#include "utility/pointer/owning_ptr.hh"

#include <protocols/jobdist/Jobs.hh>

namespace bp = boost::python;
using namespace std;
using namespace utility;

// begin blatant code duplication from _utility__by_hand.cc

/*
template <class T>
std::ostream& operator <<(std::ostream &os, utility::vector1<T> const & v)
{
	os << "[";
	for(unsigned int i=1; i<=v.size(); i++) {
		os << v[i] << ", ";
	}
	os << "]";
	return os;
} */

template <class T>
std::string vector1_repr(utility::vector1<T> const & v)
{
	std::ostringstream os;

	os << "~~~[";
	for(unsigned int i=1; i<=v.size(); i++) {
		os << v[i] << ", ";
	}
	os << "]";
	return os.str();
}

template< class TT > inline void vector1_set( vector1<TT> & v, size_t const & i, TT const & val ) { v[i] = val; }
template< class TT > inline std::size_t vector1_len( vector1<TT> & v ) { return v.size(); }

template< class TT > inline std::string vector1_str( vector1<TT> & v ) { std::ostringstream s; s<<v; return s.str(); }

template< class TT > inline typename vector1<TT>::iterator vector1_begin( vector1<TT> & v ) { return v.begin(); }
template< class TT > inline typename vector1<TT>::iterator vector1_end  ( vector1<TT> & v ) { return v.end(); }
template< class Htype, class CP, class CP_const>
void wrap_vector1(char * name) {
  typedef vector1<Htype> Ttype;
  typedef vectorL<1,Htype,allocator<Htype> > Btype;
  typedef vector<Htype> Vtype;
  bp::class_<Ttype>(name)
    .def( bp::init< size_t >() )
    .def( bp::init< vector1<Htype> const & >() )
    // .def( bp::init< size_t, TT >() )
    .def("__getitem__"
        , (Htype const & (Ttype::*)(size_t const) const)( &Ttype::at )
        , CP_const()    )
    .def("__getitem__"
        , (Htype & (Ttype::*)(size_t const))( &Ttype::at )
        , CP()        )
    .def("__setitem__"
        , &vector1_set<Htype> )
    .def("append"
        , (Btype & (Btype::*)(Htype const &))( &Btype::add_back )
        , bp::return_value_policy< bp::reference_existing_object >()        )
    .def("__len__", & vector1_len<Htype> )
    .def("__iter__", bp::range(&vector1_begin<Htype>,&vector1_end<Htype>) )

    //.def("__str__", & vector1_str<Htype> )
    .def("__str__", & vector1_repr<Htype> )
	//.def( bp::self_ns::str( bp::self ) )
  ;
}

template< class Htype, class CP, class CP_const>
void wrap_vector1_part(char * name) {
  typedef vector1<Htype> Ttype;
  typedef vectorL<1,Htype,allocator<Htype> > Btype;
  typedef vector<Htype> Vtype;
  bp::class_<Ttype>(name)
    .def("__getitem__"
        , (Htype const & (Ttype::*)(size_t const) const)( &Ttype::at )
        , CP_const()    )
    .def("__getitem__"
        , (Htype & (Ttype::*)(size_t const))( &Ttype::at )
        , CP()        )
    .def("__setitem__"
        , &vector1_set<Htype> )
    .def("append"
        , (Btype & (Btype::*)(Htype const &))( &Btype::add_back )
        , bp::return_value_policy< bp::reference_existing_object >()        )
    .def("__len__", & vector1_len<Htype> )
    .def("__iter__", bp::range(&vector1_begin<Htype>,&vector1_end<Htype>) )
  ;
}

// end copy

// begin blatant code duplication from _utility__by_hand.cc

void __jobdist_by_hand_beginning__()
{
	using namespace pointer;
	typedef bp::return_value_policy< bp::reference_existing_object > CP_REF;
	typedef bp::return_value_policy< bp::copy_const_reference >      CP_CCR;
	typedef bp::return_value_policy< bp::copy_non_const_reference >  CP_CNCR;

	wrap_vector1<protocols::jobdist::BasicJobOP,       CP_CNCR, CP_CCR>("vector1_basicjob");
}
