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

#include <core/fragment/FragID.hh>
#include <core/fragment/FragData.hh>
#include <core/types.hh>


/* Wraputils */
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>

#include <utility/stream_util.hh>

#include <ostream>
#include <set>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>


namespace bp = boost::python;
using namespace std;
using namespace utility;

typedef bp::return_value_policy< bp::reference_existing_object > CP_REF;
typedef bp::return_value_policy< bp::copy_const_reference >      CP_CCR;
typedef bp::return_value_policy< bp::copy_non_const_reference >  CP_CNCR;

// utility::vector1 --------------------------------------------------------------------------------------------------------------------
template <class T>
std::string vector1_repr(utility::vector1<T> const & v)
{
    std::ostringstream os;

    os << "[";
    for(unsigned int i=1; i<=v.size(); i++) {
        os << v[i] << ", ";
    }
    os << "]";
    return os.str();
}

template< class TT > inline void vector1_set( vector1<TT> & v, core::Size const & i, TT const & val ) { v[i] = val; }
template< class TT > inline core::Size vector1_len( vector1<TT> & v ) { return v.size(); }

template< class TT > inline std::string vector1_str( vector1<TT> & v ) { std::ostringstream s; s<<v; return s.str(); }

template< class TT > inline typename vector1<TT>::iterator vector1_begin( vector1<TT> & v ) { return v.begin(); }
template< class TT > inline typename vector1<TT>::iterator vector1_end  ( vector1<TT> & v ) { return v.end(); }

template< class TT > inline void vector1_reserve( vector1<TT> & v, core::Size n) { v.reserve(n); }
template< class TT > inline void vector1_resize( vector1<TT> & v, core::Size n) { v.resize(n); }

template< class Htype, class CP, class CP_const>
void wrap_vector1(char * name) {
  typedef vector1<Htype> Ttype;
  typedef vectorL<1,Htype,allocator<Htype> > Btype;
  typedef vector<Htype> Vtype;
  bp::class_<Ttype>(name)
    .def( bp::init< core::Size >() )
    .def( bp::init< vector1<Htype> const & >() )
    // .def( bp::init< size_t, TT >() )
    .def("__getitem__"
        , (Htype const & (Ttype::*)(core::Size const) const)( &Ttype::at )
        , CP_const()    )
    .def("__getitem__"
        , (Htype & (Ttype::*)(core::Size const))( &Ttype::at )
        , CP()        )
    .def("__setitem__"
        , &vector1_set<Htype> )
    .def("append"
        , (Btype & (Btype::*)(Htype const &))( &Btype::add_back )
        , bp::return_value_policy< bp::reference_existing_object >()        )
    .def("__len__", & vector1_len<Htype> )
    .def("__iter__", bp::range(&vector1_begin<Htype>,&vector1_end<Htype>) )

    //.def("__str__", & vector1_str<Htype> )
    //.def("__str__", & vector1_repr<Htype> )
    //.def( bp::self_ns::str( bp::self ) )

    .def("reserve", &vector1_reserve<Htype> )
    .def("resize", &vector1_resize<Htype> )

  ;
}
/*Wraputils end*/

void __fragment_by_hand_beginning__()
{
	wrap_vector1<core::fragment::FragDataOP,       CP_CNCR, CP_CCR>("vector1_FragDataOP");
	wrap_vector1<core::fragment::FragData,       CP_CNCR, CP_CCR>("vector1_FragData");

	wrap_vector1<core::fragment::Frame,       CP_CNCR, CP_CCR>("vector1_Frame");
	wrap_vector1<core::fragment::FrameOP,       CP_CNCR, CP_CCR>("vector1_FrameOP");
}
