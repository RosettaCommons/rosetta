// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Sergey Lyskov
///

#include "boost/python.hpp"

#include <utility/pointer/access_ptr.hh>

#include <core/chemical/orbitals/OrbitalTypeSet.hh>

// template< class T >
// T * getCAP( utility::pointer::access_ptr<T> rs ) {
//   T & rs_ref( *rs );
//   T * rs_ptr = &rs_ref;
//   return rs_ptr;
// }


void __orbitals_by_hand_beginning__()
{
	// boost::python::class_< utility::pointer::access_ptr< core::chemical::orbitals::OrbitalTypeSet const> >("core___chemical___orbitals__OrbitalTypeSetCAP");

	// boost::python::def("orbitals___getCAP"
	// 	   , (  core::chemical::orbitals::OrbitalTypeSet const * (*)( utility::pointer::access_ptr< core::chemical::orbitals::OrbitalTypeSet const> )  )
    //      			   ( & getCAP< core::chemical::orbitals::OrbitalTypeSet const> )
    //      , boost::python::return_value_policy< boost::python::reference_existing_object >() );

}
