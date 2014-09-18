// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/wum/MoverList.cc
/// @brief
/// @author Mike Tyka

#include <protocols/wum/MoverList.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace  wum {

static thread_local basic::Tracer TR( "MoverList" );

void MoverList::register_mover( const std::string &name, moves::MoverCOP the_mover){
	mover_list_[ name ] = the_mover;
}

moves::MoverCOP MoverList::get_mover( const std::string &name ) const{
	TR << "Getting Mover.." << std::endl;
	std::map< std::string, moves::MoverCOP >::const_iterator iter = mover_list_.find( name );
	if( iter == mover_list_.end() ){
		utility_exit_with_message( "ERROR: Cannot find Mover named '" + name + "'" );
	}
	return iter->second();
}


}
}


