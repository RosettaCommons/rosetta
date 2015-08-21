// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/WorkUnitList.cc
/// @brief
/// @author Mike Tyka

#include <protocols/wum/WorkUnitList.hh>
#include <protocols/wum/WorkUnitBase.fwd.hh>
#include <protocols/wum/WorkUnitBase.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace  wum {

static thread_local basic::Tracer TR( "WorkUnitList" );

void WorkUnitList::register_work_unit( const std::string &name, WorkUnitBaseOP the_work_unit){
	the_work_unit->set_wu_type( name );
	work_unit_list_[ name ] = the_work_unit;
}

const WorkUnitBaseCOP WorkUnitList::get_work_unit( const std::string &name ) const{
	TR.Debug << "Getting WorkUnit.." << std::endl;
	std::map< std::string, WorkUnitBaseCOP >::const_iterator iter = work_unit_list_.find( name );
	if ( iter == work_unit_list_.end() ) {
		utility_exit_with_message( "ERROR: Cannot find WorkUnit named '" + name + "'" );
	}
	return iter->second;
}

// Copy over all the WorkUnitOPs from source to the current workUnitlist.
void WorkUnitList::merge( const WorkUnitList & source ){
	for ( std::map< std::string, WorkUnitBaseCOP >::const_iterator it = source.work_unit_list_.begin(),
			end = source.work_unit_list_.end(); it != end; ++it ) {
		work_unit_list_[ it->first ] = it->second;
	}
}

}
}


