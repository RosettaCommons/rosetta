// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/doug/number_of_residuetypes.cc
/// @brief Simply print the number of ResidueTypes in the default ResidueTypeSet
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// devel headers
#include <devel/init.hh>

// core headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
//#include <core/chemical/ResidueType.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>

// c++ headers
#include <iostream>

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.doug.number_of_residuetypes" );

int
main( int argc, char * argv [] )
{
	using namespace core;
	using namespace chemical;

	try {

		devel::init( argc, argv );

		ResidueTypeSetCOP rts_fa_std( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		ResidueTypeSetCOP rts_centroid( ChemicalManager::get_instance()->residue_type_set( CENTROID ) );

		TR << "The FA_STANDARD ResidueTypeSet contains " << rts_fa_std->residue_types_DO_NOT_USE().size() << " ResidueTypes" << std::endl;
		TR << "The CENTROID ResidueTypeSet contains " << rts_centroid->residue_types_DO_NOT_USE().size() << " ResidueTypes" << std::endl;

		TR << "************************************d**o**n**e***********************************" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
