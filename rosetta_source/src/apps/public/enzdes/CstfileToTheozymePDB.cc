// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/public/enzdes/CstFileToTheozymePDB.cc
/// @brief
/// @author Florian Richter, floric@u.washington.edu, june 2010


#include <devel/init.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/match.OptionKeys.gen.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTree.hh>
#include <basic/Tracer.hh>

void
match_main();

int main( int argc, char * argv [] )
{

	devel::init( argc, argv );

	match_main();
}

void
match_main()
{

	basic::Tracer tr( "apps.public.enzdes.CstfileToTheozymePDB.cc" );

	std::string cstfile_name( basic::options::option[basic::options::OptionKeys::match::geometric_constraint_file]() );

	protocols::toolbox::match_enzdes_util::EnzConstraintIOOP enz_io = new protocols::toolbox::match_enzdes_util::EnzConstraintIO( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	enz_io->read_enzyme_cstfile( cstfile_name );

	protocols::toolbox::match_enzdes_util::InvrotTreeOP invrot_tree = new protocols::toolbox::match_enzdes_util::TheozymeInvrotTree( enz_io );
	invrot_tree->generate_targets_and_inverse_rotamers();

	//test whether the input cstfile is in the same directory or somewhere else
	std::string::size_type const slash_loc = cstfile_name.find_last_of( '/' );
	std::string outname_base;
	if( slash_loc == std::string::npos ) { // same directory
		outname_base = "PDB_Model_"+cstfile_name;
	}
	else{
		outname_base = "PDB_Model_"+cstfile_name.substr(slash_loc+1, cstfile_name.size() );
	}
	invrot_tree->dump_invrots_tree_as_multimodel_pdbs( outname_base );
}

