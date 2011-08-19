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
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/match.OptionKeys.gen.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <fstream>

void
match_main();

/// @brief small function to parse the ALGORITHM_INFO blocks of a cstfile
/// for whether they contain a SECONDARY_MATCH: UPSTREAM_CST <x> tag
/// and if so, return this information
std::pair< core::Size, core::Size >
target_downstream_res(
	protocols::toolbox::match_enzdes_util::MatchConstraintFileInfoCOP mcfi );

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
	//1. figure out what our ligand / downstream residue is
	protocols::toolbox::match_enzdes_util::EnzConstraintIOOP enz_io = new protocols::toolbox::match_enzdes_util::EnzConstraintIO( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	enz_io->read_enzyme_cstfile( cstfile_name );
	utility::vector1< core::chemical::ResidueTypeCAP > ds_restypes( enz_io->mcfi_list( 1 )->mcfi( 1 )->allowed_restypes( enz_io->mcfi_list( 1 )->mcfi( 1 )->downstream_res() ) );
	if( ds_restypes.size() == 0 ) utility_exit_with_message("Apparently there are no downstream residues for the first block in the constraint file.");
	core::conformation::ResidueCOP ligres( core::conformation::ResidueFactory::create_residue( *(ds_restypes[1]) ) );
	if( ds_restypes.size() > 1 ) tr << "More than one downstream residue types were detected in the first cstfile block. Building theozyme model according to " << ligres->type().name() << "." << std::endl;

	//2. get the inverse rotamers
	core::Size num_geomcst( enz_io->num_mcfi_lists() );
	utility::vector1< std::list< core::conformation::ResidueCOP > > theozyme_res( num_geomcst );
	utility::vector1< std::list< core::conformation::ResidueCOP >::const_iterator > res_iterators( num_geomcst );
	for( core::Size i =1; i <= num_geomcst; ++i ){
		core::conformation::ResidueCOP downstream_res( ligres );
		std::pair< core::Size, core::Size > targ_ds_res( target_downstream_res( enz_io->mcfi_list( i )->mcfi( 1 ) ) );
		if( targ_ds_res.first > i ) utility_exit_with_message("Downstream residue given for constraint block "+utility::to_string( i )+" seems to be in a later constraint block ( " +utility::to_string( targ_ds_res.first )+" ). This is not legal." );
		if( targ_ds_res.second != 1 ) downstream_res = *(theozyme_res[ targ_ds_res.first ].begin());
		std::list< core::conformation::ResidueCOP > cur_inv_rots( enz_io->mcfi_list( i )->inverse_rotamers_against_residue( enz_io->mcfi_list( i )->mcfi(1)->downstream_res(), downstream_res ) );
		theozyme_res[i].splice( theozyme_res[i].end(), cur_inv_rots );
		runtime_assert( theozyme_res[i].size() > 0 );
		res_iterators[i] = theozyme_res[i].begin();
	}

	//3. write shit
	core::Size atomcounter(0), modelcount(1);
	bool all_res_iterators_at_end( false );
	std::string outname;
	//test whether the input cstfile is in the same directory or somewhere else
	std::string::size_type const slash_loc = cstfile_name.find_last_of( '/' );
	if( slash_loc == std::string::npos ) { // same directory
		outname = "PDB_Model_"+cstfile_name+".pdb";
	}
	else{
		outname = "PDB_Model_"+cstfile_name.substr(slash_loc+1, cstfile_name.size() )+".pdb";
	}
	std::cout << "PDB model will be written to file " << outname << "... " << std::endl;

	std::ofstream file_out( outname.c_str() );
	file_out << "MODEL   1\n";
	core::io::pdb::dump_pdb_residue( *ligres, atomcounter, file_out );

	while( ! all_res_iterators_at_end ){
		all_res_iterators_at_end = true;
		for( core::Size i =1; i <= num_geomcst; ++i ){
			if( res_iterators[ i ] != theozyme_res[i].end() ){
				all_res_iterators_at_end = false;
				core::io::pdb::dump_pdb_residue( **(res_iterators[ i ]), atomcounter, file_out );
				res_iterators[i]++;
			}
		}
		if( !all_res_iterators_at_end ){
			file_out << "ENDMDL \n";
			file_out << "MODEL    "+utility::to_string( ++modelcount )+"\n";
		}
	}
	file_out.close();
}

std::pair< core::Size, core::Size >
target_downstream_res(
	protocols::toolbox::match_enzdes_util::MatchConstraintFileInfoCOP mcfi )
{
	std::pair< core::Size, core::Size> to_return(1,1);
	std::map< std::string, utility::vector1< std::string > > const &
		alg_info( mcfi->algorithm_inputs() );

	if ( alg_info.find( "match" ) == alg_info.end() ) return to_return;

	utility::vector1< std::string > const & info( alg_info.find( "match" )->second );
	for ( core::Size ll = 1; ll <= info.size(); ++ll ) {
		std::string llstr = info[ ll ];
		std::istringstream llstream( llstr );
		std::string first, second;
		llstream >> first >> second;
		if( first == "SECONDARY_MATCH:" && second == "UPSTREAM_CST" ){
			to_return.second = 2;
			llstream >> to_return.first;
			break;
		}
	}
	return to_return;
}


