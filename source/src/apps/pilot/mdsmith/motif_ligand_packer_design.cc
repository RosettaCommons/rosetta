// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief mdsmith, sthyme
///

#include <devel/init.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/io/pdb/pose_io.hh> // pose_from_pdb
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/dna/setup.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
#include <protocols/motifs/LigandMotifSearch.hh>
#include <protocols/motifs/LigandMotifSearch.fwd.hh>

//#include <protocols/motifs/MotifLigandPacker.hh>
#include <protocols/dna/PDBOutput.hh>
#include <protocols/dna/util.hh>
//#include <protocols/dna/DnaInterfacePacker.hh>
//#include <protocols/dna/DnaInterfaceFinder.hh>
#include <protocols/motifs/motif_utils.hh>
using namespace protocols::dna;

#include <basic/prof.hh>
#include <basic/Tracer.hh>
static basic::Tracer TR("apps.pilot.motif_dna_packer_design");

#include <core/import_pose/import_pose.hh> // Need since refactor

// Utility Headers
#include <utility/file/file_sys_util.hh> // file_exists
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
using utility::vector1;
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
using utility::string_split;

// c++ headers
#include <fstream>
#include <iostream>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>

using namespace core;
using namespace basic;
using namespace chemical;
using namespace pack;
	using namespace task;
using namespace scoring;

////////////////////////////////////////////////////////////////////////////////
ScoreFunctionOP
get_total_rebuild_fullatom_scorefxn()
{
	ScoreFunctionOP scorefxn( new ScoreFunction() );
//	scorefxn->set_energy_method_options( scoring::methods::EnergyMethodOptions().exclude_DNA_DNA( false ) );
	//scorefxn->add_weights_from_database_file( "my_dna.wts" );
	//if ( option[ OK::dna::specificity::score_function ].user() ) {
	//	scorefxn->add_weights_from_file( option[ OK::dna::specificity::score_function ] );
	//} else {
	//	scorefxn->add_weights_from_database_file( "my_dna.wts" );
	//}

	// need to think more about this
	scorefxn->set_weight( chainbreak, 10.0 / 3.0 ); // option[ my_options::fullatom_dna_chainbreak_weight ] );
	//scorefxn->set_weight( dna_dihedral, 0.25 );//option[ my_options::fullatom_dna_dihedral_weight ] );
	//scorefxn->set_weight( dna_dihedral_sugar, 1.0 );//option[ my_options::fullatom_dna_dihedral_weight ] );

	/// for scoring ZN binding site:
	scorefxn->set_weight( atom_pair_constraint, 1.0 ); // approx?
	scorefxn->set_weight( res_type_constraint, 1.0 ); // approx?
	return scorefxn;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
    try {
	using namespace options;
	using namespace OptionKeys;

	devel::init( argc, argv ); // reading options--name should be more descriptive

	basic::prof_reset();

	TR << "Getting input filename(s)" << '\n';

	typedef vector1< utility::file::FileName > Filenames;
	Filenames pdbnames;

	if ( option[ in::file::l ].user() ) {
		Filenames listnames( option[ in::file::l ]().vector() );
		for ( Filenames::const_iterator filename( listnames.begin() );
		      filename != listnames.end(); ++filename ) {
			std::ifstream list( (*filename).name().c_str() );
			while ( list ) {
				std::string pdbname;
				list >> pdbname;
				pdbnames.push_back( pdbname );
			}
		}

	} else if ( option[ in::file::s ].user() ) {
		pdbnames = option[ in::file::s ]().vector();

	} else {
		std::cerr << "No files given: Use either -file:s or -file:l "
		          << "to designate a single pdb or a list of pdbs"
		          << std::endl;
	}

	bool const output_pdb( ! option[ OptionKeys::dna::design::nopdb ].user() );

	for ( Filenames::const_iterator filename( pdbnames.begin() );
	      filename != pdbnames.end(); ++filename ) {
		if ( !utility::file::file_exists( *filename ) ) {
			continue;
		}
		pose::Pose pose;
		//devel::blab::motif::motif_pose_from_pdb( pose, *filename, true );
		//devel::blab::viewer::add_conformation_viewer( pose.conformation(), pose.conformation() );

	//	io::pdb::pose_from_pdb( pose, *filename ); //Doesn't work anymore since refactor
		   core::import_pose::pose_from_pdb( pose, *filename );
		std::string pdbprefix( string_split( string_split( *filename, '/' ).back(), '.' ).front() );
		bool minimize( false );
		if ( option[ OptionKeys::out::prefix ].user() ) pdbprefix = option[ OptionKeys::out::prefix ]();
		if ( option[OptionKeys::dna::design::minimize].user() ) minimize = option[ OptionKeys::dna::design::minimize ]();

		// set up scoring
		std::string const weights( option[ OptionKeys::score::weights ]() );
		// this sets BasePartner in pose cacheable data
		scoring::dna::set_base_partner( pose );
		ScoreFunctionOP scorefxn( getScoreFunction() );

		//ScoreFunctionOP scorefxn( get_total_rebuild_fullatom_scorefxn() );
		//protocols::motifs::make_dna_mutations( pose );
		//scoring::dna::set_base_partner( pose );
	//	utility::vector1< Size > design_positions;
		// need 1 to 73 for 2p09
//		for ( Size pos_count( 1); pos_count<70; ++pos_count){
	//	design_positions.push_back(pos_count);}
	  // utility::vector1< Size  > design_positions;
		/*if ( option[OptionKeys::motifs::motif_build_positions].user()) {
	  utility::vector1< Size  > design_positions( option[ OptionKeys::motifs::motif_build_positions ]().vector() );
	for (  Size position( 1); position<=design_positions.size(); ++position){
				std::cout << "Design position: " << design_positions[position] << std::endl;
				}
	  } else {
		std::cout << "No build positions given, quitting for now" << std::endl;
		break;
	  } */
		protocols::motifs::LigandMotifSearchOP motif_search = new protocols::motifs::LigandMotifSearch;

		if (  option[ OptionKeys::motifs::motif_build_positions ].user() ) {
	  utility::vector1< Size  > design_positions( option[ OptionKeys::motifs::motif_build_positions ]().vector() );
	for (  Size position( 1); position<=design_positions.size(); ++position){
				std::cout << "about to run Design position: " << design_positions[position] << std::endl;
				}
		motif_search->run( pose, design_positions );
		} //End if user design positions
		else {
		// core::Real
	core::Real ligand_motif_sphere( basic::options::option[basic::options::OptionKeys::motifs::ligand_motif_sphere] ) ;

		std::cout << "Using ligand motif sphere of " <<  ligand_motif_sphere << " angstroms." << std::endl;
		motif_search->run( pose, ligand_motif_sphere );
		}
		//protocols::motifs::LigandMotifSearch ligand_search( scorefxn, minimize, pdbprefix );
		std::cout << "SUCCESSFUL COMPLETION" << std::endl;

	}
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
    return 0;
}

