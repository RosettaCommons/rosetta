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

/// @file apps/pilot/ronj/rss_energy_calculator.cc
/// @brief functions which iterates over PDBs, chops out fragments and scores them, and then outputs residue energies
/// @author Ron Jacak (ronj@email.unc.edu)

// Unit headers
#include <devel/init.hh>

//project Headers
#include <core/chemical/AA.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <basic/options/option.hh>

#include <protocols/jobdist/standard_mains.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <basic/prof.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

// Numeric Headers
#include <numeric/random/random.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>


// C++ headers
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>


static basic::Tracer TR("rss_energy_calculator");

using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace ObjexxFCL::format;

// application specific options
namespace rss_energy_calculator {
	BooleanOptionKey const no_repack_before_fragmenting( "rss_energy_calculator::no_repack_before_fragmenting" );
	BooleanOptionKey const no_repack_before_scoring_fragment( "rss_energy_calculator::no_repack_before_scoring_fragment" );
	IntegerOptionKey const fragment_length( "rss_energy_calculator::fragment_length" );
}


std::string usage_string;

//@brief pretty prints a usage prompt if the user doesn't specify any options
void init_usage_prompt( std::string exe ) {

	// place the prompt up here so that it gets updated easily; global this way, but that's ok
	std::stringstream usage_stream;
	usage_stream << "No files given: Use either -file:s or -file:l to designate a single pdb or a list of pdbs.\n\n"
			<< "Usage: " << exe
			<< "\n\t-database path/to/minidb"
			<< "\n\t-s pdb|-l pdbs"
			<< "\n\t-ignore_unrecognized_res"
			<< "\n\t[-ex1 -ex2 -ex3 -ex4]"
			<< "\n\t[-use_input_sc]"
			<< "\n\t[-linmem_ig 10]"
			<< "\n\t[-no_his_his_pairE]"
			<< "\n\t[-mute core.io core.conformation]"
			<< "\n\t[-no_repack_before_fragmenting]"
			<< "\n\t[-no_repack_before_scoring_fragment]"
			<< "\n\n";
	usage_string = usage_stream.str();

}


///@brief Takes the input Pose and runs a fast repack protocol on it.
void repack_pose( pose::Pose & pose, scoring::ScoreFunctionOP scorefxn ) {

	pack::task::PackerTaskOP repack_task = pack::task::TaskFactory::create_packer_task( pose );
	repack_task->set_bump_check( true );
	repack_task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );

	protocols::simple_moves::PackRotamersMoverOP repack_protocol( new protocols::simple_moves::PackRotamersMover( scorefxn, repack_task, 1 /*ndruns value hardcoded*/ ) );
	repack_protocol->apply( pose );

	return;
}


//@brief the main workhorse of this protocol, it creates fragments from the passed in pdb file and stores things in the passed-in map
void
create_and_score_fragments( std::string pdb_filename, scoring::ScoreFunctionOP scorefxn,
	std::map< std::string, std::vector< std::pair< std::string, scoring::EnergyMap > > > & aa_to_vector_pair_sequence_energymap ) {

	TR << "create_and_score_fragments(): creating pose object from file " << pdb_filename << std::endl;

	pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, pdb_filename );

	if ( !( option[ rss_energy_calculator::no_repack_before_fragmenting ] ) ) {
		repack_pose( pose, scorefxn );
	}

	// need to get a count of the number of residues in the protein since n_residue (or conformation.size()) includes metal atoms
	// and other ligands as residues now. really, this should be a method in the Pose object.
	int countProteinResidues = 0;
	Size chain = pose.residue(1).chain();
	for ( Size i=1; i < pose.n_residue(); ++i ) {
		if ( pose.residue(i).type().is_protein() )
			countProteinResidues++;
		if ( chain != pose.residue(i).chain() )
			break;
	}

	// the questions we have to worry about here are how many fragments to create from a given PDB file and how to select those fragments.
	// In this first implementation, select total_residue in sequence DIV 10 as the number of fragments to create and vary
	// the fragment start locations using a MOD obtained value to select "random" fragments.
	int fragment_length;
	if ( option[ rss_energy_calculator::fragment_length ].user() )
		fragment_length = option[ rss_energy_calculator::fragment_length ];
	else
		fragment_length = 13;

	int fragments_to_create = countProteinResidues / 10;  // bigger proteins, more fragments
	if ( fragments_to_create == 0 ) { return; }       // may occur with "proteins" of size < 10

	//TR << "create_and_score_fragments(): creating " << fragments_to_create << " fragments for pdb '" << pdb_filename << "'" << std::endl;

	for ( int i=1; i <= fragments_to_create; ++i ) {

		// generate a random start location within this structure; use a random number whose range does not include the termini.
		// For example, if we want to pick a 13-mer, we can't do that at residue 57 of a 60-residue protein.  Instead, have
		// the random number be selected from between 1 and (60-13+1)=48.  That way, if 48 comes up, we get just the last 13
		// residues.  If 1 comes up, we get the first 13 residues.  Otherwise, we get some chunk in the middle of the protein.
		int random_seqpos = numeric::random::random_range( 1, countProteinResidues - fragment_length + 1 );
		if ( random_seqpos <= 0 || random_seqpos >= (int)countProteinResidues ) {
			utility_exit_with_message("Something wrong with random.");
		}

		int fragment_start = random_seqpos;
		int fragment_end = random_seqpos + fragment_length;

		// check to make sure the start and end are on the same chain; if not, regenerate the start. also check to make sure
		// both residues are of type protein.  it's still possible that we have some nonresidue types in between the ends
		// but that's a risk I'm willing to let happen.
		if ( pose.chain( fragment_start ) != pose.chain( fragment_end ) ||
				!pose.residue(fragment_start).type().is_protein() || !pose.residue(fragment_end).type().is_protein() ) {

			do {
				random_seqpos = numeric::random::random_range( 1, countProteinResidues - fragment_length + 1 );
				fragment_start = random_seqpos;
				fragment_end = random_seqpos + fragment_length;
			} while ( pose.chain( fragment_start ) != pose.chain( fragment_end ) ||
				!pose.residue(fragment_start).type().is_protein() || !pose.residue(fragment_end).type().is_protein() );

		}

		TR << "create_and_score_fragments(): considering fragment " << fragment_start << "-" << fragment_start + fragment_length << std::endl;

		/*
		// there's a situation that can come up with disulfide bonds; if a fragment has two cys near each other in primary sequence
		// they may be forming a disulfide.  if that's not recognized as such, it can lead mini to assign very high repulsive
		// energies for that residue.  hopefully, this situation occurs so rarely that it won't affect the average energies.
		// to properly deal with the issue, Doug Renfrew started on the code below to remove disulfides if they're present.

		// if 2 residues in the segment have a disulfide bond keep the bond otherwise remove it
		Energy score_orig = ( *scorefxn )( pose ); // not sure if nessesary
		scoring::disulfides::DisulfideEnergyContainerCOP dec = new scoring::disulfides::DisulfideEnergyContainer( pose );

		ResidueTypeSetCAP RTS( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		ResidueType const & RT( RTS->name_map( "CYS" ) );

		for (Size i = rand_lower; i <= rand_upper; ++i ) {
			if( dec->residue_forms_disulfide(i) ) {
				Size partner( dec->other_neighbor_id(i) );
				std::cout << "DISULF: RANGE " << rand_lower << " " << rand_upper << " " << pdb_filename << " " << i << " " << pose.residue(i).name3() << " and " << partner << " " << pose.residue(i).name3() << " in " << pdb_filename << std::endl;
				if ( partner > rand_upper || partner < rand_lower ) {
					//core::pose::remove_variant_type_from_pose_residue( pose, chemical::DISULFIDE, i );
					//core::pose::remove_variant_type_from_pose_residue( pose, chemical::DISULFIDE, partner );
					pose::replace_pose_residue_copying_existing_coordinates( pose, i, RT );
					pose::replace_pose_residue_copying_existing_coordinates( pose, partner, RT );
					std::cout << "DISULF: REMOVED varient type DISULPHIDE from residue " << i << " " << pose.residue(i).name3() << " and " << partner << " " << pose.residue(i).name3() << " in " << pdb_filename << std::endl;
				} else {
					std::cout << "DISULF: PRESERVED varient type DISULPHIDE from residue " << i << " " << pose.residue(i).name3() << " and " << partner << " " << pose.residue(i).name3() << " in " << pdb_filename << std::endl;
				}
			}
		}
		*/

		// create fragment pose and seed it with the first res by jump, then add the rest of the fragment
		pose::Pose fragment;
		fragment.append_residue_by_jump( pose.residue( fragment_start ), 1 );
		for (int jj = fragment_start + 1; jj < fragment_start + fragment_length; ++jj ) {
			fragment.append_residue_by_bond( pose.residue(jj) );
		}

		//TR << "create_and_score_fragments(): scoring fragment pose object" << std::endl;
		if ( !( option[ rss_energy_calculator::no_repack_before_scoring_fragment ] ) ) {
			repack_pose( fragment, scorefxn );
		}

		// score fragment
		Energy score = ( *scorefxn )( fragment );
		if ( score > 200 ) {
			TR << "Very high energy found for fragment '" << fragment.sequence() << "'; score: " << score << std::endl;

			//scorefxn->show( std::cout, fragment );
			//std::cout << std::endl;

			utility::file::FileName fn( pdb_filename );
			std::stringstream outputfilename;
			outputfilename << fn.base() << "(" << fragment_start << "-" << (fragment_start+fragment_length) << ").pdb";
			fragment.dump_scored_pdb( outputfilename.str() , *scorefxn );

			// utility_exit_with_message("Bad energy calculation.");
		}

		int central_residue_index = (fragment_length / 2) + 1;
		std::string central_residue_aa = fragment.residue( central_residue_index ).name3();
		std::string fragment_sequence = fragment.sequence();

		// add up the energies of the individual terms for the central residue
		scoring::EnergyMap const & central_residue_total_energies( fragment.energies().residue_total_energies( central_residue_index ) );
		/*if ( i <= 1 ) {
			utility::file::FileName fn( pdb_filename );
			std::stringstream outputfilename;
			outputfilename << fn.base() << "(" << fragment_start << "-" << (fragment_start+fragment_length) << ").pdb";
			fragment.dump_scored_pdb( outputfilename.str() , *scorefxn );
		}*/

		aa_to_vector_pair_sequence_energymap[ central_residue_aa ].push_back(
			std::pair< std::string, scoring::EnergyMap >( fragment_sequence, central_residue_total_energies ) );


	} // end loop fragments to create

} // end create_and_score_fragments()


//@brief an alternative way of calculating unfolded state energies. calculates the energy of every amino acid in all of the structures
/// and stores them in the passed-in map
void
score_folded_residues( std::string pdb_filename, scoring::ScoreFunctionOP scorefxn,
	std::map< std::string, std::vector< std::pair< std::string, scoring::EnergyMap > > > & aa_to_vector_pair_sequence_energymap ) {

	TR << "score_folded_residues(): creating pose object from file " << pdb_filename << std::endl;

	pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, pdb_filename );

	if ( !( option[ rss_energy_calculator::no_repack_before_fragmenting ] ) ) {
		TR << "score_folded_residues(): repacking pose..." << std::endl;
		repack_pose( pose, scorefxn );
	}

	// score the entire pose
	Energy score = ( *scorefxn )( pose );
	if ( score > 1000 ) {
		TR << "Very high energy found for pose '" << pose.sequence() << "'; score: " << score << std::endl;

		//scorefxn->show( std::cout, fragment );
		//std::cout << std::endl;

		utility::file::FileName fn( pdb_filename );
		std::stringstream outputfilename;
		outputfilename << fn.base() << ".scored.pdb";
		pose.dump_scored_pdb( outputfilename.str() , *scorefxn );

		// utility_exit_with_message("Bad energy calculation.");
	}

	for ( Size ii=1; ii <= pose.n_residue(); ++ii ) {
		// check to make sure this residue is of type protein.
		if ( !pose.residue(ii).type().is_protein() ) { continue; }

		std::string residue_aa = pose.residue(ii).name3();
		std::string sequence = "PROTEIN";

		// save the energies of the individual terms for this residue
		scoring::EnergyMap const & residue_total_energies( pose.energies().residue_total_energies( ii ) );

		aa_to_vector_pair_sequence_energymap[ residue_aa ].push_back( std::pair< std::string, scoring::EnergyMap >( sequence, residue_total_energies ) );
	}

} // end score_folded_residues()


//@brief main method for the solubility protocols
int
main( int argc, char * argv [] ) {

	try {


	// add application specific options to options system
	option.add( rss_energy_calculator::no_repack_before_fragmenting, "Do not repack the pose to remove any steric clashes before making fragments." );
	option.add( rss_energy_calculator::no_repack_before_scoring_fragment, "Do not repack the fragment before scoring it." );
	option.add( rss_energy_calculator::fragment_length, "Length of fragment to use for energies." );

	// initialize
	devel::init(argc, argv);

	//
	// concatenate -s and -l flags together to get total list of PDB files
	// The advantage of parsing -s and -l separately is that users can specify a list and a single structure on the
	// command line.  Not anymore.  To get an output file name that's based on the input structures, I have to use
	// either the -l or -s option.
	//
	using utility::file::file_exists;

	std::vector< utility::file::FileName > pdb_file_names;

	if ( option[ in::file::s ].user() )
		pdb_file_names = option[ in::file::s ]().vector(); // make a copy (-s)
	else {
		std::vector< utility::file::FileName > list_file_names;
		if ( option[ in::file::l ].user() ) {
			list_file_names = option[ in::file::l ]().vector(); // make a copy (-l)

			//ronj for each file input with the -l switch...
			for ( std::vector< utility::file::FileName >::iterator i = list_file_names.begin(), i_end = list_file_names.end(); i != i_end; ++i ) {
				std::string listfilename( i->name() );
				std::ifstream data( listfilename.c_str() );
				//ronj try to open that particular file first...
				if ( !data.good() ) {
					utility_exit_with_message( "Unable to open file: '" + listfilename + "'\n" );
				}
				std::string line;
				//ronj then read all the lines in that file until there are no more
				while( getline(data, line) ) {
					pdb_file_names.push_back( utility::file::FileName(line) );
				}
				data.close();
			}
		}
	}

	if ( pdb_file_names.size() == 0 ) {
		init_usage_prompt( argv[0] );
		utility_exit_with_message_status( usage_string, 1 );
	}


	// create a custom score function; basically we want the given score function without the reference energies
	// setting a ScoreType's (or is it an EnergyMethod) weight to 0 calls remove_method in the ScoreFunction, so it does remove
	// that term from the score function.
	// Don't include the solubility term or the pAA term that we are considering using here.
	TR << "Creating score function." << std::endl;
	scoring::ScoreFunctionOP scorefxn = scoring::getScoreFunction();

	scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn->energy_method_options() );
	energymethodoptions.hbond_options()->decompose_bb_hb_into_pair_energies( true );
	scorefxn->set_energy_method_options( energymethodoptions );

	scorefxn->set_weight( core::scoring::ref, 0.0 );

	scorefxn->set_weight( core::scoring::dslf_ss_dst, 0.0 );
	scorefxn->set_weight( core::scoring::dslf_cs_ang, 0.0 );
	scorefxn->set_weight( core::scoring::dslf_ss_dih, 0.0 );
	scorefxn->set_weight( core::scoring::dslf_ca_dih, 0.0 );

	// I want sequence information and energy breakdowns for each fragment.
	// The sequence information will be easy to get from the Pose object and the energies should just be kept in an EnergyMap.
	// Since there's no inherent order, store all the fragment energies in a vector.  But let's at least keep everything organized
	// using a Map that translates AA type to a vector of EnergyMaps.
	// map_aa_to_vector_pair_sequence_energymap[ D ][ 0 ] = ( PTDSGADSAIVMF, energies )
	// map_aa_to_vector_pair_sequence_energymap[ D ][ 1 ] = ( TSRSGADSAIVMF, energies )
	// map_aa_to_vector_pair_sequence_energymap[ W ][ 0 ] = ( AGDSGAWSAIPTD, energies )
	// map_aa_to_vector_pair_sequence_energymap[ S ][ 0 ] = ( IYSAIDSGASFMV, energies )

	std::map< std::string, std::vector< std::pair< std::string, scoring::EnergyMap > > > map_aa_to_vector_pair_sequence_energymap;

	// iterate through all the structures - do something to them
	for ( std::vector< utility::file::FileName >::iterator pdb = pdb_file_names.begin(), last_pdb = pdb_file_names.end(); pdb != last_pdb; ++pdb) {
		if ( !file_exists( *pdb ) ) {
			std::cerr << "Input pdb " << *pdb << " not found, skipping" << std::endl;
			continue;
		}
		//TR << "Fragmenting and scoring pdb '" << pdb->name() << "'" << std::endl;
		//create_and_score_fragments( pdb->name(), scorefxn, map_aa_to_vector_pair_sequence_energymap );

		TR << "Scoring residues within pdb '" << pdb->name() << "'" << std::endl;
		score_folded_residues( pdb->name(), scorefxn, map_aa_to_vector_pair_sequence_energymap );
	}


	// I want the output to look something like the following, since I already have a Perl script that can aggregate things for me.
	// Though, it probably wouldn't be hard to write the aggregation logic into this app.

	// aa1  aa3  frag_seq  Eatr  Erep  Eintra_rep  Esol  Epair Ehbondx4 Erama Eomega Edun Epaa_pp Etotal
	// P PRO PTDSGAPSAIVMF -1.643  0.349 -1.643  0.453  0.00000 -1.132  0.00000  0.767  0.00000 -1.18486
	// M MET APSAIVMFPVGEK -1.199  0.001 -1.199  0.463  0.00000  0.054  0.00000  1.500  0.00024  0.48003
	// K LYS GEKPNPKGAAMKP -1.069  0.078 -1.069  0.560  0.00000 -0.103  0.00000  3.213  0.00000  2.02869
	// F PHE AMKPVVFNHLIHE -3.142  0.157 -3.142  1.068  0.00000 -0.153  0.00000  0.830  0.00049 -0.60859

	//std::string central_residue_aa = fragment.residue( central_residue_index ).name3();
	//std::string fragment_sequence = fragment.sequence();

	// add up the energies of the individual terms for the central residue
	//EnergyMap const & central_residue_total_energies( fragment.energies().residue_total_energies( central_residue_index ) );

	TR << "Done scoring. Outputting values from map." << std::endl;

	// instead of placing the name of the term before the value, place them all at the top as a header line
	//ALA ETVLAEATVAPVS fa_atr: -0.821 fa_rep:  0.120 fa_sol:  0.648 fa_intra_rep:  0.159 fa_pair:  0.000 hbond_sr_bb:  0.000 hbond_lr_bb:  0.000 hbond_bb_sc:  0.000 hbond_sc:  0.000 rama:  0.653 omega:  0.095 fa_dun:  0.000 p_aa_pp:  0.515 ref:  0.000 total_score:  0.160
	// aa3 frag_seq fa_atr	fa_rep	fa_sol	fa_intra_rep	fa_pair	hbond_sr_bb	hbond_lr_bb hbond_bb_sc	hbond_sc	rama	omega	fa_dun	p_aa_pp total_score
	//ALA TPKDDEAWCATCH -2.871 0.624 2.074 0.226 0.000 0.000 0.000 0.000 0.000 -1.280 0.433 0.000 -0.325 -0.817

	utility::file::FileName outfn;
	if ( option[ in::file::s ].user() )
		outfn = pdb_file_names[0].name();
	else
		outfn = ( option[ in::file::l ]().vector() )[0];

	std::stringstream outputfilename;
	outputfilename << outfn.bare_name() << ".fragmentenergies.txt";
	std::ofstream energiesOutFile( outputfilename.str().c_str() );


	// iterate over the map and divide real by count and print shit out
	energiesOutFile << "aa frag_seq ";
	for ( int kk = 1; kk <= scoring::n_score_types; ++kk ) {
		switch( scoring::ScoreType(kk) ) {
		case scoring::fa_atr:
		case scoring::fa_rep:
		case scoring::fa_intra_rep:
		case scoring::fa_sol:
		case scoring::fa_pair:
		case scoring::hbond_sr_bb:
		case scoring::hbond_lr_bb:
		case scoring::hbond_bb_sc:
		case scoring::hbond_sc:
		case scoring::rama:
		case scoring::omega:
		case scoring::fa_dun:
		case scoring::p_aa_pp:
		case scoring::pro_close:
			energiesOutFile << scoring::ScoreType(kk) << " ";
			break;
		default:
			break;
		}
	}
	energiesOutFile << std::endl;

	std::map< std::string, std::vector< std::pair< std::string, scoring::EnergyMap > > >::const_iterator mi;
	for( mi = map_aa_to_vector_pair_sequence_energymap.begin(); mi != map_aa_to_vector_pair_sequence_energymap.end(); ++mi ) {

		for ( Size jj=0; jj < (mi->second).size(); ++jj ) {

			energiesOutFile << mi->first << " "; // prints out the 3-letter code
			energiesOutFile << mi->second[ jj ].first << " "; // prints out the sequence

			// now print out the energies, but only those we care about
			for ( int kk = 1; kk <= scoring::n_score_types; ++kk ) {
				switch( scoring::ScoreType(kk) ) {
				case scoring::fa_atr:
				case scoring::fa_rep:
				case scoring::fa_intra_rep:
				case scoring::fa_sol:
				case scoring::fa_pair:
				case scoring::hbond_sr_bb:
				case scoring::hbond_lr_bb:
				case scoring::hbond_bb_sc:
				case scoring::hbond_sc:
				case scoring::rama:
				case scoring::omega:
				case scoring::fa_dun:
				case scoring::p_aa_pp:
				case scoring::pro_close:
					energiesOutFile << ObjexxFCL::format::F(7,5,mi->second[ jj ].second[ scoring::ScoreType(kk) ]) << " ";
					break;
				default:
					break;
				}
			}

			energiesOutFile << std::endl;
		}
	}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

