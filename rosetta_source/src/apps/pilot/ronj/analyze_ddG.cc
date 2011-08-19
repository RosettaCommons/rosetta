// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file analyze_ddG_stability.cc
/// @brief A protocol which takes a mutation and computes the ddG (stability) for the mutation
/// @author Ron Jacak


// Unit headers
#include <devel/init.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/conformation/Residue.hh>
#include <protocols/moves/MinMover.hh>
#include <protocols/moves/SetReturningPackRotamersMover.hh>
#include <protocols/moves/TaskAwareMinMover.hh>
#include <protocols/toolbox/pose_metric_calculators/SurfaceCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.hh>

#include <basic/options/util.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

// Numeric Headers

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <sstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>


#undef NDEBUG

static basic::Tracer TR("analyze_ddG_stability");

using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace ObjexxFCL::fmt;


// application specific options
namespace analyze_ddG_stability {
	BooleanOptionKey const use_reweighted_score12_with_surfaceE_scorefunction( "analyze_ddG_stability::use_reweighted_score12_with_surfaceE_scorefunction" );
	BooleanOptionKey const use_reweighted_score12_scorefunction( "analyze_ddG_stability::use_reweighted_score12_scorefunction" );
	BooleanOptionKey const calculate_surface_energy( "analyze_ddG_stability::calculate_surface_energy" );
	StringOptionKey const mutation( "analyze_ddG_stability::mutation" );
}

std::string usage_string;


///
/// @begin init_usage_prompt
///
/// @brief
/// the usage prompt that gets printed when the user doesn't enter all the required command line arguments
///
void init_usage_prompt( std::string exe ) {

	// place the prompt up here so that it gets updated easily; global this way, but that's ok
	std::stringstream usage_stream;
	usage_stream
			<< "Usage: " << exe
			<< "\n\t-database path/to/minidb"
			<< "\n\t-s pdb"
			<< "\n\t-mutation ALA,16,TRP,A"

			<< "\n\t[-use_reweighted_score12_with_surfaceE_scorefunction]"
			<< "\n\t[-use_reweighted_score12_scorefunction]"
			<< "\n\t[-calculate_surface_energy]"

			<< "\n\t[-ex1 [-ex2]]"
			<< "\n\t[-ndruns 5]"

			<< "\n\t[-ignore_unrecognized_res]"
			<< "\n\t[-mute core.io core.conformation core.pack core.scoring]"

			<< "\n\n";
	usage_string = usage_stream.str();

}

///
/// @begin print_energies
///
/// @brief
/// Helper method for the main function. Takes in a pose, the scorefunction, surface energies and weights and prints everything
/// out in a pretty format.
///
void print_energies( scoring::EnergyMap e, scoring::ScoreFunctionOP scorefxn ) {

	scoring::EnergyMap const & wts( scorefxn->weights() );

	for ( int jj = 1; jj <= scoring::n_score_types; ++jj ) {
		Real const weight = wts[ scoring::ScoreType(jj) ];

		switch( scoring::ScoreType( jj ) ) {
		case scoring::fa_atr:
		case scoring::fa_rep:
		case scoring::fa_sol:
		case scoring::fa_pair:
		case scoring::hbond_lr_bb:
		case scoring::hbond_sr_bb:
		case scoring::hbond_bb_sc:
		case scoring::hbond_sc:
		case scoring::rama:
		case scoring::omega:
		case scoring::fa_dun:
		case scoring::p_aa_pp:
		case scoring::ref:
		case scoring::surface:
			std::cout << scoring::ScoreType(jj) << ": " << ObjexxFCL::fmt::F(5,2, weight * e[ scoring::ScoreType(jj) ] ) << ", ";
			break;
		default:
			break;
		}
	}
	std::cout << std::endl;

}


void tokenize_string( const std::string & str, std::vector< std::string > & tokens, const std::string & delimiters = " " ) {

	// Skip delimiters at beginning.
	std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);

	// Find first "non-delimiter".
	std::string::size_type pos = str.find_first_of(delimiters, lastPos);

	while ( std::string::npos != pos || std::string::npos != lastPos ) {
		// Found a token, add it to the vector.
		tokens.push_back( str.substr( lastPos, pos - lastPos ) );

		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of( delimiters, pos );

		// Find next "non-delimiter"
		pos = str.find_first_of( delimiters, lastPos );
	}
}

///
/// @begin main
///
/// @brief main method for the ddG protocol
///
int main( int argc, char* argv[] ) {

	//
	// add application specific options to options system
	//
	option.add( analyze_ddG_stability::use_reweighted_score12_with_surfaceE_scorefunction, "Use a score12 scorefunction that has the surface term and reference energies, with everything reweighted." );
	option.add( analyze_ddG_stability::use_reweighted_score12_scorefunction, "Use a reweighted score12 scorefunction, without the surface energy." );
	option.add( analyze_ddG_stability::mutation, "The mutation to make to the wild-type pose. Format to use is old-residue-type,PDB-res-num,new-residue-type." );
	option.add( analyze_ddG_stability::calculate_surface_energy, "Calculate the surface energies of the wild-type and mutant poses." );


	//
	// options, random initialization
	//
	devel::init( argc, argv );

	//
	// Only use the -s flag. This protocol is not going to handle list input files for structures.  Also check to make sure the file exists.
	//
	using utility::file::file_exists;

	utility::file::FileName pdb_file_name;
	if ( option[ in::file::s ].active() ) {
		//pdb_file_name = utility::file::FileName( option[ in::file::s ]()[0] );
		pdb_file_name = utility::file::FileName( basic::options::start_file() );


		// check to make sure the file exist; if not, move on to the next one
		if ( !file_exists( pdb_file_name ) ) {
			std::cerr << "Input pdb " << pdb_file_name.name() << " not found. Quitting protocol." << std::endl;
			exit(1);
		}

	} else {
		init_usage_prompt( argv[0] );
		std::cout << usage_string;
		utility_exit_with_message_status( "No files given: Use -file:s to designate a single pdb for making mutations.\n", 1 );
	}

	//
	// Error out if the user didn't specify a mutation; then parse the mutation string
	//
	if ( ! option[ analyze_ddG_stability::mutation ].active() ) {
		init_usage_prompt( argv[0] );
		std::cout << usage_string;
		utility_exit_with_message_status( "No mutation string specified.\n", 1 );
	}

	std::vector< std::string > parsed_tokens;
	std::string mutation_string = option[ analyze_ddG_stability::mutation ];

	tokenize_string( mutation_string, parsed_tokens, "," );

	std::string wt_residue_type;
	std::string mutant_residue_type;
	Size resid;
	char chain;

	wt_residue_type = parsed_tokens[0];
	mutant_residue_type = parsed_tokens[2];

	std::istringstream ss( parsed_tokens[1] );
	ss >> resid;

	// converting a string to a char is easy in C++
	chain = parsed_tokens[3][0];

	#ifndef NDEBUG
		TR << "Finished parsing mutation string. Found mutation of '" << wt_residue_type << "' to '"
			<< mutant_residue_type << "' at PDB resid: " << resid << " on chain '" << chain << "'" << std::endl;
	#endif


	//
	// create a custom score function as well as a score12 scorefunction
	//
	scoring::ScoreFunctionOP scorefxn;
	if ( option[ analyze_ddG_stability::use_reweighted_score12_with_surfaceE_scorefunction ] ) {
		TR << "Using reweighted score12 with surfaceE scorefunction." << std::endl;
		scorefxn = scoring::ScoreFunctionFactory::create_score_function( scoring::STANDARD_WTS, scoring::SCORE12_PATCH );

		utility::vector1< Real > refEs_;
		refEs_.resize( chemical::num_canonical_aas, 0.0 );
		Real const rpp_refs[20] = { 1.01553, 2.33569, 0.918812, 1.1572, 1.92699, 1.02994, 1.72094, 0.824097, 1.07864, 0.569007,
				1.19491, 0.602461, 0.883383, 0.919596, 0.878134, 0.776869, 0.596056, 0.517296, 3.53647, 1.30933 };

		for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
			refEs_[ii] = rpp_refs[ ii-1 ];
		}
		scorefxn->set_method_weights( core::scoring::ref, refEs_ );

		scorefxn->set_weight( core::scoring::fa_atr, 0.6 );
		scorefxn->set_weight( core::scoring::fa_rep, 0.807902 );
		scorefxn->set_weight( core::scoring::fa_sol, 0.458364 );
		scorefxn->set_weight( core::scoring::fa_intra_rep, 0.0588026 );
		scorefxn->set_weight( core::scoring::pro_close, 0.0690371 );
		scorefxn->set_weight( core::scoring::fa_pair, 0.127158 );
		scorefxn->set_weight( core::scoring::hbond_sr_bb, 1.7428 );
		scorefxn->set_weight( core::scoring::hbond_lr_bb, 1.76154 );
		scorefxn->set_weight( core::scoring::hbond_bb_sc, 0.816982 );
		scorefxn->set_weight( core::scoring::hbond_sc, 0.393182 );
		scorefxn->set_weight( core::scoring::dslf_ss_dst, 1 );
		scorefxn->set_weight( core::scoring::dslf_cs_ang, 1 );
		scorefxn->set_weight( core::scoring::dslf_ss_dih, 1 );
		scorefxn->set_weight( core::scoring::dslf_ca_dih, 1 );
		scorefxn->set_weight( core::scoring::rama, 0.0134665 );
		scorefxn->set_weight( core::scoring::omega, 0.5 );
		scorefxn->set_weight( core::scoring::fa_dun, 0.0603081 );
		scorefxn->set_weight( core::scoring::p_aa_pp, 0.351527 );
		scorefxn->set_weight( core::scoring::ref, 1.0 );
		scorefxn->set_weight( core::scoring::surface, 1.0 );

	} else if ( option[ analyze_ddG_stability::use_reweighted_score12_scorefunction ] ) {

		TR << "Using reweighted score12 scorefunction." << std::endl;
		scorefxn = scoring::ScoreFunctionFactory::create_score_function( scoring::STANDARD_WTS, scoring::SCORE12_PATCH );

		utility::vector1< Real > refEs_;
		refEs_.resize( chemical::num_canonical_aas, 0.0 );
		Real const rpp_refs[20] = { 1.44647, 3.20455, 1.28033, 1.50335, 2.39638, 1.47398, 2.35098, 1.24274, 1.64901, 1.11337,
				1.68414, 0.895638, 1.68502, 1.43777, 1.28772, 1.26944, 1.03373, 0.903269, 4.10334, 1.8192 };
		for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
			refEs_[ii] = rpp_refs[ ii-1 ];
		}
		scorefxn->set_method_weights( core::scoring::ref, refEs_ );

		scorefxn->set_weight( core::scoring::fa_atr, 0.6 );
		scorefxn->set_weight( core::scoring::fa_rep, 0.676039 );
		scorefxn->set_weight( core::scoring::fa_sol, 0.473486 );
		scorefxn->set_weight( core::scoring::fa_intra_rep, 0.0669174 );
		scorefxn->set_weight( core::scoring::pro_close, 0.0430282 );
		scorefxn->set_weight( core::scoring::fa_pair, 0.155594 );
		scorefxn->set_weight( core::scoring::hbond_sr_bb, 1.34807 );
		scorefxn->set_weight( core::scoring::hbond_lr_bb, 1.79619 );
		scorefxn->set_weight( core::scoring::hbond_bb_sc, 0.855615 );
		scorefxn->set_weight( core::scoring::hbond_sc, 0.34485 );
		scorefxn->set_weight( core::scoring::dslf_ss_dst, 1 );
		scorefxn->set_weight( core::scoring::dslf_cs_ang, 1 );
		scorefxn->set_weight( core::scoring::dslf_ss_dih, 1 );
		scorefxn->set_weight( core::scoring::dslf_ca_dih, 1 );
		scorefxn->set_weight( core::scoring::rama, 0.0407331 );
		scorefxn->set_weight( core::scoring::omega, 0.5 );
		scorefxn->set_weight( core::scoring::fa_dun, 0.0647276 );
		scorefxn->set_weight( core::scoring::p_aa_pp, 0.319602 );
		scorefxn->set_weight( core::scoring::ref, 1.0 );

	} else {

		TR << "Using standard score12 scorefunction." << std::endl;
		scorefxn = scoring::ScoreFunctionFactory::create_score_function( scoring::STANDARD_WTS, scoring::SCORE12_PATCH );

	}

	//EM options for bb-bb hbond scores
 	scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn->energy_method_options() );
 	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
 	scorefxn->set_energy_method_options( energymethodoptions );

	//
	// repack and score the wild type structure
	// scoring isn't really necessary here, but scoring is fast and poses seem to be happier once they've been scored
	//
	pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, pdb_file_name.name() );

	// score the wt pose
	(*scorefxn)( pose );

	//
	// create a surface metric calculator, if using surface energy
	//
	pose::metrics::PoseMetricCalculatorOP surface_calculator;

	basic::MetricValue< Real > wt_pose_surface_total;
	basic::MetricValue< utility::vector1< core::Real > > wt_pose_residue_surface;

	if ( option[ analyze_ddG_stability::calculate_surface_energy ] ) {
		surface_calculator = new protocols::toolbox::pose_metric_calculators::SurfaceCalculator;
		pose::metrics::CalculatorFactory::Instance().register_calculator( "surface", surface_calculator );

		pose.metric( "surface", "residue_surface", wt_pose_residue_surface );
	}


	//
	// now make the mutation pose and repack the designed position as well as neighboring positions
	//

	using namespace core::pack::task::operation;
	using namespace core::pack::task;
	using namespace basic::options;

	pose::Pose mutant_pose = pose;

	// convert the PDB resid to the pose resid in case the PDB numbering is screwy
	Size mutant_pose_resid = pose.pdb_info()->pdb2pose().find(chain, resid); // pdb info doesn't get copied so use the original pose

	#ifndef NDEBUG
		TR << "Converted PDB res num " << resid << " to pose res num " << mutant_pose_resid << std::endl;
	#endif


	// need to create a calculator here that I can use to identify neighbors of the mutated residue
	pose::metrics::PoseMetricCalculatorOP mutant_nb_calculator = new toolbox::pose_metric_calculators::NeighborsByDistanceCalculator( mutant_pose_resid );
	pose::metrics::CalculatorFactory::Instance().register_calculator( "mutant_nb_calculator", mutant_nb_calculator );


	// the restrict operation class (which in the end is just a TaskOperation) takes a calculator during construction. I've already
	// created that calculator above.  This operation will disable repacking and design at all positions except those in the neighborhood
	// of the mutated position.
	TaskOperationCOP nb_op = new toolbox::task_operations::RestrictToNeighborhoodOperation( "mutant_nb_calculator" );

	// extra operations we want to also include
	// the restrict residue to repacking ops are used to make sure that only repacking and not design is done to the residues in the neighborhood
	InitializeFromCommandlineOP init_op = new InitializeFromCommandline();
	IncludeCurrentOP ic_op = new IncludeCurrent();
	RestrictResidueToRepackingOP mutant_repack_op = new RestrictResidueToRepacking();
	RestrictResidueToRepackingOP wt_repack_op = new RestrictResidueToRepacking(); // will include one extra residue to repack
	// not sure how to use bump check with Task Operations
	// TODO: read resfile operation

	TaskFactoryOP wt_tf = new TaskFactory();
	TaskFactoryOP mutant_tf = new TaskFactory();

	wt_tf->push_back( init_op ); mutant_tf->push_back( init_op );
	wt_tf->push_back( ic_op ); mutant_tf->push_back( ic_op );
	wt_tf->push_back( nb_op ); mutant_tf->push_back( nb_op );

	bool found_mutant_residue = false;

	for ( Size ii = 1; ii <= mutant_pose.n_residue(); ++ii ) {
		if ( ii == mutant_pose_resid ) {

			if ( mutant_pose.residue( mutant_pose_resid ).aa() != core::chemical::aa_from_name( wt_residue_type ) ) {
				TR << "Wild-type residue differs from what was entered. wt: " << name_from_aa( mutant_pose.residue( mutant_pose_resid ).aa() ) << std::endl;
				utility_exit_with_message( "Incorrect residue type found.\n" );
			}
			found_mutant_residue = true;

			// do design on this position
			utility::vector1< bool > keep_canonical_aas( chemical::num_canonical_aas, false );
			keep_canonical_aas[ core::chemical::aa_from_name( mutant_residue_type ) ] = true;
			RestrictAbsentCanonicalAASOP restrict_op = new RestrictAbsentCanonicalAAS( ii, keep_canonical_aas );
			mutant_tf->push_back( restrict_op );

			// for the wild type, don't design on the mutant resid - but do allow repacking
			wt_repack_op->include_residue( ii );

		} else {
			// make this position repackable only; because of the commutativity of packer task ops, only the residues that are in the neighborhood
			// of the mutant will be allowed to repack. the restrict to neighborhood op will disallow packing at all positions not near the mutant.
			mutant_repack_op->include_residue( ii );
			wt_repack_op->include_residue( ii );
		}
	}

	wt_tf->push_back( wt_repack_op );
	mutant_tf->push_back( mutant_repack_op );


	// make sure we found the mutant residue. if not, stop here.
	if ( ! found_mutant_residue ) {
		utility_exit_with_message( "Mutant residue position not found. Check PDB numbering and/or chain.\n" );
	}

	#ifndef NDEBUG
		TR << "Finished creating all TaskOperation's and TaskFactory's. Creating MoveMap." << std::endl;
	#endif

	// define some variables to be used with the movers coming up
	Size pack_cycles = (Size)basic::options::option[packing::ndruns].value();
	utility::vector1< pose::Pose > repacked_mutant_poses = utility::vector1< pose::Pose >( pack_cycles );

	basic::MetricValue< std::set< core::Size > > mutant_position_neighborhood;
	mutant_pose.metric( "mutant_nb_calculator", "neighbors", mutant_position_neighborhood );

	kinematics::MoveMapOP movemap = new core::kinematics::MoveMap();
	std::set< core::Size >::iterator iter;
	for( iter = mutant_position_neighborhood.value().begin(); iter != mutant_position_neighborhood.value().end(); iter++ ) {
		// movemap_->set_bb(i, true); // don't do any backbone minimization
		movemap->set_chi( *iter, true ); // but do minimize the side chains
	}
	#ifndef NDEBUG
		TR << "Movemap created... " << std::endl;
		//movemap->show( std::cout, mutant_pose.n_residue() );
	#endif

	#ifndef NDEBUG
		TR << "Beginning repacking/minimization of mutant pose." << std::endl;
	#endif

	// the repack mover
	moves::SetReturningPackRotamersMoverOP mutant_repacker = new protocols::moves::SetReturningPackRotamersMover( pack_cycles );
	mutant_repacker->task_factory( mutant_tf );
	mutant_repacker->score_function( scorefxn );
	mutant_repacker->apply( mutant_pose );
	mutant_repacker->get_repacked_poses( repacked_mutant_poses );

	// the side-chain minimization mover
	moves::MinMoverOP min_mover = new moves::MinMover( movemap, scorefxn, option[ OptionKeys::run::min_type ].value(), 0.01, true /*use_nblist*/ );
	moves::TaskAwareMinMoverOP task_aware_min_mover = new moves::TaskAwareMinMover( min_mover, mutant_tf );

	Energy sum_repacked_scores_mutant = 0.0;
	Energy average_repacked_score_mutant = 0.0;
	Energy score = 0.0;

	// minimize the same set of side chain in all of the repacked poses
	for ( Size ii = 1; ii <= repacked_mutant_poses.size(); ++ii ) {

		// minimize the side chains
		task_aware_min_mover->apply( repacked_mutant_poses[ ii ] );

		// re-score the minimized pose and output it (with the mutation string in the filename)
		score = (*scorefxn)( repacked_mutant_poses[ ii ] );
		#ifndef NDEBUG
			TR << "mutant pose " << ii << ", after minimization total energy: " << score << std::endl;
		#endif
		sum_repacked_scores_mutant += score;

		std::stringstream out;
		out << wt_residue_type << resid << mutant_residue_type;
		std::string filename = pdb_file_name.base() + "." + out.str() + "." + I( 3, 3, ii ) + ".pdb";

		repacked_mutant_poses[ ii ].dump_scored_pdb( filename, *(scorefxn()) );
	}

	average_repacked_score_mutant = sum_repacked_scores_mutant / repacked_mutant_poses.size();

	#ifndef NDEBUG
		TR << "Beginning repacking/minimization of wt pose." << std::endl;
	#endif

	//
	// we need to repack the entire Pose so that when we go to try the desired mutation, we don't get big effects
	// from the other score terms. but, the energy we want to keep is the average energy of 20 repacked poses.  if we call
	// pack_rotamers with ndruns 20, that will return only the Pose with the lowest packer energy.  what we need is a
	// pack method that returns a vector of Pose objects representing all of the packed poses.  Then we can go through
	// and score all of them to get the average energy.  But then which one do we go on to do design with?  Might as well
	// be the best one!
	//
	// We want to repack only the residues that get repacked when designing in the mutation. We also want to do some
	// side-chain minimization on those residues.  It's very slow to repack the entire structure and then would be
	// prohibitively slow to minimize the entire structure.  So instead, what we need to do is make the desired mutation
	// first and then repack/minimize the same residues that are changed in the mutant structure.
	//

	utility::vector1< pose::Pose > repacked_poses = utility::vector1< pose::Pose >( pack_cycles );

	moves::SetReturningPackRotamersMoverOP wt_repacker = new protocols::moves::SetReturningPackRotamersMover( pack_cycles );
	wt_repacker->task_factory( wt_tf );
	wt_repacker->score_function( scorefxn );
	wt_repacker->apply( pose );
	wt_repacker->get_repacked_poses( repacked_poses );

	// the side-chain minimization mover
	moves::MinMoverOP wt_min_mover = new moves::MinMover( movemap, scorefxn, option[ OptionKeys::run::min_type ].value(), 0.01, true /*use_nblist*/ );
	moves::TaskAwareMinMoverOP wt_task_aware_min_mover = new moves::TaskAwareMinMover( wt_min_mover, wt_tf );

	Energy sum_repacked_scores_wt = 0.0;
	Energy average_repacked_score_wt = 0.0;

	for ( Size ii=1; ii <= repacked_poses.size(); ++ii ) {

		// minimize the side chains
		wt_task_aware_min_mover->apply( repacked_poses[ ii ] );

		// re-score the minimized pose and output it (with the mutation string in the filename)
		score = (*scorefxn)( repacked_poses[ ii ] );
		#ifndef NDEBUG
			TR << "wt pose " << ii << ", after minimization total energy: " << score << std::endl;
		#endif
		sum_repacked_scores_wt += score;

		std::string filename = pdb_file_name.base() + "." + I( 3, 3, ii ) + ".pdb";
		repacked_poses[ ii ].dump_scored_pdb( filename, *(scorefxn()) );
	}

	average_repacked_score_wt = sum_repacked_scores_wt / repacked_poses.size();



	//
	// Now print out all of the score/ddG information.
	//
	std::cout << "\nScore summary:" << std::endl;

	std::cout << "wt pose:" << std::endl;
	scorefxn->show( std::cout, pose );
	std::cout << std::endl;

	std::cout << "mutant pose:" << std::endl;
	scorefxn->show( std::cout, mutant_pose );
	std::cout << std::endl;

	Real ddG = average_repacked_score_mutant - average_repacked_score_wt;
	std::cout << A( "ddG total:" ) << F( 12,3, ddG ) << std::endl;

	std::cout << "ddG by energy term:" << std::endl;
	std::cout << "wt:";
	print_energies( pose.energies().total_energies(), scorefxn );
	std::cout << "mut:";
	print_energies( mutant_pose.energies().total_energies(), scorefxn );
	scoring::EnergyMap ddG_energies = mutant_pose.energies().total_energies();
	ddG_energies -= pose.energies().total_energies();
	std::cout << "diff:";
	print_energies( ddG_energies, scorefxn );

	if ( option[ analyze_ddG_stability::calculate_surface_energy ] ) {
		Real mutant_pose_surface_energy = mutant_pose.energies().total_energies()[ scoring::surface ];
		Real wt_pose_surface_energy = pose.energies().total_energies()[ scoring::surface ];
		std::cout << A( "wt pose surface energy:  " ) << F( 9,3, wt_pose_surface_energy ) << std::endl;
		std::cout << A( "mutant pose surface energy:  " ) << F( 9,3, mutant_pose_surface_energy ) << std::endl;
		Real ddG_surface = mutant_pose_surface_energy - wt_pose_surface_energy;
		std::cout << A( "ddG surface:  " ) << F( 9,3, ddG_surface ) << std::endl;
	}

}
