// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

// libRosetta headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>
#include <basic/options/option.hh>
//#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueConnection.hh>

#include <core/scoring/TenANeighborGraph.hh>

#include <protocols/simple_moves/FragmentMover.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
// Headers for RNA
#include <protocols/farna/fragments/FullAtomRNA_Fragments.hh>
#include <basic/database/open.hh>
#include <protocols/farna/movers/RNA_FragmentMover.hh>
#include <protocols/farna/movers/RNA_FragmentMover.fwd.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <core/chemical/ResidueType.hh>


#include <devel/domain_assembly/domain_assembly.hh>

#include <utility/vector1.hh>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <utility/io/izstream.hh>
#include <string>


// option key includes

#include <basic/options/keys/DomainAssembly.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector0.hh>

//Auto Headers
#include <protocols/viewer/GraphicsState.hh>


using namespace core;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::viewer;
using namespace basic::options;
using namespace farna;
using namespace farna::movers;
using namespace ObjexxFCL::format;
//using namespace basic::options::OptionKeys;
using namespace chemical;

// protocol specific options
// namespace DomainAssembly{
//  //FileOptionKey const da_setup_option_file( "DomainAssembly::da_setup_option_file" );
//  //FileOptionKey const da_setup_output_pdb( "DomainAssembly::da_setup_output_pdb" );
//   FileOptionKey const da_start_pdb( "DomainAssembly::da_start_pdb" );
//  FileOptionKey const da_linker_file( "DomainAssembly::da_linker_file" );
//  IntegerOptionKey const da_start_pdb_num( "DomainAssembly::da_start_pdb_num" );
//  IntegerOptionKey const da_nruns( "DomainAssembly::da_nruns" );
//  FileOptionKey const da_setup_output_pdb( "DomainAssembly::da_setup_output_pdb" );
// }

namespace devel {
namespace domain_assembly {

core::kinematics::MoveMapOP
read_movemap_from_da_linker_file()
{
	// read in linkers and set move map true for linker regions
	std::string filename_linkers = option[ OptionKeys::DomainAssembly::da_linker_file ]();
	utility::vector1< std::pair < Size, Size > >  linker_ranges;
	read_linker_file( filename_linkers, linker_ranges );
	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	set_movemap_for_linkers ( linker_ranges, mm );
	return mm;
}

/// @brief reads in file that specifies which regions of the protein will
///  will move during domain assembly
///  Each line of the file should have the start and end position for a linker region
bool
read_linker_file(
	std::string const & filename,
	utility::vector1< std::pair < Size, Size > > & linker_ranges
)
{
	utility::io::izstream data( filename );
	if ( !data ) {
		std::cerr << " can not open linker file: " << filename << std::endl;
		return false;
	}

	std::string line;
	while ( getline( data, line ) ) {
		std::istringstream line_stream( line );
		Size begin, end;
		line_stream >> begin >> end;
		if ( line_stream.fail() ) {
			std::cout << " can not parse line in linker file: " << line << std::endl;
			return false;
		}
		std::pair < Size, Size > range ( begin, end );
		linker_ranges.push_back( range );
	}
	data.close();
	data.clear();
	return true;
}

/// @brief function that parses the rna regions specified into a boolean array to be used by a RNA Fragment Mover Object
protocols::toolbox::AtomLevelDomainMapOP
set_moveable_rna(
	pose::Pose & full_pose,
	utility::vector1< std::pair < Size, Size > > & linker_rna
)
{
	using namespace protocols::toolbox;
	//FArray1D used to maintain RNA_FragmentMover compatability with other RNA protocols
	AtomLevelDomainMapOP atom_level_domain_map( new AtomLevelDomainMap( full_pose ) );
	for ( Size i = 1; i <= full_pose.total_residue(); ++i ) {
		atom_level_domain_map->set( i, false );
		if ( full_pose.residue_type(i).is_RNA() ) {
			//If the residue is RNA and a linker range, set it as moveable, otherwise it will remain immobile
			for ( Size ii = 1; ii <= linker_rna.size(); ++ii ) {
				if ( i >= linker_rna[ii].first && i <= linker_rna[ii].second ) {
					atom_level_domain_map->set( i, true );
					//If it's true for one linker range don't bother testing the following ranges
					break;
				}
			}
		}
	}
	return atom_level_domain_map;;
}


/// @brief sets movemap true for regions specified in linker file
void
set_movemap_for_linkers(
	utility::vector1< std::pair < Size, Size > > const & linker_ranges,
	kinematics::MoveMapOP & mm
)
{
	for ( Size i = 1; i <= linker_ranges.size(); ++i ) {
		std::pair< Size, Size > p;
		p = linker_ranges[ i ];
		for ( Size seqpos = p.first; seqpos <= p.second; ++seqpos ) {
			mm->set_bb( seqpos, true );
		}
	}
}

/// @brief centroid mode optimization of linkers
void
optimize_linkers_centroid_mode(
	kinematics::MoveMapOP & mm,
	pose::Pose & full_pose
)
{
	// Size inside_steps_stage1 ( 200 );
	//Size outside_steps_stage1 ( 10 );
	//Size inside_steps_stage2 ( 100 );
	//Size outside_steps_stage2 ( 5 );
	Size inside_steps_stage1 ( 50 );
	Size outside_steps_stage1 ( 15 );
	Size inside_steps_stage2 ( 50 );
	Size outside_steps_stage2 ( 5 );


	core::pose::Pose const saved_input_pose( full_pose ); //used to return sidechains later
	core::util::switch_to_residue_type_set( full_pose, core::chemical::CENTROID );
	scoring::ScoreFunctionOP scorefxn_centroid( scoring::ScoreFunctionFactory::create_score_function( scoring::CENTROID_WTS ) );
	MonteCarloOP mc( new MonteCarlo( full_pose, *scorefxn_centroid, 0.8 /*temperature*/ ) );

	// read fragments file
	core::fragment::ConstantLengthFragSetOP fragset3mer = NULL;
	if ( basic::options::option[ basic::options::OptionKeys::in::file::frag3].user() ) {
		fragset3mer = core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( 3 ) );
		fragset3mer->read_fragment_file( basic::options::option[ basic::options::OptionKeys::in::file::frag3 ]() );
	}

	// more traditional small moves
	protocols::simple_moves::SmallMoverOP small_mover( new protocols::simple_moves::SmallMover( mm, 0.8/*temp*/, 1/*nmoves*/ ) );
	small_mover->angle_max( 'H', 2.0 );  // max angle displacement 180 degrees
	small_mover->angle_max( 'E', 4.0 );
	small_mover->angle_max( 'L', 4.0 );

	//// STAGE 1 /////
	TrialMoverOP centroid_trial_mover = NULL;
	if ( fragset3mer ) {
		protocols::simple_moves::FragmentMoverOP frag_mover( new protocols::simple_moves::ClassicFragmentMover(fragset3mer, mm) );
		centroid_trial_mover = TrialMoverOP( new TrialMover( frag_mover, mc ) );
	} else {
		// if no fragments, use coarse small moves
		Size nmoves ( 1 );
		protocols::simple_moves::SmallMoverOP coarse_small_mover( new protocols::simple_moves::SmallMover( mm, 0.8/*temp*/, nmoves ) );
		coarse_small_mover->angle_max( 'H', 180.0 );  // max angle displacement 180 degrees
		coarse_small_mover->angle_max( 'E', 180.0 );
		coarse_small_mover->angle_max( 'L', 180.0 );
		centroid_trial_mover = TrialMoverOP( new TrialMover( coarse_small_mover, mc ) );
	}

	RepeatMoverOP inner_centroid_loop( new RepeatMover( centroid_trial_mover, inside_steps_stage1 ) );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= outside_steps_stage1; ++i ) {
		inner_centroid_loop -> apply( full_pose );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
	}
	std::cout << "end of stage one centroid refinement" << std::endl;
	//// END STAGE 1 ////

	////// STAGE 2 ///////
	TrialMoverOP stage2_trial( new TrialMover( small_mover, mc ) );
	RepeatMoverOP stage2( new RepeatMover( stage2_trial, inside_steps_stage2 ) );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= outside_steps_stage2; ++i ) {
		stage2 -> apply( full_pose );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
	}
	std::cout << "end of stage two centroid refinement" << std::endl;
	///////  END_STAGE 2 ////////

	mc -> recover_low( full_pose );
	protocols::simple_moves::ReturnSidechainMover return_sidechains( saved_input_pose );
	return_sidechains.apply( full_pose );

	// for graphics:
	// protocols::viewer::add_conformation_viewer( full_pose.conformation(), "full_pose", 400, 400 );
}

/// @brief a helper function for the domain assembly protocol. Selects
///residues near linkers and domain interfaces for repacking
void
da_residues_to_repack(
	kinematics::MoveMapOP & mm,
	utility::vector1< std::pair < Size, Size > > & nearest_movable_residues,
	pose::Pose & pose,
	utility::vector1<bool> & repack_residues
)
{
	pose.update_residue_neighbors();
	core::scoring::TenANeighborGraph const & nb_graph( pose.energies().tenA_neighbor_graph() );

	for ( Size i=1; i <= pose.total_residue(); ++i ) {

		core::graph::Node const * current_node( nb_graph.get_node(i) );  // neighbor graph node for residue i
		if ( mm->get_bb(i) ) {
			repack_residues[ i ] = true;   // movable residues
			for ( core::graph::Node::EdgeListConstIter it = current_node->const_edge_list_begin();
					it != current_node->const_edge_list_end(); ++it ) {
				repack_residues[ (*it)->get_other_ind(i) ] = true;  // neighbors of movable residues
			}
		} else {
			// residues at the interdomain interface
			for ( core::graph::Node::EdgeListConstIter it = current_node->const_edge_list_begin();
					it != current_node->const_edge_list_end(); ++it ) {
				Size nb = (*it)->get_other_ind(i);
				std::pair< Size, Size > p;
				p = nearest_movable_residues[ i ];
				if ( nb <= p.first || nb >= p.second ) {
					repack_residues[ i ] = true;
				}
			}
		}   // if mm true
	}   // for over all residues
}

/// @brief a helper function for the domain assembly protocol.  For each residue
/// it finds the closest N-terminal and C-terminal movable residue (as specified
/// in the input movemap)
void
find_nearest_movable_residues(
	kinematics::MoveMapOP & mm,
	pose::Pose & pose,
	utility::vector1< std::pair < Size, Size > > & nearest_movable_residues
)
{
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		// search for nearest movable residue towards the N-terminus
		Size low_residue = 0;
		for ( Size j = i-1; j > 0; --j ) {
			if ( mm->get_bb( j ) ) {
				low_residue = j;
				break;
			}
		}
		// search for the nearest movable residue towards the C-terminus
		Size high_residue = pose.total_residue() + 1;
		for ( Size j = i+1; j <= pose.total_residue(); ++j ) {
			if ( mm->get_bb( j ) ) {
				high_residue = j;
				break;
			}
		}
		std::pair < Size, Size > p1 ( low_residue, high_residue );
		nearest_movable_residues.push_back( p1 );
	}
}

void
optimize_linkers_fullatom_mode(
	kinematics::MoveMapOP & mm,
	pose::Pose & full_pose
)
{
	//core::util::switch_to_residue_type_set( full_pose, core::chemical::FA_STANDARD );
	scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );

	// for each residue - identify the nearest movable residue forward and backwards in sequence
	// this will be used to help determine which side chains to repack
	utility::vector1 < std::pair < Size, Size > > nearest_movable_residues;
	find_nearest_movable_residues( mm, full_pose, nearest_movable_residues );

	// global repack of the side chains
	pack::task::PackerTaskOP base_packer_task( pack::task::TaskFactory::create_packer_task( full_pose ));
	base_packer_task -> initialize_from_command_line();
	base_packer_task -> restrict_to_repacking();
	base_packer_task -> set_bump_check( true );
	base_packer_task -> or_include_current( true );
	pack::pack_rotamers( full_pose, *scorefxn, base_packer_task );

	// Add the constraints, if any, to the score function
	core::scoring::constraints::add_fa_constraints_from_cmdline(full_pose, *scorefxn);


	MonteCarloOP mc( new MonteCarlo( full_pose, *scorefxn, 0.8 /*temperature*/ ) );

	// MOVER: small moves
	protocols::simple_moves::SmallMoverOP small_mover( new protocols::simple_moves::SmallMover( mm, 0.8/*temp*/, 1 ) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 3.0 );
	small_mover->angle_max( 'L', 4.0 );

	// MOVER: rotamer trials
	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial_mover( new protocols::simple_moves::EnergyCutRotamerTrialsMover( scorefxn, *base_packer_task, mc, 0.01 /*energycut*/ ) );

	// MOVER minimization
	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( mm, scorefxn, "dfpmin", 0.001, true /*use_nblist*/ ) );
	//TrialMoverOP min_trial_mover( new TrialMover( min_mover, mc ) );

	// Initial Minimization //
	//std::cout << " Score before initial fullatom minimization " << mc->last_accepted_score() << std::endl;
	//min_mover -> apply (full_pose);
	//std::cout << " Score after initial fullatom minimization " << mc->last_accepted_score() << std::endl;

	//// STAGE1 ////
	SequenceMoverOP stage1_seq( new SequenceMover );
	stage1_seq -> add_mover( small_mover );
	stage1_seq -> add_mover( pack_rottrial_mover );
	TrialMoverOP stage1_trial( new TrialMover( stage1_seq, mc ) );
	protocols::moves::RepeatMoverOP stage1_inner_loop( new RepeatMover( stage1_trial, 15 /*100 cycles*/ ) );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= 10; ++i ) {
		stage1_inner_loop -> apply( full_pose );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
		utility::vector1<bool> repack_residues( full_pose.total_residue(), false );
		da_residues_to_repack( mm, nearest_movable_residues, full_pose, repack_residues );
		pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
		this_packer_task -> restrict_to_residues( repack_residues );
		pack::pack_rotamers( full_pose, *scorefxn, this_packer_task );
	}
	std::cout << "end of stage one fullatom refinement" << std::endl;
	//// END STAGE1 ////

	//// STAGE2 ////
	SequenceMoverOP stage2_seq( new SequenceMover );
	stage2_seq->add_mover( small_mover );
	stage2_seq->add_mover( pack_rottrial_mover );
	stage2_seq->add_mover( min_mover );
	TrialMoverOP stage2_trial( new TrialMover( stage2_seq, mc ) );
	RepeatMoverOP stage2_inner_loop( new RepeatMover( stage2_trial, 10 /*cycles*/ ) );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= 10; ++i ) {
		stage2_inner_loop -> apply( full_pose );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
		utility::vector1<bool> repack_residues( full_pose.total_residue(), false );
		da_residues_to_repack( mm, nearest_movable_residues, full_pose, repack_residues );
		pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
		this_packer_task -> restrict_to_residues( repack_residues );
		pack::pack_rotamers( full_pose, *scorefxn, this_packer_task );
	}
	std::cout << "end of stage two fullatom refinement" << std::endl;
	//// END STAGE2 ////
	mc -> recover_low( full_pose );
}


void
optimize_linkers_rna_fullatom_mode(
	kinematics::MoveMapOP & mm,
	pose::Pose & full_pose,
	protocols::farna::fragments::RNA_FragmentsOP & all_rna_fragments,
	utility::vector1< std::pair< Size, Size > > linker_rna
)
{
	// Set the score function to include energies from NA-NA interactions
	scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );
	scoring::methods::EnergyMethodOptions scorefxn_opt  = scorefxn -> energy_method_options();
	scorefxn_opt.exclude_DNA_DNA( false );
	scorefxn -> set_energy_method_options( scorefxn_opt );

	// for each residue - identify the nearest movable residue forward and backwards in sequence
	// this will be used to help determine which side chains to repack
	utility::vector1 < std::pair < Size, Size > > nearest_movable_residues;
	find_nearest_movable_residues( mm, full_pose, nearest_movable_residues );

	// global repack of the side chains
	pack::task::PackerTaskOP base_packer_task( pack::task::TaskFactory::create_packer_task( full_pose ));
	base_packer_task -> initialize_from_command_line();
	base_packer_task -> restrict_to_repacking();
	base_packer_task -> set_bump_check( true );
	base_packer_task -> or_include_current( true );
	pack::pack_rotamers( full_pose, *scorefxn, base_packer_task );
	MonteCarloOP mc_nocont( new MonteCarlo( full_pose, *scorefxn, 0.8 ) );

	Size rna_normalize_step ( 15 );
	Size inside_steps_stage1 ( 15 );
	Size outside_steps_stage1 ( 5 );
	Size inside_steps_stage2 ( 25 );
	Size outside_steps_stage2 ( 8 );
	//Size inside_steps_stage1 ( 50 );
	//Size outside_steps_stage1 ( 15 );
	//Size inside_steps_stage2 ( 50 );
	//Size outside_steps_stage2 ( 5 );

	// Create a rna fragment mover, set its fragsize, and then do an initialization move
	RNA_FragmentMoverOP rna_fragment_mover = RNA_FragmentMoverOP( new RNA_FragmentMover( *all_rna_fragments, set_moveable_rna( full_pose, linker_rna) ) );
	rna_fragment_mover -> set_frag_size( Size(1) );
	for ( Size i = 1; i <= rna_normalize_step; ++i ) {
		rna_fragment_mover -> apply(full_pose);
		mc_nocont -> boltzmann( full_pose, "frag" + SS( Size(1) ));
	}

	// Add the constraints, if any, to the score function
	scoring::ScoreFunctionOP scorefxn_cont( scoring::get_score_function() );
	core::scoring::constraints::add_fa_constraints_from_cmdline(full_pose, *scorefxn_cont);

	// Set the score function with constraints to also include NA-NA interactions
	scoring::methods::EnergyMethodOptions scorefxn_cont_opt = scorefxn_cont -> energy_method_options();
	scorefxn_cont_opt.exclude_DNA_DNA( false );
	scorefxn_cont -> set_energy_method_options( scorefxn_cont_opt );
	MonteCarloOP mc( new MonteCarlo( full_pose, *scorefxn_cont, 0.8) );

	// MOVER: small moves
	protocols::simple_moves::SmallMoverOP small_mover( new protocols::simple_moves::SmallMover( mm, 0.8, 1 ) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 4.0 );
	small_mover->angle_max( 'L', 4.0 );

	// MOVER: rotamer trials
	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial_mover( new protocols::simple_moves::EnergyCutRotamerTrialsMover( scorefxn_cont, *base_packer_task, mc, 0.01 ) );

	// MOVER: minimization
	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( mm, scorefxn_cont, "dfpmin", 0.001, true ) );

	// MOVER: coarse small moves
	protocols::simple_moves::SmallMoverOP coarse_small_mover( new protocols::simple_moves::SmallMover( mm, 0.8, 1 ) );
	coarse_small_mover->angle_max( 'H', 180.0 );  // max angle displacement 180 degrees
	coarse_small_mover->angle_max( 'E', 180.0 );
	coarse_small_mover->angle_max( 'L', 180.0 );

	// MOVER: coarse small moves
	protocols::simple_moves::SmallMoverOP semi_small_mover( new protocols::simple_moves::SmallMover( mm, 0.8, 1 ) );
	semi_small_mover->angle_max( 'H', 45.0 );  // max angle displacement 180 degrees
	semi_small_mover->angle_max( 'E', 45.0 );
	semi_small_mover->angle_max( 'L', 45.0 );


	//// STAGE 1 ////
	min_mover -> apply( full_pose );
	SequenceMoverOP stage1_seq( new SequenceMover );
	stage1_seq -> add_mover( coarse_small_mover );
	stage1_seq -> add_mover( pack_rottrial_mover );
	//stage1_seq -> add_mover( min_mover );
	TrialMoverOP stage1_trial( new TrialMover( stage1_seq, mc ) );
	protocols::moves::RepeatMoverOP stage1_inner_loop( new RepeatMover( stage1_trial, inside_steps_stage1 ) );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= outside_steps_stage1; ++i ) {
		stage1_inner_loop -> apply( full_pose );
		utility::vector1<bool> repack_residues( full_pose.total_residue(), false );
		da_residues_to_repack( mm, nearest_movable_residues, full_pose, repack_residues );
		pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
		this_packer_task -> restrict_to_residues( repack_residues );
		pack::pack_rotamers( full_pose, *scorefxn_cont, this_packer_task );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
	}
	std::cout << "end of stage one fullatom refinement" << std::endl;


	//// END STAGE 1 ////


	////// STAGE 2 ///////
	SequenceMoverOP stage2_seq( new SequenceMover );
	stage2_seq -> add_mover( semi_small_mover );
	stage2_seq -> add_mover( pack_rottrial_mover );
	stage2_seq -> add_mover( rna_fragment_mover );
	TrialMoverOP stage2_trial( new TrialMover( stage2_seq, mc ) );
	protocols::moves::RepeatMoverOP stage2_inner_loop( new RepeatMover( stage2_trial, inside_steps_stage2 ) );

	for ( Size i = 1; i <= outside_steps_stage2; ++i ) {
		stage2_inner_loop -> apply( full_pose );
		utility::vector1<bool> repack_residues( full_pose.total_residue(), false );
		da_residues_to_repack( mm, nearest_movable_residues, full_pose, repack_residues );
		pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
		this_packer_task -> restrict_to_residues( repack_residues );
		pack::pack_rotamers( full_pose, *scorefxn_cont, this_packer_task );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
	}
	std::cout << "end of stage two fullatom refinement" << std::endl;
	mc -> recover_low( full_pose );

	///////  END_STAGE 2 ////////


	//////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////

	//// STAGE3 ////
	SequenceMoverOP stage3_seq( new SequenceMover );
	stage3_seq -> add_mover( rna_fragment_mover );
	stage3_seq -> add_mover( small_mover );
	stage3_seq -> add_mover( pack_rottrial_mover );
	TrialMoverOP stage3_trial( new TrialMover( stage3_seq, mc ) );
	protocols::moves::RepeatMoverOP stage3_inner_loop( new RepeatMover( stage3_trial, inside_steps_stage2 ) );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= outside_steps_stage2; ++i ) {
		stage3_inner_loop -> apply( full_pose );
		utility::vector1<bool> repack_residues( full_pose.total_residue(), false );
		da_residues_to_repack( mm, nearest_movable_residues, full_pose, repack_residues );
		pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
		this_packer_task -> restrict_to_residues( repack_residues );
		pack::pack_rotamers( full_pose, *scorefxn_cont, this_packer_task );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
	}
	mc -> recover_low( full_pose );
	std::cout << "end of stage three fullatom refinement" << std::endl;


	//// END STAGE3 ////

	//// STAGE4 ////
	SequenceMoverOP stage4_seq( new SequenceMover );
	stage4_seq->add_mover( rna_fragment_mover );
	stage4_seq->add_mover( small_mover );
	stage4_seq->add_mover( pack_rottrial_mover );
	stage4_seq->add_mover( min_mover );
	TrialMoverOP stage4_trial( new TrialMover( stage4_seq, mc ) );
	RepeatMoverOP stage4_inner_loop( new RepeatMover( stage4_trial, inside_steps_stage1  ) );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= outside_steps_stage1; ++i ) {
		stage4_inner_loop -> apply( full_pose );
		utility::vector1<bool> repack_residues( full_pose.total_residue(), false );
		da_residues_to_repack( mm, nearest_movable_residues, full_pose, repack_residues );
		pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
		this_packer_task -> restrict_to_residues( repack_residues );
		pack::pack_rotamers( full_pose, *scorefxn_cont, this_packer_task );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
	}
	std::cout << "end of stage four fullatom refinement" << std::endl;
	//// END STAGE4 ////

	mc -> recover_low( full_pose );
	mc -> clear_poses();
}


/// @brief optimizes linkers in a multidomain protein
/// expects an input pdb with multi-domains and an
/// input file specifying which torsion angles to optimize
void
assemble_domains_optimize()
{
	pose::Pose full_pose, start_pose;
	std::string filename_start = option[ OptionKeys::DomainAssembly::da_start_pdb ]();
	core::import_pose::pose_from_pdb( start_pose, filename_start );
	RNA_FragmentsOP all_rna_fragments;
	utility::vector1< std::pair < Size, Size > > linker_ranges_rna;
	full_pose = start_pose;
	std::string filename_linkers_rna = option[ OptionKeys::DomainAssembly::da_linker_file_rna]();
	//If a rna linker file is passed, fill the fragment library and the rna linker ranges
	if ( filename_linkers_rna != "--" ) {
		all_rna_fragments = RNA_FragmentsOP( new FullAtomRNA_Fragments( basic::database::full_name("sampling/rna/1jj2.torsions"  ) ) );
		read_linker_file( filename_linkers_rna, linker_ranges_rna);
	}

	// read in linkers and set move map true for linker regions
	kinematics::MoveMapOP mm = read_movemap_from_da_linker_file();

	int nruns = option[ OptionKeys::DomainAssembly::da_nruns ];
	int start_pdb_num = option[ OptionKeys::DomainAssembly::da_start_pdb_num ];

	for ( int i = 0; i < nruns; ++i ) {
		using namespace std;
		full_pose = start_pose;
		//IMPORTANT: If your start structure contains RNA, even if it is immobile, you must specify a rna linker file.  For immobile just point to an empty file
		if ( filename_linkers_rna == "--"  ) {
			//TR << " --------- Centroid mode optimization -------------------- " << std::endl;
			optimize_linkers_centroid_mode( mm, full_pose );
			//ostringstream outputfilename;
			//outputfilename << "test_stage2_" << setw(6) << setfill('0') << start_pdb_num+i << ".pdb";
			//full_pose.dump_pdb( outputfilename.str() );
			//TR << " ---------- Fullatom mode optimization -------------------- " << std::endl;
			optimize_linkers_fullatom_mode( mm, full_pose );
		} else {
			optimize_linkers_rna_fullatom_mode( mm, full_pose, all_rna_fragments, linker_ranges_rna );
		}

		scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );
		core::scoring::constraints::add_fa_constraints_from_cmdline(full_pose, *scorefxn);
		// final output
		std::ostringstream outputfilename;
		outputfilename << "test_stage3_" << setw(6) << setfill('0') << start_pdb_num+i << ".pdb";
		std::string fname = outputfilename.str();
		//core::io::pdb::dump_pdb( full_pose, out );
		full_pose.dump_scored_pdb( fname, *scorefxn );
		//std::ofstream out( fname.c_str() );
		//out << "total_energy: " << (*scorefxn)(full_pose) << std::endl;
	}

}

}
}
