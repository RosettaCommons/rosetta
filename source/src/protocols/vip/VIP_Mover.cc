// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file VIP_Mover.cc
/// @brief Void Identification and Packing Mover. Identifies buried voids in an input structure
/// and attempts to fill these using fixed backbone packing and the GOE energy
/// @author ben bborgo@genetics.wustl.edu

#include <protocols/vip/VIP_Mover.hh>
#include <protocols/vip/VIP_Report.hh>
#include <protocols/vip/VIP_Utils.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/cp.OptionKeys.gen.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/PDBInfo.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/MiniRelax.hh>
#include <protocols/relax/ClassicRelax.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/AddCavitiesMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/ScoreMover.hh>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>

core::Size num_stored_poses( 0 );


namespace protocols {
namespace vip {


static basic::Tracer TR( "VIP" );

VIP_Mover::VIP_Mover() {}
VIP_Mover::VIP_Mover(
	core::pose::Pose p,
	core::pose::Pose cp,
	core::pose::Pose fp,
	core::pose::Pose fup,
	utility::vector1<core::conformation::ResidueOP> tr,
	utility::vector1<core::Size> tp,
	utility::vector1<core::Real> te,
	utility::vector1<core::conformation::ResidueOP> vr,
	utility::vector1<core::Size> fav_p,
	utility::vector1<core::Real> fav_e,
	core::Size nc,
	utility::vector1<core::Size> cb,
	utility::vector1<std::string> fm,
	utility::vector1<core::Size> vn,
	utility::vector1<core::Size> vm,
	core::Real fe ){
	initial_pose = p;
	cavity_pose = cp;
	final_pose = fp;
	final_unrelaxed_pose = fup;
	temp_residues = tr;
	temp_positions = tp;
	temp_energies = te;
	favorable_residues = vr;
	favorable_positions = fav_p;
	favorable_energies = fav_e;
	number_cavities = nc;
	cavity_balls = cb;
	favorable_mutations = fm;
	void_neighbors = vn;
	void_mutatables = vm;
	final_energy = fe;
	energy_to_beat = 99999.0;
	iteration_ = 0;
	use_stored_best_energy = false;
}

VIP_Mover::~VIP_Mover()= default;

void VIP_Mover::set_initial_pose( core::pose::Pose pose ){
	initial_pose = pose;
}


void VIP_Mover::minimize_conformation(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::pose::Pose pose = initial_pose;
	core::scoring::ScoreFunctionOP sf2 = core::scoring::ScoreFunctionFactory::create_score_function( option[cp::minimizer_score_fxn] );
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_jump(false);
	movemap->set_chi(true);
	movemap->set_bb(true);
	protocols::moves::MoverOP min_native( new protocols::simple_moves::MinMover( movemap, sf2, "lbfgs_armijo_nonmonotone", 1e-2, true ) );
	min_native->apply( pose );
	initial_pose = pose;
}


void VIP_Mover::compute_number_cavities(){
	core::Size numsuck = 0;
	for ( core::Size a = 1 ; a < (cavity_pose.size()+1) ; a++ ) {
		if ( cavity_pose.residue(a).name() == "SUCK" ) {
			numsuck = numsuck + 1;
		}
	}
	number_cavities = numsuck;
}


void VIP_Mover::apply_holes(){
	core::pose::Pose pose = initial_pose;
	protocols::simple_moves::AddCavitiesMover cavget;
	cavget.apply( pose );
	cavity_pose = pose;
}


void VIP_Mover::get_cavity_positions(){
	utility::vector1<core::Size> cav_positions;
	for ( core::Size i = 1; i <= cavity_pose.size(); i++ ) {
		if ( cavity_pose.residue(i).name() == "SUCK" ) {
			cav_positions.push_back(i);
		}
	}
	cavity_balls = cav_positions;
}


void VIP_Mover::dump_pdb_to_file( core::pose::Pose & posey, std::string filename ){
	posey.dump_pdb(filename);
}


core::Real
VIP_Mover::get_cav_approx( core::Size a ){
	numeric::xyzVector< core::Real > cppos = cavity_pose.residue( cavity_balls[a] ).xyz(1);
	core::Real min = 99999.9;
	for ( core::Size i = 1; i <= initial_pose.size(); i++ ) {
		for ( core::Size j = 1; j <= initial_pose.residue(i).nheavyatoms(); j++ ) {
			if ( initial_pose.residue(i).xyz(j).distance( cppos ) < min ) {
				min = initial_pose.residue(i).xyz(j).distance( cppos );
			}
		}
	}
	return min;
}

void VIP_Mover::get_neighbors(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	utility::vector1<core::Size> neighbors;

	if ( cavity_balls.size() < 1 ) {
		initial_pose.dump_pdb("intermediate.pdb");
		TR.Error << "No more cavities! Dumping structure to intermediate.pdb" << std::endl;
		return;
	}

	for ( core::Size i = 1; i <= cavity_balls.size(); i++ ) {
		numeric::xyzVector<core::Real> cav_center = cavity_pose.residue( cavity_balls[i]).xyz(1);
		core::Real min = get_cav_approx( i );
		for ( core::Size a = 1; a <= initial_pose.size(); a++ ) {
			for ( core::Size b = 1; b <= initial_pose.residue(a).nheavyatoms(); b++ ) {
				numeric::xyzVector<core::Real> test_position = initial_pose.residue(a).xyz(b);
				if ( test_position.distance(cav_center) <= (option[cp::cutoff]+min) ) {
					neighbors.push_back( a );
					neighbors.push_back( b );
				}
			}
		}
	}
	void_neighbors = neighbors;
}


void
VIP_Mover::cull_mutatable_residues(){

	utility::vector1<core::Size> mutatable_residues;
	for ( core::Size i = 1; i <= void_neighbors.size(); i+=2 ) {
		if ( (initial_pose.residue(void_neighbors[i]).is_surface() == false) &&
				(initial_pose.residue(void_neighbors[i]).is_polar() == false) &&
				(std::find( excluded_positions.begin(), excluded_positions.end(), void_neighbors[i]) == excluded_positions.end() ) &&
				(initial_pose.residue(void_neighbors[i]).aa() < core::chemical::num_canonical_aas ) && // Don't try amino acid if APOLAR will choke on it
				(initial_pose.residue(void_neighbors[i]).atom_is_backbone(void_neighbors[i+1]) == false ) ) {
			mutatable_residues.push_back(void_neighbors[i]);
		}
	}

	utility::vector1<core::Size> mm_res;
	for ( core::Size ii = 1; ii < mutatable_residues.size(); ii++ ) {
		if ( mutatable_residues[ii] != mutatable_residues[ii+1] ) {
			mm_res.push_back( mutatable_residues[ii] );
		}
	}

	std::sort( mm_res.begin(), mm_res.end() );
	auto new_end_itr = std::unique( mm_res.begin(), mm_res.end() );

	utility::vector1<core::Size> unique_mm_res;
	std::copy( mm_res.begin(), new_end_itr, std::back_inserter( unique_mm_res ) );

	void_mutatables = unique_mm_res;
}


void VIP_Mover::try_point_mutants(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	for ( core::Size i = 1; i <= void_mutatables.size(); i++ ) {
		temp_positions.push_back(void_mutatables[i]);
	}

	// TR << "Num void_mutatables is " << void_mutatables.size() << std::endl;
	// TR << "Num temp_positions  is " << temp_positions.size() << std::endl;

	core::pack::task::ResfileCommandOP command( new core::pack::task::APOLAR );

	for ( core::Size aa = 1; aa <= void_mutatables.size(); aa++ ) {
		core::pose::Pose pack_pose;
		pack_pose = initial_pose;
		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pack_pose ));
		core::scoring::ScoreFunctionOP score_fxn = core::scoring::ScoreFunctionFactory::create_score_function( option[cp::pack_sfxn] );
		core::pack::task::TaskFactoryOP main_task_factory( new core::pack::task::TaskFactory );
		for ( core::Size j = 1; j <= pack_pose.size(); j++ ) {
			if ( j != void_mutatables[aa] ) {
				task->nonconst_residue_task(j).prevent_repacking();
			} else {
				task->nonconst_residue_task(void_mutatables[aa]).or_ex1(true);
				task->nonconst_residue_task(void_mutatables[aa]).or_ex2(true);
				task->nonconst_residue_task(void_mutatables[aa]).or_ex3(true);
				command->residue_action(*task,void_mutatables[aa]);
			}
		}
		protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover(score_fxn, task) );
		pack_mover->apply( pack_pose );
		// Store this single residue
		temp_residues.push_back( pack_pose.residue(temp_positions[aa]).clone() );
		temp_energies.push_back( pack_pose.energies().total_energy() );
	}

	// TR << "Done try_point_mutants" << std::endl;

}


void VIP_Mover::print_pack_report(){
	VIP_Report();
	VIP_Report vip_report;

	vip_report.get_GOE_repack_report( initial_pose, favorable_energies, favorable_residues, favorable_positions, iteration_,  false, energy_to_beat );
}


void VIP_Mover::sort_fill_energies(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::scoring::ScoreFunctionOP score_fxn = core::scoring::ScoreFunctionFactory::create_score_function( option[cp::pack_sfxn] );
	protocols::simple_moves::ScoreMoverOP score_em( new protocols::simple_moves::ScoreMover(score_fxn) );
	score_em->apply( initial_pose );

	core::Real baseE = initial_pose.energies().total_energy();
	core::pose::Pose basePose = initial_pose;
	for ( core::Size i = 1; i <= temp_energies.size(); i++ ) {
		core::Real tempE = temp_energies[i];
		if ( tempE < baseE ) {
			if ( initial_pose.residue(temp_positions[i]).name() != temp_residues[i]->name() ) {
				//core::pose::Pose temp_pose;
				//temp_pose = initial_pose;
				//temp_pose.replace_residue( temp_positions[i], *(temp_residues[i]), true );
				favorable_residues.push_back( temp_residues[i] );
				favorable_positions.push_back( temp_positions[i] );
				favorable_energies.push_back( temp_energies[i] );
				num_stored_poses++;
			}
		}
	}

	if ( option[ cp::print_reports ] ) {
		print_pack_report();
	}
}


void VIP_Mover::print_relax_report(){
	VIP_Report();
	VIP_Report vip_report;

	vip_report.get_GOE_relaxed_report( initial_pose, favorable_energies, favorable_residues, favorable_positions, iteration_, use_stored_best_energy, energy_to_beat );
	//  vip_report.get_GOE_packstat_report( initial_pose, favorable_energies );
}

void VIP_Mover::skip_relax(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::Real bestE = 99999;
	core::Size bestP = 0;
	for ( core::Size i = 1; i <= favorable_energies.size(); i++ ) {
		if ( favorable_energies[i] < bestE ) {
			bestE = favorable_energies[i];
			bestP = i;
		}
	}

	if ( favorable_positions.size() == 0 || favorable_residues.size() == 0 ) { // Avoid segmentation fault
		TR.Warning << "No favorable positions found." << std::endl;
		return;
	}
	final_pose = initial_pose;
	final_pose.replace_residue( favorable_positions[bestP], *(favorable_residues[bestP]), true );

	dump_pdb_to_file( final_pose, "final.pdb" );
	final_energy = bestE;
}


void VIP_Mover::relax_favorable_poses(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::Real bestE( 99999.0 );
	core::Size best_index( 99999 );

	core::scoring::ScoreFunctionOP relax_score_fxn = core::scoring::ScoreFunctionFactory::create_score_function( option[cp::relax_sfxn] );

	std::string rmover = option[ cp::relax_mover ];

	for ( core::Size i = 1; i <= favorable_residues.size(); i++ ) {

		core::pose::Pose relax_pose = initial_pose;
		relax_pose.replace_residue( favorable_positions[i], *(favorable_residues[i]), true );

		if ( rmover == "relax" ) {
			protocols::relax::RelaxProtocolBaseOP relaxmover( new protocols::relax::FastRelax( relax_score_fxn, 15 ) );
			if ( option[ cp::local_relax ] ) {
				core::kinematics::MoveMapOP mmap_ptr( new core::kinematics::MoveMap );
				if ( option[ cp::local_relax ] ) {
					set_local_movemap( relax_pose, favorable_positions[i], mmap_ptr );
				}
				relaxmover->set_movemap( mmap_ptr );
			}
			relaxmover->apply(relax_pose);
		} else if ( rmover == "classic_relax" ) {
			protocols::relax::RelaxProtocolBaseOP relaxmover( new protocols::relax::ClassicRelax( relax_score_fxn ) );
			if ( option[ cp::local_relax ] ) {
				core::kinematics::MoveMapOP mmap_ptr( new core::kinematics::MoveMap );
				if ( option[ cp::local_relax ] ) {
					set_local_movemap( relax_pose, favorable_positions[i], mmap_ptr );
				}
				relaxmover->set_movemap( mmap_ptr );
			}
			relaxmover->apply(relax_pose);
		} else if ( rmover == "cst_relax" ) {
			protocols::relax::RelaxProtocolBaseOP cstrelaxmover( new protocols::relax::MiniRelax( relax_score_fxn ) );
			if ( option[ cp::local_relax ] ) {
				core::kinematics::MoveMapOP mmap_ptr( new core::kinematics::MoveMap );
				if ( option[ cp::local_relax ] ) {
					set_local_movemap( relax_pose, favorable_positions[i], mmap_ptr );
				}
				cstrelaxmover->set_movemap( mmap_ptr );
			}
			cstrelaxmover->apply(relax_pose);
		}

		core::scoring::ScoreFunctionOP sf2 =
			core::scoring::ScoreFunctionFactory::create_score_function( option[cp::relax_sfxn] );
		protocols::simple_moves::ScoreMoverOP score_em( new protocols::simple_moves::ScoreMover(sf2) );
		score_em->apply( relax_pose );
		favorable_energies[i] = relax_pose.energies().total_energy();

		if ( favorable_energies[i] < bestE ) {
			best_index = i;
			bestE = favorable_energies[i];
			final_pose = relax_pose;
			final_energy = bestE;
		}
	}


	if ( best_index < 99999 ) {
		//TR << "STASHING UNRELAXED POSE" << std::endl;
		final_unrelaxed_pose = initial_pose;
		final_unrelaxed_pose.replace_residue( favorable_positions[best_index],
			*(favorable_residues[best_index]), true );
	}

}

void VIP_Mover::sort_relaxed_poses(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ cp::print_reports ] ) {
		print_relax_report();
	}
}

void VIP_Mover::nook_finder(){
	minimize_conformation();
	apply_holes();
	get_cavity_positions();
	get_neighbors();
	cull_mutatable_residues();
}

void VIP_Mover::cranny_packer(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	try_point_mutants();
	sort_fill_energies();
	if ( option[ cp::skip_relax ] ) {
		skip_relax();
	} else {
		relax_favorable_poses();
		sort_relaxed_poses();
	}
}

void VIP_Mover::apply(){
	set_excluded_positions();
	nook_finder();
	cranny_packer();
}

void
VIP_Mover::set_excluded_positions() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	excluded_positions.clear();

	if ( option[ cp::exclude_file ].active() ) {
		std::string exclude_file_name( option[ cp::exclude_file ] );

		// Probe for file
		std::ifstream exc_file( exclude_file_name.c_str() );

		if ( !exc_file ) {
			TR << "Exclude_file " << exclude_file_name << " not found." << std::endl;
			TR << "No positions will be excluded." << std::endl;
			return;
		}

		// Process one at a time

		core::Size exc_pos;
		char exc_chain;

		exc_file >> exc_pos;
		while ( !exc_file.eof() ) {
			exc_file >> exc_chain;
			TR << "Adding position " << exc_pos << " chain " << exc_chain << " to exclude list." << std::endl;
			excluded_positions.push_back( initial_pose.pdb_info()->pdb2pose( exc_chain, exc_pos ) );
			exc_file >> exc_pos;
		}
		exc_file.close();
	}

	// Exclude any non-amino acid residues from mutation
	for ( core::Size i(1), ei( initial_pose.size() ) ; i <= ei ; ++i ) {
		if ( !initial_pose.residue( i ).is_protein() &&
				( std::find(excluded_positions.begin(), excluded_positions.end(), i ) ==
				excluded_positions.end() ) ) {
			excluded_positions.push_back( i );
		}
	}

	if ( excluded_positions.size() > 0 ) {
		TR << "Found " << excluded_positions.size() << " positions excluded from mutation" << std::endl;
	} else {
		TR << "No positions will be excluded." << std::endl;
	}
	return;
}

bool
are_seqs_different( core::pose::Pose & p1, core::pose::Pose & p2 ) {
	for ( core::Size j = 1; j <= p1.size(); j++ ) {
		if ( p1.residue(j).name() != p2.residue(j).name() ) {
			return true;
		}
	}
	return false;
}


}}
