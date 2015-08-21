// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ClashBasedRepackShellSelector.cc
/// @brief  The ClashBasedRepackShellSelector identifies all residues that clash with at least one rotamer of a design position
/// @author Noah Ollikainen (nollikai@gmail.com)
/// @author Roland A. Pache, PhD


// Unit headers
#include <core/pack/task/residue_selector/ClashBasedRepackShellSelector.hh>
#include <core/pack/task/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/types.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/packer_neighbors.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// C++ headers
#include <utility/assert.hh>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {


ClashBasedRepackShellSelector::ClashBasedRepackShellSelector()
{
	bump_overlap_factor_ = 0.5;
}

ClashBasedRepackShellSelector::~ClashBasedRepackShellSelector() {}

ClashBasedRepackShellSelector::ClashBasedRepackShellSelector(
	core::pack::task::PackerTaskOP packer_task,
	core::scoring::ScoreFunctionOP score_fxn
)
{
	packer_task_ = packer_task;
	score_fxn_ = score_fxn;
	bump_overlap_factor_ = 0.5;
}

// helper functions

bool ClashBasedRepackShellSelector::is_sc_sc_clash(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2
) const
{
	core::Size max_res1_heavyatoms = rsd1.nheavyatoms();
	core::Size max_res2_heavyatoms = rsd2.nheavyatoms();
	for ( core::Size m=5; m <= max_res1_heavyatoms; ++m ) {
		core::Vector const & atom1_coords( rsd1.xyz( m ) );
		for ( core::Size n=5; n <= max_res2_heavyatoms; ++n ) {
			core::Vector const & atom2_coords( rsd2.xyz( n ) );
			// use the sum of the LJ radii as the distance cutoff with a bump overlap factor
			core::Real lj_sum = (rsd1.atom_type( m ).lj_radius() + rsd2.atom_type( n ).lj_radius());
			core::Real distance_cutoff_squared = (lj_sum * lj_sum) * bump_overlap_factor_;
			core::Real distance_squared = ( atom1_coords - atom2_coords ).length_squared();
			if ( distance_squared < distance_cutoff_squared ) return true;
		}
	}
	return false;
}


bool ClashBasedRepackShellSelector::is_sc_bb_clash(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2
) const
{
	core::Size max_res1_heavyatoms = rsd1.nheavyatoms();
	for ( core::Size m=5; m <= max_res1_heavyatoms; ++m ) {
		core::Vector const & atom1_coords( rsd1.xyz( m ) );
		for ( core::Size n=1; n <= 4; ++n ) {
			core::Vector const & atom2_coords( rsd2.xyz( n ) );
			// use the sum of the LJ radii as the distance cutoff with a bump overlap factor
			core::Real lj_sum = (rsd1.atom_type( m ).lj_radius() + rsd2.atom_type( n ).lj_radius());
			core::Real distance_cutoff_squared = (lj_sum * lj_sum) * bump_overlap_factor_;
			core::Real distance_squared = ( atom1_coords - atom2_coords ).length_squared();
			if ( distance_squared < distance_cutoff_squared ) return true;
		}
	}
	return false;
}


utility::vector1<core::Size> ClashBasedRepackShellSelector::get_clashing_positions(
	core::pose::Pose const & pose,
	core::conformation::Residue const & rsd1,
	core::Size const resnum
) const
{
	utility::vector1<core::Size> clash_positions;
	// iterate over all neighbor residues
	for ( core::Size i=1; i <= pose.total_residue(); ++i ) {
		if ( i == resnum ) continue;
		core::conformation::Residue const & rsd2 = pose.conformation().residue(i);
		// check for backbone clash first (where repacking wouldn't help)
		if ( is_sc_bb_clash(rsd1,rsd2) ) continue;
		// check for clashes between all pairs of side-chain atoms between the query residue and the given neighbor
		bool found_sc_only_clash = false;
		if ( is_sc_sc_clash(rsd1,rsd2) ) {
			//std::cout << "sidechain clashes with " << i << std::endl;
			found_sc_only_clash=true;
			//check if this side-chain clash also involves backbone clashes with other residues in the neighborhood
			for ( core::Size j=1; j <= pose.total_residue(); ++j ) {
				if ( j == resnum ) continue;
				core::conformation::Residue const & rsd3 = pose.conformation().residue(j);
				if ( is_sc_bb_clash(rsd1,rsd3) ) {
					//ignore side-chain clash, since it also involves backbone clashes
					//std::cout << "backbone clashes with " << j << std::endl;
					found_sc_only_clash = false;
					break;
				}
			}
		}
		if ( found_sc_only_clash ) clash_positions.push_back(i);
	}
	return clash_positions;
}

ResidueSubset
ClashBasedRepackShellSelector::apply( core::pose::Pose const & pose ) const
{
	ResidueSubset subset( pose.total_residue(), false );
	core::pose::PoseOP mypose( new core::pose::Pose(pose) );

	// determine design positions based on the packer task
	std::set<core::Size> design_shell;
	for ( core::Size i = 1; i <= mypose->total_residue(); ++i ) {
		if ( packer_task_->design_residue(i) ) {
			design_shell.insert(i);
		}
	}

	// determine clash-based repack shell around the design shell
	std::set<core::Size> repack_shell;

	for ( std::set<core::Size>::const_iterator cit = design_shell.begin(); cit != design_shell.end(); ++cit ) {
		core::Size design_pos = *cit;

		// setup a bump check packer task for the given design residue
		core::pack::task::PackerTaskOP bump_check_task( packer_task_->clone() );
		bump_check_task->set_bump_check( false );//don't throw out rotamers during packer graph generation
		bump_check_task->temporarily_fix_everything();//NATRO for all residues
		bump_check_task->temporarily_set_pack_residue( design_pos, true );//set design pos to packable
		score_fxn_->setup_for_packing( *mypose, bump_check_task->repacking_residues(), bump_check_task->designing_residues() );

		// create the rotamer set for the given design residue
		core::conformation::Residue const & design_res = mypose->residue( design_pos );
		core::pack::rotamer_set::RotamerSetFactory rotamer_set_factory;
		core::pack::rotamer_set::RotamerSetOP rotamer_set = rotamer_set_factory.create_rotamer_set( design_res );
		rotamer_set->set_resid( design_pos );
		core::graph::GraphOP packer_graph = core::pack::create_packer_graph( *mypose, *score_fxn_, bump_check_task );
		rotamer_set->build_rotamers( *mypose, *score_fxn_, *bump_check_task, packer_graph );

		// determine neighboring positions where any rotamer of the given design residue would result in a clash
		for ( core::Size i = 1; i <= rotamer_set->num_rotamers(); i++ ) {
			//std::cout << "checking rotamer " << i << std::endl;
			core::conformation::ResidueOP bump_check_rotamer(  rotamer_set->rotamer( i )->clone() );
			//check if the given rotamer would result in a clash with a neighboring side chain
			utility::vector1<core::Size> clash_positions = get_clashing_positions(*mypose, *bump_check_rotamer, design_pos);
			// register all those clash-positions for repacking that are not part of the design shell
			for ( core::Size j = 1; j <= clash_positions.size(); ++j ) {
				core::Size clash_pos = clash_positions[j];
				if ( design_shell.find(clash_pos) == design_shell.end() ) {
					//std::cout << "adding " << clash_pos << std::endl;
					repack_shell.insert(clash_pos);
				}
			}
		}
	}

	for ( core::Size i = 1; i <= mypose->total_residue(); i++ ) {
		if ( repack_shell.find(i) != repack_shell.end() ) {
			subset[ i ] = true;
		}
	}

	return subset;
}

void ClashBasedRepackShellSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
)
{
	bump_overlap_factor_ = tag->getOption< core::Real >( "bump_overlap_factor", 0.5 );
}

std::string ClashBasedRepackShellSelector::get_name() const {
	return ClashBasedRepackShellSelector::class_name();
}

std::string ClashBasedRepackShellSelector::class_name() {
	return "ClashBasedRepackShell";
}

// setters
void ClashBasedRepackShellSelector::set_packer_task( core::pack::task::PackerTaskOP packer_task ) { packer_task_ = packer_task; }
void ClashBasedRepackShellSelector::set_score_fxn( core::scoring::ScoreFunctionOP score_fxn ) { score_fxn_ = score_fxn; }
void ClashBasedRepackShellSelector::set_bump_overlap_factor( core::Real bump_overlap_factor ) { bump_overlap_factor_ = bump_overlap_factor; }

// getters
core::pack::task::PackerTaskOP ClashBasedRepackShellSelector::get_packer_task() const { return packer_task_; }
core::scoring::ScoreFunctionOP ClashBasedRepackShellSelector::get_score_fxn() const { return score_fxn_; }
core::Real ClashBasedRepackShellSelector::get_bump_overlap_factor() const { return bump_overlap_factor_; }

ResidueSelectorOP
ClashBasedRepackShellSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new ClashBasedRepackShellSelector );
}

std::string
ClashBasedRepackShellSelectorCreator::keyname() const {
	return ClashBasedRepackShellSelector::class_name();
}

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core
