// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/protocols/simple_moves/RotamerizeMover.cc
/// @brief Definitions of class methods for RotamerizeMover
/// @author Jim Havranek

// Unit headers
#include <protocols/simple_moves/RotamerizeMover.hh>

#include <basic/datacache/DataMap.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/exit.hh>
//#include <utility/Tag/Tag.hh>
#include <utility/string_util.hh> // string_split

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace basic::options;
using namespace pack;
using namespace rotamer_set;
using namespace task;
using namespace operation;
using namespace scoring;
using namespace chemical;
using namespace conformation;

using basic::Warning;
using basic::t_warning;
static THREAD_LOCAL basic::Tracer TR( "protocols.moves.RotamerizeMover" );

/// RotamerizeMover

RotamerizeMover::RotamerizeMover() :
	Mover("RotamerizeMover"),
	task_(/* 0 */),
	task_factory_(/* 0 */),
	rotamer_sets_( RotamerSetsOP( new rotamer_set::RotamerSets ) )
{}

RotamerizeMover::RotamerizeMover( std::string const & type_name ) :
	Mover( type_name ),
	task_(/* 0 */),
	task_factory_(/* 0 */),
	rotamer_sets_( RotamerSetsOP( new rotamer_set::RotamerSets ) )
{}

// constructors with arguments
RotamerizeMover::RotamerizeMover(
	PackerTaskOP task
) :
	Mover("RotamerizeMover"),
	task_( task ),
	task_factory_(/* 0 */),
	rotamer_sets_( RotamerSetsOP( new rotamer_set::RotamerSets ) )
{}

RotamerizeMover::~RotamerizeMover(){}

RotamerizeMover::RotamerizeMover( RotamerizeMover const & other )
: Mover( other )
{
	task_ = other.task();
	task_factory_ = other.task_factory();
	rotamer_sets_ = RotamerSetsOP( new rotamer_set::RotamerSets );
}

void
RotamerizeMover::apply( Pose & pose )
{
	// get rotamers, energies
	// also performs lazy initialization of ScoreFunction, PackerTask
	this->setup( pose );

	// Now loop through positions and find the best geometric fit, taking into account
	// silly symmetries due to goofy renaming of atoms

	for ( utility::vector1< RotamerSetOP >::const_iterator itr = rotamer_sets_->begin() ;
			itr != rotamer_sets_->end() ; ++itr ) {

		Size this_resid( (*itr)->resid() );

		if ( pose.residue_type( this_resid ).is_NA() ) continue;

		//TR << "Working on RotamerSet for position " << this_resid << std::endl;

		Real curr_best_rmsd( 9999.0 );
		Rotamers::const_iterator best_itr;
		bool best_is_flipped( false );

		// Process the rotamers for this site
		for ( Rotamers::const_iterator r_itr = (*itr)->begin() ; r_itr != (*itr)->end() ; ++r_itr ) {

			ResidueCOP rotamer( *r_itr );

			// Skip if not the same aa type
			if ( pose.residue_type( this_resid ).aa() != rotamer->type().aa() ) {
				TR << "Rotamer residue type doesn't make pose - bad task?" << std::endl;
				continue;
			}

			// Loop over heavy atoms and get the rmsd
			Real rmsd_sum( 0.0 );
			for ( Size iatom = 1 ; iatom <= rotamer->nheavyatoms() ; ++iatom ) {
				rmsd_sum += pose.residue( this_resid ).xyz( iatom ).distance_squared( rotamer->xyz( iatom ) );
			}
			rmsd_sum = std::sqrt( rmsd_sum / (1.0 * rotamer->nheavyatoms() ) );

			if ( rmsd_sum < curr_best_rmsd ) {
				curr_best_rmsd = rmsd_sum;
				best_itr = r_itr;
				best_is_flipped = false;
			}

			// Explicitly check any needed arginine flips
			if ( rotamer->type().aa() == aa_arg ) {

				Real flip_rmsd_sum( 0.0 );
				for ( Size iatom = 1 ; iatom <= rotamer->nheavyatoms() ; ++iatom ) {
					if ( pose.residue( this_resid ).atom_name( iatom ) == "NH1" ) {
						Size other_atom( rotamer->atom_index( "NH2" ) );
						flip_rmsd_sum += pose.residue( this_resid ).xyz( iatom ).distance_squared( rotamer->xyz( other_atom ) );
					} else if ( pose.residue( this_resid ).atom_name( iatom ) == "NH2" ) {
						Size other_atom( rotamer->atom_index( "NH1" ) );
						flip_rmsd_sum += pose.residue( this_resid ).xyz( iatom ).distance_squared( rotamer->xyz( other_atom ) );
					} else {
						flip_rmsd_sum += pose.residue( this_resid ).xyz( iatom ).distance_squared( rotamer->xyz( iatom ) );
					}
				}
				flip_rmsd_sum = std::sqrt( rmsd_sum / (1.0 * rotamer->nheavyatoms() ) );

				if ( flip_rmsd_sum < curr_best_rmsd ) {
					curr_best_rmsd = flip_rmsd_sum;
					best_itr = r_itr;
					best_is_flipped = true;
				}
			} // End of Arg flip check

			if ( rotamer->type().aa() == aa_glu ) {

				// Do this by an explicit chi flip
				Residue mut_res( *rotamer ); // I'm making a copy because I'm going to mess with it

				// Flip about chi 2 by 180 degrees
				mut_res.set_chi( 3, mut_res.chi()[3] - 180.0 );

				Real flip_rmsd_sum( 0.0 );
				for ( Size iatom = 1 ; iatom <= rotamer->nheavyatoms() ; ++iatom ) {
					flip_rmsd_sum += pose.residue( this_resid ).xyz( iatom ).distance_squared( mut_res.xyz( iatom ) );
				}
				flip_rmsd_sum = std::sqrt( rmsd_sum / (1.0 * rotamer->nheavyatoms() ) );

				if ( flip_rmsd_sum < curr_best_rmsd ) {
					curr_best_rmsd = flip_rmsd_sum;
					best_itr = r_itr;
					best_is_flipped = true;
				}
			} // End of Chi 3 flip check

			if ( rotamer->type().aa() == aa_asp || rotamer->type().aa() == aa_phe || rotamer->type().aa() == aa_tyr ) {

				// Do this by an explicit chi flip
				Residue mut_res( *rotamer ); // I'm making a copy because I'm going to mess with it

				// Flip about chi 2 by 180 degrees
				mut_res.set_chi( 2, mut_res.chi()[2] - 180.0 );

				Real flip_rmsd_sum( 0.0 );
				for ( Size iatom = 1 ; iatom <= rotamer->nheavyatoms() ; ++iatom ) {
					flip_rmsd_sum += pose.residue( this_resid ).xyz( iatom ).distance_squared( mut_res.xyz( iatom ) );
				}
				flip_rmsd_sum = std::sqrt( rmsd_sum / (1.0 * rotamer->nheavyatoms() ) );

				if ( flip_rmsd_sum < curr_best_rmsd ) {
					curr_best_rmsd = flip_rmsd_sum;
					best_itr = r_itr;
					best_is_flipped = true;
				}
			} // End of Chi 2 flip check


		}

		// See what we have
		TR << "Best rmsd for this position is " << curr_best_rmsd << " aa type is " << (*best_itr)->type().name3() << std::endl;
		if ( best_is_flipped ) {
			TR << "Best rotamer was a flipped rotamer!" << std::endl;
		}

		// Slam it in there
		pose.replace_residue( this_resid, **best_itr, false );
	}
}

std::string
RotamerizeMover::get_name() const {
	return "RotamerizeMover";
}

/// @brief when the PackerTask was not generated locally, verify compatibility with pose
/// @details the pose residue types must be equivalent to the ones used to generate the ResidueLevelTasks, because of the way that prevent_repacking and its associated flags work
bool
RotamerizeMover::task_is_valid( Pose const & pose ) const
{
	if ( task_->total_residue() != pose.total_residue() ) return false;
	for ( Size i(1); i <= pose.total_residue(); ++i ) {
		chemical::ResidueTypeCOP r( pose.residue_type(i).get_self_ptr() );
		if ( ! task_->residue_task(i).is_original_type( r ) ) return false;
	}
	return true;
}

/// @brief get rotamers, energies. Also performs lazy initialization of ScoreFunction, PackerTask.
void RotamerizeMover::setup( Pose & pose )
{
	// jec update_residue_neighbors() required to update EnergyGraph (ensures graph_state == GOOD) when calling Interface.cc
	pose.update_residue_neighbors();

	core::scoring::ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "empty" );
	core::pack::interaction_graph::InteractionGraphBaseOP ig;

	// if present, task_factory_ always overrides/regenerates task_
	if ( task_factory_ != 0 ) {
		task_ = task_factory_->create_task_and_apply_taskoperations( pose );
	} else if ( task_ == 0 ) {
		Warning() << "undefined PackerTask -- creating a default one" << std::endl;
		task_ = TaskFactory::create_packer_task( pose );
	} else runtime_assert( task_is_valid( pose ) );
	// in case PackerTask was not generated locally, verify compatibility with pose

	// Make sure the bump check is off
	task_->set_bump_check( false );

	pack_rotamers_setup( pose, *scorefxn, task_, rotamer_sets_, ig );

}

// setters
void RotamerizeMover::task( task::PackerTaskOP t ) { task_ = t; }

void RotamerizeMover::task_factory( TaskFactoryCOP tf )
{
	runtime_assert( tf != 0 );
	task_factory_ = tf;
}

// accessors
PackerTaskOP RotamerizeMover::task() const { return task_; }
TaskFactoryCOP RotamerizeMover::task_factory() const { return task_factory_; }
rotamer_set::RotamerSetsCOP RotamerizeMover::rotamer_sets() const { return rotamer_sets_; }
} // moves
} // protocols

