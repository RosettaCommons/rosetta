// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_SugarRotamer.cc
/// @brief Generate sugar pucker rotamers for RNA.
/// @author Fang-Chieh Chou


// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_SugarRotamer.hh>

// Package headers
#include <core/id/TorsionID.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <core/pose/rna/RNA_IdealCoord.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/pose/Pose.hh>

// Numeric headers
#include <numeric/random/random.hh>

#include <basic/Tracer.hh>

using namespace core;
using namespace core::chemical::rna;
using namespace core::pose::rna;
static basic::Tracer TR( "protocols.rotamer_sampler.rna.RNA_SugarRotamer" );
static numeric::random::RandomGenerator RG( 256199 );  // Magic number

namespace protocols {
namespace rotamer_sampler {
namespace rna {

static basic::Tracer TR("protocols.rotamer_sampler.rna.RNA_SugarRotamer");
//////////////////////////////////////////////////////////////////////////////////////////
// Constructor
RNA_SugarRotamer::RNA_SugarRotamer(
	Size const rsd_id,
	PuckerState const pucker_state
):
	RotamerSized(),
	rsd_id_( rsd_id ),
	pucker_state_( pucker_state ),
	skip_same_pucker_( true ),
	idealize_coord_( true ),
	north_pucker_dofs_have_not_been_initialized_( true ),
	south_pucker_dofs_have_not_been_initialized_( true )
{
	runtime_assert( pucker_state_ == ANY_PUCKER || pucker_state_ == NORTH || pucker_state_ == SOUTH );
}
//////////////////////////////////////////////////////////////////////////
void RNA_SugarRotamer::init() {
	if ( pucker_state_ == ANY_PUCKER ) {
		pucker_states_.push_back( NORTH );
		pucker_states_.push_back( SOUTH );
	} else {
		pucker_states_.push_back( pucker_state_ );
	}
	set_init( true );
	reset();
}
//////////////////////////////////////////////////////////////////
void RNA_SugarRotamer::apply( pose::Pose & pose, core::Size const i ) {
	// Arvind - 10/06/2013
	// The new version of this function takes a bunch of code that was previously distributed in a bunch of places like RNA_Util.cc and util.cc, and consolidates it here. This has the additional advantage that the runtime has been dramatically reduced through pre-caching of the DOFs for ideal puckers.
	runtime_assert( is_init() );
	PuckerState pucker_state = pucker_states_[i];
	assert( pucker_state <= 2 );

	static const RNA_IdealCoord ideal_coord;
	static const RNA_FittedTorsionInfo torsion_info;
	Real delta, nu1, nu2;

	PuckerState const curr_pucker = assign_pucker( pose, rsd_id_ );
	if ( skip_same_pucker_ && pucker_state == curr_pucker ) return;

	if ( pucker_state == ANY_PUCKER ) pucker_state = curr_pucker;

	if (idealize_coord_) {
		if ( pucker_state == NORTH ) {
			if (  north_pucker_dofs_have_not_been_initialized_ ) {
				// Arvind - 10/06/2013
				// If it is the first time the RNA_SugarRotamer object has encountered a NORTH pucker, then we need to go through the whole process of running copy_dofs() from the ideal pucker. But when we do this, we are caching all of the dofs we need to alter so that we don't have to do this again
				ideal_coord.apply_pucker(pose, rsd_id_, pucker_state);
				//north_pucker_dof_key_values_ = ideal_coord.apply_and_return(pose, rsd_id_, pucker_state);
				//north_pucker_dofs_have_not_been_initialized_ = false;
			} else {
				//Record the torsions in starting pose
				utility::vector1 < core::id::TorsionID > saved_torsion_id;
				utility::vector1 < Real > saved_torsions;
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_,   id::BB,  ALPHA   ) );
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_,   id::BB,  BETA    ) );
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_,   id::BB,  GAMMA   ) );
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_,   id::BB,  EPSILON ) );
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_,   id::BB,  ZETA    ) );
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_,   id::CHI, 1       ) ); //CHI
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_,   id::CHI, 4       ) ); //O2H
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_-1, id::BB,  ZETA    ) );
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_+1, id::BB,  ALPHA   ) );

				for ( Size index = 1; index <= saved_torsion_id.size(); ++index ) {
					bool const is_exists = ideal_coord.is_torsion_exists( pose, saved_torsion_id[index] );
					if (is_exists) {
						saved_torsions.push_back( pose.torsion( saved_torsion_id[index] ) );
					} else {
						saved_torsions.push_back( -9999 );
					}
				}

				// If we've already cached all the dofs, we can just loop through them directly without using the copy_dofs() machinery.
				for ( std::map < id::DOF_ID , Real >::const_iterator
					 it=north_pucker_dof_key_values_.begin(), it_end = north_pucker_dof_key_values_.end(); it != it_end; ++it ) {
					pose.set_dof( it->first, it->second );
				}

				for ( Size index = 1; index <= saved_torsion_id.size(); ++index ) {
					if ( saved_torsions[index] > -1000 ) pose.set_torsion( saved_torsion_id[index], saved_torsions[index] );
				}
			}
		} else {
			if (  south_pucker_dofs_have_not_been_initialized_ ) {
				// Arvind - 10/06/2013
				// Same procedure for the first occurrence of a SOUTH pucker.
				ideal_coord.apply_pucker(pose, rsd_id_, pucker_state);
				//south_pucker_dof_key_values_ = ideal_coord.apply_and_return(pose, rsd_id_, pucker_state);
				//south_pucker_dofs_have_not_been_initialized_ = false;
			} else {
				//Record the torsions in starting pose
				utility::vector1 < core::id::TorsionID > saved_torsion_id;
				utility::vector1 < Real > saved_torsions;
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_,   id::BB,  ALPHA   ) );
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_,   id::BB,  BETA    ) );
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_,   id::BB,  GAMMA   ) );
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_,   id::BB,  EPSILON ) );
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_,   id::BB,  ZETA    ) );
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_,   id::CHI, 1       ) ); //CHI
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_,   id::CHI, 4       ) ); //O2H
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_-1, id::BB,  ZETA    ) );
				saved_torsion_id.push_back( core::id::TorsionID( rsd_id_+1, id::BB,  ALPHA   ) );

				for ( Size index = 1; index <= saved_torsion_id.size(); ++index ) {
					bool const is_exists = ideal_coord.is_torsion_exists( pose, saved_torsion_id[index] );
					if (is_exists) {
						saved_torsions.push_back( pose.torsion( saved_torsion_id[index] ) );
					} else {
						saved_torsions.push_back( -9999 );
					}
				}

				// If we've already cached all the dofs, we can just loop through them directly without using the copy_dofs() machinery.
				for ( std::map < id::DOF_ID , Real >::const_iterator
					 it=south_pucker_dof_key_values_.begin(), it_end = south_pucker_dof_key_values_.end(); it != it_end; ++it ) {
					pose.set_dof( it->first, it->second );
				}

				for ( Size index = 1; index <= saved_torsion_id.size(); ++index ) {
					if ( saved_torsions[index] > -1000 ) pose.set_torsion( saved_torsion_id[index], saved_torsions[index] );
				}
			}
		}

	} else {
		if (pucker_state == NORTH) {
			delta = torsion_info.delta_north();
			nu2 = torsion_info.nu2_north();
			nu1 = torsion_info.nu1_north();
		} else {
			delta = torsion_info.delta_south();
			nu2 = torsion_info.nu2_south();
			nu1 = torsion_info.nu1_south();
		}
		pose.set_torsion( id::TorsionID( rsd_id_, id::BB,  4 ), delta );
		pose.set_torsion( id::TorsionID( rsd_id_, id::CHI, 2 ), nu2 );
		pose.set_torsion( id::TorsionID( rsd_id_, id::CHI, 3 ), nu1 );
	}
}
//////////////////////////////////////////////////////////////////
}
}
}
