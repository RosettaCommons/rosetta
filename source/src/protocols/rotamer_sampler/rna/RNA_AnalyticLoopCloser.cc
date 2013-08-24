// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_AnalyticLoopCloser.cc
/// @brief Simply close the RNA loop with KIC.
/// @detailed
/// @author Rhiju Das, Fang-Chieh Chou


// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_AnalyticLoopCloser.hh>

// Package headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/pose/Pose.hh>
#include <numeric/random/random.hh>
#include <numeric/angle.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

using namespace core;
using namespace core::pose;

static numeric::random::RandomGenerator RG ( 656840 ); // <- Magic number, do not change it!
using core::id::AtomID;
using core::id::NamedAtomID;
using core::id::DOF_ID;
using core::id::BOGUS_DOF_ID;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using numeric::angle_radians;
using numeric::principal_angle;
using numeric::dihedral_radians;

static basic::Tracer TR( "protocols.rotamer_sampler.rna.RNA_AnalyticLoopCloser" );

namespace protocols {
namespace rotamer_sampler {
namespace rna {

RNA_AnalyticLoopCloser::RNA_AnalyticLoopCloser(
	Pose const & pose,
	Size const moving_suite,
	Size const chainbreak_suite
):
	RotamerSized(),
	verbose_ ( false ),
	moving_suite_ ( moving_suite ),
	chainbreak_suite_ ( chainbreak_suite ),
	ref_pose_ ( pose ),
	nsol_ ( 0 ),
	id_( 0 )
{}

////////////////////////////////////////////////////////////////
void RNA_AnalyticLoopCloser::init() {
	using namespace core::kinematics;
	using namespace numeric::kinematic_closure;
	///// kinematic loop close.
	// Following copied from, e.g., KinematicMover.cc.  Need to elaborate for terminal residues!
	// inputs to loop closure
	utility::vector1< utility::vector1< Real > > atoms;
	utility::vector1< Size > pivots ( 3 ), order ( 3 );
	utility::vector1< Real > dt_ang, db_len, db_ang;

	atom_ids_.clear();
	atom_ids_.push_back ( NamedAtomID ( " C3'", moving_suite_ ) );
	atom_ids_.push_back ( NamedAtomID ( " O3'", moving_suite_ ) );
	atom_ids_.push_back ( NamedAtomID ( " P  ", moving_suite_ + 1 ) );
	atom_ids_.push_back ( NamedAtomID ( " O5'", moving_suite_ + 1 ) );
	atom_ids_.push_back ( NamedAtomID ( " C5'", moving_suite_ + 1 ) );
	atom_ids_.push_back ( NamedAtomID ( " C4'", moving_suite_ + 1 ) );
	atom_ids_.push_back ( NamedAtomID ( " C3'", chainbreak_suite_ ) );
	atom_ids_.push_back ( NamedAtomID ( " O3'", chainbreak_suite_ ) );
	atom_ids_.push_back ( NamedAtomID ( " P  ", chainbreak_suite_ + 1 ) );
	atom_ids_.push_back ( NamedAtomID ( " O5'", chainbreak_suite_ + 1 ) );
	atom_ids_.push_back ( NamedAtomID ( " C5'", chainbreak_suite_ + 1 ) );
	atom_ids_.push_back ( NamedAtomID ( " C4'", chainbreak_suite_ + 1 ) );
	atom_ids_.push_back ( NamedAtomID ( " C3'", chainbreak_suite_ + 1 ) );
	atom_ids_.push_back ( NamedAtomID ( " O3'", chainbreak_suite_ + 1 ) );
	atom_ids_.push_back ( NamedAtomID ( " C2'", chainbreak_suite_ + 1 ) );

	for ( Size i = 1; i <= 3; i++ ) {
		order[i] = i;
		pivots[i] = 3 * i + 2;
	}

	fill_chainTORS ( ref_pose_, atom_ids_, atoms, dt_ang, db_ang, db_len );
	if ( verbose_ ) output_chainTORS ( dt_ang, db_ang, db_len );

	// These atoms_xyz are the same as computed in fill_chainTORS, but I'm
	// having some trouble sending them out (can't clear or copy vector1< Vector > ?)
	utility::vector1< Vector > atoms_xyz;

	for ( Size i = 1; i <= atoms.size(); i++ ) {
		atoms_xyz.push_back ( Vector ( atoms[i][1], atoms[i][2], atoms[i][3] ) );
	}

	//////////////////////////////////////////////
	// Parameter at chainbreak.
	// This looks a bit weird, because of a hack.
	Size cutpos_ = chainbreak_suite_;
	////////////////////////////////////////////////////////////////
	Real const d_O3star_nextP = ( ref_pose_.xyz ( NamedAtomID ( " O3'", cutpos_ ) ) -
	                              ref_pose_.xyz ( NamedAtomID ( "OVL1", cutpos_ ) ) ).length();
	db_len[ 8 ] = d_O3star_nextP;
	////////////////////////////////////////////////////////////////
	Real const theta_C3star_O3star_nextP = degrees ( angle_radians ( ref_pose_.xyz ( NamedAtomID ( " C3'", cutpos_ ) ),
	                                       ref_pose_.xyz ( NamedAtomID ( " O3'", cutpos_ ) ),
	                                       ref_pose_.xyz ( NamedAtomID ( "OVL1", cutpos_ ) ) ) );
	db_ang[ 8 ] = theta_C3star_O3star_nextP;
	Real const theta_O3star_nextP_nextO5star = degrees ( angle_radians ( ref_pose_.xyz ( NamedAtomID ( " O3'", cutpos_ ) ),
	    ref_pose_.xyz ( NamedAtomID ( "OVL1", cutpos_ ) ),
	    ref_pose_.xyz ( NamedAtomID ( "OVL2", cutpos_ ) ) ) );
	db_ang[ 9 ] = theta_O3star_nextP_nextO5star;
	////////////////////////////////////////////////////////////////
	Real const phi_C4star_C3star_O3star_nextP = degrees ( dihedral_radians ( ref_pose_.xyz ( NamedAtomID ( " C4'", moving_suite_ + 1 ) ),
	    ref_pose_.xyz ( NamedAtomID ( " C3'", cutpos_ ) ),
	    ref_pose_.xyz ( NamedAtomID ( " O3'", cutpos_ ) ),
	    ref_pose_.xyz ( NamedAtomID ( "OVL1", cutpos_ ) ) ) );
	dt_ang[ 7 ] =  phi_C4star_C3star_O3star_nextP;
	Real const phi_C3star_O3star_nextP_nextO5star = degrees ( dihedral_radians (
	      ref_pose_.xyz ( NamedAtomID ( " C3'", cutpos_ ) ),
	      ref_pose_.xyz ( NamedAtomID ( " O3'", cutpos_ ) ),
	      ref_pose_.xyz ( NamedAtomID ( "OVL1", cutpos_ ) ),
	      ref_pose_.xyz ( NamedAtomID ( "OVL2", cutpos_ ) ) ) );
	dt_ang[ 8 ] =  phi_C3star_O3star_nextP_nextO5star;
	Real const phi_O3star_nextP_nextO5star_nextC5star = degrees ( dihedral_radians ( ref_pose_.xyz ( NamedAtomID ( "OVU1", cutpos_ + 1 ) ),
	    ref_pose_.xyz ( NamedAtomID ( " P  ", cutpos_ + 1 ) ),
	    ref_pose_.xyz ( NamedAtomID ( " O5'", cutpos_ + 1 ) ),
	    ref_pose_.xyz ( NamedAtomID ( " C5'", cutpos_ + 1 ) ) ) );
	dt_ang[ 9 ] =  phi_O3star_nextP_nextO5star_nextC5star;
	if ( verbose_ ) {
		TR <<  "after chainbreak geometry fix" << std::endl;
		output_chainTORS ( dt_ang, db_ang, db_len );
	}
	///////////////////////////////////
	// Perform loop closure
	///////////////////////////////////
	t_ang_.clear();
	b_ang_.clear();
	b_len_.clear();
	nsol_ = 0;
	bridgeObjects ( atoms, dt_ang, db_ang, db_len, pivots, order, t_ang_, b_ang_, b_len_, nsol_ );
	figure_out_dof_ids_and_offsets ( ref_pose_, dt_ang );
	set_init( true );
	reset();
}
////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::figure_out_dof_ids_and_offsets(
	utility::vector1< Real > const & dt_ang
) {
	////////////////////////////////////////////////////////////////
	// Note that the torsion angles that we solved for do not directly correspond to
	// torsion angles in the atom-tree. But they are right up to an *offset*, which
	// we pre-calculate now. Also this is a good time to figure out exactly which
	// DOFs need to be changed in the atom-tree.
	//
	// Just be totally explicit.
	//
	////////////////////////////////////////////////////////////////
	offset_save_.clear();
	dof_ids_.clear();
	DOF_ID dof_id;
	AtomID id1, id2, id3, id4;
	/////////////////////////////////////////
	// pivot 1
	/////////////////////////////////////////
	id1 = AtomID ( ref_pose_.residue( moving_suite_ + 1 ).atom_index ( " P  " ), moving_suite_ + 1 );
	id2 = AtomID ( ref_pose_.residue( moving_suite_ + 1 ).atom_index ( " O5'" ), moving_suite_ + 1 );
	id3 = AtomID ( ref_pose_.residue( moving_suite_ + 1 ).atom_index ( " C5'" ), moving_suite_ + 1 );
	id4 = AtomID ( ref_pose_.residue( moving_suite_ + 1 ).atom_index ( " C4'" ), moving_suite_ + 1 );
	dof_id = ref_pose_.atom_tree().torsion_angle_dof_id ( id1, id2, id3, id4 );
	figure_out_offset ( ref_pose_, dof_id, dt_ang[ 3 * 1 + 1 ], offset_save_ );
	id1 = AtomID ( ref_pose_.residue( moving_suite_ + 1 ).atom_index ( " O5'" ), moving_suite_ + 1 );
	id2 = AtomID ( ref_pose_.residue( moving_suite_ + 1 ).atom_index ( " C5'" ), moving_suite_ + 1 );
	id3 = AtomID ( ref_pose_.residue( moving_suite_ + 1 ).atom_index ( " C4'" ), moving_suite_ + 1 );
	id4 = AtomID ( ref_pose_.residue( moving_suite_ + 1 ).atom_index ( " C3'" ), moving_suite_ + 1 );
	dof_id = ref_pose_.atom_tree().torsion_angle_dof_id ( id1, id2, id3, id4 );
	figure_out_offset ( ref_pose_, dof_id, dt_ang[ 3 * 1 + 2 ], offset_save_ );
	/////////////////////////////////////////
	// pivot 2
	/////////////////////////////////////////
	id1 = AtomID ( ref_pose_.residue( chainbreak_suite_ ).atom_index ( " C4'" ), chainbreak_suite_ );
	id2 = AtomID ( ref_pose_.residue( chainbreak_suite_ ).atom_index ( " C3'" ), chainbreak_suite_ );
	id3 = AtomID ( ref_pose_.residue( chainbreak_suite_ ).atom_index ( " O3'" ), chainbreak_suite_ );
	id4 = AtomID ( ref_pose_.residue( chainbreak_suite_ ).atom_index ( "OVL1" ), chainbreak_suite_ );
	dof_id = ref_pose_.atom_tree().torsion_angle_dof_id ( id1, id2, id3, id4 );
	figure_out_offset ( ref_pose_, dof_id, dt_ang[ 3 * 2 + 1 ], offset_save_ );
	id1 = AtomID ( ref_pose_.residue( chainbreak_suite_ ).atom_index ( " C3'" ), chainbreak_suite_ );
	id2 = AtomID ( ref_pose_.residue( chainbreak_suite_ ).atom_index ( " O3'" ), chainbreak_suite_ );
	id3 = AtomID ( ref_pose_.residue( chainbreak_suite_ ).atom_index ( "OVL1" ), chainbreak_suite_ );
	id4 = AtomID ( ref_pose_.residue( chainbreak_suite_ ).atom_index ( "OVL2" ), chainbreak_suite_ );
	dof_id = ref_pose_.atom_tree().torsion_angle_dof_id ( id1, id2, id3, id4 );
	figure_out_offset ( ref_pose_, dof_id, dt_ang[ 3 * 2 + 2 ], offset_save_ );
	/////////////////////////////////////////
	// pivot 3
	/////////////////////////////////////////
	id1 = AtomID ( ref_pose_.residue( chainbreak_suite_ + 1 ).atom_index ( " P  " ), chainbreak_suite_ + 1 );
	id2 = AtomID ( ref_pose_.residue( chainbreak_suite_ + 1 ).atom_index ( " O5'" ), chainbreak_suite_ + 1 );
	id3 = AtomID ( ref_pose_.residue( chainbreak_suite_ + 1 ).atom_index ( " C5'" ), chainbreak_suite_ + 1 );
	id4 = AtomID ( ref_pose_.residue( chainbreak_suite_ + 1 ).atom_index ( " C4'" ), chainbreak_suite_ + 1 );
	dof_id = ref_pose_.atom_tree().torsion_angle_dof_id ( id1, id2, id3, id4 );
	figure_out_offset ( ref_pose_, dof_id, dt_ang[ 3 * 3 + 1 ], offset_save_ );
	id1 = AtomID ( ref_pose_.residue( chainbreak_suite_ + 1 ).atom_index ( " O5'" ), chainbreak_suite_ + 1 );
	id2 = AtomID ( ref_pose_.residue( chainbreak_suite_ + 1 ).atom_index ( " C5'" ), chainbreak_suite_ + 1 );
	id3 = AtomID ( ref_pose_.residue( chainbreak_suite_ + 1 ).atom_index ( " C4'" ), chainbreak_suite_ + 1 );
	id4 = AtomID ( ref_pose_.residue( chainbreak_suite_ + 1 ).atom_index ( " C3'" ), chainbreak_suite_ + 1 );
	dof_id = ref_pose_.atom_tree().torsion_angle_dof_id ( id1, id2, id3, id4 );
	figure_out_offset ( ref_pose_, dof_id, dt_ang[ 3 * 3 + 2 ], offset_save_ );
}

////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::figure_out_offset (
	id::DOF_ID const & dof_id,
	Real const & original_torsion_value,
	utility::vector1< Real > & offset_save
) {
	if ( dof_id == BOGUS_DOF_ID ) { //expected at cutpoint!
		utility_exit_with_message ( "Problem with DOF_ID" );
	} else {
		offset_save.push_back ( ref_pose_.dof ( dof_id ) - radians ( original_torsion_value ) );
		dof_ids_.push_back ( dof_id );
		if ( verbose_ ) {
			TR << dof_id;
			TR << "  offset " << ref_pose_.dof ( dof_id ) << " " << radians ( original_torsion_value )
					<< " " << ref_pose_.dof ( dof_id ) - radians ( original_torsion_value ) << std::endl;
		}
	}
}
////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::fill_chainTORS(
  utility::vector1< id::NamedAtomID > const & atom_ids_,
  utility::vector1< utility::vector1< Real > > & atoms,
  utility::vector1< Real > & dt_ang,
  utility::vector1< Real > & db_ang,
  utility::vector1< Real > & db_len
) const {
	using namespace core::kinematics;
	using namespace numeric::kinematic_closure;
	utility::vector1< utility::vector1< Real > > Q0 ( 3 );
	utility::vector1< Real > R0 ( 3 );
	utility::vector1< Vector > atoms_xyz;

	for ( Size i = 1; i <= atom_ids_.size(); i++ ) {
		atoms_xyz.push_back ( ref_pose_.xyz ( atom_ids_[ i ] ) );
	}
	atoms.clear();

	for ( Size i = 1; i <= atom_ids_.size(); i++ ) {
		utility::vector1< Real > atom_xyz_vals;
		atom_xyz_vals.push_back ( atoms_xyz[i].x() );
		atom_xyz_vals.push_back ( atoms_xyz[i].y() );
		atom_xyz_vals.push_back ( atoms_xyz[i].z() );
		atoms.push_back ( atom_xyz_vals );
	}

	chainTORS ( atoms.size(), atoms, dt_ang, db_ang, db_len, R0, Q0 );
}
////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::output_chainTORS (
	utility::vector1< core::Real > const & dt_ang,
	utility::vector1< core::Real > const & db_ang,
	utility::vector1< core::Real > const & db_len )
const {
	TR << "------  chainTORS output ---- " << std::endl;

	for ( Size i = 1; i <= dt_ang.size(); i++ ) {
		TR << I ( 3, i ) << " ";
		TR << "TORSIONS: ";
		TR << F ( 8, 3, dt_ang[ i ] ) << " ";
		TR << "   BOND_ANGLES: ";
		TR << F ( 8, 3, db_ang[ i ] ) << " ";
		TR << "   BOND_LENGTHS: ";
		TR << F ( 8, 3, db_len[ i ] ) << " ";
		TR << std::endl;
	}
}

////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::apply( Pose & pose, Size const id ) const {
	runtime_assert( is_init() );
	Size count ( 0 );
	for ( Size i = 1; i <= 3; i++ ) {
		count++;
		ref_pose_.set_dof ( dof_ids_[count], principal_angle ( radians ( t_ang_[id][3 * i + 1] ) + offset_save_[count] ) );
		count++;
		ref_pose_.set_dof ( dof_ids_[count], principal_angle ( radians ( t_ang_[id][3 * i + 2] ) + offset_save_[count] ) );
	}
}
////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::operator++() {
	runtime_assert( not_end() );
	if ( is_random() )
		id_ = RG.random_range( 1, size() );
	} else {
		++id_;
	}
}
////////////////////////////////////////////////////////////////
bool
RNA_AnalyticLoopCloser::not_end() const {
	runtime_assert( is_init() );
	if ( is_random() ) return true;
	return ( id_ <= size() );
}
////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::reset() {
	runtime_assert( is_init() );
	if ( is_random() ) {
		++( *this );
	} else {
		id_ = 1;
	}
}
////////////////////////////////////////////////////////////////
}
}
}
