// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/rna/RNA_KinematicCloser.cc
/// @brief Close a RNA loop with Kinematic Closer (KIC).
/// @author Rhiju Das
/// @author Fang-Chieh Chou


// Unit headers
#include <protocols/stepwise/sampler/rna/RNA_KinematicCloser.hh>

// Package headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/kinematics/AtomTree.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/angle.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
#include <numeric/random/random.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

using namespace core;
using namespace core::pose;

using core::id::AtomID;
using core::id::NamedAtomID;
using core::id::DOF_ID;
using core::id::BOGUS_DOF_ID;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using numeric::angle_radians;
using numeric::principal_angle;
using numeric::dihedral_radians;

static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.rna.RNA_KinematicCloser" );

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

RNA_KinematicCloser::RNA_KinematicCloser(
	PoseOP const & ref_pose,
	Size const moving_suite,
	Size const chainbreak_suite
):
	StepWiseSamplerSized(),
	ref_pose_( ref_pose ),
	moving_suite_( moving_suite ),
	chainbreak_suite_( chainbreak_suite ),
	verbose_( false ),
	nsol_( 0 )
{}

RNA_KinematicCloser::~RNA_KinematicCloser() {}
////////////////////////////////////////////////////////////////
void RNA_KinematicCloser::init() {
	using namespace core::kinematics;
	using namespace numeric::kinematic_closure;
	///// kinematic loop close.
	// Following copied from, e.g., KinematicMover.cc.  Need to elaborate for terminal residues!
	// inputs to loop closure
	Pose const & pose ( *ref_pose_ );
	dt_ang_.clear();
	db_len_.clear();
	db_ang_.clear();
	t_ang_.clear();
	b_ang_.clear();
	b_len_.clear();
	offset_save_.clear();
	dof_ids_.clear();
	atom_ids_.clear();
	atoms_.clear();
	nsol_ = 0;

	utility::vector1< Size > pivots ( 3 ), order ( 3 );
	for ( Size i = 1; i <= 3; i++ ) {
		order[i] = i;
		pivots[i] = 3 * i + 2;
	}

	atom_ids_.emplace_back( NamedAtomID( " C3'", moving_suite_ ) );
	atom_ids_.emplace_back( NamedAtomID( " O3'", moving_suite_ ) );
	atom_ids_.emplace_back( NamedAtomID( " P  ", moving_suite_ + 1 ) );
	atom_ids_.emplace_back( NamedAtomID( " O5'", moving_suite_ + 1 ) );
	atom_ids_.emplace_back( NamedAtomID( " C5'", moving_suite_ + 1 ) );
	atom_ids_.emplace_back( NamedAtomID( " C4'", moving_suite_ + 1 ) );
	atom_ids_.emplace_back( NamedAtomID( " C3'", chainbreak_suite_ ) );
	atom_ids_.emplace_back( NamedAtomID( " O3'", chainbreak_suite_ ) );
	atom_ids_.emplace_back( NamedAtomID( " P  ", chainbreak_suite_ + 1 ) );
	atom_ids_.emplace_back( NamedAtomID( " O5'", chainbreak_suite_ + 1 ) );
	atom_ids_.emplace_back( NamedAtomID( " C5'", chainbreak_suite_ + 1 ) );
	atom_ids_.emplace_back( NamedAtomID( " C4'", chainbreak_suite_ + 1 ) );
	atom_ids_.emplace_back( NamedAtomID( " C3'", chainbreak_suite_ + 1 ) );
	atom_ids_.emplace_back( NamedAtomID( " O3'", chainbreak_suite_ + 1 ) );
	atom_ids_.emplace_back( NamedAtomID( " C2'", chainbreak_suite_ + 1 ) );

	fill_chainTORS();
	if ( verbose_ ) output_chainTORS();

	//////////////////////////////////////////////
	// Parameter at chainbreak.
	// This looks a bit weird, because of a hack.
	Size const cutpos = chainbreak_suite_;
	////////////////////////////////////////////////////////////////
	Real const d_O3prime_nextP = (
		pose.xyz( NamedAtomID( " O3'", cutpos ) ) -
		pose.xyz( NamedAtomID( "OVL1", cutpos ) ) ).length();
	db_len_[ 8 ] = d_O3prime_nextP;
	////////////////////////////////////////////////////////////////
	Real const theta_C3prime_O3prime_nextP = degrees( angle_radians(
		pose.xyz( NamedAtomID( " C3'", cutpos ) ),
		pose.xyz( NamedAtomID( " O3'", cutpos ) ),
		pose.xyz( NamedAtomID( "OVL1", cutpos ) ) ) );
	db_ang_[ 8 ] = theta_C3prime_O3prime_nextP;
	Real const theta_O3prime_nextP_nextO5prime = degrees( angle_radians(
		pose.xyz( NamedAtomID( " O3'", cutpos ) ),
		pose.xyz( NamedAtomID( "OVL1", cutpos ) ),
		pose.xyz( NamedAtomID( "OVL2", cutpos ) ) ) );
	db_ang_[ 9 ] = theta_O3prime_nextP_nextO5prime;
	////////////////////////////////////////////////////////////////
	Real const phi_C4prime_C3prime_O3prime_nextP = degrees( dihedral_radians(
		pose.xyz( NamedAtomID( " C4'", moving_suite_ + 1 ) ),
		pose.xyz( NamedAtomID( " C3'", cutpos ) ),
		pose.xyz( NamedAtomID( " O3'", cutpos ) ),
		pose.xyz( NamedAtomID( "OVL1", cutpos ) ) ) );
	dt_ang_[ 7 ] =  phi_C4prime_C3prime_O3prime_nextP;
	Real const phi_C3prime_O3prime_nextP_nextO5prime = degrees( dihedral_radians(
		pose.xyz( NamedAtomID( " C3'", cutpos ) ),
		pose.xyz( NamedAtomID( " O3'", cutpos ) ),
		pose.xyz( NamedAtomID( "OVL1", cutpos ) ),
		pose.xyz( NamedAtomID( "OVL2", cutpos ) ) ) );
	dt_ang_[ 8 ] =  phi_C3prime_O3prime_nextP_nextO5prime;
	Real const phi_O3prime_nextP_nextO5prime_nextC5prime =
		degrees( dihedral_radians(
		pose.xyz( NamedAtomID( "OVU1", cutpos + 1 ) ),
		pose.xyz( NamedAtomID( " P  ", cutpos + 1 ) ),
		pose.xyz( NamedAtomID( " O5'", cutpos + 1 ) ),
		pose.xyz( NamedAtomID( " C5'", cutpos + 1 ) ) ) );
	dt_ang_[ 9 ] =  phi_O3prime_nextP_nextO5prime_nextC5prime;
	if ( verbose_ ) {
		TR <<  "after chainbreak geometry fix" << std::endl;
		output_chainTORS();
	}
	///////////////////////////////////
	// Perform loop closure
	///////////////////////////////////
	int nsol = 0; //Somehow bridgeObjects takes a reference to int... weird
	bridgeObjects( atoms_, dt_ang_, db_ang_, db_len_, pivots, order,
		t_ang_, b_ang_, b_len_, nsol );
	runtime_assert( nsol >= 0 ); //I think nsol should be >= 0
	nsol_ = static_cast<Size>( nsol ); //Cast int to Size.
	figure_out_dof_ids_and_offsets();
	set_init( true );
	reset();
}
////////////////////////////////////////////////////////////////
void
RNA_KinematicCloser::figure_out_dof_ids_and_offsets() {
	////////////////////////////////////////////////////////////////
	// Note that the torsion angles that we solved for do not directly correspond to
	// torsion angles in the atom-tree. But they are right up to an *offset*, which
	// we pre-calculate now. Also this is a good time to figure out exactly which
	// DOFs need to be changed in the atom-tree.
	//
	// Just be totally explicit.
	//
	////////////////////////////////////////////////////////////////
	Pose const & pose ( *ref_pose_ );
	DOF_ID dof_id;
	AtomID id1, id2, id3, id4;
	/////////////////////////////////////////
	// pivot 1
	/////////////////////////////////////////
	id1 = AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().p_atom_index(), moving_suite_ + 1 );
	id2 = AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().o5prime_atom_index(), moving_suite_ + 1 );
	id3 = AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().c5prime_atom_index(), moving_suite_ + 1 );
	id4 = AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().c4prime_atom_index(), moving_suite_ + 1 );
	dof_id = pose.atom_tree().torsion_angle_dof_id( id1, id2, id3, id4 );
	figure_out_offset( dof_id, dt_ang_[3 * 1 + 1] );
	id1 = AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().o5prime_atom_index(), moving_suite_ + 1 );
	id2 = AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().c5prime_atom_index(), moving_suite_ + 1 );
	id3 = AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().c4prime_atom_index(), moving_suite_ + 1 );
	id4 = AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().c3prime_atom_index(), moving_suite_ + 1 );
	dof_id = pose.atom_tree().torsion_angle_dof_id( id1, id2, id3, id4 );
	figure_out_offset( dof_id, dt_ang_[3 * 1 + 2] );
	/////////////////////////////////////////
	// pivot 2
	/////////////////////////////////////////

	id1 = AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().c4prime_atom_index(), chainbreak_suite_ );
	id2 = AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().c3prime_atom_index(), chainbreak_suite_ );
	id3 = AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().o3prime_atom_index(), chainbreak_suite_ );
	id4 = AtomID( pose.residue( chainbreak_suite_ ).atom_index( "OVL1" ),
		chainbreak_suite_ );
	dof_id = pose.atom_tree().torsion_angle_dof_id( id1, id2, id3, id4 );
	figure_out_offset( dof_id, dt_ang_[3 * 2 + 1] );

	id1 = AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().c3prime_atom_index(), chainbreak_suite_ );
	id2 = AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().o3prime_atom_index(), chainbreak_suite_ );
	id3 = AtomID ( pose.residue( chainbreak_suite_ ).atom_index( "OVL1" ),
		chainbreak_suite_ );
	id4 = AtomID( pose.residue( chainbreak_suite_ ).atom_index( "OVL2" ),
		chainbreak_suite_ );
	dof_id = pose.atom_tree().torsion_angle_dof_id( id1, id2, id3, id4 );
	figure_out_offset( dof_id, dt_ang_[3 * 2 + 2] );
	/////////////////////////////////////////
	// pivot 3
	/////////////////////////////////////////
	id1 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().p_atom_index(), chainbreak_suite_ + 1 );
	id2 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().o5prime_atom_index(), chainbreak_suite_ + 1 );
	id3 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c5prime_atom_index(), chainbreak_suite_ + 1 );
	id4 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c4prime_atom_index(), chainbreak_suite_ + 1 );
	dof_id = pose.atom_tree().torsion_angle_dof_id ( id1, id2, id3, id4 );
	figure_out_offset( dof_id, dt_ang_[ 3 * 3 + 1 ] );
	id1 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().o5prime_atom_index(), chainbreak_suite_ + 1 );
	id2 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c5prime_atom_index(), chainbreak_suite_ + 1 );
	id3 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c4prime_atom_index(), chainbreak_suite_ + 1 );
	id4 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c3prime_atom_index(), chainbreak_suite_ + 1 );
	dof_id = pose.atom_tree().torsion_angle_dof_id( id1, id2, id3, id4 );
	figure_out_offset( dof_id, dt_ang_[3 * 3 + 2] );
}

////////////////////////////////////////////////////////////////
void
RNA_KinematicCloser::figure_out_offset(
	id::DOF_ID const & dof_id,
	Real const original_torsion_value
) {
	if ( dof_id == BOGUS_DOF_ID ) { //expected at cutpoint!
		utility_exit_with_message( "Problem with DOF_ID" );
	} else {
		offset_save_.push_back( ref_pose_->dof ( dof_id ) -
			radians( original_torsion_value ) );
		dof_ids_.push_back ( dof_id );
		if ( verbose_ ) {
			TR << dof_id;
			TR << "  offset " << ref_pose_->dof( dof_id ) << " " <<
				radians( original_torsion_value ) << " " <<
				ref_pose_->dof( dof_id ) - radians( original_torsion_value ) <<
				std::endl;
		}
	}
}
////////////////////////////////////////////////////////////////
void
RNA_KinematicCloser::fill_chainTORS() {
	using namespace core::kinematics;
	using namespace numeric::kinematic_closure;
	utility::vector1< utility::vector1< Real > > Q0( 3 );
	utility::vector1< Real > R0( 3 );

	for ( Size i = 1; i <= atom_ids_.size(); i++ ) {
		Vector atom_xyz( ref_pose_->xyz( atom_ids_[i] ) );
		utility::vector1< Real > atom_xyzs;
		atom_xyzs.push_back( atom_xyz.x() );
		atom_xyzs.push_back( atom_xyz.y() );
		atom_xyzs.push_back( atom_xyz.z() );
		atoms_.push_back ( atom_xyzs );
	}
	chainTORS( atoms_.size(), atoms_, dt_ang_, db_ang_, db_len_, R0, Q0 );
}
////////////////////////////////////////////////////////////////
void
RNA_KinematicCloser::output_chainTORS() const {
	using namespace ObjexxFCL::format;
	TR << "------  chainTORS output ---- " << std::endl;
	for ( Size i = 1; i <= dt_ang_.size(); i++ ) {
		TR << I( 3, i ) << " ";
		TR << "TORSIONS: ";
		TR << F( 8, 3, dt_ang_[i] ) << " ";
		TR << "   BOND_ANGLES: ";
		TR << F( 8, 3, db_ang_[i] ) << " ";
		TR << "   BOND_LENGTHS: ";
		TR << F( 8, 3, db_len_[i] ) << " ";
		TR << std::endl;
	}
}
////////////////////////////////////////////////////////////////
void
RNA_KinematicCloser::apply( Pose & pose, Size const id ) {
	runtime_assert( is_init() );
	runtime_assert( id > 0 );
	runtime_assert( id <= nsol_ );

	Size count ( 0 );
	for ( Size i = 1; i <= 3; i++ ) {
		++count;
		pose.set_dof( dof_ids_[count], principal_angle( radians (
			t_ang_[id][3 * i + 1] ) + offset_save_[count] ) );
		++count;
		pose.set_dof( dof_ids_[count], principal_angle( radians (
			t_ang_[id][3 * i + 2] ) + offset_save_[count] ) );
	}
}

} //rna
} //sampler
} //stepwise
} //protocols
