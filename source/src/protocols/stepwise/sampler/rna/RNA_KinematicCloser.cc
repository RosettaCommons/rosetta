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
/// @author Kalli Kappel


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

#include <Eigen/Dense>

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

/////////////////////////////////////////////////////////////////////////////////////////////
/// @details
///
/// How about some ASCII art to show all the atoms at play?
///
///     C1'                               C1'             C1'                               C1'
///   /    \     moving_suite           /    \          /    \    chainbreak_suite        /   \     ...
///  O4   C2'                          O4   C2'        O4   C2'                          O4     C2'
///  |     |                           |     |         |     |    OVU | OVL1 OVL2        |     |
///  C4’---C3'---O3'---P---O5'---C5'---C4'--C3'- ... --C4'---C3'---O3'---P---O5'---C5'---C4'--C3'
///   delt   eps   zeta  alph  bet  gam  deta            delt   eps  zeta  alph  bet  gam  delt
///    X     D1    D2     D3    P1--P1    X               X      P2--P2     D4    P3--P3     X
///
///  X = fixed (sugar)
///  P = pivot torsions
///  D = driver torsions
///
/// TODO: Generalize to cyclized cutpoints.
/// TODO: Test that driver torsions vs. pivot torsions agree in cyclic dinucleotide.
/// TODO: Check possibly incorrect O3' atom.
///
/////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

RNA_KinematicCloser::RNA_KinematicCloser(
	Pose const & init_pose,
	Size const moving_suite,
	Size const chainbreak_suite
):
	StepWiseSamplerSized(),
	init_pose_( init_pose ),
	moving_suite_( moving_suite ),
	chainbreak_suite_( chainbreak_suite ),
	verbose_( false ),
	nsol_( 0 ),
	calculate_jacobians_( false )
{}

RNA_KinematicCloser::~RNA_KinematicCloser() {}

////////////////////////////////////////////////////////////////
void RNA_KinematicCloser::init() {
	// the original init(), e.g., used in ERRASER & stepwise KIC moves,
	// just made use of init_pose. No separate "pose_closed",
	init( init_pose_,
				init_pose_,
				false /*use_pose_closed*/ );
}

void RNA_KinematicCloser::init( pose::Pose const & pose_current /* has CUTPOINT variants, may not be closed as long as pose_closed is provided.  */,
																pose::Pose const & pose_closed  /* no CUTPOINT variants, definitely closed,
																																	 in use in detailed_balance mode */,
																bool use_pose_closed /* = true */ )
{
	using namespace core::kinematics;
	using namespace numeric::kinematic_closure;
	///// kinematic loop close.
	// Following copied from, e.g., KinematicMover.cc.  Need to elaborate for terminal residues!
	// inputs to loop closure
	dt_ang_.clear();
	db_len_.clear();
	db_ang_.clear();
	t_ang_.clear();
	b_ang_.clear();
	b_len_.clear();
	offset_save_.clear();
	dof_ids_.clear();
	atom_ids_.clear();
	pivot_ids_.clear();
	atoms_.clear();
	nsol_ = 0;
	all_jacobians_.clear();

	utility::vector1< Size > pivots ( 3 ), order ( 3 );
	for ( Size i = 1; i <= 3; i++ ) {
		order[i] = i;
		pivots[i] = 3 * i + 2;
	}

	////////////////////////////
	// needed for bridgeObjects
	////////////////////////////
	// moving_suite
	atom_ids_.emplace_back( NamedAtomID( " C3'", moving_suite_ ) );
	atom_ids_.emplace_back( NamedAtomID( " O3'", moving_suite_ ) );
	atom_ids_.emplace_back( NamedAtomID( " P  ", moving_suite_ + 1 ) );
	atom_ids_.emplace_back( NamedAtomID( " O5'", moving_suite_ + 1 ) );
	atom_ids_.emplace_back( NamedAtomID( " C5'", moving_suite_ + 1 ) ); // pivot1
	atom_ids_.emplace_back( NamedAtomID( " C4'", moving_suite_ + 1 ) );

	// chainbreak_suite
	atom_ids_.emplace_back( NamedAtomID( " C3'", chainbreak_suite_ ) );
	atom_ids_.emplace_back( NamedAtomID( " O3'", chainbreak_suite_ ) ); // pivot2
	atom_ids_.emplace_back( NamedAtomID( " P  ", chainbreak_suite_ + 1 ) );
	atom_ids_.emplace_back( NamedAtomID( " O5'", chainbreak_suite_ + 1 ) );
	atom_ids_.emplace_back( NamedAtomID( " C5'", chainbreak_suite_ + 1 ) ); // pivot3
	atom_ids_.emplace_back( NamedAtomID( " C4'", chainbreak_suite_ + 1 ) );

	// apparently important to have one more atom
	atom_ids_.emplace_back( NamedAtomID( " C3'", chainbreak_suite_ + 1 ) );

	/////////////////////////
	// needed for jacobians
	/////////////////////////
	// pivot 1
	pivot_ids_.emplace_back( NamedAtomID( "O5'", moving_suite_ + 1 ) );
	pivot_ids_.emplace_back( NamedAtomID( "C5'", moving_suite_ + 1 ) );
	pivot_ids_.emplace_back( NamedAtomID( "C4'", moving_suite_ + 1 ) );

	// pivot 2
	pivot_ids_.emplace_back( NamedAtomID( "C3'", chainbreak_suite_ ) );
	pivot_ids_.emplace_back( NamedAtomID( "O3'", chainbreak_suite_ ) );
	pivot_ids_.emplace_back( NamedAtomID( " P  ", chainbreak_suite_ + 1 ) );

	// pivot 3
	pivot_ids_.emplace_back( NamedAtomID( "O3'", chainbreak_suite_ + 1 ) ); // HEY WAIT! SHOULDNT THIS BE O5' ??? -- rhiju 2016
	pivot_ids_.emplace_back( NamedAtomID( "C5'", chainbreak_suite_ + 1 ) );
	pivot_ids_.emplace_back( NamedAtomID( "C4'", chainbreak_suite_ + 1 ) );

	fill_chainTORS( pose_current ); // fills atoms_, including endpoints of the loop, which need to be preserved.
	if ( verbose_ ) output_chainTORS();

	//////////////////////////////////////////////
	// Parameter at chainbreak.
	// This looks a bit weird, because of a hack.
	Size const cutpos = chainbreak_suite_;

	if ( use_pose_closed ) {
		////////////////////////////////////////////////////////////////
		// Following is apparently necessary when doing detailed-
		// balance sampling, worked out by Kalli Kappel in
		// RNA_KinematicCloser_DB (now unified with RNA_KinematicCloser)
		// -- rhiju, 2016
		////////////////////////////////////////////////////////////////
		Real const d_O3prime_nextP = (
			pose_closed.xyz( NamedAtomID( " O3'", cutpos ) ) -
			pose_closed.xyz( NamedAtomID( " P  ", cutpos + 1 ) ) ).length();
		db_len_[ 8 ] = d_O3prime_nextP;
		////////////////////////////////////////////////////////////////
		Real const theta_C3prime_O3prime_nextP = degrees( angle_radians(
			pose_closed.xyz( NamedAtomID( " C3'", cutpos ) ),
			pose_closed.xyz( NamedAtomID( " O3'", cutpos ) ),
			pose_closed.xyz( NamedAtomID( " P  ", cutpos + 1 ) ) ) );
		db_ang_[ 8 ] = theta_C3prime_O3prime_nextP;
		Real const theta_O3prime_nextP_nextO5prime = degrees( angle_radians(
			pose_closed.xyz( NamedAtomID( " O3'", cutpos ) ),
			pose_closed.xyz( NamedAtomID( " P  ", cutpos + 1 ) ),
			pose_closed.xyz( NamedAtomID( " O5'", cutpos + 1 ) ) ) );
		db_ang_[ 9 ] = theta_O3prime_nextP_nextO5prime;
	} else {
		////////////////////////////////////////////////////////////////
		// Original code took geometry based on OVL atoms...
		////////////////////////////////////////////////////////////////
		Real const d_O3prime_nextP = (
		  pose_current.xyz( NamedAtomID( " O3'", cutpos ) ) -
		  pose_current.xyz( NamedAtomID( "OVL1", cutpos ) ) ).length();
		db_len_[ 8 ] = d_O3prime_nextP;
		////////////////////////////////////////////////////////////////
		Real const theta_C3prime_O3prime_nextP = degrees( angle_radians(
		  pose_current.xyz( NamedAtomID( " C3'", cutpos ) ),
		  pose_current.xyz( NamedAtomID( " O3'", cutpos ) ),
		  pose_current.xyz( NamedAtomID( "OVL1", cutpos ) ) ) );
		db_ang_[ 8 ] = theta_C3prime_O3prime_nextP;
	  Real const theta_O3prime_nextP_nextO5prime = degrees( angle_radians(
		  pose_current.xyz( NamedAtomID( " O3'", cutpos ) ),
			pose_current.xyz( NamedAtomID( "OVL1", cutpos ) ),
			pose_current.xyz( NamedAtomID( "OVL2", cutpos ) ) ) );
		db_ang_[ 9 ] = theta_O3prime_nextP_nextO5prime;
	}

	////////////////////////////////////////////////////////////////
	Real const phi_C4prime_C3prime_O3prime_nextP = degrees( dihedral_radians(
		pose_current.xyz( NamedAtomID( " C4'", moving_suite_ + 1 ) ),
		pose_current.xyz( NamedAtomID( " C3'", cutpos ) ),
		pose_current.xyz( NamedAtomID( " O3'", cutpos ) ),
		pose_current.xyz( NamedAtomID( "OVL1", cutpos ) ) ) );
	dt_ang_[ 7 ] =  phi_C4prime_C3prime_O3prime_nextP;
	Real const phi_C3prime_O3prime_nextP_nextO5prime = degrees( dihedral_radians(
		pose_current.xyz( NamedAtomID( " C3'", cutpos ) ),
		pose_current.xyz( NamedAtomID( " O3'", cutpos ) ),
		pose_current.xyz( NamedAtomID( "OVL1", cutpos ) ),
		pose_current.xyz( NamedAtomID( "OVL2", cutpos ) ) ) );
	dt_ang_[ 8 ] =  phi_C3prime_O3prime_nextP_nextO5prime;
	Real const phi_O3prime_nextP_nextO5prime_nextC5prime =
		degrees( dihedral_radians(
		pose_current.xyz( NamedAtomID( "OVU1", cutpos + 1 ) ),
		pose_current.xyz( NamedAtomID( " P  ", cutpos + 1 ) ),
		pose_current.xyz( NamedAtomID( " O5'", cutpos + 1 ) ),
		pose_current.xyz( NamedAtomID( " C5'", cutpos + 1 ) ) ) );
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
	figure_out_dof_ids_and_offsets( pose_current );
  set_init( true );

  /// Calculate the Jacobian for each of the solutions ///
  if ( calculate_jacobians_ ) {
		Pose test_pose = pose_current;
		for ( Size i = 1; i <= nsol_; i++ ) {
			apply( test_pose, i );
			all_jacobians_.push_back( get_jacobian( test_pose ) );
		}
	}

  reset();
}
////////////////////////////////////////////////////////////////
void
RNA_KinematicCloser::figure_out_dof_ids_and_offsets( pose::Pose const & pose ) {
	////////////////////////////////////////////////////////////////
	// Note that the torsion angles that we solved for do not directly correspond to
	// torsion angles in the atom-tree. But they are right up to an *offset*, which
	// we pre-calculate now. Also this is a good time to figure out exactly which
	// DOFs need to be changed in the atom-tree.
	//
	// Just be totally explicit.
	//
	////////////////////////////////////////////////////////////////

	/////////////////////////////////////////
	// pivot 1
	/////////////////////////////////////////
	figure_out_offset(
		pose.atom_tree().torsion_angle_dof_id(
		AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().p_atom_index(), moving_suite_ + 1 ),
		AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().o5prime_atom_index(), moving_suite_ + 1 ),
		AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().c5prime_atom_index(), moving_suite_ + 1 ),
		AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().c4prime_atom_index(), moving_suite_ + 1 ) ),
		dt_ang_[3 * 1 + 1] , pose);

	figure_out_offset(
		pose.atom_tree().torsion_angle_dof_id(
		AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().o5prime_atom_index(), moving_suite_ + 1 ),
		AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().c5prime_atom_index(), moving_suite_ + 1 ),
		AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().c4prime_atom_index(), moving_suite_ + 1 ),
		AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().c3prime_atom_index(), moving_suite_ + 1 ) ),
		dt_ang_[3 * 1 + 2] , pose);

	/////////////////////////////////////////
	// pivot 2
	/////////////////////////////////////////
	figure_out_offset(
		pose.atom_tree().torsion_angle_dof_id(
		AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().c4prime_atom_index(), chainbreak_suite_ ),
		AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().c3prime_atom_index(), chainbreak_suite_ ),
		AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().o3prime_atom_index(), chainbreak_suite_ ),
		AtomID( pose.residue_type( chainbreak_suite_ ).atom_index( "OVL1" ), chainbreak_suite_ ) ),
		dt_ang_[3 * 2 + 1] , pose);
	figure_out_offset(
		pose.atom_tree().torsion_angle_dof_id(
		AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().c3prime_atom_index(), chainbreak_suite_ ),
		AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().o3prime_atom_index(), chainbreak_suite_ ),
		AtomID( pose.residue_type( chainbreak_suite_ ).atom_index( "OVL1" ), chainbreak_suite_ ),
		AtomID( pose.residue_type( chainbreak_suite_ ).atom_index( "OVL2" ), chainbreak_suite_ ) ),
		dt_ang_[3 * 2 + 2], pose);
	/////////////////////////////////////////
	// pivot 3
	/////////////////////////////////////////
	figure_out_offset(
		pose.atom_tree().torsion_angle_dof_id(
		AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().p_atom_index(), chainbreak_suite_ + 1 ),
		AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().o5prime_atom_index(), chainbreak_suite_ + 1 ),
		AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c5prime_atom_index(), chainbreak_suite_ + 1 ),
		AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c4prime_atom_index(), chainbreak_suite_ + 1 ) ),
		dt_ang_[ 3 * 3 + 1 ], pose);
	figure_out_offset(
		pose.atom_tree().torsion_angle_dof_id(
		AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().o5prime_atom_index(), chainbreak_suite_ + 1 ),
		AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c5prime_atom_index(), chainbreak_suite_ + 1 ),
		AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c4prime_atom_index(), chainbreak_suite_ + 1 ),
		AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c3prime_atom_index(), chainbreak_suite_ + 1 ) ),
		dt_ang_[3 * 3 + 2], pose );
}

////////////////////////////////////////////////////////////////
void
RNA_KinematicCloser::figure_out_offset(
	id::DOF_ID const & dof_id,
	Real const original_torsion_value,
	pose::Pose const & pose
) {
	if ( dof_id == BOGUS_DOF_ID ) { //expected at cutpoint!
		utility_exit_with_message( "Problem with DOF_ID" );
	} else {
		offset_save_.push_back( pose.dof ( dof_id ) -
			radians( original_torsion_value ) );
		dof_ids_.push_back ( dof_id );
		if ( verbose_ ) {
			TR << dof_id;
			TR << "  offset " << pose.dof( dof_id ) << " " <<
				radians( original_torsion_value ) << " " <<
				pose.dof( dof_id ) - radians( original_torsion_value ) <<
				std::endl;
		}
	}
}
////////////////////////////////////////////////////////////////
void
RNA_KinematicCloser::fill_chainTORS( pose::Pose const & pose ) {
	using namespace core::kinematics;
	using namespace numeric::kinematic_closure;
	utility::vector1< utility::vector1< Real > > Q0( 3 );
	utility::vector1< Real > R0( 3 );

	for ( Size i = 1; i <= atom_ids_.size(); i++ ) {
		Vector atom_xyz( pose.xyz( atom_ids_[i] ) );
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

////////////////////////////////////////////////////////////////
Real
RNA_KinematicCloser::get_jacobian( Pose & pose ) const {
	using std::abs;
	using numeric::max;

	Eigen::Matrix<Real, 6, 3> r1, r2, s1, s2, gamma;
	Eigen::Matrix<Real, 3, 1> cross_i, cross_45, delta;
	Eigen::Matrix<Real, 4, 4> J;

	// Convert the solution from internal coordinates to  cartesian ones.
	// Since the jacobian is invariant under rigid-body transformations, the
	// origin and reference frame used in this conversion are unimportant.


	for ( Size i = 1; i <= 3; i++ ) {
		//Pivot indices, there are 3
		Size j = 2 * (i - 1);

		//x,y,z, coords of atom before pivot, pivot, and atom after pivot
		r1(j, 0) = pose.xyz( pivot_ids_[3*i-2] ).x();
		r1(j, 1) = pose.xyz( pivot_ids_[3*i-2] ).y();
		r1(j, 2) = pose.xyz( pivot_ids_[3*i-2] ).z();

		r1(j+1, 0) = pose.xyz( pivot_ids_[3*i-1] ).x();
		r1(j+1, 1) = pose.xyz( pivot_ids_[3*i-1] ).y();
		r1(j+1, 2) = pose.xyz( pivot_ids_[3*i-1] ).z();

		r2(j, 0) = pose.xyz( pivot_ids_[3*i-1] ).x();
		r2(j, 1) = pose.xyz( pivot_ids_[3*i-1] ).y();
		r2(j, 2) = pose.xyz( pivot_ids_[3*i-1] ).z();

		r2(j+1, 0) = pose.xyz( pivot_ids_[3*i] ).x();
		r2(j+1, 1) = pose.xyz( pivot_ids_[3*i] ).y();
		r2(j+1, 2) = pose.xyz( pivot_ids_[3*i] ).z();
	}

	// Calculate the jacobian following the method outlined by Nilmeier, Hua,
	// Coutsias, and Jacobson in their 2011 JCTC paper.

	for ( Size i = 0; i < 6; i++ ) {
		delta = r2.row(i) - r1.row(i);
		gamma.row(i) = delta.normalized();
	}

	cross_45 = gamma.row(4).cross(gamma.row(5));

	for ( Size i = 0; i < 4; i++ ) {
		cross_i = gamma.row(i).cross(r1.row(5) - r1.row(i));
		Real dot_i = gamma.row(i).dot(cross_45);

		J(0, i) = cross_i(0);
		J(1, i) = cross_i(1);
		J(2, i) = cross_i(2);
		J(3, i) = dot_i;
	}

	Real jacobian = 1 / max(abs(J.determinant()), 1e-100);

	return jacobian;
}

} //rna
} //sampler
} //stepwise
} //protocols
