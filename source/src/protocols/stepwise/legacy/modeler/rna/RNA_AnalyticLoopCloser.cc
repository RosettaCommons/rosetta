// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RNA_AnalyticLoopCloser
/// @brief protocols that are specific to RNA_AnalyticLoopCloser
/// @details
/// @author Rhiju Das

// Unit headers
#include <protocols/stepwise/legacy/modeler/rna/RNA_AnalyticLoopCloser.hh>

// Package headers
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/random/random.hh>
#include <numeric/angle.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

///////////////////////////////////////////////////////////////////////////
// This should be deprecated in favor of Fang's RNA_KinematicCloser.cc
// which is in protocols/stepwise/sampler/rna/
///////////////////////////////////////////////////////////////////////////

using namespace core;

using ObjexxFCL::format::I;
using ObjexxFCL::format::F;
using core::id::AtomID;
using core::id::NamedAtomID;
using core::id::DOF_ID;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using numeric::angle_radians;
using numeric::principal_angle;
using numeric::dihedral_radians;

static basic::Tracer TR( "protocols.stepwise.legacy.modeler.rna.RNA_AnalyticLoopCloser" );

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace rna {

RNA_AnalyticLoopCloser::RNA_AnalyticLoopCloser ( Size const moving_suite, Size const chainbreak_suite ) :
	moving_suite_ ( moving_suite ),
	chainbreak_suite_ ( chainbreak_suite ),
	verbose_ ( false ),
	nsol_ ( 0 ),
	choose_least_perturb_solution_ ( true ),
	choose_best_solution_ ( false ),
	choose_random_solution_ ( false ),
	save_all_solutions_ ( false ) {
	Mover::type ( "RNA_AnalyticLoopCloser" );
}

/// @brief Clone this object
protocols::moves::MoverOP RNA_AnalyticLoopCloser::clone() const {
	return protocols::moves::MoverOP( new RNA_AnalyticLoopCloser ( *this ) );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details  Apply the RNA Loop Closer -- set up for Full Atom Representation.
///
void
RNA_AnalyticLoopCloser::apply ( core::pose::Pose & pose ) {
	close_at_cutpoint ( pose );
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_AnalyticLoopCloser::close_at_cutpoint ( core::pose::Pose & pose ) {
	using namespace core::kinematics;
	using namespace numeric::kinematic_closure;
	///// kinematic loop close.
	// Following copied from, e.g., KinematicMover.cc.  Need to elaborate for terminal residues!
	// inputs to loop closure
	utility::vector1< utility::fixedsizearray1< Real,3 > > atoms;
	utility::vector1< Size > pivots ( 3 ), order ( 3 );
	utility::vector1< Real > dt_ang, db_len, db_ang;
	// doesn't matter.
	order[1] = 1;
	order[2] = 2;
	order[3] = 3;

	// Set up your NamedAtomIDs
	atom_ids_.clear();
	atom_ids_.emplace_back( " C3'", moving_suite_ );
	atom_ids_.emplace_back( " O3'", moving_suite_ );
	atom_ids_.emplace_back( " P  ", moving_suite_ + 1 );
	atom_ids_.emplace_back( " O5'", moving_suite_ + 1 );
	atom_ids_.emplace_back( " C5'", moving_suite_ + 1 );
	atom_ids_.emplace_back( " C4'", moving_suite_ + 1 );
	atom_ids_.emplace_back( " C3'", chainbreak_suite_ );
	atom_ids_.emplace_back( " O3'", chainbreak_suite_ );
	atom_ids_.emplace_back( " P  ", chainbreak_suite_ + 1 );
	atom_ids_.emplace_back( " O5'", chainbreak_suite_ + 1 );
	atom_ids_.emplace_back( " C5'", chainbreak_suite_ + 1 );
	atom_ids_.emplace_back( " C4'", chainbreak_suite_ + 1 );
	atom_ids_.emplace_back( " C3'", chainbreak_suite_ + 1 );
	atom_ids_.emplace_back( " O3'", chainbreak_suite_ + 1 );
	atom_ids_.emplace_back( " C2'", chainbreak_suite_ + 1 );

	for ( Size i = 1; i <= 3; i++ ) {
		pivots[ i ] = 3 * i + 2;
	}

	fill_chainTORS ( pose, atom_ids_, atoms, dt_ang, db_ang, db_len );

	if ( verbose_ ) output_chainTORS ( dt_ang, db_ang, db_len );

	// These atoms_xyz are the same as computed in fill_chainTORS, but I'm
	// having some trouble sending them out (can't clear or copy vector1< Vector > ?)
	utility::vector1< Vector > atoms_xyz( atoms.size() );
	std::transform( atoms.begin(), atoms.end(),
		atoms_xyz.begin(),
		[&]( utility::fixedsizearray1< Real,3 > const & atom ) {
			return Vector( atom[1], atom[2], atom[3] );
		} );

	//////////////////////////////////////////////
	// Parameter at chainbreak.
	// This looks a bit weird, because of a hack.
	Size cutpos_ = chainbreak_suite_;
	////////////////////////////////////////////////////////////////////////////////////
	Real const d_O3prime_nextP = ( pose.xyz( NamedAtomID( " O3'", cutpos_ ) ) -
		pose.xyz( NamedAtomID( "OVL1", cutpos_ ) ) ).length();
	db_len[ 8 ] = d_O3prime_nextP;
	////////////////////////////////////////////////////////////////////////////////////
	Real const theta_C3prime_O3prime_nextP = degrees( angle_radians( pose.xyz( NamedAtomID( " C3'", cutpos_ ) ),
		pose.xyz( NamedAtomID( " O3'", cutpos_ ) ),
		pose.xyz( NamedAtomID( "OVL1", cutpos_ ) ) ) );
	db_ang[ 8 ] = theta_C3prime_O3prime_nextP;
	Real const theta_O3prime_nextP_nextO5prime = degrees( angle_radians( pose.xyz( NamedAtomID( " O3'", cutpos_ ) ),
		pose.xyz( NamedAtomID( "OVL1", cutpos_ ) ),
		pose.xyz( NamedAtomID( "OVL2", cutpos_ ) ) ) );
	db_ang[ 9 ] = theta_O3prime_nextP_nextO5prime;
	////////////////////////////////////////////////////////////////////////////////////
	Real const phi_C4prime_C3prime_O3prime_nextP = degrees( dihedral_radians( pose.xyz( NamedAtomID( " C4'", moving_suite_ + 1 ) ),
		pose.xyz( NamedAtomID( " C3'", cutpos_ ) ),
		pose.xyz( NamedAtomID( " O3'", cutpos_ ) ),
		pose.xyz( NamedAtomID( "OVL1", cutpos_ ) ) ) );
	dt_ang[ 7 ] =  phi_C4prime_C3prime_O3prime_nextP;
	Real const phi_C3prime_O3prime_nextP_nextO5prime = degrees( dihedral_radians(
		pose.xyz( NamedAtomID( " C3'", cutpos_ ) ),
		pose.xyz( NamedAtomID( " O3'", cutpos_ ) ),
		pose.xyz( NamedAtomID( "OVL1", cutpos_ ) ),
		pose.xyz( NamedAtomID( "OVL2", cutpos_ ) ) ) );
	dt_ang[ 8 ] =  phi_C3prime_O3prime_nextP_nextO5prime;
	Real const phi_O3prime_nextP_nextO5prime_nextC5prime = degrees( dihedral_radians( pose.xyz( NamedAtomID( "OVU1", cutpos_ + 1 ) ),
		pose.xyz( NamedAtomID( " P  ", cutpos_ + 1 ) ),
		pose.xyz( NamedAtomID( " O5'", cutpos_ + 1 ) ),
		pose.xyz( NamedAtomID( " C5'", cutpos_ + 1 ) ) ) );
	dt_ang[ 9 ] =  phi_O3prime_nextP_nextO5prime_nextC5prime;

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
	bridgeObjects( atoms, dt_ang, db_ang, db_len, pivots, order, t_ang_, b_ang_, b_len_, nsol_ );

	if ( nsol_ == 0 ) return false;

	figure_out_dof_ids_and_offsets ( pose, dt_ang );
	apply_solutions ( pose );
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::figure_out_dof_ids_and_offsets ( pose::Pose const & pose,
	utility::vector1< Real > const & dt_ang ) {
	////////////////////////////////////////////////////////////////////////////////////
	// Note that the torsion angles that we solved for do not directly correspond to
	// torsion angles in the atom-tree. But they are right up to an *offset*, which
	// we pre-calculate now. Also this is a good time to figure out exactly which
	// DOFs need to be changed in the atom-tree.
	//
	// Just be totally explicit.
	//
	////////////////////////////////////////////////////////////////////////////////////
	offset_save_.clear();
	dof_ids_.clear();
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
	figure_out_offset( pose, dof_id, dt_ang[ 3 * 1 + 1 ], offset_save_ );
	id1 = AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().o5prime_atom_index(), moving_suite_ + 1 );
	id2 = AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().c5prime_atom_index(), moving_suite_ + 1 );
	id3 = AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().c4prime_atom_index(), moving_suite_ + 1 );
	id4 = AtomID( pose.residue_type( moving_suite_ + 1 ).RNA_info().c3prime_atom_index(), moving_suite_ + 1 );
	dof_id = pose.atom_tree().torsion_angle_dof_id( id1, id2, id3, id4 );
	figure_out_offset( pose, dof_id, dt_ang[ 3 * 1 + 2 ], offset_save_ );
	/////////////////////////////////////////
	// pivot 2
	/////////////////////////////////////////
	id1 = AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().c4prime_atom_index(), chainbreak_suite_ );
	id2 = AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().c3prime_atom_index(), chainbreak_suite_ );
	id3 = AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().o3prime_atom_index(), chainbreak_suite_ );
	id4 = AtomID( pose.residue( chainbreak_suite_ ).atom_index( "OVL1" ),
		chainbreak_suite_ );
	dof_id = pose.atom_tree().torsion_angle_dof_id ( id1, id2, id3, id4 );
	figure_out_offset ( pose, dof_id, dt_ang[ 3 * 2 + 1 ], offset_save_ );
	id1 = AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().c3prime_atom_index(), chainbreak_suite_ );
	id2 = AtomID( pose.residue_type( chainbreak_suite_ ).RNA_info().o3prime_atom_index(), chainbreak_suite_ );
	id3 = AtomID ( pose.residue( chainbreak_suite_ ).atom_index( "OVL1" ),
		chainbreak_suite_ );
	id4 = AtomID( pose.residue( chainbreak_suite_ ).atom_index( "OVL2" ),
		chainbreak_suite_ );
	dof_id = pose.atom_tree().torsion_angle_dof_id ( id1, id2, id3, id4 );
	figure_out_offset ( pose, dof_id, dt_ang[ 3 * 2 + 2 ], offset_save_ );
	/////////////////////////////////////////
	// pivot 3
	/////////////////////////////////////////
	id1 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().p_atom_index(), chainbreak_suite_ + 1 );
	id2 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().o5prime_atom_index(), chainbreak_suite_ + 1 );
	id3 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c5prime_atom_index(), chainbreak_suite_ + 1 );
	id4 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c4prime_atom_index(), chainbreak_suite_ + 1 );
	dof_id = pose.atom_tree().torsion_angle_dof_id ( id1, id2, id3, id4 );
	figure_out_offset ( pose, dof_id, dt_ang[ 3 * 3 + 1 ], offset_save_ );
	id1 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().o5prime_atom_index(), chainbreak_suite_ + 1 );
	id2 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c5prime_atom_index(), chainbreak_suite_ + 1 );
	id3 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c4prime_atom_index(), chainbreak_suite_ + 1 );
	id4 = AtomID( pose.residue_type( chainbreak_suite_ + 1 ).RNA_info().c3prime_atom_index(), chainbreak_suite_ + 1 );
	dof_id = pose.atom_tree().torsion_angle_dof_id ( id1, id2, id3, id4 );
	figure_out_offset ( pose, dof_id, dt_ang[ 3 * 3 + 2 ], offset_save_ );
}

////////////////////////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::figure_out_offset (
	core::pose::Pose const & pose,
	core::id::DOF_ID const & dof_id,
	core::Real const & original_torsion_value,
	utility::vector1< core::Real > & offset_save ) {
	if ( dof_id == DOF_ID::BOGUS_DOF_ID() ) { //expected at cutpoint!
		TR <<  "Problem with DOF_ID " << dof_id << std::endl;
		utility_exit_with_message ( "Problem with DOF_ID" );
		//  }
	} else {
		offset_save.push_back ( pose.dof ( dof_id ) - radians ( original_torsion_value ) );
		dof_ids_.push_back ( dof_id );

		if ( verbose_ ) {
			TR << dof_id;
			TR << "  offset " << pose.dof ( dof_id ) << " " << radians ( original_torsion_value )
				<< " " << pose.dof ( dof_id ) - radians ( original_torsion_value ) << std::endl;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::apply_solutions ( core::pose::Pose & pose ) {
	debug_assert ( t_ang_.size() == Size ( nsol_ ) );

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// Finally, ready to check out the solutions
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	if ( nsol_ == 0 ) return;

	if ( choose_least_perturb_solution_ ) {
		if ( verbose_ )  {
			TR << "---------------------------------- " << std::endl;
			TR << "   start pose " << std::endl;
			TR << "---------------------------------- " << std::endl;
			utility::vector1< Real > dt_ang, db_len, db_ang;
			utility::vector1< utility::fixedsizearray1< Real,3 > > atoms;
			fill_chainTORS ( pose, atom_ids_, atoms, dt_ang, db_ang, db_len );
			output_chainTORS ( dt_ang, db_ang, db_len );
			pose.dump_pdb ( "before_closed.pdb" );
		}

		Real best_deviation2 ( 0.0 );
		Size best_sol ( 0 );
		// could save time by just looking over a subset of residues. But I don't think this is rate limiting
		utility::vector1< Vector > ref_vectors;
		Size const ref_atom ( 1 );

		for ( Size i = 1; i <= pose.size(); i++ ) {
			ref_vectors.push_back( pose.xyz( id::AtomID( ref_atom, i ) ) );
		}

		for ( Size n = 1; n <= Size( nsol_ ); n++ ) {
			fill_solution( pose, n );
			Real deviation2 ( 0.0 );

			for ( Size i = 1; i <= pose.size(); i++ ) {
				deviation2 += ( pose.xyz( id::AtomID( ref_atom, i ) ) - ref_vectors[i] ).length_squared();
			}

			if ( n == 1 || deviation2 < best_deviation2 ) {
				best_deviation2 = deviation2;
				best_sol = n;
			}
		}

		fill_solution ( pose, best_sol );

		if ( verbose_ ) { //consistency check.
			TR << "---------------------------------- " << std::endl;
			TR << "   solution " << best_sol << std::endl;
			TR << "---------------------------------- " << std::endl;
			output_chainTORS ( t_ang_[best_sol], b_ang_[best_sol], b_len_[best_sol] );
		}

		//  fill_solution( pose, best_sol );

		if ( verbose_ )  {
			TR << "pose " << best_sol << ": " << std::endl;
			utility::vector1< Real > dt_ang, db_len, db_ang;
			utility::vector1< utility::fixedsizearray1< Real,3 > > atoms;
			fill_chainTORS ( pose, atom_ids_, atoms, dt_ang, db_ang, db_len );
			output_chainTORS ( dt_ang, db_ang, db_len );
			pose.dump_pdb ( "closed.pdb" );
		}
	} else if ( choose_best_solution_ ) {
		debug_assert ( scorefxn_ != 0 );
		Real best_score ( 0.0 );
		Size best_sol ( 0 );

		for ( Size n = 1; n <= Size( nsol_ ); n++ ) {
			fill_solution ( pose, n );
			Real const score = ( *scorefxn_ )( pose );

			if ( score < best_score || n == 1 ) {
				best_score = score;
				best_sol = n;
			}

			if ( verbose_ && n == 2 ) { //consistency check.
				TR << "solution " << n << ": " << std::endl;
				output_chainTORS ( t_ang_[n], b_ang_[n], b_len_[n] );
				TR << "pose " << n << ": " << std::endl;
				utility::vector1< Real > dt_ang, db_len, db_ang;
				utility::vector1< utility::fixedsizearray1< Real,3 > > atoms;
				fill_chainTORS( pose, atom_ids_, atoms, dt_ang, db_ang, db_len );
				output_chainTORS( dt_ang, db_ang, db_len );
			}
		}

		fill_solution ( pose, best_sol );
	} else {
		debug_assert ( choose_random_solution_ );
		Size const n = static_cast< int > ( nsol_ * numeric::random::rg().uniform() ) + 1;
		fill_solution ( pose, n );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::get_all_solutions ( core::pose::Pose & pose,
	utility::vector1< core::pose::PoseOP > & pose_list ) {
	pose_list.clear();

	for ( Size n = 1; n <= Size( nsol_ ); n++ ) {
		fill_solution( pose, n );
		core::pose::PoseOP pose_save( new Pose );
		*pose_save = pose;
		pose_list.push_back ( pose_save );

		if ( verbose_ ) {
			pose.dump_pdb( "KIC_" + ObjexxFCL::string_of ( n ) + ".pdb" );
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::fill_solution ( core::pose::Pose & pose,
	Size const n ) const {
	Size count ( 0 );

	for ( Size i = 1; i <= 3; i++ ) {
		count++;
		pose.set_dof( dof_ids_[count], principal_angle ( radians ( t_ang_[ n ][ 3 * i + 1 ] ) + offset_save_[count] ) );
		count++;
		pose.set_dof( dof_ids_[count], principal_angle ( radians ( t_ang_[ n ][ 3 * i + 2 ] ) + offset_save_[count] ) );
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Real >
RNA_AnalyticLoopCloser::get_torsions ( Size const n ) {
	debug_assert ( n <= t_ang_.size() );
	utility::vector1< Real > torsions;
	Size count ( 0 );

	for ( Size i = 1; i <= 3; i++ ) {
		count++;
		torsions.push_back( degrees( principal_angle( radians( t_ang_[ n ][ 3 * i + 1 ] ) + offset_save_[count] ) ) );
		count++;
		torsions.push_back( degrees( principal_angle( radians( t_ang_[ n ][ 3 * i + 2 ] ) + offset_save_[count] ) ) );
	}

	return torsions;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< utility::vector1< Real > >
RNA_AnalyticLoopCloser::get_torsions_for_all_solutions() {
	utility::vector1< utility::vector1< Real > > torsions_for_all_solutions;

	for ( Size n = 1; n <= t_ang_.size(); n++ ) torsions_for_all_solutions.push_back ( get_torsions ( n ) );

	return torsions_for_all_solutions;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::choose_best_solution_based_on_score_function ( core::scoring::ScoreFunctionOP scorefxn ) {
	choose_best_solution_ = true;
	choose_least_perturb_solution_ = false;
	scorefxn_ = scorefxn;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::choose_least_perturb_solution() {
	choose_best_solution_ = false;
	choose_least_perturb_solution_ = true;
}

///////////////////////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::output_chainTORS ( utility::vector1< core::Real > const & dt_ang,
	utility::vector1< core::Real > const & db_ang,
	utility::vector1< core::Real > const & db_len ) const {
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

///////////////////////////////////////////////////////////
void
RNA_AnalyticLoopCloser::fill_chainTORS (
	core::pose::Pose const & pose,
	utility::vector1< id::NamedAtomID > const & atom_ids_,
	utility::vector1< utility::fixedsizearray1< Real, 3 > > & atoms,
	utility::vector1< Real > & dt_ang,
	utility::vector1< Real > & db_ang,
	utility::vector1< Real > & db_len ) const {
	using namespace core::kinematics;
	using namespace numeric::kinematic_closure;
	utility::fixedsizearray1< utility::fixedsizearray1< Real,3 >,3 > Q0;
	utility::fixedsizearray1< Real,3 > R0;
	utility::vector1< Vector > atoms_xyz;

	for ( Size i = 1; i <= atom_ids_.size(); i++ ) {
		//  TR << "filling: " << atom_ids_[i].atomno() << " " << atom_ids_[i].rsd() << std::endl;
		atoms_xyz.push_back( pose.xyz ( atom_ids_[ i ] ) );
	}

	atoms.clear();

	for ( Size i = 1; i <= atom_ids_.size(); i++ ) {
		utility::fixedsizearray1< Real,3 > atom_xyz_vals;
		atom_xyz_vals[1] = atoms_xyz[i].x();
		atom_xyz_vals[2] = atoms_xyz[i].y();
		atom_xyz_vals[3] = atoms_xyz[i].z();
		atoms.push_back ( atom_xyz_vals );
	}

	chainTORS( atoms.size(), atoms, dt_ang, db_ang, db_len, R0, Q0 );
}


} //rna
} //modeler
} //legacy
} //stepwise
} //protocols

