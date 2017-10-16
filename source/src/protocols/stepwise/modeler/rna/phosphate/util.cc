// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/rna/phosphate/PhosphateUtil.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/rna/phosphate/util.hh>
#include <protocols/stepwise/modeler/rna/phosphate/PhosphateMove.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/TorsionID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.rna.phosphate.util" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace phosphate {

/////////////////////////////////////////////////////////////////
void
remove_terminal_phosphates( pose::Pose & pose ){
	utility::vector1< Size > res_list;
	for ( Size n = 1; n <= pose.size(); n++ ) res_list.push_back( n );
	remove_terminal_phosphates( pose, res_list );
}

/////////////////////////////////////////////////////////////////
void
remove_terminal_phosphates( pose::Pose & pose, utility::vector1< Size > const & res_list ){
	for ( Size const n : res_list ) {
		if ( pose.residue(n).has_variant_type( core::chemical::FIVE_PRIME_PHOSPHATE ) ) {
			remove_variant_type_from_pose_residue( pose, core::chemical::FIVE_PRIME_PHOSPHATE, n );
			add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE, n );
		}
		remove_variant_type_from_pose_residue( pose, core::chemical::THREE_PRIME_PHOSPHATE, n );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
position_five_prime_phosphate_SLOW( pose::Pose & pose, Size const res ) {
	using namespace core::chemical;
	ResidueTypeSetCOP rsd_set = pose.residue_type_set_for_pose( pose.residue_type( res ).mode() );
	conformation::ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *( rsd_set->get_representative_type_aa( aa_from_name( "RAD") ) ) ) ;
	pose.prepend_polymer_residue_before_seqpos( *new_rsd, res, true );
	pose.delete_polymer_residue( res );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
position_five_prime_phosphate( pose::Pose & pose, Size const res ) {
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::id;

	// reposition XO3' "manually".
	Residue const & rsd = pose.residue( res );
	Vector const & P_xyz = rsd.xyz( " P  " );
	Vector const & O5prime_xyz = rsd.xyz( " O5'" );
	Vector const & C5prime_xyz = rsd.xyz( " C5'" );
	Vector const & OP2_xyz     = rsd.xyz( " OP2" );
	Vector const & XO3prime_xyz  = rsd.xyz( "XO3'" );

	Real const OP2_dihedral      = dihedral_degrees( OP2_xyz,      P_xyz, O5prime_xyz, C5prime_xyz );
	Real const XO3prime_dihedral = dihedral_degrees( XO3prime_xyz, P_xyz, O5prime_xyz, C5prime_xyz );
	static Real const desired_XO3prime_OP2_offset( 114.6 ); // from rna_phenix files.
	Real const rotation_amount = ( desired_XO3prime_OP2_offset  - ( XO3prime_dihedral - OP2_dihedral ) );

	AtomID XO3prime_ID = AtomID( rsd.atom_index( "XO3'"), res );
	DOF_ID XO3prime_DOF_ID(XO3prime_ID, PHI );

	pose.set_dof( XO3prime_DOF_ID, pose.dof( XO3prime_DOF_ID ) + numeric::conversions::radians( rotation_amount ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
copy_over_phosphate_variants( pose::Pose & pose_input,
	pose::Pose const & reference_pose,
	utility::vector1< PhosphateMove > const & phosphate_move_list )
{
	using namespace core::id;
	using namespace core::chemical;
	using namespace core::chemical::rna;
	using namespace core::pose;

	Pose pose = pose_input; // to prevent some problems with graphics thread

	runtime_assert( pose.size() == reference_pose.size() );

	for ( PhosphateMove const & phosphate_move : phosphate_move_list ) {
		Size const & n = phosphate_move.rsd();
		PhosphateTerminus const & terminus = phosphate_move.terminus();

		utility::vector1< TorsionID > torsion_ids;
		if ( terminus == FIVE_PRIME_PHOSPHATE ) {
			if ( n == 1 || pose.fold_tree().is_cutpoint( n-1 ) ) { // may be copying to a pose that no longer has terminus
				bool added_new_phosphate( false );
				if ( pose.residue( n ).has_variant_type( core::chemical::FIVE_PRIME_PHOSPHATE ) ) {
					make_variants_match( pose, reference_pose, n, core::chemical::FIVE_PRIME_PHOSPHATE );
					make_variants_match( pose, reference_pose, n, core::chemical::VIRTUAL_PHOSPHATE );
					make_variants_match( pose, reference_pose, n, core::chemical::VIRTUAL_RIBOSE );
				} else { // order matters, since there's not variant with both virtual phosphate and five prime phosphate
					make_variants_match( pose, reference_pose, n, core::chemical::LOWER_TERMINUS_VARIANT );
					make_variants_match( pose, reference_pose, n, core::chemical::VIRTUAL_RIBOSE );
					make_variants_match( pose, reference_pose, n, core::chemical::VIRTUAL_PHOSPHATE );
					make_variants_match( pose, reference_pose, n, core::chemical::FIVE_PRIME_PHOSPHATE );
					added_new_phosphate = pose.residue( n ).has_variant_type( core::chemical::FIVE_PRIME_PHOSPHATE );
				}
				if ( added_new_phosphate ) {
					position_five_prime_phosphate( pose, n);
				}
				torsion_ids.push_back( TorsionID( n, id::BB, ALPHA ) );
				torsion_ids.push_back( TorsionID( n, id::BB, BETA ) );
				torsion_ids.push_back( TorsionID( n, id::BB, GAMMA ) );
			}
		} else {
			runtime_assert( terminus == THREE_PRIME_PHOSPHATE );
			if ( n == pose.size() || pose.fold_tree().is_cutpoint( n ) ) {
				make_variants_match( pose, reference_pose, n, core::chemical::UPPER_TERMINUS_VARIANT );
				make_variants_match( pose, reference_pose, n, core::chemical::VIRTUAL_RIBOSE );
				make_variants_match( pose, reference_pose, n, core::chemical::THREE_PRIME_PHOSPHATE );
				torsion_ids.push_back( TorsionID( n, id::BB, EPSILON ) );
				torsion_ids.push_back( TorsionID( n, id::BB, ZETA ) );
			}
		}

		for ( auto const & tid : torsion_ids ) {
			pose.set_torsion( tid, reference_pose.torsion( tid ) );
		}

		make_variants_match( pose, reference_pose, n, core::chemical::VIRTUAL_RIBOSE );

	}

	core::scoring::constraints::map_constraints_from_original_pose( pose_input, pose );
	pose_input = pose; // to prevent some problems with graphics thread
}

///////////////////////////////////////////////////////////////////////////////
void
copy_over_phosphate_variants( pose::Pose & pose_input,
	pose::Pose const & reference_pose,
	utility::vector1< Size > const & res_list ){

	utility::vector1< PhosphateMove > phosphate_move_list;
	for ( Size i = 1; i <= res_list.size(); i++ ) {
		phosphate_move_list.push_back( PhosphateMove( res_list[i], FIVE_PRIME_PHOSPHATE ) );
		phosphate_move_list.push_back( PhosphateMove( res_list[i], THREE_PRIME_PHOSPHATE ) );
	}
	copy_over_phosphate_variants( pose_input, reference_pose, phosphate_move_list );
}

////////////////////////////////////////////////////////////////////
core::scoring::ScoreFunctionCOP
get_phosphate_scorefxn(){
	core::scoring::methods::EnergyMethodOptions options;
	return get_phosphate_scorefxn( options );
}

////////////////////////////////////////////////////////////////////
core::scoring::ScoreFunctionCOP
get_phosphate_scorefxn( core::scoring::methods::EnergyMethodOptions const & options ) {
	using namespace scoring;
	ScoreFunctionOP phosphate_scorefxn;
	phosphate_scorefxn = ScoreFunctionOP( new ScoreFunction ); // creating anew seems wasteful, but having problems with setting from input scorefxn.
	phosphate_scorefxn->set_energy_method_options( options );
	phosphate_scorefxn->set_weight( scoring::fa_atr, 0.1); // argh seems to trigger etable?
	phosphate_scorefxn->set_weight( scoring::fa_rep, 0.1);
	phosphate_scorefxn->set_weight( scoring::free_suite, 1.0 );
	phosphate_scorefxn->set_weight( scoring::rna_torsion, 0.5 ); // trying to increase accepts
	phosphate_scorefxn->set_weight( scoring::hbond_lr_bb_sc, 3.0 ); // trying to increase accepts
	phosphate_scorefxn->set_weight( scoring::hbond_sr_bb_sc, 3.0 ); // trying to increase accepts
	return phosphate_scorefxn;
}

//////////////////////////////////////////////////////////////////////
// Most of the following functions are helpers for PhosphateMover but are also
// useful in, e.g. deciding whether or not to virtualize phosphates in
// virtualize_free_rna_moieties()
bool
check_phosphate_contacts_donor( utility::vector1< Vector > const & op_xyz_list,
	utility::vector1< Vector > const & donor_atom_xyz_list,
	utility::vector1< Vector > const & donor_base_atom_xyz_list )
{
	static Real const max_hbond_dist_cutoff2( 2.2 * 2.2 );
	static Real const max_hbond_base_dist_cutoff2( 3.2 * 3.2 );
	static Real const cos_angle_cutoff = -0.5;

	for ( Size i = 1; i <= donor_atom_xyz_list.size(); i++ ) {
		for ( auto const & op_xyz : op_xyz_list ) {
			Real const dist2 =  ( donor_atom_xyz_list[ i ] - op_xyz ).length_squared();
			if ( dist2 > max_hbond_dist_cutoff2 ) continue;

			Real const dist_to_base2 =  ( donor_base_atom_xyz_list[ i ] - op_xyz ).length_squared();
			if ( dist_to_base2 > max_hbond_base_dist_cutoff2  ) continue;

			Real const costheta =  cos_of( donor_base_atom_xyz_list[i], donor_atom_xyz_list[i], op_xyz  );
			if ( costheta > cos_angle_cutoff ) continue;
			//    TR << "dist! " << std::sqrt( dist2 ) << " to " << i << " out of " << donor_atom_xyz_list.size() << " and cos angle is " << costheta << std::endl;
			return true;
		}
	}

	return false;
}

//////////////////////////////////////////////////////////////////////
bool
check_phosphate_contacts_donor( pose::Pose const & pose, Size const n ){
	utility::vector1< Vector > op_xyz_list;
	op_xyz_list.push_back( pose.residue( n ).xyz( " OP1" ) );
	op_xyz_list.push_back( pose.residue( n ).xyz( " OP2" ) );

	utility::vector1< Vector > donor_atom_xyz_list;
	utility::vector1< Vector > donor_base_atom_xyz_list;
	utility::vector1< Size > neighbor_copy_dofs; // dummy
	get_phosphate_atom_and_neighbor_list( pose, PhosphateMove( n, FIVE_PRIME_PHOSPHATE ), donor_atom_xyz_list, donor_base_atom_xyz_list, neighbor_copy_dofs );

	return check_phosphate_contacts_donor( op_xyz_list,
		donor_atom_xyz_list,
		donor_base_atom_xyz_list );
}

//////////////////////////////////////////////////////////////////////
void
get_phosphate_atom_and_neighbor_list( core::pose::Pose const & pose,
	PhosphateMove const & phosphate_move_,
	utility::vector1< Vector > & donor_atom_xyz_list,
	utility::vector1< Vector > & donor_base_atom_xyz_list,
	utility::vector1< Size > & neighbor_copy_dofs ) {

	static Real const  phosphate_takeoff_donor_distance_cutoff2( 8.0 * 8.0 );
	static Real const  phosphate_nbr_distance_cutoff2( 12.0 * 12.0 );

	donor_atom_xyz_list.clear();
	donor_base_atom_xyz_list.clear();
	neighbor_copy_dofs.clear();

	Size const & n = phosphate_move_.rsd();
	Vector const phosphate_takeoff_xyz = (phosphate_move_.terminus() == FIVE_PRIME_PHOSPHATE) ?
		pose.residue( n ).xyz( " C5'" ) :
		pose.residue( n ).xyz( " O3'" );

	for ( Size m = 1; m <= pose.size(); m++ ) {

		if ( ( pose.residue( m ).nbr_atom_xyz()  - phosphate_takeoff_xyz ).length_squared() > phosphate_nbr_distance_cutoff2 ) continue;
		neighbor_copy_dofs.push_back( m );

		core::chemical::AtomIndices Hpos_polar = pose.residue_type( m ).Hpos_polar();
		for ( Size const q : Hpos_polar ) {
			Vector const & donor_atom_xyz = pose.residue( m ).xyz( q );
			if ( ( donor_atom_xyz - phosphate_takeoff_xyz ).length_squared() < phosphate_takeoff_donor_distance_cutoff2 ) {
				donor_atom_xyz_list.push_back( donor_atom_xyz );
				Size  q_base = pose.residue_type(m).atom_base( q );
				Vector const & donor_base_atom_xyz = pose.residue( m ).xyz( q_base );
				donor_base_atom_xyz_list.push_back( donor_base_atom_xyz );
			}
		}
	}

}

//////////////////////////////////////////////////////////////////////
utility::vector1< bool >
detect_phosphate_contacts( pose::Pose const & pose ){
	utility::vector1< bool > phosphate_makes_contact( pose.size(), false );
	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( !pose.residue_type( i ).is_RNA() ) continue;
		phosphate_makes_contact[ i ] = check_phosphate_contacts_donor( pose, i );
	}
	return phosphate_makes_contact;
}


//////////////////////////////////////////////////////////////////////
void
setup_three_prime_phosphate_based_on_next_residue( pose::Pose & pose, Size const n ) {
	add_variant_type_to_pose_residue( pose, chemical::THREE_PRIME_PHOSPHATE, n );
	runtime_assert( n < pose.size() );
	runtime_assert( pose.residue_type( n+1 ).is_RNA() );
	runtime_assert( pose.residue_type( n+1 ).has_variant_type( chemical::VIRTUAL_PHOSPHATE ) ||
		pose.residue_type( n+1 ).has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ) );
	pose.set_xyz( id::NamedAtomID( "YP  ", n ), pose.residue( n+1 ).xyz( " P  " ) );
	pose.set_xyz( id::NamedAtomID( "YOP1", n ), pose.residue( n+1 ).xyz( " OP1" ) );
	pose.set_xyz( id::NamedAtomID( "YOP2", n ), pose.residue( n+1 ).xyz( " OP2" ) );
	pose.set_xyz( id::NamedAtomID( "YO5'", n ), pose.residue( n+1 ).xyz( " O5'" ) );
}

} //phosphate
} //rna
} //modeler
} //stepwise
} //protocols
