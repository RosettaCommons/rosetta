// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/phosphate/PhosphateUtil.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/rna/phosphate/util.hh>
#include <protocols/stepwise/modeler/rna/phosphate/PhosphateMove.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.rna.phosphate.util" );

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
		for ( Size n = 1; n <= pose.total_residue(); n++ )res_list.push_back( n );
		remove_terminal_phosphates( pose, res_list );
	}

	/////////////////////////////////////////////////////////////////
	void
	remove_terminal_phosphates( pose::Pose & pose, utility::vector1< Size > const & res_list ){
		for ( Size i = 1; i <= res_list.size(); i++ ){
			Size const & n = res_list[ i ];
			if ( pose.residue(n).has_variant_type( core::chemical::FIVE_PRIME_PHOSPHATE ) ){
				remove_variant_type_from_pose_residue( pose, core::chemical::FIVE_PRIME_PHOSPHATE, n );
				add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE, n );
			}
			remove_variant_type_from_pose_residue( pose, core::chemical::THREE_PRIME_PHOSPHATE, n );
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	correctly_position_five_prime_phosphate_SLOW( pose::Pose & pose, Size const res ) {
		using namespace core::chemical;
		ResidueTypeSet const & rsd_set = pose.residue( res ).residue_type_set();
		conformation::ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *( rsd_set.aa_map( aa_from_name( "RAD") )[1] ) ) ;
		pose.prepend_polymer_residue_before_seqpos( *new_rsd, res, true );
		pose.delete_polymer_residue( res );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	correctly_position_five_prime_phosphate( pose::Pose & pose, Size const res ) {
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

		runtime_assert( pose.total_residue() == reference_pose.total_residue() );

		for ( Size i = 1; i <= phosphate_move_list.size(); i++ ){
			PhosphateMove const & phosphate_move = phosphate_move_list[ i ];
			Size const & n = phosphate_move.rsd();
			PhosphateTerminus const & terminus = phosphate_move.terminus();

			utility::vector1< TorsionID > torsion_ids;
			if ( terminus == FIVE_PRIME_PHOSPHATE ){
				if ( n == 1 || pose.fold_tree().is_cutpoint( n-1 ) ) { // may be copying to a pose that no longer has terminus
					bool added_new_phosphate( false );
					if ( pose.residue( n ).has_variant_type( core::chemical::FIVE_PRIME_PHOSPHATE ) ){
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
					if ( added_new_phosphate ){
						correctly_position_five_prime_phosphate( pose, n);
					}
					torsion_ids.push_back( TorsionID( n, id::BB, ALPHA ) );
					torsion_ids.push_back( TorsionID( n, id::BB, BETA ) );
					torsion_ids.push_back( TorsionID( n, id::BB, GAMMA ) );
				}
			} else {
				runtime_assert( terminus == THREE_PRIME_PHOSPHATE );
				if ( n == pose.total_residue() || pose.fold_tree().is_cutpoint( n ) ) {
					make_variants_match( pose, reference_pose, n, core::chemical::UPPER_TERMINUS_VARIANT );
					make_variants_match( pose, reference_pose, n, core::chemical::VIRTUAL_RIBOSE );
					make_variants_match( pose, reference_pose, n, core::chemical::THREE_PRIME_PHOSPHATE );
					torsion_ids.push_back( TorsionID( n, id::BB, EPSILON ) );
					torsion_ids.push_back( TorsionID( n, id::BB, ZETA ) );
				}
			}

			for ( Size m = 1; m <= torsion_ids.size(); m++ ) {
				pose.set_torsion( torsion_ids[m], reference_pose.torsion( torsion_ids[m] ) );
			}

			make_variants_match( pose, reference_pose, n, core::chemical::VIRTUAL_RIBOSE );

		}

		pose_input = pose; // to prevent some problems with graphics thread
	}

	///////////////////////////////////////////////////////////////////////////////
	void
	copy_over_phosphate_variants( pose::Pose & pose_input,
																pose::Pose const & reference_pose,
																utility::vector1< Size > const & res_list ){

		utility::vector1< PhosphateMove > phosphate_move_list;
		for ( Size i = 1; i <= res_list.size(); i++ ){
			phosphate_move_list.push_back( PhosphateMove( res_list[i], FIVE_PRIME_PHOSPHATE ) );
			phosphate_move_list.push_back( PhosphateMove( res_list[i], THREE_PRIME_PHOSPHATE ) );
		}
		copy_over_phosphate_variants( pose_input, reference_pose, phosphate_move_list );
	}

	////////////////////////////////////////////////////////////////////
	core::scoring::ScoreFunctionCOP
	get_phosphate_scorefxn(){
		return get_phosphate_scorefxn( 0 );
	}

	////////////////////////////////////////////////////////////////////
	// I originally wanted this to be efficient by taking etables
	//  and stuff from the input scorefxn, but couldn't get that to work..
	core::scoring::ScoreFunctionCOP
	get_phosphate_scorefxn( core::scoring::ScoreFunctionCOP /*scorefxn*/ ){
		using namespace scoring;
		ScoreFunctionOP phosphate_scorefxn;
		//		if ( scorefxn != 0 ) {
		//			phosphate_scorefxn = scorefxn->clone();
		//			//for ( Size n = 1; n <= end_of_score_type_enumeration; n++ ) phosphate_scorefxn->set_weight( ScoreType(n), 0.0 );
		//		} else {
		phosphate_scorefxn = new ScoreFunction; // creating anew seems wasteful, but having problems with setting from input scorefxn.
			//		}
		phosphate_scorefxn->set_weight( scoring::fa_atr, 0.1); // argh seems to trigger etable?
		phosphate_scorefxn->set_weight( scoring::fa_rep, 0.1);
		phosphate_scorefxn->set_weight( scoring::free_suite, 1.0 );
		phosphate_scorefxn->set_weight( scoring::rna_torsion, 0.5 ); // trying to increase accepts
		phosphate_scorefxn->set_weight( scoring::hbond_lr_bb_sc, 3.0 ); // trying to increase accepts
		phosphate_scorefxn->set_weight( scoring::hbond_sr_bb_sc, 3.0 ); // trying to increase accepts
		return phosphate_scorefxn;
	}

} //phosphate
} //rna
} //modeler
} //stepwise
} //protocols
