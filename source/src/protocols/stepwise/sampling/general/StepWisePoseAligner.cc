// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/general/StepWisePoseAligner.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/general/StepWisePoseAligner.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.general.StepWisePoseAligner" );
using ObjexxFCL::format::F;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace general {

	//Constructor
	StepWisePoseAligner::StepWisePoseAligner( pose::Pose const & reference_pose ):
		reference_pose_( reference_pose ),
		skip_bulges_( false )
	{
	}

	//Destructor
	StepWisePoseAligner::~StepWisePoseAligner()
	{}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWisePoseAligner::apply( pose::Pose & pose ){
		superimpose_at_fixed_res( pose );
	}

	//////////////////
	void
	StepWisePoseAligner::superimpose_at_fixed_res( pose::Pose & pose ){

		using namespace core::chemical;
		using namespace core::id;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace core::scoring;

		reference_pose_local_ = reference_pose_.clone(); // local working copy, mutated in cases where nucleotides have been designed ('n')

		// first need to slice up reference_pose to match residues in actual pose.
		// define atoms over which to compute RMSD, using rmsd_res.
		FullModelInfo & full_model_info = nonconst_full_model_info( pose );
		utility::vector1< Size > const & res_list = full_model_info.res_list();
		utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();
		std::string const full_sequence = full_model_info.full_sequence();

		if ( res_list.size() == 0 ) return; // special case -- blank pose.

		utility::vector1< Size > rmsd_res_in_pose, superimpose_res_in_pose;
		get_rmsd_res_and_superimpose_res_in_pose( pose, rmsd_res_in_pose, superimpose_res_in_pose );

		utility::vector1< Size > calc_rms_res;
		utility::vector1< Size > skipped_residues;

		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( rmsd_res_in_pose.has_value( n ) ) {
				// I don't see the purpose of Arvind's skipped_residues -- need to check this.
				if ( skip_bulges_ && residue_is_bulged( pose, n ) && residue_is_bulged( *reference_pose_local_, res_list[ n ] ) ) {
					skipped_residues.push_back( n );
					continue;
				}
				calc_rms_res.push_back( n );

				char const pose_nt = pose.sequence()[ n-1 ];
				if ( full_sequence[ res_list[ n ] - 1 ] == 'n' ){
					// need to generalize to protein... should be trivial.
					pose::rna::mutate_position( *reference_pose_local_, res_list[ n ], pose_nt );
				} else {
					runtime_assert( full_sequence[ res_list[ n ] - 1 ] == pose_nt);
				}
				runtime_assert( reference_pose_local_->sequence()[ res_list[ n ] - 1] == pose_nt );
			}
		}

		// super special case.
		if ( calc_rms_res.size() == 0 && pose.total_residue() == 1 ) calc_rms_res.push_back( 1 );

		calc_rms_atom_id_map_.clear();
		for ( Size k = 1; k <= calc_rms_res.size(); k++ ){
			Size const n = calc_rms_res[ k ];
			for ( Size q = 1; q <= pose.residue_type( n ).nheavyatoms(); q++ ){
				add_to_atom_id_map_after_checks( calc_rms_atom_id_map_,
																				 pose.residue_type( n ).atom_name( q ),
																				 n, res_list[ n ],
																				 pose, *reference_pose_local_ );
			}
		}

		utility::vector1< Size > calc_rms_suites;
		// additional RNA suites over which to calculate RMSD
		for ( Size n = 1; n < pose.total_residue(); n++ ){

			if ( ( pose.residue_type( n ).is_RNA() && pose.residue_type( n + 1 ).is_RNA() ) ||
					 ( pose.residue_type( n ).is_protein() && pose.residue_type( n + 1 ).is_protein() ) ){
				// Atoms at ends of rebuilt loops:
				if ( !pose.fold_tree().is_cutpoint( n ) || pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) ) {
					if ( calc_rms_res.has_value( n ) || calc_rms_res.has_value( n+1) ){
						calc_rms_suites.push_back( n ); continue;
					}
				}
				// Domain boundaries:
				if ( (res_list[ n+1 ] == res_list[ n ] + 1) &&
						 fixed_domain_map[ res_list[ n ] ] != 0 &&
						 fixed_domain_map[ res_list[ n+1 ] ] != 0 &&
						 fixed_domain_map[ res_list[ n ] ] != fixed_domain_map[ res_list[ n+1 ] ] ){
					calc_rms_suites.push_back( n );
				}
			}
		}

		utility::vector1< std::string > const extra_suite_atoms_upper = utility::tools::make_vector1( " P  ", " OP1", " OP2", " O5'" /*upper*/ );
		utility::vector1< std::string > const extra_suite_atoms_lower = utility::tools::make_vector1( " O  " /*lower*/ );
		for ( Size k = 1; k <= calc_rms_suites.size(); k++ ){
			Size const n = calc_rms_suites[ k ];
			for ( Size q = 1; q <= extra_suite_atoms_upper.size(); q++ ){
				add_to_atom_id_map_after_checks( calc_rms_atom_id_map_, extra_suite_atoms_upper[ q ],
																				 n+1, res_list[ n+1 ],
																				 pose, *reference_pose_local_ );
			}
			for ( Size q = 1; q <= extra_suite_atoms_lower.size(); q++ ){
				add_to_atom_id_map_after_checks( calc_rms_atom_id_map_, extra_suite_atoms_lower[ q ],
																				 n, res_list[ n ],
																				 pose, *reference_pose_local_ );
			}
		}
		// output_atom_id_map( calc_rms_atom_id_map_ );

		// define superposition atoms. Should be over atoms in any fixed domains. This should be
		// the 'inverse' of calc_rms atoms.
		superimpose_atom_id_map_.clear();
		for ( Size n = 1; n < pose.total_residue(); n++ ){
			if ( skipped_residues.has_value( n ) ) continue;
			if ( !superimpose_res_in_pose.has_value( n ) ) continue;
			for ( Size q = 1; q <= pose.residue_type( n ).nheavyatoms(); q++ ){
				if ( calc_rms_atom_id_map_.find( AtomID( q, n ) ) == calc_rms_atom_id_map_.end() ){
					add_to_atom_id_map_after_checks( superimpose_atom_id_map_,
													pose.residue_type( n ).atom_name( q ),
													n, res_list[ n ],
													pose, *reference_pose_local_ );
				}
			}
		}

		// What if there weren't any fixed atoms? superimpose over everything.
		if ( superimpose_atom_id_map_.size() == 0 ) superimpose_atom_id_map_ = calc_rms_atom_id_map_;
		//		output_atom_id_map( superimpose_atom_id_map_ );

		Real superimpose_rms( 0.0 );
		rmsd_ = 0.0;
		if ( superimpose_atom_id_map_.size() > 0 ) {
			superimpose_rms = scoring::superimpose_pose( pose, *reference_pose_local_, superimpose_atom_id_map_ );
			if ( calc_rms_atom_id_map_.size() > 0 ) {
				rmsd_ = rms_at_corresponding_atoms_no_super( pose, *reference_pose_local_, calc_rms_atom_id_map_ );
			}
		}

		TR << "RMSD " << F(5,3,rmsd_) <<
			" (" << natoms_rmsd() << " atoms in " << make_tag_with_dashes( sub_to_full(rmsd_res_in_pose,pose) ) << "), superimposed on " << superimpose_atom_id_map_.size() << " atoms in " <<
			make_tag_with_dashes( sub_to_full(superimpose_res_in_pose,pose) ) << " (RMSD " <<
			F(5,3,superimpose_rms) << ") " << std::endl;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Are there any fixed domains in the pose? If so, superimpose on the first on of those, and calculate rmsd over
	// everything else.
	// If no fixed domains, calculate rmsd over everything.
	Size
	StepWisePoseAligner::get_rmsd_res_and_superimpose_res_in_pose( pose::Pose const & pose,
																																 utility::vector1< Size > & rmsd_res_in_pose,
																																 utility::vector1< Size > & superimpose_res_in_pose ) const {

		utility::vector1< Size > const domain_map = core::pose::full_model_info::get_fixed_domain_from_full_model_info_const( pose );

		// figure out 'primary' domain number. Smallest number that is not zero.
		// must be drawn from root_partition (if that partition is defined)
		Size d_primary = 0;
		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( root_partition_res_.size() > 0 && !root_partition_res_.has_value( n ) ) continue;
			Size const d = domain_map[ n ];
			if ( d > 0 && ( d_primary == 0 || d < d_primary ) ) d_primary = d;
		}

		rmsd_res_in_pose.clear();
		superimpose_res_in_pose.clear();
		if ( d_primary == 0 ){ // superimpose on everything.
			for ( Size n = 1; n <= pose.total_residue(); n++ ) {
				if ( !root_partition_res_.has_value( n ) )	rmsd_res_in_pose.push_back( n );
				if ( root_partition_res_.size() == 0 || root_partition_res_.has_value( n ) )	superimpose_res_in_pose.push_back( n );
			}
		} else { // superimpose on primary domain, calculate rmsd over rest.
			for ( Size n = 1; n <= pose.total_residue(); n++ ) {
				if ( domain_map[ n ] == d_primary ) {
					superimpose_res_in_pose.push_back( n );
				} else {
					if ( !root_partition_res_.has_value( n ) )  rmsd_res_in_pose.push_back( n );
				}
			}
		}
		return d_primary;
	}

	/////////////////////////////////////////////////////////////////////////
	// adapted from arvind kannan's original hack -- rhiju, 2014
	void
	StepWisePoseAligner::add_coordinate_constraints_from_map( pose::Pose & pose, pose::Pose const & reference_pose,
																														std::map< id::AtomID, id::AtomID > const & superimpose_atom_id_map,
																														core::Real const & constraint_x0, core::Real const & constraint_tol ) const {
		using namespace core::scoring::func;
		using namespace core::scoring::constraints;

		Size const my_anchor( pose.fold_tree().root() ); //Change to use the root of the current foldtree as done by Rocco in AtomCoordinateCstMover - JAB.

		ConstraintSetOP cst_set = pose.constraint_set()->clone();
		FuncOP constraint_func = new FlatHarmonicFunc( constraint_x0, 1.0, constraint_tol );
		//FuncOP constraint_func = new FadeFunc( -0.7, 1.5, 0.8, -1.0, 0.0)) );

		for ( std::map< id::AtomID, id::AtomID >::const_iterator
			 it=superimpose_atom_id_map.begin(), it_end = superimpose_atom_id_map.end(); it != it_end; ++it ) {
			id::AtomID const mapped_atom = it->second;
			ConstraintOP constraint = new CoordinateConstraint ( it->first, id::AtomID(1, my_anchor),
																													 reference_pose.residue(mapped_atom.rsd()).xyz(mapped_atom.atomno()),
																													 constraint_func );
			cst_set->add_constraint( constraint );
		}

		pose.constraint_set( cst_set );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWisePoseAligner::create_coordinate_constraints( pose::Pose & pose,
																											Real const rmsd_screen ){

		using namespace core::pose::full_model_info;

		pose.remove_constraints();
		if ( rmsd_screen == 0.0 ) return;
		runtime_assert( reference_pose_local_ != 0 ); // needs to be setup by apply() above.

		FullModelInfo & full_model_info = nonconst_full_model_info( pose );
		utility::vector1< Size > const & res_list = full_model_info.res_list();
		std::map< id::AtomID, id::AtomID> coordinate_constraint_atom_id_map;
		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			for ( Size q = 1; q <= pose.residue_type( n ).nheavyatoms(); q++ ){
				add_to_atom_id_map_after_checks( coordinate_constraint_atom_id_map,
																				 pose.residue_type( n ).atom_name( q ),
																				 n, res_list[ n ],
																				 pose, *reference_pose_local_ );
			}
		}

		//		output_atom_id_map( coordinate_constraint_atom_id_map );

		Real const constraint_x0  = 0.0; // stay near native.
		Real const constraint_tol = rmsd_screen; // no penalty for deviations up to this amount. After that, (x - tol)^2.
		add_coordinate_constraints_from_map( pose, *reference_pose_local_, coordinate_constraint_atom_id_map,
																				 constraint_x0, constraint_tol );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// encodes different choices between RNA and proteins, including all-heavy-atom for RNA and just backbone-atoms for proteins.
	void
	StepWisePoseAligner::add_to_atom_id_map_after_checks( std::map< id::AtomID, id::AtomID> & atom_id_map,
																												std::string const & atom_name,
																												Size const & n1, Size const & n2,
																												pose::Pose const & pose1, pose::Pose const & pose2 ){

		using namespace core::id;


		runtime_assert ( n1 >= 1 && n1 <= pose1.total_residue() );
		runtime_assert ( n2 >= 1 && n2 <= pose2.total_residue() );
		runtime_assert( pose1.residue_type( n1 ).aa() == pose2.residue_type( n2 ).aa() );

		if ( ! pose1.residue_type( n1 ).has( atom_name ) ) return;
		if ( ! pose2.residue_type( n2 ).has( atom_name ) ) return;

		Size const idx1 = pose1.residue_type( n1 ).atom_index( atom_name );
		Size const idx2 = pose2.residue_type( n2 ).atom_index( atom_name );

		if ( pose1.residue_type( n1 ).is_virtual( idx1 ) ) return;
		if ( pose2.residue_type( n2 ).is_virtual( idx2 ) ) return;

		if ( pose1.residue_type( n1 ).is_protein() && ( idx1 >= pose1.residue_type( n1 ).first_sidechain_atom() ) ) return;
		if ( pose2.residue_type( n2 ).is_protein() && ( idx2 >= pose2.residue_type( n2 ).first_sidechain_atom() ) ) return;

		atom_id_map[ AtomID( idx1, n1 ) ] = AtomID( idx2, n2 );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWisePoseAligner::residue_is_bulged( pose::Pose const & pose, Size const & resid ) {
		boost::unordered_map < core::Size , core::Size > num_stacks = pose.get_stacking_map();
		if ( num_stacks[ resid ] < 2 ) {
			return true;
		}
		return false;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWisePoseAligner::output_atom_id_map( std::map< id::AtomID, id::AtomID > const & atom_id_map ){
		for ( std::map < id::AtomID, id::AtomID >::const_iterator it = atom_id_map.begin();
					it != atom_id_map.end(); it++ ){
			TR << it->first << " mapped to " << it->second << std::endl;
		}
	}


} //general
} //sampling
} //stepwise
} //protocols
