// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/align/StepWisePoseAligner.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/align/StepWisePoseAligner.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/rna/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.align.StepWisePoseAligner" );
using ObjexxFCL::format::F;
using utility::tools::make_vector1;
using namespace core::scoring;
using namespace core::pose::full_model_info;

///////////////////////////////////////////////////////////////////////////////////////////
//
// Nicely factored object that (1) can figure out which atoms to use for superposition
//  of protein, RNA, or mixed poses, (2) which atoms to use for calculating RMSDs (based
//  on information stored in fixed_domain in pose full_model_info), and (3)
//  calculate those RMSDs with or without superposition.
//
// This is the central place to do RMSD calculations, and the object is used
//  by get_rmsd(), StepWiseClusterer, NativeRMSD_Screener, and setting of
//  coordinate constaints during StepWise Monte Carlo & Assembly. So if we want
//  to expand SWA/SWM to more stuff, like ligands, DNA, or metal ions, encode
//  the choices in RMSD calculations *here*.
//
// -- rhiju, 2014
//
///////////////////////////////////////////////////////////////////////////////////////////



namespace protocols {
namespace stepwise {
namespace sampling {
namespace align {

	//Constructor
	StepWisePoseAligner::StepWisePoseAligner( pose::Pose const & reference_pose ):
		reference_pose_( reference_pose ),
		skip_bulges_( false ),
		rmsd_( 0.0 ),
		superimpose_rmsd_( 0.0 ),
		check_alignment_tolerance_( 1.0e-3 )
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


	///////////////////////////////////////////////////////////////////////////////////
	Real
	StepWisePoseAligner::superimpose_at_fixed_res( pose::Pose & pose ){
		initialize( pose );
		return do_superimposition( pose );
	}

	///////////////////////////////////////////////////////////////////////////////////
	Real
	StepWisePoseAligner::get_rmsd_no_superimpose( pose::Pose const & pose,
																								bool const check_align /* = true */ ){
		superimpose_rmsd_ = 0.0;
		rmsd_ = 0.0;

		runtime_assert( pose.annotated_sequence() == annotated_sequence_used_for_atom_id_maps_ );

		if ( superimpose_atom_id_map_.size() > 0  && check_align ) {
			Real const check_alignment_tolerance_( 1.0e-3 );
			superimpose_rmsd_ = rms_at_corresponding_atoms_no_super( pose, *reference_pose_local_, superimpose_atom_id_map_ );
			runtime_assert( superimpose_rmsd_ <= check_alignment_tolerance_ );

		}

		if ( calc_rms_atom_id_map_.size() > 0 ) {
			rmsd_ = rms_at_corresponding_atoms_no_super( pose, *reference_pose_local_, calc_rms_atom_id_map_ );
		}

		return rmsd_;
	}


	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWisePoseAligner::initialize( pose::Pose const & pose ){
		update_reference_pose_local( pose );
		get_rmsd_res_and_superimpose_res_in_pose( pose );
		update_calc_rms_atom_id_map( pose );
		update_superimpose_atom_id_map( pose );
	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWisePoseAligner::update_reference_pose_local( pose::Pose const & pose ){

		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
		utility::vector1< Size > const res_list_in_reference = get_res_list_in_reference( pose );
		std::string const & full_sequence = const_full_model_info( pose ).full_sequence();

		// local working copy, mutated in cases where nucleotides have been designed ('n')
		if ( reference_pose_local_ == 0 ) reference_pose_local_ = reference_pose_.clone();

		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			char const pose_nt = pose.sequence()[ n-1 ];
			if ( full_sequence[ res_list[ n ] - 1 ] == 'n' ) {
				if ( reference_pose_local_->sequence()[ res_list_in_reference[n] - 1 ] != pose_nt ) {
					// need to alignize to protein... should be trivial.
					pose::rna::mutate_position( *reference_pose_local_, res_list_in_reference[n], pose_nt );
				}
			} else {
				runtime_assert( full_sequence[ res_list[ n ] - 1 ] == pose_nt );
			}
			runtime_assert( reference_pose_local_->sequence()[ res_list_in_reference[n] - 1] == pose_nt );
		}
	}


	///////////////////////////////////////////////////////////////////////////////////
	// where each residue in pose ends up in the reference_pose
	utility::vector1< Size >
	StepWisePoseAligner::get_res_list_in_reference( pose::Pose const & pose ) const {

		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
		utility::vector1< Size > const & reference_pose_res_list = get_res_list_from_full_model_info_const( reference_pose_ );

		utility::vector1< Size > res_list_in_reference;
		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			Size const reference_pose_res = reference_pose_res_list.index( res_list[ n ] );
 			res_list_in_reference.push_back( reference_pose_res );
		}
		return res_list_in_reference;
	}


	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWisePoseAligner::update_calc_rms_atom_id_map( pose::Pose const & pose ){

		using namespace core::chemical;

		utility::vector1< Size > const res_list_in_reference = get_res_list_in_reference( pose );

		// first need to slice up reference_pose to match residues in actual pose.
		// define atoms over which to compute RMSD, using rmsd_res.
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const & res_list = full_model_info.res_list();
		utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();
		std::string const full_sequence = full_model_info.full_sequence();

		if ( res_list.size() == 0 ) return; // special case -- blank pose.

		utility::vector1< Size > calc_rms_res;
		skipped_res_.clear();
		for ( Size i = 1; i <= rmsd_res_in_pose_.size(); i++ ){
			Size const & n = rmsd_res_in_pose_[ i ];
			//  introduced by Arvind. Perhaps not relevant anymore.
			if ( skip_bulges_ && residue_is_bulged( pose, n ) && residue_is_bulged( *reference_pose_local_, res_list[ n ] ) ) {
				skipped_res_.push_back( n );	continue;
			}
			if ( user_defined_calc_rms_res_.size() > 0 && !user_defined_calc_rms_res_.has_value( n ) ){
			 	skipped_res_.push_back( n ); continue;
			}
			calc_rms_res.push_back( n );
		}

		// super special case.
		if ( calc_rms_res.size() == 0 && pose.total_residue() == 1 ) calc_rms_res.push_back( 1 );

		calc_rms_atom_id_map_.clear();
		for ( Size k = 1; k <= calc_rms_res.size(); k++ ){
			Size const n = calc_rms_res[ k ];
			for ( Size q = 1; q <= pose.residue_type( n ).nheavyatoms(); q++ ){
				add_to_atom_id_map_after_checks( calc_rms_atom_id_map_,
																				 pose.residue_type( n ).atom_name( q ),
																				 n, res_list_in_reference[ n ],
																				 pose, *reference_pose_local_ );
			}
		}

		// additional RNA & protein 'suites' (connections from i to i+1) over which to calculate RMSD
		utility::vector1< Size > calc_rms_suites;
		for ( Size n = 1; n < pose.total_residue(); n++ ){

			if ( ( pose.residue_type( n ).is_RNA()     && pose.residue_type( n + 1 ).is_RNA() ) ||
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

		utility::vector1< std::string > const extra_suite_atoms_upper = make_vector1( " P  ", " OP1", " OP2", " O5'" );
		utility::vector1< std::string > const extra_suite_atoms_lower = make_vector1( " O  " );
		for ( Size k = 1; k <= calc_rms_suites.size(); k++ ){
			Size const n = calc_rms_suites[ k ];
			for ( Size q = 1; q <= extra_suite_atoms_upper.size(); q++ ){
				add_to_atom_id_map_after_checks( calc_rms_atom_id_map_, extra_suite_atoms_upper[ q ],
																				 n+1, res_list_in_reference[ n+1 ],
																				 pose, *reference_pose_local_ );
			}
			for ( Size q = 1; q <= extra_suite_atoms_lower.size(); q++ ){
				add_to_atom_id_map_after_checks( calc_rms_atom_id_map_, extra_suite_atoms_lower[ q ],
																				 n, res_list_in_reference[ n ],
																				 pose, *reference_pose_local_ );
			}
		}

		annotated_sequence_used_for_atom_id_maps_ = pose.annotated_sequence();

		// output_atom_id_map( calc_rms_atom_id_map_ );
	}


	///////////////////////////////////////////////////////////////////////////////////
	// define superposition atoms. Should be over atoms in any fixed domains. This should be
	// the 'inverse' of calc_rms atoms.
	void
	StepWisePoseAligner::update_superimpose_atom_id_map( pose::Pose const & pose ) {
		using namespace core::id;
		utility::vector1< Size > const res_list_in_reference = get_res_list_in_reference( pose );

		superimpose_atom_id_map_.clear();
		for ( Size n = 1; n < pose.total_residue(); n++ ){
			if ( skipped_res_.has_value( n ) ) continue;
			if ( !superimpose_res_in_pose_.has_value( n ) ) continue;
			for ( Size q = 1; q <= pose.residue_type( n ).nheavyatoms(); q++ ){
				if ( calc_rms_atom_id_map_.find( AtomID( q, n ) ) == calc_rms_atom_id_map_.end() ){
					add_to_atom_id_map_after_checks( superimpose_atom_id_map_,
													pose.residue_type( n ).atom_name( q ),
													n, res_list_in_reference[n],
													pose, *reference_pose_local_ );
				}
			}
		}

		// What if there weren't any fixed atoms? superimpose over everything.
		// CHANGE THIS? how about superposition over just N-CA-C triad (protein) or nucleobase (RNA), which should not change.
		if ( superimpose_atom_id_map_.size() == 0 ) {
			superimpose_atom_id_map_ = get_root_triad_atom_id_map( pose );
		}

		//		output_atom_id_map( superimpose_atom_id_map_ );
  }


	/////////////////////////////////////////////////////////////////////////////////////////////
	Real
	StepWisePoseAligner::do_superimposition( pose::Pose & pose ) {

		superimpose_rmsd_ = 0.0;
		rmsd_ = 0.0;

		if ( superimpose_atom_id_map_.size() > 0 ) {
			superimpose_rmsd_ = scoring::superimpose_pose( pose, *reference_pose_local_, superimpose_atom_id_map_ );
			if ( calc_rms_atom_id_map_.size() > 0 ) {
				rmsd_ = rms_at_corresponding_atoms_no_super( pose, *reference_pose_local_, calc_rms_atom_id_map_ );
			}
		}

		TR << "RMSD " << F(5,3,rmsd_) <<
			" (" << natoms_rmsd() << " atoms in " << make_tag_with_dashes( sub_to_full(rmsd_res_in_pose_,pose) ) << "), superimposed on " << superimpose_atom_id_map_.size() << " atoms in " <<
			make_tag_with_dashes( sub_to_full(superimpose_res_in_pose_,pose) ) << " (RMSD " <<
			F(9,7,superimpose_rmsd_) << ") " << std::endl;

		return rmsd_;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Are there any fixed domains in the pose? If so, superimpose on the first on of those, and calculate rmsd over
	// everything else.
	// If no fixed domains, calculate rmsd over everything.
	Size
	StepWisePoseAligner::get_rmsd_res_and_superimpose_res_in_pose( pose::Pose const & pose ) {

		utility::vector1< Size > domain_map = get_fixed_domain_from_full_model_info_const( pose );
		utility::vector1< Size > const & extra_minimize_res = const_full_model_info( pose ).extra_minimize_res();
		utility::vector1< Size > const & res_list = const_full_model_info( pose ).res_list();

		if ( user_defined_calc_rms_res_.size() > 0 ){
			for ( Size n = 1; n <= user_defined_calc_rms_res_.size(); n++ )  domain_map[ user_defined_calc_rms_res_[n] ] = 0;
		}
		//		TR << "DOMAIN MAP " << domain_map << std::endl;

		// figure out 'primary' domain number. Smallest number that is not zero.
		// must be drawn from root_partition (if that partition is defined)
		Size d_primary = 0;
		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( root_partition_res_.size() > 0 && !root_partition_res_.has_value( n ) ) continue;
			if ( extra_minimize_res.has_value( res_list[n] ) ) continue;
			Size const d = domain_map[ n ];
			if ( d > 0 && ( d_primary == 0 || d < d_primary ) ) d_primary = d;
		}
		//		TR << "PRIMARY DOMAIN " << d_primary << std::endl;

		rmsd_res_in_pose_.clear();
		superimpose_res_in_pose_.clear();
		if ( d_primary == 0 ){ // superimpose on everything.
			for ( Size n = 1; n <= pose.total_residue(); n++ ) {
				if ( !root_partition_res_.has_value( n ) )	rmsd_res_in_pose_.push_back( n );
				if ( root_partition_res_.size() == 0 || root_partition_res_.has_value( n ) )	superimpose_res_in_pose_.push_back( n );
			}
		} else { // superimpose on primary domain, calculate rmsd over rest.
			for ( Size n = 1; n <= pose.total_residue(); n++ ) {
				if ( domain_map[ n ] == d_primary ) {
					if ( extra_minimize_res.has_value( res_list[n] ) ) continue;
					if ( root_partition_res_.size() == 0 || root_partition_res_.has_value( n ) )	superimpose_res_in_pose_.push_back( n );
				} else {
					if ( !root_partition_res_.has_value( n ) )  rmsd_res_in_pose_.push_back( n );
				}
			}
		}
		return d_primary;
	}

	/////////////////////////////////////////////////////////////////////////
	// adapted from arvind kannan's original hack -- rhiju, 2014
	void
	StepWisePoseAligner::add_coordinate_constraints_from_map( pose::Pose & pose, pose::Pose const & reference_pose,
																														std::map< id::AtomID, id::AtomID > const & atom_id_map,
																														core::Real const & constraint_x0, core::Real const & constraint_tol ) const {
		using namespace core::scoring::func;
		using namespace core::scoring::constraints;

		Size const my_anchor( pose.fold_tree().root() ); //Change to use the root of the current foldtree as done by Rocco in AtomCoordinateCstMover - JAB.

		ConstraintSetOP cst_set = pose.constraint_set()->clone();
		FuncOP constraint_func = new FlatHarmonicFunc( constraint_x0, 1.0, constraint_tol );
		//FuncOP constraint_func = new FadeFunc( -0.7, 1.5, 0.8, -1.0, 0.0)) );

		for ( std::map< id::AtomID, id::AtomID >::const_iterator
			 it=atom_id_map.begin(), it_end = atom_id_map.end(); it != it_end; ++it ) {
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
	// encodes different choices between RNA and proteins:
	//     all-heavy-atom for RNA -- no terminal phosphates.
	//     just backbone-atoms for proteins.
	void
	StepWisePoseAligner::add_to_atom_id_map_after_checks( std::map< id::AtomID, id::AtomID> & atom_id_map,
																												std::string const & atom_name,
																												Size const & n1, Size const & n2,
																												pose::Pose const & pose1, pose::Pose const & pose2 ) const {

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

		// no terminal phosphates...
		if ( pose1.residue_type( n1 ).is_RNA() && !pose1.residue_type( n1 ).has_variant_type( "CUTPOINT_UPPER" )
				 && ( n1 == 1 || pose1.fold_tree().is_cutpoint( n1 - 1 ) ) ) return;
		if ( pose2.residue_type( n2 ).is_RNA() && !pose2.residue_type( n2 ).has_variant_type( "CUTPOINT_UPPER" )
				 && ( n2 == 1 || pose2.fold_tree().is_cutpoint( n2 - 1 ) ) ) return;

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
	StepWisePoseAligner::output_atom_id_map( std::map< id::AtomID, id::AtomID > const & atom_id_map ) const {
		for ( std::map < id::AtomID, id::AtomID >::const_iterator it = atom_id_map.begin();
					it != atom_id_map.end(); it++ ){
			TR << it->first << " mapped to " << it->second << std::endl;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWisePoseAligner::output_atom_id_map( std::map< id::AtomID, id::AtomID > const & atom_id_map,
																					 pose::Pose const & pose1,
																					 pose::Pose const & pose2 ) const {
		for ( std::map < id::AtomID, id::AtomID >::const_iterator it = atom_id_map.begin();
					it != atom_id_map.end(); it++ ){
			TR << it->first << " " << pose1.residue( it->first.rsd() ).atom_name( it->first.atomno() ) <<
				" mapped to " <<
				it->second << " " << pose2.residue( it->second.rsd() ).atom_name( it->second.atomno() ) << std::endl;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::map< id::AtomID, id::AtomID >
	StepWisePoseAligner::get_root_triad_atom_id_map( pose::Pose const & pose ) const {
		Size const root_res = pose.fold_tree().root();
		Size const root_atomno = get_root_residue_root_atomno( pose.residue( root_res ), pose.fold_tree() );
		core::kinematics::tree::AtomCOP root_atom ( & pose.atom_tree().atom_dont_do_update( id::AtomID( root_atomno, root_res ) ) );
		utility::vector1< core::kinematics::tree::AtomCOP > stub_atoms =
			make_vector1( root_atom->stub_atom1(),
										root_atom->stub_atom2(),
										root_atom->stub_atom3() );

		std::map< id::AtomID, id::AtomID > root_triad_atom_id_map;
		//utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
		utility::vector1< Size > const res_list_in_reference = get_res_list_in_reference( pose );
		for ( Size i = 1; i <= 3; i++ ){
			Size const n = stub_atoms[i]->id().rsd();
			Size const q = stub_atoms[i]->id().atomno();
			add_to_atom_id_map_after_checks( root_triad_atom_id_map,
																			 pose.residue_type( n ).atom_name( q ),
																			 n, res_list_in_reference[ n ],
																			 pose, *reference_pose_local_ );
		}
		return root_triad_atom_id_map;

	}


} //align
} //sampling
} //stepwise
} //protocols
