// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/sugar/StepWiseRNA_VirtualSugarUtil.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/rna/sugar/util.hh>
#include <protocols/stepwise/modeler/rna/sugar/StepWiseRNA_VirtualSugarJustInTimeInstantiator.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/rna/sugar/SugarModeling.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeSet.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/copydofs/util.hh>
#include <core/pose/util.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/stream_util.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.rna.sugar.util" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace sugar {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is called *after* VirtualSugarSampler in two pipelines. Why not just run it
//  during the actual modeler? Well, there is a combinatorial mix/match
//  that occurs in VirtualSugarJustInTimeInstantiator -- up to four sugars might be instantiated
//  and they are applied in all available combinations and then those residues are
//  reminimized.
//
void
minimize_all_sampled_floating_bases( core::pose::Pose & viewer_pose,
																		 utility::vector1< SugarModeling > const & modeling_list,
																		 utility::vector1< PoseOP > & pose_data_list,
																		 core::scoring::ScoreFunctionCOP const & scorefxn,
																		 working_parameters::StepWiseWorkingParametersCOP const & working_parameters_,
																		 bool const virtual_sugar_is_from_prior_step_ ){

	using namespace core::optimization;
	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::io::silent;
	//	using namespace protocols::farna;
	using namespace core::id;
	using namespace ObjexxFCL;

	output_title_text( "Enter minimize_all_sampled_floating_bases", TR.Debug );

	pose::Pose const viewer_pose_copy = viewer_pose;

	if ( modeling_list.size() == 0 ) return;
	if ( pose_data_list.size() == 0 ) return;

	ScoreFunctionOP modeler_scorefxn = get_modeler_scorefxn( scorefxn );

	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.00000025 );
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true );

	core::kinematics::MoveMap mm;
	mm.set_bb( false );
	mm.set_chi( false );

	for ( Size n = 1; n <= modeling_list.size(); n++ ){

		SugarModeling const & sugar_modeling_ = modeling_list[n];

		// actually should have this on, but need to run some legacy code
		//		runtime_assert( pose_data_list[1]->residue( sugar_modeling_.bulge_res ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) );

		mm.set( TorsionID( sugar_modeling_.bulge_res - 1, id::BB,  5 ), true ); //epsilon
		mm.set( TorsionID( sugar_modeling_.bulge_res - 1, id::BB,  6 ), true ); //zeta

		mm.set( TorsionID( sugar_modeling_.bulge_res, id::BB,  1 ), true ); //alpha
		mm.set( TorsionID( sugar_modeling_.bulge_res, id::BB,  2 ), true ); //beta
		mm.set( TorsionID( sugar_modeling_.bulge_res, id::BB,  3 ), true ); //gamma
		mm.set( TorsionID( sugar_modeling_.bulge_res, id::BB,  5 ), true ); //epsilon
		mm.set( TorsionID( sugar_modeling_.bulge_res, id::BB,  6 ), true ); //zeta
		mm.set( TorsionID( sugar_modeling_.bulge_res, id::CHI, 1 ), true ); //chi (torsion between base and sugar sugar)

		mm.set( TorsionID( sugar_modeling_.bulge_res + 1, id::BB,  1 ), true ); //alpha
		mm.set( TorsionID( sugar_modeling_.bulge_res + 1, id::BB,  2 ), true ); //beta
		mm.set( TorsionID( sugar_modeling_.bulge_res + 1, id::BB,  3 ), true ); //gamma

	}

	for ( Size n = 1; n <= pose_data_list.size(); n++ ){

		TR << TR.Blue << "Minimizing pose " << n <<  " out of " << pose_data_list.size() << TR.Reset << std::endl;

		viewer_pose = ( *pose_data_list[n] );

		utility::vector1 < core::Size > const & working_moving_partition_res = working_parameters_->working_moving_partition_res();
		utility::vector1 < core::Size > already_virtualized_res_list;

		if ( virtual_sugar_is_from_prior_step_ ){ //Virtualize the other partition since it doesn't exist in prior step!
			for ( Size ii = 1; ii <= working_moving_partition_res.size(); ii++ ){
				Size const seq_num = working_moving_partition_res[ii];
				if ( viewer_pose.residue( seq_num ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ){
					already_virtualized_res_list.push_back( seq_num );
					continue;
				}
				pose::add_variant_type_to_pose_residue( viewer_pose, core::chemical::VIRTUAL_RNA_RESIDUE, seq_num );
			}
			if ( working_parameters_->gap_size() == 0 ) {
				pose::add_variant_type_to_pose_residue( viewer_pose,
						core::chemical::VIRTUAL_PHOSPHATE, working_parameters_->five_prime_chain_break_res() + 1 );
			}
		}

		TR.Debug << "minimize_all_sampled_floating_bases pose # " << n << " out of " << pose_data_list.size() << " " << std::endl;;
		minimizer.run( viewer_pose, mm, ( *modeler_scorefxn ), options );

		if ( virtual_sugar_is_from_prior_step_ ){ //Virtualize the other partition since it doesn't exist in prior step!
			for ( Size ii = 1; ii <= working_moving_partition_res.size(); ii++ ){
				Size const seq_num = working_moving_partition_res[ii];
				if ( already_virtualized_res_list.has_value( seq_num ) ) continue;
				pose::remove_variant_type_from_pose_residue( viewer_pose, core::chemical::VIRTUAL_RNA_RESIDUE, seq_num );
			}
			if ( working_parameters_->gap_size() == 0 ) {
				pose::remove_variant_type_from_pose_residue( viewer_pose,
						core::chemical::VIRTUAL_PHOSPHATE, working_parameters_->five_prime_chain_break_res() + 1 );
			}
		}

		( *pose_data_list[n] ) = viewer_pose;
	}

	viewer_pose = viewer_pose_copy;

	output_title_text( "Exit minimize_all_sampled_floating_bases", TR.Debug );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
is_sugar_virtual( core::pose::Pose const & pose, core::Size const sugar_res, core::Size const bulge_res) {
	utility::vector1< Size > bulge_residues_to_virtualize_dummy;
	return is_sugar_virtual( pose, sugar_res, bulge_res, bulge_residues_to_virtualize_dummy );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
is_sugar_virtual( core::pose::Pose const & pose, core::Size const sugar_res, core::Size const bulge_res,
									utility::vector1< Size > & bulge_residues_to_virtualize ){

	using namespace core::chemical;
	using namespace core::conformation;
	using namespace ObjexxFCL;

	Size const nres = pose.total_residue();

	if ( ( sugar_res + 1 ) != bulge_res && ( sugar_res - 1 ) != bulge_res ) {
		TR.Debug << "sugar_res = " << sugar_res << " bulge_res = " << bulge_res << std::endl;
		utility_exit_with_message( "( sugar_res + 1 ) != bulge_res && ( sugar_res - 1 ) != bulge_res )" );
	}

	if ( sugar_res < 1 || sugar_res > nres ){
		utility_exit_with_message( "sugar_res < 1 || sugar_res > nres( " + string_of( nres ) + " )!. sugar_res = " + string_of( sugar_res ) );
	}

	if ( pose.residue( sugar_res ).has_variant_type( core::chemical::VIRTUAL_RIBOSE ) ) {

		//bulge consistency checks:
		runtime_assert ( bulge_res > 1 && bulge_res < nres );
		// will virtualize any bulge residues...
		if ( !pose.residue( bulge_res ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) &&
				!bulge_residues_to_virtualize.has( bulge_res ) ) bulge_residues_to_virtualize.push_back( bulge_res );

		return true;
	} else{
		return false;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
modeler_starting_pose_data_list( utility::vector1< PoseOP > & starting_pose_data_list,
																	 utility::vector1< SugarModeling > const & sugar_modeling_list,
																	 core::pose::Pose const & pose ){

	pose::Pose const pose_copy = pose;

	utility::vector1< Size > sugar_ID_counter_list( sugar_modeling_list.size(), 1 );

	while ( true ){

		PoseOP start_pose = pose_copy.clone();

		for ( Size n = 1; n <= sugar_modeling_list.size(); n++ ){
			SugarModeling const & curr_modeling = sugar_modeling_list[n];
			Size const sugar_ID = sugar_ID_counter_list[n];
			tag_into_pose( *start_pose, tag_from_pose( *start_pose ) + tag_from_pose( *curr_modeling.pose_list[sugar_ID] ) );
			pose::remove_variant_type_from_pose_residue( *start_pose,
					core::chemical::VIRTUAL_RIBOSE, curr_modeling.moving_res );
			copy_bulge_res_and_sugar_torsion( curr_modeling, *start_pose, ( *curr_modeling.pose_list[sugar_ID] ) );
		}

		starting_pose_data_list.push_back( start_pose );

		///////////////////////////Counter/////////////////////////////
		sugar_ID_counter_list[1]++;

		for ( Size n = 1; n < sugar_modeling_list.size(); n++ ){
			if ( sugar_ID_counter_list[n] == ( sugar_modeling_list[n].pose_list.size() + 1 ) ){
				 sugar_ID_counter_list[n] = 1;
				 sugar_ID_counter_list[n + 1]++;
			}
		}

		if ( sugar_ID_counter_list[sugar_modeling_list.size()] == ( sugar_modeling_list[sugar_modeling_list.size()].pose_list.size() + 1 ) ) break;
		////////////////////////////////////////////////////////////////

	}

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::map< core::Size, core::Size >
get_res_map( SugarModeling const & sugar_modeling ){

	std::map< core::Size, core::Size > 	res_map;
	res_map[ sugar_modeling.moving_res    ] = sugar_modeling.moving_res;
	if ( sugar_modeling.bulge_res > 0 ) {
		res_map[ sugar_modeling.bulge_res     ] = sugar_modeling.bulge_res;
		res_map[ sugar_modeling.reference_res ] = sugar_modeling.reference_res;
	}
	return res_map;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_bulge_res_and_sugar_torsion( SugarModeling const & sugar_modeling, core::pose::Pose & pose, core::pose::Pose const & template_pose,
																	bool instantiate_sugar /* = false */){

	using namespace ObjexxFCL;
	using namespace core::io::silent;
	using namespace core::id;

	if ( instantiate_sugar ){
		runtime_assert( pose.residue( sugar_modeling.moving_res ).has_variant_type( core::chemical::VIRTUAL_RIBOSE ) );
		pose::remove_variant_type_from_pose_residue( pose, core::chemical::VIRTUAL_RIBOSE, sugar_modeling.moving_res );
}

	//This is map from sub numbering to input_res numbering..
	std::map< core::Size, core::Size > res_map = get_res_map( sugar_modeling );

	//Dec 24, 2011 Parin S.:Convert to Rhiju's NEW version
	core::pose::copydofs::copy_dofs_match_atom_names( pose, template_pose, res_map, false /*backbone_only*/, false /*side_chain_only*/, false /*ignore_virtual*/ );

}


	////////////////////////////////////////////////////////////////
	// Following is used in pose_setup to define jumps...
	//
 	// Virtual sugars (VIRTUAL_RIBOSE variant type) should only occur in SWA if the residue was created via a
	//  'floating base' move, with a jump to an anchor (reference) residue in the rest of the pose.
	//
	// For now, we are  also assuming that each residue can only act as an anchor for one other residue. This is the case for
	//  SWA of RNA motifs, but may no longer hold when we allow multiple strands and/or ligands.
	//
	// In the more general case, the pattern of anchors should still be a directed acyclic graph, and, I think, will
	// be parseable if the pose already has jump information. Following is written though as if the pose doesn't have
	// jumps already defined -- a separate set of functions that take advantage of that jump information are below.
	//
	//     -- rhiju, dec. 2013
	//
 	///////////////////////////////////////////////////////////////////////////////////////
	std::map< Size, Size > const
	get_reference_res_for_each_virtual_sugar( pose::Pose const & pose,
																						bool const check_for_non_jump /* = false */,
																						Size const moving_suite /*cannot place jump across partititions*/	){

		if ( check_for_non_jump ) {
			return get_reference_res_for_each_virtual_sugar_without_fold_tree( pose, moving_suite );
		}
		return get_reference_res_for_each_virtual_sugar_based_on_fold_tree( pose );
	}

 	///////////////////////////////////////////////////////////////////////////////////////
	std::map< Size, Size > const
	get_reference_res_for_each_virtual_sugar_without_fold_tree( pose::Pose const & pose, Size const moving_suite /*cannot place jump across partititions*/	){

		std::map< Size, Size > reference_res_for_each_virtual_sugar;
		utility::vector1< Size > virtual_sugar_res;
		utility::vector1< utility::vector1< Size > > possible_reference_res_lists;

		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( !pose.residue( n ).has_variant_type( core::chemical::VIRTUAL_RIBOSE ) ) continue;
			virtual_sugar_res.push_back( n );
			reference_res_for_each_virtual_sugar[ n ] = 0;
			utility::vector1< Size > possible_reference_res = get_possible_reference_res_list_from_pose_without_fold_tree( n, pose,
																																																										 moving_suite /* cannot place jump across partitions */ );
			runtime_assert( possible_reference_res.size() == 1 || possible_reference_res.size() == 2 );
			possible_reference_res_lists.push_back( possible_reference_res );
		}

		Size const num_virtual_sugar_res = virtual_sugar_res.size();
		// Look for terminal virtual sugars first, where the anchor res should be totally obvious,
		// and then work our way to an internal virtual sugars.
		for ( Size num_pass = 1; num_pass <= num_virtual_sugar_res; num_pass++ ){

			Size found_reference_res( 0 );
			Size virt_sugar_res( 0 );
			for ( Size i = 1; i <= num_virtual_sugar_res; i++ ){
				virt_sugar_res = virtual_sugar_res[i];
				if ( reference_res_for_each_virtual_sugar[ virt_sugar_res ] == 0 &&
						 possible_reference_res_lists[ i ] .size() == 1 ){
					found_reference_res = possible_reference_res_lists[ i ][ 1 ];
					reference_res_for_each_virtual_sugar[ virt_sugar_res ] = found_reference_res;
					TR.Debug << "Found anchor res " << found_reference_res << " for sugar " << virt_sugar_res << std::endl;
					break;
				}
			}
			runtime_assert( found_reference_res > 0 );

			// the sugar cannot be an anchor for any other res -- an *assumption*.
			// eliminate from the list of possibilities
			for ( Size i = 1; i <= num_virtual_sugar_res; i++ ){
			 	utility::vector1< Size > const & possible_reference_res_list = possible_reference_res_lists[ i ];
			 	utility::vector1< Size > new_list;
			 	for ( Size k = 1; k <= possible_reference_res_list.size(); k++ ){
			 		if ( /*possible_reference_res_list[k] != found_reference_res &&*/
							possible_reference_res_list[k] != virt_sugar_res )  new_list.push_back( possible_reference_res_list[k] );
			 	}
			 	possible_reference_res_lists[i] = new_list;
			}

		}

		for ( Size i = 1; i <= num_virtual_sugar_res; i++ )	runtime_assert ( reference_res_for_each_virtual_sugar[ virtual_sugar_res[i] ] > 0 );
		return reference_res_for_each_virtual_sugar;
	}

	////////////////////////////////////////////////////////////////
	Size
	look_for_non_jump_reference_to_previous( Size const virtual_sugar_res,
																					 pose::Pose const & pose,
																					 Size const moving_suite ){
		// Following was to figure out attachment location for old-school Parin fold tree.
		Size i = virtual_sugar_res - 1;
		while ( i >= 1 ){
			if ( pose.fold_tree().is_cutpoint(i) ) break;
			if ( !pose.residue( i ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ){
				if ( i != moving_suite  ){
					return i;
				}
			}
			i--;
		}
		return 0;
	}

	////////////////////////////////////////////////////////////////
	Size
	look_for_non_jump_reference_to_next( Size const virtual_sugar_res,
																					 pose::Pose const & pose,
																					 Size const moving_suite ){
		// Following was to figure out attachment location for old-school Parin fold tree.
		Size i = virtual_sugar_res + 1;
		while ( i <= pose.total_residue() ){
			if ( pose.fold_tree().is_cutpoint( i - 1 ) ) break;
			if ( !pose.residue( i ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ){
				if ( i != (moving_suite - 1) ){
					return i;
				}
			}
			i++;
		}
		return 0;
	}


	////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	get_possible_reference_res_list_from_pose_without_fold_tree( Size const virtual_sugar_res,
																															 pose::Pose const & pose,
																															 Size const  moving_suite /* for old-school poses without jumps already setup*/){
		utility::vector1< Size > possible_reference_res_list;
		Size possible_reference_res( 0 );
		possible_reference_res = look_for_non_jump_reference_to_previous( virtual_sugar_res, pose, moving_suite );
		if ( possible_reference_res == 0 ) possible_reference_res = look_for_jumps_to_previous( virtual_sugar_res, pose, true /*force_upstream*/ );
		if ( possible_reference_res > 0 ) possible_reference_res_list.push_back( possible_reference_res );

		possible_reference_res = look_for_non_jump_reference_to_next( virtual_sugar_res, pose, moving_suite );
		if ( possible_reference_res == 0 ) possible_reference_res = look_for_jumps_to_next( virtual_sugar_res, pose, true /*force_upstream*/ );
		if ( possible_reference_res > 0 ) possible_reference_res_list.push_back( possible_reference_res );

		TR.Debug << pose.fold_tree() << std::endl;
		TR.Debug << "ASSUMING MOVING_SUITE " << moving_suite << std::endl;
		TR.Debug << "REFERENCE_RES_LIST FOR " << virtual_sugar_res << " is " << possible_reference_res_list << std::endl;

		return possible_reference_res_list;
	}

 	///////////////////////////////////////////////////////////////////////////////////////
	std::map< Size, Size > const
	get_reference_res_for_each_virtual_sugar_based_on_fold_tree( pose::Pose const & pose ){

		std::map< Size, Size > reference_res_for_each_virtual_sugar;
		for ( Size n = 1; n <= pose.total_residue(); n++ ){

			if ( !pose.residue( n ).has_variant_type( core::chemical::VIRTUAL_RIBOSE ) ) continue;
			Size possible_reference_res = look_for_jumps( n, pose, true /*force_upstream*/ );

			if ( !possible_reference_res ){
				// did not find an upstream residue. as a last resort look for a downstream jump residue.
				// that can happen if fold tree gets rerooted upon a merge.
				possible_reference_res = look_for_jumps( n, pose, false /*force_upstream*/ );
			}
			runtime_assert( possible_reference_res );
			reference_res_for_each_virtual_sugar[ n ] = possible_reference_res;
		}
		return reference_res_for_each_virtual_sugar;
	}


 	///////////////////////////////////////////////////////////////////////////////////////
	Size
	look_for_jumps( Size const n, pose::Pose const & pose, bool const force_upstream ){
		Size const possible_reference_res_previous = look_for_jumps_to_previous( n, pose, force_upstream );
		Size const possible_reference_res_next = look_for_jumps_to_next( n, pose, force_upstream );
		runtime_assert( !possible_reference_res_previous || !possible_reference_res_next );
		if ( possible_reference_res_previous )	return  possible_reference_res_previous;
		if ( possible_reference_res_next ) return  possible_reference_res_next;
		return 0;
	}


	////////////////////////////////////////////////////////////////
	Size
	look_for_jumps_to_previous( Size const virtual_sugar_res,
															pose::Pose const & pose,
															bool const force_upstream ){
		Size i = virtual_sugar_res - 1;
		while ( i >= 1 ) { // look for jumps with reference residue 'upstream'
			Size const jump_nr = pose.fold_tree().jump_nr( i, virtual_sugar_res );
			if ( jump_nr > 0 && (!force_upstream || Size( pose.fold_tree().upstream_jump_residue( jump_nr ) ) == i ) ) {
				return i;
			}
			i--;
		}
		return 0;
	}


	////////////////////////////////////////////////////////////////
	Size
	look_for_jumps_to_next( Size const virtual_sugar_res,
															pose::Pose const & pose,
															bool const force_upstream ){
		Size i = virtual_sugar_res + 1;
		while ( i <= pose.total_residue() ) {  // look for jumps with reference residue 'upstream'
			Size const jump_nr = pose.fold_tree().jump_nr( i, virtual_sugar_res );
			if ( jump_nr > 0 && (!force_upstream || Size( pose.fold_tree().upstream_jump_residue( jump_nr ) ) == i ) ) {
				return i;
			}
			i++;
		}
		return 0;
	}


	sampler::copy_dofs::ResidueAlternativeSetOP
	convert_sugar_modeling_to_residue_alternative_set( SugarModeling const & sugar_modeling ){
		return sampler::copy_dofs::ResidueAlternativeSetOP( new sampler::copy_dofs::ResidueAlternativeSet( sugar_modeling.pose_list,
																																	get_res_map( sugar_modeling ),
																																	sugar_modeling.moving_res ) );
	}


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	StepWiseRNA_VirtualSugarJustInTimeInstantiatorOP
	instantiate_any_virtual_sugars( pose::Pose & pose,
																	working_parameters::StepWiseWorkingParametersCOP working_parameters,
																	core::scoring::ScoreFunctionCOP scorefxn,
																	options::StepWiseModelerOptionsCOP options ) {
		Pose pose_save = pose;
		StepWiseRNA_VirtualSugarJustInTimeInstantiatorOP virtual_sugar_just_in_time_instantiator( new StepWiseRNA_VirtualSugarJustInTimeInstantiator( working_parameters ) );
		virtual_sugar_just_in_time_instantiator->set_scorefxn( scorefxn );
		virtual_sugar_just_in_time_instantiator->set_options( options );
		virtual_sugar_just_in_time_instantiator->apply( pose );
		pose = pose_save;
		return virtual_sugar_just_in_time_instantiator;
	}

	///////////////////////////////////////////////////////////////////////////
	utility::vector1< bool >
	detect_sugar_contacts( pose::Pose const & pose ) {
		utility::vector1< bool > sugar_makes_contact( pose.total_residue(), false );
		for ( Size i = 1; i <= pose.total_residue(); i++ ) {
			if ( !pose.residue_type( i ).is_RNA() ) continue;
			sugar_makes_contact[ i ] = detect_sugar_contacts( pose, i );
		}
		return sugar_makes_contact;
	}


	///////////////////////////////////////////////////////////////////////////
	bool
	detect_sugar_contacts( pose::Pose const & pose, Size const moving_res,
												 Distance o2prime_contact_distance_cutoff_ ) {

		bool found_contact( false );
		Vector const & moving_O2prime_xyz = pose.residue( moving_res ).xyz( " O2'" );
		core::pose::PDBInfoCOP pdb_info = pose.pdb_info();
		for ( Size i = 1; i <= pose.total_residue(); i++ ){
			if ( i == moving_res ) continue;

			// sometimes see artifactual H-bonds of O2' to next phosphate or sugar.
			if ( i > 1 &&
					 pdb_info->number( i )-1 == pdb_info->number( moving_res ) &&
					 pdb_info->chain( i ) == pdb_info->chain( moving_res ) ) continue;

			if ( pose.residue( i ).is_virtual_residue() ) continue;
			for ( Size j = 1; j <= pose.residue( i ).nheavyatoms(); j++ ){
				if ( pose.residue_type( i ).is_virtual( j ) ) continue;
				if ( pose.residue_type( i ).heavyatom_is_an_acceptor( j ) ||
						 pose.residue_type( i ).heavyatom_has_polar_hydrogens( j )  ){
					Distance dist = ( pose.residue(i).xyz( j ) - moving_O2prime_xyz ).length();
					//std::cout << "Distance to rsd " << i << "  atom " << pose.residue_type( i ).atom_name( j ) << ": " << dist << std::endl;
					if ( dist < o2prime_contact_distance_cutoff_ ) {
						found_contact = true; break;
					}
				}
			} // atoms j
			if ( found_contact ) break;
		} // residues i

		return found_contact;
	}

} //sugar
} //rna
} //modeler
} //stepwise
} //protocols
