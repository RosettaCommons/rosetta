// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler/util.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/recces/sampler/util.hh>
#include <protocols/recces/sampler/MC_Any.hh>
#include <protocols/recces/sampler/MC_Loop.hh>
#include <protocols/recces/sampler/MC_OneTorsion.hh>
#include <protocols/recces/sampler/rna/MC_RNA_Suite.hh>
#include <protocols/recces/sampler/rna/MC_RNA_MultiSuite.hh>
#include <protocols/recces/sampler/rna/MC_RNA_KIC_Sampler.hh>
#include <protocols/recces/sampler/rna/MC_RNA_OneJump.hh>
#include <protocols/recces/params/RECCES_Parameters.hh>
#include <protocols/recces/options/RECCES_Options.hh>
#include <protocols/recces/util.hh>
#include <protocols/farna/secstruct/RNA_SecStructLegacyInfo.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.recces.sampler.util" );

using namespace core;
using namespace core::pose;
using namespace protocols::recces::sampler::rna;

namespace protocols {
namespace recces {
namespace sampler {

//////////////////////////////////////////////////////////////////////////////
///  TODO -- *unify* recces_turner & thermal_sampler setup by making
///           HELIX and DANGLE full_model_parameters specify recces_turner-like
///           residues. May also want rna_secstruct setup in full_model_parameters!
///
///          Tests: rb_recces, recces_turner, & thermal_sampler integration tests:
///                    must not change.
///
///          More tests: 1. handle long helices setup with -seq & -jump_res in the middle.
///                      2. long helices + KIC to improve sampling.
///                      3. check compatibility of sample_jump with helix+dangling ends.
///                      4. multiple dangles hanging off helix.
///
protocols::recces::sampler::MC_CombOP
initialize_sampler( pose::Pose const & pose,
										options::RECCES_Options const & options,
										params::RECCES_Parameters const & params )
{
	using namespace protocols::recces::sampler;
	using namespace protocols::recces::sampler::rna;

	MC_CombOP sampler( new MC_Comb );

	// TODO: Instead, make this *the* 'standard_bb' sampler used inside thermal_sampler-style sampler.
	//       do it in such a way that all tests remain unchanged, i.e. chi_sampler will be unfilled.
	if ( options.legacy_turner_mode() ) {
		sampler = get_recces_turner_sampler( pose, options.a_form_range(), options.rna_secstruct(), params );
	}

	// TODO: deprecate this? Its a legacy of rb_recces, and could be handled separately.
	// make sure base pair sampling can be carried out after legacy_turner_mode. Again, handle
	// as one of the samplers inside the thermal_sampler-style MC_Loop sampler.
	if (  !options.thermal_sampler_mode() &&
				( ( pose.size() == 2 && pose.fold_tree().is_cutpoint( 1 ) ) || options.sample_jump() ) )  {
		MC_RNA_OneJumpOP jump_sampler = initialize_jump_sampler( pose, 1, options );
		sampler->add_rotamer( jump_sampler );
		print_base_centroid_atoms_for_rb_entropy( pose.residue( pose.fold_tree().downstream_jump_residue( 1 ) ), options.xyz_file() );
	}

	// sweet: thermal sampler mode instead. this should probably be *default* style.
	if ( sampler->num_rotamers() == 0 ) {
		runtime_assert( options.thermal_sampler_mode() );
		runtime_assert( options.sample_residues().size() > 0 ); // the standard signature for thermal_sampler
		sampler = initialize_thermal_sampler( pose, options );
	}

	sampler->init();

	if ( !options.suppress_sampler_display() ) sampler->show( TR, 0 );

	return sampler;
}


////////////////////////////////////////////////////////////////////////
///  @details pretty satisfactory sampler setup function. Developed
///       to handle single-nucleotide bulge & other motifs more
///       complex than helices.
///
///  TODO: This remains a little cryptic due to use of MC_RNA_Suite wrapper
///         which has convenient set_sample_near_a_form functions ...
///         would be better (and pretty easy) to create more 'atomic' MC_RNA_Sugar,
///         MC_RNA_Chi, and MC_RNA_SuiteBB functions with set_sample_near_a_form()
///         function. -- rhiju, jan2017
///
protocols::recces::sampler::MC_CombOP
get_recces_turner_sampler_from_secstruct( pose::Pose const & pose,
																					core::Real const & a_form_range,
																					core::pose::rna::RNA_SecStruct const & rna_secstruct )
{
	MC_CombOP sampler( new MC_Comb );
	Size n_rsd( pose.size() );
	/// Chi & sugar samplers. Use MC_RNA_Suite due to its useful 'aform_range' helper functions.
	for ( Size i = 1; i <= n_rsd; ++ i ) {
		if ( pose.residue( i ).has_variant_type( chemical::VIRTUAL_RIBOSE ) ) continue;
		MC_RNA_SuiteOP suite_sampler( new MC_RNA_Suite( i ) );
		suite_sampler->set_sample_bb( false );
		suite_sampler->set_sample_lower_nucleoside( true ); // sugar and chi of i
		suite_sampler->set_sample_upper_nucleoside( false );
		suite_sampler->set_sample_near_a_form( rna_secstruct.in_helix( i ) );
		suite_sampler->set_a_form_range( a_form_range );
		sampler->add_rotamer( suite_sampler );
	}

	/// Suite samplers. Use MC_RNA_Suite due to its useful 'aform_range' helper functions.
	for ( Size i = 1; i < n_rsd; ++ i ) {
		if ( pose.fold_tree().is_cutpoint( i ) ) continue;
		MC_RNA_SuiteOP suite_sampler( new MC_RNA_Suite( i ) );
		suite_sampler->set_sample_bb( true ); // 5 backbone torsions between i and i+1.
		suite_sampler->set_sample_lower_nucleoside( false );
		suite_sampler->set_sample_upper_nucleoside( false );
		suite_sampler->set_sample_near_a_form( rna_secstruct.in_same_helix( i, i+1 ) );
		suite_sampler->set_a_form_range( a_form_range );
		sampler->add_rotamer( suite_sampler );
	}
	return sampler;
}

////////////////////////////////////////////////////////////////////////
/// @detailed Legacy sampler setup from Fang
///     + goes through a lot of cryptic book-keeping in order to fit
///        samplers into "MultiSuite" framework.
///
///
/// TODO: Could unify with get_recces_turner_sampler_from_secstruct() if
///        we had a silly function that could set up order of chi & bb samplers
///        to *exactly* match the order created by this legacy order -- and
///        then check *exact* numerical trajectory with recces_turner integration test.
///       Use of recces_turner app would trigger that reordering (or perhaps recognition
///        that sequence is a perfect helix, perhaps with dangle in recces_turner style).
///       Note that we should also be able to deprecate bp_res and dangle_res --
///        and actually then deprecate RECCES_Parameters (since those are its only variables).
///                              -- rhiju, 2016
protocols::recces::sampler::MC_CombOP
get_recces_turner_sampler_legacy( pose::Pose const & pose,
													 core::Real const & a_form_range,
													 params::RECCES_Parameters const & params )
{
	MC_RNA_MultiSuiteOP sampler( new MC_RNA_MultiSuite );
	Size total_len( pose.size() );
	runtime_assert( total_len == params.bp_res().size() + params.dangling_res().size() );
	Size len1( 1 );
	while ( !pose.fold_tree().is_cutpoint( len1 ) && len1 < total_len ) len1++; // length of first strand.
	for ( Size i = 1; i <= total_len; ++i ) {
		bool const sample_near_a_form( params.bp_res().has_value( i ) );
		if ( i == 1 || ( i > len1 && i != total_len ) ) {
			MC_RNA_SuiteOP suite_sampler( new MC_RNA_Suite( i ) );
			suite_sampler->set_sample_bb( i != 1 );
			suite_sampler->set_sample_lower_nucleoside( true );
			suite_sampler->set_sample_upper_nucleoside( false );
			suite_sampler->set_sample_near_a_form( sample_near_a_form );
			suite_sampler->set_a_form_range( a_form_range );
			sampler->add_rotamer( suite_sampler );
		} else {
			MC_RNA_SuiteOP suite_sampler( new MC_RNA_Suite( i - 1 ) );
			suite_sampler->set_sample_bb( len1 == total_len || i != total_len );
			suite_sampler->set_sample_lower_nucleoside( false );
			suite_sampler->set_sample_upper_nucleoside( true );
			suite_sampler->set_sample_near_a_form( sample_near_a_form );
			suite_sampler->set_a_form_range( a_form_range );
			sampler->add_rotamer( suite_sampler );
		}
	}
	return sampler;
}

////////////////////////////////////////////////////////////////////////
protocols::recces::sampler::MC_CombOP
get_recces_turner_sampler( pose::Pose const & pose,
													 core::Real const & a_form_range,
													 core::pose::rna::RNA_SecStruct const & secstruct,
													 params::RECCES_Parameters const & params )
{
	if ( secstruct.blank() ) {
		// this one is legacy
		return get_recces_turner_sampler_legacy( pose, a_form_range, params);
	} else {
		// user-defined secstruct.
		return get_recces_turner_sampler_from_secstruct( pose, a_form_range, secstruct );
	}
}

////////////////////////////////////////////////////////////////////////
MC_CombOP
initialize_thermal_sampler( pose::Pose const & pose,
	options::RECCES_Options const & options )
{
	using namespace sampler;
	using namespace sampler::rna;
	using namespace core::id;

	utility::vector1< Size > const & sample_res( options.sample_residues() );
	utility::vector1< Size > const & free_res( options.free_residues() );
	std::cout << "Sample residues: " << sample_res << std::endl;

	PoseCOP pose_cop( pose.get_self_ptr() ); // oh please mercy mercy

	//Set up the internal move samplers
	MC_AnyOP internal_bb_sampler( new MC_Any );
	for ( Size i = 1; i<= sample_res.size(); ++i ) {
		if ( sample_res[ i ] == 1 ) continue;
		if ( pose.fold_tree().is_cutpoint( sample_res[ i ] - 1 ) ) continue;
		if ( pose.fold_tree().is_cutpoint( sample_res[ i ] ) ) continue; // already a cutpoint there
		MC_RNA_KIC_SamplerOP suite_sampler( new MC_RNA_KIC_Sampler( pose_cop, sample_res[i]-1, sample_res[i]) );
		suite_sampler->init();
		if ( free_res.has_value( sample_res[i] ) || ( i > 1 && free_res.has_value( sample_res[i-1] ) ) ) {
			suite_sampler->set_angle_range_from_init_torsions( options.angle_range_free_bb() );
		} else {
			suite_sampler->set_angle_range_from_init_torsions( options.angle_range_bb() );
		}
		internal_bb_sampler->add_rotamer( suite_sampler );
	}
	internal_bb_sampler->set_name( "Backbone (KIC)" );
	internal_bb_sampler->set_update_pose( pose_cop );
	internal_bb_sampler->init();

	//Set up the chi samplers
	MC_AnyOP chi_sampler( new MC_Any );
	core::Real init_torsion;
	utility::vector1<TorsionID> chi_torsion_ids;
	for ( Size i = 1; i<= sample_res.size(); ++i ) {
		if ( pose.residue( sample_res[i] ).has_variant_type( chemical::VIRTUAL_RIBOSE ) ) continue;
		TorsionID chi_ID (TorsionID( sample_res[i] , id::CHI, 1));
		chi_torsion_ids.push_back( chi_ID );
		init_torsion = pose.torsion( chi_ID );
		MC_OneTorsionOP chi_torsion( new MC_OneTorsion( chi_ID, pose.torsion( chi_ID)));
		if ( free_res.has_value( sample_res[i] ) ) {
			chi_torsion->set_angle_range( init_torsion - options.angle_range_free_chi() , init_torsion + options.angle_range_free_chi() );
		} else {
			chi_torsion->set_angle_range( init_torsion - options.angle_range_chi() , init_torsion + options.angle_range_chi() );
		}
		// TODO -- need to include sugar pucker sampler, perhaps even with KIC closure.
		// TODO -- or do that pucker sampling in MultiSuite below. --  rhiju, 2017
		chi_torsion->set_update_tolerance( 1.0e-5 );
		chi_sampler->add_rotamer( chi_torsion );
	}
	chi_sampler->set_name( "Chi" );
	chi_sampler->set_update_pose( pose_cop );
	chi_sampler->init();

	//Set up the regular torsion samplers, necessary for free energy comparisons
	//Don't sample chi torsions here because there is a separate chi sampler
	MC_RNA_MultiSuiteOP standard_bb_sampler( new MC_RNA_MultiSuite );
	for ( Size i = 1; i <= sample_res.size(); ++i ) {
		if ( sample_res[i] > 1 && ( i == 1 || ( i > 1 && sample_res[i] != sample_res[i-1]+1 ) ) ) {
			// make sure we create sampler for [i]-1, i.e. "takeoff" for any loop
			MC_RNA_SuiteOP suite_sampler_1( new MC_RNA_Suite( sample_res[i] - 1 ) );
			suite_sampler_1->set_init_from_pose( pose );
			suite_sampler_1->set_sample_bb( true );
			suite_sampler_1->set_sample_lower_nucleoside( false );
			suite_sampler_1->set_sample_upper_nucleoside( false );
			suite_sampler_1->set_pucker_flip_rate( 0 );
			suite_sampler_1->set_sample_near_a_form( false );
			if ( free_res.has_value( sample_res[i] ) ) {
				suite_sampler_1->set_angle_range_from_init( options.angle_range_free_bb() );
			} else {
				suite_sampler_1->set_angle_range_from_init( options.angle_range_bb() );
			}
			standard_bb_sampler->add_rotamer( suite_sampler_1 );
		}

		if ( sample_res[i] >= pose.total_residue() ) continue;
		MC_RNA_SuiteOP suite_sampler( new MC_RNA_Suite( sample_res[i] ) );
		suite_sampler->set_init_from_pose( pose );
		suite_sampler->set_sample_bb( true );
		suite_sampler->set_sample_lower_nucleoside( false );
		suite_sampler->set_sample_upper_nucleoside( false );
		suite_sampler->set_pucker_flip_rate( 0 );
		suite_sampler->set_sample_near_a_form( false );
		// double-check this
		if ( free_res.has_value( sample_res[i] ) || ( i > 1 && free_res.has_value( sample_res[i-1] ) ) ) {
			suite_sampler->set_angle_range_from_init( options.angle_range_free_bb() );
		} else {
			suite_sampler->set_angle_range_from_init( options.angle_range_bb() );
		}
		standard_bb_sampler->add_rotamer( suite_sampler );
	}
	standard_bb_sampler->set_name( "Standard backbone" );
	standard_bb_sampler->set_do_no_op_random( true );
	standard_bb_sampler->set_update_pose( pose_cop );
	standard_bb_sampler->init();

	MC_RNA_OneJumpOP jump_sampler;
	if ( options.sample_jump() ) {
		jump_sampler = initialize_jump_sampler( pose, 1, options );
		jump_sampler->set_name( "Jump" );
	}

	MC_LoopOP loop_sampler( new MC_Loop );
	for ( Size n = 1; n <= 10; n++ ) {
	 	if ( (n % 10) == 0 ) {
	 		loop_sampler->add_rotamer( standard_bb_sampler );
	 	} else if ( (n % 2) == 0 ) {
	 		loop_sampler->add_rotamer( internal_bb_sampler );
		} else if ( jump_sampler != 0 && ( n % 5 ) == 0 ) {
			loop_sampler->add_rotamer( jump_sampler );
		} else {
			loop_sampler->add_rotamer( chi_sampler );
		}
	}

	return loop_sampler;
}

////////////////////////////////////////////////////////////////////////////
MC_RNA_OneJumpOP
initialize_jump_sampler( core::pose::Pose const & pose,
												 core::Size const & num_jump,
												 options::RECCES_Options const & options )
{

	MC_RNA_OneJumpOP jump_sampler( new MC_RNA_OneJump( pose, num_jump ) );
	jump_sampler->set_translation_mag( options.base_pair_translation_mag() );
	jump_sampler->set_rotation_mag( options.base_pair_rotation_mag() );
	jump_sampler->set_rmsd_cutoff( options.base_pair_rmsd_cutoff() );
	return jump_sampler;
}

} //sampler
} //recces
} //protocols
