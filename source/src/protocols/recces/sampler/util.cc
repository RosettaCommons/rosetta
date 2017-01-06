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
#include <protocols/recces/params/RECCES_Parameters.hh>
#include <protocols/recces/options/RECCES_Options.hh>
#include <protocols/farna/secstruct/RNA_SecStructLegacyInfo.hh>
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

protocols::recces::sampler::MC_CombOP
get_recces_turner_sampler( pose::Pose const & pose,
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
		TorsionID chi_ID (TorsionID( sample_res[i] , id::CHI, 1));
		chi_torsion_ids.push_back( chi_ID );
		init_torsion = pose.torsion( chi_ID );
		MC_OneTorsionOP chi_torsion( new MC_OneTorsion( chi_ID, pose.torsion( chi_ID)));

		if ( free_res.has_value( sample_res[i] ) ) {
			chi_torsion->set_angle_range( init_torsion - options.angle_range_free_chi() , init_torsion + options.angle_range_free_chi() );
		} else {
			chi_torsion->set_angle_range( init_torsion - options.angle_range_chi() , init_torsion + options.angle_range_chi() );
		}
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
		if ( i == 1 || ( i > 1 && sample_res[i] != sample_res[i-1]+1 ) ) {
			//create samplers for [i]-1 and [i]
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

		MC_RNA_SuiteOP suite_sampler( new MC_RNA_Suite( sample_res[i] ) );
		suite_sampler->set_init_from_pose( pose );
		suite_sampler->set_sample_bb( true );
		suite_sampler->set_sample_lower_nucleoside( false );
		suite_sampler->set_sample_upper_nucleoside( false );
		suite_sampler->set_pucker_flip_rate( 0 );
		suite_sampler->set_sample_near_a_form( false );
		if ( free_res.has_value( sample_res[i] ) ) {
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

	MC_LoopOP loop_sampler( new MC_Loop );
	for ( Size n = 1; n <= 10; n++ ) {
	 	if ( (n % 10) == 0 ) {
	 		loop_sampler->add_rotamer( standard_bb_sampler );
	 	} else if ( (n % 2) == 0 ) {
	 		loop_sampler->add_rotamer( internal_bb_sampler );
		} else {
	 		loop_sampler->add_rotamer( chi_sampler );
		}
	}

	return loop_sampler;
}

} //sampler
} //recces
} //protocols
