// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/rna/RNA_StepWiseSamplerSamplerUtil.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampler/rna/util.hh>
#include <protocols/stepwise/sampler/rna/RNA_KIC_Sampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_SuiteStepWiseSampler.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerBase.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/rna/util.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.sampler.rna.RNA_StepWiseSamplerSamplerUtil" );

using namespace core;
using namespace protocols::stepwise::modeler::rna;
using namespace protocols::stepwise::modeler;

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

//////////////////////////////////////////////////////////////////
// Fang-Chieh Chou nicely refactored the rotamer modeler
//////////////////////////////////////////////////////////////////
StepWiseSamplerBaseOP
setup_sampler( pose::Pose const & pose,
							 modeler::options::StepWiseModelerOptionsCOP options,
							 working_parameters::StepWiseWorkingParametersCOP working_parameters,
							 bool const build_pose_from_scratch,
							 bool const kic_modeler,
							 bool const close_chain ) {

	using namespace chemical::rna;
	using namespace pose::rna;

	/////Load in constants being used/////
	utility::vector1<Size> const & working_moving_suite_list(	working_parameters->working_moving_suite_list() );
	utility::vector1<Size> const & syn_chi_res(	working_parameters->working_force_syn_chi_res_list() );
	utility::vector1<Size> const & north_puckers(	working_parameters->working_force_north_sugar_list() );
	utility::vector1<Size> const & south_puckers(	working_parameters->working_force_south_sugar_list() );
	bool const is_prepend_ = working_parameters->is_prepend();
	bool const is_internal_ = working_parameters->is_internal();

	runtime_assert( working_moving_suite_list.size() == 1 );
	Size const moving_suite_( working_moving_suite_list[1] );

	/////Get the base and pucker state/////
	utility::vector1<bool> sample_sugar( 2, false );
	utility::vector1<ChiState> base_state( 2, ANY_CHI );
	utility::vector1<PuckerState> pucker_state( 2, ANY_PUCKER );

	if ( build_pose_from_scratch || working_parameters->sample_both_sugar_base_rotamer()) {
		sample_sugar[1] = true;
		sample_sugar[2] = true;
	} else if ( !is_internal_  ) {
		if ( is_prepend_ ) {
			sample_sugar[1] = true;
		} else {
			sample_sugar[2] = true;
		}
	} else {
		runtime_assert( is_internal_ );
	}

	// new ... should be more robust
	if ( modeler_sugar_at_five_prime(  pose, moving_suite_ ) ) sample_sugar[ 1 ] = true;
	if ( modeler_sugar_at_three_prime( pose, moving_suite_ ) ) sample_sugar[ 2 ] = true;
	if ( working_parameters->sample_both_sugar_base_rotamer() ){
		runtime_assert( sample_sugar[1] );
		runtime_assert( sample_sugar[2] );
	}

	for ( Size i = 1; i <= 2; ++i ) {
		Size const curr_rsd( moving_suite_ + i - 1 );
		if ( sample_sugar[i] ) {
			bool is_north ( north_puckers.has_value( curr_rsd ) );
			bool is_south ( south_puckers.has_value( curr_rsd ) );
			bool is_syn ( syn_chi_res.has_value( curr_rsd ) );
			runtime_assert( !is_north || !is_south );
			if ( is_north ) pucker_state[i] = NORTH;
			if ( is_south ) pucker_state[i] = SOUTH;
			if (  is_purine( pose.residue( curr_rsd ) ) ) {
				if ( is_syn ) base_state[i] = SYN;
			} else { // pyrimidine
				if ( !options->allow_syn_pyrimidine() ){
					runtime_assert( !is_syn );
					base_state[i] = ANTI;
				}
			}
		} else {
			pucker_state[i] = assign_pucker( pose, curr_rsd );
			base_state[i] = NO_CHI;
		}
	}


	/////Set up the sampler/////
	if ( kic_modeler ) {
		runtime_assert( close_chain );

		Size const chainbreak_suite( working_parameters->five_prime_chain_break_res() );
		runtime_assert( chainbreak_suite > 0 );

		pose::PoseOP new_pose = new pose::Pose( pose ); //hard copy
		RNA_KIC_SamplerOP sampler = new RNA_KIC_Sampler(
				new_pose, moving_suite_, chainbreak_suite );
		//		runtime_assert( (moving_suite_ == chainbreak_suite + 1) || (moving_suite_ == chainbreak_suite - 1) );
		Size const which_nucleoside_to_sample = ( moving_suite_ < chainbreak_suite ) ? 2 : 1;
		Size const sample_nucleoside_res = ( moving_suite_ < chainbreak_suite ) ? (moving_suite_ + 1) : moving_suite_;
		sampler->set_sample_nucleoside( sample_nucleoside_res );
		if ( !sample_sugar[ which_nucleoside_to_sample ] ) {
			sampler->set_base_state( NO_CHI );
			sampler->set_pucker_state( NO_PUCKER );
		} else {
			sampler->set_base_state( base_state[ which_nucleoside_to_sample ] );
			sampler->set_pucker_state( pucker_state[ which_nucleoside_to_sample ] );
		}

		sampler->set_verbose( options->verbose() );
		sampler->set_skip_same_pucker( options->use_phenix_geo() );
		sampler->set_idealize_coord( options->use_phenix_geo() );
		sampler->set_extra_epsilon( options->sampler_extra_epsilon_rotamer() );
		sampler->set_extra_chi(	options->extra_chi() );
		sampler->set_random( options->choose_random() );
		sampler->set_fast( options->integration_test_mode() ); // overrules extra_chi, extra_epsilon; and sets bin size to 40!
		if ( options->finer_sampling_at_chain_closure() ) sampler->set_bin_size( 10 );
		sampler->init();
		return sampler;
	}

	RNA_SuiteStepWiseSamplerOP sampler = new RNA_SuiteStepWiseSampler( moving_suite_,
			pucker_state[1], pucker_state[2], base_state[1], base_state[2] );
	sampler->set_skip_same_pucker( options->use_phenix_geo() );
	sampler->set_idealize_coord( options->use_phenix_geo() );
	sampler->set_sample_nucleoside_lower( sample_sugar[1] );
	sampler->set_sample_nucleoside_upper( sample_sugar[2] );
	sampler->set_fast( options->integration_test_mode() );
	sampler->set_extra_epsilon( options->sampler_extra_epsilon_rotamer() );
	sampler->set_extra_beta( options->sampler_extra_beta_rotamer() );
	sampler->set_extra_chi(	options->extra_chi() );
	sampler->set_random( options->choose_random() );
	if ( close_chain && options->finer_sampling_at_chain_closure()  )	sampler->set_bin_size( 10 );
	sampler->init();

	return sampler;
}

/////////////////////////////////////////////////////////////////////////
bool
modeler_sugar_at_five_prime( pose::Pose const & pose,
															Size const moving_suite ) {
	using namespace core::pose::full_model_info;
	if ( moving_suite == 1 ||
			 ( pose.fold_tree().is_cutpoint( moving_suite - 1 ) &&
				 !pose.residue_type( moving_suite ).has_variant_type( chemical::CUTPOINT_UPPER ) ) ) {
		if ( full_model_info_defined( pose ) ) {
				utility::vector1< Size > const & sample_res = const_full_model_info( pose ).sample_res();
				utility::vector1< Size > const & res_list   = const_full_model_info( pose ).res_list();
				return sample_res.has_value( res_list[ moving_suite ] );
		} else {
			return true;
		}
	}
	return false;
}

/////////////////////////////////////////////////////////////////////////
bool
modeler_sugar_at_three_prime( pose::Pose const & pose,
															Size const moving_suite ) {
	using namespace core::pose::full_model_info;
	if ( (moving_suite + 1) == pose.total_residue() ||
			 ( pose.fold_tree().is_cutpoint( moving_suite + 1 ) &&
				 !pose.residue_type( moving_suite + 1 ).has_variant_type( chemical::CUTPOINT_LOWER ) ) ) {
		if ( full_model_info_defined( pose ) ) {
				utility::vector1< Size > const & sample_res = const_full_model_info( pose ).sample_res();
				utility::vector1< Size > const & res_list   = const_full_model_info( pose ).res_list();
				return sample_res.has_value( res_list[ moving_suite + 1 ] );
		} else {
			return true;
		}
	}
	return false;
}

} //rna
} //sampler
} //stepwise
} //protocols
