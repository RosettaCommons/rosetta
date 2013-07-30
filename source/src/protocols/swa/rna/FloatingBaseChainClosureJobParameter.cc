// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/FloatingBaseChainClosureJobParameter.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/swa/rna/FloatingBaseChainClosureJobParameter.hh>

#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>

#include <core/pose/Pose.hh>>
#include <core/scoring/rna/RNA_Util.hh>

#include <core/types.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.swa.FloatingBaseChainClosureJobParameter" );

namespace protocols {
namespace swa {
namespace rna {

	FloatingBaseChainClosureJobParameter::FloatingBaseChainClosureJobParameter(core::Size const input_moving_res, core::Size const input_reference_res):
			sample_sugar(true),
			moving_res ( input_moving_res ),
			reference_res( input_reference_res ),
			moving_res_pucker_state( ALL ),
			bulge_res_pucker_state( ALL ),
			moving_res_base_state( BOTH ),
			bulge_res_base_state( BOTH )
		{

			PDL.clear(); //pose_data_list
			Is_prepend= (moving_res>reference_res) ? false: true;
			bulge_res= (Is_prepend) ? reference_res-1: reference_res+1;
			bulge_suite = (Is_prepend) ? bulge_res : bulge_res-1;
			five_prime_chain_break= (Is_prepend) ? moving_res: moving_res-1;

			////////////////////////Consistency check!////////////////////////////////////////
			if(Is_prepend){
				if(moving_res+2!=reference_res){
					std::cout << "moving_res= " << moving_res << " reference_res= " << reference_res;
					utility_exit_with_message("prepend, but moving_res+2!=reference_res!");
				}
			}else{
				if(moving_res-2!=reference_res){
					std::cout << "moving_res= " << moving_res << " reference_res= " << reference_res;
					utility_exit_with_message("append, but moving_res-2!=reference_res!");
				}
			}

		}

	FloatingBaseChainClosureJobParameter::FloatingBaseChainClosureJobParameter():
			sample_sugar( false)
		{
			PDL.clear(); //pose_data_list
		};

	FloatingBaseChainClosureJobParameter::~FloatingBaseChainClosureJobParameter(){};


	void
	FloatingBaseChainClosureJobParameter::check_compatibility(core::Size const nres) const{
			using namespace ObjexxFCL;

			if(moving_res<1 || moving_res> nres) utility_exit_with_message( "moving_res<1 || moving_res> nres. moving_res= " + string_of(moving_res) );
			if(bulge_res <1 || bulge_res> nres) utility_exit_with_message( "bulge_res <1 || bulge_res> nres. bulge_res= " + string_of(bulge_res) );
			if(reference_res<1 || reference_res> nres) utility_exit_with_message( "reference_res<1 || reference_res> nres. reference_res= " + string_of(reference_res) );

			//Should check here that moving_res contain virtual ribose and bulge_res is a virtual_rna_residue?
	}

	void
	FloatingBaseChainClosureJobParameter::set_base_and_pucker_state(core::pose::Pose const & pose, StepWiseRNA_JobParametersCOP const & JP){

			////////////////////////June 02, 2011 Add BaseState and PuckerState information///
			moving_res_pucker_state=ALL;
			if(  (JP->working_force_north_ribose_list()).has_value(moving_res) ) moving_res_pucker_state=NORTH;
			if(  (JP->working_force_south_ribose_list()).has_value(moving_res) ) moving_res_pucker_state=SOUTH;

			moving_res_base_state = (core::scoring::rna::is_purine( pose.residue( moving_res ) ) ) ? BOTH: ANTI;
			if(  (JP->working_force_syn_chi_res_list()).has_value(moving_res) ) moving_res_base_state=SYN;

			bulge_res_pucker_state=ALL;
			if(  (JP->working_force_north_ribose_list()).has_value(bulge_res) ) bulge_res_pucker_state=NORTH;
			if(  (JP->working_force_south_ribose_list()).has_value(bulge_res) ) bulge_res_pucker_state=SOUTH;

			bulge_res_base_state = (core::scoring::rna::is_purine( pose.residue( bulge_res ) ) ) ? BOTH: ANTI;
			if(  (JP->working_force_syn_chi_res_list()).has_value(bulge_res) ) bulge_res_base_state=SYN;

			////////////////////////Print data!////////////////////////////////////////
			std::cout << "FloatingBaseChainClosureJobParameter: " << std::endl;
			Output_boolean(" Is_prepend= " , Is_prepend, TR );
			std::cout << " reference_res=" << reference_res;

			std::cout <<" moving_res=" << moving_res;
			print_base_state("|base_state=", moving_res_base_state );
			print_ribose_pucker_state("|pucker_state=", moving_res_pucker_state);

			std::cout << " bulge_res= " << bulge_res;
			print_base_state("|base_state=", bulge_res_base_state );
			print_ribose_pucker_state("|pucker_state=", bulge_res_pucker_state);

			std::cout << " bulge_suite= " << bulge_suite << " five_prime_chain_break= " << five_prime_chain_break;

	}

} //rna
} //swa
} //protocols
