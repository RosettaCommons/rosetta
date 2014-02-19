// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/enumerate/rna/sugar/SugarModeling.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/enumerate/rna/sugar/SugarModeling.hh>

#include <protocols/stepwise/enumerate/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/types.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.rna.SugarModeling" );

using namespace core::chemical::rna;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {
namespace sugar {

	SugarModeling::SugarModeling( core::Size const input_moving_res, core::Size const input_reference_res ):
			sample_sugar( true ),
			moving_res ( input_moving_res ),
			reference_res( input_reference_res ),
			moving_res_pucker_state( WHATEVER ),
			bulge_res_pucker_state( WHATEVER ),
			moving_res_base_state( WHATEVER ),
			bulge_res_base_state( WHATEVER )
		{

			pose_list.clear();
			is_prepend = ( moving_res <= reference_res );
			bulge_res  = ( is_prepend ) ? reference_res - 1: reference_res + 1;
			bulge_suite = ( is_prepend ) ? bulge_res : bulge_res - 1;
			five_prime_chain_break = ( is_prepend ) ? moving_res: moving_res - 1;
		}

	SugarModeling::SugarModeling():
		sample_sugar( false )
	{
		pose_list.clear(); //pose_data_list
	}

	SugarModeling::~SugarModeling(){}


	void
	SugarModeling::check_compatibility( core::Size const nres ) const{
		using namespace ObjexxFCL;
		runtime_assert( moving_res >= 1 && moving_res <= nres );
		runtime_assert( reference_res >= 1 && reference_res <= nres );
	}

	void
	SugarModeling::set_base_and_pucker_state( core::pose::Pose const & pose, StepWiseRNA_JobParametersCOP const & JP ){

			////////////////////////June 02, 2011 Add core::Size and core::Size information///
			moving_res_pucker_state = WHATEVER;
			if (  ( JP->working_force_north_sugar_list() ).has_value( moving_res ) ) moving_res_pucker_state = NORTH;
			if (  ( JP->working_force_south_sugar_list() ).has_value( moving_res ) ) moving_res_pucker_state = SOUTH;

			moving_res_base_state = ( core::chemical::rna::is_purine( pose.residue( moving_res ) ) ) ? WHATEVER: ANTI;
			if (  ( JP->working_force_syn_chi_res_list() ).has_value( moving_res ) ) moving_res_base_state = SYN;

			bulge_res_pucker_state = WHATEVER;
			if (  ( JP->working_force_north_sugar_list() ).has_value( bulge_res ) ) bulge_res_pucker_state = NORTH;
			if (  ( JP->working_force_south_sugar_list() ).has_value( bulge_res ) ) bulge_res_pucker_state = SOUTH;

			bulge_res_base_state = ( core::chemical::rna::is_purine( pose.residue( bulge_res ) ) ) ? WHATEVER: ANTI;
			if (  ( JP->working_force_syn_chi_res_list() ).has_value( bulge_res ) ) bulge_res_base_state = SYN;

			////////////////////////Print data!////////////////////////////////////////
			TR.Debug << "SugarModeling: " << std::endl;
			output_boolean( " is_prepend = ", is_prepend, TR.Debug );
			TR.Debug << " reference_res = " << reference_res;

			TR.Debug << " moving_res = " << moving_res;
			print_base_state( "|base_state = ", moving_res_base_state, TR.Debug );
			print_sugar_pucker_state( "|pucker_state = ", moving_res_pucker_state, TR.Debug );

			TR.Debug << " bulge_res = " << bulge_res;
			print_base_state( "|base_state = ", bulge_res_base_state, TR.Debug );
			print_sugar_pucker_state( "|pucker_state = ", bulge_res_pucker_state, TR.Debug );

			TR.Debug << " bulge_suite = " << bulge_suite << " five_prime_chain_break = " << five_prime_chain_break;
	}

} //sugar
} //rna
} //enumerate
} //stepwise
} //protocols
