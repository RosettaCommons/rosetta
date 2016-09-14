// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/rna/sugar/SugarModeling.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/rna/sugar/SugarModeling.hh>

#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/output_util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/rna/util.hh>
#include <core/types.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.rna.sugar.SugarModeling" );

using namespace core::chemical::rna;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace sugar {

SugarModeling::SugarModeling( core::Size const input_moving_res, core::Size const input_reference_res ):
	sample_sugar( true ),
	moving_res ( input_moving_res ),
	reference_res( input_reference_res ),
	moving_res_pucker_state( ANY_PUCKER ),
	bulge_res_pucker_state( ANY_PUCKER ),
	moving_res_base_state( ANY_CHI ),
	bulge_res_base_state( ANY_CHI )
{
	pose_list.clear();
	is_prepend = ( moving_res <= reference_res );
	bulge_res  = ( is_prepend ) ? reference_res - 1: reference_res + 1;
	bulge_suite = ( is_prepend ) ? bulge_res : bulge_res - 1;
	five_prime_chain_break = ( is_prepend ) ? moving_res: moving_res - 1;
}

SugarModeling::SugarModeling():
	sample_sugar( false ),
	moving_res( 0 ),
	reference_res( 0 ),
	is_prepend( false ),
	bulge_res( 0 ),
	bulge_suite( 0 ),
	five_prime_chain_break( 0 ),
	moving_res_pucker_state( ANY_PUCKER ),
	bulge_res_pucker_state( ANY_PUCKER ),
	moving_res_base_state( ANY_CHI ),
	bulge_res_base_state( ANY_CHI )
{
	pose_list.clear(); //pose_data_list
}

SugarModeling &
SugarModeling::operator = ( SugarModeling const & src )
{
	sample_sugar = src.sample_sugar;
	moving_res = src.moving_res;
	reference_res = src.reference_res;
	is_prepend = src.is_prepend;
	bulge_res = src.bulge_res;
	bulge_suite = src.bulge_suite;
	five_prime_chain_break = src.five_prime_chain_break;
	moving_res_pucker_state = src.moving_res_pucker_state;
	bulge_res_pucker_state = src.bulge_res_pucker_state;
	moving_res_base_state = src.moving_res_base_state;
	bulge_res_base_state = src.bulge_res_base_state;
	pose_list = src.pose_list;

	return *this;
}

SugarModeling::~SugarModeling(){}


void
SugarModeling::check_compatibility( core::Size const nres ) const{
	using namespace ObjexxFCL;
	runtime_assert( moving_res >= 1 && moving_res <= nres );
	//runtime_assert( reference_res >= 1 && reference_res <= nres );
}

void
SugarModeling::set_base_and_pucker_state( core::pose::Pose const & pose, working_parameters::StepWiseWorkingParametersCOP const & WP ){

	////////////////////////June 02, 2011 Add core::Size and core::Size information///
	moving_res_pucker_state = ANY_PUCKER;
	if ( ( WP->working_force_north_sugar_list()  ).has_value( moving_res ) ) moving_res_pucker_state = NORTH;
	if ( ( WP->working_force_south_sugar_list()  ).has_value( moving_res ) ) moving_res_pucker_state = SOUTH;

	// note: hard-wired anti for pyrimidine
	moving_res_base_state = ( pose.residue_type( moving_res ).is_purine() ) ? ANY_CHI: ANTI;
	if ( ( WP->working_force_syn_chi_res_list()  ).has_value( moving_res ) ) moving_res_base_state = SYN;
	if ( ( WP->working_force_anti_chi_res_list() ).has_value( moving_res ) ) moving_res_base_state = ANTI;

	bulge_res_pucker_state = ANY_PUCKER;
	if ( ( WP->working_force_north_sugar_list()  ).has_value(  bulge_res ) ) bulge_res_pucker_state = NORTH;
	if ( ( WP->working_force_south_sugar_list()  ).has_value(  bulge_res ) ) bulge_res_pucker_state = SOUTH;

	// note: hard-wired anti for pyrimidine
	bulge_res_base_state = ( pose.residue_type( bulge_res ).is_purine() ) ? ANY_CHI: ANTI;
	if ( ( WP->working_force_syn_chi_res_list()  ).has_value(  bulge_res ) ) bulge_res_base_state = SYN;
	if ( ( WP->working_force_anti_chi_res_list() ).has_value(  bulge_res ) ) bulge_res_base_state = ANTI;

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
} //modeler
} //stepwise
} //protocols
