// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/full_model_info/FullModelParameterType.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/full_model_info/FullModelParameterType.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "core.pose.full_model_info.FullModelParameterType" );

namespace core {
namespace pose {
namespace full_model_info {

std::map< FullModelParameterType, std::string> full_model_parameter_type_name;

///////////////////////////////////////////////////////////////////////////////////////////
void
initialize_parameters( FullModelParameters & full_model_parameters ){
	for ( Size n = 1; n < LAST_TYPE; n++ ) {
		FullModelParameterType type = static_cast< FullModelParameterType >( n );
		//   if ( type == CHAIN ) continue; // special case, must match CUTPOINT_OPEN
		full_model_parameters.set_parameter_as_res_list( type,           utility::vector1< Size >() );
	}
	utility::vector1< Size > working_res;
	for ( Size n = 1; n <= full_model_parameters.full_sequence().size(); n++ ) working_res.push_back( n );
	full_model_parameters.set_parameter_as_res_list( WORKING,            working_res );
}

///////////////////////////////////////////////////////////////////////////////////////////
void
initialize_full_model_parameter_type_name(){
	static bool init( false );
	if ( !init ) {
		full_model_parameter_type_name[ CALC_RMS ]        = "CALC_RMS";
		full_model_parameter_type_name[ CUTPOINT_OPEN ]   = "CUTPOINT_OPEN";
		full_model_parameter_type_name[ DOCK_DOMAIN ]     = "DOCK_DOMAIN";
		full_model_parameter_type_name[ EXTRA_MINIMIZE ]  = "EXTRA_MINIMIZE";
		full_model_parameter_type_name[ EXTRA_MINIMIZE_JUMP ] = "EXTRA_MINIMIZE_JUMP";
		full_model_parameter_type_name[ FIXED_DOMAIN ]    = "FIXED_DOMAIN";
		full_model_parameter_type_name[ INPUT_DOMAIN ]    = "INPUT_DOMAIN";
		full_model_parameter_type_name[ JUMP ]            = "JUMP";
		full_model_parameter_type_name[ PREFERRED_ROOT ]  = "PREFERRED_ROOT";
		full_model_parameter_type_name[ SAMPLE ]          = "SAMPLE";
		full_model_parameter_type_name[ WORKING ]         = "WORKING";
		full_model_parameter_type_name[ RNA_SAMPLE_SUGAR ]     = "RNA_SAMPLE_SUGAR";
		full_model_parameter_type_name[ RNA_SYN_CHI ]     = "RNA_SYN_CHI";
		full_model_parameter_type_name[ RNA_ANTI_CHI ]     = "RNA_ANTI_CHI";
		full_model_parameter_type_name[ RNA_BULGE ]       = "RNA_BULGE";
		full_model_parameter_type_name[ RNA_NORTH_SUGAR ] = "RNA_NORTH_SUGAR";
		full_model_parameter_type_name[ RNA_SOUTH_SUGAR ] = "RNA_SOUTH_SUGAR";
		full_model_parameter_type_name[ RNA_TERMINAL ]    = "RNA_TERMINAL";
		full_model_parameter_type_name[ DISULFIDE ]      = "DISULFIDE";
		init = true;
	}
}
///////////////////////////////////////////////////////////////////////////////////////////
std::string
to_string( FullModelParameterType const & full_model_parameter_type ){
	initialize_full_model_parameter_type_name();
	return full_model_parameter_type_name[ full_model_parameter_type ];
}
///////////////////////////////////////////////////////////////////////////////////////////
FullModelParameterType
full_model_parameter_type_from_string( std::string const & name ){
	initialize_full_model_parameter_type_name();
	FullModelParameterType full_model_parameter_type( NO_TYPE );
	for ( std::map< FullModelParameterType, std::string>::const_iterator it = full_model_parameter_type_name.begin(),
			end = full_model_parameter_type_name.end(); it != end; ++it ) {
		if ( it->second == name ) {
			full_model_parameter_type = it->first;
		}
	}
	runtime_assert( full_model_parameter_type != NO_TYPE);
	return full_model_parameter_type;
}

} //full_model_info
} //pose
} //core
