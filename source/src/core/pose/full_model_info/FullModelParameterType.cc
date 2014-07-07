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
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/full_model_info/FullModelParameterType.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.pose.full_model_info.FullModelParameterType" );

namespace core {
namespace pose {
namespace full_model_info {


	void
	initialize_parameters( FullModelParameters & full_model_parameters ){
		full_model_parameters.set_parameter_as_res_list( CALC_RMS,           utility::vector1< Size >() );
		full_model_parameters.set_parameter_as_res_list( CUTPOINT_OPEN,      utility::vector1< Size >() );
		full_model_parameters.set_parameter_as_res_list( FIXED_DOMAIN,       utility::vector1< Size >() );
		full_model_parameters.set_parameter_as_res_list( EXTRA_MINIMIZE,     utility::vector1< Size >() );
		full_model_parameters.set_parameter_as_res_list( SAMPLE,             utility::vector1< Size >() );
		full_model_parameters.set_parameter_as_res_list( WORKING,            utility::vector1< Size >() );
		full_model_parameters.set_parameter_as_res_list( RNA_SYN_CHI,        utility::vector1< Size >() );
		full_model_parameters.set_parameter_as_res_list( RNA_TERMINAL,       utility::vector1< Size >() );
	}

} //full_model_info
} //pose
} //core
