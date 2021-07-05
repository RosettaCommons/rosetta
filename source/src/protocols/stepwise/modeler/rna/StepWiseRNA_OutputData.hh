// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file StepWiseRNA_OutputData.hh
/// @brief
/// @details
///
/// @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_OutputData_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_OutputData_HH


#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <core/io/silent/BinarySilentStruct.hh> //Sept 26, 2011 Parin S.

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <numeric/xyzMatrix.hh>
#include <string>
#include <core/io/silent/SilentFileData.fwd.hh>


typedef  numeric::xyzMatrix< core::Real > Matrix;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {


core::io::silent::BinarySilentStruct
get_binary_rna_silent_struct_safe( core::pose::Pose const & const_pose, std::string const & tag, std::string const & silent_file );

core::io::silent::BinarySilentStruct
get_binary_rna_silent_struct_safe_wrapper( core::pose::Pose const & const_pose, std::string const & tag, std::string const & silent_file, bool const write_score_only );

void
output_data( std::string const & silent_file, std::string const & tag, bool const write_score_only, core::pose::Pose const & pose, core::pose::PoseCOP native_poseCOP, working_parameters::StepWiseWorkingParametersCOP working_parameters_, bool const NAT_rmsd = true );

void
output_data( core::io::silent::SilentFileData& silent_file_data, std::string const & silent_file, std::string const & tag, bool const write_score_only, core::pose::Pose const & pose, core::pose::PoseCOP native_poseCOP, working_parameters::StepWiseWorkingParametersCOP working_parameters_, bool const NAT_rmsd = true );


} //rna
} //modeler
} //stepwise
} //protocols

#endif
