// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/mmtf/ExtraDataEnum.hh
/// @brief Enum definitions for pose extra data we will be storing/retrieving.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_core_io_mmtf_ExtraDataEnum_hh
#define INCLUDED_core_io_mmtf_ExtraDataEnum_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace core {
namespace io {
namespace mmtf {




/// @brief Enum definitions for pose extra data we will be storing/retrieving.
///   ALL Enums must also have their string equivalents defined in mmtf/ExtraDataEnumManager.cc
enum ExtraDataEnum{
	pose_cache_string_data = 1,
	pose_cache_real_data,
	pdb_comments,
	simple_metric_string_data,
	simple_metric_real_data,
	simple_metric_composite_string_data,
	simple_metric_composite_real_data,
	simple_metric_per_residue_string_data,
	simple_metric_per_residue_real_data,
	simple_metric_per_residue_string_output,
	simple_metric_per_residue_real_output,

	ExtraDataEnum_total=simple_metric_per_residue_real_output
};


} //core
} //io
} //mmtf



#endif //INCLUDED_core_io_mmtf_ExtraDataEnum_hh





