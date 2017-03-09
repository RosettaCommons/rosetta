// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/util.hh
/// @brief Util functions for Carbohydrates.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Jason W. Labonte (JWLabonte@jhu.edu)

#ifndef INCLUDED_protocols_carbohydrates_util_hh
#define INCLUDED_protocols_carbohydrates_util_hh

#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace carbohydrates {

///@brief Get a TaskFactory of all residues in the subset and neighboring residues.
///@details
///
/// Operations:
///  InitializeFromCommandline
///  ReadResFile? If option given on cmd-line, returns TF up to this.
///  NeighborhoodResidueSelector/OperateOnResidueSubset
///  RestrictRepacking/PreventRepacking
core::pack::task::TaskFactoryOP
get_all_glycans_and_neighbor_res_task_factory(utility::vector1< bool > const & subset, core::Real pack_distance = 6.0, bool read_resfile = true);


} //protocols
} //carbohydrates


#endif	//protocols/carbohydrates_util_hh

