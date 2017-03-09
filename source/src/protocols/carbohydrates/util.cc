// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/util.cc
/// @brief Util functions for Carbohydrates.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/carbohydrates.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>


namespace protocols {
namespace carbohydrates {
	using namespace core::pack::task;
	using namespace basic::options;
	
core::pack::task::TaskFactoryOP
get_all_glycans_and_neighbor_res_task_factory(utility::vector1< bool > const & glycan_positions, core::Real pack_distance, bool read_resfile) {

	using namespace core::pack::task::operation;
	using namespace core::select::residue_selector;

	TaskFactoryOP tf = TaskFactoryOP( new TaskFactory());
	tf->push_back(InitializeFromCommandlineOP( new InitializeFromCommandline));

	//If a resfile is provided, we just use that and get out.
	if ( read_resfile && option[ OptionKeys::packing::resfile ].user() ) {
		tf->push_back( ReadResfileOP( new ReadResfile()) );
	} else {
		NeighborhoodResidueSelectorOP neighbor_selector = NeighborhoodResidueSelectorOP( new NeighborhoodResidueSelector(glycan_positions, pack_distance, true /* include focus */));
		PreventRepackingRLTOP prevent_repacking = PreventRepackingRLTOP( new PreventRepackingRLT());

		OperateOnResidueSubsetOP subset_op = OperateOnResidueSubsetOP( new OperateOnResidueSubset( prevent_repacking, neighbor_selector, true /* flip */));
		tf->push_back( subset_op );
		tf->push_back( RestrictToRepackingOP( new RestrictToRepacking()));

	}
	return tf;

}



} //protocols
} //carbohydrates


