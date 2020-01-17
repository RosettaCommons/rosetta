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

#include <core/chemical/ResidueProperty.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/select/residue_selector/ResiduePropertySelector.hh>
#include <core/select/residue_selector/GlycanLayerSelector.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>
#include <core/select/residue_selector/AndResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/Pose.hh>

#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/carbohydrates.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/string_util.hh>
#include <cmath>

static basic::Tracer TR( "protocols.carbohydrates.util" );

namespace protocols {
namespace carbohydrates {

using namespace core::pack::task;
using namespace core::select::residue_selector;
using namespace basic::options;
using namespace utility::pointer;

core::pack::task::TaskFactoryOP
get_all_glycans_and_neighbor_res_task_factory(utility::vector1< bool > const & glycan_positions, core::Real pack_distance, bool read_resfile) {

	using namespace core::pack::task::operation;
	using namespace core::select::residue_selector;

	TaskFactoryOP tf = utility::pointer::make_shared< TaskFactory >();
	tf->push_back(utility::pointer::make_shared< InitializeFromCommandline >());

	//If a resfile is provided, we just use that and get out.
	if ( read_resfile && option[ OptionKeys::packing::resfile ].user() ) {
		tf->push_back( utility::pointer::make_shared< ReadResfile >() );
	} else {
		NeighborhoodResidueSelectorOP neighbor_selector = utility::pointer::make_shared< NeighborhoodResidueSelector >(glycan_positions, pack_distance, true /* include focus */);
		PreventRepackingRLTOP prevent_repacking = make_shared< PreventRepackingRLT >();

		OperateOnResidueSubsetOP subset_op = make_shared< OperateOnResidueSubset >( prevent_repacking, neighbor_selector, true /* flip */);
		tf->push_back( subset_op );

		//Skip virts
		ResiduePropertySelectorOP virt_selector = utility::pointer::make_shared<ResiduePropertySelector>();
		virt_selector->set_property( core::chemical::VIRTUAL_RESIDUE);

		OperateOnResidueSubsetOP virt_subset_op = make_shared< OperateOnResidueSubset >( prevent_repacking, virt_selector, false /* flip */);

		tf->push_back(virt_subset_op);
		tf->push_back( utility::pointer::make_shared< RestrictToRepacking >());

	}
	return tf;

}

void
run_shear_min_pack(
	minimization_packing::MinMover & min_mover,
	minimization_packing::PackRotamersMover & packer,
	simple_moves::ShearMover & shear,
	moves::MonteCarlo & mc,
	core::Size n_glycan_residues,
	core::pose::Pose & pose,
	bool use_shear)
{

	core::Size accepts = 0;
	core::Size shear_trials = 40 * n_glycan_residues; //About 3 sec for 14 glycan residues
	if ( use_shear ) {
		for ( core::Size i = 1; i <= shear_trials; ++i ) {
			shear.apply(pose);
			bool accepted = mc.boltzmann(pose);
			if ( accepted ) accepts+=1;
		}
		TR << "Shear accepts: "<< accepts << "/"<<shear_trials << std::endl;
	}

	//I've had odd times when energy is increased here.
	min_mover.apply( pose );
	bool accepted = mc.boltzmann( pose );
	TR << "min accept: " << accepted << std::endl;

	packer.apply( pose );
	accepted = mc.boltzmann( pose );
	TR << "pack accepted: " << accepted << std::endl;
}

///@brief Used for benchmarking to test even sampling of different kinematic protocols
///  Get the total number of sampling rounds for the GlycanTreeModeler protocol with previous default settings.
core::Size
get_total_rounds_for_overlap_one_layer_two(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & residue_subset,
	core::Size sampler_rounds
){

	//Matches sampling to that of the GlycanTreeSampler with default settings (window_size 2, overlap 1)
	// Used for benchmarking

	//R*nL0 + R*nLm1 + 2*R*(N-(nL0+nL1))
	//
	//R = Number of set rounds
	//N = Number of total glycans
	//nL0 = Number of glycan residues in layer 0
	//nLm1 = Number of glycan residues in last Layer (-1 index)

	//Setup
	GlycanLayerSelector layer_selector = GlycanLayerSelector();
	layer_selector.set_layer(0, 0);
	utility::vector1< bool > const layer_0 = layer_selector.apply(pose);

	core::Size max_end_layer = pose.glycan_tree_set()->get_largest_glycan_tree_layer( residue_subset );
	layer_selector.set_layer(max_end_layer, max_end_layer);
	utility::vector1< bool > const layer_m1  = layer_selector.apply(pose);

	//Variables
	core::Size const R = sampler_rounds;
	core::Size const N = count_selected(residue_subset);
	core::Size const nL0  = count_selected( AND_combine( residue_subset, layer_0 ) );
	core::Size const nLm1 = count_selected( AND_combine( residue_subset, layer_m1) );

	//Calculation
	core::Size total_rounds = (R * nL0) + (R * nLm1) + (2 * R * (N - (nL0+nLm1) ));
	return total_rounds;

}


} //protocols
} //carbohydrates


