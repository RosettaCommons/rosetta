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
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/select/residue_selector/ResiduePropertySelector.hh>

#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/carbohydrates.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/pointer/owning_ptr.hh>

static basic::Tracer TR( "protocols.carbohydrates.util" );

namespace protocols {
namespace carbohydrates {

using namespace core::pack::task;
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

} //protocols
} //carbohydrates


