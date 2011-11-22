// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file MembraneTopologyClaimer
/// @brief membrane topology
/// @author Yifan Song

// Unit Headers
#include <protocols/topology_broker/MembraneTopologyClaimer.hh>
#include <core/scoring/MembraneTopology.hh>
#include <protocols/moves/PoseMembraneRigidBodyMover.hh>

// Project headers
#include <basic/options/option.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>


// Package Headers
#include <protocols/topology_broker/DofClaim.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/moves/MoverContainer.hh>

// option key includes
#include <basic/options/keys/membrane.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>

#include <protocols/jobdist/Jobs.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


static basic::Tracer TR("protocols.topo_broker.membrane_topology",basic::t_info);
static numeric::random::RandomGenerator RG(786137);

namespace protocols {
namespace topology_broker {

using namespace core;


MembraneTopologyClaimer::MembraneTopologyClaimer() {}

MembraneTopologyClaimer::MembraneTopologyClaimer( pose::Pose const& input_pose ) :
	input_pose_(input_pose)
{}

void
MembraneTopologyClaimer::add_mover(
    moves::RandomMover& random_mover,
		core::pose::Pose const& pose,
		abinitio::StageID /*stageID,  abinitio sampler stage */,
		core::scoring::ScoreFunction const& /*scorefxn*/,
		core::Real /*progress  progress within stage */
)
{
	if ( basic::options::option[basic::options::OptionKeys::membrane::fixed_membrane] ) {
		//moves::MoverOP move_pose_to_membrane  =	new moves::MovePoseToMembraneCenterMover;
		//core::Real move_pose_to_membrane_weight(1.0);
		//random_mover.add_mover( move_pose_to_membrane, move_pose_to_membrane_weight);

		moves::MoverOP membrane_center_perturbation_mover =	new moves::MembraneCenterPerturbationMover;
		core::Real membrane_center_perturbation_weight(2.0);
		random_mover.add_mover( membrane_center_perturbation_mover, membrane_center_perturbation_weight );

		// TR << pose.fold_tree();
		// if pose is symmetric or option has symmetry, no rotation trial
		if ( !core::pose::symmetry::is_symmetric( pose ) ) {
			moves::MoverOP membrane_normal_perturbation_mover =	new moves::MembraneNormalPerturbationMover;
			core::Real membrane_normal_perturbation_weight(5.0);
			random_mover.add_mover( membrane_normal_perturbation_mover, membrane_normal_perturbation_weight);
		}
	}
}

/*
void MembraneTopologyClaimer::new_decoy( core::pose::Pose const& pose ) {
	//get rigid-body orientation of jump and safe for later use in initialize_dofs();

}
*/

void MembraneTopologyClaimer::initialize_dofs(
											  core::pose::Pose& pose,
											  DofClaims const& /*init_dofs*/,
											  DofClaims& /*failed_to_init*/ ) {
	using basic::options::option;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;
	runtime_assert ( option[in::file::spanfile].user() );

	std::string const spanfile = option[ in::file::spanfile ]();
	core::scoring::MembraneTopologyOP topologyOP = new core::scoring::MembraneTopology;
	pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY, topologyOP );
	core::scoring::MembraneTopology & topology=*( static_cast< core::scoring::MembraneTopology * >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY )() ));
	topology.initialize(spanfile);

	if ( basic::options::option[basic::options::OptionKeys::membrane::fixed_membrane] ) {
		if ( !core::pose::symmetry::is_symmetric( pose ) ) {
			addVirtualResAsRootMembrane(pose);
		}
	}
}

// void MembraneTopologyClaimer::generate_sequence_claims( DofClaims& /* new_claims */) {
// 	new_claims.push_back( new SequenceClaim( this, 1, 1, label(), DofClaim::INIT /* for now... eventually CAN_INIT ? */ ) );
// }

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Adds a virtual residue to a pose as the root.
/// Jump is to a residue in the middle of a transmembrane segment.
void MembraneTopologyClaimer::addVirtualResAsRootMembrane( core::pose::Pose & pose ) {
	int nres = pose.total_residue();

	core::scoring::MembraneTopology const & topology( core::scoring::MembraneTopology_from_pose(pose) );

	// return if the pose is empty (otherwise will segfault)
	if (nres == 0) {
		TR.Warning << "addVirtualResAsRootMembrane() called with empty pose!" << std::endl;
		return;
	}

	// find residue in the middle of transmembrane region
	// pick a random transmembrane segment
	Size ihelix = int(RG.uniform() * topology.tmhelix() + 1.);
	Size jump_res (int(0.5 * (topology.span_begin(ihelix) + topology.span_end(ihelix))));

	// create a virtual residue, fullatom or centroid
	bool fullatom = pose.is_fullatom();
	core::chemical::ResidueTypeSetCAP const &residue_set(
														 core::chemical::ChemicalManager::get_instance()->residue_type_set
														 ( fullatom ? core::chemical::FA_STANDARD : core::chemical::CENTROID )
														 );
	core::chemical::ResidueTypeCAPs const & rsd_type_list( residue_set->name3_map("VRT") );
	core::conformation::ResidueOP new_res( core::conformation::ResidueFactory::create_residue( *rsd_type_list[1] ) );

	// move to membrane_center if it's defined
	if ( basic::options::option[basic::options::OptionKeys::membrane::membrane_center].user() ) {
		Vector mem_center;
		mem_center.x() = basic::options::option[basic::options::OptionKeys::membrane::membrane_center]()[1];
		mem_center.y() = basic::options::option[basic::options::OptionKeys::membrane::membrane_center]()[2];
		mem_center.z() = basic::options::option[basic::options::OptionKeys::membrane::membrane_center]()[3];

		for ( Size j=1; j<= new_res->natoms(); ++j ) {
			new_res->atom(j).xyz( new_res->atom(j).xyz() + mem_center );
		}
	}

	core::pose::Pose closed_loops_pose( pose );
	closed_loops_pose.fold_tree( broker().final_fold_tree() );

	pose.append_residue_by_jump( *new_res , jump_res );
	closed_loops_pose.append_residue_by_jump( *new_res, jump_res );

	// make the virt atom the root
	kinematics::FoldTree newF( pose.fold_tree() );
	newF.reorder( nres+1 );
	TR << "addVirtualResAsRoot() setting new fold tree to " << newF << std::endl;
	TR << "   jump_res = " << jump_res << std::endl;
	pose.fold_tree( newF );

	kinematics::FoldTree finalF( closed_loops_pose.fold_tree() );
	finalF.reorder( nres+1 );
	broker().final_fold_tree() = finalF;

}

} //topology_broker
} //protocols
