// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MembraneTopologyClaimer
/// @brief membrane topology
/// @author Yifan Song

// Unit Headers
#include <protocols/topology_broker/MembraneTopologyClaimer.hh>
#include <core/scoring/MembraneTopology.hh>
#include <protocols/rigid/PoseMembraneRigidBodyMover.hh>

// Project headers
#include <basic/options/option.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/symmetry/util.hh>


// Package Headers
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/BBClaim.hh>
#include <protocols/topology_broker/claims/CutClaim.hh>
#include <protocols/topology_broker/claims/JumpClaim.hh>
#include <protocols/topology_broker/claims/LegacyRootClaim.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/moves/MoverContainer.hh>

// option key includes
#include <basic/options/keys/membrane.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/conformation/Residue.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.topo_broker.membrane_topology", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace core;


MembraneTopologyClaimer::MembraneTopologyClaimer(){}

MembraneTopologyClaimer::MembraneTopologyClaimer( pose::Pose const& input_pose ) :
	input_pose_(input_pose)
{}

bool
MembraneTopologyClaimer::claimer_builds_own_fold_tree()
{
	return true;
}

void
MembraneTopologyClaimer::set_pose_from_broker(core::pose::Pose& pose)
{
	input_pose_ = pose;
}

core::kinematics::FoldTreeOP
MembraneTopologyClaimer::get_fold_tree(core::pose::Pose& pose)
{
	return core::kinematics::FoldTreeOP( new core::kinematics::FoldTree(pose.fold_tree()) );
}


//@brief reads in spanfile, sets up membrane topology, adds virtual residue to pose as the root residue
void
MembraneTopologyClaimer::pre_process(core::pose::Pose& pose)
{
	using basic::options::option;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;
	runtime_assert ( option[in::file::spanfile].user() );

	std::string const spanfile = option[ in::file::spanfile ]();
	core::scoring::MembraneTopologyOP topologyOP( new core::scoring::MembraneTopology );
	pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY, topologyOP );
	core::scoring::MembraneTopology & topology=*( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ));
	topology.initialize(spanfile);

	if ( basic::options::option[basic::options::OptionKeys::membrane::fixed_membrane] ) {
		if ( !core::pose::symmetry::is_symmetric( pose ) ) {
			addVirtualResAsRootMembrane(pose);
		}
	}
}

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
		//moves::MoverOP move_pose_to_membrane  = new moves::MovePoseToMembraneCenterMover;
		//core::Real move_pose_to_membrane_weight(1.0);
		//random_mover.add_mover( move_pose_to_membrane, move_pose_to_membrane_weight);

		moves::MoverOP membrane_center_perturbation_mover( new rigid::MembraneCenterPerturbationMover );
		core::Real membrane_center_perturbation_weight(2.0);
		random_mover.add_mover( membrane_center_perturbation_mover, membrane_center_perturbation_weight );

		// TR << pose.fold_tree();
		// if pose is symmetric or option has symmetry, no rotation trial
		if ( !core::pose::symmetry::is_symmetric( pose ) ) {
			moves::MoverOP membrane_normal_perturbation_mover( new rigid::MembraneNormalPerturbationMover );
			core::Real membrane_normal_perturbation_weight(5.0);
			random_mover.add_mover( membrane_normal_perturbation_mover, membrane_normal_perturbation_weight);
		}
	}
}

void
MembraneTopologyClaimer::generate_claims( claims::DofClaims& dof_claims )
{
	core::pose::Pose pose = broker().current_pose();
	ObjexxFCL::FArray2D_int jump_array(2,pose.fold_tree().num_jump());

	for ( core::Size jump_array_index = 1; jump_array_index <= pose.fold_tree().num_jump(); ++jump_array_index ) {
		jump_array(1,jump_array_index) = pose.fold_tree().jump_edge(jump_array_index).start();
		jump_array(2,jump_array_index) = pose.fold_tree().jump_edge(jump_array_index).stop();
	}

	if ( TR.Trace.visible() ) { pose.fold_tree().show(TR.Trace);}
	for ( core::Size i=1; i<=pose.fold_tree().nres(); ++i ) {
		if ( i == pose.fold_tree().root() ) {
			dof_claims.push_back(claims::DofClaimOP( new claims::LegacyRootClaim(get_self_weak_ptr(),i,claims::DofClaim::CAN_INIT) ));
			dof_claims.push_back(claims::DofClaimOP( new claims::BBClaim(get_self_weak_ptr(),i,claims::DofClaim::CAN_INIT) ));
		} else if ( pose.residue(i).is_virtual_residue() || pose.residue(i).name3() == "XXX" || pose.residue(i).name3() == "VRT" ) {
			dof_claims.push_back(claims::DofClaimOP( new claims::BBClaim(get_self_weak_ptr(),i,claims::DofClaim::CAN_INIT) ));
		}
	}
	for ( Size jump_num = 1; jump_num <= pose.fold_tree().num_jump(); ++jump_num ) {
		dof_claims.push_back(claims::DofClaimOP( new claims::JumpClaim(get_self_weak_ptr(),jump_array(1,jump_num),jump_array(2,jump_num),claims::DofClaim::CAN_INIT) ));
	}
}

void MembraneTopologyClaimer::initialize_dofs( core::pose::Pose& pose, claims::DofClaims const& /*init_dofs*/, claims::DofClaims& /*failed_to_init*/) {
	using basic::options::option;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;
	runtime_assert ( option[in::file::spanfile].user() );

	std::string const spanfile = option[ in::file::spanfile ]();
	core::scoring::MembraneTopologyOP topologyOP( new core::scoring::MembraneTopology );
	pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY, topologyOP );
	core::scoring::MembraneTopology & topology=*( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ));
	topology.initialize(spanfile);

	if ( basic::options::option[basic::options::OptionKeys::membrane::fixed_membrane] ) {
		if ( !core::pose::symmetry::is_symmetric( pose ) ) {
			addVirtualResAsRootMembrane(pose);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Adds a virtual residue to a pose as the root.
/// Jump is to a residue in the middle of a transmembrane segment.
void MembraneTopologyClaimer::addVirtualResAsRootMembrane( core::pose::Pose & pose ) {
	int nres = pose.size();
	core::scoring::MembraneTopology const & topology( core::scoring::MembraneTopology_from_pose(pose) );

	// return if the pose is empty (otherwise will segfault)
	if ( nres == 0 ) {
		TR.Warning << "addVirtualResAsRootMembrane() called with empty pose!" << std::endl;
		return;
	}

	//only add the virtual residue as the root once
	if ( pose.residue(pose.fold_tree().root()).name3() == "XXX" ) {
		return;
	}

	// find residue in the middle of transmembrane region
	// pick a random transmembrane segment
	Size ihelix = int(numeric::random::rg().uniform() * topology.tmhelix() + 1.);
	core::Size jump_res = (int(0.5 * (topology.span_begin(ihelix) + topology.span_end(ihelix))));

	// create a virtual residue, of the appropriate type
	core::chemical::ResidueTypeSetCOP residue_set( pose.residue_type_set_for_pose() );
	core::chemical::ResidueTypeCOP rsd_type( residue_set->get_representative_type_name3("VRT") );
	core::conformation::ResidueOP new_virt_res( core::conformation::ResidueFactory::create_residue( *rsd_type ) );

	// move to membrane_center if it's defined
	if ( basic::options::option[basic::options::OptionKeys::membrane::membrane_center].user() ) {
		Vector mem_center;
		mem_center.x() = basic::options::option[basic::options::OptionKeys::membrane::membrane_center]()[1];
		mem_center.y() = basic::options::option[basic::options::OptionKeys::membrane::membrane_center]()[2];
		mem_center.z() = basic::options::option[basic::options::OptionKeys::membrane::membrane_center]()[3];

		for ( Size j=1; j<= new_virt_res->natoms(); ++j ) {
			new_virt_res->atom(j).xyz( new_virt_res->atom(j).xyz() + mem_center );
		}
	}

	core::pose::Pose closed_loops_pose( pose );
	//TR << "closed_loops_pose FoldTree is:  ";
	//closed_loops_pose.fold_tree().show(TR);
	if ( broker().does_final_fold_tree_exist() ) {
		closed_loops_pose.fold_tree( broker().final_fold_tree() );
	}
	pose.append_residue_by_jump( *new_virt_res , jump_res );
	closed_loops_pose.append_residue_by_jump( *new_virt_res, jump_res );

	if ( TR.Trace.visible() ) {
		TR.Trace << "before fold tree...jump_res:  " << jump_res << " virt res:  " << new_virt_res->seqpos() << " " << new_virt_res->name() << std::endl;
		TR.Trace << "before reordering...beginning of fold tree: " << pose.fold_tree().begin()->start() << std::endl;
	}


	// make the virt atom the root
	kinematics::FoldTree newF( pose.fold_tree() );
	newF.reorder( nres+1 );
	TR << "addVirtualResAsRoot() setting new fold tree to " << newF << std::endl;
	TR << "   jump_res = " << jump_res << std::endl;
	pose.fold_tree( newF );

	kinematics::FoldTree finalF( closed_loops_pose.fold_tree() );
	finalF.reorder( nres+1 );
	if ( broker().does_final_fold_tree_exist()==true ) {
		broker().final_fold_tree() = finalF;
	} else {
		TR << "Broker says final fold tree does not exist!" << std::endl;
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "after reordering...beginning of fold tree: " << pose.fold_tree().begin()->start() << std::endl;
		TR.Debug << "position of virtual res:  " << new_virt_res->seqpos() << std::endl;

		TR.Debug << "root of fold tree:  " << pose.fold_tree().root() << " is virtual? " <<
			pose.residue(pose.fold_tree().root()).is_virtual_residue() << " is residue type:  " <<
			pose.residue(pose.fold_tree().root()).name() << std::endl;
		core::Size root(pose.fold_tree().root());
		for ( core::Size i = 1; i <= pose.residue(root).natoms(); i++ ) {
			std::string atom_name = pose.residue(root).atom_name(i);
			TR.Debug << pose.residue(root).name() << " " << atom_name << " " << pose.residue(root).atom(atom_name) << std::endl;
		}
	}
}

void
MembraneTopologyClaimer::build_fold_tree(core::pose::Pose& pose, core::kinematics::FoldTree& fold_tree_in)
{
	pose.fold_tree(fold_tree_in);
}

} //topology_broker
} //protocols
