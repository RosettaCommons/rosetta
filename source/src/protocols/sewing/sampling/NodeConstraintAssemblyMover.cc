// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file NodeConstraintAssemblyMover.cc
///
/// @brief
/// @author Tim Jacobs

#ifdef NOT_IN_SCONS_DEPRECATED

// Unit Headers
#include <devel/sewing/sampling/NodeConstraintAssemblyMover.hh>
#include <devel/sewing/sampling/NodeConstraintAssemblyMoverCreator.hh>

// Package Headers
#include <devel/sewing/conformation/AssemblyFactory.hh>
#include <devel/sewing/conformation/ContinuousAssembly.hh>
#include <devel/sewing/hashing/Hasher.hh>
#include <devel/sewing/util/io.hh>
#include <devel/sewing/util/util.hh>
#include <devel/sewing/sampling/SewGraph.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <utility/tag/Tag.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <core/scoring/TwelveANeighborGraph.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/TaskOperationCreators.hh>
#include <core/pack/task/TaskFactory.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/features/SecondaryStructureSegmentFeatures.hh>
#include <protocols/features/SmotifFeatures.hh>

namespace devel {
namespace sewing {

static basic::Tracer TR( "devel.sewing.NodeConstraintAssemblyMover" );

////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Mover Creator Functions   ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

protocols::moves::MoverOP
NodeConstraintAssemblyMoverCreator::create_mover() const
{
	return new NodeConstraintAssemblyMover;
}

std::string
NodeConstraintAssemblyMoverCreator::keyname() const
{
	return NodeConstraintAssemblyMoverCreator::mover_name();
}

std::string
NodeConstraintAssemblyMoverCreator::mover_name()
{
	return "NodeConstraintAssemblyMover";
}

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Mover  Functions   ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

protocols::moves::MoverOP
NodeConstraintAssemblyMover::clone() const {
	return( protocols::moves::MoverOP( new NodeConstraintAssemblyMover( *this ) ) );
}
protocols::moves::MoverOP
NodeConstraintAssemblyMover::fresh_instance() const {
	return protocols::moves::MoverOP( new NodeConstraintAssemblyMover );
}

std::string
NodeConstraintAssemblyMover::get_name() const {
	return "NodeConstraintAssemblyMover";
}

NodeConstraintAssemblyMover::NodeConstraintAssemblyMover():
	parent()
{}


bool
NodeConstraintAssemblyMover::check_constraints(
	core::Size index,
	ModelNode const * node
) const {

	if(index > node_constraints_.size()) {
		utility_exit_with_message("Attempting to access an invalid constraint! You need a constraint for every node!");
	}
	NodeConstraint constraint = node_constraints_[index];
	Model const & model = node->model();

	if(
		model.segments_[1].dssp_ == constraint.ss1_ && model.segments_.back().dssp_ == constraint.ss2_ &&
		model.distance_ >= constraint.ideal_distance_-constraint.distance_range_ && model.distance_ <= constraint.ideal_distance_+constraint.distance_range_ &&
		model.hoist_angle_degrees_ >= constraint.ideal_hoist_-constraint.hoist_range_ && model.hoist_angle_degrees_ <= constraint.ideal_hoist_+constraint.hoist_range_ &&
		model.packing_angle_degrees_ >= constraint.ideal_packing_-constraint.packing_range_ && model.packing_angle_degrees_ <= constraint.ideal_packing_+constraint.packing_range_ &&
		model.meridian_angle_degrees_ >= constraint.ideal_meridian_-constraint.meridian_range_ && model.meridian_angle_degrees_ <= constraint.ideal_meridian_+constraint.meridian_range_
	) {
		return true;
	}
	return false;
}



///@brief check all constraints that will be evaluated
bool
NodeConstraintAssemblyMover::check_all_nodes() const {
	core::Size n_nodes = graph_->num_nodes();
	for(core::Size ci = 1; ci <= num_edges_to_follow_+1; ++ci) {
		core::Size counter = 0;
		for(core::Size i=1; i<=n_nodes; ++i) {
			ModelNode const * const node = graph_->get_model_node(i);
			if(check_constraints(ci, node)) {
				++counter;
			}
		}
		TR << "Total valid nodes for constraint " << ci << ": " << counter << std::endl;
		if(counter == 0) { return false; }
	}
	return true;
}



///@details get a starting node that conforms to the
///first node definition. No sense in starting with something
///that won't work
ModelNode const *
NodeConstraintAssemblyMover::get_starting_model_node() {

	//randomly permute nodes until you find one that fits with
	//the first node definition
	core::Size n_nodes = graph_->num_nodes();
	utility::vector1<core::Size> order;
	for(core::Size i=1; i<=n_nodes; ++i) {
		order.push_back(i);
	}
	numeric::random::random_permutation(order, numeric::random::rg());

	for(core::Size i=1; i<=n_nodes; ++i) {
		ModelNode const * const node = graph_->get_model_node(order[i]);
		if(check_constraints(1, node)) {
			graph_->add_edges_from_binary(edge_file_, node->get_node_index());
			return node;
		}
	}

	utility_exit_with_message("No model found that satisfies the first node constraint!");
	return 0;
}



///@details setup the node constraints and then apply the base
///class method
void
NodeConstraintAssemblyMover::apply(
	core::pose::Pose & pose
) {

	//Read in:file:native if it exits
	//calculate the smotif angles from these values

	using basic::options::option;
	using namespace basic::options::OptionKeys;

	core::pose::PoseOP template_pose;
	if(option[in::file::native].user()) {
		core::chemical::ResidueTypeSetCOP res_type_set =
			core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		template_pose = core::import_pose::pose_from_pdb( *res_type_set, option[ in::file::native ]() );
	}
	else {
		utility_exit_with_message("You must provide a template with the in:file:native option");
	}

	core::scoring::dssp::Dssp dssp(*template_pose);
	dssp.dssp_reduced();


	//Find segments of secondary structure
	char seg_dssp = 'N'; //invalid starting value
	utility::vector1< std::pair < std::pair<core::Size, core::Size>, char > > segments;
	std::pair<core::Size, core::Size> cur_segment;
	for(core::Size i=1; i<= template_pose->total_residue(); ++i) {
		if(template_pose->residue(i).is_protein()) {
			if(seg_dssp == 'N' || dssp.get_dssp_secstruct(i) != seg_dssp) {
				if(seg_dssp != 'N') {
					cur_segment.second = i-1;
					segments.push_back(std::make_pair(cur_segment, seg_dssp));
				}
				cur_segment.first = i;
				seg_dssp = dssp.get_dssp_secstruct(i);
			}
		}
	}
	cur_segment.second = template_pose->total_residue();
	segments.push_back(std::make_pair(cur_segment, seg_dssp));

	if(segments.size() < 2) {
		utility_exit_with_message("Template PDB has fewer than two secondary structure segments and is not suitable for mimicry!");
	}

	//If we are starting with a loop, fold it into the first segment
	if(segments[1].second == 'L') {
		segments[2].first.first = segments[1].first.first;
		segments.erase(segments.begin());
	}

	protocols::features::SmotifFeatures smotifs;
	for(core::Size i=1; i<=segments.size()-2; i+=2) {
		if(segments[i].second != 'L' && segments[i+1].second == 'L' && segments[i+2].second != 'L') {
			core::Size ss1_begin = segments[i].first.first;
			core::Size ss1_end = segments[i].first.second;
			utility::vector1< numeric::xyzVector< core::Real > > ss1_coords;
			for(core::Size j = ss1_begin; j <= ss1_end; ++j) {
				ss1_coords.push_back(template_pose->residue(j).atom("CA").xyz());
			}

			core::Size ss2_begin = segments[i+2].first.first;
			core::Size ss2_end = segments[i+2].first.second;
			utility::vector1< numeric::xyzVector< core::Real > > ss2_coords;
			for(core::Size j = ss2_begin; j <= ss2_end; ++j) {
				ss2_coords.push_back(template_pose->residue(j).atom("CA").xyz());
			}

			core::Real distance, hoist, packing, meridian;
			smotifs.calculate_angles(ss1_coords, ss2_coords, distance, hoist, packing, meridian);

			NodeConstraint constraint;
			constraint.ss1_ = segments[i].second;
			constraint.ss2_ = segments[i+2].second;
			constraint.ideal_distance_ = distance;
			constraint.ideal_hoist_ = hoist;
			constraint.ideal_packing_ = packing;
			constraint.ideal_meridian_ = meridian;

			if(TR.Debug.visible()) {
				TR.Debug << "Constraining smotif starting with segment " << i << ": " << std::endl;
				TR.Debug << "\tSS1 DSSP: " << segments[i].second << std::endl;
				TR.Debug << "\tSS2 DSSP: " << segments[i+2].second << std::endl;
				TR.Debug << "\tdistance: " << distance << std::endl;
				TR.Debug << "\thoist: " << hoist << std::endl;
				TR.Debug << "\tpacking: " << packing << std::endl;
				TR.Debug << "\tmeridian: " << meridian << std::endl;
			}

			node_constraints_.push_back(constraint);
		}
		else {
			TR.Warning << "Template protein doesn't follow canonical smotif patterning" << std::endl;
		}
	}

	if(num_edges_to_follow_+1 > node_constraints_.size()) {
		utility_exit_with_message("Requested " + utility::to_string(num_edges_to_follow_) + " nodes, but only have " + utility::to_string(node_constraints_.size()) + " constraints");
	}
	TR << "Added " << node_constraints_.size() << " node constraints" << std::endl;

	if(!check_all_nodes()) {
		utility_exit_with_message("Impossible to satisfy all node constraints");
	}

	//call base-class apply
	parent::apply(pose);
}


bool
NodeConstraintAssemblyMover::check_edge(
	AssemblyOP assembly,
	ModelNode const * const reference_node,
	HashEdge const * const new_edge,
	bool allow_repeat
){
	ModelNode const * const new_model =
		graph_->get_model_node(new_edge->get_other_ind(reference_node->get_node_index()));
	core::Size added_edges = assembly->path().size();

	if(!check_constraints(added_edges+2, new_model)) {
		return false;
	}

	return parent::check_edge(assembly, reference_node, new_edge, allow_repeat);
}

} //sewing
} //devel

#endif
