// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file score_interface.cc
/// @brief simple pilot app to score an interface
/// @author Noah Ollikainen

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <core/id/types.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/methods/Methods.hh> //for long range energies
#include <core/scoring/LREnergyContainer.hh> //long range energies
#include <core/scoring/constraints/util.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

basic::Tracer TR("apps.score_interface");

void register_metrics() {

	core::pose::metrics::PoseMetricCalculatorOP interface_neighbor_calculator = new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator((core::Size)1,(core::Size)2);
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "interface_neighbor", interface_neighbor_calculator );
	
}

core::Real compute_interface_energy(
	core::pose::Pose & this_pose) {
	
	basic::MetricValue<core::Size> mv_size;
	this_pose.metric("interface_neighbor","first_chain_first_resnum",mv_size);
	core::Size ch1_begin_num = mv_size.value();
	this_pose.metric("interface_neighbor","first_chain_last_resnum",mv_size);
	core::Size ch1_end_num = mv_size.value();
	this_pose.metric("interface_neighbor","second_chain_first_resnum",mv_size);
	core::Size ch2_begin_num = mv_size.value();
	this_pose.metric("interface_neighbor","second_chain_last_resnum",mv_size);
	core::Size ch2_end_num = mv_size.value();

	// Clear the energy-holders, get the (unweighted) energies from the pose
	core::scoring::EnergyMap delta_energies_unweighted_;
	delta_energies_unweighted_.clear();
	
	core::scoring::EnergyGraph const & energy_graph( this_pose.energies().energy_graph() );
	
	std::vector<core::Size> interface_positions;
	
	// Loop over interactions across the interface
	for ( core::Size i = ch1_begin_num; i <= ch1_end_num; ++i ) {
		for ( core::graph::Graph::EdgeListConstIter
						iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
						irue = energy_graph.get_node(i)->const_upper_edge_list_end();
					iru != irue; ++iru ) {
			const core::scoring::EnergyEdge * edge( static_cast< const core::scoring::EnergyEdge *> (*iru) );
			core::Size const j( edge->get_second_node_ind() );
						
			if ( ( j >= ch2_begin_num ) && ( j <= ch2_end_num ) ) {
				delta_energies_unweighted_ += edge->fill_energy_map();
			}
		}
	}

	// Graph is asymmetric, so switch i/j and redo
	for ( core::Size i = ch2_begin_num; i <= ch2_end_num; ++i ) {
		for ( core::graph::Graph::EdgeListConstIter
						iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
						irue = energy_graph.get_node(i)->const_upper_edge_list_end();
					iru != irue; ++iru ) {
			const core::scoring::EnergyEdge * edge( static_cast< const core::scoring::EnergyEdge *> (*iru) );
			core::Size const j( edge->get_second_node_ind() );
						
			if ( ( j >= ch1_begin_num ) && ( j <= ch1_end_num ) ) {
				delta_energies_unweighted_ += edge->fill_energy_map();
			}
		}
	}
	
	core::scoring::EnergyMap weights_ = this_pose.energies().weights();
	core::Real weighted_total_ = delta_energies_unweighted_.dot(weights_);
	return weighted_total_;
	
}

int main( int argc, char * argv [] )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::task;
		
	// initialize Rosetta
	devel::init(argc, argv);
	
	register_metrics();
	
	core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();
	
	utility::vector1<std::string> pdbs = basic::options::option[ in::file::s ]();
	
	core::Real min_interface_score = std::numeric_limits<core::Real>::infinity();
	core::Real min_total_score = std::numeric_limits<core::Real>::infinity();
	std::string min_interface_filename;
	std::string min_total_filename;
	
	for(core::Size i = 1; i <= pdbs.size(); i++) {
	
		core::pose::Pose pose;
		core::import_pose::pose_from_pdb(pose, pdbs[i]);
		(*score_fxn)(pose);

		core::Real interface_score = compute_interface_energy(pose);
		if (min_interface_score > interface_score) {
			min_interface_score = interface_score;
			min_interface_filename = pdbs[i];
		}
		
		core::Real total_score = pose.energies().total_energy();
		if (min_total_score > total_score) {
			min_total_score = total_score;
			min_total_filename = pdbs[i];
		}
		TR << pdbs[i] << "\t" << "\t" << total_score << "\t" << interface_score << std::endl;

	}
	
	TR << "MINIMUM TOTAL" << "\t" << min_total_filename << "\t" << min_total_score << std::endl;
	TR << "MINIMUM INTERFACE" << "\t" << min_interface_filename << "\t" << min_interface_score << std::endl;
	
	return 0;
}
