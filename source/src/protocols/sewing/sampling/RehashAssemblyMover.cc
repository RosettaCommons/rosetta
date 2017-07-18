// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RehashAssemblyMover.cc
///
/// @brief Derived from AssemblyMover, this mover generates an assembly and then uses Loophash segments
/// to connect any broken segments resulting from the assembly process
/// @author Tim Jacobs

// Unit Headers
#include <devel/sewing/sampling/RehashAssemblyMover.hh>
#include <devel/sewing/sampling/RehashAssemblyMoverCreator.hh>

//Package headers
#include <devel/sewing/util/io.hh>

//Protocol headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>

#include <core/pose/util.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/util.hh>

#include <protocols/analysis/LoopAnalyzerMover.hh>

#include <protocols/simple_moves/MinMover.hh>

#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <protocols/sic_dock/loophash_util.hh>
#include <protocols/relax/AtomCoordinateCstMover.hh>

//Utility headers
#include <basic/Tracer.hh>

#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#include <devel/denovo_design/ConnectJumps.hh>

namespace devel {
namespace sewing {

static basic::Tracer TR( "devel.sewing.RehashAssemblyMover" );

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Boiler Plate Code   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
protocols::moves::MoverOP
RehashAssemblyMoverCreator::create_mover() const
{
	return new RehashAssemblyMover;
}

std::string
RehashAssemblyMoverCreator::keyname() const
{
	return RehashAssemblyMoverCreator::mover_name();
}

std::string
RehashAssemblyMoverCreator::mover_name()
{
	return "RehashAssemblyMover";
}

protocols::moves::MoverOP
RehashAssemblyMover::clone() const {
	return( protocols::moves::MoverOP( new RehashAssemblyMover( *this ) ) );
}
protocols::moves::MoverOP
RehashAssemblyMover::fresh_instance() const {
	return protocols::moves::MoverOP( new RehashAssemblyMover );
}

std::string
RehashAssemblyMover::get_name() const {
	return "RehashAssemblyMover";
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  RehashAssemblyMover functions   //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

RehashAssemblyMover::RehashAssemblyMover():
	AssemblyMover(),
	max_loop_distance_(14)
{}

bool
RehashAssemblyMover::generate_assembly(
	AssemblyOP & assembly
){
	using namespace basic::options;

	if ( !AssemblyMover::generate_assembly(assembly) ) {
		return false;
	}

	bool succesfully_rearranged = rearrange_assembly(assembly);
	if ( !succesfully_rearranged ) { return false; }

	utility::vector1<core::Size> disconnected_segments = assembly->disconnected_segments();

	Model test_model = assembly->add_bridge_model(disconnected_segments[1], disconnected_segments[1]+1);

	Hasher hasher;
	for ( std::map< int, Model >::const_iterator it = bridge_models_.begin(); it != bridge_models_.end(); ++it ) {
		hasher.insert(it->second);
	}
	core::Size num_segments_to_match = 2;
	core::Size min_hash_score = option[OptionKeys::sewing::min_hash_score].value();
	core::Size max_clash_score = option[OptionKeys::sewing::max_clash_score].value();

	TR << "Begin scoring" << std::endl;
	core::Size starttime = time(NULL);
	ScoreResults scores = hasher.score(test_model, num_segments_to_match, min_hash_score, max_clash_score, true);
	core::Size endtime = time(NULL);
	TR << "Scoring complete in time " << endtime-starttime << " seconds, found " << scores.size() << " unique alignments." << std::endl;
	if ( scores.size() == 0 ) {
		return false;
	}
	ScoreResults::const_iterator it = scores.begin();
	ScoreResults::const_iterator it_end = scores.end();
	core::Size counter=0;
	for ( ; it != it_end; ++it ) {
		++counter;
		TR << "Alignment " << counter << std::endl;
		TR << "basis 1 " << it->first.first.model_id << " " << it->first.first.resnum << std::endl;
		TR << "basis 2 " << it->first.second.model_id << " " << it->first.second.resnum << std::endl;
	}
	bridge_models_.insert(std::make_pair(0, test_model));

	SewGraphOP bridge_graph = new SewGraph(bridge_models_, 1);
	//, scores_to_alignments(scores));
	// ModelNode const * reference_node = bridge_graph->get_model_node_from_model_id(0);
	ModelNode const * reference_node = 0;
	utility::graph::EdgeListConstIterator edge_it = reference_node->const_edge_list_begin();
	//for(core::Size j=0; j<edge_order[cur_edge_ind]; ++j) {
	// ++edge_it;
	//}
	HashEdge const * const cur_edge = static_cast< HashEdge const * >(*edge_it);

	//Check the current edge to see if it can possible be added to the assembly
	assembly->to_multichain_pose(core::chemical::CENTROID).dump_pdb("before_bridge.pdb");
	if ( assembly->check_edge(bridge_graph, cur_edge, reference_node) ) {
		assembly->follow_edge(bridge_graph, cur_edge, reference_node->get_node_index());
		assembly->to_multichain_pose(core::chemical::CENTROID).dump_pdb("after_bridge.pdb");
	}

	return true;
}

bool
RehashAssemblyMover::rearrange_assembly(
	AssemblyOP & assembly
) const {

	utility::vector1< utility::vector1<core::Size> > orders = assembly->find_possible_orders(max_loop_distance_);
	numeric::random::random_permutation(orders, numeric::random::rg());
	TR << "Permuting " << orders.size() << " valid assembly orders" << std::endl;

	AssemblyOP reordered_assembly = assembly->clone();
	for ( core::Size i=1; i<=orders.size(); ++i ) {
		reordered_assembly = assembly->clone();
		reordered_assembly->reorder(orders[i]);
		core::pose::Pose assembly_pose = reordered_assembly->to_pose(core::chemical::CENTROID);
		//core::Size num_lh_frags = count_loophash_fragments(reordered_assembly, assembly_pose);
		//if((int)num_lh_frags >= basic::options::option[ basic::options::OptionKeys::sewing::min_lh_fragments].value()) {
		assembly = reordered_assembly;
		return true;
		//}
	}
	TR << "Assembly could not be rearranged in a connectable way." << std::endl;
	return false;
}

void
RehashAssemblyMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
){
	using namespace basic::options;

	AssemblyMover::parse_my_tag(tag, data, filters, movers, pose);

	if ( tag->hasOption("max_loop_distance") ) {
		max_loop_distance_ = tag->getOption<core::Real>("max_loop_distance");
	}

	if ( tag->hasOption("bridge_model_file") ) {
		std::string model_file = tag->getOption<std::string>("bridge_model_file");
		bridge_models_ = read_model_file(model_file);
	}
}

} //sewing
} //devel
