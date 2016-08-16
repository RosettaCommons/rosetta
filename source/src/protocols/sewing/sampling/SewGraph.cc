// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SewGraph.cc
///
/// @brief
/// @author   Tim Jacobs

//Unit headers
#include <protocols/sewing/sampling/SewGraph.hh>

//Package headers
#include <protocols/sewing/util/io.hh>

//Protocol headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.functions.hh>

//Utility headers
#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>

#include <numeric/random/random.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/LexicographicalIterator.hh>

#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

// Boost Headers
#include <boost/cstdint.hpp>
#include <boost/pool/pool.hpp>
#include <core/graph/unordered_object_pool.hpp>

namespace protocols {
namespace sewing  {

static basic::Tracer TR("protocols.sewing.SewGraph");

//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////   ModelNode Class   //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

ModelNode::ModelNode( core::graph::Graph * owner, core::Size index ) :
	core::graph::Node( owner, index )
{}

Model const &
ModelNode::model() const{
	return model_;
}

void
ModelNode::model(
	Model const & model
){
	model_ = model;
}

std::set<core::Size> const &
ModelNode::segment_ids() const {
	return segment_ids_;
}

void
ModelNode::segment_ids(
	std::set<core::Size> const & segment_ids
){
	segment_ids_ = segment_ids;
}

void
ModelNode::copy_from(
	core::graph::Node const * source
){
	ModelNode const * mn_source = static_cast< ModelNode const * > ( source );
	model_ = mn_source->model_;
	segment_ids_ = mn_source->segment_ids_;
}

std::ostream &
operator<<(std::ostream& out, ModelNode const & node ) {
	out << "Node " << node.get_node_index() << ": model_id ("
		<< node.model().model_id_ << ") - segments (";
	std::set< core::Size > const & segment_ids = node.segment_ids();
	std::set< core::Size >::const_iterator it = segment_ids.begin();
	std::set< core::Size >::const_iterator it_end = segment_ids.end();
	for ( ; it != it_end; ++it ) {
		out << *it;
		if ( it != --segment_ids.end() ) {
			out << ", ";
		}
	}
	out << ")";
	return out;
}

//////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  HashEdge Class   //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

HashEdge::HashEdge(
	SewGraph * owner,
	core::Size n1,
	core::Size n2
) :
	core::graph::Edge( owner, n1, n2 )
{}

core::Size
HashEdge::model_resnum(
	int model_id
) const {
	if ( model_id == basis_pair_.first.model_id ) {
		return basis_pair_.first.resnum;
	} else if ( model_id == basis_pair_.second.model_id ) {
		return basis_pair_.second.resnum;
	} else {
		utility_exit_with_message("No model with id " +utility::to_string(model_id)+ " in BasisPair");
		return 0;
	}
}

std::ostream &
operator<<(std::ostream& out, HashEdge const & edge ) {
	out << "HashEdge: " <<
		"\t" << edge.basis_pair_.first.model_id << " " << edge.basis_pair_.first.resnum <<
		"\t" << edge.basis_pair_.second.model_id << " " << edge.basis_pair_.second.resnum << std::endl;
	return out;
}

//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////   SewGraph Class   ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

SewGraph::~SewGraph() {
	delete_everything();
	delete hash_edge_pool_;
	hash_edge_pool_ = 0;
}

/// @details Notice that this does not call the parent( src ) copy constructor.
/// This is because the copy constructor relies on polymorphic functions which
/// are unavailable during the Graph constructor.  Instead, this function waits
/// until parent construction is complete, and relies on the assigmnent operator.
SewGraph::SewGraph( SewGraph const & src ):
	core::graph::Graph(),
	hash_edge_pool_( new boost::unordered_object_pool< HashEdge > ( 256 ) ),
	model_indices_(src.model_indices_),
	last_node_added_(src.last_node_added_),
	special_edge_data_(src.special_edge_data_)
{
	core::graph::Graph::operator = ( src );
}

SewGraph::SewGraph():
	core::graph::Graph(),
	hash_edge_pool_( new boost::unordered_object_pool< HashEdge > ( 256 ) )
{}

///@details iterate through the models and add all the nodes to the graph
SewGraph::SewGraph(
	std::map< int, Model > const & models,
	core::Size segment_matches_per_edge
):
	core::graph::Graph(),
	hash_edge_pool_( new boost::unordered_object_pool< HashEdge > ( 256 ) )
{

	core::Size starttime = time(NULL);

	//first, figure out how many nodes we have (this is silly to do twice)
	core::Size n_nodes = 0;
	std::map< int, Model >::const_iterator models_it = models.begin();
	std::map< int, Model >::const_iterator models_end = models.end();
	for ( ; models_it != models_end; ++models_it ) {
		utility::vector1<core::Size> dim_sizes(segment_matches_per_edge, models_it->second.segments_.size());
		for ( utility::LexicographicalIterator lex( dim_sizes ); ! lex.at_end(); ++lex ) {
			bool invalid = false;
			for ( core::Size i=1; i<=segment_matches_per_edge; ++i ) {
				for ( core::Size j=i+1; j<=segment_matches_per_edge; ++j ) {
					if ( lex[j] >= lex[i] ) {
						invalid = true;
					}
				}
			}
			if ( invalid ) {
				continue;
			}
			++n_nodes;
		}
	}

	//Next, create the graph and map a model index to the index used by the graph
	set_num_nodes(n_nodes);//this clears the graph
	model_indices_.clear();

	models_it = models.begin();
	core::Size counter = 1;
	for ( ; models_it != models_end; ++models_it ) {
		if ( models_it->first <= 0 ) {
			continue;
		}
		if ( model_indices_.find(models_it->first) != model_indices_.end() ) {
			std::stringstream err;
			err << "Duplicate model found with ID: " << models_it->first << ", Exiting!" << std::endl;
			utility_exit_with_message(err.str());
		}

		//Iterate all combinations of 'num_matches_per_edge' segments. Each of these will be a node in the graph
		utility::vector1<core::Size> dim_sizes(segment_matches_per_edge, models_it->second.segments_.size());
		for ( utility::LexicographicalIterator lex( dim_sizes ); ! lex.at_end(); ++lex ) {

			//Don't allow multiple combinations of the same combinations of segments
			bool invalid = false;
			for ( core::Size i=1; i<=segment_matches_per_edge; ++i ) {
				for ( core::Size j=i+1; j<=segment_matches_per_edge; ++j ) {
					if ( lex[j] >= lex[i] ) {
						invalid = true;
					}
				}
			}
			if ( invalid ) {
				continue;
			}

			model_indices_[models_it->first].insert(counter);
			ModelNode * node = static_cast< ModelNode * >( get_node( counter ) );
			node->model(models_it->second);

			std::set<core::Size> segments;
			for ( core::Size i=1; i<=segment_matches_per_edge; ++i ) {
				segments.insert(models_it->second.segments_[lex[i]].segment_id_);
			}
			node->segment_ids(segments);
			++counter;
		}
	}

	//Add negative (special) models at the end so as not to mess with indexing of the binary
	//edge files
	models_it = models.begin();
	for ( ; models_it != models_end; ++models_it ) {
		if ( models_it->first > 0 ) {
			break;
		}

		if ( model_indices_.find(models_it->first) != model_indices_.end() ) {
			std::stringstream err;
			err << "Duplicate model found with ID: " << models_it->first << ", Exiting!" << std::endl;
			utility_exit_with_message(err.str());
		}

		//Iterate all combinations of 'num_matches_per_edge' segments. Each of these will be a node in the graph
		utility::vector1<core::Size> dim_sizes(segment_matches_per_edge, models_it->second.segments_.size());
		for ( utility::LexicographicalIterator lex( dim_sizes ); ! lex.at_end(); ++lex ) {

			//Don't allow multiple combinations of the same combinations of segments
			bool invalid = false;
			for ( core::Size i=1; i<=segment_matches_per_edge; ++i ) {
				for ( core::Size j=i+1; j<=segment_matches_per_edge; ++j ) {
					if ( lex[j] >= lex[i] ) {
						invalid = true;
					}
				}
			}
			if ( invalid ) {
				continue;
			}

			model_indices_[models_it->first].insert(counter);
			ModelNode * node = static_cast< ModelNode * >( get_node( counter ) );
			node->model(models_it->second);

			std::set<core::Size> segments;
			for ( core::Size i=1; i<=segment_matches_per_edge; ++i ) {
				segments.insert(models_it->second.segments_[lex[i]].segment_id_);
			}
			node->segment_ids(segments);
			++counter;
		}
	}

	core::Size endtime = time(NULL);
	TR << "SewGraph initialized with " << num_nodes() << " nodes. (" << endtime-starttime << " seconds)" << std::endl;
}

void
SewGraph::set_special_edges(
	ScoreResults const & scores
){
	special_edge_data_ = scores;
}

void
SewGraph::add_special_edges(){
	ScoreResults::const_iterator it = special_edge_data_.begin();
	ScoreResults::const_iterator end = special_edge_data_.end();
	for ( ; it != end; ++it ) {

		int model_id_1 = it->first.first.model_id;
		int model_id_2 = it->first.second.model_id;

		if ( model_id_1 == model_id_2 ) {
			continue;
		}

		core::Size resnum_1 = it->first.first.resnum;
		core::Size resnum_2 = it->first.second.resnum;

		std::set<core::Size> model_1_segs;
		std::set<core::Size> model_2_segs;

		std::map< SegmentPair, core::Size > const & segment_pairs = it->second.segment_match_counts;
		std::map< SegmentPair, core::Size >::const_iterator seg_it = segment_pairs.begin();
		std::map< SegmentPair, core::Size >::const_iterator seg_it_end = segment_pairs.end();
		for ( ; seg_it != seg_it_end; ++seg_it ) {
			model_1_segs.insert(seg_it->first.first);
			model_2_segs.insert(seg_it->first.second);
		}

		ModelNode const * const n1 = get_model_node(model_id_1, model_1_segs);
		ModelNode const * const n2 = get_model_node(model_id_2, model_2_segs);

		//For continuous models, remove edges that aren't between first elements and last elements
		if ( basic::options::option[basic::options::OptionKeys::sewing::assembly_type].value() == "continuous" ) {
			if ( !
					((resnum_1 >= n1->model().segments_[1].residues_[1].resnum_ && resnum_1 <= n1->model().segments_[1].residues_.back().resnum_ &&
					resnum_2 >= n2->model().segments_.back().residues_[1].resnum_ && resnum_2 <= n2->model().segments_.back().residues_.back().resnum_) ||
					(resnum_1 >= n1->model().segments_.back().residues_[1].resnum_ && resnum_1 <= n1->model().segments_.back().residues_.back().resnum_ &&
					resnum_2 >= n2->model().segments_[1].residues_[1].resnum_ && resnum_2 <= n2->model().segments_[1].residues_.back().resnum_))
					) {
				continue;
			}
		}

		BasisPair bp;
		bp.first.model_id = model_id_1;
		bp.first.resnum = resnum_1;
		bp.second.model_id = model_id_2;
		bp.second.resnum = resnum_2;
		HashEdge * const e = static_cast< HashEdge * >( add_edge(n1->get_node_index(), n2->get_node_index()));
		e->basis_pair(bp);
	}
}//add_special_edges


/// @brief Factory method for node creation
core::graph::Node*
SewGraph::create_new_node( Size index )
{
	return new ModelNode( this, index );
}

/// @brief Factory method for edge creation
core::graph::Edge*
SewGraph::create_new_edge( Size index1, Size index2 )
{
	return hash_edge_pool_->construct( this, index1, index2 );
}

void SewGraph::delete_edge( core::graph::Edge * edge )
{
	assert( dynamic_cast< HashEdge* > (edge) );
	hash_edge_pool_->destroy( static_cast< HashEdge* > (edge) );
}


HashEdge *
SewGraph::find_hash_edge(
	Size n1,
	Size n2
){
	core::graph::Edge * edge( find_edge( n1, n2 ) );
	if ( edge ) {
		return static_cast< HashEdge * > ( edge );
	} else {
		return 0;
	}
}

HashEdge const *
SewGraph::find_hash_edge(
	Size n1,
	Size n2
) const {
	core::graph::Edge const * edge( find_edge( n1, n2 ) );
	if ( edge ) {
		return static_cast< HashEdge const * > ( edge );
	} else {
		return 0;
	}
}

// possible example of usage, SewGraph.cc:line 791, ModelNode const * const node = get_model_node(node_id);
ModelNode const *
SewGraph::get_model_node(
	Size n
) const {
	core::graph::Node const * node( get_node(n) );
	if ( node ) {
		return static_cast< ModelNode const * > ( node );
	} else {
		return 0;
	}
}

ModelNode const *
SewGraph::get_model_node(
	int model_id,
	std::set<core::Size> segment_ids
) const {

	std::map< int, std::set<core::Size> >::const_iterator it =
		model_indices_.find(model_id);
	if ( it == model_indices_.end() ) {
		utility_exit_with_message("No nodes found for model id " + utility::to_string(model_id));
	}
	std::set<core::Size>::const_iterator node_it = it->second.begin();
	std::set<core::Size>::const_iterator node_it_end = it->second.end();
	for ( ; node_it != node_it_end; ++node_it ) {
		ModelNode const * node( static_cast< ModelNode const * >(get_node(*node_it)) );
		if ( node->segment_ids() == segment_ids ) {
			return node;
		}
	}
	utility_exit_with_message("No node found with given model id and segment ids");
	return 0;

}

ModelNode const *
SewGraph::get_random_node() const {
	core::Size rand_index = numeric::random::random_range(1, (int)num_nodes());
	return get_model_node(rand_index);
}

ModelNode const *
SewGraph::get_random_node_with_edges() const {
	//No infinite loops
	runtime_assert(num_edges() > 0);

	while ( true ) {
		ModelNode const * node = get_random_node();
		if ( node->num_edges() > 0 ) { return node; }
	}
}



void
SewGraph::generate_binary_score_file(
	std::string score_filename,
	std::string binary_filename
){
	utility::io::izstream file;
	file.open(score_filename);

	if ( !file.good() ) {
		utility_exit_with_message("Could not find Hasher score file with name: " + score_filename);
	}

	//model id/segment set -> vector of other pairs of other model id/segment sets, pair of resnums
	typedef std::pair<int, std::set<core::Size> > score_node;
	typedef std::map< score_node, utility::vector1< std::pair< score_node, std::pair<core::Size,core::Size> > > > EdgeList;

	EdgeList edges;

	boost::uint32_t n_segs_per_edge = 0;
	std::string line;
	while ( getline( file, line ) ) {
		utility::vector1<std::string> tokens = utility::string_split(line);
		assert(tokens.size() >= 7);

		if ( n_segs_per_edge == 0 ) {
			n_segs_per_edge = (boost::int32_t)((tokens.size() - 5)/2);
		} else {
			assert(n_segs_per_edge == (tokens.size() - 5)/2);
		}

		score_node node1;
		score_node node2;

		core::Size token = 1;
		node1.first = utility::string2int(tokens[token++]);
		core::Size resnum_1 = utility::string2int(tokens[token++]);
		for ( core::Size i=1; i<=n_segs_per_edge; ++i ) {
			node1.second.insert(utility::string2int(tokens[token++]));
		}

		node2.first = utility::string2int(tokens[token++]);
		core::Size resnum_2 = utility::string2int(tokens[token++]);
		for ( core::Size i=1; i<=n_segs_per_edge; ++i ) {
			node2.second.insert(utility::string2int(tokens[token++]));
		}

		edges[node1].push_back(std::make_pair(node2, std::make_pair(resnum_1, resnum_2)));
		edges[node2].push_back(std::make_pair(node1, std::make_pair(resnum_2, resnum_1)));
	}
	file.close();
	TR << "Done reading plain-text score file - " << edges.size() << " nodes have edges." << std::endl;

	utility::io::ozstream binary_file;
	binary_file.open(binary_filename, std::ios::out | std::ios::binary);
	if ( !binary_file ) {
		utility_exit_with_message("Couldn't open " + binary_filename + " for writing!");
	}

	boost::uint32_t n_nodes_with_edges = (boost::uint32_t)edges.size();

	// if(num_models > num_nodes()) {
	//  utility_exit_with_message("Error converting scores to binary! There are more nodes in the score file than contained in the graph!");
	// }

	if ( n_nodes_with_edges < num_nodes() ) {
		TR.Warning << "There are " << num_nodes() << " model nodes in your model file, but only " << n_nodes_with_edges << " nodes have edges" << std::endl;

		//We need elements for all nodes, even if there are no edges, so make sure we have them here
		for ( core::Size i=1; i<=num_nodes(); ++i ) {
			score_node query_node;
			query_node.first= get_model_node(i)->model().model_id_;
			query_node.second = get_model_node(i)->segment_ids();
			if ( edges.find(query_node) == edges.end() ) {
				edges[query_node].clear();
			}
		}
	}
	n_nodes_with_edges = (boost::uint32_t)edges.size();

	if ( n_nodes_with_edges != num_nodes() ) {
		TR << "n_nodes_with_edges: " << n_nodes_with_edges << std::endl;
		TR << "num_nodes(): " << num_nodes() << std::endl;
	}
	runtime_assert(n_nodes_with_edges == num_nodes());

	//First write the number of models for validation later
	binary_file.write( (char*) & n_nodes_with_edges, sizeof( boost::uint32_t ));

	//Next write the number of segments per edge
	binary_file.write( (char*) & n_segs_per_edge, sizeof( boost::uint32_t ));

	//Next write the number edges for each score node
	EdgeList::const_iterator it = edges.begin();
	EdgeList::const_iterator it_end = edges.end();
	for ( ; it != it_end; ++it ) {
		boost::uint32_t num = (boost::uint32_t)it->second.size();
		if ( TR.Debug.visible() ) {
			TR << "model " << it->first.first << " segments ";
			std::set<core::Size>::const_iterator seg_it = it->first.second.begin();
			std::set<core::Size>::const_iterator seg_it_end = it->first.second.end();
			for ( ; seg_it != seg_it_end; ++seg_it ) {
				TR << *seg_it << " ";
			}
			TR << "has " << num << " edges" << std::endl;
		}
		binary_file.write( (char*) & num , sizeof( boost::uint32_t ));
	}

	//Finally, write all the edges
	it = edges.begin();
	for ( ; it != it_end; ++it ) {
		utility::vector1< std::pair< score_node, std::pair<core::Size,core::Size> > >::const_iterator edge_it = it->second.begin();
		utility::vector1< std::pair< score_node, std::pair<core::Size,core::Size> > >::const_iterator edge_it_end = it->second.end();
		for ( ; edge_it != edge_it_end; ++edge_it ) {

			std::set<core::Size>::const_iterator seg_it = edge_it->first.second.begin();
			std::set<core::Size>::const_iterator seg_it_end = edge_it->first.second.end();
			for ( ; seg_it != seg_it_end; ++seg_it ) {
				boost::int32_t seg_id = (boost::int32_t)*seg_it;
				binary_file.write( (char*) & seg_id , sizeof( boost::int32_t ));
			}

			boost::int32_t other_model_id = edge_it->first.first;
			boost::int32_t resnum_1 = (boost::int32_t)edge_it->second.first;
			boost::int32_t resnum_2 = (boost::int32_t)edge_it->second.second;
			binary_file.write( (char*) & other_model_id , sizeof( boost::int32_t ));
			binary_file.write( (char*) & resnum_1 , sizeof( boost::int32_t ));
			binary_file.write( (char*) & resnum_2 , sizeof( boost::int32_t ));
		}
	}

	binary_file.close();
}

void
SewGraph::report_binary_stats(
	std::map< int, Model > const & models,
	std::string filename
){
	std::ifstream binary_file;
	std::string binary_filename = filename + ".bin";
	binary_file.open(binary_filename.c_str(), std::ios::in | std::ios::binary);
	if ( !binary_file ) {
		utility_exit_with_message("Couldn't open " + utility::to_string(binary_filename) + " for reading!");
	}

	//First, read the number of nodes for this binary file
	boost::uint32_t n_nodes;
	binary_file.read( (char*) &n_nodes, sizeof( boost::uint32_t )  );
	TR << "There are " << models.size() << " models in the model file and " << n_nodes << " nodes for the binary file" << std::endl;

	boost::uint32_t n_segments_per_node;
	binary_file.read( (char*) &n_segments_per_node, sizeof( boost::uint32_t )  );
	TR << "There are " << n_segments_per_node << " matching segments for each edge" << std::endl;

	boost::uint32_t n_edges;
	std::map< int, Model >::const_iterator it = models.begin();
	std::map< int, Model >::const_iterator it_end = models.end();
	for ( ; it != it_end; ++it ) {
		std::set<core::Size> const & nodes = model_indices_[it->first];
		std::set<core::Size>::const_iterator n_it = nodes.begin();
		std::set<core::Size>::const_iterator n_it_end = nodes.end();
		for ( ; n_it != n_it_end; ++n_it ) {
			binary_file.read( (char*) &n_edges, sizeof( boost::uint32_t )  );
			ModelNode const * const node = get_model_node(*n_it);
			TR << "Model " << it->first << ", Segments ";
			std::set<core::Size> const & segs = node->segment_ids();
			std::set<core::Size>::const_iterator seg_it = segs.begin();
			std::set<core::Size>::const_iterator seg_it_end = segs.end();
			for ( ; seg_it != seg_it_end; ++seg_it ) {
				TR << *seg_it << " ";
			}
			TR << "(node " << *n_it << ") has " << n_edges << " edges" << std::endl;
		}
	}

	//core::Size node_id = num_nodes();
	core::Size node_id = 3;
	while ( num_edges() == 0 ) {
		add_edges_from_binary(filename, node_id);
		ModelNode const * const node = get_model_node(node_id);
		TR << "Just added edges for model " << node->model().model_id_ << ", segments: ";
		std::set<core::Size> const & segments = node->segment_ids();
		std::set<core::Size>::const_iterator it = segments.begin();
		std::set<core::Size>::const_iterator it_end = segments.end();
		for ( ; it != it_end; ++it ) {
			TR << *it << std::endl;
		}
		++node_id;
	}
}



std::set<core::Size>
SewGraph::get_node_indices_from_model_id(
	int model_id
) const{

	std::map< int, std::set<core::Size> >::const_iterator find_it = model_indices_.find(model_id);
	if ( find_it == model_indices_.end() ) {
		utility_exit_with_message("Error: no model with model id " + utility::to_string(model_id) + " found in SewGraph");
	}
	return find_it->second;
}

void
SewGraph::add_all_model_edges_from_binary(
	std::string filename,
	int model_id
){
	std::set<core::Size> const & nodes = model_indices_[model_id];
	std::set<core::Size>::const_iterator it = nodes.begin();
	std::set<core::Size>::const_iterator it_end = nodes.end();
	for ( ; it != it_end; ++it ) {
		add_edges_from_binary(filename, *it);
	}
}

void
SewGraph::add_edges_from_binary(
	std::string filename,
	core::Size node_id
){

	//Don't add nodes that we just added
	if ( last_node_added_ == node_id ) { return; }
	last_node_added_ = node_id;

	//Drop all edges before every addition. This is MUCH faster than checking to see
	//which edges exist already and only create the new ones
	drop_all_edges();

	//Don't look for edges between 'special' models in the binary file
	if ( get_model_node(node_id)->model().model_id_ < 0 ) {
		add_special_edges();
		return;
	}

	std::ifstream binary_file;
	std::string binary_filename = filename + ".bin";
	binary_file.open(binary_filename.c_str(), std::ios::in | std::ios::binary);
	if ( !binary_file ) {
		utility_exit_with_message("Couldn't open " + utility::to_string(binary_filename) + " for writing!");
	}

	//First, read the number of nodes for this binary file
	boost::uint32_t n_nodes;
	binary_file.read( (char*) &n_nodes, sizeof( boost::uint32_t )  );
	//if(TR.Debug.visible()) {
	// TR << "There are " << n_nodes << " nodes in this file " << std::endl;
	//}
	if ( num_nodes() > n_nodes ) {
		TR.Warning << "There are more nodes in the graph than in the edge file. This should only happen if you have added PDB nodes!" << std::endl;
	}
	if ( TR.Debug.visible() ) {
		TR << "n_nodes: " << n_nodes << std::endl;
		TR << "num_nodes(): " << num_nodes() << std::endl;
	}
	runtime_assert(n_nodes <= num_nodes());

	//Next, read the number of segments per node
	boost::uint32_t n_segments_per_node;
	binary_file.read( (char*) &n_segments_per_node, sizeof( boost::uint32_t )  );
	//TR << "There are " << n_segments_per_node << " matching segments for each edge" << std::endl;

	//Advance past the number of nodes, n_segments/node and the counts for previous nodes
	long long advance = sizeof(boost::uint32_t) + sizeof(boost::uint32_t);
	advance += (node_id-1) * sizeof(boost::uint32_t);

	//Now, read how many edges our target node has
	boost::uint32_t num_scores;
	binary_file.seekg(advance, std::ios::beg);
	binary_file.read( (char*) &num_scores, sizeof( boost::uint32_t )  );

	if ( TR.Debug.visible() ) {
		TR << "Reading " << num_scores << " edges" <<  std::endl;
	}

	//If this node has no edges, abort
	if ( num_scores == 0 ) {
		return;
	}

	//Now, count up the edges for nodes prior to this node
	binary_file.seekg(sizeof(boost::uint32_t) + sizeof(boost::uint32_t));
	core::Size sum = 0;
	for ( core::Size i=1; i<node_id; ++i ) {
		boost::uint32_t cur_node_count;
		binary_file.read( (char*) &cur_node_count, sizeof( boost::uint32_t )  );
		sum += cur_node_count;
	}

	if ( TR.Debug.visible() ) {
		TR << sum << " edges before this edge" << std::endl;
	}

	//Now, set advance past nodes count, and all the edge counts
	advance = sizeof(boost::uint32_t) + sizeof(boost::uint32_t);
	advance += n_nodes * sizeof(boost::uint32_t);

	//Now, advance past all the actual edges prior to the edges we are interested in
	core::Size edge_size = (3 * sizeof(boost::int32_t)) + (n_segments_per_node * sizeof(boost::int32_t));
	advance += sum * edge_size;
	binary_file.seekg(advance);

	for ( core::Size i=1; i<= num_scores; ++i ) {
		std::set<core::Size> segments;
		for ( core::Size j=1; j<= n_segments_per_node; ++j ) {
			boost::int32_t seg_id;
			binary_file.read( (char*) & seg_id, sizeof(boost::int32_t) );
			segments.insert(seg_id);
		}

		boost::int32_t other_model_id;
		boost::int32_t resnum_1;
		boost::int32_t resnum_2;
		binary_file.read( (char*) & other_model_id, sizeof(boost::int32_t) );
		binary_file.read( (char*) & resnum_1, sizeof(boost::int32_t) );
		binary_file.read( (char*) & resnum_2, sizeof(boost::int32_t) );

		ModelNode const * const node = get_model_node(node_id);

		BasisPair bp;
		bp.first.model_id = node->model().model_id_;
		bp.first.resnum = resnum_1;
		bp.second.model_id = other_model_id;
		bp.second.resnum = resnum_2;

		ModelNode const * const other_node = get_model_node(other_model_id, segments);
		HashEdge * const e = static_cast< HashEdge * >( add_edge(node_id, other_node->get_node_index()));
		e->basis_pair(bp);
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Added " << num_scores << " edges for node " << node_id << std::endl;
	}
}

utility::vector1<BasisPair>
scores_to_alignments(
	ScoreResults const & scores
){
	utility::vector1<BasisPair> basis_pairs;
	for ( ScoreResults::const_iterator it = scores.begin(); it != scores.end(); ++it ) {
		basis_pairs.push_back(it->first);
	}
	return basis_pairs;
}


//utility::vector1<BasisPair>
//read_hashing_scores_from_file(
// std::string filename
//) {
//
// utility::io::izstream file(filename);
// if(!file.good()) {
//  utility_exit_with_message("Could not find Hasher score file with name: " + filename);
// }
//
// utility::vector1<BasisPair> basis_pairs;
// std::string line;
// while( getline( file, line ) ) {
//
//  utility::vector1<std::string> tokens = utility::string_split(line);
//  assert(tokens.size() > 0);
//  //runtime_assert(tokens.size() == 5);
//  Basis b1(utility::string2int(tokens[1]), utility::string2int(tokens[2]));
//  Basis b2(utility::string2int(tokens[3]), utility::string2int(tokens[4]));
//  basis_pairs.push_back(std::make_pair(b1, b2));
// }
// TR.Debug << "Added " << basis_pairs.size() << " alignments from score file: " << filename << std::endl;
// file.close();
// return basis_pairs;
//}

std::string
serialize_graph_json(
	SewGraphCOP graph,
	core::Size max_nodes
){

	std::stringstream nodes;
	nodes << "\"nodes\": [" << std::endl;

	std::stringstream edges;
	edges << "\"edges\": [" << std::endl;

	for ( core::Size i=1; i<=graph->num_nodes() && i<=max_nodes; ++i ) {

		ModelNode const * model_node = graph->get_model_node(i);
		Model const & model = model_node->model();

		nodes << "  { \"data\": {" << std::endl;
		nodes << "    \"id\": \"" << i << "\"," << std::endl;
		nodes << "    \"model_id\": \"" << model.model_id_ << "\"," << std::endl;
		nodes << "    \"pdb_code\": \"" << model.pdb_code_ << "\"" << std::endl;
		nodes << "  } }," << std::endl;

		TR << "Node " << i << " " << model.model_id_ << std::endl;

		core::graph::EdgeListConstIterator edge_it = model_node->const_upper_edge_list_begin();
		core::graph::EdgeListConstIterator edge_it_end = model_node->const_upper_edge_list_end();
		for ( ; edge_it != edge_it_end; ++edge_it ) {

			HashEdge const * const cur_edge = static_cast< HashEdge const * >(*edge_it);
			core::Size other_index = cur_edge->get_other_ind(i);
			Model const & target_model = graph->get_model_node(other_index)->model();
			TR << "Edge " << other_index << " " << target_model.model_id_ << std::endl;

			core::Size source_resnum;
			core::Size target_resnum;
			if ( cur_edge->basis_pair().first.model_id == model.model_id_ ) {
				source_resnum = model.pose_number(cur_edge->basis_pair().first.resnum);
				target_resnum = target_model.pose_number(cur_edge->basis_pair().second.resnum);
			} else {
				source_resnum = model.pose_number(cur_edge->basis_pair().second.resnum);
				target_resnum = target_model.pose_number(cur_edge->basis_pair().first.resnum);
			}
			if ( other_index > max_nodes ) { continue; }

			edges << "  { \"data\": {" << std::endl;
			edges << "    \"source\": \"" << i << "\"," << std::endl;
			edges << "    \"source_resnum\": \"" << source_resnum << "\"," << std::endl;
			edges << "    \"target\": \"" << other_index << "\"," << std::endl;
			edges << "    \"target_resnum\": \"" << target_resnum << "\"" << std::endl;
			edges << "  } }," << std::endl;
		}
	}
	nodes << "]," << std::endl;
	edges << "]" << std::endl;

	std::stringstream elements;
	elements << "{ " << nodes.str() << edges.str() << "}" << std::endl;
	return elements.str();
}

} //sewing namespace
} //protocols namespace
