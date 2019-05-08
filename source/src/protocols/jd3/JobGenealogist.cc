// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobGenealogist.cc
/// @brief  class method definitions for JobGenealogist
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Unit headers
#include <protocols/jd3/JobGenealogist.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/vector1.functions.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.jd3.JobGenealogist" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
//#include <cereal/types/hash.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/unordered_set.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {

//////////
//JGJobNode

JGJobNode::JGJobNode() :
	global_job_id_( 0 ),
	node_( 0 ),
	input_source_id_( 0 ),
	parents_( 0 ),
	children_( 0 )
{}

JGJobNode::JGJobNode(
	core::Size global_job,
	unsigned int job_dag_node,
	JGResultNodeAP parent_node,
	unsigned int input_source
) :
	global_job_id_( global_job ),
	node_( job_dag_node ),
	input_source_id_( input_source ),
	parents_( 1, parent_node ),
	children_( 0 )
{}

//JGJobNode::~JGJobNode(){};

void
JGJobNode::add_parent( JGResultNodeAP p, bool tell_parent_to_add_child ){
	parents_.push_back( p );
	if ( tell_parent_to_add_child ) {
		JGResultNodeOP const p_op = p.lock();
		debug_assert( p_op != nullptr );
		p_op->add_child( shared_from_this(), false );
	}
}

bool
JGJobNode::remove_parent( JGResultNodeAP p, bool tell_parent_to_remove_child ){
	JGResultNodeOP p_op = p.lock();

	core::Size zero_indexed_element_location = 0;
	for ( JGResultNodeAP const & parent_ap : parents_ ) {
		JGResultNodeOP const parent = parent_ap.lock();

		if ( parent != p_op ) {
			++zero_indexed_element_location;
			continue;
		}

		if ( tell_parent_to_remove_child ) {
			parent->remove_child( shared_from_this(), false );
		}

		//update input_source_id to new primary parent
		if ( zero_indexed_element_location == 0 && parents_.size() > 1 ) {
			input_source_id_ = parents_[ 2 ].lock()->parent().lock()->input_source_id_;
		}
		parents_.erase( std::next( parents_.begin(), zero_indexed_element_location ) );
		return true;
	}

	return false;
}

JGResultNode::JGResultNode() :
	result_id_( 0 ),
	parent_(),
	children_( 0 )
{}

JGResultNode::JGResultNode( unsigned int result, JGJobNodeAP par ) :
	result_id_( result ),
	parent_( par ),
	children_( 0 )
{}

void
JGResultNode::add_child( JGJobNodeAP c, bool tell_child_to_add_parent ){
	children_.push_back( c );
	if ( tell_child_to_add_parent ) {
		c.lock()->add_parent( shared_from_this(), false );
	}
}

bool
JGResultNode::remove_child( JGJobNodeAP c, bool tell_child_to_remove_parent ){
	JGJobNodeOP c_op = c.lock();

	core::Size zero_indexed_element_location = 0;
	for ( JGJobNodeAP const & child_ap : children_ ) {
		JGJobNodeOP const child = child_ap.lock();

		if ( child != c_op ) {
			++zero_indexed_element_location;
			continue;
		}

		if ( tell_child_to_remove_parent ) {
			child->remove_parent( shared_from_this(), false );
		}

		children_.erase( std::next( children_.begin(), zero_indexed_element_location ) );
		return true;
	}

	return false;
}



JobGenealogist::JobGenealogist(
	core::Size num_nodes,
	core::Size num_input_sources
) :
	num_input_sources_( num_input_sources )
{
	job_nodes_for_dag_node_.resize( num_nodes );
	all_result_nodes_.max_load_factor( 0.5 );
	all_job_nodes_.max_load_factor( 0.5 );
}

JobGenealogist::~JobGenealogist() {}

JGJobNodeOP
JobGenealogist::register_new_job(
	core::Size const job_dag_node_id,
	core::Size const global_job_id,
	core::Size const input_source_id
){

	debug_assert( num_input_sources_ >= input_source_id );

	JGJobNodeOP new_node = utility::pointer::make_shared< JGJobNode >();
	all_job_nodes_.insert( new_node );

	new_node->node( job_dag_node_id );
	new_node->global_job_id( global_job_id );
	new_node->input_source_id( input_source_id );

	auto iter = std::lower_bound(
		job_nodes_for_dag_node_[ job_dag_node_id ].begin(),
		job_nodes_for_dag_node_[ job_dag_node_id ].end(),
		new_node,
		sorter_
	);
	job_nodes_for_dag_node_[ job_dag_node_id ].insert( iter, new_node );

	return new_node;
}

JGJobNodeOP
JobGenealogist::register_new_job(
	core::Size const job_dag_node_id,
	core::Size const global_job_id,
	core::Size const job_dag_node_id_of_parent,
	core::Size const global_job_id_of_parent,
	core::Size const result_id_of_parent
){
	JGJobNodeOP new_node = utility::pointer::make_shared< JGJobNode >();
	all_job_nodes_.insert( new_node );

	new_node->node( job_dag_node_id );
	new_node->global_job_id( global_job_id );

	debug_assert( job_dag_node_id_of_parent && global_job_id_of_parent && result_id_of_parent );

	JGResultNodeOP const parent =
		get_result_node( job_dag_node_id_of_parent, global_job_id_of_parent, result_id_of_parent ).lock();
	debug_assert( parent != nullptr );
	new_node->parents().push_back( parent );
	parent->children().push_back( new_node );

	/*debug_assert(
	std::find(
	parent->children().begin(),
	parent->children().end(),
	new_node
	) != parent->children().end()
	);*/

	new_node->input_source_id( parent->parent().lock()->input_source_id() );

	auto iter = std::lower_bound(
		job_nodes_for_dag_node_[ job_dag_node_id ].begin(),
		job_nodes_for_dag_node_[ job_dag_node_id ].end(),
		new_node,
		sorter_
	);
	job_nodes_for_dag_node_[ job_dag_node_id ].insert( iter, new_node );

	return new_node;
}

JGJobNodeOP
JobGenealogist::register_new_job(
	core::Size const job_dag_node_id,
	core::Size const global_job_id,
	JGResultNodeAP const parent
){
	JGJobNodeOP new_node = utility::pointer::make_shared< JGJobNode >();
	all_job_nodes_.insert( new_node );

	new_node->node( job_dag_node_id );
	new_node->global_job_id( global_job_id );
	new_node->input_source_id( parent.lock()->parent().lock()->input_source_id() );

	new_node->add_parent( parent, true );

	auto iter = std::lower_bound(
		job_nodes_for_dag_node_[ job_dag_node_id ].begin(),
		job_nodes_for_dag_node_[ job_dag_node_id ].end(),
		new_node,
		sorter_
	);
	job_nodes_for_dag_node_[ job_dag_node_id ].insert( iter, new_node );

	return new_node;
}

JGJobNodeOP
JobGenealogist::register_new_job(
	core::Size const job_dag_node_id,
	core::Size const global_job_id,
	utility::vector1< JGResultNodeAP > const & parents
){
	JGJobNodeOP new_node = utility::pointer::make_shared< JGJobNode >();
	all_job_nodes_.insert( new_node );

	new_node->node( job_dag_node_id );
	new_node->global_job_id( global_job_id );

	new_node->parents().reserve( parents.size() );
	for ( JGResultNodeAP ii : parents ) {
		debug_assert( ! ii.expired() );
		new_node->add_parent( ii, true );
	}

	//Note: not all parent are created equal! new_node gets its identity from it's first parent
	new_node->input_source_id( new_node->parents().front().lock()->parent().lock()->input_source_id() );

	auto iter = std::lower_bound(
		job_nodes_for_dag_node_[ job_dag_node_id ].begin(),
		job_nodes_for_dag_node_[ job_dag_node_id ].end(),
		new_node,
		sorter_
	);
	job_nodes_for_dag_node_[ job_dag_node_id ].insert( iter, new_node );

	return new_node;
}


void JobGenealogist::note_job_completed(
	JGJobNodeAP const job_node,
	core::Size const nresults
){
	for ( core::Size ii = 1; ii <= nresults; ++ii ) {
		JGResultNodeOP new_node = utility::pointer::make_shared< JGResultNode >();
		all_result_nodes_.insert( new_node );

		new_node->result_id( ii );
		new_node->parent( job_node );
		job_node.lock()->children().push_back( new_node );
	}
}

void JobGenealogist::discard_job_result(
	core::Size const node,
	core::Size const global_job_id,
	core::Size const result_id
){
	JGResultNodeAP const result_node = get_result_node( node, global_job_id, result_id );
	debug_assert( ! result_node.expired() );
	JGResultNodeOP const result_node_op = result_node.lock();

	JGJobNodeOP const parent = result_node.lock()->parent().lock();
	debug_assert( parent != nullptr );

	for ( core::Size ii = 1; ii <= parent->children().size(); ++ii ) {
		if ( parent->children()[ ii ].lock() == result_node_op ) {
			parent->children().erase( std::next( parent->children().begin(), ii - 1 ) );
			break;
		}
	}

	delete_node( result_node );
}


void JobGenealogist::garbage_collection(
	core::Size const node,
	bool const delete_downstream_job_if_it_has_no_results,
	std::list< jd3::JobResultID > & container_for_discarded_result_ids
) {
	for ( core::Size aa = 1; aa <= job_nodes_for_dag_node_[ node ].size(); ) {

		auto const & job_node = job_nodes_for_dag_node_[ node ][ aa ];

		//prune children
		for ( core::Size ii = 1; ii <= job_node->children().size(); ) {
			JGResultNodeOP result_node = job_node->children()[ ii ].lock();

			if ( delete_downstream_job_if_it_has_no_results ) {
				for ( core::Size jj = 1; jj <= result_node->children().size(); ) {
					JGJobNodeOP grandchild = result_node->children()[ jj ].lock();
					if ( grandchild->children().empty() ) {
						//stack overflow told me this is faster than erase:
						std::swap( result_node->children()[ jj ], result_node->children().back() );
						result_node->children().pop_back();
						delete_node( grandchild );
					} else {
						++jj;
					}
				}
			}

			if ( result_node->children().empty() ) {
				container_for_discarded_result_ids.emplace_back( job_node->global_job_id(), result_node->result_id() );

				//stack overflow told me this is faster than erase:
				std::swap( job_node->children()[ ii ], job_node->children().back() );
				job_node->children().pop_back();
				delete_node( result_node );
			} else {
				++ii;
			}
		}//ii

		if ( job_node->children().empty() ) {
			if ( ! job_node->parents().empty() ) {
				for ( JGResultNodeAP parent : job_node->parents() ) {
					utility::vector1< JGJobNodeAP > & siblings = parent.lock()->children();
					for ( core::Size ii = 1; ii <= siblings.size(); ++ii ) {
						if ( siblings[ ii ].lock() == job_node ) {
							siblings.erase( std::next( siblings.begin(), ii - 1 ) );
							break;
						}
					}
				}//for parent
			} // if ! job_node->parents().empty()

			delete_node( job_node, 0, false );
			job_nodes_for_dag_node_[ node ].erase( job_nodes_for_dag_node_[ node ].begin() + aa - 1 );
		} else {
			++aa;
		}
	}

#ifndef NDEBUG
	TR << "printing nodes after garbage collection on node " << node << std::endl;
	print_all_nodes();
#endif

}

void JobGenealogist::all_job_results_for_node(
	core::Size const job_dag_node,
	std::list< jd3::JobResultID > & container_for_output
) const {
	for ( JGJobNodeCOP const & job_node : job_nodes_for_dag_node_[ job_dag_node ] ) {
		for ( JGResultNodeCAP const & result_node : job_node->children() ) {
			container_for_output.emplace_back( job_node->global_job_id(), result_node.lock()->result_id() );
		}
	}
}

void JobGenealogist::add_newick_tree_for_node(
	JGResultNodeCAP const result_node,
	std::stringstream & ss
) const {

	JGResultNodeCOP result_node_cop = result_node.lock();

	std::list< JGResultNodeCAP > spawned_results;
	for ( JGJobNodeCAP child : result_node_cop->children() ) {
		debug_assert( ! child.expired() );
		for ( JGResultNodeCAP grandchild : child.lock()->children() ) {
			spawned_results.push_back( grandchild );
		}
	}

	if ( spawned_results.empty() ) {
		ss << "JR_" << result_node_cop->parent().lock()->global_job_id() << "_" << result_node_cop->result_id();
		return;
	} else {
		ss << "(";
		for ( JGResultNodeCAP spawned_result : spawned_results ) {
			add_newick_tree_for_node( spawned_result, ss );
			ss << ",";
		}
		ss.seekp( -1, ss.cur );//remove last comma https://stackoverflow.com/questions/4546021/remove-char-from-stringstream-and-append-some-data/26492431#26492431
		ss << ")JR_" << result_node_cop->parent().lock()->global_job_id() << "_" << result_node_cop->result_id();
	}
}

std::string
JobGenealogist::newick_tree() const {

	//collect starting nodes
	utility::vector1< std::list< JGResultNodeCAP > > parentless_results_for_input_source;
	parentless_results_for_input_source.resize( num_input_sources_ );

	for ( unsigned int i = 1; i < job_nodes_for_dag_node_.size(); ++i ) {
		for ( JGJobNodeCAP const & job_node : job_nodes_for_dag_node_[ i ] ) {
			JGJobNodeCOP job_node_cop = job_node.lock();
			debug_assert( job_node_cop != nullptr );
			if ( job_node_cop->parents().empty() ) {
				unsigned int const input_src = job_node_cop->input_source_id();
				for ( JGResultNodeCAP const & result_node : job_node_cop->children() ) {
					parentless_results_for_input_source[ input_src ].push_back( result_node );
				}
			}
		}
	}

	//traverse
	std::stringstream ss;
	ss << "(";

	for ( unsigned int input_src = 1; input_src <= parentless_results_for_input_source.size(); ++input_src ) {
		if ( parentless_results_for_input_source[ input_src ].empty() ) {
			ss << "input_source_" << input_src;
			if ( input_src != parentless_results_for_input_source.size() ) ss << ",";
			continue;
		}

		ss << "(";

		unsigned int counter = 0;
		for ( JGResultNodeCAP result_node : parentless_results_for_input_source[ input_src ] ) {
			add_newick_tree_for_node( result_node, ss );
			if ( ++counter != parentless_results_for_input_source[ input_src ].size() ) {
				ss << ",";
			}
		}
		//ss.seekp( -1, ss.cur );//remove last comma https://stackoverflow.com/questions/4546021/remove-char-from-stringstream-and-append-some-data/26492431#26492431
		ss << ")input_source_" << input_src;
		if ( input_src != parentless_results_for_input_source.size() ) ss << ",";
	}

	ss << ")all";
	return ss.str();
}

void
JobGenealogist::print_all_nodes(){
	for ( core::Size node = 1; node <= job_nodes_for_dag_node_.size(); ++node ) {
		TR << "node " << node << " : ";
		for ( JGJobNodeCAP const & job_node : job_nodes_for_dag_node_[ node ] ) {
			TR << " " << job_node.lock()->global_job_id();
		}
		TR << std::endl;
	}
}


JGJobNodeOP
JobGenealogist::get_job_node(
	core::Size job_dag_node,
	core::Size global_job_id
){
	JGJobNode dummy;
	dummy.global_job_id( global_job_id );
	auto iter = std::lower_bound(
		job_nodes_for_dag_node_[ job_dag_node ].begin(),
		job_nodes_for_dag_node_[ job_dag_node ].end(),
		& dummy,
		sorter_
	);

	if ( iter == job_nodes_for_dag_node_[ job_dag_node ].end() ) {
		return nullptr;
	} else if ( (*iter)->global_job_id() != global_job_id ) {
		return nullptr;
	} else {
		return *iter;
	}
}

JGJobNodeCOP
JobGenealogist::get_const_job_node(
	core::Size job_dag_node,
	core::Size global_job_id
) const {
	JGJobNode dummy;
	dummy.global_job_id( global_job_id );
	auto iter = std::lower_bound(
		job_nodes_for_dag_node_[ job_dag_node ].begin(),
		job_nodes_for_dag_node_[ job_dag_node ].end(),
		& dummy,
		sorter_
	);

	if ( iter == job_nodes_for_dag_node_[ job_dag_node ].end() ) {
		return nullptr;
	} else if ( (*iter)->global_job_id() != global_job_id ) {
		return nullptr;
	} else {
		return *iter;
	}
}


JGResultNodeAP
JobGenealogist::get_result_node(
	core::Size const node,
	core::Size const global_job_id,
	core::Size const result_id
) {
	JGJobNodeAP job_node = get_job_node( node, global_job_id );
	debug_assert( ! job_node.expired() );
	for ( JGResultNodeAP child : job_node.lock()->children() ) {
		if ( child.lock()->result_id() == result_id ) return child;
	}
	return JGResultNodeAP();//blank
}

JGResultNodeCAP
JobGenealogist::get_const_result_node(
	core::Size const node,
	core::Size const global_job_id,
	core::Size const result_id
) const {
	JGJobNodeCAP job_node = get_const_job_node( node, global_job_id );
	for ( JGResultNodeCAP child : job_node.lock()->children() ) {
		if ( child.lock()->result_id() == result_id ) return child;
	}
	return JGResultNodeCAP();
}

void
JobGenealogist::delete_node(
	JGJobNodeAP const job_node,
	unsigned int const job_dag_node,
	bool const delete_from_vec
){
	JGJobNodeOP job_node_op = job_node.lock();

	if ( delete_from_vec ) {
		//job_nodes_for_dag_node_[ job_dag_node ].erase( job_node_op );
		JGJobNode dummy;
		dummy.global_job_id( job_node_op->global_job_id() );
		auto iter = std::lower_bound(
			job_nodes_for_dag_node_[ job_dag_node ].begin(),
			job_nodes_for_dag_node_[ job_dag_node ].end(),
			& dummy,
			sorter_
		);
		if ( *iter == job_node_op ) {
			job_nodes_for_dag_node_[ job_dag_node ].erase( iter );
		}
	}

	all_job_nodes_.erase( job_node_op );
}

} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::JGResultNode::save( Archive & arc ) const {
	arc( CEREAL_NVP( result_id_ ) ); // unsigned int
	arc( CEREAL_NVP( parent_ ) ); // JGJobNodeAP
	arc( CEREAL_NVP( children_ ) ); // utility::vector1<JGJobNodeAP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::JGResultNode::load( Archive & arc ) {
	arc( result_id_ ); // unsigned int
	arc( parent_ ); // JGJobNodeAP
	arc( children_ ); // utility::vector1<JGJobNodeAP>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::JGResultNode );

/// @brief Default constructor required by cereal to deserialize this class
protocols::jd3::JobGenealogist::JobGenealogist() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::JobGenealogist::save( Archive & arc ) const {
	arc( CEREAL_NVP( num_input_sources_ ) ); // core::Size
	arc( CEREAL_NVP( job_nodes_for_dag_node_ ) ); // utility::vector1<utility::vector1<JGJobNodeOP> >
	arc( CEREAL_NVP( all_job_nodes_ ) ); // std::unordered_set<JGJobNodeOP>
	arc( CEREAL_NVP( all_result_nodes_ ) ); // std::unordered_set<JGResultNodeOP>
	arc( CEREAL_NVP( sorter_ ) ); // struct protocols::jd3::compare_job_nodes
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::JobGenealogist::load( Archive & arc ) {
	arc( num_input_sources_ ); // core::Size
	arc( job_nodes_for_dag_node_ ); // utility::vector1<utility::vector1<JGJobNodeOP> >
	arc( all_job_nodes_ ); // std::unordered_set<JGJobNodeOP>
	arc( all_result_nodes_ ); // std::unordered_set<JGResultNodeOP>
	arc( sorter_ ); // struct protocols::jd3::compare_job_nodes
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::JobGenealogist );
CEREAL_REGISTER_TYPE( protocols::jd3::JobGenealogist )


/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::compare_job_nodes::save( Archive & arc ) const {
	arc( cereal::base_class< struct std::binary_function<const class protocols::jd3::JGJobNode *, const class protocols::jd3::JGJobNode *, bool> >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::compare_job_nodes::load( Archive & arc ) {
	arc( cereal::base_class< struct std::binary_function<const class protocols::jd3::JGJobNode *, const class protocols::jd3::JGJobNode *, bool> >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::compare_job_nodes );

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::JGJobNode::save( Archive & arc ) const {
	arc( CEREAL_NVP( global_job_id_ ) ); // core::Size
	arc( CEREAL_NVP( node_ ) ); // unsigned int
	arc( CEREAL_NVP( input_source_id_ ) ); // unsigned int
	arc( CEREAL_NVP( parents_ ) ); // utility::vector1<JGResultNodeAP>
	arc( CEREAL_NVP( children_ ) ); // utility::vector1<JGResultNodeAP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::JGJobNode::load( Archive & arc ) {
	arc( global_job_id_ ); // core::Size
	arc( node_ ); // unsigned int
	arc( input_source_id_ ); // unsigned int
	arc( parents_ ); // utility::vector1<JGResultNodeAP>
	arc( children_ ); // utility::vector1<JGResultNodeAP>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::JGJobNode );
CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_JobGenealogist )
#endif // SERIALIZATION
