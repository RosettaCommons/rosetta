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
/// @author Jack Maguire, jack@med.unc.edu

// Unit headers
#include <protocols/jd3/JobGenealogist.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/vector1.functions.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.jd3.JobGenealogist" );

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
	JGResultNode * parent_node,
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
JGJobNode::add_parent( JGResultNode * p, bool tell_parent_to_add_child ){
	parents_.push_back( p );
	if ( tell_parent_to_add_child ) {
		p->add_child( this, false );
	}
}

bool
JGJobNode::remove_parent( JGResultNode * p, bool tell_parent_to_remove_child ){
	auto iter = std::find( parents_.begin(), parents_.end(), p );
	if ( iter == parents_.end() ) {
		return false;
	} else {
		if ( tell_parent_to_remove_child ) {
			(*iter)->remove_child( this, false );
		}
		//update input_source_id to new primary parent
		if ( iter == parents_.begin() && parents_.size() > 1 ) {
			input_source_id_ = parents_[ 2 ]->parent()->input_source_id_;
		}
		parents_.erase( iter );
		return true;
	}
}

JGResultNode::JGResultNode() :
	result_id_( 0 ),
	parent_( 0 ),
	children_( 0 )
{}

JGResultNode::JGResultNode( unsigned int result, JGJobNode * par ) :
	result_id_( result ),
	parent_( par ),
	children_( 0 )
{}

void JGResultNode::add_child( JGJobNode * c, bool tell_child_to_add_parent ){
	children_.push_back( c );
	if ( tell_child_to_add_parent ) {
		c->add_parent( this, false );
	}
}

bool JGResultNode::remove_child( JGJobNode * c, bool tell_child_to_remove_parent ){
	auto const iter = std::find( children_.begin(), children_.end(), c );
	if ( iter == children_.end() ) {
		return false;
	} else {
		if ( tell_child_to_remove_parent ) {
			(*iter)->remove_parent( this, false );
		}
		children_.erase( iter );
		return true;
	}
}



JobGenealogist::JobGenealogist(
	core::Size num_nodes,
	core::Size num_input_sources
) :
	num_input_sources_( num_input_sources )
{
	job_nodes_for_dag_node_.resize( num_nodes );
}

JobGenealogist::~JobGenealogist() {}

JGJobNode * JobGenealogist::register_new_job(
	core::Size job_dag_node_id,
	core::Size global_job_id,
	core::Size input_source_id
){

	debug_assert( num_input_sources_ >= input_source_id );

	JGJobNode * new_node = job_node_pool_.construct();
	new_node->node( job_dag_node_id );
	new_node->global_job_id( global_job_id );
	new_node->input_source_id( input_source_id );

	auto iter = std::lower_bound(
		job_nodes_for_dag_node_[ job_dag_node_id ].begin(),
		job_nodes_for_dag_node_[ job_dag_node_id ].end(),
		new_node,
		sorter
	);
	job_nodes_for_dag_node_[ job_dag_node_id ].insert( iter, new_node );

	return new_node;
}

JGJobNode * JobGenealogist::register_new_job(
	core::Size job_dag_node_id,
	core::Size global_job_id,
	core::Size job_dag_node_id_of_parent,
	core::Size global_job_id_of_parent,
	core::Size result_id_of_parent
){
	JGJobNode * new_node = job_node_pool_.construct();
	new_node->node( job_dag_node_id );
	new_node->global_job_id( global_job_id );

	debug_assert( job_dag_node_id_of_parent && global_job_id_of_parent && result_id_of_parent );

	JGResultNode * const parent = get_result_node( job_dag_node_id_of_parent, global_job_id_of_parent, result_id_of_parent );
	debug_assert( parent );
	new_node->parents().push_back( parent );
	parent->children().push_back( new_node );

	debug_assert(
		std::find(
		parent->children().begin(),
		parent->children().end(),
		new_node
		) != parent->children().end()
	);

	new_node->input_source_id( parent->parent()->input_source_id() );

	auto iter = std::lower_bound(
		job_nodes_for_dag_node_[ job_dag_node_id ].begin(),
		job_nodes_for_dag_node_[ job_dag_node_id ].end(),
		new_node,
		sorter
	);
	job_nodes_for_dag_node_[ job_dag_node_id ].insert( iter, new_node );

	return new_node;
}

JGJobNode * JobGenealogist::register_new_job(
	core::Size job_dag_node_id,
	core::Size global_job_id,
	JGResultNode * parent
){
	JGJobNode * new_node = job_node_pool_.construct();
	new_node->node( job_dag_node_id );
	new_node->global_job_id( global_job_id );
	new_node->input_source_id( parent->parent()->input_source_id() );

	new_node->add_parent( parent, true );

	auto iter = std::lower_bound(
		job_nodes_for_dag_node_[ job_dag_node_id ].begin(),
		job_nodes_for_dag_node_[ job_dag_node_id ].end(),
		new_node,
		sorter
	);
	job_nodes_for_dag_node_[ job_dag_node_id ].insert( iter, new_node );

	return new_node;
}

JGJobNode * JobGenealogist::register_new_job(
	core::Size job_dag_node_id,
	core::Size global_job_id,
	utility::vector1< JGResultNode * > const & parents
){
	JGJobNode * new_node = job_node_pool_.construct();
	new_node->node( job_dag_node_id );
	new_node->global_job_id( global_job_id );

	new_node->parents().reserve( parents.size() );
	for ( JGResultNode * ii : parents ) {
		debug_assert( ii );
		new_node->add_parent( ii, true );
	}

	//Note: not all parent are created equal! new_node gets its identity from it's first parent
	new_node->input_source_id( new_node->parents().front()->parent()->input_source_id() );

	auto iter = std::lower_bound(
		job_nodes_for_dag_node_[ job_dag_node_id ].begin(),
		job_nodes_for_dag_node_[ job_dag_node_id ].end(),
		new_node,
		sorter
	);
	job_nodes_for_dag_node_[ job_dag_node_id ].insert( iter, new_node );

	return new_node;
}


void JobGenealogist::note_job_completed(
	JGJobNode * job_node,
	core::Size nresults
){
	for ( core::Size ii = 1; ii <= nresults; ++ii ) {
		JGResultNode * new_node = result_node_pool_.construct();
		debug_assert( new_node );

		new_node->result_id( ii );
		new_node->parent( job_node );
		job_node->children().push_back( new_node );
	}
}

void JobGenealogist::discard_job_result(
	core::Size node,
	core::Size global_job_id,
	core::Size result_id
){
	JGResultNode * const result_node = get_result_node( node, global_job_id, result_id );
	debug_assert( result_node );

	JGJobNode * const parent = result_node->parent();
	debug_assert( parent );

	auto iter = std::find( parent->children().begin(), parent->children().end(), result_node );
	debug_assert( iter != parent->children().end() );
	parent->children().erase( iter );

	delete_node( result_node );
}


void JobGenealogist::garbage_collection(
	core::Size const node,
	bool const delete_downstream_job_if_it_has_no_results,
	std::list< jd3::JobResultID > & container_for_discarded_result_ids
) {
	for ( core::Size aa = 1; aa <= job_nodes_for_dag_node_[ node ].size(); ) {
		//for ( auto it = job_nodes_for_dag_node_[ node ].begin(); it != job_nodes_for_dag_node_[ node ].end(); ) {
		JGJobNode * const job_node = job_nodes_for_dag_node_[ node ][ aa ];

		//prune children
		for ( core::Size ii = 1; ii <= job_node->children().size(); ) {
			JGResultNode * result_node = job_node->children()[ ii ];

			if ( delete_downstream_job_if_it_has_no_results ) {
				for ( core::Size jj = 1; jj <= result_node->children().size(); ) {
					JGJobNode * grandchild = result_node->children()[ jj ];
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
				for ( JGResultNode * parent : job_node->parents() ) {
					utility::vector1< JGJobNode * > & siblings = parent->children();
					auto pos = std::find( siblings.begin(), siblings.end(), job_node );
					runtime_assert( pos != siblings.end() );
					siblings.erase( pos );
				}
			}

			job_nodes_for_dag_node_[ node ].erase( job_nodes_for_dag_node_[ node ].begin() + aa - 1 );
			delete_node( job_node, 0, false );
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
	core::Size job_dag_node,
	std::list< jd3::JobResultID > & container_for_output
) const {
	for ( auto it = job_nodes_for_dag_node_[ job_dag_node ].begin(); it != job_nodes_for_dag_node_[ job_dag_node ].end(); ++it ) {
		JGJobNode const * const job_node = *it;
		for ( core::Size ii = 1; ii <= job_node->children().size(); ++ii ) {
			JGResultNode const * const result_node = job_node->children()[ ii ];
			container_for_output.push_back( std::make_pair( job_node->global_job_id(), result_node->result_id() ) );
		}
	}
}

void JobGenealogist::add_newick_tree_for_node(
	JGResultNode const * result_node,
	std::stringstream & ss
) const {

	std::list< JGResultNode const * > spawned_results;
	for ( JGJobNode const * child : result_node->children() ) {
		if ( child ) {
			for ( JGResultNode const * grandchild : child->children() ) {
				spawned_results.push_back( grandchild );
			}
		}
	}

	if ( spawned_results.empty() ) {
		ss << "JR_" << result_node->parent()->global_job_id() << "_" << result_node->result_id();
		return;
	} else {
		ss << "(";
		for ( JGResultNode const * spawned_result : spawned_results ) {
			add_newick_tree_for_node( spawned_result, ss );
			ss << ",";
		}
		ss.seekp( -1, ss.cur );//remove last comma https://stackoverflow.com/questions/4546021/remove-char-from-stringstream-and-append-some-data/26492431#26492431
		ss << ")JR_" << result_node->parent()->global_job_id() << "_" << result_node->result_id();
	}
}

std::string JobGenealogist::newick_tree() const {

	//collect starting nodes
	utility::vector1< std::list< JGResultNode const * > > parentless_results_for_input_source;
	parentless_results_for_input_source.resize( num_input_sources_ );

	for ( unsigned int i = 1; i < job_nodes_for_dag_node_.size(); ++i ) {
		for ( JGJobNode const * const job_node : job_nodes_for_dag_node_[ i ] ) {
			debug_assert( job_node );
			if ( job_node->parents().empty() ) {
				unsigned int const input_src = job_node->input_source_id();
				for ( JGResultNode const * result_node : job_node->children() ) {
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
		for ( JGResultNode const * result_node : parentless_results_for_input_source[ input_src ] ) {
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

void JobGenealogist::print_all_nodes(){
	for ( core::Size node = 1; node <= job_nodes_for_dag_node_.size(); ++node ) {
		TR << "node " << node << " : ";
		for ( JGJobNode const * const job_node : job_nodes_for_dag_node_[ node ] ) {
			TR << " " << job_node->global_job_id();
		}
		TR << std::endl;
	}
}


JGJobNode * JobGenealogist::get_job_node( core::Size job_dag_node, core::Size global_job_id ){
	JGJobNode dummy;
	dummy.global_job_id( global_job_id );
	auto iter = std::lower_bound(
		job_nodes_for_dag_node_[ job_dag_node ].begin(),
		job_nodes_for_dag_node_[ job_dag_node ].end(),
		& dummy,
		sorter
	);

	if ( iter == job_nodes_for_dag_node_[ job_dag_node ].end() ) {
		return 0;
	} else if ( (*iter)->global_job_id() != global_job_id ) {
		return 0;
	} else {
		return *iter;
	}
}


JGJobNode const * JobGenealogist::get_job_node( core::Size job_dag_node, core::Size global_job_id ) const {
	JGJobNode dummy;
	dummy.global_job_id( global_job_id );
	auto iter = std::lower_bound(
		job_nodes_for_dag_node_[ job_dag_node ].begin(),
		job_nodes_for_dag_node_[ job_dag_node ].end(),
		& dummy,
		sorter
	);

	if ( iter == job_nodes_for_dag_node_[ job_dag_node ].end() ) {
		return 0;
	} else if ( (*iter)->global_job_id() != global_job_id ) {
		return 0;
	} else {
		return *iter;
	}
}

JGResultNode * JobGenealogist::get_result_node( core::Size node, core::Size global_job_id, core::Size result_id ) {
	JGJobNode * job_node = get_job_node( node, global_job_id );
	debug_assert( job_node );
	//TODO std::find_if() ?
	for ( JGResultNode * child : job_node->children() ) {
		if ( child->result_id() == result_id ) return child;
	}
	return 0;
}

JGResultNode const * JobGenealogist::get_result_node( core::Size node, core::Size global_job_id, core::Size result_id ) const {
	JGJobNode const * job_node = get_job_node( node, global_job_id );
	//TODO std::find_if() ?
	for ( JGResultNode const * child : job_node->children() ) {
		if ( child->result_id() == result_id ) return child;
	}
	return 0;
}

void JobGenealogist::delete_node( JGJobNode * job_node, unsigned int job_dag_node, bool delete_from_vec ){
	if ( delete_from_vec ) {
		JGJobNode dummy;
		dummy.global_job_id( job_node->global_job_id() );
		auto iter = std::lower_bound(
			job_nodes_for_dag_node_[ job_dag_node ].begin(),
			job_nodes_for_dag_node_[ job_dag_node ].end(),
			& dummy,
			sorter
		);
		if ( *iter == job_node ) {
			job_nodes_for_dag_node_[ job_dag_node ].erase( iter );
		}
	}
	job_node_pool_.destroy( job_node );
}

} // namespace jd3
} // namespace protocols

