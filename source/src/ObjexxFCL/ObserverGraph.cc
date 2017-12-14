// ObserverGraph: Observer Graph Representation
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


// ObjexxFCL Headers
#include <ObjexxFCL/ObserverGraph.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/SetWrapper.hh>

// C++ Headers
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <utility>


namespace ObjexxFCL {
namespace internal {


// ObserverGraph: Observer Graph Representation


/// @brief Subject Constructor
ObserverGraph::ObserverGraph( Subject const & s )
{
	// Construct the graph with zero in-degree counts
	if ( ! push( s, s ) ) {
		std::cerr << "\n*** ObjexxFCL Error: " <<
			"Cyclic FArray/Dimension dependency detected" << std::endl;
		std::exit( EXIT_FAILURE );
	}

	// Set the in-degree counts
	for ( auto ig = graph_.begin(), eg = graph_.end(); ig != eg; ++ig ) {
		Observer * const observer_p( ig->first );
		if ( auto * const os_p = dynamic_cast< ObserverSingle * >( observer_p ) ) { // Single Observer
			Observer * const oso_p( os_p->observer_p() );
			if ( oso_p ) {
				assert( graph_.find( oso_p ) != graph_.end() );
				++graph_[ oso_p ]; // Increment the in-degree count
			}
		} else if ( auto * const om_p = dynamic_cast< ObserverMulti * >( observer_p ) ) { // Multi Observer
			if ( om_p->observers_p() ) {
				ObserverMulti::Observers const & observers( om_p->observers() );
				for ( auto io : observers() ) {
					assert( graph_.find( io ) != graph_.end() );
					++graph_[ io ]; // Increment the in-degree count
				}
			}
		}
	}

	// Set the sources
	for ( auto ig = graph_.begin(), eg = graph_.end(); ig != eg; ++ig ) {
		if ( ig->second == 0 ) sources_.push_back( ig ); // In-degree == zero => Source Observer
	}
	assert( ( ! sources_.empty() ) || ( graph_.empty() ) );
}


/// @brief Push a Subject's Transitive Observers onto Graph and Return Acyclicity
bool
ObserverGraph::push( Subject const & s_root, Subject const & s )
{
	if ( auto const * const ss_p = dynamic_cast< SubjectSingle const * >( &s ) ) { // Single Observer
		Observer * const o_p( ss_p->observer_p() );
		if ( o_p ) { // Subject has an Observer
			if ( graph_.find( o_p ) == graph_.end() ) { // New Observer
				graph_.insert( std::make_pair( o_p, static_cast< size_type >( 0 ) ) ); // Add it
				if ( ( o_p == &s ) || ( o_p == &s_root ) ) return false; // Cyclic
				if ( ! push( s_root, *o_p ) ) return false; // Recurse: Abort if cyclic
			}
		}
	} else if ( auto const * const sm_p = dynamic_cast< SubjectMulti const * >( &s ) ) { // Multi Observer
		if ( sm_p->observers_p() ) { // Subject has Observers
			ObserverMulti::Observers const & observers( sm_p->observers() );
			for ( auto o_p : observers() ) {
				if ( graph_.find( o_p ) == graph_.end() ) { // New Observer
					graph_.insert( std::make_pair( o_p, static_cast< size_type >( 0 ) ) ); // Add it
					if ( ( o_p == &s ) || ( o_p == &s_root ) ) return false; // Cyclic
					if ( ! push( s_root, *o_p ) ) return false; // Recurse: Abort if cyclic
				}
			}
		}
	}
	return true;
}


/// @brief Pop a Source Observer from Graph
Observer *
ObserverGraph::pop()
{
	if ( sources_.empty() ) { // No more sources
		assert( graph_.empty() );
		return nullptr;
	} else { // Pop Last source
		Graph::iterator const ig( sources_.back() ); // Last source
		sources_.pop_back(); // Remove the last source

		// Decrement in-degree counts of its Observers
		Observer * const observer_p( ig->first );
		if ( auto * const os_p = dynamic_cast< ObserverSingle * >( observer_p ) ) { // Single Observer
			Observer * const oso_p( os_p->observer_p() );
			if ( oso_p ) {
				Graph::iterator const igo( graph_.find( oso_p ) );
				assert( igo != graph_.end() );
				size_type & in_degree( igo->second );
				assert( in_degree > 0 );
				if ( --in_degree == 0 ) sources_.push_back( igo ); // Decrement the in-degree count / Add to sources if zero
			}
		} else if ( auto * const om_p = dynamic_cast< ObserverMulti * >( observer_p ) ) { // Multi Observer
			if ( om_p->observers_p() ) {
				ObserverMulti::Observers const & observers( om_p->observers() );
				for ( auto io : observers() ) {
					Graph::iterator const igo( graph_.find( io ) );
					assert( igo != graph_.end() );
					size_type & in_degree( igo->second );
					assert( in_degree > 0 );
					if ( --in_degree == 0 ) sources_.push_back( igo ); // Decrement the in-degree count / Add to sources if zero
				}
			}
		}

		graph_.erase( ig ); // Remove the source Observer from the graph

		return observer_p;
	}
}


// ObserverGraph


} // namespace internal
} // namespace ObjexxFCL
