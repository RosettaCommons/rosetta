// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/LoopGraph.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/loop_graph/LoopGraph.hh>
#include <core/scoring/loop_graph/Loop.hh>
#include <core/scoring/loop_graph/LoopCycle.hh>
#include <core/scoring/loop_graph/evaluator/LoopClosePotentialEvaluator.hh>
#include <core/scoring/loop_graph/util.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/GaussianChainFunc.hh>
#include <core/chemical/VariantType.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// for figuring out elementary cycles
#include <boost/graph/directed_graph.hpp>
#include <boost/graph/tiernan_all_cycles.hpp>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/numbers.hh> // for is_finite
#include <utility/tools/make_vector1.hh>

static basic::Tracer TR( "core.scoring.loop_graph.LoopGraph" );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// @details
//
// For Stepwise modeling, full model has several instantiated poses whose covalent loop interconnections
//  have not been instantiated. This LoopGraph keeps track of those loops, and any that form cycles
//  lead to a 'loop_close' term (see LoopCloseEnergy in this namespace).
//
// Full models like this are OK:
//
//            Loop
//       ---- ~ ~ ~ -- ~ ~ ~
//       ||||  cycA || cycB  ~  Loop
//       ---- ~ ~ ~ -- ~ ~ ~
//     POSE 1  |     POSE 2
//            Loop
//
//         ~ ~ = loops that are not created yet.
//
//
// Need to specify loop_fixed_cost, which is the 'typical' restriction volume for the nucleotides at the ends of the loop, in Angstrom^3.
//   [For the 6D potentials, interpret this loop_fixed_cost as the restricted 6D volume = translational volume x [fraction of SO(3)].
//
// Note: models like this are currently a problem:
//
//     POSE 1      POSE 2
//       ---- ~ ~ ~  -- ~ ~
//       |||| cycA   ||     ~
//       ---- ~ -- ~ -- ~   ~
//   POSE 3 --> ||  cycB  ~  ~
//            ~ -- ~ ~ ~ ~   ~
//           ~     cycC     ~
//            ~ ~ ~ ~ ~ ~ ~
//
//  Note that some loops are shared between cycles -- that's the issue.
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace core::scoring::loop_graph::evaluator;

namespace core {
namespace scoring {
namespace loop_graph {

//Constructor
LoopGraph::LoopGraph():
	total_energy_( 0.0 ),
	loop_fixed_cost_( basic::options::option[ basic::options::OptionKeys::score::loop_close::loop_fixed_cost ]() ), // -0.29 default, in Rosetta energy units
	error_out_on_complex_cycles_( !basic::options::option[ basic::options::OptionKeys::score::loop_close::allow_complex_loop_graph ]() ),
	has_just_simple_cycles_( true ),
	use_6D_potential_( basic::options::option[ basic::options::OptionKeys::score::loop_close::use_6D_potential ]() )
{
}

//Destructor
LoopGraph::~LoopGraph() = default;

/////////////////////////////////////////////////////////////////////
void
LoopGraph::update( pose::Pose & pose, bool const verbose /* = false */ ){

	using namespace core::pose::full_model_info;
	utility::vector1< Size > const pose_domain_map = figure_out_pose_domain_map( pose );
	utility::vector1< Size > const & cutpoint_open = const_full_model_info( pose ).cutpoint_open_in_full_model();

	update_loops_and_cycles( pose_domain_map, cutpoint_open );

	current_pose_loop_score_evaluators_.clear();
	total_energy_ = 0.0;
	check_for_unexpected_cutpoints( pose );
	for ( Size k = 1; k <= loop_cycles_.size(); ++k ) {
		LoopCycle & loop_cycle = loop_cycles_[ k ];
		LoopClosePotentialEvaluatorCOP potential_evaluator( get_loop_close_potential( pose, loop_cycle, loop_fixed_cost_, use_6D_potential_ ) );
		Real const & loop_closure_energy = potential_evaluator->loop_closure_energy();
		if ( verbose ) TR << TR.Blue << "Cycle " << k << " " << loop_cycle << " ==> " << loop_closure_energy << TR.Reset << std::endl;
		total_energy_ += loop_closure_energy;
		if ( potential_evaluator->involves_current_pose() ) { // save for derivative calculations
			current_pose_loop_score_evaluators_.push_back( potential_evaluator );
		}

	}

	if ( verbose ) TR << TR.Blue << "Total loop close energy  ==> " << total_energy_ << TR.Reset << std::endl;
	debug_assert( utility::isfinite( total_energy_ ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
LoopGraph::update_loops_and_cycles( utility::vector1< Size > const & pose_domain_map,
	utility::vector1< Size > const & cutpoint_open ){

	update_loops( pose_domain_map, cutpoint_open );
	TR.Debug << "NUMBER OF LOOPS " << loops_.size() << std::endl;

	figure_out_loop_cycles();
	TR.Debug << "NUMBER OF CYCLES " << loop_cycles_.size() << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// let's figure out cycles, for which we can use GaussianChainFunc term developed by Rhiju.
// pose_domains are vertices.
// loops are edges (directed!)
void
LoopGraph::figure_out_loop_cycles(){
	figure_out_loop_cycles_tiernan();
	// figure_out_loop_cycles_legacy();
	check_loop_cycles_are_disjoint();
}

// stolen from http://www.boost.org/doc/libs/1_55_0/libs/graph/example/tiernan_print_cycles.cpp
using namespace std;
using namespace boost;

struct cycle_printer
{
	explicit cycle_printer( utility::vector1< utility::vector1< Size > > & cycles ):
		cycles_( cycles )
	{
	}

	template <typename Path, typename Graph>
	void cycle(const Path& p, const Graph& g)
	{
		// Get the property map containing the vertex indices
		// so we can print them.
		typedef typename property_map<Graph, vertex_index_t>::const_type IndexMap;
		IndexMap indices = get(vertex_index, g);

		// Iterate over path printing each vertex that forms the cycle.
		typename Path::const_iterator i, end = p.end();
		utility::vector1< Size > cycle;
		for ( i = p.begin(); i != end; ++i ) {
			// add 1, to convert from 0-indexed vertices to 1-indexed domains
			cycle.push_back( get(indices, *i)+1 );
		}
		cycles_.push_back( cycle );
	}

	utility::vector1< utility::vector1< Size > > & cycles_;
};

//////////////////////////////////////////////////////////////////////////////////////////////
// using canned function in boost library. there is apparently a more efficient one
// by donaldson, but not available in library.
void
LoopGraph::figure_out_loop_cycles_tiernan() {
	loop_cycles_.clear();

	using Graph = boost::directed_graph<>;
	Size const num_domains = loops_from_domain_.size();
	Graph g( num_domains );
	for ( Size n = 1; n <= loops_.size(); n++ ) {
		Size const domain1 = loops_[n].takeoff_domain();
		if ( domain1 == 0 ) continue; // terminal loop
		Size const domain2 = loops_[n].landing_domain();
		if ( domain2 == 0 ) continue; // terminal loop
		// subtract 1, to convert to 0-indexed vertices.
		add_edge( vertex( domain1-1, g ), vertex( domain2-1,g ), g  );
		// internal loops won't be caught by elementary cycles decomposition below -- save now
		if ( domain1 == domain2 ) loop_cycles_.push_back( utility::tools::make_vector1( loops_[n] ) );
	}

	// Instantiate the visitor for saving cycles
	utility::vector1< utility::vector1< Size > > elementary_cycles;
	cycle_printer vis( elementary_cycles );

	// Use the Tiernan algorithm to visit all cycles, printing them
	// as they are found.
	tiernan_all_cycles(g, vis);
	// TR << "Elementary cycles: " << elementary_cycles << std::endl;

	// The cycles above are in terms of domain numbers.
	// For calculating loop scores, we need to work out cycles in terms of loops.
	// Set them up via a recursion, since each domain edge may actually
	// correspond to multiple connections.
	for ( Size n = 1; n <= elementary_cycles.size(); n++ ) {
		utility::vector1< Loop > loops_for_cycle;
		record_loop_cycle( elementary_cycles[ n ], 1, loops_for_cycle );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
LoopGraph::record_loop_cycle(
	utility::vector1< Size > const & elementary_cycle,
	Size const & idx,
	utility::vector1< Loop > const & loops_for_cycle_in )
{
	Size const & current_domain = elementary_cycle[ idx ];
	Size const next_idx = ( idx < elementary_cycle.size() ) ? (idx + 1) : 1;
	Size const & next_domain = elementary_cycle[ next_idx ];

	utility::vector1< Size > const & loops_from_current_domain = loops_from_domain_[ current_domain ];
	bool found_loop( false );
	for ( Size const loop_idx : loops_from_current_domain ) {
		Loop const & loop = loops_[ loop_idx ];
		if ( loop.landing_domain() != next_domain ) continue;

		found_loop = true;

		utility::vector1< Loop > loops_for_cycle = loops_for_cycle_in;
		loops_for_cycle.push_back( loop );

		if ( idx == elementary_cycle.size() ) { // done with cycle
			loop_cycles_.push_back( LoopCycle( loops_for_cycle ) );
		} else {
			runtime_assert( idx < elementary_cycle.size() );
			record_loop_cycle( elementary_cycle, idx + 1, loops_for_cycle );
		}
	}
	// verbiage for debugging -- DEPRECATE in 2016 if not in use.
	if ( !found_loop ) {
		TR << "All loops " << std::endl;
		TR << loops_ << std::endl;
		TR << "loops_from_current_domain (by idx)" << loops_from_current_domain << std::endl;
		TR << "checking elementary_cycle" << elementary_cycle << " at idx " << idx << " looking for loop from current_domain " << current_domain << " to next_domain " << next_domain << "  [ next_idx: " << next_idx << " ] " << std::endl;
	}
	runtime_assert( found_loop );
}

// does not handle cases with complex cycles (multiple cycles sharing same loop)
// DEPRECATE! remove from code in 2016 if boost replacement (tiernan algorithm) does the trick.
void
LoopGraph::figure_out_loop_cycles_legacy()
{
	loop_cycles_.clear();
	Size const num_domains = loops_from_domain_.size();
	if ( num_domains == 0 ) return; // no domains means there's no graph...

	TR.Debug << "NUMBER OF DOMAINS " << num_domains << std::endl;

	// these are used for consistency checks.
	loop_visited_   = utility::vector1< bool >( loops_.size(), false );
	domain_visited_ = utility::vector1< bool >( num_domains, false );

	// there may be several unconnected subgraphs of the graph of domains and loop interconnections.
	bool all_domains_visited( false );
	while ( !all_domains_visited ) {
		all_domains_visited = true;
		for ( Size k = 1; k <= num_domains; k++ ) {
			TR.Debug << "Checking visited: domain " << k << " " << domain_visited_[ k ] << std::endl;
			if ( !domain_visited_[ k ] ) {
				utility::vector1< Size > parent_domains;
				utility::vector1< Loop > loops_so_far;
				look_for_cycles_recursively( k, parent_domains, loops_so_far );
				all_domains_visited = false;
				break; // do another scan for unvisited domains.
			}
		}
	}

}
//////////////////////////////////////////////////////////////////////////////////////////////
// does not handle cases with complex cycles (multiple cycles sharing same loop)
// remove from code in 2016 if boost replacement (tiernan algorithm) does the trick.
void
LoopGraph::look_for_cycles_recursively( Size const current_domain,
	utility::vector1< Size > const parent_domains_in,
	utility::vector1< Loop > const loops_so_far_in ){

	TR.Debug << "Entering domain " << current_domain << std::endl;
	runtime_assert( !domain_visited_[ current_domain ] );
	domain_visited_[ current_domain ] = true;

	utility::vector1< Size > parent_domains = parent_domains_in;
	parent_domains.push_back( current_domain );

	utility::vector1< Size > const & loops_from_current_domain = loops_from_domain_[ current_domain ];
	TR.Debug << " num loops from domain " << current_domain << ": " << loops_from_current_domain.size() << std::endl;
	for ( Size const loop_idx : loops_from_current_domain ) {
		Loop const & loop = loops_[ loop_idx ];
		runtime_assert( !loop_visited_[ loop_idx ] );
		loop_visited_[ loop_idx ] = true;

		TR.Debug << "Going into loop " << loop_idx << std::endl;

		utility::vector1< Loop > loops_so_far = loops_so_far_in;
		loops_so_far.push_back( loop );

		runtime_assert( loop.takeoff_domain() == current_domain );
		Size const & next_domain = loop.landing_domain();

		if ( next_domain  == 0 ) continue; // terminal.

		if ( parent_domains.has_value( loop.landing_domain() ) ) {

			// found a cycle!
			Size const cycle_idx = parent_domains.index( loop.landing_domain() );
			utility::vector1< Loop > loops_for_cycle;
			for ( Size n = cycle_idx; n <= loops_so_far.size(); n++ ) loops_for_cycle.push_back( loops_so_far[ n ] );
			loop_cycles_.push_back( LoopCycle( loops_for_cycle ) );
			// don't keep traversing the graph -- we'd actually get an infinite cycle!

		} else { // keep traversing the graph
			if ( domain_visited_[ next_domain ] ) continue; // may have visited on a previous pass through subgraphs.
			look_for_cycles_recursively( next_domain, parent_domains, loops_so_far );
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
LoopGraph::check_loop_cycles_are_disjoint( bool const verbose /* = false */) {

	has_just_simple_cycles_ = true;

	if ( verbose ) {
		TR << "Looking through loop cycles " << std::endl;
		for ( Size i = 1; i <= loop_cycles_.size(); i++ ) {
			TR << "CYCLE " << i << loop_cycles_[ i ].loops() << std::endl;
		}
	}

	// check that cycles are disjoint from each other...
	for ( Size i = 1; i <= loop_cycles_.size(); i++ ) {
		for ( Size j = i+1; j <= loop_cycles_.size(); j++ ) {
			if ( !check_disjoint( loop_cycles_[i], loop_cycles_[j] ) ) {
				if ( verbose ) TR << "NOT DISJOINT! CYCLE " << i << " CYCLE " << j << std::endl;
				has_just_simple_cycles_ = false;
				return;
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////
bool
LoopGraph::check_disjoint( LoopCycle loop_cycle1, LoopCycle loop_cycle2 ) const {

	for ( Size i = 1; i <= loop_cycle1.size(); i++ ) {
		for ( Size j = 1; j <= loop_cycle2.size(); j++ ) {
			// Don't change to != - that's not defined
			if ( ! ( loop_cycle1.loop(i) == loop_cycle2.loop( j ) ) ) continue;

			if ( error_out_on_complex_cycles_ ) {
				TR << "loop # " << loop_cycle1.loop(i) << " shared between different cycles: "  << std::endl;
				TR << "Cycle1: " << std::endl;
				TR << loop_cycle1 << std::endl;
				TR << "Cycle2: " << std::endl;
				TR << loop_cycle2 << std::endl;
				utility_exit_with_message( "The term loop_close cannot handle multiloops beyond simple cycles! You can run again with approximate handling of nested cycles using the flag: -allow_complex_loop_graph." );
			}
			return false;
		}
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
LoopGraph::update_loops( pose::Pose const & pose ){
	using namespace pose::full_model_info;
	utility::vector1< Size > const pose_domain_map = figure_out_pose_domain_map_const( pose );
	utility::vector1< Size > const & cutpoint_open = const_full_model_info( pose ).cutpoint_open_in_full_model();
	update_loops( pose_domain_map, cutpoint_open );
}


//////////////////////////////////////////////////////////////////////////////////////////////
void
LoopGraph::update_loops( utility::vector1< Size > const & pose_domain_map,
	utility::vector1< Size > const & cutpoint_open ){
	loops_.clear();


	Size takeoff_pos( 0 ), takeoff_domain( 0 ), landing_pos( 0 ), landing_domain( 0 );
	Size nres = pose_domain_map.size();

	for ( Size n = 1; n <= nres; n++ ) {

		bool const at_cutpoint = ( n == nres || cutpoint_open.has_value( n ) );
		bool const at_domain_boundary = ( ( n < nres ) && pose_domain_map[ n ] != pose_domain_map[ n+1 ]);

		if ( at_cutpoint || at_domain_boundary ) {
			landing_pos = n + 1;
			landing_domain = at_cutpoint ? 0 : pose_domain_map[ landing_pos ];
			if ( pose_domain_map[ n ] == 0 ) {
				// loops with size of one residue or greater:
				loops_.push_back( Loop( takeoff_pos, landing_pos, takeoff_domain, landing_domain ) );
			} else if ( ( n < nres ) && pose_domain_map[ n+1 ] > 0 && !at_cutpoint ) {
				// at boundary between fixed domains -- counts as loop
				// sort of weird -- do an early update of takeoff_pos & takeoff_domain.
				takeoff_pos = n;
				takeoff_domain = pose_domain_map[ n ];
				loops_.push_back( Loop( takeoff_pos, landing_pos, takeoff_domain, landing_domain ) );
			}
			takeoff_pos = n;
			takeoff_domain = at_cutpoint ? 0 : pose_domain_map[ n ];
		}
	}

	Size num_domains( 0 );
	for ( Size n = 1; n <= nres; n++ ) {
		if ( n >= 1 && n <= nres && num_domains < pose_domain_map[ n ] ) num_domains = pose_domain_map[ n ];
	}

	loops_from_domain_.clear();
	utility::vector1< Size > blank_vector;
	for ( Size k = 1; k <= num_domains; k++ ) {
		loops_from_domain_[ k ] = blank_vector; // empty
	}

	for ( Size k = 1; k <= loops_.size(); k++ ) {
		Size const & takeoff_domain = loops_[ k ].takeoff_domain();
		if ( takeoff_domain == 0 ) continue;
		loops_from_domain_[ takeoff_domain ].push_back( k );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
// just a consistency check.
void
LoopGraph::check_for_unexpected_cutpoints( pose::Pose const & pose ) const {
	using namespace pose::full_model_info;

	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
	for ( Size n = 1; n < pose.size(); n++ ) {
		if ( pose.fold_tree().is_cutpoint( n ) && ( res_list[ n+1 ] != res_list[n] + 1 ) ) {
			// better not be a closed chainbreak!
			//TR << "Unexpected nonlocal cutpoint-style connection, possibly O2' based: " << res_list[n] << " " << res_list[n+1] << std::endl;
			runtime_assert( ! pose.residue( n   ).has_variant_type( chemical::CUTPOINT_LOWER )
				&& ! pose.residue( n   ).has_variant_type( chemical::C2_BRANCH_POINT ) );
			runtime_assert( ! pose.residue( n+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) );
		}
	}
}

/////////////////////////////////////////////
LoopClosePotentialEvaluatorCOP
LoopGraph::loop_score_evaluator( Size const n ) const {
	return current_pose_loop_score_evaluators_[ n ];
}

//////////////////////////////////////////////////////////////////
// metric for # missing residues + # missing 'connections' (as occur between
//  contiguous fixed domains). Was used historically to gauge how far
//  a stepwise run was getting towards filling in a full model.
//////////////////////////////////////////////////////////////////
Size
LoopGraph::nmissing( pose::Pose const & pose ) const {

	utility::vector1< Size > const & input_domain_map = pose::full_model_info::const_full_model_info( pose ).input_domain_map();
	utility::vector1< Size > const & working_res = pose::full_model_info::const_full_model_info( pose ).working_res();

	Size nmissing_total( 0 );

	for ( Loop const & loop : loops_ ) {
		for ( Size q = loop.takeoff_pos() + 1; q < loop.landing_pos(); q++ ) {
			if ( input_domain_map[ q ] == 0 && working_res.has_value( q ) /*loop residue*/ ) {
				nmissing_total++;
			}
		}

		Size const nmissing_loop = loop.landing_pos() - loop.takeoff_pos() - 1;
		if ( nmissing_loop == 0 ) {
			// super-special -- connection between separate domains.
			// historically used for, e.g., four-way junctions with predefined helices.
			// nmissing was incremented if the helices hadn't been connected yet.
			Size const n = loop.takeoff_pos();
			runtime_assert( loop.landing_pos() == n+1 );
			if ( input_domain_map.size() > 0 && /* needs to be specified by user */
					input_domain_map[ n   ] > 0 &&
					input_domain_map[ n+1 ] > 0 &&
					input_domain_map[ n+1 ] != input_domain_map[ n ] ) {
				if ( loop.takeoff_domain() != loop.landing_domain() ) {
					nmissing_total++;
				}
			}
		}
	}
	return nmissing_total;
}

//////////////////////////////////////////////////////////////////
// actual identities (a, c, g, u, etc.) of missing residues.
//////////////////////////////////////////////////////////////////
utility::vector1< char >
LoopGraph::missing_residues( pose::Pose const & pose ) const{
	utility::vector1< Size > missing_pos;
	utility::vector1< Size > const & input_domain_map = pose::full_model_info::const_full_model_info( pose ).input_domain_map();
	utility::vector1< Size > const & working_res = pose::full_model_info::const_full_model_info( pose ).working_res();
	for ( Loop const & loop : loops_ ) {
		for ( Size k = loop.takeoff_pos() + 1; k < loop.landing_pos(); k++ ) {
			if ( input_domain_map[ k ] == 0 && working_res.has_value( k ) ) missing_pos.push_back( k );
		}
	}
	std::sort( missing_pos.begin(), missing_pos.end() );

	std::string const & full_sequence = pose::full_model_info::const_full_model_info( pose ).full_sequence();
	utility::vector1< char > missing_residues;
	for ( Size n = 1; n <= missing_pos.size(); n++ ) {
		missing_residues.push_back( full_sequence[ missing_pos[n] - 1 ] );
	}
	return missing_residues;
}

//////////////////////////////////////////////////////////////////
/// @details Returns a vector of loop_suites
/// @brief   If include_free_loops is turned off, loops with
/// takeoff_domain == 0 or landing_domain == 0 will not be included
utility::vector1< utility::vector1< Size > >
LoopGraph::loop_suites( bool include_free_loops /* = true */ ) const {
	utility::vector1< utility::vector1< Size > > loop_suites;
	for ( Loop const & loop : loops_ ) {
		if ( !include_free_loops && ( loop.takeoff_domain() == 0 || loop.landing_domain() == 0 ) ) continue;
		utility::vector1< Size > loop_suite_set;
		for ( Size k = loop.takeoff_pos()+1 ; k < loop.landing_pos(); k++ ) loop_suite_set.push_back( k );
		loop_suites.push_back( loop_suite_set );
	}
	return loop_suites;
}


} //loop_graph
} //scoring
} //core
