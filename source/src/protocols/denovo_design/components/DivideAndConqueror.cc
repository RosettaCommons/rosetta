// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/components/DivideAndConqueror.cc
/// @brief Splits a denovo structure into pieces, and devises a strategy for folding the structure piece-by-piece
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/components/DivideAndConqueror.hh>

// Protocol headers
#include <protocols/denovo_design/architects/StrandArchitect.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/util.hh>

// Basic headers
#include <basic/Tracer.hh>

// Boost headers
#include <boost/assign.hpp>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.components.DivideAndConqueror" );


namespace protocols {
namespace denovo_design {
namespace components {

std::ostream &
operator<<( std::ostream & os, SegmentPairSet const & pairs )
{
	os << "[ ";
	for ( SegmentPairSet::const_iterator p=pairs.begin(); p!=pairs.end(); ++p ) {
		if ( p != pairs.begin() ) os << ", ";
		os << p->first << "<-->" << p->second;
	}
	os << "]";
	return os;
}

DivideAndConqueror::DivideAndConqueror():
	utility::pointer::ReferenceCount(),
	start_segments_(),
	stop_segments_()
{
}

DivideAndConqueror::~DivideAndConqueror(){}

void
DivideAndConqueror::set_start_segments( SegmentNameSet const & start_seg )
{
	start_segments_ = start_seg;
}

void
DivideAndConqueror::set_stop_segments( SegmentNameSet const & stop_seg )
{
	stop_segments_ = stop_seg;
}

std::ostream &
operator<<( std::ostream & os, StructureSlice const & slice )
{
	SegmentNameList const names( slice.first, slice.second );
	os << names;
	return os;
}

BuildPhases
DivideAndConqueror::divide_and_conquer( StructureData const & sd ) const
{
	if ( sd.segments_begin() == sd.segments_end() ) {
		std::stringstream msg;
		msg << "DivideAndConqueror::divide_and_conquer(): No segments are present in the StructureData!"
			<< std::endl << " SD = " << sd << std::endl;
		utility_exit_with_message( msg.str() );
	}
	// determine things that need to be paired
	SegmentPairSet const pairs( sd );
	TR << "Segment pairs are: " << pairs << std::endl;

	// Compute pieces of the structure that bridge the gaps between pairs segments
	FoldingSolutions const unordered_solutions = compute_all_units( sd, pairs.segment_name_set(), sd.segments_begin() );
	TR << "Computed " << unordered_solutions.size() << " structure units with satisfied pairings." << std::endl;
	TR.Debug << unordered_solutions << std::endl;

	// catch case where we are only building one segment
	SegmentNameList::const_iterator next = ++sd.segments_begin();
	if ( next == sd.segments_end() ) {
		StructureSlices const slices = boost::assign::list_of (StructureSlice( sd.segments_begin(), sd.segments_end()));
		return BuildPhases( slices );
	}

	// enumerate all possible build orders
	FoldingSolutions solutions;
	for ( FoldingSolutions::const_iterator solution=unordered_solutions.begin(); solution!=unordered_solutions.end(); ++solution ) {
		FoldingSolutions const sol_list = all_permutations( *solution, solution->end() );
		solutions.insert( solutions.end(), sol_list.begin(), sol_list.end() );
	}

	TR.Debug << "Possible solutions: " << solutions << std::endl;

	// collect valid solutions
	SolutionPredicate const is_valid( sd, pairs, start_segments_, stop_segments_ );
	utility::vector1< FoldingSolutions::const_iterator > valid;
	for ( FoldingSolutions::const_iterator sol=solutions.begin(); sol!=solutions.end(); ++sol ) {
		if ( is_valid( *sol ) ) {
			TR << "Valid: " << *sol << std::endl;
			valid.push_back( sol );
		}
	}

	if ( valid.empty() ) {
		std::stringstream msg;
		msg << "DivideAndConqueror: Could not find any valid build solutions.  Check your inputs and try again" << std::endl;
		utility_exit_with_message( msg.str() );
	}

	FoldingSolutions::const_iterator const selected_valid = select_best_solution( valid );

	TR << "Found " << valid.size() << " build solutions.  Selected "
		<< *selected_valid << std::endl;
	return BuildPhases( *selected_valid );
}

DivideAndConquerorOP
DivideAndConqueror::clone() const
{
	return DivideAndConquerorOP( new DivideAndConqueror( *this ) );
}

FoldingSolutions::const_iterator
DivideAndConqueror::select_best_solution( utility::vector1< FoldingSolutions::const_iterator > const & valid ) const
{
	return *valid.begin();
}

FoldingSolutions
DivideAndConqueror::all_permutations(
	StructureSlices const & slices,
	StructureSlices::const_iterator chosen_item ) const
{
	if ( slices.size() == 1 ) return boost::assign::list_of (slices);

	// make a copy of list without the chosen item
	StructureSlices slices_copy;
	for ( StructureSlices::const_iterator slice=slices.begin(); slice!=slices.end(); ++slice ) {
		if ( ( chosen_item != slices.end() ) && ( slice == chosen_item ) ) continue;
		slices_copy.push_back( *slice );
	}

	FoldingSolutions retval;
	for ( StructureSlices::const_iterator slice=slices_copy.begin(); slice!=slices_copy.end(); ++slice ) {
		// find all permutations of copy list
		FoldingSolutions const permutations = all_permutations( slices_copy, slice );

		// Add x to the beginning of each permutation of copy list
		for ( FoldingSolutions::const_iterator p=permutations.begin(); p!=permutations.end(); ++p ) {
			StructureSlices new_list;
			if ( chosen_item != slices.end() ) new_list.push_back( *chosen_item );
			for ( StructureSlices::const_iterator t=p->begin(); t!=p->end(); ++t ) {
				new_list.push_back( *t );
			}
			retval.push_back( new_list );
		}
	}
	return retval;
}

class StructureSlicePredicate {
public:
	StructureSlicePredicate() {};

	bool
	operator()( StructureData const & sd, StructureSlice const & slice ) const;
};

bool
StructureSlicePredicate::operator()(
	StructureData const & sd,
	StructureSlice const & slice ) const
{
	// start iterator == stop iterator
	if ( slice.first == slice.second ) return false;

	// bad slice if all segments have template
	for ( SegmentNameList::const_iterator s=slice.first; s!=slice.second; ++s ) {
		if ( ! sd.segment( *s ).has_template_pose() ) return true;
	}
	return false;
}

/// @brief compute all possible slices
FoldingSolutions
DivideAndConqueror::compute_all_units(
	StructureData const & sd,
	SegmentNameSet const & paired_segments,
	SegmentNameList::const_iterator start_seg ) const
{
	SegmentNameList::const_iterator next_segment = start_seg;
	++next_segment;
	if ( next_segment == sd.segments_end() ) return FoldingSolutions();

	// Create list of slices starting at start_seg at which we might stop building
	StructureSlicePredicate const slice_ok;

	StructureSlices slices;
	for ( SegmentNameList::const_iterator stop_seg=start_seg; stop_seg!=sd.segments_end(); ++stop_seg ) {
		if ( stop_seg == start_seg ) continue;
		if ( paired_segments.find( *stop_seg ) == paired_segments.end() ) continue;
		SegmentNameList::const_iterator next = stop_seg;
		if ( next != sd.segments_end() ) ++next;
		StructureSlice const newslice( start_seg, next );
		if ( slice_ok( sd, newslice ) ) {
			slices.push_back( newslice );
		}
	}
	if ( slices.empty() || ( slices.rbegin()->second != sd.segments_end() ) ) {
		slices.push_back( StructureSlice( start_seg, sd.segments_end() ) );
	}

	FoldingSolutions retval;
	for ( StructureSlices::const_iterator slice=slices.begin(); slice!=slices.end(); ++slice ) {
		SegmentNameList::const_iterator next_start_segment = slice->second;
		--next_start_segment;
		FoldingSolutions solutions_from_me = compute_all_units( sd, paired_segments, next_start_segment );
		if ( solutions_from_me.empty() ) {
			retval.push_back( boost::assign::list_of (*slice) );
		} else {
			for ( FoldingSolutions::const_iterator sol=solutions_from_me.begin(); sol!=solutions_from_me.end(); ++sol ) {
				StructureSlices my_solution = boost::assign::list_of (*slice);
				my_solution.insert( my_solution.end(), sol->begin(), sol->end() );
				retval.push_back( my_solution );
			}
		}
	}

	TR.Debug << "Units starting with segment " << *start_seg << ": " << retval << std::endl;
	return retval;
}


BuildPhases::BuildPhases( StructureSlices const & slices ):
	utility::vector1< SegmentNames >()
{
	for ( StructureSlices::const_iterator slice=slices.begin(); slice!=slices.end(); ++slice ) {
		push_back( SegmentNames( slice->first, slice->second ) );
	}
}

BuildPhases::BuildPhases( BuildPhases::const_iterator const & begin_it, BuildPhases::const_iterator const & end_it ):
	utility::vector1< SegmentNames >( begin_it, end_it )
{}

SegmentPairSet::SegmentPairSet( components::StructureData const & sd )
{
	for ( SegmentPairingCOPs::const_iterator p=sd.pairings_begin(); p!=sd.pairings_end(); ++p ) {
		if ( (*p)->type() != SegmentPairing::STRAND ) continue;
		SegmentPair const pair = std::make_pair( *(*p)->segments().begin(), *( (*p)->segments().begin() + 1 ) );
		insert( pair );
	}
}

SegmentNameSet
SegmentPairSet::segment_name_set() const
{
	SegmentNameSet visited;
	for ( SegmentPairSet::const_iterator pair=begin(); pair!=end(); ++pair ) {
		visited.insert( pair->first );
		visited.insert( pair->second );
	}
	return visited;
}

SolutionPredicate::SolutionPredicate(
	StructureData const & sd,
	SegmentPairSet const & pairs,
	SegmentNameSet const & start_seg,
	SegmentNameSet const & stop_seg ):
	pairs_( pairs ),
	start_segments_( start_seg ),
	stop_segments_( stop_seg ),
	sd_( StructureDataOP( new StructureData( sd ) ) )
{}

bool
SolutionPredicate::operator()( StructureSlices const & slices ) const
{
	if ( !check_basic_connectivity( slices ) ) return false;
	if ( !check_template_poses( slices ) ) return false;
	if ( !check_for_unpaired_segments( slices ) ) return false;
	if ( !check_start_segment( slices ) ) return false;
	if ( !check_stop_segment( slices ) ) return false;
	return true;
}

bool
SolutionPredicate::check_start_segment( StructureSlices const & slices ) const
{
	if ( start_segments_.empty() ) return true;

	// start segment must be present in first StructureSegment
	StructureSlice const & seg = *slices.begin();
	SegmentNameSet phase1_segments( seg.first, seg.second );

	SegmentNames const diff =
		protocols::denovo_design::set_difference< SegmentNames, SegmentName >( start_segments_.begin(), start_segments_.end(), phase1_segments.begin(), phase1_segments.end() );

	TR.Debug << "Difference: " << diff << std::endl;
	if ( !diff.empty() ) {
		TR.Debug << slices << " is bad because the first phase does not contain the user-set start_segments ("
			<< start_segments_ << ")" << std::endl;
		return false;
	}
	return true;
}

bool
SolutionPredicate::check_stop_segment( StructureSlices const & slices ) const
{
	if ( stop_segments_.empty() ) return true;

	// start segment must be present in first StructureSegment
	StructureSlice const & seg = *slices.rbegin();
	SegmentNameSet phaseN_segments( seg.first, seg.second );

	SegmentNames const diff =
		protocols::denovo_design::set_difference< SegmentNames, SegmentName >( stop_segments_.begin(), stop_segments_.end(), phaseN_segments.begin(), phaseN_segments.end() );

	TR.Debug << "Difference: " << diff << std::endl;
	if ( !diff.empty() ) {
		TR.Debug << slices << " is bad because the first phase does not contain the user-set stop_segments ("
			<< stop_segments_ << ")" << std::endl;
		return false;
	}
	return true;
}

bool
SolutionPredicate::check_template_poses( StructureSlices const & slices ) const
{
	StructureSlicePredicate const slice_ok;
	for ( StructureSlices::const_iterator slice=slices.begin(); slice!=slices.end(); ++slice ) {
		if ( slice_ok( *sd_, *slice ) ) continue;

		TR.Debug << slices << " contains a slice that doesn't build any segments." << std::endl;
		return false;
	}
	return true;
}

bool
SolutionPredicate::check_basic_connectivity( StructureSlices const & slices ) const
{
	if ( slices.empty() ) {
		TR.Debug << slices << " is bad because it is empty" << std::endl;
		return false;
	}

	SegmentNameSet visited( slices.begin()->first, slices.begin()->second );

	for ( StructureSlices::const_iterator slice=++slices.begin(); slice!=slices.end(); ++slice ) {
		if ( slice->first == slice->second ) {
			TR.Debug << slices << " is bad because it contains an empty build segment" << std::endl;
			return false;
		}
		SegmentNameList::const_iterator last = slice->second;
		--last;
		if ( ( visited.find( *last ) == visited.end() ) && ( visited.find( *slice->first ) == visited.end() ) ) {
			// neither the start nor end of this segment has been visited yet -- bad connectivity!
			TR.Debug << slices << " is bad because it contains an orphaned segment ("
				<< *slice->first << " or " << *last << " visited=" << visited << std::endl;
			return false;
		}
		if ( ( visited.find( *last ) != visited.end() ) && ( visited.find( *slice->first ) != visited.end() ) ) {
			// both slices are visited.  This segment is useless
			TR.Debug << slices << " is bad because both ends have already been visited" << std::endl;
			return false;
		}
		visited.insert( slice->first, slice->second );
	}
	return true;
}

bool
SolutionPredicate::check_for_unpaired_segments( StructureSlices const & slices ) const
{
	SegmentNameSet const in_pair = pairs_.segment_name_set();
	SegmentNameSet visited;

	for ( StructureSlices::const_iterator slice=slices.begin(); slice!=slices.end(); ++slice ) {
		TR << "Looking at slice " << *slice << std::endl;

		// 1. add segement names to visited list
		visited.insert( slice->first, slice->second );

		// 2. all slices in a pair must have at least one adjacent segment in "visited"
		for ( SegmentNameList::const_iterator seg=slice->first; seg!=slice->second; ++seg ) {
			if ( in_pair.find( *seg ) == in_pair.end() ) continue;
			// beyond here, we are in a pair
			// here, we could all pairs this segment is associated with
			core::Size const visited_pairs = count_visited_pairs_with_segment( *seg, visited );
			if ( visited_pairs == 0 ) {
				TR.Debug << "Solution " << slices
					<< " is bad because there it contains a build phase with an orphaned pairing" << std::endl;
				return false;
			}
		}
	}

	return true;
}

core::Size
SolutionPredicate::count_visited_pairs_with_segment(
	SegmentName const & segment_name,
	SegmentNameSet const & visited ) const
{
	core::Size count = 0;
	for ( SegmentPairSet::const_iterator p=pairs_.begin(); p!=pairs_.end(); ++p ) {
		// if this segment has nothing to do with the pairing, skip over it
		if ( ( p->first != segment_name ) && ( p->second != segment_name ) ) continue;

		if ( visited.find( segment_name ) == visited.end() ) continue;

		// Find pairing partner
		SegmentName partner;
		if ( p->first == segment_name ) partner = p->second;
		else partner = p->first;

		if ( visited.find( partner ) == visited.end() ) continue;
		count += 1;
	}
	return count;
}

} //protocols
} //denovo_design
} //components

