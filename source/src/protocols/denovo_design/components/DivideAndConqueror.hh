// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/components/DivideAndConqueror.hh
/// @brief Splits a denovo structure into pieces, and devises a strategy for folding the structure piece-by-piece
/// @author Tom Linsky (tlinsky@uw.edu)


#ifndef INCLUDED_protocols_denovo_design_components_DivideAndConqueror_hh
#define INCLUDED_protocols_denovo_design_components_DivideAndConqueror_hh

#include <protocols/denovo_design/components/DivideAndConqueror.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/denovo_design/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

namespace protocols {
namespace denovo_design {
namespace components {

class SegmentPairSet;
class StructureSlice;
typedef utility::vector1< StructureSlice > StructureSlices;
typedef utility::vector1< StructureSlices > FoldingSolutions;

class StructureSlice : public std::pair< SegmentNameList::const_iterator, SegmentNameList::const_iterator > {
public:
	StructureSlice( SegmentNameList::const_iterator const & start, SegmentNameList::const_iterator const & stop ):
		std::pair< SegmentNameList::const_iterator, SegmentNameList::const_iterator >( start, stop ) {}
};

class BuildPhases : public utility::vector1< SegmentNames > {
public:
	BuildPhases( StructureSlices const & slices );

	BuildPhases( BuildPhases::const_iterator const & begin_it, BuildPhases::const_iterator const & end_it );

	virtual
	~BuildPhases() {};

private:
	BuildPhases() {};
};

typedef std::pair< SegmentName, SegmentName > SegmentPair;

class SegmentPairSet : public std::set< SegmentPair > {
public:
	SegmentPairSet( components::StructureData const & sd );
	SegmentNameSet
	segment_name_set() const;

	friend std::ostream &
	operator<<( std::ostream & os, SegmentPairSet const & pairs );

private:
	SegmentPairSet() {};
};

///@brief Splits a denovo structure into pieces, and devises a strategy for folding the structure piece-by-piece
class DivideAndConqueror : public utility::pointer::ReferenceCount {
public:
	DivideAndConqueror();

	virtual ~DivideAndConqueror();

	DivideAndConquerorOP
	clone() const;

	BuildPhases
	divide_and_conquer( StructureData const & sd ) const;

	void
	set_start_segments( SegmentNameSet const & start_seg_names );

	void
	set_stop_segments( SegmentNameSet const & stop_seg_names );

private:
	/// @brief compute all possible segments
	FoldingSolutions
	compute_all_units(
		StructureData const & sd,
		SegmentNameSet const & paired_segments,
		SegmentNameList::const_iterator start_seg ) const;

	FoldingSolutions
	all_permutations(
		StructureSlices const & segments,
		StructureSlices::const_iterator chosen_item ) const;

	FoldingSolutions::const_iterator
	select_best_solution( utility::vector1< FoldingSolutions::const_iterator > const & valid ) const;

private:
	SegmentNameSet start_segments_;
	SegmentNameSet stop_segments_;
};

class SolutionPredicate {
public:
	SolutionPredicate(
		StructureData const & sd,
		SegmentPairSet const & pairs,
		SegmentNameSet const & start_seg,
		SegmentNameSet const & stop_seg	);

	bool
	operator()( StructureSlices const & slices ) const;

private:
	bool
	check_basic_connectivity( StructureSlices const & slices ) const;

	bool
	check_for_unpaired_segments( StructureSlices const & slices ) const;

	bool
	check_start_segment( StructureSlices const & slices ) const;

	bool
	check_stop_segment( StructureSlices const & slices ) const;

	bool
	check_template_poses( StructureSlices const & slices ) const;

	core::Size
	count_visited_pairs_with_segment(
		SegmentName const & segment_name,
		SegmentNameSet const & visited ) const;

private:
	SegmentPairSet pairs_;
	SegmentNameSet start_segments_;
	SegmentNameSet stop_segments_;
	StructureDataCOP sd_;
};

} //protocols
} //denovo_design
} //components

#endif //INCLUDED_protocols_denovo_design_components_DivideAndConqueror_hh

