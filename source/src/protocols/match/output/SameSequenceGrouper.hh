// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/SameSequenceGrouper.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_SameSequenceGrouper_hh
#define INCLUDED_protocols_match_output_SameSequenceGrouper_hh

// Unit headers
#include <protocols/match/output/SameSequenceGrouper.fwd.hh>

// Package headers
#include <protocols/match/downstream/DownstreamBuilder.fwd.hh>
#include <protocols/match/output/MatchGrouper.hh>
#include <protocols/match/output/UpstreamHitCacher.fwd.hh>
#include <protocols/match/Hit.fwd.hh>

//Project headers
#include <core/id/AtomID.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <map>

#include <utility/OrderedTuple.fwd.hh>


namespace protocols {
namespace match {
namespace output {

/// @brief Class to group matches that represent the same amino acids at
/// the same launch points. E.g. Two matches that both put a cys at launch
/// point #33, a ser a launch point #42 and another ser at launch point
/// #80 would be grouped together -- even if they are different rotamers.
class SameSequenceGrouper : public MatchGrouper {
public:
	typedef core::Real Real;
	typedef core::Size Size;
	typedef std::map< utility::OrderedTuple< utility::vector1< Size > >, Size > SequenceMap;

public:
	SameSequenceGrouper();
	SameSequenceGrouper( Size ncst );

	virtual
	~SameSequenceGrouper();

	virtual
	Size
	assign_group_for_match(
		match const & m
	);

	virtual
	Size
	assign_group_for_match(
		match_dspos1 const & m
	);

	virtual
	void
	reset();

	virtual
	void
	set_n_geometric_constraints( Size n_csts );

	void
	set_hit_cacher( UpstreamHitCacherOP cacher );

private:

	Size n_geometric_constraints_;
	UpstreamHitCacherOP hit_cacher_;
	SequenceMap sequence_indexer_;

};


/// @brief class that groups based on same sequence and proximity of the downstream object (based on rms )
/// NOTE: right now only the downstream position according to the first geomcst id upstream residue is
/// taken into account
class SameSequenceAndDSPositionGrouper : public SameSequenceGrouper {

public:
	typedef std::map< std::pair< Size, Size >, Size > SequenceLigPosMap;
	typedef core::Vector Vector;
	typedef SameSequenceGrouper parent;

	SameSequenceAndDSPositionGrouper();
	SameSequenceAndDSPositionGrouper( Size ncst );

	virtual
	~SameSequenceAndDSPositionGrouper();

	//virtual
	//Size
	//assign_group_for_match(
	// match const & m
	//);

	virtual
	Size
	assign_group_for_match(
		match_dspos1 const & m
	);

	virtual
	void
	reset();

	virtual
	void
	set_n_geometric_constraints( Size n_csts );

	void
	set_rms_group_cutoff( Real cutoff );

	void
	set_downstream_builder(
		Size geomcst_id,
		downstream::DownstreamBuilderCOP dsbuilder
	);

	void
	set_relevant_atom_ids(
		utility::vector1< core::id::AtomID > const & relevant_atom_ids
	);

private:

	Size
	assign_downstream_position_group_for_match(
		match_dspos1 const & m
	);

	SequenceLigPosMap sequence_pos_map_;

	Real rms_group_cutoff_;
	utility::vector1< utility::vector1< Vector> > representative_dspos_;
	utility::vector1< core::id::AtomID > relevant_atom_ids_;
	//needed to build the downstream target for rms comparison
	utility::vector1< downstream::DownstreamBuilderCOP > dsbuilders_;
};

}
}
}

#endif
