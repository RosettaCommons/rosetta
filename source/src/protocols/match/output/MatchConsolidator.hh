// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/MatchConsolidator.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_MatchConsolidator_hh
#define INCLUDED_protocols_match_output_MatchConsolidator_hh

// Unit headers
#include <protocols/match/output/MatchConsolidator.fwd.hh>

// Package headers
#ifdef WIN32
#include <protocols/match/Hit.hh>
#endif
#include <protocols/match/output/MatchProcessor.hh>
#include <protocols/match/output/MatchGrouper.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/heap.hh>

// C++ headers
#include <list>

#include <utility/vector0.hh>
#include <utility/vector1_bool.hh>


namespace protocols {
namespace match {
namespace output {

class MatchConsolidator : public MatchProcessor {
public:
	typedef core::Real Real;
	typedef core::Size Size;
public:
	MatchConsolidator();

	virtual
	~MatchConsolidator();

	virtual
	void
	begin_processing();

	virtual
	void
	process_match(
		match const & m
	);

	virtual
	void
	process_match(
		match_dspos1 const & m
	);

	virtual
	void
	end_processing();

	void
	set_n_to_output_per_group( Size setting );

	void
	set_grouper( MatchGrouperOP grouper );

	void
	reset_grouper();

private:
	void end_processing_of_regular_match_groups();
	void end_processing_of_match_dspos1_groups();

private:

	MatchGrouperOP   grouper_;

	Size n_to_output_per_group_;

	utility::vector1< BestMatchesCollectionOP > match_groups_;

};


class BestMatchesCollection : public utility::pointer::ReferenceCount
{
public:
	typedef core::Size Size;
	typedef core::Real Real;

public:
	BestMatchesCollection( Size n_to_keep, bool dspos1_mode = false );
	virtual ~BestMatchesCollection();

	void
	add_match( match const & m, Real score );

	void
	add_match_dspos1( match_dspos1 const & m, Real score );

	bool dspos1_mode() const { return dspos1_mode_; }

	Size
	n_kept_matches() const;

	match const & 
	kept_match( Size which_match ) const;

	match_dspos1 const & 
	kept_match_dspos1( Size which_match ) const;

private:

	Size
	index_for_new_match( Real score );

private:
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// Size const n_top_matches_to_keep_;
	bool const dspos1_mode_; // store matches or match_dspos1s?


	utility::heap scores_heap_;

	utility::vector1< match > best_matches_;
	utility::vector1< match_dspos1 > best_match_dspos1s_;

};


}
}
}

#endif
