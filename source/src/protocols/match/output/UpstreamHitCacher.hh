// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/UpstreamHitCacher.hh
/// @brief  Declaration for class to cache and recall the conformations of upstream hits.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_UpstreamHitCacher_hh
#define INCLUDED_protocols_match_output_UpstreamHitCacher_hh

// Unit headers
#include <protocols/match/output/UpstreamHitCacher.fwd.hh>

// Package headers
#include <protocols/match/Hit.fwd.hh>
#include <protocols/match/Matcher.fwd.hh>
#include <protocols/match/upstream/UpstreamBuilder.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#ifdef WIN32
#include <utility/OrderedTuple.hh>
#include <utility/vector1.hh>
#endif

// C++ headers
#include <map>

#include <utility/OrderedTuple.fwd.hh>
#include <utility/vector1_bool.hh>

#ifdef PYROSETTA
#include <protocols/match/Matcher.hh>
#include <utility/OrderedTuple.hh>
#endif


namespace protocols {
namespace match {
namespace output {

class UpstreamHitCacher : public upstream::UpstreamResidueProcessor {
public:

	typedef utility::fixedsizearray1< Size, 2 >          ScaffoldRotamerPair;
	typedef utility::OrderedTuple< ScaffoldRotamerPair > ScaffoldRotamerTuple;

public:
	UpstreamHitCacher( MatcherCOP matcher );
	virtual ~UpstreamHitCacher();

	void
	set_cache_size( Size n_rotamers_to_cache );

	core::conformation::ResidueCOP
	upstream_conformation_for_hit( Size geometric_constraint_id, Hit const & hit );

public:

	/// @brief The method by which the UpstreamHitCacher communicates with the UpstreamBuilders.
	/// CAUTION: this class should not be handed to an UpstreamBuilder directly.  It should only
	/// hand itself to an UpstreamBuilder. Ask the UpstreamHitCacher for a particular hit, and it
	/// will call recover_rotamer() on the upstream builder; it's doing bookkeeping in the background
	/// and that bookkeeping is important for the success of this function.
	virtual
	void
	process_hit(
		Hit const & hit,
		core::conformation::Residue const & upstream_conformation
	);


private:

	/// @brief Allocate space in the arrays dependent on n_geometric_constraints_ and
	/// n_confs_to_cache_.
	void
	resize_arrays();

	/// @brief Returns 0 if the scaffold/rotamer pair is not already in the queue, and
	/// the non-zero index in the queue if it is.
	Size
	already_in_queue( Size cst_id, ScaffoldRotamerTuple const & rotid ) const;

	/// @brief Construct the rotamer for the requested scaffold/rotamer pair and
	/// put it into the queue, evicting the previous queue resident if necessary.
	Size
	fetch( Size cst_id, ScaffoldRotamerTuple const & rotid );

private:
	MatcherCOP matcher_;

	Size n_geometric_constraints_;
	Size n_confs_to_cache_;

	utility::vector1< std::map< ScaffoldRotamerTuple, Size > > index_for_rotamer_;

	Size which_cst_being_processed_;

	utility::vector1< Size > queue_head_;
	utility::vector1< utility::vector1< ScaffoldRotamerTuple > >           scafrot_pair_for_conf_;
	utility::vector1< utility::vector1< core::conformation::ResidueCOP > > upstream_confs_;

};

}
}
}

#endif
