// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/output/DownstreamRMSEvaluator.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_DownstreamRMSEvaluator_hh
#define INCLUDED_protocols_match_output_DownstreamRMSEvaluator_hh

// Unit headers
#include <protocols/match/output/DownstreamRMSEvaluator.fwd.hh>

// Package headers
#include <protocols/match/output/MatchEvaluator.hh>
#include <protocols/match/downstream/DownstreamBuilder.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers

// Numeric headers

#include <core/id/AtomID.fwd.hh>
#include <utility/vector1.hh>

#ifdef PYROSETTA
#include <core/id/AtomID.hh>
#endif

namespace protocols {
namespace match {
namespace output {

/// @brief In the best of all possible worlds, this class would be sufficiently
/// generic such that I could compare RMS for arbitrary subsets of atoms on
/// the downstream partner, but in my first pass implementation, I'm writing
/// while aiming at the RigidLigandBuilder -- 1 residue -- and I'll compare
/// all heavy atoms.
class DownstreamRMSEvaluator : public MatchEvaluator
{
public:
	typedef core::Size       Size;
	typedef core::Vector     Vector;
	typedef core::id::AtomID AtomID;
public:
	DownstreamRMSEvaluator();
	virtual ~DownstreamRMSEvaluator();

	void
	set_downstream_pose( core::pose::PoseCOP dspose );

	void
	set_n_geometric_constraints( Size setting );

	void
	set_downstream_builder(
		Size which_geom_cst,
		downstream::DownstreamBuilderCOP ds_builder
	);

	virtual
	Real
	score( match const & m ) const;

	/// @brief Causes a graceful exit.  The DownstreamRMSEvaluator is incompatible with
	/// the match_dspos1 match definition!
	virtual
	Real
	score( match_dspos1 const & m ) const;

private:
	core::pose::PoseCOP dspose_;

	Size n_geometric_constraints_;

	utility::vector1< AtomID > atoms_to_compare_;
	utility::vector1< downstream::DownstreamBuilderCOP > ds_builders_;

};

}
}
}

#endif
