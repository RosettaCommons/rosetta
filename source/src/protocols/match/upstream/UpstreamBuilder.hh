// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/upstream/UpstreamBuilder.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_upstream_UpstreamBuilder_hh
#define INCLUDED_protocols_match_upstream_UpstreamBuilder_hh

// Unit headers
#include <protocols/match/upstream/UpstreamBuilder.fwd.hh>

// Package headers
#include <protocols/match/BumpGrid.fwd.hh>
// AUTO-REMOVED #include <protocols/match/Hit.hh>
#include <protocols/match/upstream/ScaffoldBuildPoint.fwd.hh>
#include <protocols/match/upstream/ProteinUpstreamBuilder.fwd.hh>

// Project headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>

#include <protocols/match/Hit.fwd.hh>


namespace protocols {
namespace match {
namespace upstream {

class UpstreamBuilder : public utility::pointer::ReferenceCount {
public:
	typedef core::Size   Size;
	typedef core::Vector Vector;

public:
	virtual ~UpstreamBuilder();

	virtual
	UpstreamBuilderOP
	clone() const = 0;

	/// @brief Iterate across possible conformations for the upstream
	/// half of the hit, and for each (non-colliding) conformation,
	/// invoke build on the downstream algorithm.
	/// Return a list of hits.
	virtual
	std::list< Hit >
	build(
		ScaffoldBuildPoint const & build_point
	) const = 0;

	/// @brief Reconstruct the upstream conformation for a hit and pass that conformation to
	/// an upstream residue processor.
	virtual
	void
	recover_hit(
		Hit const & hit,
		ScaffoldBuildPoint const & build_point,
		UpstreamResidueProcessor & processor
	) const = 0;

	/// @brief Reconstruct the upstream conformation for a set of hits and pass their conformations to
	/// an upstream residue processor.
	virtual
	void
	recover_hits(
		std::list< Hit >::const_iterator hits_begin,
		std::list< Hit >::const_iterator hits_end,
		ScaffoldBuildPoint const & build_point,
		UpstreamResidueProcessor & processor
	) const = 0;

	virtual
	Size
	n_restypes_to_build() const = 0;

	virtual
	core::chemical::ResidueTypeCOP
	restype( Size which_restype ) const = 0;

	virtual bool compatible(
		Hit const & my_hit,
		ScaffoldBuildPoint const & build_point_mine,
		UpstreamBuilder const & other,
		Hit const & other_hit,
		ScaffoldBuildPoint const & build_point_other,
		bool first_dispatch = true
	) const;


	virtual bool compatible(
		Hit const & my_hit,
		ScaffoldBuildPoint const & build_point_mine,
		ProteinUpstreamBuilder const & other,
		Hit const & other_hit,
		ScaffoldBuildPoint const & build_point_other,
		bool first_dispatch = true
	) const;


	void
	set_bb_grid( BumpGridCOP bbgrid );

protected:

	BumpGrid const &
	bbgrid() const {
		return *bbgrid_;
	}

private:
	BumpGridCOP bbgrid_;

};

class UpstreamResidueProcessor : public utility::pointer::ReferenceCount
{
public:
	virtual ~UpstreamResidueProcessor();

	virtual
	void
	process_hit(
		Hit const & hit,
		core::conformation::Residue const & upstream_conformation
	) = 0;

};

}
}
}

#endif
