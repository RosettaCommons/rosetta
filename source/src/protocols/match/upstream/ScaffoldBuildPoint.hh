// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/ScaffoldBuildPoint.hh
/// @brief  Class declarations for the launch point geometry on the Scaffold.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_upstream_ScaffoldBuildPoint_hh
#define INCLUDED_protocols_match_upstream_ScaffoldBuildPoint_hh

// Unit headers
#include <protocols/match/upstream/ScaffoldBuildPoint.fwd.hh>

// Package headers
// AUTO-REMOVED #include <protocols/match/Hit.hh>
#include <protocols/match/upstream/OriginalScaffoldBuildPoint.fwd.hh>
#include <protocols/match/upstream/UpstreamBuilder.fwd.hh>

// Project headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers
// AUTO-REMOVED #include <numeric/xyzVector.hh>

#include <protocols/match/Hit.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace upstream {

class ScaffoldBuildPoint : public utility::pointer::ReferenceCount
{
public:
	typedef core::Size Size;

public:
	ScaffoldBuildPoint();
	ScaffoldBuildPoint( Size index );
	virtual ~ScaffoldBuildPoint();

	virtual bool compatible( ScaffoldBuildPoint const &, bool first_dispatch = true ) const;

	virtual bool compatible( OriginalBackboneBuildPoint const &, bool first_dispatch = true ) const;

	/// @brief Inform the calling function where in the original scaffold
	/// this build point should be inserted.  If the output pose from a matching
	/// has a different number of residues than the original scaffold, then
	/// the calling function must determine where the hit from this build point
	/// should be inserted.
	virtual Size original_insertion_point() const = 0;

	virtual
	void
	insert(
		Size seqpos_to_insert_at,
		Hit const & hit,
		UpstreamBuilderCOP builder,
		core::pose::Pose & pose
	) const = 0;


	Size index() const { return index_; }
	void index( Size setting );

private:
	Size index_;

};


}
}
}

#endif
