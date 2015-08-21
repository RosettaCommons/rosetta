// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/InvrotTreeNodeBase.hh
/// @brief  Forward declaration for inverse rotamer tree node base
/// @author Florian Richter, flosopher@gmail.com, mar 2012

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_AllowedSeqposForGeomCst_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_AllowedSeqposForGeomCst_hh

/// unit headeers
#include <protocols/toolbox/match_enzdes_util/AllowedSeqposForGeomCst.fwd.hh>

//project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

//c++ heades

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

/// @brief a simple helper class that holds a list of what sequence positions
/// each geom cst is allowd to be at
/// not sure about the ideal home of this class yet, the matcher task
/// could use it too
class AllowedSeqposForGeomCst : public utility::pointer::ReferenceCount {

public:
	typedef core::Size Size;

	AllowedSeqposForGeomCst( utility::vector1< utility::vector1< Size > > const & seqpos_for_geomcst );

	AllowedSeqposForGeomCst();

	virtual ~AllowedSeqposForGeomCst();

	core::Size
	num_seqpos_lists() const {
		return seqpos_for_geomcst_.size(); }

	utility::vector1< Size > const &
	seqpos_for_geomcst( Size geomcst ) const;

	/// @brief this function used to live in the matcher task
	/// pose can be passed in optionally to support the ALL tag
	/// in the pos file
	void
	initialize_from_command_line( core::pose::PoseCOP pose = NULL );

private:
	//dimension of this vector is num_geomcst
	utility::vector1< utility::vector1< Size > > seqpos_for_geomcst_;
};


}
}
}

#endif
