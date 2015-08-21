// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_kinematic_closure_pivot_pickers_PivotPicker_HH
#define INCLUDED_protocols_kinematic_closure_pivot_pickers_PivotPicker_HH

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <boost/noncopyable.hpp>

namespace protocols {
namespace kinematic_closure {
namespace pivot_pickers {

/// @brief Base class for all the pivot picking algorithms.
class PivotPicker
	: public utility::pointer::ReferenceCount, private boost::noncopyable {

public:

	/// @brief Return the name of this pivot picker.
	virtual std::string get_name() const = 0;

	/// @brief Return a loop object.  The pivots will be taken to be the start,
	/// cut and stop residues.
	virtual Loop pick(Pose const & pose, Loop const & loop) = 0;

};

}
}
}

#endif

