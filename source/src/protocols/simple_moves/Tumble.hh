// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   Tumble.hh
///
/// @brief
/// @author Yifan Song

#ifndef INCLUDED_protocols_docking_Tumble_hh
#define INCLUDED_protocols_docking_Tumble_hh

// Unit Headers
//#include <protocols/docking/Tumble.fwd.hh>

// Package Headers

// Project Headers
#include <core/types.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <protocols/toolbox/task_operations/InterfaceTaskOperation.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <protocols/moves/Mover.hh>

#include <utility/tag/Tag.fwd.hh>


// Utility Headers
#include <utility/vector1.hh>

// Numeric Headers

// ObjexxFCL Headers

// C++ headers

namespace protocols {
namespace simple_moves {

class Tumble: public protocols::moves::Mover
{
public:
	/// @brief
	/// 	empty constructor fills values with the expected defaults
	Tumble();

	Tumble(
		Size const rb_jump_in,
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionOP(NULL)
	);

	//destructor
	~Tumble();

    numeric::xyzVector<core::Real> center_of_mass(core::pose::Pose const & pose);

    void apply( core::pose::Pose & pose );

    void
    parse_my_tag(
                 TagCOP const tag,
                 basic::datacache::DataMap & datamap,
                 Filters_map const & filters,
                 moves::Movers_map const & movers,
                 Pose const & pose
                 );

	virtual protocols::moves::MoverOP clone() const;

	virtual protocols::moves::MoverOP fresh_instance() const;

    std::string get_name() const;


private:
    utility::vector1< core::Size > residue_list_;
};
} // simple_moves
} // protocols

#endif

