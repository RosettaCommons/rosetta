// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_moves/ModifyVariantTypeMover.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_simple_moves_ModifyVariantTypeMover_HH
#define INCLUDED_protocols_simple_moves_ModifyVariantTypeMover_HH

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/ModifyVariantTypeMover.fwd.hh>

#include <core/pack/task/TaskFactory.fwd.hh>

// Project headers

#include <core/pose/Pose.fwd.hh>


// ObjexxFCL Headers

// C++ Headers
#include <string>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

/// @brief Adds variant types to selected residues
class ModifyVariantTypeMover : public protocols::moves::Mover
{
public:
	// default constructor (nmoves=1)
	ModifyVariantTypeMover();

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

private:
	core::pack::task::TaskFactoryOP task_factory_;
	utility::vector1<std::string> add_target_types_;
	utility::vector1<std::string> remove_target_types_;
};

} // moves
} // protocols


#endif
