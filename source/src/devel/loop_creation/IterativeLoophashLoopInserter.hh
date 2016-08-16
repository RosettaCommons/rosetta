// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LoophashIterativeLoophashLoopInserter.hh
///
/// @brief
/// @author Tim Jacobs


#ifndef INCLUDED_protocols_loophash_IterativeLoophashLoopInserter_HH
#define INCLUDED_protocols_loophash_IterativeLoophashLoopInserter_HH

//Unit
#include <devel/loop_creation/LoophashLoopInserter.hh>

//Protocols
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loops/Loop.fwd.hh>


//namespace protocols {
//namespace loophash {
namespace devel {
namespace loop_creation {

class IterativeLoophashLoopInserter : public devel::loop_creation::LoophashLoopInserter
{
public:

	IterativeLoophashLoopInserter();

	protocols::moves::MoverOP
	clone() const;

	protocols::moves::MoverOP
	fresh_instance() const;

	std::string
	get_name() const;

	void
	apply(
		core::pose::Pose & pose
	);

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

private:
	core::Real max_closure_deviation_;
	core::Real max_insertions_;
};

} //loop creation
} //devel
//} //protocols
//} //loophash

#endif
