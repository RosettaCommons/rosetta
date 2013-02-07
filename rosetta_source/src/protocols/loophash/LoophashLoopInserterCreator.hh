// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file LoophashLoopInserterCreator.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_loophash_LoophashLoopInserterCreator_HH
#define INCLUDED_protocols_loophash_LoophashLoopInserterCreator_HH

// Unit Headers
#include <devel/loop_creation/LoopInserterCreator.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace loophash {

/// @brief creator for the LoophashLoopInserter class
class LoophashLoopInserterCreator : public protocols::loops::loop_creation::LoopInserterCreator
{
public:
	LoophashLoopInserterCreator();
	virtual ~LoophashLoopInserterCreator();

	virtual protocols::loops::loop_creation::LoopInserterOP create_loop_inserter() const;
	virtual std::string inserter_name() const;
};
	
} //protocols
} //loophash

#endif
