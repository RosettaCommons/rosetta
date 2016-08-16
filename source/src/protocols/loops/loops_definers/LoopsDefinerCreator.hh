// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loops_definers/LoopsDefinerCreator.hh
/// @brief  Base class LoopsDefinerCreator for the load-time factory registration scheme
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_loops_loops_definers_LoopsDefinerCreator_hh
#define INCLUDED_protocols_loops_loops_definers_LoopsDefinerCreator_hh

// Unit Headers
#include <protocols/loops/loops_definers/LoopsDefiner.fwd.hh>
#include <protocols/loops/loops_definers/LoopsDefinerCreator.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>


#include <core/types.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace loops {
namespace loops_definers {

/// @brief creator for the LoopsDefiner class
class LoopsDefinerCreator : public utility::pointer::ReferenceCount {
public:
	LoopsDefinerCreator() {}
	virtual ~LoopsDefinerCreator() {}

	virtual LoopsDefinerOP create_loops_definer() const = 0;
	virtual std::string type_name() const = 0;
};

} //namespace
} //namespace
} //namespace

#endif // include guard
