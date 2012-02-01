// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file 
/// @brief
/// @author Javier Castellanos	(javiercv@uw.edu)



#ifndef INCLUDED_devel_constrained_sequence_design_SequenceConstraint_fwd_HH
#define INCLUDED_devel_constrained_sequence_design_SequenceConstraint_fwd_HH

#include <utility/pointer/owning_ptr.hh>

// Boost
#include <boost/unordered_map.hpp>

namespace devel {
namespace constrained_sequence_design {

class SequenceConstraint;

typedef utility::pointer::owning_ptr< SequenceConstraint > SequenceConstraintOP;
typedef utility::pointer::owning_ptr< SequenceConstraint const > SequenceConstraintCOP;

typedef boost::unordered_map<std::string, SequenceConstraintOP> SequenceConstraintMap;
} // devel 
} // constrained_sequence_design

#endif
