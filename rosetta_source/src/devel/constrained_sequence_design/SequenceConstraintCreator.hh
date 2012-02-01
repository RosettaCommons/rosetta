// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   
/// @brief  
/// @author  Javier Castellanos	( javiercv@uw.edu )

#ifndef INCLUDED_devel_constrained_sequence_design_SequenceConstraintCreator_hh
#define INCLUDED_devel_constrained_sequence_design_SequenceConstraintCreator_hh

// Unit Headers
#include <devel/constrained_sequence_design/SequenceConstraint.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers

namespace devel {
namespace constrained_sequence_design {

/// @brief Abstract base class for a SequenceConstraint factory; the Creator class is responsible for
/// creating a particular SequenceConstraint derived class.
class SequenceConstraintCreator : public utility::pointer::ReferenceCount
{
public:
	SequenceConstraintCreator();
	virtual ~SequenceConstraintCreator();

	virtual SequenceConstraintOP create() const = 0;
	virtual std::string keyname() const = 0;
};

typedef utility::pointer::owning_ptr< SequenceConstraintCreator > SequenceConstraintCreatorOP;
typedef utility::pointer::owning_ptr< SequenceConstraintCreator const > SequenceConstraintCreatorCOP;

} // namespace devel
} // constrained_sequence_design

#endif         
