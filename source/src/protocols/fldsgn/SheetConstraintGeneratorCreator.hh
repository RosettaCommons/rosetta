// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/fldsgn/SheetConstraintGeneratorCreator.hh
/// @brief This class will create instances of the SheetConstraintGenerator
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_fldsgn_SheetConstraintGeneratorCreator_hh
#define INCLUDED_protocols_fldsgn_SheetConstraintGeneratorCreator_hh

#include <protocols/constraint_generator/ConstraintGeneratorCreator.hh>

namespace protocols {
namespace fldsgn {

class SheetConstraintGeneratorCreator : public protocols::constraint_generator::ConstraintGeneratorCreator {
public:
	virtual protocols::constraint_generator::ConstraintGeneratorOP create_constraint_generator() const;
	virtual std::string keyname() const;
	static std::string constraint_generator_name();
};

}
}

#endif
