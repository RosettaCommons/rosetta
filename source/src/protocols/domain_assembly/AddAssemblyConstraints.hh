// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/domain_assembly/AddAssemblyConstraints.fwd.hh
/// @brief  Adds constraints for the assembling process
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_domain_assembly_AddAssemblyConstraints_hh
#define INCLUDED_protocols_domain_assembly_AddAssemblyConstraints_hh


#include <protocols/domain_assembly/AddAssemblyConstraints.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <core/types.hh>
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace domain_assembly {

class AddAssemblyConstraints : public protocols::moves::Mover {
public:

	AddAssemblyConstraints() : Mover("AddAssemblyConstraints") {}

	std::string get_name() const override {
		return "AddAssemblyConstraints";
	}

	void apply( core::pose::Pose & pose ) override;
};


}//namespace domain_assembly
}//namespace protocols

#endif

