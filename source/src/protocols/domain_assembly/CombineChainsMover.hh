// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/domain_assembly/CombineChainsMover.fwd.hh
/// @brief  A mover for assembling domains
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_domain_assembly_CombineChainsMover_hh
#define INCLUDED_protocols_domain_assembly_CombineChainsMover_hh

#include <protocols/domain_assembly/CombineChainsMover.fwd.hh>

#include <protocols/moves/Mover.hh>

#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace domain_assembly {

class CombineChainsMover : public protocols::moves::Mover {
public:

	CombineChainsMover() : Mover("CombineChainsMover") {}

	virtual std::string get_name() const {
		return "CombineChainsMover";
	}

	void apply( core::pose::Pose & pose );
};

}//namespace domain_assembly
}//namespace protocols

#endif

