// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/splice/RBInMoverCreator.hh
/// @brief This class will create instances of protocols::moves::Mover RBInMover for the protocols::moves::MoverFactory
/// @author Christoffer Norn

#ifndef INCLUDED_devel_splice_RBInMoverCreator_hh
#define INCLUDED_devel_splice_RBInMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace devel {
namespace splice {

class RBInMoverCreator : public protocols::moves::MoverCreator {
public:
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
};

}
}

#endif // INCLUDED_protocols_simple_moves_RBInMoverCreator_hh

