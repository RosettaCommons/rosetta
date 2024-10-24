// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/minimization_packing/PackRotamersMoverPackRotamersMoverCreator.hh
/// @brief This class will create instances of protocols::moves::Mover PackRotamersMover for the protocols::moves::MoverFactory
/// @author Andrew Leaver-Fay via code_writer.py (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_minimization_packing_PackRotamersMoverLazyCreator_hh
#define INCLUDED_protocols_minimization_packing_PackRotamersMoverLazyCreator_hh

#include <protocols/minimization_packing/PackRotamersMoverCreator.hh>

namespace protocols {
namespace minimization_packing {

class PackRotamersMoverLazyCreator : public protocols::minimization_packing::PackRotamersMoverCreator {
public:
	static std::string mover_name();

};

}
}

#endif

