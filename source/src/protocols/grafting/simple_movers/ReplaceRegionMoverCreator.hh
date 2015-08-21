// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/grafting/simple_movers/DeleteRegionMover.hh
/// @brief
/// @author  Jared Adolf-Bryfogle

#ifndef INCLUDED_protocols_grafting_simple_movers_REPLACEREGIONMOVERCREATOR_HH
#define INCLUDED_protocols_grafting_simple_movers_REPLACEREGIONMOVERCREATOR_HH

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace grafting {
namespace simple_movers {

class ReplaceRegionMoverCreator : public protocols::moves::MoverCreator {
public:
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();



};

}
}
}

#endif //#ifndef INCLUDED_protocols/antibody_design_REPLACEREGIONMOVERCREATOR_HH

