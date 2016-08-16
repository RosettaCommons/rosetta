// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/scoring_grid/GridInitMoverCreator.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_scoring_grid_GridInitMoverCreator_HH_
#define INCLUDED_protocols_qsar_scoring_grid_GridInitMoverCreator_HH_

#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class GridInitMoverCreator : public protocols::moves::MoverCreator {
public:
	virtual moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
};

}
}
}


#endif /* GRIDINITMOVERCREATOR_HH_ */
