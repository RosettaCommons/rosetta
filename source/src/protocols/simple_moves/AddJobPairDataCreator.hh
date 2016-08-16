// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/AddJobPairDataCreator.hh
/// @brief MoverCreator definitions for the AddJobPairData mover
/// @author Sam DeLuca <Samuel.l.deluca@vanderbilt.edu)

#ifndef INCLUDED_protocols_simple_moves_AddJobPairDataCreator_hh
#define INCLUDED_protocols_simple_moves_AddJobPairDataCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace simple_moves {

class AddJobPairDataCreator : public protocols::moves::MoverCreator {
public:
	virtual moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
};

}
}


#endif /* ADDJOBSTRINGPAIR_HH_ */
