// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols//ParatopeEpitopeSiteConstraintMoverCreator.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_constraints_ParatopeEpitopeSiteConstraintMoverCreator_hh
#define INCLUDED_protocols_antibody_constraints_ParatopeEpitopeSiteConstraintMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace antibody {
namespace constraints {

class ParatopeEpitopeSiteConstraintMoverCreator : public protocols::moves::MoverCreator {
public:

	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();



};

} //constraints
} //antibody
} //protocols

#endif //constraints_ParatopeEpitopeSiteConstraintMoverCreator_hh
