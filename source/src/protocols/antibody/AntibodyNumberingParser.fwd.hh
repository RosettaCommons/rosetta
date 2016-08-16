// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/AntibodyNumberingParser.fwd.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_AntibodyNumberingParser_fwd_hh
#define INCLUDED_protocols_antibody_AntibodyNumberingParser_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace antibody {

class AntibodyNumberingParser;
typedef utility::pointer::shared_ptr< AntibodyNumberingParser > AntibodyNumberingParserOP;
typedef utility::pointer::shared_ptr< AntibodyNumberingParser const > AntibodyNumberingParserCOP;

class PDBLandmark;
typedef utility::pointer::shared_ptr< PDBLandmark > PDBLandmarkOP;
typedef utility::pointer::shared_ptr< PDBLandmark const > PDBLandmarkCOP;
}
}

#endif //INCLUDED_protocols_antibody_design_AntibodyNumberingParser.fwd.hh


