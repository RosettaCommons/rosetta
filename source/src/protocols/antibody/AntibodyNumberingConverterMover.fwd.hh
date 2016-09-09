// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/AntibodyNumberingConverterMover.fwd.hh
/// @brief Converts numbering schemes of an antibody.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_AntibodyNumberingConverterMover_fwd_hh
#define INCLUDED_protocols_antibody_AntibodyNumberingConverterMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace antibody {

class AntibodyNumberingConverterMover;

typedef utility::pointer::shared_ptr< AntibodyNumberingConverterMover > AntibodyNumberingConverterMoverOP;
typedef utility::pointer::shared_ptr< AntibodyNumberingConverterMover const > AntibodyNumberingConverterMoverCOP;



} //protocols
} //antibody


#endif //INCLUDED_protocols_antibody_AntibodyNumberingConverterMover_fwd_hh





