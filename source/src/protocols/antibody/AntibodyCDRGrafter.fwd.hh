// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/AntibodyCDRGrafter.fwd.hh
/// @brief Class to graft CDR loops from an antibody to a new antibody or from a CDR pose into a different antibody.  
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_AntibodyCDRGrafter_fwd_hh
#define INCLUDED_protocols_antibody_AntibodyCDRGrafter_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace antibody {

class AntibodyCDRGrafter;

typedef utility::pointer::shared_ptr< AntibodyCDRGrafter > AntibodyCDRGrafterOP;
typedef utility::pointer::shared_ptr< AntibodyCDRGrafter const > AntibodyCDRGrafterCOP;


}//protocols
}//antibody

#endif	//INCLUDED_protocols/antibody_AntibodyCDRGrafter_fwd_hh





