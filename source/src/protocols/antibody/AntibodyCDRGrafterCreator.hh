// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/AntibodyCDRGrafterCreator.hh
/// @brief Class to graft CDR loops from an antibody to a new antibody or from a CDR pose into a different antibody.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_AntibodyCDRGrafterCreator_hh
#define INCLUDED_protocols_antibody_AntibodyCDRGrafterCreator_hh

#include <protocols/moves/MoverCreator.hh>


namespace protocols {
namespace antibody {

class AntibodyCDRGrafterCreator : public protocols::moves::MoverCreator {

public:

	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();



};


}//protocols
}//antibody


#endif //INCLUDED_protocols_antibody_AntibodyCDRGrafterCreator_hh

