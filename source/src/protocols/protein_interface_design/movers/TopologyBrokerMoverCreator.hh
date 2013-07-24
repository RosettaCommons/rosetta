// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///@file protocols/protein_interface_design/movers/TopologyBrokerMoverCreator.hh
///@brief This class will create instances of protocols::moves::Mover TopologyBrokerMover for the protocols::moves::MoverFactory
///@author Andrew Leaver-Fay via code_writer.py (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_protein_interface_design_movers_TopologyBrokerMoverCreator_hh
#define INCLUDED_protocols_protein_interface_design_movers_TopologyBrokerMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

class TopologyBrokerMoverCreator : public protocols::moves::MoverCreator {
public:
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
};

} //movers
} //protein_interface_design
} //protocols
#endif

