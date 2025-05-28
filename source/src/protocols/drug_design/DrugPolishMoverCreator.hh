// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///@file protocols/drug_design/DrugPolishMoverDrugPolishMoverCreator.hh
///@brief This class will create instances of Mover DrugPolishMover for the MoverFactory
///@author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_DrugPolishMoverCreator_hh
#define INCLUDED_protocols_drug_design_DrugPolishMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace drug_design {

class DrugPolishMoverCreator : public protocols::moves::MoverCreator {
public:
	moves::MoverOP create_mover() const;
	std::string keyname() const;
	static std::string mover_name();
};

}
}

#endif

