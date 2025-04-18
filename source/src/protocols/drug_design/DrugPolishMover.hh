// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/DrugPolishMover.hh
/// @brief Exhaustively enumerate reactions on a ligand
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_DrugPolishMover_hh
#define INCLUDED_protocols_drug_design_DrugPolishMover_hh

// Unit header
#include <protocols/drug_design/DrugPolishMover.fwd.hh>

#include <protocols/drug_design/DrugDesignMover.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/AtomRefMapping.hh>

// Project Headers

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

// External headers

// C/C++ headers
#include <string>

namespace protocols {
namespace drug_design {

class DrugPolishMover : public protocols::drug_design::DrugDesignMover {
	// TODO: Checkpointing
public:
	/// @brief default constructor
	DrugPolishMover();

	/// @brief destructor
	~DrugPolishMover() override;

	/// @brief create copy constructor
	protocols::moves::MoverOP clone() const override;

	/// @brief create this type of objectt
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief apply DrugPolishMover
	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override;

public: // accessors

	/// @brief The bonus given to scoring for applying (another) reaction
	core::Real bonus() const { return bonus_; }

public: // setters

	/// @brief The bonus given to scoring for applying (another) reaction
	void bonus(core::Real setting) { bonus_ = setting; }

public: // Other functions

	/// @brief parse xml file
	void
	parse_my_tag( TagCOP tag, basic::datacache::DataMap & data ) override;

private: // Data

	/// @brief The bonus given to scoring for applying (another) reaction
	core::Real bonus_;
};

} // namespace drug_design
} // namespace protocols

#endif
