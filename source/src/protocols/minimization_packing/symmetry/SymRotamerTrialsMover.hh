// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author  Ingemar Andre

#ifndef INCLUDED_protocols_minimization_packing_symmetry_SymRotamerTrialsMover_hh
#define INCLUDED_protocols_minimization_packing_symmetry_SymRotamerTrialsMover_hh

// Unit headers
#include <protocols/minimization_packing/symmetry/SymRotamerTrialsMover.fwd.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.fwd.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>

#ifdef PYROSETTA
#include <protocols/moves/MonteCarlo.hh>
#endif

namespace protocols {
namespace minimization_packing {
namespace symmetry {

/// @brief This class is provided for historical support only
/// The RotamerTrialsMover should handle symmetry transparently and is preferred.
class SymRotamerTrialsMover : public protocols::minimization_packing::RotamerTrialsMover {
public:

public:

	// Use the parent Mover
	using RotamerTrialsMover::RotamerTrialsMover;

	~SymRotamerTrialsMover();

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
};

/// @brief This class is provided for historical support only
/// The EnergyCutRotamerTrialsMover should handle symmetry transparently and is preferred.
class SymEnergyCutRotamerTrialsMover : public EnergyCutRotamerTrialsMover {
public:

	using EnergyCutRotamerTrialsMover::EnergyCutRotamerTrialsMover;

	~SymEnergyCutRotamerTrialsMover();
public:

	std::string get_name() const override;

};

} // symmetry
} // moves
} // protocols

#endif
