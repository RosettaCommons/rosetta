// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/jumping/PairingTemplate
/// @brief header file for ClassicAbinitio protocol
/// @details
///  from converting jumping_pairings.cc of rosetta++ into mini
///
///
///
/// @author Oliver Lange

#ifndef INCLUDED_protocols_jumping_RandomSheetBuilder_hh
#define INCLUDED_protocols_jumping_RandomSheetBuilder_hh

// Unit Headers
#include <protocols/jumping/SheetBuilder.hh>

// Package Headers
#include <protocols/jumping/SameStrand.fwd.hh>
#include <core/scoring/dssp/PairingsList.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2A.fwd.hh>

#include <utility/vector1.hh>


//// C++ headers
//#include <cstdlib>
//#include <string>
//#include <vector>

namespace protocols {
namespace jumping {


/// @brief select jumps to build a given topology
/// @detail this class encapsulates the functionality of choose_random_pairings in jumping_pairings.cc of Rosetta++
class RandomSheetBuilder : public SheetBuilder {
public:
	RandomSheetBuilder( core::fragment::SecondaryStructureOP, core::scoring::dssp::PairingsList const&, SheetTopology const& );

protected:
	//default do nothing always use input_sheet_sizes_ as sheet_sizes_.
	virtual SheetTopology create_new_random_topol() const;
	std::string type_name() const {
		return "RandomSheetBuilder";
	}
private:
	SheetTopology input_sheet_sizes_;

};

} //protocols
} //jumping

#endif

