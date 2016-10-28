// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/enzymes/EnzymeData.hh
/// @brief   Definitions for EnzymeData.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_enzymes_EnzymeData_HH
#define INCLUDED_core_enzymes_EnzymeData_HH

// Project header
#include <core/types.hh>

// Utility header
#include <utility/vector1.hh>

// C++ headers
#include <string>


namespace core {
namespace enzymes {

/// @brief  The type of consensus sequence stored in an instance of EnzymeData.
enum ConsensusSequenceType {
	AA = 1,
	NA,
	SACCHARIDE
};


/// @brief  A structure for storing reaction information for specific virtual enzymes.
struct EnzymeData {
	// Data read from the database
	std::string consensus_sequence;
	ConsensusSequenceType cs_type;
	core::uint cs_resnum;  // the position in the consensus sequence that is potentially modified
	std::string atom_to_modify;
	core::Real efficiency;  // ratio of times this enzyme acts on a recognized site
	utility::vector1< std::string > second_substrates_or_byproducts;  // We do not know the kind of reaction yet.

	// Derived data
	utility::vector1< utility::vector1< std::string > > consensus_residues;  // A position may have options.
};  // struct EnzymeData

}  // namespace enzymes
}  // namespace core

#endif  // INCLUDED_core_enzymes_EnzymeData_HH
