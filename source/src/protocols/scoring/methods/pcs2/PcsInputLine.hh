// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs2/PcsInputLine.hh
///
/// @brief Class that hold a line data of the input file
///
/// @details This class hold the information of a pcs value:
/// atom name, residue name, PCS value, PCS experimental
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsInputLine_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsInputLine_hh

// Package headers

// Project headers
#include <core/types.hh>

// Utility headers

// Numeric headers

// Objexx headers

// C++ headers
#include <string>

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

///////////////////////////////////////////////////////////////////////////
/// @brief PcsInputLine class: hold a line of the input file information (.npc format)
/// One PcsInputLine per line in the input file
class PcsInputLine {

public:
	PcsInputLine(); //Construc

	~PcsInputLine(); //Destruct

	PcsInputLine(PcsInputLine const & other); //Copy

	PcsInputLine &
	operator=( PcsInputLine const & other ); //=

	PcsInputLine(core::Size const residue_num,
		std::string const atom_name,
		core::Real const PCS_experimental,
		core::Real const PCS_tolerance
	);

	/// @brief Give me the residue number of this PcsInputLine
	core::Size
	get_residue_num() const;

	/// @brief Give me the atom name of this PcsInputLine
	std::string
	get_atom_name() const;

	/// @brief Give me the experimental PCS of this PcsInputLine
	core::Real
	get_PCS_experimental() const;

	/// @brief Give me the tolerance for the PCS of this PcsInputLine
	core::Real
	get_PCS_tolerance() const;

	/// @brief Output myself on the stream
	friend std::ostream &
	operator<<(std::ostream& out, const PcsInputLine &me);

private:
	core::Size const residue_num_;
	std::string const atom_name_;
	core::Real const PCS_experimental_;
	core::Real const PCS_tolerance_;
};

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols

#endif
