// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

 //////////////////////////////////////////////
 ///
 /// @file protocols/scoring/methods/pcs2/PcsInputFile.hh
 ///
 /// @brief Read all input from a .npc input file, and hold the data in the class.
 /// One file per lanthanide data.
 ///
 /// @details
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

#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsInputFile_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsInputFile_hh

// Package headers
#include <protocols/scoring/methods/pcs2/PcsInputLine.hh>

// Project headers

// Utility headers
#include <utility/vector1.hh>

// Numeric headers

// Objexx headers

// C++ headers

namespace protocols{
namespace scoring{
namespace methods{
namespace pcs2{

//////////////////////////////////////////////////////////
/// @brief PcsInputFile contain all the information of a .npc file
/// one per lanthanide.
class PcsInputFile {
private:
	std::string const filename_;
	utility::vector1<PcsInputLine> PcsInputLine_all_;
	core::Real const weight_;

	/// @brief read the file containing the PCS. Private and called by the constructor
	void
	read_PCS_file();

public:
	PcsInputFile(); //construct

	~PcsInputFile(); //destruct

	PcsInputFile(PcsInputFile const & other); //copy

	PcsInputFile &
	operator=( PcsInputFile const & other ); // =

	PcsInputFile(std::string const & filename, core::Real const my_weight);

	/// @brief Give me the name of the file
	std::string
	get_filename() const;

	/// @brief Give me the weight associated
	core::Real
	get_weight() const;

	/// @brief Give me the vector of all the line of the file
	utility::vector1<PcsInputLine> &
	get_PcsInputLine_all();

	/// @brief Print me on the stream
	friend std::ostream &
	operator<<(std::ostream& out, const PcsInputFile &me);
};


}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols

#endif
