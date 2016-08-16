// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief   Declarations and simple accessor/mutator definitions for MatchScoreWriter.
/// @details A class to save and output matcher data, if matches are scored
/// @author  Kyle Barlow (kb@kylebarlow.com)

#ifndef INCLUDED_protocols_match_output_MatchScoreWriter_HH
#define INCLUDED_protocols_match_output_MatchScoreWriter_HH

// Unit header
#include <protocols/match/output/MatchScoreWriter.fwd.hh>

// Package headers

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Numeric headers

// C++ headers
#include <iostream>

namespace protocols {
namespace match {
namespace output {

/// @details Writes out a score file based on the sum of RMS of geometric constrainst satisfaction with the desired ligan positioning
class MatchScoreWriter : public utility::pointer::ReferenceCount {
public:
	// Standard methods ////////////////////////////////////////////////////////
	/// @brief  Default constructor
	MatchScoreWriter();

	/// @brief  Copy constructor
	MatchScoreWriter(MatchScoreWriter const & object_to_copy);

	// Assignment operator
	MatchScoreWriter & operator=(MatchScoreWriter const & object_to_copy);

	// Destructor
	~MatchScoreWriter();

	// Standard Rosetta methods ////////////////////////////////////////////////
	/// @brief  Generate string representation of MatchScoreWriter for debugging purposes.
	virtual void show(std::ostream & output=std::cout) const;

	// Insertion operator (overloaded so that MatchScoreWriter can be "printed" in PyRosetta).
	friend std::ostream & operator<<(std::ostream & output, MatchScoreWriter const & object_to_output);


	// Accessors/Mutators
	/// @brief Add a match name (along with score) to be outputed in score file
	virtual void add_match(std::string match_name, core::Real score);

	/// @brief Add an unnamed match (along with score) to be outputed in score file
	/// @details Each additional match will have an incremented name
	virtual void add_match(core::Real score);

	/// @brief Write out the match scores once all scores have been added
	virtual void write_match_scores();

	/// @brief Set output filename, which also turns on output
	virtual void set_output_filename(std::string output_filename);

	/// @brief Set whether or not to output
	virtual void set_output(bool output_on);

private:
	// Private methods /////////////////////////////////////////////////////////
	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	virtual void copy_data(MatchScoreWriter object_to_copy_to, MatchScoreWriter object_to_copy_from);

	// Initialize data members.
	virtual void init();


	// Private data ////////////////////////////////////////////////////////////
	std::string score_output_filename_;
	utility::vector1<std::string> match_names_;
	utility::vector1<core::Real> match_scores_;
	std::string output_header_;
	bool write_file_;
	// Used to increment output names if no name is given
	core::Real generic_output_number_;

};  // class MatchScoreWriter

}  // namespace output
}  // namespace match
}  // namespace protocols

#endif  // INCLUDED_protocols_match_output_MatchScoreWriter_HH
