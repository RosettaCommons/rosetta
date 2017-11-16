// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief   Method definitions for MatchScoreWriter.
/// @author  Kyle Barlow (kb@kylebarlow.com)

// Unit headers
#include <protocols/match/output/MatchScoreWriter.hh>

// Package headers

// Project headers

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers

// C++ headers
#include <iostream>
#include <fstream>

// Construct tracer.
static basic::Tracer TR( "protocols.match.output.MatchScoreWriter" );

namespace protocols {
namespace match {
namespace output {

using namespace core;

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////
// Default constructor
/// @details  default is to not really do much of anything, unless filename output is set
MatchScoreWriter::MatchScoreWriter() : utility::pointer::ReferenceCount()
{
	init();
}

// Copy constructor
MatchScoreWriter::MatchScoreWriter(MatchScoreWriter const & object_to_copy) : utility::pointer::ReferenceCount()
{
	copy_data(*this, object_to_copy);
}

// Assignment operator
MatchScoreWriter &
MatchScoreWriter::operator=(MatchScoreWriter const & object_to_copy)
{
	// Abort self-assignment.
	if ( this == &object_to_copy ) {
		return *this;
	}

	copy_data(*this, object_to_copy);
	return *this;
}

// Destructor
MatchScoreWriter::~MatchScoreWriter() {}


// Standard Rosetta methods ////////////////////////////////////////////////////
// General methods
void
MatchScoreWriter::show(std::ostream & output) const
{
	output << output_header_ << std::endl;

	for ( core::Size i = 1 ; i<=match_names_.size() ; i++ ) {
		output << match_names_[i] << "\t" << match_scores_[i] << std::endl;
	}
}


// Accessors/Mutators
/// @brief Add a match (along with score) to be outputted in score file
void MatchScoreWriter::add_match(std::string match_name, core::Real score)
{
	match_names_.push_back(match_name);
	match_scores_.push_back(score);
}

/// @brief Add an unnamed match (along with score) to be outputed in score file
/// @details Each additional match will have an incremented name
void MatchScoreWriter::add_match(core::Real score)
{
	std::ostringstream match_name;
	match_name << "MatchOutput_" << generic_output_number_++;
	match_names_.push_back(match_name.str());
	match_scores_.push_back(score);
}

/// @brief Write out the match scores once all scores have been added
void MatchScoreWriter::write_match_scores()
{
	// Check to see if there is anything to do
	if ( (!write_file_) || (match_names_.size() == 0) ) {
		return;
	}

	std::ofstream fileout;
	fileout.open(score_output_filename_.c_str());
	fileout << output_header_ << std::endl;

	for ( core::Size i = 1 ; i<=match_names_.size() ; i++ ) {
		fileout << match_names_[i] << "\t" << match_scores_[i] << std::endl;
	}

	TR << "Match scores written to file: " << score_output_filename_ << std::endl;
	fileout.close();
}


/// @brief Set output filename, which also turns on output
void
MatchScoreWriter::set_output_filename(std::string output_filename) {
	score_output_filename_ = output_filename;
	set_output( true );
}

/// @brief Set whether or not to output
void
MatchScoreWriter::set_output(bool output_on) {
	write_file_ = output_on;
}

// Private methods /////////////////////////////////////////////////////////////
// Initialize data members.
void
MatchScoreWriter::init()
{
	score_output_filename_ = "score.out";
	output_header_="match_name\tmatch_score";
	write_file_ = false;
	generic_output_number_ = 1;
}

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
MatchScoreWriter::copy_data(
	MatchScoreWriter object_to_copy_to,
	MatchScoreWriter object_to_copy_from)
{
	object_to_copy_to.score_output_filename_ = object_to_copy_from.score_output_filename_;
	object_to_copy_to.output_header_ = object_to_copy_from.output_header_;
	object_to_copy_to.write_file_ = object_to_copy_from.write_file_;
	object_to_copy_to.generic_output_number_ = object_to_copy_from.generic_output_number_;
}


// Friend methods //////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that MatchScoreWriter can be "printed" in PyRosetta).
std::ostream &
operator<<(std::ostream & output, MatchScoreWriter const & object_to_output)
{
	object_to_output.show(output);
	return output;
}

}  // namespace output
}  // namespace match
}  // namespace protocols
