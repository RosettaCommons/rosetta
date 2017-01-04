// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/silent/SilentFileOptions.hh
/// @brief  Options for constructing or writing a pose from / to a silent file
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_io_silent_SilentFileOptions_hh
#define INCLUDED_core_io_silent_SilentFileOptions_hh


// Unit Headers
#include <core/io/silent/SilentFileOptions.fwd.hh>

// Platform Headers
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ Headers
#include <string>

namespace core {
namespace io {
namespace silent {

class SilentFileOptions : public utility::pointer::ReferenceCount
{
public:
	SilentFileOptions();
	SilentFileOptions( utility::options::OptionCollection const & options );
	void read_from_global_options();
	void read_from_options( utility::options::OptionCollection const & options );
	void read_from_tag( utility::tag::TagCOP tag );

	static void list_read_options( utility::options::OptionKeyList & read_options );
	static void append_attributes_for_tag_parsing(
		utility::tag::XMLSchemaDefinition & xsd,
		utility::tag::AttributeList & attributes );

	bool keep_input_scores() const;
	bool out_user_tag_set() const;
	std::string out_user_tag() const;
	bool out_weight_silent_scores() const;
	std::string in_silent_score_prefix() const;
	bool in_silent_scores_wanted_set() const;
	utility::vector1< std::string > in_silent_scores_wanted() const;
	bool in_fullatom() const;
	bool write_failures_set() const;
	bool write_failures() const;
	bool in_silent_struct_type_set() const;
	std::string in_silent_struct_type() const;
	bool out_silent_struct_type_set() const;
	std::string out_silent_struct_type() const;
	bool read_through_errors() const;
	bool select_random_set() const;
	int  select_random() const;
	int  select_range_start() const;
	int  select_range_len() const;
	int  select_range_mul() const;
	bool force_silent_bitflip_on_read_set() const;
	bool force_silent_bitflip_on_read() const;
	bool print_all_score_headers() const;

	void keep_input_scores( bool setting );
	void out_user_tag( std::string setting );
	void out_weight_silent_scores( bool setting );
	void in_silent_score_prefix( std::string setting );
	void in_silent_scores_wanted( utility::vector1< std::string > const & setting );
	void in_fullatom( bool setting );
	void write_failures( bool setting );
	void in_silent_struct_type( std::string setting );
	void out_silent_struct_type( std::string setting );
	void read_through_errors( bool setting );
	void select_random( int setting );
	void select_range_start( int setting );
	void select_range_len( int setting );
	void select_range_mul( int setting );
	void force_silent_bitflip_on_read( bool setting );
	void print_all_score_headers( bool setting );

private:

	// Read by SilentStruct
	bool keep_input_scores_;
	bool out_user_tag_set_;
	std::string out_user_tag_;
	bool out_weight_silent_scores_;
	std::string in_silent_score_prefix_;
	bool in_silent_scores_wanted_set_;
	utility::vector1< std::string > in_silent_scores_wanted_;

	// Read by BinarySecStruct
	bool in_fullatom_;
	bool write_failures_set_;
	bool write_failures_;

	// Read by SSFactory
	bool in_silent_struct_type_set_;
	std::string in_silent_struct_type_;
	bool out_silent_struct_type_set_;
	std::string out_silent_struct_type_;

	// read by SilentFileData
	bool read_through_errors_;
	bool select_random_set_;
	int  select_random_;
	int  select_range_start_;
	int  select_range_len_;
	int  select_range_mul_;
	bool force_silent_bitflip_on_read_set_;
	bool force_silent_bitflip_on_read_;
	bool print_all_score_headers_;
};


} // namespace
} // namespace
} // namespace


#endif
