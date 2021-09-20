// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rosetta_scripts/RosettaScriptsSchemaValidator.hh
/// @brief  A singleton class for generating the schema accepted by the
///         RosettaScriptsParser and for holding the libxml2-library
///         schema validation object.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_rosetta_scripts_RosettaScriptsSchemaValidator_hh
#define INCLUDED_protocols_rosetta_scripts_RosettaScriptsSchemaValidator_hh

// Unit Headers
#include <protocols/rosetta_scripts/RosettaScriptsSchemaValidator.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <utility/SingletonBase.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/XMLSchemaValidation.hh>

// c++ headers
#include <map>
#include <set>

namespace protocols {
namespace rosetta_scripts {


class RosettaScriptsSchemaValidator : public utility::SingletonBase< RosettaScriptsSchemaValidator >
{
public:
	friend class utility::SingletonBase< RosettaScriptsSchemaValidator >;

	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;

public:
	virtual ~RosettaScriptsSchemaValidator();

	static std::string xsd_for_rosetta_scripts();
	static void write_ROSETTASCRIPTS_complex_type( utility::tag::XMLSchemaDefinition & xsd );
	static std::string rosetta_scripts_element_name(); // returns "ROSETTASCRIPTS"
	static std::string rosetta_scripts_complex_type_naming_func( std::string const & element_name );

	utility::tag::XMLSchemaValidatorCOP validator() const;

private:
	RosettaScriptsSchemaValidator();

	RosettaScriptsSchemaValidator( RosettaScriptsSchemaValidator const & ) = delete;
	RosettaScriptsSchemaValidator const & operator = ( RosettaScriptsSchemaValidator const & ) = delete;

private:
	/// @brief The object used to validate XML input files against the schema
	utility::tag::XMLSchemaValidatorOP validator_;

};

} //namespace rosetta_scripts
} //namespace protocols

#endif
