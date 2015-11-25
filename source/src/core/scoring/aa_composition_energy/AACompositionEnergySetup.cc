// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/aa_composition_energy/AACompositionEnergySetup.cc
/// @brief A helper for the EnergyMethod that penalizes deviation from a desired amino acid composition.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

// Unit headers
#include <core/scoring/aa_composition_energy/AACompositionEnergySetup.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>

// Options system
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// File I/O
#include <basic/database/open.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

// Other Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <math.h>

namespace core {
namespace scoring {
namespace aa_composition_energy {

static THREAD_LOCAL basic::Tracer TR("core.scoring.aa_composition_energy.AACompositionEnergySetup");

/**************************************************
AACompositionPropertiesSet class:
A helper class that stores:
- names of residue types that WILL be counted (TYPE)
- names of residue types that WILL NOT be counted (NOT_TYPE)
- properties that MUST be present in order for a residue to be counted (PROPERTIES).
- properties, any of which is sufficient for a residue to be counted (OR_PROPERTIES).
- properties that MUST NOT be present in order for a residue to be counted (NOT_PROPERTIES).
The logic is as follows: a residue is counted if and only if [any TYPE matches] OR [ (no NOT_TYPE matches) AND
( {all PROPERTIES match} OR {any OR_PROPERTIES match} OR {no PROPERTIES defined AND no OR_PROPERTIES defined } ) AND
( no NOT_PROPERTIES match) ]
**************************************************/


/// @brief Default constructor for AACompositionPropertiesSet.
///
AACompositionPropertiesSet::AACompositionPropertiesSet() :
	included_types_(),
	excluded_types_(),
	included_properties_(),
	or_properties_(),
	excluded_properties_()
{}

/// @brief Constructor for AACompositionPropertiesSet that takes lists of
/// included and excluded properties.
AACompositionPropertiesSet::AACompositionPropertiesSet(
	utility::vector1< std::string > const &included_type_strings,
	utility::vector1< std::string > const &excluded_type_strings,
	utility::vector1< std::string > const &included_properties_strings,
	utility::vector1< std::string > const &or_properties_strings,
	utility::vector1< std::string > const &excluded_properties_strings
) :
	included_types_(),
	excluded_types_(),
	included_properties_(),
	or_properties_(),
	excluded_properties_()
{
	if ( !included_type_strings.empty()       ) parse_included_types( included_type_strings );
	if ( !excluded_type_strings.empty()       ) parse_excluded_types( excluded_type_strings );
	if ( !included_properties_strings.empty() ) parse_included_properites( included_properties_strings );
	if ( !or_properties_strings.empty()       ) parse_or_properties( or_properties_strings );
	if ( !excluded_properties_strings.empty() ) parse_excluded_properites( excluded_properties_strings );
}

/// @brief Copy constructor for AACompositionPropertiesSet.
///
AACompositionPropertiesSet::AACompositionPropertiesSet( AACompositionPropertiesSet const &src ) :
	utility::pointer::ReferenceCount(),
	included_types_( src.included_types_ ),
	excluded_types_( src.excluded_types_ ),
	included_properties_( src.included_properties_ ),
	or_properties_( src.or_properties_ ),
	excluded_properties_( src.excluded_properties_ )
{}

/// @brief Default destructor for AACompositionPropertiesSet.
///
AACompositionPropertiesSet::~AACompositionPropertiesSet()
{}

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
AACompositionPropertiesSetOP AACompositionPropertiesSet::clone() const {
	return AACompositionPropertiesSetOP( new AACompositionPropertiesSet(*this) );
}

/// @brief Add a type to the list of types that are always counted.
/// @details Checks that it hasn't yet been added to any list.
void AACompositionPropertiesSet::add_included_type( std::string const &type ) {
	runtime_assert_string_msg( !is_in_list(type, included_types_), "Error in core::scoring::aa_composition_energy::AACompositionPropertiesSet::add_included_type(): This type has already been added ot the list of included types!" );
	runtime_assert_string_msg( !is_in_list(type, excluded_types_), "Error in core::scoring::aa_composition_energy::AACompositionPropertiesSet::add_included_type(): This type has already been added ot the list of excluded types!" );
	included_types_.push_back( type );
	return;
}

/// @brief Add a type to the list of types that are never counted.
/// @details Checks that it hasn't yet been added to any list.
void AACompositionPropertiesSet::add_excluded_type( std::string const &type ) {
	runtime_assert_string_msg( !is_in_list(type, included_types_), "Error in core::scoring::aa_composition_energy::AACompositionPropertiesSet::add_included_type(): This type has already been added ot the list of included types!" );
	runtime_assert_string_msg( !is_in_list(type, excluded_types_), "Error in core::scoring::aa_composition_energy::AACompositionPropertiesSet::add_included_type(): This type has already been added ot the list of excluded types!" );
	excluded_types_.push_back( type );
	return;
}

/// @brief Add a property to the list of properties that must be present.
///
void AACompositionPropertiesSet::add_included_property( core::chemical::ResidueProperty const property ) {
	runtime_assert_string_msg( !is_in_list(property, included_properties_), "Error in core::scoring::aa_composition_energy::AACompositionPropertiesSet::add_included_property(): Property has already been added to the included property list!" );
	runtime_assert_string_msg( !is_in_list(property, or_properties_), "Error in core::scoring::aa_composition_energy::AACompositionPropertiesSet::add_included_property(): Property has already been added to the or properties list!" );
	runtime_assert_string_msg( !is_in_list(property, excluded_properties_), "Error in core::scoring::aa_composition_energy::AACompositionPropertiesSet::add_included_property(): Property has already been added to the excluded property list!" );
	included_properties_.push_back(property);
	return;
}

/// @brief Add a property to the list of properties that, if present, result in the residue being counted
/// if it's not in the excluded_types_ or excluded_properties_ lists.
void AACompositionPropertiesSet::add_or_property( core::chemical::ResidueProperty const property ) {
	runtime_assert_string_msg( !is_in_list(property, included_properties_), "Error in core::scoring::aa_composition_energy::AACompositionPropertiesSet::add_or_property(): Property has already been added to the included property list!" );
	runtime_assert_string_msg( !is_in_list(property, or_properties_), "Error in core::scoring::aa_composition_energy::AACompositionPropertiesSet::add_or_property(): Property has already been added to the or properties list!" );
	runtime_assert_string_msg( !is_in_list(property, excluded_properties_), "Error in core::scoring::aa_composition_energy::AACompositionPropertiesSet::add_or_property(): Property has already been added to the excluded property list!" );
	or_properties_.push_back(property);
	return;
}


/// @brief Add a property to the list of properties that must not be present.
///
void AACompositionPropertiesSet::add_excluded_property( core::chemical::ResidueProperty const property ) {
	runtime_assert_string_msg( !is_in_list(property, included_properties_), "Error in core::scoring::aa_composition_energy::AACompositionPropertiesSet::add_excluded_property(): Property has already been added to the included property list!" );
	runtime_assert_string_msg( !is_in_list(property, or_properties_), "Error in core::scoring::aa_composition_energy::AACompositionPropertiesSet::add_excluded_property(): Property has already been added to the or properties list!" );
	runtime_assert_string_msg( !is_in_list(property, excluded_properties_), "Error in core::scoring::aa_composition_energy::AACompositionPropertiesSet::add_excluded_property(): Property has already been added to the excluded properties list!" );
	excluded_properties_.push_back(property);
	return;
}

/// @brief Take a list of included type strings and add it to the list of included types, checking that
/// none of the types has already been added.
/// @details Populates the included_types_ vector based on the types named in the list.
void AACompositionPropertiesSet::parse_included_types( utility::vector1< std::string > const &typelist ) {
	included_types_.clear();
	core::Size const ntype(typelist.size());
	if ( ntype==0 ) return; //Do nothing if passed no types.
	for ( core::Size i=1; i<=ntype; ++i ) { //Loop through the type list.
		if ( typelist[i] != "" ) add_included_type( typelist[i] ); //Store the string.
	}
	return;
}

/// @brief Take a list of excluded type strings and add it to the list of excluded types, checking that
/// none of the types has already been added.
/// @details Populates the excluded_types_ vector based on the types named in the list.
void AACompositionPropertiesSet::parse_excluded_types( utility::vector1< std::string > const &typelist ) {
	excluded_types_.clear();
	core::Size const ntype(typelist.size());
	if ( ntype==0 ) return; //Do nothing if passed no types.
	for ( core::Size i=1; i<=ntype; ++i ) { //Loop through the type list.
		if ( typelist[i] != "" ) add_excluded_type( typelist[i] ); //Store the string.
	}
	return;
}

/// @brief Take a list of included property strings and parse it.
/// @details Populates the included_properties_ vector based on the properties named in
/// the list.
void AACompositionPropertiesSet::parse_included_properites( utility::vector1< std::string > const &proplist ) {
	included_properties_.clear();
	core::Size const nprop(proplist.size());
	if ( nprop==0 ) return; //Do nothing if passed no properties.
	for ( core::Size i=1; i<=nprop; ++i ) { //Loop through the property list.
		if ( proplist[i] != "" ) add_included_property( parse_property( proplist[i] ) ); //Convert the string to a ResidueProperty and store it.
	}
	return;
}

/// @brief Take a list of included property strings and parse it.
/// @details Populates the included_properties_ vector based on the properties named in
/// the list.
void AACompositionPropertiesSet::parse_or_properties( utility::vector1< std::string > const &proplist ) {
	or_properties_.clear();
	core::Size const nprop(proplist.size());
	if ( nprop==0 ) return; //Do nothing if passed no properties.
	for ( core::Size i=1; i<=nprop; ++i ) { //Loop through the property list.
		if ( proplist[i] != "" ) add_or_property( parse_property( proplist[i] ) ); //Convert the string to a ResidueProperty and store it.
	}
	return;
}

/// @brief Take a list of excluded property strings and parse it.
/// @details Populates the excluded_properties_ vector based on the properties named in
/// the list.
void AACompositionPropertiesSet::parse_excluded_properites( utility::vector1< std::string > const &proplist ) {
	excluded_properties_.clear();
	core::Size const nprop(proplist.size());
	if ( nprop==0 ) return; //Do nothing if passed no properties.
	for ( core::Size i=1; i<=nprop; ++i ) { //Loop through the property list.
		if ( proplist[i] != "" ) add_excluded_property( parse_property( proplist[i] ) ); //Convert the string to a ResidueProperty and store it.
	}
	return;
}

/// @brief Generate a one-line summary of the properties stored in this AACompositionPropertySet
///
std::string AACompositionPropertiesSet::one_line_report() const {
	std::stringstream outstream("");

	outstream << "PROPERTIES: {";

	for ( core::Size i=1, imax=included_properties_.size(); i<=imax; ++i ) {
		outstream << core::chemical::ResidueProperties::get_string_from_property(included_properties_[i]);
		if ( i<imax ) outstream << ",";
	}

	outstream << "} OR_PROPERTIES: {";
	for ( core::Size i=1, imax=or_properties_.size(); i<=imax; ++i ) {
		outstream << core::chemical::ResidueProperties::get_string_from_property(or_properties_[i]);
		if ( i<imax ) outstream << ",";
	}

	outstream << "} NOT_PROPERTIES: {";
	for ( core::Size i=1, imax=excluded_properties_.size(); i<=imax; ++i ) {
		outstream << core::chemical::ResidueProperties::get_string_from_property(excluded_properties_[i]);
		if ( i<imax ) outstream << ",";
	}

	outstream << "} TYPES: {";
	for ( core::Size i=1, imax=included_types_.size(); i<=imax; ++i ) {
		outstream << included_types_[i];
		if ( i<imax ) outstream << ",";
	}

	outstream << "} NOT_TYPES: {";
	for ( core::Size i=1, imax=excluded_types_.size(); i<=imax; ++i ) {
		outstream << excluded_types_[i];
		if ( i<imax ) outstream << ",";
	}

	outstream << "}" << std::endl;

	return outstream.str();
}

/**************************************************
AACompositionEnergySetup class:
**************************************************/

/// @brief Default constructor for AACompositionEnergySetup.
///
AACompositionEnergySetup::AACompositionEnergySetup() :
	utility::pointer::ReferenceCount(),
	property_sets_(),
	expected_by_properties_fraction_(),
	expected_by_properties_absolute_(),
	property_penalties_(),
	property_deviation_ranges_(),
	property_tailfunctions_()
{}

/// @brief Copy constructor for AACompositionEnergySetup.
///
AACompositionEnergySetup::AACompositionEnergySetup( AACompositionEnergySetup const &src ) :
	utility::pointer::ReferenceCount(),
	property_sets_(), //Cloned below
	expected_by_properties_fraction_( src.expected_by_properties_fraction_ ),
	expected_by_properties_absolute_( src.expected_by_properties_absolute_ ),
	property_penalties_( src.property_penalties_ ),
	property_deviation_ranges_( src.property_deviation_ranges_ ),
	property_tailfunctions_( src.property_tailfunctions_ )
{
	property_sets_.clear();
	for ( core::Size i=1, imax=src.property_sets_.size(); i<=imax; ++i ) {
		property_sets_.push_back( src.property_sets_[i]->clone() );
	}

	check_data(); //Double-check that the object has been copied properly.
}

/// @brief Default destructor for AACompositionEnergySetup.
///
AACompositionEnergySetup::~AACompositionEnergySetup() {}

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
AACompositionEnergySetupOP AACompositionEnergySetup::clone() const {
	return AACompositionEnergySetupOP( new AACompositionEnergySetup(*this) );
}

/// @brief Reset all data in this data storage object.
///
void AACompositionEnergySetup::reset() {
	property_sets_.clear();
	expected_by_properties_fraction_.clear();
	expected_by_properties_absolute_.clear();
	property_penalties_.clear();
	property_deviation_ranges_.clear();
	property_tailfunctions_.clear();
	return;
}

/// @brief Initialize from a .comp file.
///
void AACompositionEnergySetup::initialize_from_file( std::string const &filename ) {
	using namespace utility::io;

	//First, reset all data:
	reset();

	std::string filename_formatted = filename;
	if ( utility::file::file_extension(filename_formatted)!="comp" ) filename_formatted+= ".comp";

	izstream infile;
	infile.open( filename_formatted );
	if ( !infile.good() ) {
		filename_formatted = "scoring/score_functions/aa_composition/" + utility::file::file_basename(filename_formatted) + ".comp";
		basic::database::open( infile, filename_formatted );
		runtime_assert_string_msg( infile.good(), "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::initialize_from_file():  Unable to open .comp file for read!" );
	}

	if ( TR.Debug.visible() ) TR.Debug << "Reading aa_composition scoring term setup data from " << filename_formatted << "." << std::endl;
	std::string curline(""); //Buffer for current line.
	utility::vector1< std::string > lines; //Storing all lines

	//Read the file:
	while ( getline(infile, curline) ) {
		if ( curline.size() < 1 ) continue; //Ignore blank lines.
		//Find and process comments:
		std::string::size_type pound = curline.find('#', 0);
		if ( pound == std::string::npos ) {
			lines.push_back( ObjexxFCL::strip( curline, " \t\n" ) );
		} else {
			curline = curline.substr(0,pound);
			lines.push_back( ObjexxFCL::strip( curline, " \t\n" ) );
		}
	}
	infile.close();

	if ( TR.Debug.visible() ) TR.Debug << "Read complete.  Parsing penalty definitions." << std::endl;
	parse_penalty_definitions( lines );

	check_data();

	return;
}

/// @brief Get tail function name from enum.
///
std::string AACompositionEnergySetup::get_tailfunction_name( TailFunction const tf ) const
{
	switch(tf) {
	case tf_linear :
		return "LINEAR";
	case tf_quadratic :
		return "QUADRATIC";
	case tf_constant :
		return "CONSTANT";
	default :
		break;
	}

	return "UNKNOWN";
}

/// @brief Get tail function enum from name.
/// @details This is slow; it calls get_tailfunction_name repeatedly.  Intended only for use during setup.
TailFunction AACompositionEnergySetup::get_tailfunction_from_name( std::string const &name ) const {
	for ( core::Size i=1; i<tf_end_of_list; ++i ) {
		if ( get_tailfunction_name( static_cast<TailFunction>(i) ) == name ) return static_cast<TailFunction>(i);
	}
	return tf_unknown;
}

/// @brief Get a summary of the data stored in this object
///
std::string AACompositionEnergySetup::report() const {
	std::stringstream output("");

	output << "property_sets_:" << std::endl;
	for ( core::Size i=1, imax=property_sets_.size(); i<=imax; ++i ) {
		output << "\t" << i << ":\t" << property_sets_[i]->one_line_report();
	}

	output << "property_penalties_:" << std::endl;
	for ( core::Size i=1, imax=property_penalties_.size(); i<=imax; ++i ) {
		output << "\t" << i << ":" << "\t";
		for ( core::Size j=1, jmax=property_penalties_[i].size(); j<=jmax; ++j ) {
			output << property_penalties_[i][j];
			if ( j<jmax ) output << " ";
			else output << std::endl;
		}
	}

	output << "property_deviation_ranges_:" << std::endl;
	for ( core::Size i=1, imax=property_deviation_ranges_.size(); i<=imax; ++i ) {
		output << "\t" << i << ":\t" << property_deviation_ranges_[i].first << ", " << property_deviation_ranges_[i].second << std::endl;
	}

	output << "expected_by_properties_fraction_ / expected_by_properties_absolute_:" << std::endl;
	for ( core::Size i=1, imax=expected_by_properties_fraction_.size(); i<=imax; ++i ) {
		output << "\t" << i << ":\t" << expected_by_properties_fraction_[i] << " / " << expected_by_properties_absolute_[i] << std::endl;
	}

	return output.str();
}


/// @brief Parse out penalty definition blocks from the data read from file.
///
void AACompositionEnergySetup::parse_penalty_definitions( utility::vector1 < std::string > const &lines ) {
	core::Size const nlines( lines.size() );
	runtime_assert_string_msg( nlines>0, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_penalty_definitions():  No lines were read from file!" );

	bool in_a_block(false);
	utility::vector1 < std::string > lines_subset;

	for ( core::Size i=1; i<=nlines; ++i ) { //Loop through all lines
		std::istringstream curline(lines[i]);
		std::string oneword("");

		curline >> oneword;
		if ( !in_a_block ) { //If we're not already in a PENALTY_DEFINITION block.
			if ( oneword == "PENALTY_DEFINITION" ) {
				in_a_block=true;
				lines_subset.clear(); //Starting with the next line, we'll start filling in this vector
			}
		} else { //If we are already in a PENALTY_DEFINITION block
			runtime_assert_string_msg( oneword != "PENALTY_DEFINITION", "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_penalty_definitions(): One \"PENALTY_DEFINITION\" line was found inside a PENALTY_DEFINITION...END_PENALTY_DEFINITION block." );
			if ( oneword == "END_PENALTY_DEFINITION" ) {
				parse_a_penalty_definition( lines_subset ); //Once we reach the end, we parse the set of lines.
				in_a_block=false;
			} else {
				lines_subset.push_back( lines[i] ); //If we haven't yet reached the end, add this line to the set of lines to parse.
			}
		} //if(!in_a_block)
	} //Loop through all lines

	return;
}

/// @brief Parse out penalty definition from a single block of lines from file.
/// @details The lines vector should be the lines BETWEEN the PENALTY_DEFINITION and END_PENALTY_DEFINITION lines.
void AACompositionEnergySetup::parse_a_penalty_definition( utility::vector1 < std::string > const &lines ) {

	// Bools for whether we've found all the lines we're looking for:
	bool typefound(false);
	bool nottypefound(false);
	bool propertiesfound(false);
	bool notpropertiesfound(false);
	bool orpropertiesfound(false);
	bool deltastartfound(false);
	bool deltaendfound(false);
	bool penaltiesfound(false);
	bool fractionfound(false);
	bool absolutefound(false);
	signed long deltastart(0);
	signed long deltaend(0);
	TailFunction beforefxn( tf_unknown ); //TailFunction is an enum defined in AACompositionEnergySetup.hh
	TailFunction afterfxn( tf_unknown );

	core::Size const nlines( lines.size() ); //Number of lines we'll be going through.

	//Temporary storage for lists of types.
	utility::vector1 < std::string > types_list;
	utility::vector1 < std::string > not_types_list;

	//Temporary storage for lists of properties.
	utility::vector1 < std::string > properties_list;
	utility::vector1 < std::string > or_properties_list;
	utility::vector1 < std::string > not_properties_list;

	//Temporary storage for the penalties.
	utility::vector1 < core::Real > penalties_vector;

	//Temporary storage for the fraction
	core::Real fraction(0.0);

	//Temporary storage for the absolute value
	signed long absolute(0);

	for ( core::Size i=1; i<=nlines; ++i ) { //Loop through all lines
		std::istringstream curline(lines[i]);
		std::string oneword("");

		curline >> oneword;
		if ( curline.fail() ) continue; //A blank line.

		if ( oneword == "TYPE" ) {
			runtime_assert_string_msg( !typefound, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  A \"TYPE\" line can only be present only once per \"PENALTY_DEFINITION\" block." );
			bool at_least_one(false);
			while ( !curline.eof() ) {
				curline >> oneword;
				runtime_assert_string_msg( !curline.fail() && oneword.length()<=3 && oneword.length()>=1, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): One or more three-letter residue type codes must be present after \"TYPE\" on a \"TYPE\" line." );
				at_least_one=true;
				types_list.push_back( oneword );
			}
			runtime_assert_string_msg( at_least_one, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  A \"TYPE\" line needs at least one three-letter code following \"TYPE\"." );
			typefound=true;
		} else if ( oneword == "NOT_TYPE" ) {
			runtime_assert_string_msg( !nottypefound, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  A \"NOT_TYPE\" line can only be present only once per \"PENALTY_DEFINITION\" block." );
			bool at_least_one(false);
			while ( !curline.eof() ) {
				curline >> oneword;
				runtime_assert_string_msg( !curline.fail() && oneword.length()<=3 && oneword.length()>=1, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): One or more three-letter residue type codes must be present after \"NOT_TYPE\" on a \"NOT_TYPE\" line." );
				at_least_one=true;
				not_types_list.push_back( oneword );
			}
			runtime_assert_string_msg( at_least_one, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  A \"NOT_TYPE\" line needs at least one three-letter code following \"NOT_TYPE\"." );
			typefound=true;
		} else if ( oneword == "PROPERTIES" ) {
			runtime_assert_string_msg( !propertiesfound, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  A \"PROPERTIES\" line can only be present only once per \"PENALTY_DEFINITION\" block." );
			bool at_least_one(false);
			while ( !curline.eof() ) {
				curline >> oneword;
				runtime_assert_string_msg( !curline.fail(), "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Could not parse \"PROPERTIES\" line." );
				at_least_one=true;
				properties_list.push_back( oneword );
			}
			runtime_assert_string_msg( at_least_one, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  A \"PROPERTIES\" line needs at least one property following \"PROPERTIES\"." );
			propertiesfound=true;
		} else if ( oneword == "OR_PROPERTIES" ) {
			runtime_assert_string_msg( !orpropertiesfound, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  An \"OR_PROPERTIES\" line can only be present only once per \"PENALTY_DEFINITION\" block." );
			bool at_least_one(false);
			while ( !curline.eof() ) {
				curline >> oneword;
				runtime_assert_string_msg( !curline.fail(), "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Could not parse \"OR_PROPERTIES\" line." );
				at_least_one=true;
				or_properties_list.push_back( oneword );
			}
			runtime_assert_string_msg( at_least_one, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  An \"OR_PROPERTIES\" line needs at least one property following \"OR_PROPERTIES\"." );
			orpropertiesfound=true;
		} else if ( oneword == "NOT_PROPERTIES" ) {
			runtime_assert_string_msg( !notpropertiesfound, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  A \"NOT_PROPERTIES\" line can only be present only once per \"PENALTY_DEFINITION\" block." );
			bool at_least_one(false);
			while ( !curline.eof() ) {
				curline >> oneword;
				runtime_assert_string_msg( !curline.fail(), "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Could not parse \"NOT_PROPERTIES\" line." );
				at_least_one=true;
				not_properties_list.push_back( oneword );
			}
			runtime_assert_string_msg( at_least_one, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  A \"NOT_PROPERTIES\" line needs at least one property following \"NOT_PROPERTIES\"." );
			notpropertiesfound=true;
		} else if ( oneword == "DELTA_START" ) {
			runtime_assert_string_msg( !deltastartfound, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  A \"DELTA_START\" line can only be present only once per \"PENALTY_DEFINITION\" block." );
			curline >> deltastart;
			runtime_assert_string_msg( !curline.fail(), "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Could not parse \"DELTA_START\" line." );
			deltastartfound=true;
		} else if ( oneword == "DELTA_END" ) {
			runtime_assert_string_msg( !deltaendfound, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  A \"DELTA_END\" line can only be present only once per \"PENALTY_DEFINITION\" block." );
			curline >> deltaend;
			runtime_assert_string_msg( !curline.fail(), "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Could not parse \"DELTA_END\" line." );
			deltaendfound=true;
		} else if ( oneword == "PENALTIES" ) {
			runtime_assert_string_msg( !penaltiesfound, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  A \"PENALTIES\" line can only be present only once per \"PENALTY_DEFINITION\" block." );
			core::Real curval;
			bool at_least_one(false);
			while ( !curline.eof() ) {
				curline >> curval;
				runtime_assert_string_msg( !curline.fail(), "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Could not parse \"PENALTIES\" line." );
				penalties_vector.push_back( curval );
				at_least_one=true;
			}
			runtime_assert_string_msg( at_least_one, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): a \"PENALTIES\" line needs at leas one value." );
			penaltiesfound=true;
		} else if ( oneword == "FRACTION" ) {
			runtime_assert_string_msg( !absolutefound, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  A \"FRACTION\" line cannot be present alongside an \"ABSOLUTE\" line in a \"PENALTY_DEFINITION\" block." );
			runtime_assert_string_msg( !fractionfound, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  A \"FRACTION\" line can only be present only once per \"PENALTY_DEFINITION\" block." );
			curline >> fraction;
			runtime_assert_string_msg( !curline.fail() && fraction >= 0.0, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Could not parse \"FRACTION\" line.  (Note that negative values are not permitted)." );
			absolute=-1; //Negative indicates that fraction should be used.
			fractionfound=true;
		} else if ( oneword == "ABSOLUTE" ) {
			runtime_assert_string_msg( !fractionfound, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  A \"FRACTION\" line cannot be present alongside an \"ABSOLUTE\" line in a \"PENALTY_DEFINITION\" block." );
			runtime_assert_string_msg( !absolutefound, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  An \"ABSOLUTE\" line can only be present only once per \"PENALTY_DEFINITION\" block." );
			curline >> absolute; //Positive value means that absolute should be used
			runtime_assert_string_msg( !curline.fail() && absolute >= 0, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Could not parse \"ABSOLUTE\" line.  (Note that negative values are not permitted)." );
			fraction=0.0;
			absolutefound=true;
		} else if ( oneword == "BEFORE_FUNCTION" ) {
			runtime_assert_string_msg( beforefxn == tf_unknown, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  More than one \"BEFORE_FUNCTION\" statement was found in a \"PENALTY_DEFINITION\" block." );
			std::string namestring(""); //Temporary storage for the name of the before function.
			curline >> namestring;
			runtime_assert_string_msg( !curline.fail() && namestring!="", "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Could not parse \"BEFORE_FUNCTION\" line." );
			beforefxn = get_tailfunction_from_name( namestring );
			runtime_assert_string_msg( beforefxn != tf_unknown, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Could not parse \"BEFORE_FUNCTION\" line." );
		} else if ( oneword == "AFTER_FUNCTION" ) {
			runtime_assert_string_msg( afterfxn == tf_unknown, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition():  More than one \"AFTER_FUNCTION\" statement was found in a \"PENALTY_DEFINITION\" block." );
			std::string namestring(""); //Temporary storage for the name of the after function.
			curline >> namestring;
			runtime_assert_string_msg( !curline.fail() && namestring!="", "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Could not parse \"AFTER_FUNCTION\" line." );
			afterfxn = get_tailfunction_from_name( namestring );
			runtime_assert_string_msg( afterfxn != tf_unknown, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Could not parse \"AFTER_FUNCTION\" line." );
		}
	} //End loop through all lines

	runtime_assert_string_msg( deltastartfound && deltaendfound, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Each \"PENALTY_DEFINITION\" block needs to have a \"DELTA_START\" line and a \"DELTA_END\" line." );
	runtime_assert_string_msg( penaltiesfound, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Each \"PENALTY_DEFINITION\" block needs to have a \"PENALTIES\" line." );
	runtime_assert_string_msg( (fractionfound || absolutefound) && !(fractionfound && absolutefound), "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::parse_a_penalty_definition(): Each \"PENALTY_DEFINITION\" block needs to have a \"FRACTION\" line OR an \"ABSOLUTE\" line (but not both)." );

	// If a BEFORE_FUNCTION line was not provided, default to quadratic:
	if ( beforefxn == tf_unknown ) { beforefxn = tf_quadratic;}
	// If an AFTER_FUNCTION line was not provided, default to quadratic:
	if ( afterfxn == tf_unknown ) { afterfxn = tf_quadratic;}

	//Store everything that we've parsed:
	property_deviation_ranges_.push_back( std::pair<signed long, signed long>(deltastart, deltaend) );
	AACompositionPropertiesSetOP new_property_set( new AACompositionPropertiesSet( types_list, not_types_list, properties_list, or_properties_list, not_properties_list ) ); //Create a new properties set for the properties.
	property_sets_.push_back( new_property_set );
	property_penalties_.push_back( penalties_vector );
	expected_by_properties_fraction_.push_back( fraction );
	expected_by_properties_absolute_.push_back( absolute );
	property_tailfunctions_.push_back( std::pair< TailFunction, TailFunction >(beforefxn, afterfxn) );

	return;
}

/// @brief Do some final checks to ensure that data were loaded properly.
///
void AACompositionEnergySetup::check_data() const {
	if ( TR.Debug.visible() ) TR.Debug << "Checking loaded data." << std::endl;

	core::Size const nproperties( property_sets_.size() );
	runtime_assert_string_msg( property_penalties_.size() == nproperties, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::check_data(): Penalty data were not found for all property sets." );
	runtime_assert_string_msg( property_deviation_ranges_.size() == nproperties, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::check_data(): Deviation range data were not found for all property sets." );
	runtime_assert_string_msg( expected_by_properties_fraction_.size() == nproperties, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::check_data(): Expected fraction data were not found for all property sets." );
	runtime_assert_string_msg( expected_by_properties_absolute_.size() == nproperties, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::check_data(): Absolute expected number data were not found for all property sets." );
	runtime_assert_string_msg( property_tailfunctions_.size() == nproperties, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::check_data(): Tail function data were not found for all property sets." );

	for ( core::Size i=1; i<=nproperties; ++i ) {
		runtime_assert_string_msg( property_deviation_ranges_[i].first <= property_deviation_ranges_[i].second, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::check_data():  The delta range min must be less than the delta range max for each property set.");
		runtime_assert_string_msg( static_cast< signed long >( property_penalties_[i].size() ) == property_deviation_ranges_[i].second-property_deviation_ranges_[i].first+1, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::check_data():  Penalties must be provided for every delta in the range from DELTA_START to DELTA_END.  Too many or too few were found.");
		runtime_assert_string_msg( property_penalties_[i].size() >=2, "Error in core::scoring::aa_composition_energy::AACompositionEnergySetup::check_data():  At least two penalty values must be specified.  Too few were found." );
	}

	if ( TR.Debug.visible() ) TR.Debug << "Data checks passed!" << std::endl;
	return;
}

} // aa_composition_energy
} // scoring
} // core
