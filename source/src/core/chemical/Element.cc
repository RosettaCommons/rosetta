// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
// This file is based off of code from the Biochemistry Library (BCL).
// The BCL is copyright Vanderbilt University (Meiler Lab), a RosettaCommons member

/// @file   core/chemical/Element.cc
/// @brief  The data for the element types
/// @author Rosetta conversion: Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/Element.hh>

#include <core/chemical/gasteiger/util.hh>

#include <utility/numbers.hh>

#include <utility/exit.hh>

namespace core {
namespace chemical {


///////////
// Enums //
///////////

//! @brief element type property as string
//! @param PROPERTY the property desired
//! @return the property as string
const std::string &Element::get_property_name( const Element::Properties &PROPERTY)
{
	static const std::string s_properties[] =
	{
			"Mass",
			"CovalentRadius",
			"VDWaalsRadius",
			"Properties" // GetStaticClassName< Properties>()
	};
	return s_properties[ PROPERTY];
}

//////////////////////////////////
// construction and destruction //
//////////////////////////////////

//! @brief construct undefined element
Element::Element() :
		element_( core::chemical::element::UnknownElement ),
		atomic_number_( utility::get_undefined_size() ),
		period_( utility::get_undefined_size() ),
		main_group_( utility::get_undefined_size() ),
		chemical_symbol_( "X"),
		chemical_name_( "UNDEFINED_ELEMENT"),
		electron_configuration_(),
		properties_(NumberOfProperties, utility::get_undefined_real() ) // set all properties to undefined
{}

//! @brief construct element from all its data
//! @param ATOMIC_NUMBER           - number in the PSE
//! @param PERIOD                  - in which period is the element
//! @param MAIN_GROUP              - in which of the main groups does the element belong (0 for transition metals)
//! @param CHEMICAL_SYMBOL         - one or two letters as in international PSE, first letter capital
//! @param CHEMICAL_NAME           - full international name (first letter capital)
//! @param MASS                    - the atomic mass as a weighted avergage of all isotopes
//! @param COVALENT_RADIUS         - radius of atom with electrons
//! @param VDW_RADIUS              - vdw radius
//! @param ELECTRON_CONFIGURATION  - the electron configuration
Element::Element
(
		const core::Size ATOMIC_NUMBER,
		const core::Size PERIOD,
		const core::Size MAIN_GROUP,
		const std::string &CHEMICAL_SYMBOL,
		const std::string &CHEMICAL_NAME,
		const core::Real MASS,
		const core::Real COVALENT_RADIUS,
		const core::Real VDW_RADIUS,
		const ElectronConfiguration &ELECTRON_CONFIGURATION
) :
atomic_number_( ATOMIC_NUMBER),
period_( PERIOD),
main_group_( MAIN_GROUP),
chemical_symbol_( CHEMICAL_SYMBOL),
chemical_name_( CHEMICAL_NAME),
electron_configuration_( ELECTRON_CONFIGURATION),
properties_(NumberOfProperties, utility::get_undefined_real() )
{
		element_ = core::chemical::element::elements_from_name(chemical_symbol_);
		properties_[ Mass] = MASS;
		properties_[ CovalentRadius]      = COVALENT_RADIUS;
		properties_[ VDWaalsRadius]       = VDW_RADIUS;

}

//! @brief virtual copy constructor
ElementOP Element::Clone() const
{
	return new Element( *this );
}


//////////////////////
// input and output //
//////////////////////

//! @brief read from std::istream
//! @param ISTREAM input stream
//! @return istream which was read from
std::istream &Element::read( std::istream &ISTREAM)
{
	// read member
	std::string tag;
	ISTREAM >> tag;
	if ( tag != "Element:" ) {
		utility_exit_with_message( "Malformated elements file. 'Element:' tag expected. '" + tag + "' found.");
	}
	ISTREAM >> chemical_symbol_; //io::Serialize::Read( chemical_symbol_, ISTREAM);
	element_ = core::chemical::element::elements_from_name(chemical_symbol_);
	ISTREAM >> chemical_name_; //io::Serialize::Read( chemical_name_, ISTREAM);
	ISTREAM >> atomic_number_; //io::Serialize::Read( atomic_number_, ISTREAM);
	gasteiger::safe_read( ISTREAM, period_); //io::Serialize::Read( period_, ISTREAM);
	gasteiger::safe_read( ISTREAM, main_group_);

	// Ensure that the number of properties is the same as when the file was written
	core::Size properties_in_files;
	ISTREAM >> properties_in_files; //io::Serialize::Read( properties_in_files, ISTREAM);
	if( properties_in_files != NumberOfProperties ) {
		utility_exit_with_message("Number of properties in Element file was incorrect");
	}

	for( core::Size a = 0; a < core::Size( NumberOfProperties); a++)
	{
		gasteiger::safe_read( ISTREAM, properties_[ a] );
	}

	ISTREAM >> electron_configuration_; //io::Serialize::Read( m_ElectronConfiguration, ISTREAM)

	// end
	return ISTREAM;
}

//! @brief write to std::ostream
//! @param OSTREAM output stream
//! @param INDENT number of indentations
//! @return ostream which was written to
std::ostream &Element::write( std::ostream &OSTREAM) const
{
	// write member
	OSTREAM << "Element:" << ' ';
	OSTREAM << chemical_symbol_ << ' '; //io::Serialize::Write( chemical_symbol_, OSTREAM, INDENT) << '\n';
	OSTREAM << chemical_name_ << ' '; //io::Serialize::Write( chemical_name_, OSTREAM, INDENT) << '\n';
	OSTREAM << atomic_number_ << ' '; //io::Serialize::Write( atomic_number_, OSTREAM, INDENT) << '\n';
	gasteiger::safe_write( OSTREAM, period_ ); //io::Serialize::Write( period_, OSTREAM, INDENT) << '\n';
	gasteiger::safe_write( OSTREAM, main_group_);

	// Write out the number of properties, if this changes, the read function will fail
	OSTREAM << NumberOfProperties << ' '; //io::Serialize::Write( core::Size( NumberOfProperties), OSTREAM, INDENT) << '\n';

	for( core::Size a = 0; a < core::Size( NumberOfProperties); a++)
	{
		gasteiger::safe_write( OSTREAM, properties_[ a]); //io::Serialize::Write( properties_[ a], OSTREAM, INDENT) << '\n';
	}


	OSTREAM << std::endl;
	OSTREAM << electron_configuration_; //io::Serialize::Write( electron_configuration_, OSTREAM, INDENT) << '\n';


	// end
	return OSTREAM;
}


} // namespace core
} // namespace chemical

