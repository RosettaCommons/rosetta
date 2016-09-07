// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
// This file is based off of code from the Biochemistry Library (BCL).
// The BCL is copyright Vanderbilt University (Meiler Lab), a RosettaCommons member

/// @file   core/chemical/Element.hh
/// @brief  The data for the element types
/// @author Rosetta conversion: Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_Element_hh
#define INCLUDED_core_chemical_Element_hh

#include <core/chemical/Element.fwd.hh>
#include <core/chemical/ElectronConfiguration.hh>
#include <core/chemical/Elements.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector0.hh>

#include <iostream>
#include <string>

namespace core {
namespace chemical {


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!
//! @class Element
//! @brief stores element properties
//! @details This is a low level helper class to store element properties
//!
//! @see @link example_chemistry_element_type_data.cpp @endlink
//! @author meilerj, woetzen, mendenjl
//! @date 08/31/2005
//!
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Element : public utility::pointer::ReferenceCount
{

public:

	///////////
	// Enums //
	///////////

	//! enum  properties for element types
	enum Properties
	{
		Mass,                //!< Mass
		CovalentRadius,      //!< CovalentRadius
		VDWaalsRadius,       //!< VdWaalsRadius
		NumberOfProperties   //!< Number of properties
	};

	//! @brief element type property as string
	//! @param PROPERTY the property desired
	//! @return the property as string
	static const std::string &get_property_name( const Properties &PROPERTY);


public:


	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	//! @brief construct undefined element
	Element();

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
	Element
	(
		const core::Size ATOMIC_NUMBER,
		const core::Size PERIOD,
		const core::Size MAIN_GROUP,
		std::string CHEMICAL_SYMBOL, // moving
		std::string CHEMICAL_NAME, // moving
		const core::Real MASS,
		const core::Real COVALENT_RADIUS,
		const core::Real VDW_RADIUS,
		ElectronConfiguration ELECTRON_CONFIGURATION
	);

	//! @brief virtual copy constructor
	ElementOP Clone() const;

	/////////////////
	// data access //
	/////////////////

	//! @brief The element enumeration
	core::chemical::element::Elements element() const {
		return element_;
	}

	//! @brief atomic number
	//! @return atomic number
	core::Size get_atomic_number() const
	{
		return atomic_number_;
	}

	//! @return Period
	core::Size get_period() const
	{
		return period_;
	}

	//! @return main Group #
	core::Size get_main_group() const
	{
		return main_group_;
	}

	//! @brief GetChemicalSymbol
	//! @return chemical symbol one or two letters as AtomName
	const std::string &get_chemical_symbol() const
	{
		return chemical_symbol_;
	}

	//! @brief GetChemicalName
	//! @return full chemical name
	const std::string &get_chemical_name() const
	{
		return chemical_name_;
	}

	//! @brief element type property as core::Real
	//! @param PROPERTY the property desired
	//! @return the property as core::Real
	core::Real get_property( const Element::Properties &PROPERTY) const
	{
		return properties_[ PROPERTY];
	}

	//! @brief electron configuration
	//! @return the ElectronConfiguration
	const ElectronConfiguration &get_electron_configuration() const
	{
		return electron_configuration_;
	}


	//! @brief tell whether this element type can participate in a conjugated system
	//! @return true if this element can participate in a common conjugated system
	//! Specifically tests if the element has 1-4 valence electrons in P orbitals
	bool is_conjugatable() const
	{
		const core::Size n_valence_p( electron_configuration_.valence_electrons_p());

		return n_valence_p > 0 && n_valence_p < 5;
	}


	///This is legacy code from old element set
	/// @brief Return the full name of the Element
	Real weight() const { return properties_[Mass]; }


	/// @brief Return true unless the element actually exists in the periodic table.
	bool is_fake() const { return atomic_number_ == 0; }

	//////////////////////
	// input and output //
	//////////////////////

public:


	//! @brief read from std::istream
	//! @param ISTREAM input stream
	//! @return istream which was read from
	std::istream &read( std::istream &ISTREAM);

	//! @brief write to std::ostream
	//! @param OSTREAM output stream
	//! @return ostream which was written to
	std::ostream &write( std::ostream &OSTREAM) const;


private:

	//////////
	// data //
	//////////

	core::chemical::element::Elements element_;                       //!< Element enum
	core::Size             atomic_number_;                            //!< atomic number
	core::Size             period_;                                   //!< Period
	core::Size             main_group_;                               //!< Group # in the main group (1-8) such as transistion metals
	std::string            chemical_symbol_;                          //!< ChemicalSymbol
	std::string            chemical_name_;                            //!< ChemicalName
	ElectronConfiguration  electron_configuration_;                   //!< electron configuration
	utility::vector0< core::Real > properties_; //!< real-valued properties


}; // class Element

inline std::ostream &
operator<< (std::ostream & out, Element const & obj ){
	return obj.write( out );
}

inline std::istream &
operator>> (std::istream & in, Element & obj ){
	return obj.read( in );
}


} // namespace core
} // namespace chemical

#endif
