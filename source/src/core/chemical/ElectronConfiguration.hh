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

/// @file   core/chemical/gasteiger/ElectronConfiguration.hh
/// @brief  The data for the BCL electron configuration
/// @author To Rosetta transitioning: Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_ElectronConfiguration_hh
#define INCLUDED_core_chemical_ElectronConfiguration_hh

#include <core/types.hh>

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

namespace core {
namespace chemical {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!
//! @class ElectronConfiguration
//! @brief describes the electron configuration of atoms
//! @details Describes the electron configuration of an atom on quantum chemistry level.
//!
//! @see @link example_chemistry_electron_configuration.cpp @endlink
//! @author woetzen, karakam
//! @date 10/28/2007
//!
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ElectronConfiguration {
public:

	enum PrincipalQuantumNumber
	{
		e_1,
		e_2,
		e_3,
		e_4,
		e_5,
		e_6,
		e_7,
		MaxPrincipleQuantumNumber
	};

	static std::vector<std::string> const & PrincipalQuantumNumber_strings();

	enum AngularMomentumQuantumNumber
	{
		e_S,
		e_P,
		e_D,
		e_F,
		MaxAngularMomentumQuantumNumber
	};

	static std::vector<std::string> const & AngularMomentumQuantumNumber_strings();

	//! @brief PrincipalQuantumNumber as string
	//! @param NUM the PrincipalQuantumNumber desired
	//! @return the PrincipalQuantumNumber as string
	static const std::string &get_descriptor( const PrincipalQuantumNumber &NUM);

	//! @brief PrincipalQuantumNumber from string
	//! @return NUM the PrincipalQuantumNumber desired
	//! @param the PrincipalQuantumNumber as string
	static PrincipalQuantumNumber get_principal_quantum_number( std::string const &STR );

	//! @brief AngularMomentumQuantumNumber as string
	//! @param NUM the AngularMomentumQuantumNumber desired
	//! @return the AngularMomentumQuantumNumber as string
	static const std::string &get_descriptor( const AngularMomentumQuantumNumber &NUM);

	//! @brief AngularMomentumQuantumNumber as string
	//! @return NUM the AngularMomentumQuantumNumber desired
	//! @param the AngularMomentumQuantumNumber as string
	static AngularMomentumQuantumNumber get_angular_momentum_quantum_number( std::string const &STR );

	//#      //! PrincipalQuantumNumberEnum is used for I/O of PrincipalQuantumNumber
	//#      typedef util::WrapperEnum< PrincipalQuantumNumber, &get_descriptor, MaxPrincipleQuantumNumber>
	//#                PrincipalQuantumNumberEnum;
	//#
	//#      //! AngularMomentumQuantumNumberEnum is used for I/O of AngularMomentumQuantumNumbers
	//#      typedef util::WrapperEnum< AngularMomentumQuantumNumber, &get_descriptor, MaxAngularMomentumQuantumNumber>
	//#                AngularMomentumQuantumNumberEnum;


public:

	//#      //! single instance of that class
	//#      static const util::SiPtr< const util::ObjectInterface> s_Instance;

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	//! default constructor
	ElectronConfiguration();

	//! construct from actual number of electrons
	ElectronConfiguration
	(
		const core::Size VALENCE_ELECTRONS_SP,
		const core::Size VALENCE_ELECTRONS_SPD,
		const core::Size NUMBER_ELECTRONS_1S,
		const core::Size NUMBER_ELECTRONS_1P,
		const core::Size NUMBER_ELECTRONS_1D,
		const core::Size NUMBER_ELECTRONS_1F,
		const core::Size NUMBER_ELECTRONS_2S,
		const core::Size NUMBER_ELECTRONS_2P,
		const core::Size NUMBER_ELECTRONS_2D,
		const core::Size NUMBER_ELECTRONS_2F,
		const core::Size NUMBER_ELECTRONS_3S,
		const core::Size NUMBER_ELECTRONS_3P,
		const core::Size NUMBER_ELECTRONS_3D,
		const core::Size NUMBER_ELECTRONS_3F,
		const core::Size NUMBER_ELECTRONS_4S,
		const core::Size NUMBER_ELECTRONS_4P,
		const core::Size NUMBER_ELECTRONS_4D,
		const core::Size NUMBER_ELECTRONS_4F,
		const core::Size NUMBER_ELECTRONS_5S,
		const core::Size NUMBER_ELECTRONS_5P,
		const core::Size NUMBER_ELECTRONS_5D,
		const core::Size NUMBER_ELECTRONS_5F,
		const core::Size NUMBER_ELECTRONS_6S,
		const core::Size NUMBER_ELECTRONS_6P,
		const core::Size NUMBER_ELECTRONS_6D,
		const core::Size NUMBER_ELECTRONS_6F,
		const core::Size NUMBER_ELECTRONS_7S,
		const core::Size NUMBER_ELECTRONS_7P,
		const core::Size NUMBER_ELECTRONS_7D,
		const core::Size NUMBER_ELECTRONS_7F
	);

	//#     //! @brief virtual copy constructor
	//#     ElectronConfiguration *Clone() const;

	/////////////////
	// data access //
	/////////////////

	//#      //! @brief returns class name
	//#      //! @return the class name as const ref std::string
	//#      const std::string &GetClassIdentifier() const;

	//! return number ValenceElectrons in the sigma valence orbitals
	core::Size valence_electrons_s() const
	{
		return electrons_[ valence_quantum_number_][ e_S];
	}

	//! @brief return the number of valence electrons in SP orbitals
	core::Size unpaired_valence_electrons_sp() const
	{
		return std::min( valence_electrons_sp_, max_valence_electrons_sp() - valence_electrons_sp_);
	}

	//! @return number ValenceElectrons in the pi valence orbitals
	core::Size valence_electrons_p() const
	{
		return electrons_[ valence_quantum_number_][ e_P];
	}

	//! @return number valence_electrons_sp
	core::Size valence_electrons_sp() const;

	//! @return the maximum number of electrons in SP orbitals for the noble gas in this period
	core::Size max_valence_electrons_sp() const;

	//! @return number valence_electrons_spd
	core::Size valence_electrons_spd() const;

	///////////////
	// operators //
	///////////////

	//! @brief number of electrons in that orbital
	//! @param PRINCIPAL_QUANTUM_NUMBER 1, 2, 3, 4, 5, 6. or 7
	//! @param ANGULAR_MOMENTUM_QUANTUM_NUMBER S, P, D, or F
	//! @return number of electrons in that particular orbital indicated by PRINCIPAL_QUANTUM_NUMBER and ANGULAR_MOMENTUM_QUANTUM_NUMBER
	core::Size operator ()
	(
		const PrincipalQuantumNumber PRINCIPAL_QUANTUM_NUMBER,
		const AngularMomentumQuantumNumber ANGULAR_MOMENTUM_QUANTUM_NUMBER
	) const
	{
		return electrons_[ PRINCIPAL_QUANTUM_NUMBER][ ANGULAR_MOMENTUM_QUANTUM_NUMBER];
	}

	//! @brief number of electrons in the valence orbital
	//! @param ANGULAR_MOMENTUM_QUANTUM_NUMBER S, P, D, or F
	//! @return number of electrons in the valence orbital indicated by ANGULAR_MOMENTUM_QUANTUM_NUMBER
	core::Size operator ()
	(
		const AngularMomentumQuantumNumber ANGULAR_MOMENTUM_QUANTUM_NUMBER
	) const
	{
		return electrons_[ valence_quantum_number_][ ANGULAR_MOMENTUM_QUANTUM_NUMBER];
	}

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
	//! @param INDENT number of indentations
	//! @return ostream which was written to
	std::ostream &write( std::ostream &OSTREAM ) const;


private:

	//////////
	// data //
	//////////

	core::Size valence_electrons_sp_;     //!< number S and P valence electrons (main group 1,...,8)
	core::Size valence_electrons_spd_;    //!< number S, P, and D valence electrons (group 1,...,18)
	core::Size electrons_[ MaxPrincipleQuantumNumber][ 4];      //!< ElectronConfiguration
	core::Size valence_quantum_number_;   //!< last quantum number with a non-zero number of electrons


	// core::Size [ 7][ 4];
	static std::vector< std::vector<core::Size> > const & s_MaxElectronsInOrbital(); //!< maximum electron for each orbital


}; // ElectronConfiguration

inline std::ostream &
operator<< (std::ostream & out, ElectronConfiguration const & obj ){
	return obj.write( out );
}

inline std::istream &
operator>> (std::istream & in, ElectronConfiguration & obj ){
	return obj.read( in );
}

} // namespace core
} // namespace chemical

#endif
