// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A class that contains orbital parameters
///
/// @details
/// This class contains the "chemical" information for orbitals. This does not contain the actual
/// xyz coordinates of the class which is managed by the conformation/Residue.hh. The orbital_type properties
/// are assigned by the class OrbitalTypeSet which is initiated from the ChemicalManager. Orbital type properties
/// are currently are read in from the file located chemical/orbital_type_sets/fa_standard/orbital_properties.txt.
/// These properties contain the the parameters of distance, but can be modified. Currently this is a very small
/// class that will be added on as more and more properties are identified and added.
/// Note that information about the atomtype is stored along with the orbital type. This may or may not be useful
/// later. Just adding the functionality for shits and giggles.
///
/// Orbital type name: the orbital type name contains the hybridization, orbital name, and element associated with the orbital
///
/// Hybridization: the hybridiztion of the atom that the orbital is bonded to (sp, sp2, sp3)
///
/// Orbital Name: the name of the orbital. This usually is p, d, pi, sigma. The orbital name is different than the orbital type name
///
/// Atom Type Name: the type of atom associated with an orbital
///
/// Distance: distance the orbital comes off of the atom. Currently, for residues, the distance is the Bohr radius of Hydrogen+element
///
///
///
///
///
/// @author Steven Combs
///
///
/////////////////////////////////////////////////////////////////////////
#ifndef INCLUDED_core_chemical_orbitals_OrbitalType_hh
#define INCLUDED_core_chemical_orbitals_OrbitalType_hh

#include <utility/vector1.hh>
#include <core/types.hh>
#include <string>
#include <core/chemical/orbitals/OrbitalTypeMapper.fwd.hh>

namespace core{
namespace chemical{
namespace orbitals{

class OrbitalType {
public:


	OrbitalType(std::string & orbital_name, std::string & atom_type_name);

	/// @brief The parameters are the actual headings in the orbital_properties.txt. If you want to add more paramters,
	/// you must edit orbital_properties.txt and add another heading. You also need to edit AtomTypeSet.txt so that
	/// it recognizes that parameter and parses it. The parameters are different form the properties in that they
	/// are Reals/Size and properties are strings.
	void
	set_parameter(
		std::string const & param,
		core::Real const setting
	);

	/// @brief Currently, these properties are not actually in the orbital_properties.txt. I have them here
	/// as an example on how to add properties. This is also a place holder as the ligand code will
	/// soon be using these properties. The Acceptor/Donor could refer to orbitals that have a lone pair
	/// and are donating to a hydrogen, or an electron defficient region. In order to add properties, one
	/// must add the properties to the last line of orbital_properties.txt and make a private member variable
	/// for that property in the header file. Then do a string match comparision, like seen below. These properties
	/// are set via OrbitalTypeSet.hh
	void
	set_property(
		std::string const & property,
		bool const setting
	);

	/// @brief returns the name of the orbital type. defined in orbital_properties.txt
	std::string name() const;

	/// @brief returns the distance from the atom the orbital comes off. defined in orbital_properties.txt
	Real distance() const;

	/// @brief returns the atom_types associated with the orbital type. defined in orbital_properties.txt
	utility::vector1<std::string> atom_type_name() const;

	/// @brief returns hybrdiziation of atom the orbital is attached to
	std::string hybridization() const;

	/// @brief returns the orbital associated with the type
	std::string orbital_name() const;

	orbital_type_enum orbital_enum() const;


private:
	/// @brief is the orbital an acceptor?
	bool is_orbital_acceptor_;

	/// @brief is the orbital a donor?
	bool is_orbital_donor_;

	/// @brief the orbital type name
	std::string orbital_type_name_;

	/// @brief hybridization of atom the orbital is attached to
	std::string hybridization_;

	/// @brief orbital name associated with the orbital type
	std::string orbital_name_;

	/// @brief atom type associated with a given orbital type
	utility::vector1<std::string> atom_type_name_;


	/// @brief distance of orbital from center of atom
	core::Real distance_;

	orbital_type_enum orbital_type_enum_;


};

}
}
}


#endif /* ORBITALTYPE_HH_ */
