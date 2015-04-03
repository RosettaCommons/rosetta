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
/// Donor: does the orbital donate electrons? currently not implemented
///
/// Acceptor: is the orbital accept electrons? currently not implemented
///
/// @author Steven Combs
///
///
/////////////////////////////////////////////////////////////////////////
#include <core/chemical/orbitals/OrbitalType.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <vector>
#include <core/chemical/orbitals/OrbitalTypeMapper.hh>
#include <core/types.hh>

#include <utility/vector1.hh>


namespace core{
namespace chemical{
namespace orbitals{

//mjo commenting out 'atom_type_name' because it is not used and cause a warnings
/// @brief Constructor that is generally initialized in OrbitalTypeSet.hh. If you add a property,
/// you must initialize that property as false!!!!!!!!
OrbitalType::OrbitalType(std::string & orbital_name, std::string & /*atom_type_name*/):
			is_orbital_acceptor_(false),
			is_orbital_donor_(false),
			distance_(0.0)

{


	utility::vector1<std::string> name_split(utility::string_split(orbital_name, '.'));
	orbital_type_name_= orbital_name;
	orbital_type_enum_ = OrbitalTypeMapper::get_instance()->get_orbital_enum(orbital_name);
	orbital_name_= name_split[2];
	hybridization_= name_split[3];

	name_split=utility::string_split(orbital_name, '_');
	for(core::Size i=1; i <= name_split.size(); ++i){
		atom_type_name_.push_back(name_split[i]);
	}

}


/// @brief The parameters are the actual headings in the orbital_properties.txt. If you want to add more paramters,
/// you must edit orbital_properties.txt and add another heading. You also need to edit AtomTypeSet.txt so that
/// it recognizes that parameter and parses it. The parameters are different form the properties in that they
/// are Reals/Size and properties are strings.
void
OrbitalType::set_parameter(
	std::string const & param,
	core::Real const setting
)
{
	if ( param == "distance" ) {
		distance_ = setting;
	}  else {
		utility_exit_with_message( "unrecognized orbital_type parameter "+param );
	}
}


/// @brief Currently, these properties are not actually in the orbital_properties.txt. I have them here
/// as an example on how to add properties. This is also a place holder as the ligand code will
/// soon be using these properties. The Acceptor/Donor could refer to orbitals that have a lone pair
/// and are donating to a hydrogen, or an electron defficient region. In order to add properties, one
/// must add the properties to the last line of orbital_properties.txt and make a private member variable
/// for that property in the header file. Then do a string match comparision, like seen below. These properties
/// are set via OrbitalTypeSet.hh
void
OrbitalType::set_property(
	std::string const & property,
	bool const setting
)
{
	if(property == "DONOR"){
		is_orbital_donor_=setting;
	}else if(property == "ACCEPTOR"){
		is_orbital_acceptor_=setting;
	} else {
	utility_exit_with_message( "unrecognized atomtype property "+property );
}


}


/// @brief returns the name of the orbital type. defined in orbital_properties.txt
std::string OrbitalType::name() const
{
	return orbital_type_name_;
}

orbital_type_enum OrbitalType::orbital_enum() const
{
	return orbital_type_enum_;
}


/// @brief returns the distance from the atom the orbital comes off. defined in orbital_properties.txt
Real OrbitalType::distance() const
{
	return distance_;
}

/// @brief returns the atom_types associated with the orbital type. defined in orbital_properties.txt
utility::vector1<std::string> OrbitalType::atom_type_name() const
{
	return atom_type_name_;
}

/// @brief returns hybrdiziation of atom the orbital is attached to
std::string OrbitalType::hybridization() const
{
	return hybridization_;
}

/// @brief returns the orbital associated with the type
std::string OrbitalType::orbital_name() const
{
	return orbital_name_;
}


}
}
}

