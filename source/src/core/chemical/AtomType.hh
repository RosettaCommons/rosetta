// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
/// @file AtomType.hh
///
/// @brief
/// A class for defining atom parameters, known as atom_types
///
/// @details
/// This class contains the "chemical" information for atoms. This does not contain the actual
/// xyz coordinates of the class (xyz found in core/conformation/Atom.hh. The atom_type properties
/// are assigned by the class AtomTypeSet which is initiated from the ChemicalManager. Atom type properties
/// are currently are read in from the file located chemical/atom_type_sets/fa_standard/atom_properties.txt.
/// These properties contain the the properties of LJ_RADIUS, LJ_WDEPTH, LK_DGRFREE, LK_LAMBDA, LK_VOLUME.
/// These properties are used in the scoring function fa_atr, fa_rep, fa_sol, which is located in the Etable
/// (core/scoring/etable/Etable.hh)
/// Additional parameters are acceptor/donor, hybridization, and orbital parameters.
///
/// @author Phil Bradley
/// @author Steven Combs - comments
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_core_chemical_AtomType_hh
#define INCLUDED_core_chemical_AtomType_hh


// Unit headers
#include <core/chemical/AtomType.fwd.hh>
#include <core/types.hh>

// Package headers
#include <core/chemical/types.hh>

// Utility headers
#include <utility>
#include <utility/vector1.hh>

// C++ headers
#include <string>

#ifdef    SERIALIZATION
#include <utility/serialization/serialization.hh>
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

/// @brief Maximum distance between a heavy atom and a hydrogen atom
/// to which it is chemically bound Set in .cc file.
extern Real const MAX_CHEMICAL_BOND_TO_HYDROGEN_LENGTH;


/// @brief basic atom type
///
/// @details name, element, certain properties and parameters
class AtomType {

public:

	// Can't make this private because then the ResidueType -- which
	// holds a vector of AtomTypes -- can't be serialized.
	AtomType() {} // for serialization


	/// @brief Construct a new atom type with its name and element.
	///
	/// @details All its properties are unset by default.
	///
	AtomType( std::string const & name_in, std::string  element_in ):
		name_( name_in ),
		element_(std::move( element_in )),
		lj_radius_( 0 ),
		lj_wdepth_( 0 ),
		lk_lambda_( 0 ),
		lk_volume_( 0 ),
		lk_dgfree_( 0 ),
		is_acceptor_( false ),
		is_donor_( false ),
		is_polar_hydrogen_( false ),
		is_h2o_( name_in == "HOH" ),
		is_aromatic_( false ),
		atom_has_orbitals_(false),
		atom_is_virtual_(false),
		atom_is_repulsive_(false),
		hybridization_( UNKNOWN_HYBRID )
	{}

	AtomType(AtomType const & ) = default;

	void
	print( std::ostream & out ) const;

	friend
	std::ostream &
	operator<< ( std::ostream & out, const AtomType & atom_type );

	/// @brief Lazaridis and Karplus solvation parameter -- lambda
	Real
	lk_lambda() const
	{
		return lk_lambda_;
	}

	/// @brief Lazaridis and Karplus solvation parameter -- dgfree
	Real
	lk_dgfree() const
	{
		return lk_dgfree_;
	}

	/// @brief Lazaridis and Karplus solvation parameter -- volume
	Real
	lk_volume() const
	{
		return lk_volume_;
	}

	/// @brief Lennard-Jones 6-12 potential parameter -- atom radius
	///
	/// @details
	///   There are two functionally identical versions of the Lennard-Jones potential:
	///
	///   E ~ 4eps(sigma1/d)^12 - (sigma1/d)^6 and
	///      and
	///   E ~ eps(sigma2/d)^12 - 2*(sigma2/d)^6
	///
	/// where sigma1 and sigma2 represent two different interpretations of the radius. eps represents the depth of the potential well.
	/// Sigma1 represents the distance between the two atoms where the Lennard-Jones energy is 0, i.e. where a collision is just forming/resolving.
	/// Sigma2 represents the distance between the two atoms where the derivative of the Lennard-Jones energy is 0, i.e. the minimum of the well depth.
	///
	/// In rosetta, we mean sigma2 when we talk about radii, but PyMOL usually sets the radii to sigma1.
	///  If you see two atoms overlapping using the Rosetta radii, they're not necessarily in collision. They are just not
	///  at their minimum value.
	///
	/// The distances are related as sigma2 = 2^(1.0/6)*sigma1;  sigma2 ~= 1.22*sigma1
	Real
	lj_radius() const
	{
		return lj_radius_;
	}

	/// @brief Lennard-Jones 6-12 potential parameter -- well depth
	Real
	lj_wdepth() const
	{
		return lj_wdepth_;
	}


	/// @brief whether atom is a hydrogen bond acceptor
	bool
	is_acceptor() const
	{
		return is_acceptor_;
	}

	/// @brief whether atom is a hydrogen bond donor
	bool
	is_donor() const
	{
		return is_donor_;
	}

	/// @brief whether atom is a polar hydrogen atom
	bool
	is_polar_hydrogen() const
	{
		return is_polar_hydrogen_;
	}

	/// @brief whether atom is a hydrogen atom
	bool
	is_hydrogen() const
	{
		return ( element_ == "H" );
	}

	// this is a little silly
	/// @brief whether atom is a heavy atom
	bool
	is_heavyatom() const
	{
		return ( element_ != "H" );
	}

	/// @brief is atom type virtual?
	bool is_virtual() const {
		debug_assert( (name_ == "VIRT") ? atom_is_virtual_ : true ); // Raise an error if an atom type named VIRT is not virtual.
		return atom_is_virtual_;
	}

	/// @brief is atom type repulsive (REPL, REPLS, HREPS)
	bool is_repulsive() const {
		return atom_is_repulsive_;
	}

	/// @brief whether atom is a water
	bool
	is_h2o() const
	{
		return is_h2o_;
	}

	/// @brief whether atom is aromatic
	bool
	is_aromatic() const
	{
		return is_aromatic_;
	}

	/// @brief atom has an orbital attached
	bool
	atom_has_orbital() const
	{
		return atom_has_orbitals_;
	}

	/// @brief is the H atom aromatic?
	bool is_haro() const
	{
		return ( name_ == "Haro" );
	}


	/// @brief set LJ and LK solvation parameter for this atom type
	void
	set_parameter(
		std::string const & param,
		Real const setting
	);


	/// @brief set relevant properties for this atom type hh
	void
	set_property(
		std::string const & property,
		bool const setting
	);

	/// @brief retrieve an atom's hybridization status.
	Hybridization const &
	hybridization() const
	{
		return hybridization_;
	}

	/// @brief set all standard properties to false, set hybridization to
	///UNKNOWN_HYBRID, and clear extra properties
	///
	void
	clear_properties();

	/// @brief set standard property to true, or set the specified hybridization
	void
	add_property(std::string const & property);

	utility::vector1< std::string >
	get_all_properties() const;

	/// @brief returns the one- or two-letter element type
	std::string
	element() const
	{
		return element_;
	}

	///@brief Get the name of this AtomType.  Note that this is NOT the same as the atom name that is written out to PDB.
	/// You want pose.residue(i).atom_name(j) for that.
	std::string
	atom_type_name() const
	{
		return name_;
	}

	/// @brief return an additional, non-hardcoded property
	Real
	extra_parameter( Size const index ) const
	{
		return extra_parameters_[ index ];
	}

	/// @brief return an additional, non-hardcoded property
	void
	set_extra_parameter( Size const index, Real const setting )
	{
		if ( extra_parameters_.size() < index ) extra_parameters_.resize( index, 0.0 );
		extra_parameters_[ index ] = setting;
	}

	/// @set all the extra parameters at once
	void
	set_all_extra_parameters(utility::vector1< Real > const & extra_parameters);

	// used after we clone an atomtype
	void
	name( std::string const & setting )
	{
		name_ = setting;
	}

	///@brief Get the name of this AtomType.  Note that this is NOT the same as the atom name that is written out to PDB.
	/// You want pose.residue(i).atom_name(j) for that.
	std::string const& name() const { return name_; };

	// data
private:

	// name
	std::string name_; // non-const for use in cloning at time of atomtypeset creation...

	// element
	//std::string const element_;
	std::string element_;  // changing this to non-const so copy assignment operator of 'AtomType' could be used in C++11

	// lj
	Real lj_radius_;
	Real lj_wdepth_;

	// lk
	Real lk_lambda_;
	Real lk_volume_;
	Real lk_dgfree_;

	// extras
	utility::vector1< Real > extra_parameters_;

	// props -- false by default
	// doh! gccdebug does not default construct bools to false!
	// so if you add anything here and you want it to be false by default,
	// set it in the constructor!!!!
	bool is_acceptor_;
	bool is_donor_;
	bool is_polar_hydrogen_;
	bool is_h2o_;
	bool is_aromatic_;
	bool atom_has_orbitals_;
	bool atom_is_virtual_;
	bool atom_is_repulsive_;

	Hybridization hybridization_;

#ifdef    SERIALIZATION
public:
	friend class cereal::access;
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

#ifdef    SERIALIZATION
/// @brief Serialize an AtomType
template < class Archive >
void serialize_atom_type( Archive & arc, AtomType const & ptr );

/// @brief Deserialize an AtomType
template < class Archive >
void deserialize_atom_type( Archive & arc, AtomType & ptr );
#endif // SERIALIZATION

} // chemical
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_AtomType )
#endif // SERIALIZATION

#endif // INCLUDED_core_chemical_AtomType_HH
