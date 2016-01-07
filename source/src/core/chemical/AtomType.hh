// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
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
#include <utility/vector1.hh>

// C++ headers
#include <string>

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

	/// @brief Construct a new atom type with its name and element.
	///
	/// @details All its properties are unset by default.
	///
	AtomType( std::string const & name_in, std::string const & element_in ):
		name_( name_in ),
		element_( element_in ),
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

	AtomType(AtomType const & src) :
		name_(src.name_ ),
		element_(src.element_),
		lj_radius_(src.lj_radius_),
		lj_wdepth_(src.lj_wdepth_),
		lk_lambda_(src.lk_lambda_),
		lk_volume_(src.lk_volume_),
		lk_dgfree_(src.lk_dgfree_),
		extra_parameters_(src.extra_parameters_),
		is_acceptor_(src.is_acceptor_),
		is_donor_(src.is_donor_),
		is_polar_hydrogen_(src.is_polar_hydrogen_),
		is_h2o_(src.is_h2o_),
		is_aromatic_(src.is_aromatic_),
		atom_has_orbitals_(src.atom_has_orbitals_),
		atom_is_virtual_(src.atom_is_virtual_),
		atom_is_repulsive_(src.atom_is_repulsive_),
		hybridization_(src.hybridization_)
	{}

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

	std::string const& name() const { return name_; };

	// data
private:

	// name
	std::string name_; // non-const for use in cloning at time of atomtypeset creation...

	// element
	std::string const element_;

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
};

} // chemical
} // core

#endif // INCLUDED_core_chemical_AtomType_HH
