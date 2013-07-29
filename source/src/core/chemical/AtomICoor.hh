// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief declaration of implementation class for abstract class Residue
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_AtomICoor_hh
#define INCLUDED_core_chemical_AtomICoor_hh


// Unit headers
#include <core/chemical/AtomICoor.fwd.hh>

// Project headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
// AUTO-REMOVED #include <core/id/AtomID.hh>

// Utility headers
#include <utility/exit.hh>

#include <core/id/AtomID.fwd.hh>


// Commented by inclean daemon #include <ObjexxFCL/string.functions.hh>

// C++ headers
// Commented by inclean daemon #include <string>


namespace core {
namespace chemical {

/// @brief Atom 's ID in internal coordinates in a ResidueType
class ICoorAtomID {
public:
	typedef conformation::Residue Residue;
	typedef conformation::Conformation Conformation;

public:
	/// ICoordAtomID type
	/**
		 - INTERNAL: atoms which inherently belong to this ResidueType
		 - POLYMER_LOWER: atom at the polymer lower connection, such as backbone C in
		 the previous residue (N-term side)
		 - POLYMER_UPPER: atom at the polymer upper connection, such as backbone N in
		 the next residue (C-term side)
		 - CONNECT: atoms from a non-adjacent residue which connect to this residue
		 by non-polymer connection, such as disulfide
	*/
	enum Type {
		INTERNAL = 1,
		POLYMER_LOWER,
		POLYMER_UPPER,
		CONNECT
	};


public:
	/// @brief default constructor
	ICoorAtomID():
		type_( INTERNAL ),
		atomno_( 0 )
	{}

	/// @brief construct ICoorAtomID by atom name and its ResidueType
	ICoorAtomID(
		std::string name,
		ResidueType const & rsd_type
	);

public:
	/// @brief get ICoorAtomID atomno
	Size
	atomno() const
	{
		return atomno_;
	}

	/// @brief set ICoorAtomID atomno
	void
	atomno( int const atomno_in )
	{
		atomno_ = atomno_in;
	}

	/// @brief get ICoordAtomID type
	Type const &
	type() const
	{
		return type_;
	}

	///
	bool
	is_internal() const
	{
		return ( type_ == INTERNAL );
	}

	///
	bool
	is_polymer_lower() const
	{
		return ( type_ == POLYMER_LOWER );
	}

	///
	bool
	is_polymer_upper() const
	{
		return ( type_ == POLYMER_UPPER );
	}

	///
	bool
	is_connect( Size const connid ) const
	{
		return ( type_ == CONNECT && atomno_ == connid );
	}

public:

	///
	Vector const &
	xyz( Residue const & rsd, Conformation const & conformation ) const;

	///
	Vector // const &
	xyz( ResidueType const & rsd_type ) const;

	/// @brief WARNING: Slightly dangerous function intended for black magic use only.
	///    Only to be used for situations where you *know* the ICoorAtomID can't be anything but
	///    a real atom on the given residue, and where a conformation is absolutely not availible.
	///    If you /can/ use ICoorAtomID::xyz( Residue const &, Conformation const &), you /should/.
	Vector const &
	xyz( conformation::Residue const & rsd ) const;

	///
	id::AtomID
	atom_id( Size const seqpos, Conformation const & conformation ) const;


private:

#ifdef USEBOOSTSERIALIZE
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
			ar & type_;
			ar & atomno_;
	}
	
#endif

	/// atom's "connection" type
	Type type_;
	/// atom's index number
	Size atomno_;
};

/// @brief A basic class containing info of internal coordinates needed for building an atom within a ResidueType
/**
	 In atom tree, each atom is defined by its internal coordinates, which include a bond distance,
	 a bond angle and a torsion angle. Of course, all these internal coordinates are only meaningful
	 in the context of three reference (stub) atoms. AtomICoor information is stored in the residue param
	 files and some terms are defined as following:\n
	 - bond distance d_ is that between the atom to be built (child) and stub_atom1 (parent)
	 - bond angle theta_ is that defined by child-parent-stub2(angle)
	 - torsion angle phi_ is that defined by child-parent-stub2-stub3(torsion)
*/
class AtomICoor {
public:
	/// @brief default constructor
	AtomICoor():
		index_(0),
		phi_(0.0),
		theta_(0.0),
		d_(0.0),
		stub_atom1_(),
		stub_atom2_(),
		stub_atom3_()
	{}

	/// @brief constructor
	AtomICoor(
		Real const phi_in,
		Real const theta_in,
		Real const d_in,
		std::string const & stub_atom1_name,
		std::string const & stub_atom2_name,
		std::string const & stub_atom3_name,
		ResidueType const & rsd_type
	):
		index_(0),
		phi_( phi_in ),
		theta_( theta_in ),
		d_( d_in ),
		stub_atom1_( stub_atom1_name, rsd_type ),
		stub_atom2_( stub_atom2_name, rsd_type ),
		stub_atom3_( stub_atom3_name, rsd_type )
	{}

	AtomICoor(
		Size const index,
		Real const phi_in,
		Real const theta_in,
		Real const d_in,
		std::string const & stub_atom1_name,
		std::string const & stub_atom2_name,
		std::string const & stub_atom3_name,
		ResidueType const & rsd_type
	):
		index_(index),
		phi_( phi_in ),
		theta_( theta_in ),
		d_( d_in ),
		stub_atom1_( stub_atom1_name, rsd_type ),
		stub_atom2_( stub_atom2_name, rsd_type ),
		stub_atom3_( stub_atom3_name, rsd_type )
	{}

public:
	/// @brief accessor to stub_atom1 ICoorAtomID
	Real
	phi() const
	{
		return phi_;
	}

	///
	Real
	theta() const
	{
		return theta_;
	}

	///
	Real
	d() const
	{
		return d_;
	}

	///
	ICoorAtomID const &
	stub_atom1() const
	{
		return stub_atom1_;
	}

	/// @brief accessor to stub_atom2 ICoorAtomID
	ICoorAtomID const &
	stub_atom2() const
	{
		return stub_atom2_;
	}

	/// accessor to stub_atom3 ICoorAtomID
	ICoorAtomID const &
	stub_atom3() const
	{
		return stub_atom3_;
	}

	///
	bool
	is_internal() const
	{
		return ( stub_atom1_.is_internal() && stub_atom2_.is_internal() && stub_atom3_.is_internal() );
	}

	///
	bool
	depends_on_polymer_lower() const
	{
		return ( stub_atom1_.is_polymer_lower() || stub_atom2_.is_polymer_lower() || stub_atom3_.is_polymer_lower() );
	}

	///
	bool
	depends_on_polymer_upper() const
	{
		return ( stub_atom1_.is_polymer_upper() || stub_atom2_.is_polymer_upper() || stub_atom3_.is_polymer_upper() );
	}

	///
	bool
	depends_on_residue_connection( Size const connid ) const
	{
		return ( stub_atom1_.is_connect( connid ) ||
						 stub_atom2_.is_connect( connid ) ||
						 stub_atom3_.is_connect( connid ) );
	}

	/// @brief accessor to stub_atom ICoorAtomID
	ICoorAtomID &
	stub_atom( int const atm )
	{
		switch( atm ) {
		case 1: return stub_atom1_;
		case 2: return stub_atom2_;
		case 3: return stub_atom3_;
		}
		utility_exit_with_message( "ICoorAtomID::stub_atom should be 1--3" );
		return stub_atom1_;
	}

	/// @brief constant accessor to stub_atom ICoorAtomID
	ICoorAtomID const &
	stub_atom( int const atm ) const
	{
		switch( atm ) {
		case 1: return stub_atom1_;
		case 2: return stub_atom2_;
		case 3: return stub_atom3_;
		}
		utility_exit_with_message( "ICoorAtomID::stub_atom should be 1--3" );
		return stub_atom1_;
	}

	void index(core::Size index)
	{
		index_ = index;
	}

	Size index()
	{
		return index_;
	}

public:

	Vector
	build(
		conformation::Residue const & rsd,
		conformation::Conformation const & conformation
	) const;

	Vector
	build(
		ResidueType const & rsd_type
	) const;

	/// @brief WARNING: Slightly dangerous function intended for black magic use only.
	///    Only to be used for situations where you *know* the Atom /and all it's stub atoms/ can't be
	///    anything but real atoms on the given residue, and where a conformation is absolutely not availible.
	///    If you /can/ use AtomICoor::build( Residue const &, Conformation const &), you /should/.
	Vector
	build( conformation::Residue const & rsd ) const;

private:
	
	
#ifdef USEBOOSTSERIALIZE
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
			ar & index_;
			ar & phi_;
			ar & theta_;
			ar & d_;
			ar & stub_atom1_;
			ar & stub_atom2_;
			ar & stub_atom3_;
	}
#endif


	Size index_;
	Real phi_;
	Real theta_;
	Real d_;
	ICoorAtomID stub_atom1_;
	ICoorAtomID stub_atom2_;
	ICoorAtomID stub_atom3_;
};

} // chemical
} // core



#endif
