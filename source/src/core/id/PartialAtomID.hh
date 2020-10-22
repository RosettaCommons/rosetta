// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/id/PartialAtomID.hh
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_id_PartialAtomID_hh
#define INCLUDED_core_id_PartialAtomID_hh

// Unit headers
#include <core/id/PartialAtomID.fwd.hh>

#include <utility/exit.hh>
#include <utility/json_spirit/json_spirit_value.h>
#include <utility/tools/make_vector.hh>

// C++ headers

#include <core/types.hh>


namespace core {
namespace id {


/////////////////////////////////////////////////////
// PartialAtomID
// Specify either a particular atom on a residue or
// an atom near an inter-residue connection atom on
// a residue
/////////////////////////////////////////////////////


/// @brief Partial atom identifier class. Defined by an atom index
/// and a residue number, or by an inter-residue-connection index,
/// the number of bonds away from the inter-residue-connection atom,
/// and a residue number.
class PartialAtomID
{

public: // Creation

	/// @brief Default constructor
	inline
	constexpr
	PartialAtomID() :
		atomno_( 0 ),
		resconnid_( 0 ),
		bonds_from_resconn_( 0 ),
		rsd_( 0 )
	{}

	/// @brief Copy constructor
	inline
	PartialAtomID( PartialAtomID const & ) = default;

	/// @brief AtomID-like constructor; this ID is fully (not partially) resolved
	inline
	PartialAtomID(
		Size const atomno_in,
		Size const rsd_in
	) :
		atomno_( atomno_in ),
		resconnid_( 0 ),
		bonds_from_resconn_( 0 ),
		rsd_( rsd_in )
	{}

	/// @brief Constructor to specify an inter-residue connection atom
	///
	/// @details The third argument specifies the number of bonds
	/// away from the residue-connection atom that the indicated
	/// atom is. The three-argument constructor also disambiguates
	/// this partial ID from the two-argument constructor that is
	/// essentially a fully-resolved AtomID
	inline
	PartialAtomID(
		Size const resconnid_in,
		Size const rsd_in,
		Size const bonds_from_resconn_in
	) :
		atomno_( 0 ),
		resconnid_( resconnid_in ),
		bonds_from_resconn_( bonds_from_resconn_in ),
		rsd_( rsd_in )
	{}

	static constexpr PartialAtomID BOGUS_PARTIAL_ATOM_ID() { return PartialAtomID(); }
	static constexpr PartialAtomID CHAINBREAK_BOGUS_PARTIAL_ATOM_ID() { return PartialAtomID(); }

	///@brief Set the value of atom and residue.
	void
	set_complete( core::Size atomno_in, core::Size rsd_in);

	///@brief Set the value of the inter-residue connection index and residue.
	void
	set_partial( core::Size resconnid_in, core::Size bonds_from_resconn_in, core::Size rsd_in);

public: // Properties

	/// @brief Returns the residue number
	inline
	Size
	rsd() const { return rsd_; }

	/// @brief Returns the atom index
	inline
	Size
	atomno() const { return atomno_; }

	/// @brief Returns the inter-residue connection index
	inline
	Size
	resconnid() const { return resconnid_; }

	/// @brief Returns the number of chemical bonds away from the
	/// residue-connection atom the target atom is located
	inline
	Size
	bonds_from_resconn() const { return bonds_from_resconn_; }

	/// @brief Returns true if the PartialAtomID is valid
	/// @note must return false for BOGUS_ATOM_ID
	inline
	bool
	valid() const { return rsd_ > 0 && ( atomno_ > 0 || resconnid_ > 0 ); }

	inline
	bool
	complete() const { return rsd_ > 0 && atomno_ > 0; }

	inline
	bool
	partial() const { return rsd_ > 0 && resconnid_ > 0; }


	/// @brief serialize an AtomID to a json_spirit object
	inline
	utility::json_spirit::Value serialize() const {
		utility::json_spirit::Pair atomno("atomno", utility::json_spirit::Value(static_cast<uint64_t>(atomno_)));
		utility::json_spirit::Pair resconnid("resconnid", utility::json_spirit::Value(static_cast<uint64_t>(resconnid_)));
		utility::json_spirit::Pair bonds_from_resconn("bonds_from_resconn", utility::json_spirit::Value(static_cast<uint64_t>(bonds_from_resconn_)));
		utility::json_spirit::Pair rsd("rsd", utility::json_spirit::Value(static_cast<uint64_t>(rsd_)));
		return utility::json_spirit::Value(utility::tools::make_vector(
			atomno,resconnid, bonds_from_resconn,rsd));
	}

	/// @brief deserialize a json_spirit object to an AtomID
	inline
	void deserialize(utility::json_spirit::mObject data) {
		atomno_ = static_cast<Size>(data["atomno"].get_uint64());
		resconnid_ = static_cast<Size>(data["resconnid"].get_uint64());
		bonds_from_resconn_ = static_cast<Size>(data["bonds_from_resconn"].get_uint64());
		rsd_    = static_cast<Size>(data["rsd"].get_uint64());
	}

public: // Friends

	friend
	std::ostream &
	operator <<(
		std::ostream & os,
		PartialAtomID const & a
	);

	/// @brief a and b are the same atom
	friend
	inline
	bool
	operator ==(
		PartialAtomID const & a,
		PartialAtomID const & b
	) {
		return a.atomno_ == b.atomno_ &&
			a.resconnid_ == b.resconnid_ &&
			a.bonds_from_resconn_ == b.bonds_from_resconn_ &&
			a.rsd_ == b.rsd_;
	}

	/// @brief a and b are different atom
	friend
	inline
	bool
	operator !=(
		PartialAtomID const & a,
		PartialAtomID const & b
	) { return a.atomno_ != b.atomno_ || a.resconnid_ != b.resconnid_ || a.bonds_from_resconn_ != b.bonds_from_resconn_ || a.rsd_ != b.rsd_; }

	/// @brief a is LOWER than b (e.g., first by smaller residue index number then by smaller atom index number)
	friend
	inline
	bool
	operator <(
		PartialAtomID const & a,
		PartialAtomID const & b
	)
	{
		return std::make_tuple(a.rsd_, a.atomno_, a.resconnid_, a.bonds_from_resconn_) <
			std::make_tuple(b.rsd_, b.atomno_, b.resconnid_, b.bonds_from_resconn_);
	}

private: // Fields


	/// @brief Atom number within the Residue
	Size atomno_;

	/// @brief The index of the inter-residue connection
	/// for the atom of interest; the index of that atom
	/// may vary between different residue types and so
	/// when this information is used, it only partially
	/// specifies the identity of the atom on a residue
	Size resconnid_;

	/// @brief The number of chemical bonds from the inter-
	/// residue connection atom and the atom of interest.
	/// 0 if the atom of interest is the inter-residue-
	/// connection atom.
	Size bonds_from_resconn_;

	/// @brief Residue number within the complex
	Size rsd_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // PartialAtomID

/// Globals -- may not be used until after core::init is called.
extern PartialAtomID const GLOBAL_BOGUS_PARTIAL_ATOM_ID;
extern PartialAtomID const GLOBAL_CHAINBREAK_BOGUS_PARTIAL_ATOM_ID;


} // namespace id
} // namespace core


#endif // INCLUDED_core_id_PartialAtomID_HH
