// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableResRotPairFloatMap.hh
/// @brief  A CacheableData that stores a std::map<ResRotPair,float>
/// @author Brian Coventry
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitute.org) -- Added caching of floats indexed by rotamer memory address,
/// used _only_ for self-interactions in the symmetric case.


#ifndef INCLUDED_basic_datacache_CacheableResRotPairFloatMap_hh
#define INCLUDED_basic_datacache_CacheableResRotPairFloatMap_hh

// unit headers
#include <basic/datacache/CacheableResRotPairFloatMap.fwd.hh>

// package headers
#include <basic/datacache/CacheableData.hh>
#include <numeric/MathMatrix.hh>
#include <utility/VirtualBase.fwd.hh>

// C++ headers
#include <map>
#include <unordered_map>
#include <string>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <basic/datacache/CacheableData.fwd.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {

// Using uint32_t because there's no way someone will have more than 4B residues or rotamers
struct ResRotPair {
	uint32_t first_res;
	uint32_t first_rot;
	uint32_t second_res;
	uint32_t second_rot;

	ResRotPair() {}

	ResRotPair( platform::Size first_res_in, platform::Size first_rot_in, platform::Size second_res_in, platform::Size second_rot_in )
	:
		first_res(uint32_t(first_res_in)),
		first_rot(uint32_t(first_rot_in)),
		second_res(uint32_t(second_res_in)),
		second_rot(uint32_t(second_rot_in))
	{}

	bool operator==( ResRotPair const & ot ) const {
		return (
			(first_res == ot.first_res) &&
			(first_rot == ot.first_rot) &&
			(second_res == ot.second_res) &&
			(second_rot == ot.second_rot)
		);
	}

	bool operator<( ResRotPair const & ot ) const {
		if ( first_res != ot.first_res ) return first_res < ot.first_res;
		if ( first_rot != ot.first_rot ) return first_rot < ot.first_rot;
		if ( second_res != ot.second_res ) return second_res < ot.second_res;
		if ( second_rot != ot.second_rot ) return second_rot < ot.second_rot;
		return false;
	}

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION
};

struct ResRotPairHasher {
	std::size_t operator()( ResRotPair const & resrot ) const {
		return size_t(resrot.first_res) << 48 |
			size_t(resrot.first_rot) << 32 |
			size_t(resrot.second_res) << 16 |
			size_t(resrot.second_rot);
	}

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & ) const {}
	template< class Archive > void load( Archive & ) {}
#endif // SERIALIZATION
};


/// @brief Wrapper for std::map< uint64_t, MathMatrix<float> >
class CacheableResRotPairFloatMap : public CacheableData
{
public:
	CacheableResRotPairFloatMap();

	CacheableResRotPairFloatMap( CacheableResRotPairFloatMap const & ot );

	CacheableResRotPairFloatMap & operator=( CacheableResRotPairFloatMap const & ot );

	~CacheableResRotPairFloatMap() override;

	CacheableDataOP
	clone() const override;

	CacheableResRotPairFloatMapOP
	shared_from_this();

	virtual std::unordered_map< ResRotPair, float, ResRotPairHasher > &
	map();

	virtual const std::unordered_map< ResRotPair, float, ResRotPairHasher > &
	map() const;

	/// @brief Nonconst access to map of four ints --> float value.
	/// @details Used only to store energies of a rotamer interacting with its own symmetric copies in the symmetric case.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	virtual std::map< std::tuple< platform::Size, platform::Size, platform::Size, platform::Size>, float > &
	four_int_indexed_map();

	/// @brief Const access to map of four ints --> float value.
	/// @details Used only to store energies of a rotamer interacting with its own symmetric copies in the symmetric case.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	virtual std::map< std::tuple< platform::Size, platform::Size, platform::Size, platform::Size>, float > const &
	four_int_indexed_map() const;

	/// @brief Only perform a shallow copy when clone() is called
	void
	set_shallow_copy( bool shallow );

private:
	/// @brief Map of seqpos1/rotamer1/seqpos2/rotamer2 --> float value.
	utility::pointer::shared_ptr< std::unordered_map< ResRotPair, float, ResRotPairHasher > > mapOP_;

	/// @brief Map of four ints --> float value.
	/// @details Used only to store energies of a rotamer interacting with its own symmetric copies in the symmetric case.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	utility::pointer::shared_ptr< std::map< std::tuple< platform::Size, platform::Size, platform::Size, platform::Size>, float > > four_int_indexed_mapOP_;

	/// @brief Only perform a shallow copy when clone() is called
	bool shallow_copy_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace datacache
} // namespace basic

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( basic_datacache_CacheableResRotPairFloatMap )
#endif // SERIALIZATION


#endif /* INCLUDED_basic_datacache_CacheableStringFloatMap_HH */
