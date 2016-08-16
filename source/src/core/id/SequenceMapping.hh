// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_id_SequenceMapping_hh
#define INCLUDED_core_id_SequenceMapping_hh

// Unit headers
#include <core/id/SequenceMapping.fwd.hh>

// Project headers
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace id {

class SequenceMapping : public utility::pointer::ReferenceCount {

public:
	// constructors, destructors and assignment operator
	/// @brief ctor
	SequenceMapping() :
		size2_(0)
	{}

	/// @brief ctor
	SequenceMapping( Size const s1, Size const s2 ) :
		size2_(s2),
		mapping_(s1,0)
	{}

	/// @brief ctor
	SequenceMapping( utility::vector1< Size > const & mapping );

	/// @brief convenience constructor from LengthEvent
	SequenceMapping( conformation::signals::LengthEvent const & event );

	/// @brief dtor
	virtual ~SequenceMapping();

	/// @brief copy constructor
	SequenceMapping( SequenceMapping const & src );

	SequenceMapping &
	operator = ( SequenceMapping const & src );

public:

	/// @brief resize
	void resize( Size const s1, Size const s2 );

	/// @brief go from an A->B mapping to a B->A mapping
	void reverse();

	/// @brief Apply a B->C mapping to the current A->B mapping to get an A->C mapping
	/// i.e. smap[j] becomes smap_to_add[ smap[j] ]
	void downstream_combine( core::id::SequenceMapping const & smap_to_add );

	/// @brief Apply a C->A mapping to the current A->B mapping to get a C->B mapping
	/// i.e. smap[j] becomes smap[ smap_to_add[ j ] ]
	void upstream_combine( core::id::SequenceMapping const & smap_to_add );

	/// @brief size of target sequence
	Size size1() const;

	/// @brief size of aligned sequence ???
	Size size2() const;

	//@brief the pose might have a missing N-terminal offset the sequence mapping accordingly
	//NOTE: set_offset( 5) followed by set_offset( -5 ) might not be an identity operation:
	//if residues are pushed below zero... there is no coming back.
	void set_offset( int setting );

	bool all_aligned() const;

	bool is_identity() const;
	bool is_identity_ignore_gaps() const;

	void size2( Size const s2 );

	void push_back( Size const al );

	void delete_source_residue( Size const pos1 );

	void show() const;

	void show( std::ostream & output ) const;

	void insert_source_residue( Size const pos1 );

	void insert_aligned_residue( Size const pos1, Size const pos2 );

	/// @brief same as insert_aligned_residue, but a couple of extra checks on size1 and size2.
	void insert_aligned_residue_safe( Size const pos1, Size const pos2 );

	void insert_target_residue( Size const pos );

	void delete_target_residue( Size const pos );

	void clear();

	// residue of aligned sequence at target position pos1
	Size   operator[]( Size const pos1 ) const;
	Size & operator[]( Size const pos1 );

	/// @brief Equality operator.
	virtual bool operator == ( SequenceMapping const & rhs ) const;

	/// @brief virtual function for ensuring both sequence mappings are the same type;
	/// essential for a valid equality operator.
	virtual bool same_type_as_me( SequenceMapping const & other ) const;


	/// @brief Construct an identity mapping, which means that for all positions,
	/// i in seq1 maps to i in seq2.
	inline
	static
	SequenceMapping
	identity( Size const size )
	{
		SequenceMapping id( size, size );
		for ( Size i=1; i<= size; ++i ) {
			id[i] = i;
		}
		runtime_assert( id.is_identity() && id.size1() == size );
		return id;
	}

	utility::vector1< core::Size > const & mapping() const {
		return mapping_;
	}

	void mapping( utility::vector1< core::Size > const & mapping ) {
		mapping_ = mapping;
	}

	std::string to_string() const;

private:
	Size size2_; // length of second sequence
	utility::vector1< Size > mapping_; // mapping_[ residue_i ] = residue_j
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class SequenceMapping

inline std::ostream & operator<< (
	std::ostream & out,
	SequenceMapping const & map
) {
	map.show( out );
	return out;
}

/// @brief make one sequence mapping out of all input ones
/// utility function added by flo, feb 2011
core::id::SequenceMappingOP
combine_sequence_mappings(
	utility::vector1< core::id::SequenceMapping > const & smaps
);

/// @brief combine the input sequence mappings into one
/// utility function added by flo, feb 2011
void
combine_sequence_mappings(
	core::id::SequenceMapping & smap,
	core::id::SequenceMapping const & smap_to_add
);

} // id
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_id_SequenceMapping )
#endif // SERIALIZATION


#endif
