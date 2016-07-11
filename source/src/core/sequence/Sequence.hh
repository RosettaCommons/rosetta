// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Sequence.hh
/// @brief class definition for a sequence
/// @author James Thompson

#ifndef INCLUDED_core_sequence_Sequence_hh
#define INCLUDED_core_sequence_Sequence_hh

// Unit headers
#include <core/sequence/Sequence.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/file/FileName.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers

#include <string>

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace sequence {

class Sequence : public utility::pointer::ReferenceCount {

public:

	/// @brief ctor
	Sequence() : id_( "blank" ), start_( 1 ), gap_char_('-') {}
	Sequence( std::string seq, std::string id, core::Size start = 1 ) :
		id_( id ),
		start_( start ),
		gap_char_('-')
	{
		sequence( seq );
	}

	Sequence( core::pose::Pose const & pose );

	/// @brief copy constructor.
	Sequence( Sequence const & src ) :
		ReferenceCount()
	{
		*this = src;
	}

public: // virtual functions
	/// @brief dtor
	virtual ~Sequence();

	/// @brief Returns an owning pointer to a copy of this sequence.
	virtual SequenceOP clone() const;

	/// @brief initializes this sequence object from a file.
	virtual void read_from_file( utility::file::FileName const & /*fn*/ );

	/// @brief Returns the number of characters in this object.
	virtual core::Size length() const;
	/// @brief Inserts a character at the given position.

	virtual void insert_char( core::Size pos, char new_char );

	/// @brief Deletes the given position from the Sequence and shifts
	/// everything else back by one.
	virtual void delete_position( core::Size pos );

	virtual std::string to_string() const;

	virtual std::string type() const;

public: // non-virtual functions

	/// @brief sets sequence to the given value.
	void sequence( std::string sequence );

	/// @brief sets id to the given value.
	void id      ( std::string new_id );

	/// @brief sets starting index to the given value.
	void start   ( core::Size new_start );

	/// @brief sets gap_char to the given value.
	void gap_char( char gap_char );

	/// @brief sets spacer_positions to the given value.
	void spacer_positions( utility::vector1< Size > const & spacer_positions );

	/// @brief Returns the number of characters in this object, ignoring gaps.
	core::Size ungapped_length() const;

	/// @brief Returns the start of this object.
	core::Size start() const;

	/// @brief Returns the id of this object.
	std::string id() const;

	/// @brief Returns the character used to represent a gap for this object.
	char gap_char() const;

	/// @brief Returns the spacer positions (chain breaks, e.g., from spaces or commas in original sequence text)
	utility::vector1< Size > spacer_positions() const;

	/// @brief Returns the string representing this sequence without gaps.
	std::string ungapped_sequence() const;

	/// @brief Returns the full sequence, which may include gaps.
	std::string sequence() const;

	/// @brief assignment operator.
	Sequence & operator = ( Sequence const & rhs ) {
		if ( this == &rhs ) return *this;

		start_    = rhs.start();
		id_       = rhs.id();
		gap_char_ = rhs.gap_char();
		seq_      = rhs.sequence();
		spacer_positions_ = rhs.spacer_positions();

		return *this;
	}

	/// @brief Returns true if this Sequence object's id is lexicographically
	/// less than the given Sequence object's id, returns false otherwise. Uses
	/// C++ string < operator to compare this Sequence object's id() to the
	/// given Sequence object's id().
	inline bool operator < ( const Sequence & s ) const {
		return id() < s.id();
	}

	/// @brief Returns true if the given Sequence object is equal to this
	/// Sequence object. Tests for string equality of id(), start(), and
	/// sequence(), and returns false if any of these are not equal.
	inline bool operator == ( const Sequence & s ) const {
		return (
			id()       == s.id()    &&
			start()    == s.start() &&
			sequence() == s.sequence() &&
			spacer_positions() == s.spacer_positions()
		);
	}

	/// @brief Returns the character at the given sequence position.
	char operator[]( core::Size pos ) const;
	char at( core::Size pos ) const;


	/// @brief Inserts a gap at the given position, where insert_gap( 0 )
	/// inserts the character at the beginning of the sequence, and
	/// insert_gap( length() ) inserts the character at the end of the
	/// sequence.
	void insert_gap( core::Size pos );

	/// @brief Append a character
	void append_char( char new_char );

	/// @brief Append a gap
	void append_gap( );

	/// @brief Returns true if this position in the sequence represents a gap,
	/// returns false otherwise.
	bool is_gap( core::Size pos ) const;

	/// @brief Initializes the information in this sequence from the given
	/// std::istream.  The istream should yield three pieces of information in
	/// the following order:
	/// - id
	/// - start
	/// - sequence

	void read_data( std::istream & in );

	/// @brief Returns the index of the given sequence position, which is the
	/// position in the sequence minus any gaps that occur earlier in the
	/// sequence. For example, if the sequence is ---AT, resnum(5) will return 2.
	/// Returns 0 for unaligned positions.
	core::Size resnum( core::Size idx ) const;

	/// @brief Prints the information a given Sequence object to the given
	/// std::ostream.
	friend std::ostream & operator<<( std::ostream & out, const Sequence & seq );

	/// @brief Prints the information a given Sequence object to the given
	/// std::ostream.
	friend std::ostream & operator<<( std::istream & out, Sequence & seq );
	friend std::istream & operator>>( std::istream & in,  Sequence & seq );

private:
	std::string id_;
	core::Size start_;
	char gap_char_;
	utility::vector1< Size > spacer_positions_;

	std::string seq_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class Sequence

} // sequence
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_sequence_Sequence )
#endif // SERIALIZATION


#endif
