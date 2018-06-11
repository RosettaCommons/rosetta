// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SequenceProfile.hh
/// @brief class definition for a sequence profile that represents each
/// position in a sequence as a probability distribution over the allowed amino
/// acids at that position.
/// @author James Thompson

#ifndef INCLUDED_core_sequence_SequenceProfile_hh
#define INCLUDED_core_sequence_SequenceProfile_hh

// Unit headers
#include <core/sequence/SequenceProfile.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/sequence/Sequence.hh>
#include <core/chemical/AA.hh>


// Utility headers
#include <utility/file/FileName.fwd.hh>

// C++ headers
#include <utility/vector1.hh>
#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace sequence {

class SequenceProfile : public Sequence {
	typedef std::string string;
	typedef utility::file::FileName FileName;
public:

	/// @brief ctors
	SequenceProfile() :
		temp_( 1.0 ),
		negative_better_(false)
	{}

	SequenceProfile( FileName const & fn ) :
		temp_( 1.0 ),
		negative_better_(false)
	{
		read_from_file( fn );
	}

	SequenceProfile(
		utility::vector1< utility::vector1< core::Real > > prof,
		std::string const & sequence,
		std::string const & id,
		Size start = 1,
		bool negative_better = false
	) :
		Sequence( sequence, id, start ),
		temp_( 1.0 ),
		negative_better_(negative_better)
	{
		profile( prof );
		debug_assert( profile().size() == length() );
	}

	/// @brief copy ctor
	SequenceProfile( SequenceProfile const & src ):
		Sequence()
	{
		*this = src;
	}

	/// @brief assignment operator.
	SequenceProfile & operator = ( SequenceProfile const & rhs ) {
		if ( this == &rhs ) return *this;

		id      ( rhs.id() );
		start   ( rhs.start()    );
		gap_char( rhs.gap_char() );
		sequence( rhs.sequence() );

		profile ( rhs.profile() );
		alphabet( rhs.alphabet() );
		temp_ = rhs.temp_;
		negative_better_ = rhs.negative_better_;

		return *this;
	}

	/// @brief dtor
	~SequenceProfile() override = default;

	/// @brief Returns an owning pointer to a new SequenceProfile object,
	/// with data that is a deep copy of the information in this object.
	SequenceOP clone() const override {
		SequenceOP new_seq_op( new SequenceProfile( *this ) );
		return new_seq_op;
	}

	/// @brief Read an profile matrix from the given filename using the NCBI
	/// PSSM format for a position-specific scoring matrix.
	void read_from_file( FileName const & fn ) override;

	/// @brief Generate the profile matrix from a sequence and a given substitution matrix
	virtual void generate_from_sequence( Sequence const & seq, std::string matrix="BLOSUM62" );

	/// @brief Multiply all profile weights by factor
	void rescale(core::Real factor=1);

	/// @brief Use boltzman scaling on a per-residue basis to convert the current profile values to probabilities ( range 0.0-1.0 )
	void convert_profile_to_probs( core::Real temp = 1.0 );

	/// @brief Use linear rescaling (with a fixed zero) to fit the values within the range -1.0 to 1.0
	void global_auto_rescale();

	/// @brief Read profile matrix from the given filename using the NNMAKE
	/// .checkpoint format.
	/// For compatability, negative_better defaults to true. Set manually if necessary.
	void read_from_checkpoint( FileName const & fn, bool negative_better = true );

	/// @brief Read profile matrix from the given filename in the legacy BLAST binary format
	void read_from_binary_chk(FileName const & fn);

	/// @brief Returns the 2D vector1 of Real values representing this profile.
	utility::vector1< utility::vector1< Real > > const & profile() const;

	/// @brief Returns the 2D vector1 of Real values of the probabilties of each aa.
	utility::vector1< utility::vector1< Real > > const & occurrence_data() const;

	/// @brief Sets the 2D vector1 of Real values representing this profile.
	void profile(
		utility::vector1< utility::vector1< Real > > const & new_profile
	);
	/// @brief Sets the 2D vector1 of Real values of the probabilties of each aa.
	void occurrence_data(
		utility::vector1< utility::vector1< Real > >  const & new_occurrence_data
	);

	/// @brief Inserts a character at the given position.
	void insert_char( core::Size pos, char new_char ) override;

	/// @brief Deletes the given position from the Sequence and shifts
	/// everything else back by one.
	void delete_position( core::Size pos ) override;

	std::string type() const override {
		return "sequence_profile";
	}

	/// @brief Returns the number of distinct values at each position in this
	/// profile.
	Size width() const;

	/// @brief Returns the vector1 of values at this position.
	//utility::vector1< Real >
	utility::vector1< Real > const &
	prof_row( Size pos ) const;

	utility::vector1< Real > const &
	probability_row( Size pos ) const;

	/// @brief Sets the 1D vector1 of Real values representing this profile at pos X.
	void prof_row(
		utility::vector1< Real > const & new_prof_row, core::Size pos
	);

	void probabilty_row(
		utility::vector1< Real > const & new_prob_row, core::Size pos
	);

	Size size() const {
		return profile_.size();
	}

	/// @brief Returns true if negative values are better identities.
	/// @details The "default" use case is for storing log likelihood values
	/// where positive is better. If you're using this class to store energy-like
	/// values, set negative_better to true.
	bool negative_better() const { return negative_better_; }

	/// @brief Set whether negative identities are better.
	void negative_better( bool negbet ) { negative_better_ = negbet; }

	/// @brief returns the temperature used in computing profile probabilities
	core::Real temp() const {
		return temp_;
	}

	/// @brief Return the alphabet used by this sequence profile. This is an
	/// N-dimensional vector1 where N is the width of the profile, and the ith
	/// entry of any row in the profile represents the probability of ith
	/// character at that row in the sequence.
	utility::vector1< std::string > alphabet() const {
		return alphabet_;
	}

	void alphabet( utility::vector1< std::string > new_alphabet ) {
		alphabet_ = new_alphabet;
	}

	/// @brief Print this SequenceProfile object to the given std::ostream.
	friend std::ostream & operator<<(
		std::ostream & out, const SequenceProfile & p
	);

	bool operator==( SequenceProfile const & other ) const;

private:

	/// @brief converts a vector1 of arbitrary scores to values using Bolztmann
	/// averaging at a given kT. Scores should follow the convention that more positive -> better score.
	/// If not, set negative_better to true.
	void scores_to_probs_(
		utility::vector1< core::Real > & scores,
		core::Real kT,
		bool negative_better = false
	) const;

	/// @brief Internal consistency check. Returns true if passed, causes a runtime_assertion failure if not.
	bool check_internals_() const;

	utility::vector1< std::string > alphabet_;
	utility::vector1< utility::vector1< Real > > profile_;
	utility::vector1< utility::vector1< Real > > occurrence_data_;//This matrix holds the the % of occurrences of an amino acid in the data base used to construct the pssm matrix.
	// This data is genrated by psi-blast program as an additional matrix along side the pssm matrix


	/// @brief temp used to convert arbitrary scores to/from probabilities
	core::Real temp_;
	/// @brief The orientation of the values. Are negative values better than zero/positive ones?
	bool negative_better_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class SequenceProfile

} // sequence
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_sequence_SequenceProfile )
#endif // SERIALIZATION


#endif
