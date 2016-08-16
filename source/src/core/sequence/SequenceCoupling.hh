// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SequenceCoupling.hh
/// @brief class definition for a sequence coupling profile that
//   represents a  probability distribution over the entire protein.
/// @author Hetu Kamisetty

#ifndef INCLUDED_core_sequence_SequenceCoupling_hh
#define INCLUDED_core_sequence_SequenceCoupling_hh

// Unit headers
#include <core/sequence/SequenceCoupling.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>

// Utility headers
#include <utility/file/FileName.fwd.hh>

// C++ headers
#include <utility/vector1.hh>

namespace core {
namespace sequence {

class SequenceCoupling: public SequenceProfile{
	typedef std::string string;
	typedef utility::file::FileName FileName;
public:

	/// @brief ctors
	SequenceCoupling() : temp_( 1.0 ) {}
	/// @brief copy ctor
	SequenceCoupling( SequenceCoupling const & src ):
		SequenceProfile()
	{
		*this = src;
	}

	/*
	void profile(
	utility::vector1< utility::vector1< Real > > const & new_profile
	);
	utility::vector1< utility::vector1< Real > > const & profile() const;
	*/

	/// @brief assignment operator.
	SequenceCoupling& operator = ( SequenceCoupling const & rhs ) {
		if ( this == &rhs ) return *this;

		id      ( rhs.id() );
		start   ( rhs.start()    );
		gap_char( rhs.gap_char() );
		sequence( rhs.sequence() );

		edgePots( rhs.edgePots() );
		edgeList( rhs.edgeList() );
		profile( rhs.profile() );
		//alphabet( rhs.alphabet() );
		temp    ( rhs.temp() );

		return *this;
	}

	/// @brief dtor
	virtual ~SequenceCoupling() {}

	/// @brief Read a SequenceCoupling model in GREMLIN format .
	virtual void read_from_file( FileName const & fn );

	std::string type() const {
		return "sequence_coupling";
	}

	/// @brief Print this SequenceCoupling object to the given std::ostream.
	/// Commenting out to fix PyRosetta build  friend std::ostream & operator<<( std::ostream & out, const SequenceCoupling & p );

	utility::vector1< utility::vector1< Size> > edgeList() const{
		return edgeList_;
	}

	utility::vector1< utility::vector1< utility::vector1< Real > > > edgePots() const{
		return edgePots_;
	}
	core::Real temp() const{
		return temp_;
	}
	Size npos() const{
		return numVerts_;
	}

	Size findEdgeId(Size vert1, Size vert2) const;

	utility::vector1< utility::vector1 < Real > > const & edgePotBetween(Size vert1, Size vert2) const;
	utility::vector1< utility::vector1 < Real > >  const & edgePotBetween(Size edgeId) const;

	void temp(core::Real const t){
		temp_ = t;
	}
	void edgeList(utility::vector1< utility::vector1< Size > > const & eList){
		edgeList_=eList;
	}
	void edgePots(utility::vector1< utility::vector1< utility::vector1< Real> > > const & ePots){
		edgePots_ = ePots;
	}


	///virtual void generate_from_sequence( Sequence const & seq, std::string matrix="BLOSUM62" );
	//
	///void rescale(core::Real factor=1);

	///void convert_profile_to_probs();

	/// @brief Read profile matrix from the given filename using the NNMAKE
	/// .checkpoint format.
	///void read_from_checkpoint( FileName const & fn );
	//
	virtual void insert_char( core::Size pos, char new_char );

	/// @brief Deletes the given position from the Sequence and shifts
	/// everything else back by one.
	virtual void delete_position( core::Size pos );

	///Size width() const;

	///virtual SequenceOP clone() const {
	/// SequenceOP new_seq_op( new SequenceProfile( *this ) );
	/// return new_seq_op;
	///}


	virtual void read_from_file(utility::file::FileName const & fn, Real temp);

private:
	// temp used to convert arbitrary scores to/from probabilities
	core::Real temp_;
	Size numVerts_, numEdges_;
	//utility::vector1< utility::vector1< Real > > vertexPots_;
	utility::vector1< utility::vector1< Size > > edgeList_;
	utility::vector1< utility::vector1< utility::vector1< Real> > > edgePots_;

}; // class SequenceCoupling

} // sequence
} // core

#endif
