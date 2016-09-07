// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/sequence/DP_Matrix.hh
/// @brief class for holding information on a dynamic programming matrix.
/// @author James Thompson

#ifndef INCLUDED_core_sequence_DP_Matrix_hh
#define INCLUDED_core_sequence_DP_Matrix_hh

#include <core/types.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1_bool.hh>

///////////////////////////////////////////////////////////////////////////////

namespace core {
namespace sequence {

enum AlignMove {
	diagonal = 1,
	left,
	above,
	end
};

// simple little class to hold scores and store pointers to other Cells. forward declaration
// for holding onto CellOPs
class Cell;
typedef utility::pointer::shared_ptr< Cell > CellOP;
class Cell : public utility::pointer::ReferenceCount {

public:
	Cell() :
		score_( 0.0 ), type_( diagonal ), coords_( 2, 0 )
	{}

	Cell(
		Real const & sc,
		AlignMove const & ty = end
	) :
		score_( 0.0 ), type_( ty ), coords_( 2, 0 )
	{
		score( sc );
		came_from( ty );
		backptr_.reset();
	}

	~Cell() override;

	void score( const Real & score ) {
		score_ = score;
	}

	Real score() const {
		return score_;
	}

	Cell & operator =( Cell const & c ) {
		score( c.score() );
		return *this ;
	}

	void next( CellOP n ) {
		backptr_ = n;
	}

	CellOP next() {
		return backptr_;
	}

	AlignMove came_from() {
		return type_;
	}

	void came_from( AlignMove const & type ) {
		type_ = type;
	}

	void x( Size const & s ) {
		coords_[1] = s;
	}

	void y( Size const & s ) {
		coords_[2] = s;
	}

	Size x() {
		return coords_[1];
	}

	Size y() {
		return coords_[2];
	}

private:
	Real score_;
	CellOP backptr_;
	AlignMove type_;
	utility::vector1< Size > coords_;
}; // class Cell

class DP_Matrix : public utility::pointer::ReferenceCount {

	typedef utility::pointer::shared_ptr< Cell > CellOP;
	typedef utility::vector1< CellOP > Row;
	typedef utility::vector1< utility::vector1< CellOP > > Matrix;

public:
	DP_Matrix( Size rs, Size cs ) {
		for ( Size i = 1; i <= rs; ++i ) {
			Row temp;
			for ( Size j = 1; j <= cs; ++j ) {
				CellOP c( new Cell( 0.0, end ) );
				c->x( i );
				c->y( j );
				temp.push_back( c );
			}
			scoring_matrix_.push_back( temp );
		}
	}
	~DP_Matrix() override;

	void clear();

	CellOP
	operator ()( Size row, Size col ) const;

	Size rows() const;

	Size cols() const;

	void xlab( utility::vector1< char > const & xs );

	void ylab( utility::vector1< char > const & ys );

	utility::vector1< char > xlab() const;

	utility::vector1< char > ylab() const;

	friend std::ostream & operator<<( std::ostream & out, const DP_Matrix & m );

private:
	Matrix scoring_matrix_;
	utility::vector1< char > xlabels_, ylabels_;
}; // DP_Matrix


} // sequence
} // core

#endif
