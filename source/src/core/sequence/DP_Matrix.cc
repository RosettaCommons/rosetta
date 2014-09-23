// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/sequence/DP_Matrix.hh
/// @brief class for holding information on a dynamic programming matrix.
/// @author James Thompson


#include <core/sequence/DP_Matrix.hh>

// C++ headers
// AUTO-REMOVED #include <iostream>
#include <string>

#include <ObjexxFCL/format.hh>


///////////////////////////////////////////////////////////////////////////////
namespace core {
namespace sequence {

using core::Real;
using core::Size;
using std::string;
using utility::vector1;
using ObjexxFCL::format::A;
using ObjexxFCL::format::F;

Cell::~Cell() {}

DP_Matrix::~DP_Matrix() {}

void
DP_Matrix::clear() {
		scoring_matrix_.erase( scoring_matrix_.begin(), scoring_matrix_.end() );
}

CellOP
DP_Matrix::operator ()( Size row, Size col ) const {
	// matrix is implicitly surrounded by 0's
	if ( row < 1 || row > rows() || col < 1 || col > cols() ) {
		return CellOP( new Cell( 0 ) );
	}

	return scoring_matrix_[ row ][ col ];
}

Size DP_Matrix::rows() const {
	return scoring_matrix_.size();
}

Size DP_Matrix::cols() const {
	return scoring_matrix_[1].size();
}

void
DP_Matrix::xlab( vector1< char > const & xs ) {
	xlabels_ = xs;
}

void
DP_Matrix::ylab( vector1< char > const & ys ) {
	ylabels_ = ys;
}

vector1< char >
DP_Matrix::xlab() const {
	return xlabels_;
}

vector1< char >
DP_Matrix::ylab() const {
	return ylabels_;
}

std::ostream & operator<<( std::ostream & out, const DP_Matrix & m ) {
	Size width = 8;
	Size precision = 3;

	vector1< char > xlabs( m.xlab() ), ylabs( m.ylab() );

	for ( Size i = 1; i <= m.rows(); ++i ) {
		if ( ylabs.size() > 1 && i == 1 ) {
			out << A( width, ' ' );

			for ( Size k = 1; k <= m.cols(); ++k )
				out << A( width, ylabs[k] );
			out << std::endl;
		}

		for ( Size j = 1; j <= m.cols(); ++j ) {
			if ( xlabs.size() > 1 && j == 1 ) {
				out << A( width, xlabs[i] );
			}

			out << F( width, precision, m(i,j)->score() );
		}
		out << std::endl;
	}

	return out;
}

} // sequence
} // core
