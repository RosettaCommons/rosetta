// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/GaussianFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author James Thompson


// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/GaussianFunc.hh>

#include <core/scoring/constraints/util.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

// AUTO-REMOVED #include <numeric/angle.functions.hh>
// AUTO-REMOVED #include <ObjexxFCL/format.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>

// C++ Headers

#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace constraints {

	void
	GaussianFunc::read_data( std::istream& in ) {
		in 	>> mean_ >> sd_;

		//I'm not sure what in.good() is meant to do here; it does NOT prevent GaussianFunc from chomping the next line of the constraint file as its tag if no tag is present.  SML 08.21.12
		if ( in.good() ) {
			std::string tag;
			in >> tag;

			if (tag == "NOLOG") {
				use_log_score_ = false;
			}
		}

	}

	Real
	GaussianFunc::func( Real const x ) const	{
		if ( use_log_score_ ) return - logdgaussian( x, mean_, sd_, 1 );
		else                  return dgaussian( x, mean_, sd_, 1 );
	} // func

	Real
	GaussianFunc::dfunc( Real const x ) const {
		if ( use_log_score_ ) return - logdgaussian_deriv( x, mean_, sd_, 1 );
		else                  return gaussian_deriv( x, mean_, sd_, 1 );
	} // dfunc

	void GaussianFunc::show_definition( std::ostream& out ) const {
		out << "GAUSSIANFUNC " << mean_ << ' ' << sd_ << "\n";
	}

} // namespace constraints
} // namespace scoring
} // namespace core

