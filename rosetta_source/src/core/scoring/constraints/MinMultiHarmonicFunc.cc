// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/HarmonicFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


#include <core/scoring/constraints/MinMultiHarmonicFunc.hh>

#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

// AUTO-REMOVED #include <numeric/angle.functions.hh>
// AUTO-REMOVED #include <numeric/random.functions.hh>
// AUTO-REMOVED #include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>




// C++ Headers

namespace core {
namespace scoring {
namespace constraints {

static basic::Tracer trMinMultiHarmonicFunc(
                "fragment.picking.scores.MinMultiHarmonicFunc");

MinMultiHarmonicFunc::MinMultiHarmonicFunc( utility::vector1<Real> const & x0_in, utility::vector1<Real> const & sd_in ) {

    assert(x0_in.size() == sd_in.size());
    x0_.clear();
    sd_.clear();
    for(Size i=1;i<=x0_in.size();++i) {
      x0_.push_back( x0_in[i] );
      sd_.push_back( sd_in[i] );
    }
}


Real
MinMultiHarmonicFunc::func( Real const x )  const {

        Real min = ( x-x0_[1] );
	min *= min;
	Real z2 = min;
	which_component_ = 1;
	for(Size i=2;i<=n_;++i) {
	  Real const z = ( x-x0_[i] );
	  if( (z2=z*z) < min ) {
	    which_component_ = i;
	    min = z2;
	  }
	}

	return min / sd_[which_component_];
}

Real
MinMultiHarmonicFunc::dfunc( Real const x ) const {
        Real min = ( x-x0_[1] )/sd_[1];
	which_component_ = 1;
	for(Size i=2;i<=n_;++i) {
	  Real const z = ( x-x0_[i] )/sd_[i];
	  if( z < min ) {
	    which_component_ = i;
	    min = z;
	  }
	}

	return 2 * min;
}

void
MinMultiHarmonicFunc::read_data( std::istream& in ) {

  x0_.clear();
  sd_.clear();
  Real r;
  std::string line;
  getline( in, line );
  std::istringstream line_stream( line );

  utility::vector1<Real> entries;
  do {
    line_stream >> r;
    entries.push_back( r );
  }
  while( !line_stream.fail() );
  n_ = entries.size() /2;
  if(n_*2 != entries.size() ) {
    trMinMultiHarmonicFunc.Warning<< "Expected even number of parameters but got"<<entries.size()<<
        "; the last value seems to be useless or one is missing"<<std::endl;
  }
  for(Size i=1;i<=n_;++i) {
    x0_.push_back( entries[i*2-1] );
    sd_.push_back( entries[i*2] );
  }
}

void
MinMultiHarmonicFunc::show_definition( std::ostream &out ) const {

	out << "MINMULTIHARMONIC";
	for(Size i=1;i<=n_;++i)
	   out<< " "<<x0_[i] << " " << sd_[i];
	out << std::endl;
}

Size
MinMultiHarmonicFunc::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {
	if (verbose_level > 100 ) {
		out << "HARM " <<  this->func(x) << std::endl;
	} else if (verbose_level > 70 ) {
		if ( ( x < x0_[which_component_])  && ( this->func(x) > threshold ) ) out << "-";
		else if ( ( x > x0_[which_component_]) && ( this->func(x) > threshold )) out << "+";
		else out << ".";
	}
	return Func::show_violations( out, x, verbose_level, threshold);
}

} // namespace constraints
} // namespace scoring
} // namespace core

