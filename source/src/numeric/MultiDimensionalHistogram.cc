// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/MultiDimensionalHistogram.cc
///
/// @brief  a class for accumulating a histogram of one or more numeric variables
/// @author Colin A. Smith <colin.smith@ucsf.edu>


#include <numeric/MultiDimensionalHistogram.hh>

// C++ headers
#include <iomanip>

namespace numeric {

std::ostream & operator << (
	std::ostream & os,
	MultiDimensionalHistogram const & mdhist
)
{
	std::streamsize oldprecision = os.precision();
	os << std::setprecision(16);

	os << mdhist.num_dimensions() << " " << mdhist.label() << std::endl;

	for ( numeric::Size i = 1; i <= mdhist.num_dimensions(); ++i ) {

		os << mdhist.num_bins()[i] << " " << mdhist.start()[i] << " " << mdhist.end()[i] << " "
			<< mdhist.dim_labels()[i] << std::endl;
	}

	for ( numeric::Size i = 0; i < mdhist.counts().size(); ++i ) {

		if ( i ) os << " ";
		os << mdhist.counts()[i];
	}
	os << std::endl;

	os << std::setprecision(oldprecision);

	return os;
}

} // numeric
