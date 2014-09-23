// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/interpolation/Histogram.fwd.hh
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date   1/26/09


#ifndef INCLUDED_numeric_interpolation_Histogram_fwd_hh
#define INCLUDED_numeric_interpolation_Histogram_fwd_hh

//Utility Headers
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>

namespace numeric {
namespace interpolation {

template< typename, typename>
class Histogram;

/*
 * Templated typedef hack
 * To call, declare like
 *    HistogramOP<X,Y>::Type hist_ptr = new Histogram<X,Y>()
 */
template<class X, class Y>
struct HistogramOP {
	typedef utility::pointer::shared_ptr< Histogram<X,Y> > Type;
};

template<class X, class Y>
struct HistogramCOP {
	typedef utility::pointer::shared_ptr< Histogram<X,Y> const > Type;
};

template<class X, class Y>
struct HistogramAP {
	typedef utility::pointer::weak_ptr< Histogram<X,Y> > Type;
};

template<class X, class Y>
struct HistogramCAP {
	typedef utility::pointer::weak_ptr< Histogram<X,Y> const > Type;
};

} //interpolation
} //numeric


#endif //INCLUDED_numeric_interpolation_sHistogram_FWD_HH
