// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/util.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_loop_graph_evaluator_util_HH
#define INCLUDED_core_scoring_loop_graph_evaluator_util_HH

#include <core/types.hh>
#include <utility/json_spirit/json_spirit.h>
#include <utility/json_spirit/json_spirit_tools.hh>
#include <numeric/MathNTensor.hh>

namespace core {
namespace scoring {
namespace loop_graph {
namespace evaluator {

	void
	inline
	get_minval_binwidth( numeric::MathNTensor< double, 6 > const & T,
											 utility::json_spirit::mObject const & json,
											 utility::fixedsizearray1< numeric::Real, 6 > & minval,
											 utility::fixedsizearray1< numeric::Real, 6 > & binwidth)
	{
		using namespace utility::json_spirit;
		mArray json_binwidth( get_mArray( json, "binwidth" ) );
		mArray json_minval( get_mArray( json, "minval" ) );
		mArray json_maxval( get_mArray( json, "maxval" ) );
		for ( numeric::Size i=1; i<=6; ++i ) {
			runtime_assert( numeric::Size( ( json_maxval[i-1].get_real() - json_minval[i-1].get_real() )/json_binwidth[i-1].get_real() + 1 ) == T.n_bins( i ) );
			minval[ i ]   = json_minval[i-1].get_real();
			binwidth[ i ] = json_binwidth[i-1].get_real();
		}
	}

} //evaluator
} //loop_graph
} //scoring
} //core

#endif
