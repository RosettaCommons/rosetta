// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/Histogram.hh
/// @brief Light-weight histogram class used to save data for simulated tempering.
/// @detailed intentionally uses floats to save memory!
///           TODO: replace with one of the numeric::histograms classes?
/// @author Fang-Chieh Chou
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_recces_Histogram_HH
#define INCLUDED_protocols_recces_Histogram_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/recces/Histogram.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace recces {


class Histogram: public utility::pointer::ReferenceCount {

public:

	//constructor
	Histogram( core::Real const min,
		core::Real const max,
		core::Real const spacing );

	//destructor
	~Histogram();

public:

	void add( float const value, core::Size const n_items );

	void clear();

	utility::vector1<core::Real> get_scores() const;

	utility::vector1<core::Size> const & get_hist() const { return hist_; }

private:

	core::Real min_, max_, spacing_;
	core::Size n_elem_;
	utility::vector1<core::Size> hist_;

};

} //recces
} //protocols

#endif
