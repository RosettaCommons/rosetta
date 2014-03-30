// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/Interface_/DDPlookup.hh
/// @author Hermann Zellner (hermann1.zellner@biologie.uni-regensburg.de)

#ifndef INCLUDED_core_scoring_interface_DDPlookup_hh
#define INCLUDED_core_scoring_interface_DDPlookup_hh

// Unit headers
#include <core/scoring/interface_/DDPlookup.fwd.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/types.hh>

// Numeric headers
#include <numeric/interpolation/spline/SplineGenerator.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <string>

namespace core {
namespace scoring {
namespace interface_ {

class DDPlookup : public utility::pointer::ReferenceCount {
public:
	DDPlookup(std::string filename);

	core::Real
	get_potentials(
			const core::chemical::AA & aa1, const core::chemical::AA & aa2, core::Real distance
			) const;

private:
	numeric::interpolation::spline::SplineGenerator*** lookup_table_;
	utility::vector1< utility::vector1< core::Real > > left_;
	utility::vector1< utility::vector1< core::Real > > right_;
};


} //Interface_
} //scoring
} //core

#endif /* INCLUDED_core_scoring_Interface_DDPLOOKUP_HH_ */
