// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/Interface/DDPlookup.hh
/// @author Hermann Zellner (hermann1.zellner@biologie.uni-regensburg.de)

#ifndef INCLUDED_core_scoring_interface_DDPlookup_hh
#define INCLUDED_core_scoring_interface_DDPlookup_hh

#include <core/chemical/AA.hh>

//#include <ObjexxFCL/FArray3D.hh>

//#include <utility/io/izstream.hh>
//#include <numeric/interpolation/spline/Interpolator.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>

#include <core/chemical/AA.hh>
#include <core/types.hh>



namespace core {
namespace scoring {
namespace interface {

class DDPlookup {
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


} //Interface
} //scoring
} //core

#endif /* INCLUDED_core_scoring_Interface_DDPLOOKUP_HH_ */

