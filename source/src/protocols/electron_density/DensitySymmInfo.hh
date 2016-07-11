// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief This represents the symmetry info needed for dealing with symmetric density maps.
///        Symmetry is detected from the map directly, without a model.
/// @author fpd

#ifndef INCLUDED_protocols_electron_density_DensitySymmInfo_hh
#define INCLUDED_protocols_electron_density_DensitySymmInfo_hh

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/constants.hh>
#include <numeric/xyzVector.io.hh>
#include <core/types.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>

namespace protocols {
namespace electron_density {


/// scale density map intensities to match a pose's
// helper class: find & store symmetry axes
class DensitySymmInfo {
public:
	DensitySymmInfo() :
		type('C'), count_primary(1), symm_center(0,0,0), axis_primary(0,0,0), axis_secondary(0,0,0) {}

	// initialize from a string, e.g., 'C4'
	DensitySymmInfo(std::string tag) : symm_center(0,0,0), axis_primary(0,0,0), axis_secondary(0,0,0) {
		type = tag[0];
		runtime_assert( type == 'C' || type == 'D' );
		count_primary = (core::Size)std::atoi( tag.substr(1).c_str() );
	}

	// detect symm operators from an electron density map
	void
	detect_axes( core::scoring::electron_density::ElectronDensity const &e );

	// mask a volume to the asymmetric unit
	void
	mask_asu(
		ObjexxFCL::FArray3D< double > &vol,
		core::scoring::electron_density::ElectronDensity const &e,
		double value );

	// min dist between X and any symm copy of Y
	core::Real
	min_symm_dist2(
		numeric::xyzVector< core::Real > const &X,
		numeric::xyzVector< core::Real > const &Y
	) const;

	// is symmetry enabled?
	bool
	enabled() const { return (count_primary!=1); }

private:
	char type; // one of 'C' or 'D' (currently)

	core::Size count_primary;
	numeric::xyzVector< core::Real > symm_center;
	numeric::xyzVector< core::Real > axis_primary;
	numeric::xyzVector< core::Real > axis_secondary; // 2-fold axis for D symmetries

	// helper function
	// autocorrelate with  rotation
	// find: symm_center _and_ autocorrelation
	void
	autocorrelate(
		core::scoring::electron_density::ElectronDensity const &e,
		core::Real angle,
		numeric::xyzVector< core::Real > const &axis,
		numeric::xyzVector< core::Real >  &symm_center,
		core::Real &autocorrelation
	);

	// helper function
	// resolve two symm axes to find center
	numeric::xyzVector< core::Real >
	resolve_symm_axes(
		numeric::xyzVector< core::Real > const &c1, numeric::xyzVector< core::Real > const &a1,
		numeric::xyzVector< core::Real > const &c2, numeric::xyzVector< core::Real > const &a2,
		core::Real &error );

};


} // moves
} // protocols

#endif
