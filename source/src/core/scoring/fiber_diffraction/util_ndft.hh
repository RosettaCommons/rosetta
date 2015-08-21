// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffraction/FiberDiffractionEnergyDens.hh
/// @brief  FiberDiffractionDens utility scoring
/// @author Wojciech Potrzebowski and Ingemar Andre


#ifndef INCLUDED_core_scoring_fiber_diffraction_util_ndft_hh
#define INCLUDED_core_scoring_fiber_diffraction_util_ndft_hh

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <complex>

namespace core {
namespace scoring {
namespace fiber_diffraction {

void
ft_nfft(
	ObjexxFCL::FArray3D< float > & fourier_in,
	utility::vector0< utility::vector0 < int > >::iterator & nvals,
	core::Size & lmax,
	core::Size & total_rvals,
	ObjexxFCL::FArray3D< std::complex<float> > & fourier_out
);

void
ndft(
	ObjexxFCL::FArray1D< std::complex <float> > & fin,
	ObjexxFCL::FArray1D< float > & phi,
	ObjexxFCL::FArray1D< int > & n,
	ObjexxFCL::FArray1D< std::complex <float> > & fout
);

void
FBnl_R(
	ObjexxFCL::FArray3D< float > & fourier,
	utility::vector0< utility::vector0 < int > >::iterator & nvals,
	core::Size & lmax,
	core::Size & rindex,
	ObjexxFCL::FArray3D< std::complex<float> > & fourier_out
);

void
ndft_1d(
	ObjexxFCL::FArray1D < std::complex< float> >  & f1,
	ObjexxFCL::FArray1D< float > & x0,
	ObjexxFCL::FArray1D< std::complex< float> > & f2
);

}
}
}

#endif
