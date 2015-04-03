// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/Minimizer.cc
/// @brief  Minimizer class
/// @author Phil Bradley


// Unit headers
#include <core/optimization/LineMinimizer.hh>
#include <core/optimization/NelderMeadSimplex.hh>
#include <core/optimization/GA_Minimizer.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <utility/exit.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>

// C++ headers
#include <cmath>
#include <iostream>
#include <algorithm>

#include <utility/vector1.hh>

#ifdef WIN32
#include <functional>
#endif


namespace core {
namespace optimization {

using namespace ObjexxFCL;
static thread_local basic::Tracer TR( "core.optimization.NelderMeadSimplex" );

// sort our vector list min->max
struct sort_pred {
    bool operator()(const std::pair< core::Real, Multivec > &left, const std::pair< core::Real, Multivec > &right) {
        return left.first < right.first;
    }
};


// set the function and the options
NelderMeadSimplex::NelderMeadSimplex(
	Multifunc & func_in,
	MinimizerOptions const & options_in
) : func_( func_in ), options_( options_in ) {}


Real
NelderMeadSimplex::run(
	Multivec & phipsi_inout,
	Multivec const& upperbound
) {
	int const ITMAX( options_.max_iter() );
	core::Real ALPHA=1, GAMMA=2, RHO=-0.5, SIGMA=0.5;

	core::Size ndim = phipsi_inout.size();
	Size npoints = ndim+1;

	utility::vector1< std::pair< core::Real, Multivec > > vertices(npoints, std::make_pair( 0.0, Multivec(ndim,0.0)) );

	// populate and evaluate vertex list
	for (Size i=1; i<=npoints; ++i) {
		Multivec &m_i = vertices[i].second;
		for (Size j=1; j<=ndim; ++j) {
			if (i==j) {
				m_i[j] = upperbound[j];
			} else {
				m_i[j] = phipsi_inout[j];
			}
		}
		//std::cerr << "bound (dim " << i << ")" << std::endl;
		vertices[i].first = func_(m_i);
	}

	for (int i=1; i<=ITMAX; ++i) {
		std::sort(vertices.begin(), vertices.end(), sort_pred() );

		// 1 - median of all pts but the worst
		Multivec x0(ndim,0.0);
		for (Size j=1; j<=npoints-1; ++j) {
			Multivec &m_j = vertices[j].second;
			for (Size k=1; k<=ndim; ++k) {
				x0[k]+=m_j[j];
			}
		}
		for (Size k=1; k<=ndim; ++k) {
			x0[k]/=(npoints-1);
		}

		// 2 reflected point
		Multivec xR(ndim);
		Multivec &xN = vertices[npoints].second;
		for (Size k=1; k<=ndim; ++k)
			xR[k] = x0[k] + ALPHA*(x0[k]-xN[k]);
		//std::cerr << "reflect worst=" << vertices[npoints].first << std::endl;
		core::Real yR = func_(xR);

		if (yR < vertices[1].first) {
			// 3A new best? expansion
			Multivec xE(ndim);
			for (Size k=1; k<=ndim; ++k)
				xE[k] = x0[k] + GAMMA*(x0[k]-xN[k]);
			//std::cerr << "expand new best =" << yR << std::endl;
			core::Real yE = func_(xE);
			if (yE<yR) {
				vertices[npoints] = std::make_pair( yE,xE );
			} else {
				vertices[npoints] = std::make_pair( yR,xR );
			}
		} else if (yR < vertices[npoints-1].first) {
			// 3B better than second worst? replacement
			vertices[npoints] = std::make_pair( yR,xR );
		} else {
			// 3C no improvement --- contraction
			Multivec xC(ndim);
			for (Size k=1; k<=ndim; ++k)
				xC[k] = x0[k] + RHO*(x0[k]-xN[k]);
			//std::cerr << "contract new worst =" << yR << std::endl;
			core::Real yC = func_(xC);
			if (yC<vertices[npoints].first) {
				vertices[npoints] = std::make_pair( yC,xC );
			} else {
				// 4 collapse on best
				for (Size j=2; j<=npoints; ++j) {
					Multivec &m_1 = vertices[1].second;
					for (Size k=1; k<=ndim; ++k) {
						vertices[j].second[k] = m_1[k] + SIGMA*(vertices[j].second[k]-m_1[k]);
					}
					//std::cerr << "collapse point " << j << " " << vertices[j].first << std::endl;
					vertices[j].first = func_(vertices[j].second);
				}
			}
		}

		// to do: convergence conditions
	}

	phipsi_inout = vertices[1].second;
	return vertices[1].first;
}

} // namespace optimization
} // namespace core
