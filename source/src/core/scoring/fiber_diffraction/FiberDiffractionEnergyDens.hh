// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffraction/FiberDiffractionEnergyDens.hh
/// @brief  FiberDiffractionDens scoring (low-resolution, electron density based)
/// @author Wojciech Potrzebowski and Ingemar Andre


#ifndef INCLUDED_core_scoring_fiber_diffraction_FiberDiffractionEnergyDens_hh
#define INCLUDED_core_scoring_fiber_diffraction_FiberDiffractionEnergyDens_hh

#include <core/scoring/fiber_diffraction/FiberDiffractionEnergyDensCreator.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>
#include <core/id/AtomID.hh>

#include <utility/vector1.hh>
#include <utility/vector0.hh>

#include <numeric/xyzVector.hh>
#include <numeric/interpolation/spline/Interpolator.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>

#include <string>
#include <complex>

namespace core {
namespace scoring {
namespace fiber_diffraction {

class FiberDiffractionEnergyDens : public methods::WholeStructureEnergy  {
public:
	typedef methods::WholeStructureEnergy  parent;

public:

	FiberDiffractionEnergyDens();
	FiberDiffractionEnergyDens( FiberDiffractionEnergyDens const & src ) :
		parent(src) {
		chi2_ = src.chi2_;
		scale_factor_ = src.scale_factor_;
		dchi2_d = src.dchi2_d;
		dchi2_d_cross_R = src.dchi2_d_cross_R;
	}

	virtual methods::EnergyMethodOP clone() const;

	virtual void finalize_total_energy(pose::Pose & pose,ScoreFunction const &,EnergyMap & totals) const;

	virtual void setup_for_scoring( pose::Pose & pose, ScoreFunction const & scorefxn ) const;

	/*virtual void setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sf) const;

	virtual void
	eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
	) const;*/

	virtual void indicate_required_context_graphs(utility::vector1< bool > & /*context_graphs_required*/
	) const {}


	void
	calculate_rho_fast2(
		pose::Pose & pose,
		ObjexxFCL::FArray1D< float > & rc,
		ObjexxFCL::FArray1D< float > & phic,
		ObjexxFCL::FArray1D< float > & zc,
		Real const c_
	) const;

private:
	//////
	// PRIVATE DATA becuase the functions in this class are const we need to make variables mutable
	// store calculated intensities
	mutable utility::vector0< utility::vector1 < Real > > I;
	mutable utility::vector1< numeric::xyzVector< core::Real > > dchi2_d, dchi2_d_cross_R;
	mutable std::map<  core::id::AtomID, core::Size > AtomID_to_atomnbr_;
	mutable core::Real chi2_;
	mutable core::Real scale_factor_;
	mutable core::Real square_obs_, sum_obs_;
	mutable core::Real p_, c_;
	mutable core::Real a_, b_;
	mutable ObjexxFCL::FArray3D< float > rho_cylindrical_;

	virtual
	core::Size version() const;
};


void
gnl_R_qfht(
	ObjexxFCL::FArray3D< float > & fourier,
	utility::vector0 < utility::vector0 < int > >::iterator & nvals,
	core::Size & lmax,
	core::Real & max_r_value,
	double & qfht_K1,
	double & qfht_K2,
	ObjexxFCL::FArray3D< std::complex<float> > & Gnl
);


void
find_max_r(
	pose::Pose const & pose,
	core::Real & maxR
);

utility::pointer::shared_ptr< numeric::interpolation::spline::Interpolator >
fit_layer_lines_with_splines(
	ObjexxFCL::FArray1D < float > xvals,
	ObjexxFCL::FArray1D < float > yvals
);

}
}
}

#endif
