// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_core_scoring_orbitals_OrbitalsLookup_hh
#define INCLUDED_core_scoring_orbitals_OrbitalsLookup_hh

#include <utility/vector1.hh>

#include <core/types.hh>

//project headers
#include <core/chemical/orbitals/OrbitalTypeMapper.fwd.hh>

#include <numeric/interpolation/spline/Bicubic_spline.hh>
#include <numeric/interpolation/spline/Cubic_spline.hh>
// AUTO-REMOVED #include <numeric/interpolation/spline/Interpolator.hh>
//#include <map>

//STL header

//#include <map>
#include <vector>

namespace core {
namespace scoring {
namespace orbitals {

class OrbitalsLookup {
public:
	enum h_type { hpol, sc_orb_bb_H, haro, bb };


	OrbitalsLookup( utility::vector1< std::string > const & filename, utility::vector1< std::string > const & cubic_files );

	void get_dist_AOH_angle_energy (
		const h_type h_enum,
		const core::chemical::orbitals::orbital_type_enum orb_type_name,
		const core::Real distance,
		const core::Real AOH_angle,
		core::Real & energy,
		core::Real & distance_derivative,
		core::Real & angle_derivative,
		bool check_derivative
	) const;


	bool check_distance(
		const core::Real distance,
		const h_type h_enum,
		const core::chemical::orbitals::orbital_type_enum orb_type_name
	) const;


	bool check_AOH_angle(
		const core::Real AOH_angle,
		const h_type h_enum,
		const core::chemical::orbitals::orbital_type_enum orb_type_name
	) const;


	void get_DHO_angle_energy(
			const h_type h_enum,
			const core::chemical::orbitals::orbital_type_enum orb_type_name,
			const core::Real DHO_angle,
			core::Real & energy,
			core::Real & angle_derivative,
			bool check_derivative
	) const;



private:
	///@brief number of statistics to put into matrix
	core::chemical::orbitals::orbital_type_enum number_stats_;
	///@brief number of elements in the KBP
	core::Size number_elements_;

	utility::vector1< numeric::interpolation::spline::BicubicSpline  > HPOL_sc_H_sc_orb_splines_;
	utility::vector1< numeric::interpolation::spline::BicubicSpline  > HPOL_bb_H_sc_orb_splines_;
	utility::vector1< numeric::interpolation::spline::BicubicSpline  > HPOL_sc_H_bb_orb_splines_;
	utility::vector1< numeric::interpolation::spline::BicubicSpline  > HARO_sc_H_sc_orb_splines_;


	utility::vector1<numeric::interpolation::spline::CubicSpline> cubic_HPOL_sc_H_sc_orb_splines_;
	utility::vector1<numeric::interpolation::spline::CubicSpline> cubic_HPOL_sc_H_bb_orb_splines_;
	utility::vector1<numeric::interpolation::spline::CubicSpline> cubic_HPOL_bb_H_sc_orb_splines_;
	utility::vector1<numeric::interpolation::spline::CubicSpline> cubic_HARO_sc_H_sc_orb_splines_;



/*	std::map< core::Size, numeric::interpolation::spline::BicubicSpline > hpol_splines_; //a map that contains all
	std::map< core::Size, numeric::interpolation::spline::BicubicSpline > haro_splines_; //a map that contains all
	std::map< core::Size, numeric::interpolation::spline::BicubicSpline > hpol_sc_orb_bb_splines_; //a map that contains all

	std::map< core::Size, numeric::interpolation::spline::BicubicSpline > orb_bb_splines_; //a map that contains all


	utility::vector1<numeric::interpolation::spline::BicubicSpline> HPOL_sc_H_sc_orb_splines_;
	utility::vector1<numeric::interpolation::spline::BicubicSpline> haro_spline_vector_;*/


	core::Real Cpisp2_min2_;
	core::Real Npisp2_min2_;
	core::Real Npsp2_min2_;
	core::Real Opisp2_min2_;
	core::Real Opsp2_min2_;
	core::Real Opsp3_min2_;
	core::Real Spsp3_min2_;


};



}
}
}



#endif /* ORBITALSLOOKUP_HH_ */
