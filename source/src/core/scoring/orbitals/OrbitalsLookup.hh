// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_core_scoring_orbitals_OrbitalsLookup_hh
#define INCLUDED_core_scoring_orbitals_OrbitalsLookup_hh

// Unit headers
#include <core/scoring/orbitals/OrbitalsLookup.fwd.hh>

#include <utility/vector1.hh>

#include <core/types.hh>

#include <core/chemical/orbitals/OrbitalTypeMapper.fwd.hh>

#include <numeric/interpolation/spline/BicubicSpline.hh>
#include <numeric/interpolation/spline/CubicSpline.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <map>
#include <string>

#include <vector>
#include <core/scoring/ScoreFunction.fwd.hh>

namespace core {
namespace scoring {
namespace orbitals {

class OrbitalsLookup : public utility::pointer::ReferenceCount {
public:
	enum h_type { Hpol_scOrbH, Haro_scOrbH, Hpol_bbOrbH };


	OrbitalsLookup(
		utility::vector1< std::string > const & DHO_energies,
		utility::vector1< std::string > const & AOH_energies,
		utility::vector1< std::string > const & AOO_energies,
		utility::vector1< std::string > const & DOO_energies,
		utility::vector1< std::string > const & ACO_AOH_orb_Hpol_energies
	);

	utility::vector1< utility::vector1< core::Real > > parse_files(
		std::string const & file,
		std::map<core::Size, std::pair<core::Size, core::Size> > & orbital_angle_dist_map
	)const;

	void OrbHdist_cosDHO_energy (
		const h_type h_enum,
		const core::Size orb_type_name,
		const core::Real distance,
		const core::Real AOH_angle,
		core::Real & energy,
		core::Real & distance_derivative,
		core::Real & angle_derivative,
		bool check_derivative
	) const;

	void OrbHdist_cosAOH_energy
	(
		const h_type h_enum,
		const core::Size orb_type_name,
		const core::Real distance,
		const core::Real AOH_angle,
		core::Real & energy,
		core::Real & distance_derivative,
		core::Real & angle_derivative,
		bool check_derivative,
		bool ACO
	) const;

	void OrbOrbDist_cosAOD_energy(
		const core::Size orb_type_name1,
		const core::Size orb_type_name2,
		const core::Real distance,
		const core::Real AOO_angle,
		core::Real & energy,
		core::Real & distance_derivative,
		core::Real & angle_derivative,
		bool check_derivative
	)const;
	void OrbOrbDist_cosDOA_energy(
		const core::Size orb_type_name1,
		const core::Size orb_type_name2,
		const core::Real distance,
		const core::Real DOO_angle,
		core::Real & energy,
		core::Real & distance_derivative,
		core::Real & angle_derivative,
		bool check_derivative

	)const;


private:
	/// @brief number of statistics to put into matrix
	core::Size number_stats_;
	/// @brief number of elements in the KBP
	core::Size number_elements_;

	utility::vector1< numeric::interpolation::spline::BicubicSpline  > DHO_Hpol_scOrbH_splines_;
	utility::vector1< numeric::interpolation::spline::BicubicSpline  > DHO_Haro_scOrbH_splines_;
	utility::vector1< numeric::interpolation::spline::BicubicSpline  > DHO_Hpol_bbOrbH_splines_;

	utility::vector1< numeric::interpolation::spline::BicubicSpline  > AOH_Hpol_scOrbH_splines_;
	utility::vector1< numeric::interpolation::spline::BicubicSpline  > AOH_Haro_scOrbH_splines_;
	utility::vector1< numeric::interpolation::spline::BicubicSpline  > AOH_Hpol_bbOrbH_splines_;


	utility::vector1< numeric::interpolation::spline::BicubicSpline  > AOD_orb_orb_splines_;
	utility::vector1< numeric::interpolation::spline::BicubicSpline  > DOA_orb_orb_splines_;
	utility::vector1< numeric::interpolation::spline::BicubicSpline  > ACO_AOH_orb_Hpol_splines_;

	// unused mutable core::Real scOrb_scHpol_weight_;
	// unused mutable core::Real scOrb_scOrb_weight_;

};


}
}
}


#endif /* ORBITALSLOOKUP_HH_ */
