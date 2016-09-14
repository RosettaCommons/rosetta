// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/electron_density_atomwise/ElectronDensityAtomwise.hh
/// @brief  Scoring a structure against an electron density map
/// @author Fang-Chieh Chou

#ifndef INCLUDED_core_scoring_electron_density_atomwise_ElectronDensityAtomwise_HH
#define INCLUDED_core_scoring_electron_density_atomwise_ElectronDensityAtomwise_HH

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

// Utility headers
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.hh>
#include <utility/exit.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray3D.hh>

// C++ headers
#include <string>
#include <map>

//Auto Headers
#include <core/kinematics/RT.hh>
#include <utility/vector1_bool.hh>
#include <numeric/xyzMatrix.hh>


namespace core {
namespace scoring {
namespace electron_density_atomwise {

const core::Real MAX_FLT = 1e37;

class ElectronDensityAtomwise {
public:
	/// @brief constructor
	ElectronDensityAtomwise();

	/// @brief Is a map loaded?
	bool isMapLoaded() {
		return is_map_loaded;
	};

	/// @brief Load an MRC (="new-CCP4") density map
	void readMRCandResize();

	//Compute Normalization factor given a pose
	void compute_normalization( pose::Pose const & pose );

	//Pre-compute the unweighted score of atom in a gride
	void precompute_unweighted_score();

	//Return the score of a given residue
	core::Real
	residue_score( core::conformation::Residue const & rsd );

	//Return the gradient of an atom
	numeric::xyzVector< core::Real > atom_gradient( core::pose::Pose const &
		pose, core::Size const & rsd_id, core::Size const & atm_id );

private:

	numeric::xyzMatrix<core::Real> f2c, c2f, i2c, c2i;

	bool is_map_loaded, is_score_precomputed;

	core::Real map_reso, cell_volume, r_cell_volume, max_del_grid, grid_spacing, gaussian_max_d, normalization, avg_rho_obs;

	//Stored 1d gaussian
	utility::vector1< core::Real > atom_gaussian_value;

	// the density data and precomputed score array
	ObjexxFCL::FArray3D< float > density;
	ObjexxFCL::FArray3D< double > unweighted_score_coeff;
	//Stored atom weight
	utility::vector1< utility::vector1 < core::Size > > atom_weight_stored;

	// map info
	numeric::xyzVector< int > grid;
	numeric::xyzVector< core::Real > orig;
	numeric::xyzVector< float > cell_angles, cell_dimensions, r_cell_angles,
		cos_r_cell_angles, r_cell_dimensions;
	//symmetry
	utility::vector1< core::kinematics::RT > symmOps;

	//helper function for symmetry
	void computeCrystParams();
	void expandToUnitCell();
	void initializeSymmOps( utility::vector1< std::string > const & symList );

	//resize the map
	void resize( core::Real approxGridSpacing );

	//////////////////////////////////////////////////////////////////////
	//compute index to cartesian transformation
	void calculate_index2cart();

	//generate 1d gaussian function and store it.
	void generate_gaussian_1d( core::Real const & sigma );
	//Return the value of 1D gaussian given the distance using stored values
	core::Real gaussian_1d( core::Real const & dist );

	//return the weight of atom given its element type
	core::Size get_atom_weight( std::string const & elt );

	//Spline interpolation
	core::Real spline_interpolation( ObjexxFCL::FArray3D < double >
		& coeffs, numeric::xyzVector< core::Real > const & idxX ) const;

	void spline_coeffs( ObjexxFCL::FArray3D< double > & data,
		ObjexxFCL::FArray3D< double > & coeffs );

	//Trilinear Interpolation
	core::Real trilinear_interpolation( ObjexxFCL::FArray3D< double > & score,
		numeric::xyzVector< core::Real > const & index );

	numeric::xyzVector<core::Real> trilinear_gradient( ObjexxFCL::FArray3D
		<double> & score, numeric::xyzVector< core::Real > const & index );

	//Convert a vector from xyz coordinate to index coordinate, shift w/ respect to the origin and fold into the unit cell
	numeric::xyzVector< core::Real > xyz2index_in_cell( numeric::xyzVector< core::Real > const & xyz_vector );

	//Convert a vector from fractional coordinate to index coordinate
	numeric::xyzVector< core::Real > frac2index( numeric::xyzVector< core::Real > const &frac_vector );

	//Convert a vector from index coordinate to fractional coordinate
	numeric::xyzVector< core::Real > index2frac( numeric::xyzVector< core::Real > const & frac_vector );
	numeric::xyzVector< core::Real > index2frac( numeric::xyzVector< int > const & frac_vector );
};

//@brief The EDM instance
ElectronDensityAtomwise& get_density_map();

} // electron_density_atomwise
} // scoring
} // core


#endif
