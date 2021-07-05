// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/EPRSpinLabel.hh
/// @brief  This is a container class specific to the use of double electron-electron resonance data
/// @details This container manages the simulated distance distributions for the deer_decay and
///      deer_distance energy method. It also manages individual electron coordinate ensembles for
///      a given protein, although a reference to that pose is not stored in this object.
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_EPRSpinLabel_hh
#define INCLUDED_core_scoring_epr_deer_EPRSpinLabel_hh

// Unit headers
#include <core/scoring/epr_deer/EPRSpinLabel.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/AtomVDW.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>

#include <map>

namespace core {
namespace scoring {
namespace epr_deer {

/// @brief Alias to save space throughout core/scoring/epr_deer
using PseudoSL = std::pair< numeric::xyzVector< Real >, Real >;

/// @brief Alias to asve space when defining residues and spin labels
using PairSizeString = std::pair< Size, std::string >;

/// @brief Object use to calculate VDW radii of centroid atoms
/// @details Global to avoid re-instantiating with every new object
static scoring::AtomVDWOP atom_vdw_ = scoring::AtomVDWOP( nullptr );

/// @brief Default coordinates map
/// @details Global to avoid re-reading from database with every new object
static std::map< std::string, utility::vector1< PseudoSL > > deflt_coords_;

/// @brief Virtual, idealized CB atom for center-of-mass calculation
/// @details Global to avoid reinstantiation with every new object
static const numeric::xyzVector< Real > vrt_cb_( -1.3116, -0.5715, 0.5398 );

class EPRSpinLabel {

public:

	/// @brief Constructor
	EPRSpinLabel();

	/// @brief Destructor
	~EPRSpinLabel();

	/// @brief Initialize object used to calculate centroid clashes
	void
	init_vdw() const;

	/// @brief Read DB file for a given type of spin label
	/// @param Name of file/SL type
	/// @return PseudoSLs used for simulation of DEER distributions
	utility::vector1< PseudoSL >
	read_db_file(
		std::string const & name
	);

	/// @brief Operator to return nnn-const PseudoSL from specific residue
	/// @param res: Residue info (number and spin label type)
	/// @return Coords corresponding to residue
	utility::vector1< PseudoSL > &
	operator[]( PairSizeString const & res );

	/// @brief Return const PseudoSL from specific residue
	/// @param res: Residue info (number and spin label type)
	/// @return Coords corresponding to residue
	utility::vector1< PseudoSL > const &
	at( PairSizeString const & res ) const;

	/// @brief Returns histogram between coordinate sets for residues
	/// @param residues: Residues contributing to histogram
	/// @param bins_per_a: Granularity of histogram (bins per angstrom)
	/// @param mod: What to add to the X-axis
	/// @param dist_ids: If a custom X-axis is used (default: empty)
	/// @return Histogram with X- and Y-values being keys and values
	/// @detail Note: the X-axis = bins_per_a * distance, rounded to an int
	/// @detail Note: Equal labeling assumed for all sites
	/// @detail This matters if residues.size() > 2
	std::map< Size, Real >
	histogram(
		utility::vector1< PairSizeString > const & residues,
		Size const & bins_per_a,
		int const & mod = 0,
		Real const & stdev = 1.0,
		std::map< Size, Real > const & dist_ids = {}
	);

	/// @brief Return histogram for pair of coordinate sets
	/// @param res1_coords: Pair of PseudoSL coords for res1
	/// @param res2_coords: Pair of PseudoSL coords for res2
	/// @param bins_per_a: Bins per angstrom for distribution
	/// @param mod: How much to shift X-axis
	/// @param stdev: St deviation of gauss used to convolute pairwise dists
	/// @param dist_ids: If a custom X-axis is used (default: empty)
	/// @return Histogram with X- and Y-values being keys and values
	/// @detail Note: the X-axis = bins_per_a * distance, rounded to an int
	std::map< Size, Real >
	histogram(
		utility::vector1< PseudoSL > const & res1_coords,
		utility::vector1< PseudoSL > const & res2_coords,
		Size const & bins_per_a,
		int const & mod = 0,
		Real const & stdev = 1.0,
		std::map< Size, Real > const & dist_ids = {}
	);

	/// @brief Returns histogram between coordinate sets for residues
	/// @param res1: Residue 1
	/// @param res2: Residue 2
	/// @param bins_per_a: Granularity of histogram (bins per angstrom)
	/// @param mod: What to add to the X-axis
	/// @param dist_ids: If a custom X-axis is used (default: empty)
	/// @return Histogram with X- and Y-values being keys and values
	/// @detail Note: the X-axis = bins_per_a * distance, rounded to an int
	std::map< Size, Real >
	histogram(
		PairSizeString const & res1,
		PairSizeString const & res2,
		Size const & bins_per_a,
		int const & mod = 0,
		Real const & stdev = 1.0,
		std::map< Size, Real > const & dist_ids = {}
	);

	/// @brief Label a residue with a certain spin label
	/// @param res: Residue index
	/// @param label: SL type
	/// @param pose: Pose used for superimposition
	/// @param skip_clash_eval: Whether clash evaluation is skipped
	/// @detail Clashes are skipped for custom coordinates for reasons
	/// @detail  relating to the way they have been calculated
	void
	label(
		PairSizeString const & res_label,
		pose::Pose const & pose,
		bool const & skip_clash_eval = false
	);

	/// @brief Normalize distribution so that the sum is equal to 1.0
	/// @param sim_map: Simulated DEER distribution
	/// @result Identical std::map except values add up to 1.0
	std::map< Size, Real >
	normalize(
		std::map< Size, Real > sim_map
	) const;

	/// @brief Get positions of unpaired electrons at a residue
	/// @param res: Residue number
	/// @param pose: Pose for clash eval
	/// @param sl_vec: Vector of PseudoSLs, which have positions of unpaired e
	/// @param skip_clash_eval: Exactly what it suggested by the title
	/// @param min_rad: Lowest radius for clash eval to check
	/// @return Vector of PseudoSLs in local coordinate frame of residue
	utility::vector1< PseudoSL >
	calc_sl_for_res(
		Size const & res,
		pose::Pose const & pose,
		utility::vector1< PseudoSL > const & sl_vec,
		bool const & skip_clash_eval = false,
		Real const & min_rad = 0.0
	);

	/// @brief retrieve weight for given coordinate
	/// @param res1: Residue over which the coordinate is being superimposed
	/// @param clash_xyz: Coordinate used for clash calculation
	/// @param w: Weight, passed by value since we need a new obj to modify
	/// @param pose: Pose with all the residues we check for clash evaluation
	/// @param vdw_rad: Radius of the clash_xyz atom to consider
	/// @return Weight of PseudoSL at position given local environment of pose
	Real
	weight(
		Size const & res1,
		numeric::xyzVector< core::Real > const & clash_xyz,
		Real w,
		pose::Pose const & pose,
		Real const & rad
	);

	/// @brief Getter for cutoff for weights
	/// @return Weight cutoff
	Real
	cutoff() const;

	/// @brief Allows a custom set of electrons to be read without superimposition
	/// @param all_coords: Custom residue-specific PseudoSLs
	void
	load_custom_electrons(
		std::map< Size, utility::vector1< PseudoSL > > const & all_coords
	);

	/// @brief Goes through every residue in provided list and calculates
	/// @param pose: Pose to label
	/// @param residues: Residues that need to be labeled
	/// @param skip_clash_eval: Whether clash evaluation should be skipped
	void
	label_everything(
		pose::Pose & pose,
		utility::vector1< PairSizeString > const & residues,
		bool const & skip_clash_eval
	);

private:

	/// @brief Returns the center of mass between coordinate and CB
	/// @param coord: Coordinate xyz
	/// @param cb: CB coordinate
	/// @return XYZ of bulk / center of mass of nitroxide ring
	/// @detail See del Alamo et al 2020 Biophysical Journal for details on
	///   why 0.875 was computed/chosen
	numeric::xyzVector< Real >
	bulk(
		numeric::xyzVector< Real > const & coord,
		numeric::xyzVector< Real > const & cb
	) const;

private: // data

	/// @brief Coordinates for each residue being stored
	std::map< PairSizeString, utility::vector1< PseudoSL > > coords_;

	/// @brief Custom coordinates (residue-specific)
	std::map< Size, utility::vector1< PseudoSL > > custom_coords_;

	/// @brief Cutoff for clash evaluation
	Real cutoff_ = 0.001;

	/// @brief Path for files describing default PseudoSL objects to read
	std::string path_ = "scoring/epr_deer/";

};

} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
