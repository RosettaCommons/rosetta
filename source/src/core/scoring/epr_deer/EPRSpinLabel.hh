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
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>

#include <set>

namespace core {
namespace scoring {
namespace epr_deer {

typedef std::pair< numeric::xyzVector< Real >, Real > PseudoElectron;

static
	std::map< Size, utility::vector1< PseudoElectron > > custom_coords_;

static const
numeric::xyzVector< Real > cb_coord_( -1.3116, -0.5715, 0.5398 );

class EPRSpinLabel {

public:

	/// @brief Constructor
	EPRSpinLabel();

	/// @brief Destructor
	~EPRSpinLabel();

	/// @brief Operator to return electrons from specific residue
	utility::vector1< PseudoElectron > &
	operator[]( std::pair< Size, std::string > const & res );

	/// @brief Allows const spin label data to be accessed
	utility::vector1< PseudoElectron > const &
	at( std::pair< Size, std::string > const & res ) const;

	/// @brief Returns a histogram between all coordinates for all residues - assumes complete labeling
	std::map< Size, Real >
	histogram(
		utility::vector1< std::pair< Size, std::string > > const & residues,
		Size const & bins_per_a,
		Real const & modifer = 0.0
	);

	/// @brief Returns a histogram between all coordinates between two residue/SL combinations
	std::map< Size, Real >
	histogram(
		std::pair< Size, std::string > const & res1,
		std::pair< Size, std::string > const & res2,
		Size const & bins_per_a,
		Real const & modifer = 0.0,
		Real const & stdev = 1.0
	);

	/// @brief Return a value for a gaussian distribution at a particular value, given a average and standard deviation
	Real
	gauss(
		Real const & dist,
		Real const & avg,
		Real const & stdev
	) const;

	/// @brief Label a residue with a certain spin label
	void
	label(
		Size const & res,
		std::string const & label,
		pose::Pose const & pose,
		bool const & skip_clash_eval = false
	);

	/// @brief Normalize distribution so that the sum is equal to one
	std::map< Size, Real >
	normalize_distribution(
		std::map< Size, Real > sim_map
	) const;

	/// @brief Given a set of electrons, a pose, and a residue of interest, find viable coords
	utility::vector1< PseudoElectron >
	get_coords(
		core::Size const & res,
		pose::Pose const & pose,
		utility::vector1< PseudoElectron > const & electrons,
		bool const & skip_clash_eval = false
	);

	/// @brief retrieve weight for given coordinate
	Real
	get_weight(
		core::Size const & source_res,
		numeric::xyzVector< core::Real > const & center_of_mass,
		pose::Pose const & pose,
		Real const & forgive_factor
	);

	/// @brief Retrieve cuttof for weights (if the weight is less than this, it is set to zero)
	Real
	cutoff() const;

	/// @brief Allows a custom set of electrons to be read without superimposition
	void
	load_custom_electrons(
		std::map< Size, utility::vector1< PseudoElectron > > const & all_coords
	);

private:

	/// @brief   Returns the center of mass between a given electron coordinate and its CB
	/// @details  The precise value 0.875 descrbed the average value that is crystallographically observed
	///      in MTSSL rotamers for the distance between the CB and the electron vs the distance between
	///      the CB and the nitroxide ring center of mass, from which clashes are evaluated here.
	///      Described in detail in a forthcoming publication
	numeric::xyzVector< Real >
	center_of_mass(
		numeric::xyzVector< Real > const & electron,
		numeric::xyzVector< Real > const & cb
	) const;

private: // data

	std::map< std::pair< Size, std::string >, utility::vector1< PseudoElectron > > mapped_coords_;
	utility::vector1< PseudoElectron > custom_electrons_;
	std::map< Size, utility::vector1< PseudoElectron > > custom_coords_;
	Real cutoff_ = 0.001;

};

} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
