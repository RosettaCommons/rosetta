// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRSpinlabel.hh
/// @brief   Class that contains data to describe an NMR spinlabel e.g. its chemical structure.
///          The NMRSpinlabel contains the spinlabel's ResidueType, name, code and other related data.
///          An NMRDummySpinlabelEnsemble object is used to represent the spinlabel's conformations.
///          Filtering of NMRDummySpinlabel conformers is done either by neighbor count or bump energy
///          calculation
/// @details last Modified: 08/31/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_NMRSpinlabel_HH
#define INCLUDED_core_scoring_nmr_NMRSpinlabel_HH

// Unit headers
#include <core/scoring/nmr/NMRSpinlabel.fwd.hh>

// Package headers
#include <core/scoring/nmr/NMRDummySpinlabelEnsemble.fwd.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// C++ headers
#include <iostream>
#include <string>

namespace core {
namespace scoring {
namespace nmr {

class NMRSpinlabel : public utility::pointer::ReferenceCount {

public: // Types

	typedef core::chemical::ResidueTypeOP ResidueTypeOP;
	typedef core::chemical::ResidueTypeCOP ResidueTypeCOP;
	typedef utility::vector1< std::pair< Real, Vector > > WeightCoordVector;

public:

	/// @brief type of spinlabel conformer filtering
	///        DISTANCE    = measure the distance of spinlabel side chain heavy atoms to the
	///                      NBR_ATOM in every neighborhood residue in the spinlabel
	///                      environment and count how many atoms are within the NBR_RADIUS
	///        BUMP_ENERGY = filter based on packer energy of a spinlabel conformer with all
	///                      neighboring residues and keep conformers with energy lower than
	///                      bump energy threshold
	enum CONFORMER_FILTER {
		DISTANCE = 1,
		BUMP_ENERGY = 2
	};

public: // Methods

	/// @brief construct from strings of residue type set and residue type;
	///        the other properties (e.g. radical atom) are looked up in the database
	NMRSpinlabel(
		std::string const & residue_type_set,
		std::string const & residue_type
	);

	/// @brief construct from ResidueType, radical atom name
	///        and NMRDummySpinlabelEnsemble
	NMRSpinlabel(
		ResidueTypeCOP residue_type,
		std::string const & radical_atom,
		NMRDummySpinlabelEnsembleCOP dummy_ensemble
	);

	/// @brief copy constructor
	NMRSpinlabel(NMRSpinlabel const & other);

	/// @brief assignment operator
	NMRSpinlabel &
	operator=(NMRSpinlabel const & rhs);

	/// @brief destructor
	~NMRSpinlabel();

	/// Getters
	std::string get_name() const { return name_; }
	std::string get_code() const { return three_letter_code_; }
	ResidueTypeCOP get_residue_type() const { return residue_type_; }
	std::string get_radical_atom() const { return radical_atom_; }
	std::string get_distance_potential_histogram_file() const { return distance_potential_histogram_file_; }
	NMRDummySpinlabelEnsembleCOP get_dummy_ensemble() const { return dummy_ensemble_; }
	NMRDummySpinlabelEnsembleOP get_dummy_ensemble() { return dummy_ensemble_; }
	WeightCoordVector const & get_radical_atom_coordinates() const { return weights_coordinates_table_; }
	Size get_max_ensemble_size() const { return max_ensemble_size_; }
	Size get_current_ensemble_size() const { return weights_coordinates_table_.size(); }
	void set_highres_conformer_filter_type(std::string const & filter_type);
	CONFORMER_FILTER get_highres_conformer_filter_type() const { return highres_conformer_filter_type_; }
	Real get_boltzman_factor() const { return boltzmann_factor_; }
	void set_boltzmann_factor(Real const kt) { boltzmann_factor_ = kt; }
	void set_path_to_distance_potential_histogram_file(std::string const & filename) { distance_potential_histogram_file_ = filename; }


	/// These are the two functions through which the user interfaces with the class
	/// to create an ensemble of valid spinlabel conformers (i.e. those that are expected
	/// to be seen given the current environment). In the first function, a simple distance
	/// check is performed and conformers that have heavy atoms within the NBR_RADIUS of
	/// other residues are marked as clashing. In the second function, the user provides a
	/// boolean mask for spinlabels which should be kept.

	/// @brief filter dummy spinlabel ensemble given the neighborhood
	///        of a particular target residue in the pose. Return a vector of
	///        each spinlabel conformer's weight and radical atom coordinates.
	///        Performs also clustering of coordinates internally such that the
	///        vector size does not exceed the maximal number of spinlabel conformers.
	WeightCoordVector filter_spinlabel_ensemble_by_distance_check(
		pose::Pose const & pose,
		Size const & target_resid
	);

	/// @brief filter dummy ensemble and keep those spinlabels provided in boolean mask
	WeightCoordVector filter_spinlabel_ensemble_by_mask(
		pose::Pose const & pose,
		Size const & target_resid,
		utility::vector1<bool> const & mask,
		utility::vector1<Real> const & scores
	);

	void show(std::ostream & TR);

private: // Methods

	/// @brief default constructor
	NMRSpinlabel();

	void init_spinlabel_name();

	/// @brief utility function to convert string to class specific enum
	void convert_string_to_conformer_filter_type(std::string const & filter_type);

	/// @brief Performs clustering of the spinlabel coordinates vector before setting the class member.
	///        The input vector doesn't store all spinlabel atoms but only the coordinates
	///        of the radical atom. Clustering is done based on input RMSD matrix.
	void
	cluster_conformers_and_set_weights_and_coordinates(
		WeightCoordVector & weights_coords,
		utility::vector1<utility::vector1<Real>> & rmsd_mat
	);

	void
	set_weights_and_coordinates(
		WeightCoordVector const & weights_coords
	);

	void init_radical_atom_weights_and_coordinates();

	/// @brief register command line options
	void register_options();
	void init_from_cml();

private: // Data

	std::string name_;
	std::string three_letter_code_;
	ResidueTypeCOP residue_type_;
	std::string radical_atom_;
	std::string distance_potential_histogram_file_;

	// A container of a minimal, dummy spinlabel conformer representation
	NMRDummySpinlabelEnsembleOP dummy_ensemble_;

	// Cartesian coordinates and weight of each spinlabel conformer.
	// The spinlabel may exist as an ensemble, so we store a vector here.
	// We need only the coordinates of the radical atom to calculate the PRE/PCS.
	// Thus, a one-dimensional vector is sufficient.
	// Some conformers may have a higher probability (i.e. are observed more frequently).
	// Thus, we also store a weight here.
	WeightCoordVector weights_coordinates_table_;

	// Maximal number of spinlabels. If ensemble size is bigger, conformers in the
	// counts_coordinates_table_ are binned based on the spinlabel ensemble's RMSD matrix.
	Size max_ensemble_size_;

	CONFORMER_FILTER highres_conformer_filter_type_;

	Real boltzmann_factor_;
};

} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_NMRSpinlabel_HH
