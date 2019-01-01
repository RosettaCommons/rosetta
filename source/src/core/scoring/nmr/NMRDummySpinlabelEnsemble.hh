// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRDummySpinlabelEnsemble.hh
/// @brief   NMRDummySpinlabelEnsemble stores a set of conformers of a particular spin-label residue,
///          and other related data such as the similarity (Rmsd) matrix between spin-label conformers,
///          the RB transformation needed to align the spin-label onto a selected target site, and a
///          voxel grid to look up clashes of a spin-label conformer with the protein environment.
/// @details last Modified: 11/28/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_NMRDummySpinlabelEnsemble_HH
#define INCLUDED_core_scoring_nmr_NMRDummySpinlabelEnsemble_HH

// Unit headers
#include <core/scoring/nmr/NMRDummySpinlabelEnsemble.fwd.hh>

// Package headers
#include <core/scoring/nmr/NMRDummySpinlabelVoxelGrid.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/HomogeneousTransform.hh>

// boost headers
#include <boost/unordered/unordered_map.hpp>
#include <boost/functional/hash.hpp>

// C++ headers
#include <iostream>
#include <string>
#include <map>

namespace core {
namespace scoring {
namespace nmr {

typedef boost::unordered_map< std::string, Vector > AtomNamePosMap;
typedef boost::unordered_map< std::string, NMRDummySpinlabelAtom > NMRDummySpinlabelAtomTable;
typedef boost::unordered_map< std::string, NMRDummySpinlabelAtom >::iterator NMRDummySpinlabelAtomTableIter;
typedef boost::unordered_map< std::string, NMRDummySpinlabelAtom >::const_iterator NMRDummySpinlabelAtomTableCOIter;
typedef std::pair< std::string, std::string > AtomPairKey;
typedef core::conformation::ResidueOP ResidueOP;
typedef core::conformation::ResidueCOP ResidueCOP;

/// @brief Class that stores information about one dummy spinlabel conformer (e.g. atom names and xyz
///        coordinates, and whether it clashes with neighboring amino acid residues or not). A vector of
///        NMRDummySpinlabelConformer objects is member of class NMRDummySpinlabelEnsemble.
class NMRDummySpinlabelConformer : public utility::pointer::ReferenceCount {

private: //Methods

	/// @brief default constructor
	NMRDummySpinlabelConformer();

public: // Methods

	/// @brief Construct from ID, number of observations, frequency and residue
	NMRDummySpinlabelConformer(
		Size const id,
		Size const nobs,
		Real const freq,
		conformation::Residue const & residue
	);

	/// @brief Construct from ID, number of observations (frequency),
	///        atom names, xyz coordinates and residue type
	NMRDummySpinlabelConformer(
		Size const id,
		Size const nobs,
		Real const freq,
		chemical::ResidueType const & restype,
		utility::vector1< std::string > const & names,
		utility::vector1< Vector > const & coords
	);

	NMRDummySpinlabelConformer(
		Size const id,
		Size const nobs,
		Real const freq,
		chemical::ResidueType const & restype,
		AtomNamePosMap const & names_coords
	);

	/// @brief copy constructor
	NMRDummySpinlabelConformer(NMRDummySpinlabelConformer const & other);

	/// @brief assignment operator
	NMRDummySpinlabelConformer &
	operator=(NMRDummySpinlabelConformer const & rhs);

	/// @brief destructor
	~NMRDummySpinlabelConformer();

	/// Getters and Setters
	Size get_id() const { return id_; }
	Size get_nobs() const { return nobs_; }
	Real get_frequency() const { return frequency_; }
	Real & frequency() { return frequency_; }
	bool has_clash() const { return clash_; }
	Real get_clash_score() const { return clash_score_; }
	Real & clash_score() { return clash_score_; }
	NMRDummySpinlabelAtomTable & get_atom_table() { return atom_table_; }
	NMRDummySpinlabelAtomTable const & get_atom_table() const { return atom_table_; }
	ResidueOP get_residue() { return residue_; }
	ResidueCOP get_residue() const { return residue_; }

	void clash_on() { clash_ = true; }
	void clash_off() { clash_ = false; }
	void set_clash_score(Real score) { clash_score_ = score; }

private: // Data

	// The ID number of this NMRDummySpinlabelConformer.
	// This is used to identify it at certain places, e.g. during the
	// clash check. If an NMRDummySpinlabelConformer has been marked
	// as clashing, it is not necessary to test another atom of the
	// same conformer.
	Size id_;

	// Number of observations of this NMRDummySpinlabelConformer.
	// We can use this to give certain conformers higher weight (see below).
	// At the moment, this information is not contained in the
	// conformer file yet, but this might change in the future
	// and we could read this information from there.
	Size nobs_;

	// Weight of this NMRDummySpinlabelConformer.
	// Initially this is nobs/Sum(nobs) where the sum runs over all
	// NMRSpinlabelConformers. When clustering NMRSpinlabelConformers,
	// the cluster centroids have a frequency variable that is proportional
	// to the cluster size
	Real frequency_;

	// Does this NMRDummySpinlabelConformer make a clash with any
	// neighboring protein residue?
	bool clash_;

	// How many clashes does this NMRDummySpinlabelConformer make?
	// If all NMRDummySpinlabelConformers make a clash, we return that
	// one with the lowest clash score.
	Real clash_score_;

	// A Residue object of this spin-label residue. We use this to do an
	// FA energy-based bump check with the packer instead of the distance-based
	// clash check with the voxel grid.
	ResidueOP residue_;

	// NMRDummySpinlabelAtom objects are used for clash check with voxel grid
	// and are stored in an unordered map
	NMRDummySpinlabelAtomTable atom_table_;

};

/// @brief Class that represents a "dummy" spinlabel ensemble. It holds information
///        about each individual conformer member (atom coordinates, number of observations)
///        and data to score those based on atom clashes with the neighborhood of the ensemble.
///        To speed up the clash score calculation, this class holds also an object of type
///        NMRDummySpinlabelVoxelGrid which gets filled with constant pointers to NMRDummySpinlabelAtoms.
class NMRDummySpinlabelEnsemble : public utility::pointer::ReferenceCount {

public: // typedefs

	typedef utility::vector1< NMRDummySpinlabelConformerOP >::iterator NMRDummySpinlabelConformerTableIter;
	typedef utility::vector1< NMRDummySpinlabelConformerOP >::const_iterator NMRDummySpinlabelConformerTableCOIter;

private: // Methods

	/// @brief default constructor
	NMRDummySpinlabelEnsemble();

public: // Methods

	/// @brief Construct from database file and residue type
	NMRDummySpinlabelEnsemble(
		std::string const & database_file,
		chemical::ResidueType const & restype
	) /* throw(utility::excn::EXCN_FileNotFound, utility::excn::EXCN_Base) */ ;

	/// @brief copy constructor
	NMRDummySpinlabelEnsemble(NMRDummySpinlabelEnsemble const & other);

	/// @brief assignment operator
	NMRDummySpinlabelEnsemble &
	operator=(NMRDummySpinlabelEnsemble const & rhs);

	/// @brief destructor
	~NMRDummySpinlabelEnsemble();

	/// @brief Perform a clash score calculation for every spinlabel conformer in the ensemble by
	///        calculating pairwise distances to neighborhood residues which are within a given radius.
	///        Specifically, the distance between every side-chain heavy atom in the
	///        spinlabel conformer and the neighbor atom of a pose residue is calculated.
	///        If the distance is smaller than the pose residue's neighbor radius, is is assumed
	///        that the spinlabel conformer will fall within the side-chain radius of the
	///        neighborhood residue and thus it will be labeled as clashing otherwise it will not.
	void
	clash_check(
		pose::Pose const & pose,
		Size const target_resid,
		Real const radius
	);

	/// @brief returns the transformation matrix which transforms a xyz coordinate from
	///        the frame of the spinlabel ensemble into the frame of the target site residue.
	HT
	coordinate_transform_onto_target_site(
		pose::Pose const & pose,
		Size const resid
	) const;

	/// Getters
	utility::vector1< NMRDummySpinlabelConformerOP > & get_conformer_table() { return conformer_table_; }
	Size get_ensemble_size() const { return ensemble_size_; }
	utility::vector1<utility::vector1<Real>> const & get_rmsd_mat() { return rmsd_mat_; }
	NMRDummySpinlabelVoxelGridOP get_voxel_grid() { return grid_; } // non-constant method, in case we want to modify (e.g. translate) and reuse the grid
	NMRDummySpinlabelVoxelGridCOP get_voxel_grid() const { return grid_; }

private: // Method

	/// @brief utility function used for initialization from database files
	void init_from_database_file(
		std::string const & database_file,
		chemical::ResidueType const & restype
	);

	/// @brief define the origin of the xyz coordinate frame of the ensemble.
	///        this needs to be done only once when we initialize the object.
	void define_ensemble_frame();

	/// @brief calculate the homogeneous transform that relates the target residue
	///        (i.e. in our case the residue that represents the spinlabel site)
	///        with the coordinate frame of the dummy spinlabel ensemble. In this
	///        way we don't need to move the spinlabel ensemble which contains many
	///        more coordinates but can leave it static.
	HT coordinate_transform_from_target_site(
		pose::Pose const & pose,
		Size const resid
	) const;

	void register_options();
	void init_from_cml();

private: // Data

	// Number of NMRDummySpinlabelConformers in NMRDummySpinlabelEnsemble
	Size ensemble_size_;

	// Vector of NMRDummySpinlabelConformers
	utility::vector1< NMRDummySpinlabelConformerOP > conformer_table_;

	// For convenience, I store the local frame of the NMRDummySpinlabelEnsemble here.
	// This is used for aligning it to the target site in the protein.
	HT ensemble_origin_;

	// RMSD matrix between NMRDummySpinlabelConformers used for clustering
	utility::vector1<utility::vector1<Real>> rmsd_mat_;

	// Voxel grid for clash check
	NMRDummySpinlabelVoxelGridOP grid_;

	// Count number of clashes for every NMRDummySpinlabelConformer.
	// If all conformers have a clash, return that one with the smallest
	// number of clashes.
	bool elaborate_clash_check_;

};

} // nmr
} // scoring
} // core

#endif // INCLUDED_core_scoring_nmr_NMRDummySpinlabelEnsemble_HH
