// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRDummySpinlabelVoxelGrid.hh
/// @brief   Child class of VoxelGrid. Stores references to atom coordinates of every spinlabel
///          ensemble member. For fast lookup of relevant spinlabel atoms within a specified
///          neighborhood radius that should be considered for clash score calculation.
/// @details last Modified: 02/11/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_NMRDummySpinlabelVoxelGrid_HH
#define INCLUDED_core_scoring_nmr_NMRDummySpinlabelVoxelGrid_HH

// Unit headers
#include <core/scoring/nmr/NMRDummySpinlabelVoxelGrid.fwd.hh>

// Package headers
#include <numeric/VoxelGrid.impl.hh>
#include <core/scoring/nmr/NMRDummySpinlabelEnsemble.fwd.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

namespace core {
namespace scoring {
namespace nmr {

/// @brief Base class of a point object from which the voxel grid is created.
class VoxelGridPoint : public utility::pointer::ReferenceCount {

public: // Methods

	/// @brief Default constructor
	VoxelGridPoint();

	/// @brief Construct from atom name, ID and 3D coordinates
	VoxelGridPoint(
		std::string const & atom,
		Size const id,
		Vector const & coords
	);

	/// @brief Destructor
	virtual ~VoxelGridPoint();

	/// @brief Type name of this voxel grid point
	virtual std::string type() const;

	/// @brief Is relevant for neighbor search
	virtual bool is_relevant_neighbor() const;

	/// @brief Integer key to identify membership of this point to a group (e.g. a conformer or residue).
	///        Points belonging to the same group will have the same key. A key is a fast way in retrieving
	///        information about the group. Locking up the derived classes' member pointer and eventual
	///        casts would take too long.
	Size id() const { return id_; }

	// Getters
	Vector const & get_coordinates() const { return coords_; }
	std::string atom_name() const { return atom_name_; }

private: // Data

	std::string atom_name_;
	Size id_;
	Vector coords_;
};

/// @brief Derived class of VoxelGridPoint to store the 3D coordinates of one atom
///        of a NMRDummySpinlabelConformer. NMRDummySpinlabelAtom is also a member of
///        NMRDummySpinlabelConformer and holds a const pointer to its outer class.
class NMRDummySpinlabelAtom : public VoxelGridPoint {

public: // Methods

	/// @brief Default constructor
	NMRDummySpinlabelAtom();

	/// @brief Construct from atom name, 3D coordinates and owning pointer to NMRDummySpinlabelConformer
	///        Convert owning pointer into access pointer to NMRDummySpinlabelConformer
	NMRDummySpinlabelAtom(
		std::string const & name,
		Size const id,
		Vector const & coords,
		NMRDummySpinlabelConformer const & conformer
	);

	/// @brief Destructor
	~NMRDummySpinlabelAtom() override;

	/// @brief Type name of this voxel grid point
	std::string type() const override;

	/// @brief Is relevant for neighbor search
	bool is_relevant_neighbor() const override;

	// Getters
	NMRDummySpinlabelConformer const * conformer() const { return conformer_; }

private: // Data

	NMRDummySpinlabelConformer const * conformer_;
};

/// @brief Derived class of VoxelGridPoint to store the 3D coordinates of one atom
///        of an amino acid residue. In addition, it holds a const pointer to the residue
///        object to which it refers to.
class VoxelGridPoint_AA : public VoxelGridPoint {

public: // Methods

	/// @brief Default constructor
	VoxelGridPoint_AA();

	/// @brief Construct from 3D coordinates and owning pointer to Residue
	VoxelGridPoint_AA(
		std::string const & name,
		Size const id,
		Vector const & coords,
		conformation::Residue const & resi
	);

	/// @brief Destructor
	~VoxelGridPoint_AA() override;

	/// @brief Type name of this voxel grid point
	std::string type() const override;

	// Getters
	conformation::Residue const * residue() const { return residue_; }

private: // Data

	conformation::Residue const * residue_;
};

/// @brief Derived class of VoxelGrid to create a grid of VoxelGridPoint objects.
class NMRDummySpinlabelVoxelGrid : public numeric::VoxelGrid< VoxelGridPoint > {

public: // Methods

	/// @brief Construct from a vector of points
	NMRDummySpinlabelVoxelGrid(
		Real const & resolution,
		utility::vector1< VoxelGridPoint > const & points,
		bool const & cache_edges = false
	);

	/// @brief Construct from a vector of points
	NMRDummySpinlabelVoxelGrid(
		Real const & resolution,
		utility::vector1< VoxelGridPoint const * > const & points,
		bool const & cache_edges = false
	);

	/// @brief Destructor
	~NMRDummySpinlabelVoxelGrid() override;

	/// @brief Extract the 3D coordinate of a given object of type VoxelGridPoint.
	Vector const *
	ExtractPosition(VoxelGridPoint const & point) const override;

	/// @brief Check if two grid points are the same.
	bool
	IsSameItem(
		VoxelGridPoint const & point1,
		VoxelGridPoint const & point2
	) const override;

	/// @brief Check if this item is relevant for neighbor search
	bool
	IsRelevantItem(VoxelGridPoint const & point) const override;
};

} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_NMRDummySpinlabelVoxelGrid_HH
