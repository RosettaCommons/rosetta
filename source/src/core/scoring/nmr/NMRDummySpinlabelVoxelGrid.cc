// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRDummySpinlabelVoxelGrid.hh
/// @brief   Implementation of classes NMRDummySpinlabelVoxelGrid and VoxelGridPoint
/// @details last Modified: 11/02/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/NMRDummySpinlabelVoxelGrid.hh>

// Package headers
#include <core/scoring/nmr/NMRDummySpinlabelEnsemble.hh>

// Project headers

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// C++ headers

namespace core {
namespace scoring {
namespace nmr {

static basic::Tracer TR( "core.scoring.nmr.NMRDummySpinlabelVoxelGrid" );

////////////////////////////////////////
/// Class VoxelGridPoint
////////////////////////////////////////

/// @brief Default constructor
VoxelGridPoint::VoxelGridPoint() :
	utility::pointer::ReferenceCount(),
	atom_name_(""),
	id_(),
	coords_()
{}

/// @brief Construct from atom name and 3D coordinates
VoxelGridPoint::VoxelGridPoint(
	std::string const & name,
	Size const id,
	Vector const & coords
) :
	utility::pointer::ReferenceCount(),
	atom_name_(name),
	id_(id),
	coords_(coords)
{}

/// @brief Destructor
VoxelGridPoint::~VoxelGridPoint() {}

/// @brief Type name of this voxel grid point
std::string
VoxelGridPoint::type() const {
	return "";
}

/// @brief Is relevant for neighbor search
bool
VoxelGridPoint::is_relevant_neighbor() const {
	return true;
}

////////////////////////////////////////
/// Class NMRDummySpinlabelAtom
////////////////////////////////////////

/// @brief Default constructor
NMRDummySpinlabelAtom::NMRDummySpinlabelAtom() :
	VoxelGridPoint(),
	conformer_(nullptr)
{}

/// @brief Construct from 3D coordinates and owning pointer to NMRDummySpinlabelConformer
///        Convert owning pointer into access pointer to NMRDummySpinlabelConformer
NMRDummySpinlabelAtom::NMRDummySpinlabelAtom(
	std::string const & name,
	Size const id,
	Vector const & coords,
	NMRDummySpinlabelConformer const & conformer
) :
	VoxelGridPoint(name, id, coords),
	conformer_(&conformer)
{}

/// @brief Destructor
NMRDummySpinlabelAtom::~NMRDummySpinlabelAtom() {}

/// @brief Type name of this voxel grid point
std::string
NMRDummySpinlabelAtom::type() const {
	return "SL";
}

/// @brief Is relevant for neighbor search
bool
NMRDummySpinlabelAtom::is_relevant_neighbor() const {
	// atom is relevant only if the conformer to which it belongs
	// isn't marked as clashing yet.
	return !conformer_->has_clash();
}

////////////////////////////////////////
/// Class VoxelGridPoint_SL
////////////////////////////////////////

/// @brief Default constructor
VoxelGridPoint_AA::VoxelGridPoint_AA() :
	VoxelGridPoint(),
	residue_(nullptr)
{}

/// @brief Construct from 3D coordinates and owning pointer to Residue
VoxelGridPoint_AA::VoxelGridPoint_AA(
	std::string const & name,
	Size const id,
	Vector const & coords,
	conformation::Residue const & resi
) :
	VoxelGridPoint(name, id, coords),
	residue_(&resi)
{}

/// @brief Destructor
VoxelGridPoint_AA::~VoxelGridPoint_AA() {}

/// @brief Type name of this voxel grid point
std::string
VoxelGridPoint_AA::type() const {
	return "AA";
}

////////////////////////////////////////
/// Class NMRDummySpinlabelVoxelGrid
////////////////////////////////////////

/// @brief Construct from a vector of points
NMRDummySpinlabelVoxelGrid::NMRDummySpinlabelVoxelGrid(
	Real const & resolution,
	utility::vector1< VoxelGridPoint > const & points,
	bool const & cache_edges
) :
	VoxelGrid< VoxelGridPoint >(resolution, cache_edges)
{
	utility::vector1< VoxelGridPoint const * > points_ptr;
	points_ptr.reserve(points.size());
	for ( auto & p : points ) {
		points_ptr.push_back(&p);
	}
	SetObjects(points_ptr);
}

/// @brief Construct from a vector of points
NMRDummySpinlabelVoxelGrid::NMRDummySpinlabelVoxelGrid(
	Real const & resolution,
	utility::vector1< VoxelGridPoint const * > const & points,
	bool const & cache_edges
) :
	VoxelGrid< VoxelGridPoint >(resolution, cache_edges)
{
	SetObjects(points);
}

/// @brief Destructor
NMRDummySpinlabelVoxelGrid::~NMRDummySpinlabelVoxelGrid() {}

/// @brief Extract the 3D coordinate of a given object of type VoxelGridPoint.
Vector const *
NMRDummySpinlabelVoxelGrid::ExtractPosition(VoxelGridPoint const & point) const {
	return &(point.get_coordinates());
}

/// @brief Check if two grid points are the same.
///        Note, that we could also lock() the derived classes' member pointer and
///        test them for equality, to know if the two points belong to the same group,
///        but this operation would probably take too long.
bool
NMRDummySpinlabelVoxelGrid::IsSameItem(
	VoxelGridPoint const & point1,
	VoxelGridPoint const & point2
) const
{
	return point1.type() == point2.type() && point1.atom_name() == point2.atom_name()
		&& point1.id() == point2.id();
}

/// @brief Check if this item is relevant for neighbor search
bool
NMRDummySpinlabelVoxelGrid::IsRelevantItem(VoxelGridPoint const & point) const {
	return point.is_relevant_neighbor();
}

} // namespace nmr
} // namespace scoring
} // namespace core
