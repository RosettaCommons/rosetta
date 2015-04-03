// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/solid_surface/SurfaceEnergies.hh
/// @brief  SurfaceEnergies class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Mike Pacella (mpacella@gmail.com)


#ifndef INCLUDED_core_scoring_solid_surface_SurfaceEnergies_hh
#define INCLUDED_core_scoring_solid_surface_SurfaceEnergies_hh

// Package Headers
#include <core/scoring/solid_surface/SurfaceEnergies.fwd.hh>

// Project Headers
#include <core/scoring/Energies.hh>
#include <core/conformation/find_neighbors.fwd.hh>

// Numeric Headers
//#include <numeric/geometry/BoundingBox.hh>

namespace core {
namespace scoring {
namespace solid_surface {

	/// A derived class from class Energies for efficiently representing the
	/// interactions between a protein and a fixed surface

class SurfaceEnergies : public Energies
{
public:
	typedef Energies parent;

public:
	/// ctor -- ensure correct initial state
	SurfaceEnergies();

	// Explicit copy ctor to avoid #include .hh's
	SurfaceEnergies( SurfaceEnergies const & src );

	// Explicit assignmnet operator to avoid #include .hh's
	virtual Energies const & operator = ( Energies const & rhs );

	/// @details determine whether my type is the same as another Conformation's
	virtual
	bool
	same_type_as_me( Energies const & other, bool recurse = true ) const;

	/// dtor
	virtual
	~SurfaceEnergies();

	virtual
	EnergiesOP
	clone() const;

	/// @brief The SurfaceEnergies object needs to know how many residues are in its pose;
	/// it also has to be told which residues are considered part of the surface and which
	/// residues are not part of the surface.
	void
	set_total_residue( Size total_residue );

	/// @brief Wipe away all the information in the SurfacEnergies object describing which
	/// residues are considered part of the surface and which residues are not part of the
	/// surface.  (Afterwards, all residues are going to be considered part of the surface).
	void
	reset_surface_residue_information();
	
	/// @brief Tell the SurfacEnergies that the following residues are considered not
	/// part of the surface.
	void
	set_residue_range_not_surface( Size seqpos_begin, Size seqpos_end );

	/// @brief Does the SurfaceEnergies object consider a particular residue to be part of the surface?
	bool
	residue_is_surface( Size seqpos ) const;
	
	
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	// private methods
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////

private:

	/// @brief Create a point graph representing the xyz coordinates of the "neighbor
	/// atoms."  (Each residue_type indicates one of its atoms to be used for
	/// neighbor detection -- e.g. CB for most amino acids, CA for GLY.)  The point
	/// graph will then be used in a call to find_neighbors to add upper-edges to
	/// neighboring residues.
	virtual
	void fill_point_graph( pose::Pose const & pose, conformation::PointGraphOP pg ) const;

	/// @brief Create a new surface grid, bounding box and dimension for the surface residues.
	/// This might get called if the range on the score function changes or if the surface 
	/// moves.  (NOTE: currently not prepared to handle the case where the surface moves!)
	//void
	//prepare_surface_grid( conformation::PointGraphOP pg, Real neighbor_cutoff ) const;

	/// @brief Wipe everything held in the surface grid, ensuring that in the next score function
	/// evaluation, a new grid will be computed.
	//void
	//reset_surface_grid() const;

private:
	Size total_residue_; // The total number of residues in the pose. NOTE: if this disagrees with parent::size(), we're in trouble 
	utility::vector1< std::pair< Size, Size > > non_surface_ranges_; // The residue beginning and ending indices which are considered not part of the surface
	utility::vector1< bool > is_surface_; // The residues considered part of the surface


	// Data for representing the surface that's used in neighbor detection
	mutable Real neighbor_cutoff_;
	//mutable std::map< core::conformation::CubeKey, utility::vector1< Size > > surface_grid_;
	//mutable numeric::geometry::BoundingBox< numeric::xyzVector< Real > > surf_bb_;
	//mutable core::conformation::CubeKey surf_dim_;
	
};


} // namespace solid_surface
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_solid_surface_Energies_HH
