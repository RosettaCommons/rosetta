// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/BumpGrid.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_BumpGrid_hh
#define INCLUDED_protocols_match_BumpGrid_hh

// Unit headers
#include <protocols/match/BumpGrid.fwd.hh>

// Package headers

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/fixedsizearray1.hh>
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/geometry/BoundingBox.hh>
#include <numeric/geometry/hashing/SixDHasher.fwd.hh> // typedefs for Bin3D.

// C++ headers
#include <string>
#include <iosfwd>

#include <utility/vector0_bool.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {

class Bool3DGrid : public utility::pointer::ReferenceCount
{
public:
	typedef core::Real   Real;
	typedef core::Vector Vector;
	typedef numeric::geometry::BoundingBox< core::Vector > BoundingBox;
	typedef numeric::geometry::hashing::Bin3D Bin3D;
	typedef utility::fixedsizearray1< Vector, 8 > CornerPoints;
	typedef std::pair< Bool3DGrid::Size, unsigned char > index_mask_pair;

public:
	/// Creation and initialization
	Bool3DGrid();
	virtual ~Bool3DGrid();

	void set_bounding_box( BoundingBox const & bb );
	void set_bin_width( Real width );

	/// @brief create a grid for the input sphere that aligns to this grid,
	/// such that it is large enough to hold a particular sphere.
	Bool3DGrid
	create_grid_for_sphere( Vector const & center, Real radius ) const;

	/// @brief create a grid for the input bounding box that aligns to this grid
	Bool3DGrid
	create_grid_for_bb( BoundingBox const & bb );

public:
	/// Accessors
	Bin3D dimsizes() const;

	Real bin_width() const { return bin_width_; }

	CornerPoints corners( Bin3D const & bin ) const;

	/// @brief bounding box points for a bin.
	BoundingBox  bin_extrema( Bin3D const & bin ) const;

	Vector
	bin_center( Bin3D const & bin ) const;

	bool occupied( Vector const & ) const;

	bool occupied( Bin3D const & bin ) const;

	BoundingBox actual_bb() const { return bb_extended_; }

public:
	// Mutators

	/// @brief set the boolean value for a particular bin.
	/// The bin dimensions are indexed from 0 to nbins-1.
	void set_value_for_bin( Bin3D const & bin, bool setting );

	/// @brief Set the value to true for any voxel that is wholy contained
	/// by the given sphere.  A voxel is wholy contained iff all 8 corners
	/// of the voxel are within the sphere.
	void
	or_by_sphere_conservative( Vector const & center, Real radius );


	/// @brief Set the value to true for any voxel that is partially contained
	/// by a given sphere.  A voxel is partially contained iff any of the 8 corners
	/// of the voxel are within the sphere.
	void
	or_by_sphere_liberal( Vector const & center, Real radius );

	/// @brief Consider a voxel filled if each of its corners are covered,
	/// even if they are covered by seperate spheres.  Dangerous, in that some
	/// voxels will be counted as being fully occupied when they are only partially occupied.
	void
	or_by_spheres_conservative(
		utility::vector1< std::pair< Vector, Real > > const & spheres
	);

	//// @brief Turn the values of all the bins that overlap with the
	/// volume in this bounding box to true.
	void
	or_by_box_liberal(
		BoundingBox const & bb
	);

	/// @brief Performs a boolean OR on the voxels shared by the two grids,
	/// and ignores all voxels that are not shared by the two grids.  These grids
	/// must be "compatible", in that the lower corner of other must lie on a grid
	/// point of this.
	void or_with( Bool3DGrid const & other );

	/// @brief Performs a boolean AND on the voxels shared by the two grids,
	/// and ignores all voxels that are not shared by the two grids.  These grids
	/// must be "compatible", in that the lower corner of other must lie on a grid
	/// point of this.
	void and_with( Bool3DGrid const & other );

	/// @brief Sets any voxel on this grid to "false" which is true on the other grid.
	/// but does not set any voxel to "true" on this grid -- even if the other voxel is "false".
	void subtract( Bool3DGrid const & other );

	/// @brief Set all values in all bins to false
	void clear();

private:
	index_mask_pair
	index_and_mask_for_point( Vector const & point ) const;

	index_mask_pair
	index_and_mask_for_bin( Bin3D const & bin ) const;

	Bin3D bin_for_point( Vector const & point ) const;

	void
	reset_grid();

	unsigned char
	mask_from_offsets( Size xmod2, Size ymod2, Size zmod2 ) const;

	unsigned char
	negmask_from_offsets( Size xmod2, Size ymod2, Size zmod2 ) const;

	Size
	byte_index_from_doublebin( Bin3D const & halfbin ) const;


private:
	BoundingBox bb_;          // user defined volume of space covered by this grid.
	BoundingBox bb_extended_; // the actual volume of space covered by this grid.

	Real bin_width_;
	Real bin_width_2x_;
	Real half_bin_width_root_three_;
	Bin3D dimsizes_;
	Bin3D dimprods_;

	Bin3D halfdimsizes_; // the double-bins are twice as large as the bins in each dimension
	Bin3D halfdimprods_;

	Bin3D supervoxel_dimsizes_; // super voxels hold several double-bins, as defined by the
	Bin3D supervoxel_dimprods_;

	/// creates super voxels of size 8x8x8 = 512 bits -> 64 bytes = cache line size on Core2 Duo,
	/// though this size should not have a very large effect for caches of larger or smaller size.
	/// The main point of the super voxel is to allow data locality as it's very often the same
	/// region of space that's being queried repeatedly during rotamer building.
	static const Size n_doublebins_per_supervoxel = 4;

	/// Each byte represents the boolean "covered" status for 8 voxels
	/// in a double-voxel.  The double voxel
	utility::vector0< unsigned char > grid_;

};

class BumpGrid : public utility::pointer::ReferenceCount
{
public:
	typedef Bool3DGrid::BoundingBox          BoundingBox;
	typedef Bool3DGrid::Bin3D                Bin3D;
	typedef Bool3DGrid::Vector               Vector;
	typedef std::pair< Vector, core::Real >  Sphere;
	typedef core::Real                       Real;
	typedef utility::pointer::ReferenceCount parent;

public:
	BumpGrid();
	BumpGrid( BumpGrid const & );
	virtual ~BumpGrid();

	BumpGrid & operator = ( BumpGrid const & );


public:
	/// Initialization
	void set_bounding_box( BoundingBox const & bb );

	//// @brief Set the tolerance for sphere overlap between all sphere types.
	/// All tolerances must be set *before* any spheres are ORed into the grid.
	void set_general_overlap_tolerance( Real tolerated_overlap );

	//// @brief Set the overlap tolerance for a specific pair of sphere types (e.g. OXY/NIT)
	/// The specific tolerance is combined with the general overlap tolerance.
	/// All tolerances must be set *before* any spheres are ORed into the grid.
	void set_pair_overlap_tolerance( ProbeRadius rad1, ProbeRadius rad2, Real tolerated_overlap );

public:
	BumpGrid
	create_bump_grid_for_bb( BoundingBox const & bb ) const;

	BumpGridOP
	create_new_bump_grid_for_bb( BoundingBox const & bb ) const;

	void or_with( BumpGrid const & );
	void and_with( BumpGrid const & );

	void or_by_sphere( Sphere const & );
	void or_by_sphere( Vector const & center, ProbeRadius radius_type );

	bool overlaps( BumpGrid const & ) const;

	/// @brief Collision detection by a point p with a given (fixed) radius type.
	/// Collision detection is performed in "configuration space", where the obstacles
	/// have been convolved with a spherical probe of a given radius.  This collision
	/// detection is conservative: it will not report a collision to exist that does not,
	/// but may miss a collision that does exist.
	bool occupied( ProbeRadius radius_type, Vector const & p ) const;

	static Real probe_radius( ProbeRadius radtype ) {
		switch ( radtype ) {
		case ZERO  : return 0.0;
		case H_ARO : return 1.0;
		case H_ALA : return 1.17;
		case OXY   : return 1.40;
		case NIT   : return 1.55;
		case C_CAR : return 1.65;
		case C_ALA : return 1.75;
		case SULPH : return 1.85;
		}
		return 0.0;
	}

	Real pair_overlap_tolerance( ProbeRadius rad1, ProbeRadius rad2 ) const;

	Bool3DGrid const &
	grid( ProbeRadius radtype ) const {
		return * grids_[ radtype ];
	}

	inline
	Real
	required_separation_distance( ProbeRadius rad1, ProbeRadius rad2 ) const {
		if ( rad1 == ZERO || rad2 == ZERO ) return 0.0;
		return std::max( probe_radius( rad1 ) + probe_radius( rad2 )
			- pair_permit_overlap_[ rad1 ][ rad2 ]
			- general_permit_overlap_, 0.0 );
	}

private:

	utility::vector1< Bool3DGridOP > grids_;

	/// allow some pairs of atom types to overlap by a certain amount -- e.g. O and N
	utility::vector1< utility::vector1< Real > > pair_permit_overlap_;

	/// allow all pairs of atom types to overlap by a certain amount
	Real general_permit_overlap_;

};

ProbeRadius
probe_radius_for_atom_type( core::Size atomtype );

BumpGridOP
bump_grid_to_enclose_pose( core::pose::Pose const & pose );


/// @brief Construct a BumpGrid that encloses a single residue.  Use the
/// original_grid as a starting point, copying over all pertinent data
/// such that the two grids can later be merged together.
BumpGridOP
bump_grid_to_enclose_residue(
	core::conformation::Residue const & residue,
	BumpGrid const & original_grid
);

/// @brief Construct a BumpGrid that encloses the backbone atoms for a single residue.
/// Use the original_grid as a starting point, copying over all pertinent data
/// such that the two grids can later be merged together.
BumpGridOP
bump_grid_to_enclose_residue_backbone(
	core::conformation::Residue const & residue,
	BumpGrid const & original_grid
);


void
fill_grid_with_residue_spheres(
	core::conformation::Residue const & residue,
	BumpGrid & grid
);

void
fill_grid_with_residue_heavyatom_spheres(
	core::conformation::Residue const & residue,
	BumpGrid & grid
);

void
fill_grid_with_backbone_heavyatom_spheres(
	core::conformation::Residue const & residue,
	BumpGrid & grid
);


class Bool3DGridKinemageWriter
{
public:
	typedef core::Real   Real;
	typedef core::Size   Size;
	typedef core::Vector Vector;
	typedef numeric::geometry::hashing::Bin3D Bin3D;

public:
	Bool3DGridKinemageWriter();
	~Bool3DGridKinemageWriter();

	void set_unselectable( bool unselectable );
	void set_line_color( std::string const & line_color );
	void set_master( std::string const & master );
	void set_shrink_factor( Real shrink_factor );
	void set_skip_completely_buried_positions( bool skip_completely_buried_positions );

	void set_write_empty_voxels( bool write_empty_voxels_ );
	void set_empty_voxel_color( std::string const & empty_voxel_color_ );

	void set_write_facets( bool write_facets );
	void set_facet_master( std::string const & facet_master );
	void set_facet_color( std::string const & facet_color );
	void set_transparent_facets( bool transparent_facets );
	void set_facet_alpha( Real facet_alpha );

	void write_grid_to_kinemage(
		std::ostream & ostr,
		std::string const & group_name,
		Bool3DGrid const & grid ) const;

	void write_grid_to_file(
		std::string const & fname,
		std::string const & group_name,
		Bool3DGrid const & grid
	) const;

private:
	bool unselectable_;
	std::string line_color_;
	std::string master_;
	Real shrink_factor_;
	bool skip_completely_buried_positions_;

	bool write_empty_voxels_;
	std::string empty_voxel_color_;

	bool write_facets_;
	std::string facet_master_;
	std::string facet_color_;
	bool transparent_facets_;
	Real facet_alpha_;

};

}
}

#endif
