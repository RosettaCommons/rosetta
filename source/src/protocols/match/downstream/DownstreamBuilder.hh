// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/downstream/DownstreamBuilder.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini
/// @author Florian Richter (floric@u.washington.edu)

#ifndef INCLUDED_protocols_match_downstream_DownstreamBuilder_hh
#define INCLUDED_protocols_match_downstream_DownstreamBuilder_hh

// Unit headers
#include <protocols/match/downstream/DownstreamBuilder.fwd.hh>

// Package headers
#include <protocols/match/BumpGrid.fwd.hh>
#include <protocols/match/OccupiedSpaceHash.fwd.hh>
#include <protocols/match/downstream/ActiveSiteGrid.fwd.hh>
#include <protocols/match/downstream/RigidLigandBuilder.fwd.hh>
#include <protocols/match/downstream/LigandConformerBuilder.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/LigandConformer.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Numeric headers
#include <numeric/HomogeneousTransform.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>

#include <core/id/AtomID.fwd.hh>
#include <protocols/match/Hit.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace downstream {

class DownstreamBuilder : public utility::pointer::ReferenceCount {
public:
	typedef utility::pointer::ReferenceCount      parent;
	typedef core::id::AtomID                      AtomID;
	typedef core::Vector                          Vector;
	typedef core::Real                            Real;
	typedef numeric::HomogeneousTransform< Real > HTReal;

public:
	DownstreamBuilder();

	DownstreamBuilder( DownstreamBuilder const & );

	virtual ~DownstreamBuilder();

	virtual
	DownstreamBuilderOP
	clone() const = 0;

	virtual
	std::list< Hit >
	build(
		HTReal const & atom6_frame,
		Size const scaffold_build_point_id,
		Size const upstream_conf_id,
		Size const external_geometry_id,
		core::conformation::Residue const & upstream_residue
	) const = 0;

	virtual
	void
	set_bb_grid(
		BumpGridCOP bbgrid
	);

	void
	set_occupied_space_hash(
		OccupiedSpaceHashCOP occ_space
	);

	void
	set_active_site_grid(
		ActiveSiteGridCOP active_site_grid
	);

	/// @brief In case downstream builders can return hits that
	/// are incompatible with each other (e.g. different ligand
	/// conformations ) the matcher needs to know about this to
	/// allow for speedy match enumeration
	virtual bool
	hits_potentially_incompatible() const = 0;

	virtual bool compatible(
		Hit const & my_hit,
		DownstreamBuilder const & other,
		Hit const & other_hit,
		bool first_dispatch = true
	) const;

	virtual bool compatible(
		Hit const & my_hit,
		RigidLigandBuilder const & other,
		Hit const & other_hit,
		bool first_dispatch = true
	) const;

	virtual bool compatible(
		Hit const & my_hit,
		LigandConformerBuilder const & other,
		Hit const & other_hit,
		bool first_dispatch = true
	) const;

	virtual
	void
	require_atom_to_reside_in_active_site(
		core::id::AtomID const & id
	) = 0;

	virtual
	ProbeRadius
	atom1_radius() const = 0;

	virtual
	ProbeRadius
	atom2_radius() const = 0;

	virtual
	ProbeRadius
	atom3_radius() const = 0;

	virtual
	bool
	atom1_belongs_in_active_site() const = 0;

	virtual
	bool
	atom2_belongs_in_active_site() const = 0;

	virtual
	bool
	atom3_belongs_in_active_site() const = 0;


	virtual
	Real
	atom1_atom2_distance() const = 0;

	virtual
	Real
	atom2_atom3_distance() const = 0;

	/// @brief Returns an angle in degrees between the three downstream atoms.
	virtual
	Real
	atom1_atom2_atom3_angle() const = 0;

	virtual
	void
	coordinates_from_hit(
		Hit const & hit,
		utility::vector1< AtomID > const & atom_indices, // these are the indices of the atoms whose coordinates are requested
		utility::vector1< Vector > & atom_coords       // ouput coordinates are placed in this array.
	) const = 0;

	virtual
	toolbox::match_enzdes_util::LigandConformerOP
	get_lig_conformers(core::Size conf_id) const = 0;

	//  virtual
	//  utility::vector1< utility::vector1< std::pair< core::Size, core::Real > > >
	//  get_min_sep_d2_from_upstream_atoms() const = 0;

	virtual
	core::chemical::ResidueTypeCOP
	get_upstream_restype() const = 0;

	virtual
	core::pose::PoseCOP
	downstream_pose_from_hit(
		Hit const & hit
	) const = 0;

	virtual
	Size
	n_possible_hits_per_at3frame() const = 0;


protected:
	bool
	bbgrid_set() const {
		return bb_grid_.get() != 0;
	}

	BumpGrid const &
	bbgrid() const {
		return *bb_grid_;
	}

	bool
	occ_space_set() const {
		return space_ != 0;
	}

	OccupiedSpaceHash const &
	occ_space() const {
		return *space_;
	}

	bool
	active_site_grid_set() const {
		return active_site_grid_ != 0;
	}

	ActiveSiteGrid const &
	active_site_grid() const {
		return *active_site_grid_;
	}

private:
	BumpGridCOP bb_grid_;
	OccupiedSpaceHashCOP space_;
	ActiveSiteGridCOP active_site_grid_;

};

}
}
}

#endif
