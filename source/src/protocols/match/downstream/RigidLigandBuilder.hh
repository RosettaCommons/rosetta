// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/downstream/RigidLigandBuilder.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini
/// @author Florian Richter (floric@u.washington.edu)

#ifndef INCLUDED_protocols_match_downstream_RigidLigandBuilder_hh
#define INCLUDED_protocols_match_downstream_RigidLigandBuilder_hh

// Unit headers
#include <protocols/match/downstream/RigidLigandBuilder.fwd.hh>

// Package headers
#include <protocols/match/downstream/DownstreamBuilder.hh>
#include <protocols/toolbox/match_enzdes_util/LigandConformer.fwd.hh>
#include <protocols/match/BumpGrid.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#ifdef __clang__
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#endif

// Numeric headers

// Utility headers
#include <utility/VirtualBase.hh>

// C++ headers
#include <list>

#include <utility/fixedsizearray1.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace downstream {

class RigidLigandBuilder : public DownstreamBuilder {
public:
	typedef DownstreamBuilder   parent;

public:
	RigidLigandBuilder();

	RigidLigandBuilder( RigidLigandBuilder const & );

	//initiates collision check in the Secondary matcher
	// RigidLigandBuilder( RigidLigandBuilder const & , core::chemical::ResidueTypeCOP upstream_restype );

	~RigidLigandBuilder() override;

	DownstreamBuilderOP
	clone() const override;

	std::list< Hit >
	build(
		HTReal const & atom3_frame,
		core::Size const scaffold_build_point_id,
		core::Size const upstream_conf_id,
		core::Size const external_geometry_id,
		core::conformation::Residue const & upstream_residue
	) const override;

	void
	set_bb_grid(
		BumpGridCOP bbgrid
	) override;

	bool
	hits_potentially_incompatible() const override{
		return false; }

	bool compatible(
		Hit const & my_hit,
		DownstreamBuilder const & other,
		Hit const & other_hit,
		bool first_dispatch = true
	) const override;

	bool compatible(
		Hit const & my_hit,
		RigidLigandBuilder const & other,
		Hit const & other_hit,
		bool first_dispatch = true
	) const override;

	void
	require_atom_to_reside_in_active_site(
		core::id::AtomID const & id
	) override;

	ProbeRadius
	atom1_radius() const override;

	ProbeRadius
	atom2_radius() const override;

	ProbeRadius
	atom3_radius() const override;

	bool
	atom1_belongs_in_active_site() const override;

	bool
	atom2_belongs_in_active_site() const override;

	bool
	atom3_belongs_in_active_site() const override;


	Real
	atom1_atom2_distance() const override;

	Real
	atom2_atom3_distance() const override;

	/// @brief Returns an angle in degrees between the three downstream atoms.
	Real
	atom1_atom2_atom3_angle() const override;

	void
	coordinates_from_hit(
		Hit const & hit,
		utility::vector1< AtomID > const & atom_indices,
		utility::vector1< Vector > & atom_coords
	) const override;

	core::pose::PoseCOP
	downstream_pose_from_hit(
		Hit const & hit
	) const override;

	core::Size
	n_possible_hits_per_at3frame() const override;

	toolbox::match_enzdes_util::LigandConformerOP
	get_lig_conformers(core::Size conf_id) const override;

	// virtual
	//  utility::vector1< utility::vector1< std::pair< core::Size, core::Real > > >
	//  get_min_sep_d2_from_upstream_atoms() const;

	core::chemical::ResidueTypeCOP
	get_upstream_restype() const override;

public:
	// Initialization

	/// @brief Specify the residue, with coordinates, that's being used as the downstream
	/// partner.  This class is meant to be used in conjuction with the ClassicMatchAglrotihm,
	/// and therefore the initialization routines are specific for that algorithm.  In this
	/// initialization function, one must list atoms "D1, D2 and D3" in the convention of
	/// describing the rigid-body orientation between three atoms of the upstream partner
	/// (atoms U3, U2 & U1) and three atoms of the downstream partner (atoms D1, D2 & D3) in terms
	/// of 2 angles, 1 distance, and 3 dihedrals.  The user must also list the 3 atoms used to
	/// define the orientation frame of the downstream ligand.  It is essential to the
	/// matching algorithm that the same three orientation atoms are used for all RigidLigandBuilders.
	void
	initialize_from_residue(
		core::Size atom1,
		core::Size atom2,
		core::Size atom3,
		core::Size orientation_atom1,
		core::Size orientation_atom2,
		core::Size orientation_atom3,
		core::conformation::Residue const & residue
	);

	void
	initialize_upstream_residue(
		core::chemical::ResidueTypeCOP upstream_res,
		core::scoring::etable::count_pair::CountPairFunctionCOP count_pair = nullptr
	);

	void ignore_h_collisions( bool setting );

	/*Real6
	global_orientation_from_frame3(
	HTReal const & frame3
	) const;

	HTReal
	frame_from_global_orientation(
	Real6 orientation
	) const;*/

private:

	void
	initialize_upstream_nonbonded_min_separation_d2();

private:

	core::chemical::ResidueTypeCOP downstream_restype_;
	core::chemical::ResidueTypeCOP upstream_restype_;

	bool ignore_h_collisions_;

	/// The indices of the three atoms defining the orientation of the
	/// ligand in the global coordinate frame
	/// These indices are in the restype indexing of atoms.
	utility::fixedsizearray1< core::Size, 3 >         orientation_atoms_;
	utility::fixedsizearray1< core::Size, 3 >         atoms_123_;
	utility::fixedsizearray1< ProbeRadius, 3 >  radii_123_;
	utility::fixedsizearray1< bool, 3 >   ats123_reqd_in_active_site_;

	utility::vector1< ProbeRadius >  atom_radii_;
	utility::vector1< bool >         atom_required_in_active_site_;

	utility::vector1< core::Size >  non_collision_detection_atoms_reqd_in_active_site_;

	//LigandConformerOP lig_conformer_;
	utility::vector1< toolbox::match_enzdes_util::LigandConformerOP > lig_conformers_;

	/// Detect collision between the upstream residue (sidechain?!) conformation and the atoms of the
	/// downstream residue
	utility::vector1< utility::vector1< std::pair< core::Size, Real > > > min_sep_d2_from_upstream_atoms_;

};

}
}
}

#endif
