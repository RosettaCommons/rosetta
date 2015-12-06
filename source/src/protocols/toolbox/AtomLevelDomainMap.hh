// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#ifndef INCLUDED_protocols_toolbox_AtomLevelDomainMap_HH
#define INCLUDED_protocols_toolbox_AtomLevelDomainMap_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh>
#include <protocols/toolbox/AtomID_Mapper.fwd.hh>
#include <core/pose/copydofs/CopyDofs.hh> // for FIXED_DOMAIN (999)
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <map>

#ifdef WIN32
#include <core/id/AtomID.hh>
#endif

namespace protocols {
namespace toolbox {

using core::pose::copydofs::FIXED_DOMAIN;

//////////////////////////////////////////////////////////////////////////////////////////////
//
// Contains a map of each atom to a "domain":
//   0      = no domain, can move this atom if you want
//   1,2... = domain occupied by a predefined "chunk". Don't shift relative atom
//              positions within each domain.
//   999    = (FIXED) part of a user-defined fixed domain.

class AtomLevelDomainMap : public utility::pointer::ReferenceCount {
public:

	//constructor
	AtomLevelDomainMap( core::pose::Pose const & pose,
											bool const map_to_vanilla_pose = false,
											utility::vector1< Size > const & allow_insert_res = utility::vector1< Size >() /* do not have to specify this*/ );

	virtual ~AtomLevelDomainMap();

	AtomLevelDomainMap( AtomLevelDomainMap const & src );

	virtual AtomLevelDomainMapOP clone() const;

public:

	bool
	has_domain( core::id::AtomID const & atom_id ) const;

	bool
	get( Size const & i ) const;

	bool
	get( core::id::AtomID const & atom_id  ) const;

	bool
	get( core::id::TorsionID const & torsion_id, core::conformation::Conformation const & conformation ) const;

	bool
	get_jump( Size const & jump_number, core::conformation::Conformation const & conformation ) const;

	Size
	get_domain( Size const & i ) const;

	Size
	get_domain( core::id::AtomID const & atom_id  ) const;

	Size
	get_domain( core::id::NamedAtomID const & named_atom_id, core::pose::Pose const & pose  ) const;

	void
	set_domain( Size const & i, Size const & setting  );

	void
	set_domain( core::id::AtomID const & atom_id, Size const & setting  );

	void
	set_domain( core::id::NamedAtomID const & atom_id, core::pose::Pose const & pose,  Size const & setting  );

	void
	set_domain( Size const & setting  );

	void
	set_phosphate_domain( core::Size const & i,
		core::pose::Pose const & pose,
		Size const & setting );
	void
	set_phosphate( core::Size const & i,
		core::pose::Pose const & pose,
		bool const & setting );

	void
	set_sugar_domain( core::Size const & i,
		core::pose::Pose const & pose,
		Size const & setting );
	void
	set_sugar( core::Size const & i,
		core::pose::Pose const & pose,
		bool const & setting );

	void
	set( Size const & i, bool const & setting  );

	void
	set( core::id::NamedAtomID const & named_atom_id, core::pose::Pose const & pose, bool const & setting  );

	void
	set( core::id::AtomID const & atom_id, bool const & setting  );

	void
	set( bool const & setting  );

	void
	set_fixed_if_moving( Size const & i );

	void
	set_fixed_if_moving( core::id::AtomID const & atom_id );

	void
	show() const;

	std::map< core::id::AtomID, Size >
	calculate_atom_id_domain_map( core::pose::Pose const & pose ) const;

	void
	renumber_after_variant_changes( core::pose::Pose const & pose );

	void
	setup_movemap( core::kinematics::MoveMap & mm,
								 core::pose::Pose const & pose );

	AtomID_MapperCOP
	atom_id_mapper() const { return atom_id_mapper_; }

private:

	void
	initialize( core::pose::Pose const & pose,
							bool const map_to_vanilla_pose,
							utility::vector1< Size > const & allow_insert_res );

	void
	update_to_move_internal_phosphates( pose::Pose const & pose );

	void
	update_to_not_move_virtual_phosphates( pose::Pose const & pose );

	void
	update_to_not_move_last_virtual_residue( pose::Pose const & pose );

	void
	apply_allow_insert_res( utility::vector1< Size > const & allow_insert_res );

private:

	std::map< core::id::AtomID, Size > domain_map_;
	AtomID_MapperCOP atom_id_mapper_;

};


} //rna
} //protocols

#endif
