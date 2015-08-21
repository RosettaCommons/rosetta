// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#ifndef INCLUDED_protocols_toolbox_AllowInsert_HH
#define INCLUDED_protocols_toolbox_AllowInsert_HH

#include <protocols/toolbox/AllowInsert.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

//XRW2 suggestion: refactor to protocols/toolbox

// C++ Headers
#include <map>

#ifdef WIN32
#include <core/id/AtomID.hh>
#endif


namespace protocols {
namespace toolbox {

//////////////////////////////////////////////////////////////////////////////////////////////
//
// Contains a map of each atom to a "domain":
//   0      = no domain, can move this atom if you want
//   1,2... = domain occupied by a predefined "chunk". Don't shift relative atom
//              positions within each domain.
//   999    = (FIXED) part of a user-defined fixed domain.

extern core::Size const FIXED_DOMAIN;

class AllowInsert : public utility::pointer::ReferenceCount {
public:

	//constructor
	AllowInsert( core::pose::Pose const & pose );
	virtual ~AllowInsert();

	AllowInsert &
	operator=( AllowInsert const & src );

	AllowInsert( AllowInsert const & src );

	virtual AllowInsertOP clone() const;

public:

	bool
	has_domain( core::id::AtomID const & atom_id  ) const;

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
	and_allow_insert(AllowInsertOP allow_insert_in );

	core::Size nres() const{ return nres_;}

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
	show();

	std::map< core::id::AtomID, Size > const &
	calculated_atom_id_domain_map();

	std::map< core::id::AtomID, Size > const &
	calculate_atom_id_domain_map( core::pose::Pose const & pose );

	void
	renumber_after_variant_changes( core::pose::Pose const & pose );

	void
	calculate_atom_id_map(
		core::pose::Pose const & pose,
		std::map< core::Size, core::Size > const & res_map,
		core::kinematics::FoldTree const & scratch_fold_tree,
		std::map< core::id::AtomID, core::id::AtomID > & atom_id_map );

	void
	set_force_ideal_chainbreak( bool const & setting ){ force_ideal_chainbreak_ = setting; }

	void
	append_residue( core::pose::Pose const & pose,
		Size const & i,
		bool const & setting );

	void
	setup_movemap( core::kinematics::MoveMap & mm,
		core::pose::Pose const & pose );

private:
	void
	initialize( core::pose::Pose const & pose );

private:

	std::map< core::id::AtomID, Size > allow_insert_;
	std::map< core::id::NamedAtomID, core::id::AtomID > named_atom_id_map_;
	utility::vector1< utility::vector1< core::id::AtomID > > atom_ids_in_res_;

	std::map< core::id::AtomID, Size > calculated_atom_id_domain_map_;

	std::map< core::id::AtomID, core::id::AtomID > map_to_original_;

	core::Size nres_;

	bool force_ideal_chainbreak_;
};


} //rna
} //protocols

#endif
