// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/AtomID_Mapper.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_toolbox_AtomID_Mapper_HH
#define INCLUDED_protocols_toolbox_AtomID_Mapper_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/toolbox/AtomID_Mapper.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <map>

namespace protocols {
namespace toolbox {

class AtomID_Mapper: public utility::pointer::ReferenceCount {

public:

	//constructor
	AtomID_Mapper( core::pose::Pose const & pose,
		bool const map_to_vanilla_pose = false );

	//destructor
	~AtomID_Mapper() override;

	//constructor
	AtomID_Mapper( AtomID_Mapper const & src );

	AtomID_MapperOP
	clone() const;

public:

	bool
	has_atom_id( core::id::AtomID const & atom_id ) const;

	void
	renumber_after_variant_changes( core::pose::Pose const & pose );

	utility::vector1< core::id::AtomID > const &
	atom_ids_in_res( Size const & n ) const {
		return atom_ids_in_res_[ n ];
	}

	// @brief In use by RNA_ChunkLibrary.
	std::map< core::id::AtomID, core::id::AtomID >
	calculate_atom_id_map( core::pose::Pose const & target_pose,
		std::map< core::Size, core::Size > const & res_map /* from target to source */,
		core::kinematics::FoldTree const & source_fold_tree,
		AtomID_MapperCOP source_mapper_to_vanilla = nullptr ) const;

	core::id::AtomID const &
	map_to_reference( core::id::AtomID const & atom_id ) const;

	core::id::AtomID const &
	map_from_reference( core::id::AtomID const & atom_id ) const;

	core::Size nres() const { return atom_ids_in_res_.size();}

private:

	void
	initialize( core::pose::Pose const & pose, bool const map_to_vanilla_pose );

	void
	initialize_from_pose( core::pose::Pose const & pose );

private:

	// @brief The named_atom_id_map never changes after initialization.
	std::map< core::id::NamedAtomID, core::id::AtomID > named_atom_id_map_;
	std::string sequence_; // used to confirm viability of any renumbering requests.

	// @brief These maps change around with calls to renumber_after_variant_changes()
	std::map< core::id::AtomID, core::id::AtomID > map_to_reference_;
	std::map< core::id::AtomID, core::id::AtomID > map_from_reference_;
	utility::vector1< utility::vector1< core::id::AtomID > > atom_ids_in_res_;

};

} //toolbox
} //protocols

#endif
