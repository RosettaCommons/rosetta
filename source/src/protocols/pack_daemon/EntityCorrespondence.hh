// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pack_daemon/EntityCorrespondence.hh
/// @brief  declaration for class EntityCorrespondence
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_pack_daemon_EntityCorrespondence_hh
#define INCLUDED_protocols_pack_daemon_EntityCorrespondence_hh

// Unit headers
#include <protocols/pack_daemon/EntityCorrespondence.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/PDBPoseMap.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/VirtualBase.hh>

// C++ headers
#include <list>
#include <iosfwd>

namespace protocols {
namespace pack_daemon {

class EntityCorrespondence : public utility::VirtualBase {
public:
	typedef core::Size                Size;
	typedef std::list< core::Size >         ResIDList;
	typedef ResIDList::const_iterator ResIDListConstIter;
	typedef utility::VirtualBase parent;
public:
	EntityCorrespondence();
	~EntityCorrespondence() override;
	EntityCorrespondence(EntityCorrespondence const &);
	EntityCorrespondence & operator = ( EntityCorrespondence const & );

	void set_pose( core::pose::PoseCOP pose );
	void set_num_entities( core::Size num_entities );
	void initialize_from_correspondence_file( std::istream & );
	void add_resid_to_entity_list( core::Size EntityID, core::Size ResID );

	core::Size num_entities() const;
	core::Size num_residues() const;
	core::Size entity_for_residue( core::Size resid ) const;
	core::Size n_residues_for_entity( core::Size entity_id ) const;
	ResIDListConstIter residues_for_entity_begin( core::Size entity_id ) const;
	ResIDListConstIter residues_for_entity_end( core::Size entity_id ) const;


private:
	void bounds_check_entity( std::string const & funcname, core::Size entity_id ) const;
	void bounds_check_residue( std::string const & funcname, core::Size resid ) const;

private:
	core::pose::PoseCOP           pose_;
	core::pose::PDBPoseMapCOP     pdb_pose_map_;
	utility::vector1< ResIDList > entity_id_2_resids_; // one entity may correspond to several residues
	utility::vector1< core::Size >      resid_2_entity_id_;  // each residue corresponds to at most one entity.

	utility::vector1< std::string > funcnames_;

};


}
}

#endif
