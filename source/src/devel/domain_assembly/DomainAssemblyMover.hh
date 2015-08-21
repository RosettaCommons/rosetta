// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/domain_assembly/DomainAssemblyMover
/// @brief


#ifndef INCLUDED_devel_domain_assembly_DomainAssemblyMover_hh
#define INCLUDED_devel_domain_assembly_DomainAssemblyMover_hh

// Unit Headers
#include <devel/domain_assembly/DomainAssemblyMover.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/fragment/ConstantLengthFragSet.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.hh>

#include <protocols/moves/Mover.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// STL Headers
#include <iosfwd>
#include <map>

#include <utility/vector1.hh>


namespace devel {
namespace domain_assembly {

class DomainAssemblyMover : public protocols::moves::Mover
{
public:

	DomainAssemblyMover();
	virtual ~DomainAssemblyMover();
	virtual void apply( core::pose::Pose & pose );
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual std::string get_name() const;

protected: // protocol stages
	virtual void run_fullatom_stage( core::pose::Pose & pose );
	void run_centroid_stage( core::pose::Pose & pose );
	void run_abinitio_centroid_stage( core::pose::Pose & pose );
	virtual void run_fullatom_relax( core::pose::Pose & pose );
	void evaluate_pose( core::pose::Pose const & pose ) const;
protected: // helper functions
	void get_domain_definition( core::pose::Pose const & pose, std::vector< std::string > & domains ) const;
	std::string get_linker_definition( core::pose::Pose const & pose ) const;
	core::Real target_rmsd( core::pose::Pose const & pose ) const;
	core::Real target_rmsd_no_align( core::pose::Pose const & pose ) const;
	void recover_sidechains( core::pose::Pose & pose ) const;
protected: // accessors and such
	core::kinematics::MoveMap const & move_map() const { return *movemap_; }
	core::pose::Pose const & target_pose() const { return target_pose_; }
	core::pose::Pose const & starting_pose() const { return starting_pose_; }
	std::string const & buried_residues() const { return buried_residues_; }
	std::map< core::Size, core::Size > const & target_pose_map() const { return target_pose_map_; }
	bool movemap_set() const { return movemap_set_; }
	bool fragsets_set() const { return (fragset3mer_ && fragset9mer_); }
protected: // initializers
	void initialize();
	void initialize_movemap_from_commandline();
	void initialize_fragments_from_commandline();
	void initialize_start_pose_from_commandline();
	void initialize_target_pose();
	void initialize_pose_map_from_commandline();
	void initialize_buried_from_commandline();
private: // data
	core::kinematics::MoveMapOP movemap_;
	bool movemap_set_;
	core::pose::Pose target_pose_;
	core::pose::Pose starting_pose_;
	std::map< core::Size, core::Size > target_pose_map_;
	std::string buried_residues_; // user-specified residues to be buried by the assembly (?caging)
	core::fragment::ConstantLengthFragSetOP fragset3mer_;
	core::fragment::ConstantLengthFragSetOP fragset9mer_;
};

}
}


#endif
