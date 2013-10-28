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
#include <core/pose/Pose.fwd.hh>

#include <protocols/moves/Mover.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// STL Headers
#include <iosfwd>
#include <map>

#include <utility/vector1.hh>



namespace devel{
namespace domain_assembly{

class DomainAssemblyMover : public protocols::moves::Mover
{
public:

	DomainAssemblyMover();

	virtual void apply( core::pose::Pose & pose );
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP	fresh_instance() const;
	virtual std::string get_name() const;

	void move_map( core::kinematics::MoveMap const & setting );

private:
	void run_fullatom_stage( core::pose::Pose & pose );
	void run_centroid_stage( core::pose::Pose & pose );

private:
	void initialize();
	void initialize_movemap_from_commandline();
	void initialize_fragments_from_commandline();

private:

	core::kinematics::MoveMapOP movemap_;
	bool movemap_set_;

	core::fragment::ConstantLengthFragSetOP fragset3mer_;

};

}
}


#endif
