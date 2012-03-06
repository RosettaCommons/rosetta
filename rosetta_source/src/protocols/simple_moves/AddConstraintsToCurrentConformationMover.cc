// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#include <protocols/simple_moves/AddConstraintsToCurrentConformationMover.hh>
#include <protocols/simple_moves/AddConstraintsToCurrentConformationMoverCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>

// task operation
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <protocols/rosetta_scripts/util.hh>

// utility
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

// option
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

static basic::Tracer TR( "protocols.simple_moves.AddConstraintsToCurrentConformationMover" );

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace basic::options;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;
using namespace constraints;
	
AddConstraintsToCurrentConformationMover::AddConstraintsToCurrentConformationMover(){}

AddConstraintsToCurrentConformationMover::~AddConstraintsToCurrentConformationMover(){}
		
void AddConstraintsToCurrentConformationMover::apply( core::pose::Pose & pose )
{
	using namespace conformation;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose::datacache;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::id;
	using namespace protocols::moves;
	using namespace core::scoring;
	

	if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
		core::pose::addVirtualResAsRoot(pose);
	}
	
	TR << pose.fold_tree() << std::endl;
	
	core::pose::Pose constraint_target_pose = pose;
	
	Size nres = pose.total_residue();
	Real const coord_sdev( option[ OptionKeys::relax::coord_cst_stdev ] );
	for ( Size ires = 1; ires <= nres; ++ires ) {
		if ( !pose.residue_type(ires).has("CA") ) continue;
		Size iatom = pose.residue_type(ires).atom_index("CA");
		if (!option[ OptionKeys::relax::coord_cst_width ].user() ) {
			pose.add_constraint( new CoordinateConstraint(
								  AtomID(iatom,ires), AtomID(1,pose.fold_tree().root()), pose.residue(ires).xyz(iatom),
								  new HarmonicFunc( 0.0, coord_sdev ) ) );
		} else {
			Real const cst_width( option[ OptionKeys::relax::coord_cst_width ]() );
			pose.add_constraint( new CoordinateConstraint(
								  AtomID(iatom,ires), AtomID(1,pose.fold_tree().root()), pose.residue(ires).xyz(iatom),
								  new BoundFunc( 0, cst_width, coord_sdev, "xyz" )) );
		}
	}
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
AddConstraintsToCurrentConformationMover::parse_my_tag(
								TagPtr const tag,
								moves::DataMap & datamap,
								Filters_map const & filters,
								moves::Movers_map const & movers,
								Pose const & pose
								)
{
	if ( tag->hasOption("cst_width") ) {
		cst_width_ = tag->getOption<core::Real>("cst_width");
	}
	if ( tag->hasOption("cst_stdev") ) {
		cst_stdev_ = tag->getOption<core::Real>("cst_stdev");
	}
	parse_task_operations( tag, datamap, filters, movers, pose );
}

void
AddConstraintsToCurrentConformationMover::parse_task_operations(
										 TagPtr const tag,
										 moves::DataMap const & datamap,
										 Filters_map const &,
										 moves::Movers_map const &,
										 Pose const &
										 )
{
	TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == 0) return;
	task_factory( new_task_factory );
}

void AddConstraintsToCurrentConformationMover::task_factory( TaskFactoryCOP tf )
{
	runtime_assert( tf );
	task_factory_ = tf;
}
	
moves::MoverOP AddConstraintsToCurrentConformationMover::clone() const { return new AddConstraintsToCurrentConformationMover( *this ); }
moves::MoverOP AddConstraintsToCurrentConformationMover::fresh_instance() const { return new AddConstraintsToCurrentConformationMover; }

protocols::moves::MoverOP
AddConstraintsToCurrentConformationMoverCreator::create_mover() const {
	return new AddConstraintsToCurrentConformationMover;
}
	
std::string
AddConstraintsToCurrentConformationMoverCreator::keyname() const
{
	return AddConstraintsToCurrentConformationMoverCreator::mover_name();
}
	
std::string
AddConstraintsToCurrentConformationMoverCreator::mover_name()
{
	return "AddConstraintsToCurrentConformationMover";
}
	
std::string
AddConstraintsToCurrentConformationMover::get_name() const {
	return "AddConstraintsToCurrentConformationMover";
}
	
} // moves
} // protocols

