// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/TaskAwareCsts.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/TaskAwareCsts.hh>
#include <protocols/protein_interface_design/movers/TaskAwareCstsCreator.hh>
// Package headers
#include <core/pose/Pose.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
//Auto Headers
#include <core/conformation/Residue.hh>
#include <numeric/xyzVector.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/id/AtomID.hh>
#include <core/pose/selection.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.movers.TaskAwareCsts" );
std::string
TaskAwareCstsCreator::keyname() const
{
	return TaskAwareCstsCreator::mover_name();
}

protocols::moves::MoverOP
TaskAwareCstsCreator::create_mover() const {
	return protocols::moves::MoverOP( new TaskAwareCsts );
}

std::string
TaskAwareCstsCreator::mover_name()
{
	return "TaskAwareCsts";
}

TaskAwareCsts::TaskAwareCsts() :
	Mover( TaskAwareCstsCreator::mover_name() ),
	task_factory_( /* NULL */ ),
	cst_type_( "coordinate" ),
	anchor_resnum_( "" )
{
}


TaskAwareCsts::~TaskAwareCsts() {}

void
TaskAwareCsts::apply( core::pose::Pose & pose )
{
	using namespace protocols::rosetta_scripts;
	using namespace core::scoring::constraints;
	using core::id::AtomID;

	ConstraintCOPs cst;
	utility::vector1< core::Size > const designable( residue_packer_states( pose, task_factory(), true/*designable*/, false/*packable*/ ) );
	runtime_assert( designable.size() );
	AtomID const anchor_atom( AtomID( pose.residue( designable[ 1 ] ).atom_index( "CA" ),
		( anchor_resnum_ == "" ? designable[ 1 ] /*anchor to first designable residue*/ : core::pose::parse_resnum( anchor_resnum_, pose ) ) ) );
	core::scoring::func::HarmonicFuncOP coord_cst_func( new core::scoring::func::HarmonicFunc( 0.0, 1.0/*sd*/ ) ); // hardwired for now
	TR<<"Adding constraints to pose at positions: ";
	for ( core::Size const resid : designable ) {
		core::conformation::Residue const rsd_i( pose.residue( resid ) );
		if ( cst_type_ == "coordinate" ) {
			cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( AtomID( rsd_i.atom_index( "CA" ), resid ), anchor_atom, rsd_i.xyz( "CA" ), coord_cst_func ) ) );
			TR<<resid<<',';
		}
	}
	TR<<std::endl;
	pose.add_constraints( cst );
}

std::string
TaskAwareCsts::get_name() const {
	return TaskAwareCstsCreator::mover_name();
}

void
TaskAwareCsts::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	using namespace protocols::rosetta_scripts;
	task_factory( parse_task_operations( tag, data ) );
	cst_type( tag->getOption< std::string >( "cst_type", "coordinate" ) );
	anchor_resnum( tag->getOption< std::string > ( "anchor_resnum", "" ) );

	TR<<"cst_type: "<<cst_type()<<" anchor_resnum: "<<anchor_resnum()<<std::endl;
}

protocols::moves::MoverOP
TaskAwareCsts::clone() const {
	return( protocols::moves::MoverOP( new TaskAwareCsts( *this ) ));
}

core::pack::task::TaskFactoryOP
TaskAwareCsts::task_factory() const{ return task_factory_; }

void
TaskAwareCsts::task_factory( core::pack::task::TaskFactoryOP tf ){ task_factory_ = tf; }

} //movers
} //protein_interface_design
} //protocols
