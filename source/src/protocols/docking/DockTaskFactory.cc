// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DockTaskFactory
/// @brief sets up the task factory for docking movers
/// @details
///  This contains the functions that set up the taskfactory for docking
///  movers depending on the command line options
/// @author Krishna Kilambi

#include <basic/options/option.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/docking/DockTaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh> // trans-clude <core/pack/task/operation/TaskOperation.hh>
#include <protocols/simple_task_operations/InterfaceTaskOperation.hh>
#include <protocols/simple_task_operations/RestrictToInterface.hh>
#include <protocols/task_operations/RestrictChainToRepackingOperation.hh>
#include <core/conformation/Residue.hh> // for design() flag
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

//for resfile reading
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

// option key includes
#include <basic/options/keys/docking.OptionKeys.gen.hh>

// C++ headers
#include <string>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/docking/DockingHighRes.hh>
#include <utility/vector0.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.docking.DockTaskFactory" );

using namespace core;

namespace protocols {
namespace docking {

// default constructor
DockTaskFactory::DockTaskFactory() : utility::pointer::ReferenceCount()
{
	set_default();
	register_options();
	init_from_options();
}

DockTaskFactory::DockTaskFactory( DockTaskFactory const & old_instance ) : utility::pointer::ReferenceCount( old_instance )
{
	resfile_ = old_instance.resfile_;
	norepack1_ = old_instance.norepack1_;
	norepack2_ = old_instance.norepack2_;
	design_chains_ = old_instance.design_chains_;
	prepack_only_ = old_instance.prepack_only_;
	restrict_to_interface_ = old_instance.restrict_to_interface_;
}
//destructor
DockTaskFactory::~DockTaskFactory() = default;

void
DockTaskFactory::set_default()
{
	prepack_only_ = false;
	norepack1_ = false;
	norepack2_ = false;
	resfile_ = false;
	//design_chains_ = utility::tools::make_vector1( NULL );
	design_chains_.clear();
	additional_task_operations_.clear();

	restrict_to_interface_ = simple_task_operations::InterfaceTaskOperationOP( new simple_task_operations::RestrictToInterface() );

	// @TODO needs to change so that these options can be set through setters and
	// do not have to be called from the commandline
	// the following should default to true always for docking
	// ex1_ = true;
	// ex2aro_ = true;
	// include_input_sc_ = true;
}

void
DockTaskFactory::init_from_options()
{
	using namespace basic::options;

	if ( option[ OptionKeys::docking::norepack1 ].user() ) {
		set_norepack1(option[ OptionKeys::docking::norepack1 ]());
	} else {
		set_norepack1(false);
	}
	if ( option[ OptionKeys::docking::norepack2 ].user() ) {
		set_norepack2(option[ OptionKeys::docking::norepack2 ]());
	} else {
		set_norepack2(false);
	}
	// why is it like this? should true be the default value?
	if ( option[ OptionKeys::packing::resfile].user() ) {
		resfile_ = true;
	}

	// @TODO should be packing level, not docking level
	if ( option[ OptionKeys::docking::design_chains ].user() ) {
		utility::vector1< std::string > chains = option[ OptionKeys::docking::design_chains ]();
		for ( core::Size i = 1; i <= chains.size(); ++i ) {
			design_chains_.push_back( chains[i][0] );
		}
	}
}

void
DockTaskFactory::register_options()
{
	using namespace basic::options;

	option.add_relevant( OptionKeys::docking::norepack1 );
	option.add_relevant( OptionKeys::docking::norepack2 );
	option.add_relevant( OptionKeys::packing::resfile );
}
void DockTaskFactory::set_additional_task_operarations( utility::vector1< core::pack::task::operation::TaskOperationOP > additional_task_operations )
{
	additional_task_operations_ = additional_task_operations;
}

void DockTaskFactory::add_additional_task_operaration( core::pack::task::operation::TaskOperationOP task_operation )
{
	additional_task_operations_.push_back( task_operation );
}

utility::vector1< core::pack::task::operation::TaskOperationOP > DockTaskFactory::get_additional_task_operarations()
{
	return additional_task_operations_;
}

void DockTaskFactory::set_interface_definition_task_operation(  protocols::simple_task_operations::InterfaceTaskOperationOP /*interface_definition*/ )
{
	//restrict_to_interface_ = interface_definition;
	return;
}
void
DockTaskFactory::create_and_attach_task_factory(
	DockingHighRes * docker,
	core::pose::Pose const & pose
) const
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace protocols::simple_task_operations;
	using namespace protocols::task_operations;

	TaskFactoryOP tf( new TaskFactory() );

	// if ( init_tf_ ) tf = new TaskFactory( *init_tf_ );

	// check to see which specific chains are being redesigned
	// set restrict_to_repacking on all other chains
	if ( design_chains_.size() > 0 ) {
		// check for movable jumps
		for ( core::Size i = 1; i <= pose.num_jump(); ++i ) {
			core::Size const cutpoint = pose.fold_tree().cutpoint_by_jump( i );
			char chain = pose.pdb_info()->chain( pose.pdb_info()->number( cutpoint ) );
			if ( find( design_chains_.begin(), design_chains_.end(), chain ) != design_chains_.end() ) {
				tf->push_back( TaskOperationCOP( new protocols::task_operations::RestrictChainToRepackingOperation( chain ) ) );
			}
		}
	} else {
		tf->push_back( TaskOperationCOP( new RestrictToRepacking ) );
	}

	tf->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
	tf->push_back( TaskOperationCOP( new IncludeCurrent ) );
	tf->push_back( TaskOperationCOP( new NoRepackDisulfides ) );
	if ( resfile_ ) tf->push_back( TaskOperationCOP( new ReadResfile ) );

	// DockingNoRepack only works over the first rb_jump in movable_jumps
	// In a 2-body case this separates 1 & 2 based on the only cutpoint
	// In a multibody case, this separates 1 & 2 based on the first cutpoint
	if ( norepack1_ ) tf->push_back( TaskOperationCOP( new DockingNoRepack1( docker->movable_jumps()[1] ) ) );
	if ( norepack2_ ) tf->push_back( TaskOperationCOP( new DockingNoRepack2( docker->movable_jumps()[1] ) ) );

	// incorporating Ian's UnboundRotamer operation.
	// note that nothing happens if unboundrot option is inactive!
	core::pack::rotamer_set::UnboundRotamersOperationOP unboundrot( new core::pack::rotamer_set::UnboundRotamersOperation() );
	unboundrot->initialize_from_command_line();
	operation::AppendRotamerSetOP unboundrot_operation( new operation::AppendRotamerSet( unboundrot ) );
	tf->push_back( unboundrot_operation );

	//@TODO pose needed for this, set it up somewhere else
	// core::pack::dunbrack::load_unboundrot(pose); // adds scoring bonuses for the "unbound" rotamers, if any

	// note that RestrictToInterfaceOperation is added during set_dock_mcm_protocol
	if ( !prepack_only_ ) {
		restrict_to_interface_->set_movable_jumps( docker->movable_jumps() );
		tf->push_back( restrict_to_interface_ );  //JQX: add restrict to interface, the Legacy code used this in the initial packing as well
	}

	// Add user specified task operations
	for ( auto const & additional_task_operation : additional_task_operations_ ) {
		tf->push_back( additional_task_operation );
	}

	docker->set_task_factory( tf );
}


} // namespace docking
} // namespace protocols
