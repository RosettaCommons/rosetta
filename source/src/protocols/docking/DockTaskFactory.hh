// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DockTaskFactory
/// @brief sets up the task factory for docking movers
/// @detailed
///		This contains the functions that set up the taskfactory for docking
///		movers depending on the command line options
/// @author Krishna Kilambi


#ifndef INCLUDED_protocols_docking_DockTaskFactory_hh
#define INCLUDED_protocols_docking_DockTaskFactory_hh

#include <protocols/docking/DockTaskFactory.fwd.hh>

// Package headers
#include <protocols/docking/types.hh>

#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <protocols/toolbox/task_operations/InterfaceTaskOperation.fwd.hh>
// AUTO-REMOVED #include <protocols/toolbox/task_operations/RestrictToInterface.fwd.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>
#include <protocols/docking/DockingHighRes.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

//Auto Headers
// option key includes
#ifdef WIN32
	#include <core/pack/task/operation/TaskOperation.hh>
#endif


namespace protocols {
namespace docking {

class DockTaskFactory : public utility::pointer::ReferenceCount
{
public:

	/// @brief Default constructor
	DockTaskFactory();
	// Undefinded, commenting out to fix PyRosetta build  DockTaskFactory(bool restrict_2_interface); //JQX

	//copy constructor
	DockTaskFactory( DockTaskFactory const & old_instance );

	// destructor
	virtual ~DockTaskFactory();

	/// @brief Creates an appropriate TaskFactory based on command line options and attach it to a DockingHighResOP
	///		Two arguments: DockingHighResOP and a pose.
	void create_and_attach_task_factory(
		DockingHighRes * docker,
		core::pose::Pose const & pose
	) const;
	/// @brief Sets booleans to default values
	void set_default();

	/// @brief Associates relevant options with the DockTaskFactory class
	void register_options();

	/// @brief Assigns user specified values to members using command line options
	void init_from_options();

	void set_norepack1( bool norepack1 ) { norepack1_=norepack1; }
	void set_norepack2( bool norepack2 ) { norepack2_=norepack2; }
	void set_design_chains( utility::vector1< char > design_chains ) { design_chains_ = design_chains; }
    void set_additional_task_operarations( utility::vector1< core::pack::task::operation::TaskOperationOP > additional_task_operations );
    void add_additional_task_operaration( core::pack::task::operation::TaskOperationOP task_operation );
    utility::vector1< core::pack::task::operation::TaskOperationOP > get_additional_task_operarations();

    void set_interface_definition_task_operation( protocols::toolbox::task_operations::InterfaceTaskOperationOP interface_definition );
   	bool get_norepack1() const { return norepack1_; }
	bool get_norepack2() const { return norepack2_; }
	void set_prepack_only( bool prepack_only ) { prepack_only_ = prepack_only;} //JQX: add this function, one can decide to do restrict2interface

private:
	// commandline args
	bool resfile_;
	bool norepack1_;
	bool norepack2_;
	utility::vector1< char > design_chains_;
    utility::vector1< core::pack::task::operation::TaskOperationOP > additional_task_operations_;
	bool prepack_only_;  //JQX

    toolbox::task_operations::InterfaceTaskOperationOP restrict_to_interface_;

//	core::pack::task::TaskFactoryOP init_tf_;
};

} // docking
} // protocols

#endif
