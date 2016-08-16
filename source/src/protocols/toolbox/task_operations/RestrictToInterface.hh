// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RestrictToInterface.hh
/// @brief When passed to a PackerTask, pack/design is limited to the interface
/// @author ashworth

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictToInterface_hh
#define INCLUDED_protocols_toolbox_task_operations_RestrictToInterface_hh

#include <protocols/toolbox/task_operations/RestrictToInterface.fwd.hh>
#include <protocols/toolbox/task_operations/InterfaceTaskOperation.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// for parsing
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

//#include <core/conformation/Interface.hh>
#include <core/types.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/tools/make_vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <ObjexxFCL/FArray1D.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

//class RestrictTaskForDocking : public core::pack::task::operation::TaskOperation
//{
//public:
// typedef core::pack::task::operation::TaskOperation TaskOperation;
// typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
// typedef TaskOperation parent;
//public:
// RestrictTaskForDocking();
//
// RestrictTaskForDocking( core::scoring::ScoreFunctionCOP scorefxn, core::Size rb_jump, bool include_current, core::Real distance_ = 8 );
//
// virtual ~RestrictTaskForDocking();
//
// virtual TaskOperationOP clone() const;
//
// virtual void apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const;
//
//private:
// core::scoring::ScoreFunctionCOP scorefxn_;
// core::Size rb_jump_;
// bool include_current_;
// core::Real distance_;
//};

class DockingNoRepack1 : public core::pack::task::operation::TaskOperation
{
public:
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef TaskOperation parent;
public:
	DockingNoRepack1();
	DockingNoRepack1( int rb_jump_in );

	virtual ~DockingNoRepack1();

	virtual TaskOperationOP clone() const;

	virtual void apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "DockingNoRepack1"; }

private:
	int rb_jump_;

};

class DockingNoRepack2 : public core::pack::task::operation::TaskOperation
{
public:
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef TaskOperation parent;
public:
	DockingNoRepack2();
	DockingNoRepack2( int rb_jump_in );

	virtual ~DockingNoRepack2();

	virtual TaskOperationOP clone() const;

	virtual void apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "DockingNoRepack2"; }

private:
	int rb_jump_;

};

class RestrictToInterface : public InterfaceTaskOperation
{
public:
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef InterfaceTaskOperation parent;
public:
	RestrictToInterface() : parent(), distance_( 8 ), loopy_interface_( false )
	{
		set_movable_jumps( utility::tools::make_vector1< int >( 1 ) );
	}

	RestrictToInterface( int rb_jump_in, core::Real distance_in = 8 ) :
		parent(), distance_ ( distance_in ), loopy_interface_( false ) {
		set_movable_jumps( utility::tools::make_vector1< int >( rb_jump_in ) );
	}

	/// @brief Constructor with arguments for multiple jumps
	RestrictToInterface( utility::vector1_int rb_jump_in, core::Real distance_in
		= 8 ) : parent(), distance_ ( distance_in ), loopy_interface_( false ) {
		set_movable_jumps( rb_jump_in );
	}

	RestrictToInterface( utility::vector1_int rb_jump_in,
		ObjexxFCL::FArray1D_bool loop_residues ) : parent(), distance_( 8 ),
		loopy_interface_( true ) {
		loop_residues_ = loop_residues;
		set_movable_jumps( rb_jump_in );
	}

	RestrictToInterface( ObjexxFCL::FArray1D_bool loop_residues ) :
		parent(), distance_( 8 ), loopy_interface_( true ) {
		loop_residues_ = loop_residues;
	}

	RestrictToInterface( utility::vector1<bool> loop_residues );


	virtual ~RestrictToInterface();

	void add_jump( int rb_jump_in ) {
		add_movable_jump( rb_jump_in );
	}

	virtual TaskOperationOP clone() const;
	void rb_jump( int jump_in );
	//void set_movable_jumps( utility::vector1_int const movable_jumps );
	void distance( core::Real const distance_in );
	void symmetric_task( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const;

	virtual void apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const;
	virtual void parse_tag( TagCOP, DataMap & );
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "RestrictToInterface"; }

private:
	//utility::vector1_int rb_jump_;
	core::Real distance_;
	bool loopy_interface_;
	ObjexxFCL::FArray1D_bool loop_residues_;
};

}
}
}

#endif
