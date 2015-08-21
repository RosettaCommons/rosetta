// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file protocols/flxbb/DesignTask.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


// Unit Headers
#include <protocols/flxbb/DesignTask.hh>

// Package Headers
//#include <protocols/flxbb/DesignLayerOperation.hh>
#include <protocols/flxbb/LayerDesignOperation.hh>
#include <protocols/flxbb/FilterStructs.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/options/option.hh>
#include <protocols/moves/Mover.hh>
#include <basic/Tracer.hh>

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "protocols.flxbb.DesignTask" );

namespace protocols {
namespace flxbb {

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
DesignTask::DesignTask():
	ncycle_( 1 ),
	scorefxn_( /* NULL */ ),
	mover_( /* NULL */ ),
	filter_structs_( /* NULL */ ),
	task_( /* NULL */ ),
	resfile_( "" )
{
	if ( option[ packing::resfile ].user() ) resfile_ = option[ packing::resfile ].value().at(1);
}

/// @brief value constructor
DesignTask::DesignTask(
	Size const ncycle,
	ScoreFunctionOP const sfxn,
	MoverOP const mover,
	FilterStructsOP const filter_structs,
	PackerTaskOP const taskf,
	String const & resfile
):
	ncycle_( ncycle ),
	scorefxn_( sfxn ),
	mover_( mover ),
	filter_structs_( filter_structs ),
	task_( taskf ),
	resfile_( resfile )
{
	if ( option[ packing::resfile ].user() ) resfile_ = option[ packing::resfile ].value().at(1);
}

/// @brief value constructor
DesignTask::DesignTask( DesignTask const & rval ) :
	utility::pointer::ReferenceCount(),
	ncycle_( rval.ncycle_ ),
	scorefxn_( rval.scorefxn_ ),
	mover_( rval.mover_ ),
	filter_structs_( rval.filter_structs_ ),
	task_( rval.task_ ),
	resfile_( rval.resfile_ )
{}

/// @brief destructor
DesignTask::~DesignTask() {}

/// @brief the number of cycles of fixbb design and mover
Size
DesignTask::ncycle() const
{
	return ncycle_;
}

/// @brief scorefxn for fixbb design
DesignTask::ScoreFunctionOP
DesignTask::scorefxn() const
{
	return scorefxn_;
}

/// @brief mover after fixbb design
DesignTask::MoverOP
DesignTask::mover() const
{
	return mover_;
}

/// @brief filter during fixbb design
DesignTask::FilterStructsOP
DesignTask::filter_structs() const
{
	return filter_structs_;
}

/// @brief packer task for fixbb design
DesignTask::PackerTaskOP
DesignTask::packertask() const
{
	return task_;
}

/// @brief resfile
DesignTask::String
DesignTask::resfile() const
{
	return resfile_;
}

/// @brief the number of cycles of fixbb design and mover
void
DesignTask::set_ncycle( Size const & ncycle )
{
	ncycle_  = ncycle;
}

/// @brief filter during fixbb design
void
DesignTask::set_scorefxn( ScoreFunctionOP const sfxn )
{
	scorefxn_  = sfxn;
}

/// @brief mover after fixbb design
void
DesignTask::set_mover( MoverOP const value )
{
	mover_  = value;
}

/// @brief filter during fixbb design
void
DesignTask::set_filter_structs( FilterStructsOP const value )
{
	filter_structs_ = value;
}

/// @brief packer task for fixbb design
void
DesignTask::set_packertask( PackerTaskOP const taskf )
{
	task_ = taskf;
}

/// @brief set resfile
void
DesignTask::set_resfile( String const & resfile )
{
	resfile_ = resfile;
}

void
DesignTask::dump_packertask( std::ostream & os )
{
	os << *task_;
}


void
DesignTask::add_task_operations( utility::vector1< TaskOperationOP > const tops )
{
	for ( utility::vector1< TaskOperationOP >::const_iterator it=tops.begin(); it!=tops.end(); ++it ) {
		add_task_operation( *it );
	}
}

void
DesignTask::add_task_operation( TaskOperationOP const top )
{
	task_operations_.push_back( top );
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
DesignTask_Layer::DesignTask_Layer() {}

/// @brief value constructor
DesignTask_Layer::DesignTask_Layer(
	bool dsgn_core,
	bool dsgn_boundary,
	bool dsgn_surface,
	bool use_original_seq,
	Size ncycle,
	ScoreFunctionOP sfxn,
	MoverOP mover,
	FilterStructsOP filter_structs
):
	DesignTask( ncycle, sfxn, mover, filter_structs ),
	dsgn_core_(dsgn_core),
	dsgn_boundary_(dsgn_boundary),
	dsgn_surface_(dsgn_surface),
	use_original_seq_( use_original_seq )
{}

/// @brief destructor
DesignTask_Layer::~DesignTask_Layer() {}

/// @brief setup PackerTask
void DesignTask_Layer::setup( pose::Pose const & pose, pack::task::PackerTaskOP const task )
{

	using namespace core::pack;
	operation::InitializeFromCommandlineOP cmop( new operation::InitializeFromCommandline );
	cmop->apply( pose, *task );
	//DesignLayerOperationOP op = new DesignLayerOperation( dsgn_core_, dsgn_boundary_, dsgn_surface_ );
	protocols::flxbb::LayerDesignOperationOP op( new protocols::flxbb::LayerDesignOperation( dsgn_core_, dsgn_boundary_, dsgn_surface_ ) );

	if ( resfile() != "" ) {
		TR << "Resfile is applied, except for the positions of AUTO " << std::endl;
		operation::ReadResfileOP rrop( new operation::ReadResfile( resfile() ) );
		rrop->apply( pose, *task );
	}

	TR << "Designed layer: core = " << dsgn_core_ << ", boundary = " << dsgn_boundary_
		<< ", sufrace = " << dsgn_surface_ << std::endl;

	if ( use_original_seq_ ) {
		op->use_original_seq();
		TR << "Original sequences are preserved for the layer you don't design."  << std::endl;
	} else {
		TR << "The region you don't design turned into Ala."  << std::endl;
	}

	op->apply( pose, *task );

	// apply additional task_opertaions
	for ( utility::vector1< TaskOperationOP >::const_iterator it=task_operations_.begin(); it!=task_operations_.end(); ++it ) {
		(*it)->apply( pose, *task );
	}

	this->set_packertask( task );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
DesignTask_Normal::DesignTask_Normal() {}

/// @brief value constructor
DesignTask_Normal::DesignTask_Normal(
	Size ncycle,
	scoring::ScoreFunctionOP sfxn,
	MoverOP mover,
	FilterStructsOP filter_structs
):
	DesignTask( ncycle, sfxn, mover, filter_structs )
{}

/// @brief destructor
DesignTask_Normal::~DesignTask_Normal() {}

/// @brief set up packer task
void DesignTask_Normal::setup( pose::Pose const & pose, pack::task::PackerTaskOP const task ){

	using namespace core::pack;
	operation::InitializeFromCommandlineOP cmop( new operation::InitializeFromCommandline );

	// set packertask based on resfile
	if ( resfile() != "" ) {
		TR << "Resfile " << resfile() << " is applied." << std::endl;
		operation::ReadResfileOP rrop( new operation::ReadResfile( resfile() ) );
		rrop->apply( pose, *task );
	}
	// initialize from command line
	cmop->apply( pose, *task );

	// apply additonal taskoperations
	for ( utility::vector1< TaskOperationOP >::const_iterator it=task_operations_.begin(); it!=task_operations_.end(); ++it ) {
		(*it)->apply( pose, *task );
	}

	this->set_packertask( task );
}

} // namespace flxbb
} //namespace protocols

