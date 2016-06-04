// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file protocols/flxbb/DesignTask.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_flxbb_DesignTask_hh
#define INCLUDED_protocols_flxbb_DesignTask_hh

// Unit header
#include <protocols/flxbb/DesignTask.fwd.hh>

// Package headers
#include <protocols/flxbb/FilterStructs.fwd.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#ifdef __clang__
#include <core/pack/task/PackerTask.fwd.hh>
#endif
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <core/types.hh>


namespace protocols {
namespace flxbb {

///////////////////////////////////////////////////////////////////////////////////////////////////////
class DesignTask : public utility::pointer::ReferenceCount {
public:

	typedef std::string String;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef protocols::moves::MoverOP MoverOP;
	typedef protocols::flxbb::FilterStructsOP FilterStructsOP;
	typedef core::chemical::AA AA;

public:


	/// @brief default constructor
	DesignTask();

	/// @brief value constructor
	DesignTask(
		Size const ncycle,
		ScoreFunctionOP const sfxn,
		MoverOP const mover,
		FilterStructsOP const filter_structs=0,
		PackerTaskOP const taskf=0,
		String const & resfile=""
	);

	/// @brief copy constructor
	DesignTask( DesignTask const & rval );

	/// @brief destructor
	virtual ~DesignTask();

	/// @brief setup packer task
	virtual void setup( Pose const &, PackerTaskOP const ) = 0;


public: // accessors


	/// @brief the number of cycles of fixbb design and mover
	Size ncycle() const;

	/// @brief scorefxn for fixbb design
	ScoreFunctionOP scorefxn() const;

	/// @brief mover after fixbb design
	MoverOP mover() const;

	/// @brief filter during fixbb design
	FilterStructsOP filter_structs() const;

	/// @brief packer task for fixbb design
	PackerTaskOP packertask() const;

	/// @brief resfile
	String resfile() const;


public: // mutators


	/// @brief the number of cycles of design and mover
	void set_ncycle( Size const & ncycle );

	/// @brief scorefxn for fixbb design
	void set_scorefxn( ScoreFunctionOP const sfxn );

	/// @brief mover after fixbb design
	void set_mover( MoverOP const value );

	/// @brief filter during fixbb design
	void set_filter_structs( FilterStructsOP const value );

	/// @brief packer task
	void set_packertask( PackerTaskOP const taskf );

	/// @brief set resfile
	void set_resfile( String const & resfile );

	/// @brief add task operations
	void add_task_operations( utility::vector1< TaskOperationOP > const & top );

	/// @brief add task operation
	void add_task_operation( TaskOperationOP const top );


public:


	/// @brief output packertask
	void dump_packertask( std::ostream & os );


protected:


	/// @brief task operations
	utility::vector1< TaskOperationOP > task_operations_;


private:


	/// @brief the number of cycles of design and mover
	Size ncycle_;

	/// @brief scorefxn for fixbb design
	ScoreFunctionOP scorefxn_;

	/// @brief mover after fixbb design
	MoverOP mover_;

	/// @brief filter during fixbb design
	FilterStructsOP filter_structs_;

	/// @brief packertask used for fixbb design
	PackerTaskOP task_;

	/// @brief resfile name
	String resfile_;


};


///////////////////////////////////////////////////////////////////////////////////////////////////////
class DesignTask_Normal: public DesignTask {
public:

	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef protocols::flxbb::FilterStructsOP FilterStructsOP;

public:
	DesignTask_Normal();
	DesignTask_Normal(
		Size ncycle,
		ScoreFunctionOP sfxn,
		MoverOP mover,
		FilterStructsOP filter_structs=0 );

	virtual ~DesignTask_Normal();

	virtual void setup( Pose const & pose, PackerTaskOP const task );

};

///////////////////////////////////////////////////////////////////////////////////////////////////////
class DesignTask_Layer: public DesignTask {
public:

	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef protocols::flxbb::FilterStructsOP FilterStructsOP;

public:
	DesignTask_Layer();

	DesignTask_Layer(
		bool dsgn_core,
		bool dsgn_boundary,
		bool dsgn_surface,
		bool use_original_seq,
		Size ncycle,
		ScoreFunctionOP sfxn,
		MoverOP mover,
		FilterStructsOP filter_structs=0 );

	virtual ~DesignTask_Layer();

	virtual void setup( Pose const & pose, PackerTaskOP const task );

private:

	bool dsgn_core_, dsgn_boundary_, dsgn_surface_;
	bool use_original_seq_;

};


}
}

#endif
