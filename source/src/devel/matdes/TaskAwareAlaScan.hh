// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/matdes/TaskAwareAlaScan.hh
/// @brief definition of filter class TaskAwareAlaScan.
/// @author Neil King (neilking@uw.edu)

#ifndef INCLUDED_devel_matdes_TaskAwareAlaScan_hh
#define INCLUDED_devel_matdes_TaskAwareAlaScan_hh

#include <devel/matdes/TaskAwareAlaScan.fwd.hh>


// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
//#include <utility/exit.hh>

//#include <utility/vector1.hh>

#include <set>

namespace devel {
namespace matdes {
	
class TaskAwareAlaScan : public protocols::filters::Filter
{
public :
	// Constructors, destructor, virtual constructors
	TaskAwareAlaScan();
	TaskAwareAlaScan(
			core::pack::task::TaskFactoryOP task_factory,
			core::pack::task::TaskFactoryOP ddG_task_factory,
			bool use_custom_task,
			core::Size const jump,
			std::string const sym_dof_name,
			core::Size const repeats,
			core::scoring::ScoreFunctionCOP scorefxn,
			bool repack,
			bool report_diffs );
	TaskAwareAlaScan( TaskAwareAlaScan const & rval );
	virtual protocols::filters::FilterOP clone() const;
	virtual protocols::filters::FilterOP fresh_instance() const;
	virtual ~TaskAwareAlaScan();

  // @brief get name of this filter
  virtual std::string name() const { return "TaskAwareAlaScan"; }

	// Setters
	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	void ddG_task_factory( core::pack::task::TaskFactoryOP ddG_task_factory );
	void use_custom_task( bool const uct );
	void jump( core::Size const j );
	void sym_dof_name( std::string const j );
	void repeats( core::Size const r );
	void scorefxn( core::scoring::ScoreFunctionOP const scorefxn );
	void repack( bool const repack );
	void report_diffs( bool const report_diffs );
	void write2pdb( bool const write );

	// Getters
	core::pack::task::TaskFactoryOP task_factory() const;
	core::pack::task::TaskFactoryOP ddG_task_factory() const;
	bool use_custom_task() const;
	core::Size jump() const;
	std::string sym_dof_name() const;
	core::Size repeats() const;
	bool repack() const;
	bool report_diffs() const;
	bool write2pdb() const;

	// Dummy apply function
	virtual bool apply( core::pose::Pose const & ) const;

	// Parse xml
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	void parse_ddG_task_operations( utility::tag::TagCOP tag, basic::datacache::DataMap const & data );
	void parse_def( utility::lua::LuaObject const & def,
				utility::lua::LuaObject const & score_fxns,
				utility::lua::LuaObject const & tasks );

	// Effector functions
	core::Real ddG_for_single_residue( core::pose::Pose const & const_pose, core::Size const resi ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	void write_to_pdb( core::Size const & residue, std::string const & residue_name, core::Real const & ddG ) const;

private:

	core::pack::task::TaskFactoryOP task_factory_;
	core::pack::task::TaskFactoryOP ddG_task_factory_;
	bool use_custom_task_;
	core::Size jump_;
	std::string sym_dof_name_;
	core::Size repeats_;
	core::scoring::ScoreFunctionOP scorefxn_;
	bool repack_;
	bool report_diffs_;
	std::set< std::string > exempt_identities_;
	bool write2pdb_;

};

} // matdes
} // devel
#endif
