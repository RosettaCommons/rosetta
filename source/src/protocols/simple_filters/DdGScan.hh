// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/DdGScan.hh
/// @brief definition of filter class DdGScan.
/// @author Neil King (neilking@uw.edu)
/// @author Kyle Barlow (kb@kylebarlow.com)

#ifndef INCLUDED_protocols_simple_filters_DdGScan_hh
#define INCLUDED_protocols_simple_filters_DdGScan_hh

#include <protocols/simple_filters/DdGScan.fwd.hh>

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
#include <protocols/simple_moves/ddG.fwd.hh>

#include <utility/vector1.hh>
#include <boost/tuple/tuple.hpp>

#include <set>

typedef boost::tuples::tuple< core::Size, std::string, core::Real > ddG_data_tuple;

namespace protocols {
namespace simple_filters {

class DdGScan : public protocols::filters::Filter
{
public :
	// Constructors, destructor, virtual constructors
	DdGScan();
	DdGScan(
		core::pack::task::TaskFactoryOP task_factory,
		core::Size const repeats,
		core::scoring::ScoreFunctionCOP scorefxn,
		bool report_diffs,
		bool write2pdb
	);
	DdGScan( DdGScan const & rval );
	virtual protocols::filters::FilterOP clone() const;
	virtual protocols::filters::FilterOP fresh_instance() const;
	virtual ~DdGScan();
	void initialize();

	// @brief get name of this filter
	virtual std::string name() const { return "DdGScan"; }

	// Setters
	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	void repeats( core::Size const r );
	void scorefxn( core::scoring::ScoreFunctionOP const scorefxn );
	void report_diffs( bool const report_diffs );
	void write2pdb( bool const write );
	void ddG_mover( protocols::simple_moves::ddGOP const ddG_mover_op );

	// Getters
	core::pack::task::TaskFactoryOP task_factory() const;
	core::Size repeats() const;
	bool report_diffs() const;
	bool write2pdb() const;
	protocols::simple_moves::ddGOP ddG_mover() const;

	// Dummy apply function
	virtual bool apply( core::pose::Pose const & ) const;

	// Parse xml
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const &
	);

	void parse_def( utility::lua::LuaObject const & def,
		utility::lua::LuaObject const & score_fxns,
		utility::lua::LuaObject const & tasks );

	// Effector functions
	core::Real ddG_for_single_residue(
		core::pose::Pose const & const_pose,
		core::Size const resi,
		core::pack::task::PackerTaskOP const general_task,
		core::pose::Pose & pose_to_mutate
	) const;

	void write_to_pdb(
		core::pose::Pose const & pose,
		core::Size const & residue,
		core::Size const & output_resi,
		std::string const & residue_name,
		core::Real const & ddG
	) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	utility::vector1< ddG_data_tuple > calculate( std::ostream & out, core::pose::Pose const & pose ) const;

private:

	core::pack::task::TaskFactoryOP task_factory_;
	core::Size repeats_;
	core::scoring::ScoreFunctionOP scorefxn_;
	bool report_diffs_;
	bool write2pdb_;
	protocols::simple_moves::ddGOP ddG_mover_;

};

} // simple_filters
} // protocols
#endif
