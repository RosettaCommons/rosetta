// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/TaskAwareScoreTypeFilter.hh
/// @brief  header file for TaskAwareScoreTypeFilter
/// @author Jacob Bale (balej@uw.edu), Neil King (neilking@uw.edu)

#ifndef INCLUDED_protocols_simple_filters_TaskAwareScoreTypeFilter_hh
#define INCLUDED_protocols_simple_filters_TaskAwareScoreTypeFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/simple_filters/TaskAwareScoreTypeFilter.fwd.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/vector1.hh>
#include <string>

// Unit headers

namespace protocols {
namespace simple_filters {

class TaskAwareScoreTypeFilter : public protocols::filters::Filter

{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	TaskAwareScoreTypeFilter();
	/// @brief Constructor with a single target residue
	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual protocols::filters::FilterOP clone() const;
	virtual protocols::filters::FilterOP fresh_instance() const;
	core::Real compute( core::pose::Pose const & pose, bool const & write ) const;
	virtual ~TaskAwareScoreTypeFilter();
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	core::pack::task::TaskFactoryOP task_factory() const;
	core::scoring::ScoreFunctionOP scorefxn() const;
	core::scoring::ScoreType score_type() const;
	core::Real threshold() const;
	bool bb_bb() const;
	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	void score_type( core::scoring::ScoreType const st );
	std::string score_type_name() const;
	void score_type_name( std::string const name );
	void threshold( core::Real const thresh );
	void bb_bb( bool const bb );
	bool unbound() const;
	void unbound( bool const unbound );
	std::string sym_dof_names() const;
	void sym_dof_names( std::string const sym_dofs );
	core::Size jump() const;
	void jump( core::Size const jump );
	std::string mode() const;
	void mode( std::string const mode );
	void write2pdb( bool const write );
	bool write2pdb() const;
	void write_to_pdb( core::pose::Pose const & pose, core::Size const residue, std::string const residue_name, core::Real const score ) const;
	bool individual_hbonds() const;
	void individual_hbonds( bool const individual_hbonds );
private:
	core::pack::task::TaskFactoryOP task_factory_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreType score_type_; // scoretype of interest, defaults to total_energy
	std::string score_type_name_;
	core::Real threshold_; // If mode is individual, then this threshold is for the maximum number of failing residues. If mode is average, then this is the maximum average score over all considered positions. If mode is total, then this is the maximum total score for the given scoretype over all considered positions.
	bool bb_bb_; // if one wants to evaluate backbone - backbone hydrogen bonding energies (short and long range)
	bool unbound_; // For evaluating energies in the unbound state
	bool write2pdb_; //Output per residue information to the PDB file for residues that fail the user-defined threshold
	bool individual_hbonds_; //Caution. Under development. Calculate energies for individual hbonds between the sidechains of the residues specified by the task operations and backbone polars.
	std::string sym_dof_names_; // For generating the unbound pose
	core::Size jump_; // Alternative method of generating the unbound pose
	std::string mode_; // This filter can be set to operate in three different modes: total, average, or individual
};

} // simple_filters
} // protocols

#endif

