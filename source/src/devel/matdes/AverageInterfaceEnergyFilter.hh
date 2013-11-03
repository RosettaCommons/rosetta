// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/AverageInterfaceEnergyFilter.hh
/// @brief Reports to Tracer which residues are designable in a taskfactory
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_devel_matdes_AverageInterfaceEnergyFilter_hh
#define INCLUDED_devel_matdes_AverageInterfaceEnergyFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <devel/matdes/AverageInterfaceEnergyFilter.fwd.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/vector1.hh>
#include <string>

// Unit headers

namespace devel {
namespace matdes {

class AverageInterfaceEnergyFilter : public protocols::filters::Filter

{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
		AverageInterfaceEnergyFilter();
	///@brief Constructor with a single target residue
	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual protocols::filters::FilterOP clone() const;
	virtual protocols::filters::FilterOP fresh_instance() const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~AverageInterfaceEnergyFilter();
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	void parse_def( utility::lua::LuaObject const & def,
					utility::lua::LuaObject const & score_fxns,
					utility::lua::LuaObject const & tasks );
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
	bool return_total() const;
	void return_total( bool const total );
private:
	core::pack::task::TaskFactoryOP task_factory_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreType score_type_; // scoretype of interest, defaults to total_energy 
	std::string score_type_name_;
	core::Real threshold_; // upper threshold for average interface energy of a given scoretype 
	bool bb_bb_; // if one wants to evaluate backbone - backbone hydrogen bonding energies (short and long range)
	bool unbound_; // For evaluating energies in the unbound state
	std::string sym_dof_names_; // For generating the unbound pose
	core::Size jump_; // Alternative method of generating the unbound pose
	bool return_total_;
};

} // matdes
} // devel

#endif

