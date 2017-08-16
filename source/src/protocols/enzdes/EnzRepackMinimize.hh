// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/enzdes/EnzRepackMinimize.hh
/// @author Sagar Khare (khares@uw.edu)

#ifndef INCLUDED_protocols_enzdes_EnzRepackMinimize_hh
#define INCLUDED_protocols_enzdes_EnzRepackMinimize_hh

// Unit headers
#include <protocols/enzdes/EnzRepackMinimize.fwd.hh>

// package headers

// Project Headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

#include <protocols/filters/Filter.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <utility/tag/Tag.fwd.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace enzdes {

class EnzRepackMinimize : public protocols::moves::Mover
{
public:
	EnzRepackMinimize();
	EnzRepackMinimize( std::string const & name );
	~EnzRepackMinimize() override;

	/// @brief turns on constraint score terms in case the pose
	/// has an enzdes cst cache
	void ensure_scofx_cstfication( core::pose::Pose const & pose);

	void minimize_in_stages(
		core::pose::Pose & pose,
		core::pack::task::PackerTaskCOP task,
		bool const & cstopt,
		bool const & min_sc,
		bool const & min_rb,
		bool const & min_lig);

	//parser stuff
	void apply( core::pose::Pose &pose ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	utility::vector1< core::Size > rb_min_jumps() const;
	void rb_min_jumps( utility::vector1< core::Size > const );


	//managing minimization options
	void set_min_rb( bool const min_rb ) { min_rb_ = min_rb; }
	bool get_min_rb() const { return min_rb_; }
	void set_min_sc( bool const min_sc ) { min_sc_ = min_sc; }
	bool get_min_sc() const {  return min_sc_; }
	void set_min_lig( bool const min_lig ) { min_lig_ = min_lig; }
	bool get_min_lig() const {  return min_lig_; }
	void set_min_bb(  bool const min_bb ) { min_bb_ = min_bb; }
	bool get_min_bb() const { return min_bb_; }

	//design?
	bool get_design() const { return design_; }
	void set_design( bool const des ) { design_ = des; }

	//scorefunctions
	void set_scorefxn_repack( core::scoring::ScoreFunctionCOP scorefxn );
	void set_scorefxn_minimize( core::scoring::ScoreFunctionCOP scorefxn );
	core::scoring::ScoreFunctionOP scorefxn_repack() const;
	core::scoring::ScoreFunctionOP scorefxn_minimize() const;

	//task stuff
	core::pack::task::PackerTaskOP create_ptask( core::pose::Pose & pose );
	void task_factory( core::pack::task::TaskFactoryOP p );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::scoring::ScoreFunctionOP scorefxn_repack_;
	core::scoring::ScoreFunctionOP scorefxn_minimize_;
	bool cst_opt_; //do polyala active-site optimization?
	bool design_, repack_, fix_catalytic_; // design or only repack?
	bool minimize_in_stages_; // minimize non-bb dofs followed by bb+rb-dofs?
	bool min_rb_, min_sc_, min_bb_, min_lig_, minimize_prot_jumps_, backrub_;
	core::pack::task::PackerTaskOP task_;
	core::pack::task::TaskFactoryOP task_factory_; // sequence positions and residue-level tasks to apply when setup_packer_task is done
	core::Size n_cycles_;
	utility::vector1< core::Size > rb_min_jumps_; //dflt empty vector; this supersedes rb_min_ above if it is defined, allowing the user to handpick which jumps to minimize. SJF 15Oct13
};

} // enzdes
} // protocols

#endif /*INCLUDED_protocols_enzdes_EnzRepackMinimize_HH*/
