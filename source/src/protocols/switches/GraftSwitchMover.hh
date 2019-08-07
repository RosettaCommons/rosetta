// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Bobby Langan (robert.langan@gmail.com)

#ifndef INCLUDED_protocols_switches_GraftSwitchMover_hh
#define INCLUDED_protocols_switches_GraftSwitchMover_hh

// Unit headers
#include <protocols/switches/GraftSwitchMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>

// Utility headers
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <set>


namespace protocols {
namespace switches {

/// @brief A protocols::moves::Mover that grafts a sequence onto the terminal
/// helix of a helical-bundle monomer.  Utilizes SimpleThreadingMover alongside
///
///
class GraftSwitchMover : public protocols::moves::Mover {
public:
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef core::pack::task::TaskFactoryCOP TaskFactoryCOP;
	typedef core::pack::task::TaskFactoryOP TaskFactoryOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;

public:
	/// @brief default constructor
	GraftSwitchMover();

	// destructor (important for properly forward-declaring smart-pointer members)
	~GraftSwitchMover() override;

	// copy constructor
	GraftSwitchMover( GraftSwitchMover const & other );

	// methods

	/// @brief Grafts, minimizes, and scores functional sequences grafted
	/// onto LOCKR type de novo protien switches
	void apply( Pose & pose ) override;

	std::string get_name() const override;

	/// @brief parse XML (specifically in the context of the parser/scripting scheme)
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & ) override;

	/// @brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
	virtual void parse_score_function(
		TagCOP,
		basic::datacache::DataMap const &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	/// @brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
	virtual void parse_task_operations(
		TagCOP,
		basic::datacache::DataMap const &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP clone() const override;

	// setters

	/// @brief Sets the TaskFactory to  <tf>
	///
	/// example(s):
	///     graftmover.task_factory(task_design)
	void task_factory( TaskFactoryCOP tf );

	/// @brief Sets the ScoreFunction to  <sf>
	///
	/// example(s):
	///     graftmover.score_function(scorefxn)
	void score_function( ScoreFunctionCOP sf );

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string mover_name();

	core::pose::PoseOP get_additional_output() override;

	void add_sequence( std::string seq ) {
		sequences_.push_back(seq);
	}

	void set_start_graft( core::Size start ) {
		start_graft_ = start;
	}

	void set_end_graft( core::Size end ) {
		end_graft_ = end;
	}

	void set_important_residues( utility::vector1< core::Size > residue_list ) {
		important_residues_.push_back(residue_list);
	}

	void set_pack_neighbors(bool input) {
		pack_neighbors_ = input;
	}

	void set_pack_min(bool input){
		pack_min_ = input;
	}

	std::list< utility::vector1<core::Size> > get_output_list() {
		return output_list_;
	}

private:
	core::Size get_c_terminal_helix( core::pose::PoseOP const& pose );
	core::Size get_n_terminal_helix( core::pose::PoseOP const& pose );
	bool init_burial_filter( core::pose::Pose const& pose, core::Size index, core::Size thread_position);
	std::list< utility::vector1< core::Size > > all_threading_start_combinations();

private:
	// pointers to data that are passed in
	ScoreFunctionCOP scorefxn_;
	TaskFactoryCOP task_factory_;

	// 'really private:' packer data, actually created and owned by this class
	core::Size start_graft_;
	core::Size end_graft_;
	bool n_terminus_;
	bool pack_neighbors_;
	bool graft_loop_;
	bool find_helix_;
	bool seq_in_any_order_;
	bool pack_min_;
	utility::vector1< std::string > sequences_;
	utility::vector1< utility::vector1< core::Size > > important_residues_;
	std::set< core::Size > threadable_residues_;
	core::Real burial_cutoff_;
	std::list< utility::vector1<core::Size> > output_list_;
	core::pose::PoseOP orig_pose_;                          //OP to original pose
	core::select::residue_selector::ResidueSelectorOP threadable_residues_selector_;
};

} // switches
} // protocols

#endif
