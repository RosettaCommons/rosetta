// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/FindConsensusSequence.hh
/// @brief Takes in multiple poses from the MSDJobDistributor and finds the consensus sequence that optimizes energy of all input poses
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_FindConsensusSequence_hh
#define INCLUDED_protocols_simple_moves_FindConsensusSequence_hh

#include <protocols/simple_moves/FindConsensusSequence.fwd.hh>
#include <protocols/moves/VectorPoseMover.hh>
#include <protocols/moves/Mover.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {

class FindConsensusSequence : public moves::VectorPoseMover {

public:

	FindConsensusSequence();

	~FindConsensusSequence() override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;


	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		moves::Movers_map const &,
		Pose const & ) override;


	void apply( Pose & pose ) override;

	void parse_resfiles();

	utility::vector1< core::Size >
	parse_resfile ( core::pack::task::PackerTaskCOP design_task );

	std::string
	resfile_at ( core::Size index );

	// XRW TEMP  std::string get_name() const override;

	core::scoring::ScoreFunctionOP score_function() const;

	core::pack::task::TaskFactoryOP task_factory() const;

	void score_function( core::scoring::ScoreFunctionOP );

	void task_factory( core::pack::task::TaskFactoryOP );

	void resfiles ( utility::vector1< std::string > & resfiles );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	utility::vector1< utility::vector1< core::Size > > res_links ();

private:
	core::scoring::ScoreFunctionOP sfxn_;
	core::pack::task::TaskFactoryOP task_factory_;
	utility::vector1< std::string > resfiles_;
	utility::vector1< utility::vector1< core::Size > > res_links_;

};

} //simple_moves
} //protocols

#endif
