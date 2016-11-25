
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file
/// @brief
/// @author will sheffler

#ifndef INCLUDED_protocols_matdes_SchemePlaceMotifsMover_HH
#define INCLUDED_protocols_matdes_SchemePlaceMotifsMover_HH

// Unit headers
#include <protocols/matdes/SchemePlaceMotifsMover.fwd.hh>

// Package headers

// project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <string>

namespace protocols {
namespace matdes {

class SchemePlaceMotifsMover : public protocols::moves::Mover {


public:
	SchemePlaceMotifsMover();
	// XRW TEMP  std::string get_name() const override { return "SchemePlaceMotifs"; }
	void apply(Pose& pose) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap,
		protocols::filters::Filters_map const & filtermap,
		protocols::moves::Movers_map const & movermap,
		core::pose::Pose const & pose
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	utility::vector1<std::string> motif_sets_;
	std::string tag_name_,dumpfile_,mode_,allowed_aas_;
	bool halt_on_error_;
	core::scoring::ScoreFunctionCOP scorefxn_ ;
	core::pack::task::TaskFactoryOP task_factory_ ;
};

} // matdes
} // protocols
#endif
