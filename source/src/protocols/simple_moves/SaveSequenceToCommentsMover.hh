// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#ifndef INCLUDED_protocols_simple_moves_SaveSequenceToCommentsMover_hh
#define INCLUDED_protocols_simple_moves_SaveSequenceToCommentsMover_hh

#include <protocols/moves/Mover.hh>

namespace protocols {
namespace simple_moves {

class SaveSequenceToCommentsMover : public moves::Mover {
public:
	SaveSequenceToCommentsMover();
	~SaveSequenceToCommentsMover() override;

	void apply( Pose & ) override;
	// XRW TEMP  std::string get_name() const override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void
	parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, Filters_map const &, moves::Movers_map const &, Pose const & ) override;

	std::string save_seq_name() const{ return save_seq_name_; }
	void save_seq_name( std::string const & s ){ save_seq_name_ = s; };

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	std::string save_seq_name_;

};

} // moves
} // protocols

#endif
