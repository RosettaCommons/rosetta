// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/StructureRestrictor.hh
///
/// @brief  lookup relevant chains for a structure in a table.
/// @author Matthew O'Meara

/// This should probably be a pilot app, but the way Rosetta Scripts
/// is set up, it can't be in the pilot apps


#ifndef INCLUDED_protocols_moves_StructureRestrictor_hh
#define INCLUDED_protocols_moves_StructureRestrictor_hh

#include <protocols/moves/StructureRestrictor.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <utility/tag/Tag.fwd.hh>

#include <map>
#include <string>

#include <utility/vector1.hh>


// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
using namespace pose;

namespace protocols {
namespace moves {


class StructureRestrictor : public protocols::moves::Mover {

public:
	StructureRestrictor();

	StructureRestrictor( std::string const & name);

	StructureRestrictor(StructureRestrictor const & src);

	~StructureRestrictor() override;

	MoverOP fresh_instance() const override;

	MoverOP clone() const override;

	std::string get_name() const override { return "StructureRestrictor"; }


	// So this this can be called from RosettaScripts
	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & /*data*/,
		Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		Pose const & /*pose*/ ) override;

	void
	setup_relevant_chains(
		std::string const & relevant_chains_fname,
		std::map<std::string, std::string> & chain_map
	);

	// this is a hack because poses do not have canonical names!
	std::string pose_name(Pose const & pose);

	void apply( Pose& pose ) override;

private:
	std::map<std::string, std::string> chain_map;
	std::string relevant_chains_fname;
	bool initialized;
};

}//moves
}//protocols

#endif //INCLUDED_protocols_moves_StructureRestrictor_FWD_HH
