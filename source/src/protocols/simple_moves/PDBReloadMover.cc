// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   PDBReloadMover.cc
///
/// @brief
/// @author Javier Castellanos (javiercv@uw.edu)

// unit headers
#include <protocols/simple_moves/PDBReloadMover.hh>
#include <protocols/simple_moves/PDBReloadMoverCreator.hh>

// type headers
#include <core/types.hh>
#include <core/id/types.hh>

// project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

// utility header
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <utility/tag/Tag.hh>

#include <sstream>

//option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>

namespace protocols {
namespace simple_moves {

static basic::Tracer TR("protocols.simple_moves.PDBReloadMover");

std::string
PDBReloadMoverCreator::keyname() const
{
	return PDBReloadMoverCreator::mover_name();
}

protocols::moves::MoverOP
PDBReloadMoverCreator::create_mover() const {
	return new PDBReloadMover;
}

std::string
PDBReloadMoverCreator::mover_name()
{
	return "PDBReload";
}

PDBReloadMover::PDBReloadMover() :
	protocols::moves::Mover("PDBReloadMover")
{}

PDBReloadMover::~PDBReloadMover() {}

protocols::moves::MoverOP
PDBReloadMover::clone() const
{
	return new PDBReloadMover( *this );
}

protocols::moves::MoverOP
PDBReloadMover::fresh_instance() const
{
	return new PDBReloadMover();
}

void
PDBReloadMover::apply( Pose & pose ) {
	std::ostringstream ss;	
	pose.dump_pdb(ss);
	core::import_pose::pose_from_pdbstring(pose, ss.str());
}

std::string
PDBReloadMover::get_name() const {
	return "PDBReloadMover";
}

void
PDBReloadMover::parse_my_tag( utility::tag::TagCOP const /*tag*/,
		basic::datacache::DataMap & /*data_map*/,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & )
{
}

} // moves
} // protocols
