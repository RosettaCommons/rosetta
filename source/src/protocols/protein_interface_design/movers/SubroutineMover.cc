// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/SubroutineMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/SubroutineMover.hh>
#include <protocols/protein_interface_design/movers/SubroutineMoverCreator.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Parser.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>

#include <protocols/moves/Mover.hh>
#include <basic/Tracer.hh>

#include <protocols/jd2/Job.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


using namespace protocols::protein_interface_design;

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.movers.SubroutineMover" );

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace protocols::moves;
using namespace core;

std::string
SubroutineMoverCreator::keyname() const
{
	return SubroutineMoverCreator::mover_name();
}

protocols::moves::MoverOP
SubroutineMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SubroutineMover );
}

std::string
SubroutineMoverCreator::mover_name()
{
	return "Subroutine";
}

protocols::moves::MoverOP
SubroutineMover::clone() const {
	return( protocols::moves::MoverOP( new SubroutineMover( *this ) ) );
}


void
SubroutineMover::apply( core::pose::Pose & pose )
{
	mover_->apply( pose );
	set_last_move_status( mover_->get_last_move_status() );
}

void
SubroutineMover::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	Movers_map const &,
	core::pose::Pose const & pose )
{
	using namespace protocols::jd2;

	std::string const xml_fname( tag->getOption< std::string >( "xml_fname" ) );

	JobOP job( JobDistributor::get_instance()->current_job() );
	ParserOP rsparser( new protocols::rosetta_scripts::RosettaScriptsParser );
	TR<<"Parsing a subroutine xml_file"<<std::endl;
	TR<<"*************WARNING: AT THIS POINT, CONSTRAINTS ADDED TO THE POSE IN A SUBROUTINE WILL BE IGNORED***********"<<std::endl;
	core::pose::Pose nonconst_pose( pose );
	rsparser->generate_mover_from_pose( job, nonconst_pose, mover_, true /*new input*/, xml_fname );
}

protocols::moves::MoverOP
SubroutineMover::fresh_instance() const {
	return protocols::moves::MoverOP( new SubroutineMover );
}

SubroutineMover::~SubroutineMover(){}

SubroutineMover::SubroutineMover() :
	Mover( SubroutineMoverCreator::mover_name() ),
	mover_( /* NULL */ )
{}

std::string
SubroutineMover::get_name() const {
	return SubroutineMoverCreator::mover_name();
}

} //movers
} //protein_interface_design
} //protocols

