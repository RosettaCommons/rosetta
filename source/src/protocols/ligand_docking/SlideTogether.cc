// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/SlideTogether.hh>
#include <protocols/ligand_docking/SlideTogetherCreator.hh>

#include <protocols/docking/DockingInitialPerturbation.hh>

// Utility Headers

#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>
#include <boost/foreach.hpp>

namespace protocols {
namespace ligand_docking {

static THREAD_LOCAL basic::Tracer slide_together_tracer( "protocols.ligand_docking.ligand_options.slide_together" );

std::string
SlideTogetherCreator::keyname() const
{
	return SlideTogetherCreator::mover_name();
}

protocols::moves::MoverOP
SlideTogetherCreator::create_mover() const {
	return protocols::moves::MoverOP( new SlideTogether );
}

std::string
SlideTogetherCreator::mover_name()
{
	return "SlideTogether";
}

SlideTogether::SlideTogether(){}

SlideTogether::SlideTogether(std::string const & chain): chain_(chain), jumps_(){}

SlideTogether::SlideTogether(SlideTogether const & that):
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( that ),
	chain_(that.chain_),
	jumps_(that.jumps_)
{}

SlideTogether::~SlideTogether() {}

protocols::moves::MoverOP SlideTogether::clone() const {
	return protocols::moves::MoverOP( new SlideTogether( *this ) );
}

protocols::moves::MoverOP SlideTogether::fresh_instance() const {
	return protocols::moves::MoverOP( new SlideTogether );
}

std::string SlideTogether::get_name() const{
	return "SlideTogether";
}

//@brief parse XML (specifically in the context of the parser/scripting scheme)
void
SlideTogether::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*datamap*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & pose
)
{
	if ( tag->getName() != "SlideTogether" ) throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");
	if ( ! tag->hasOption("chains") ) throw utility::excn::EXCN_RosettaScriptsOption("'SlideTogether' mover requires chains tag");

	std::string const chains_str = tag->getOption<std::string>("chains");
	utility::vector1<std::string> chain_strs= utility::string_split(chains_str, ',');
	BOOST_FOREACH ( std::string chain_str, chain_strs ) {
		utility::vector1<core::Size> chain_ids= get_chain_ids_from_chain(chain_str, pose);
		BOOST_FOREACH ( core::Size chain_id, chain_ids ) {
			core::Size jump_id= core::pose::get_jump_id_from_chain_id(chain_id, pose);
			jumps_.push_back(jump_id);
		}
	}
}

void
SlideTogether::apply( core::pose::Pose & pose ){
	slide_together_tracer<< "Applying slide_together"<< std::endl;

	// If we are not using parse_my_tags...
	if ( ! chain_.empty() ) {
		if ( jumps_.size() > 0 ) {
			utility_exit_with_message("This should be impossible");// jumps are set through parse_my_tags, chain through the 1-arg constructor
		}
		utility::vector1<core::Size> chain_ids= get_chain_ids_from_chain(chain_, pose);
		BOOST_FOREACH ( core::Size chain_id, chain_ids ) {
			core::Size jump_id= core::pose::get_jump_id_from_chain_id(chain_id, pose);
			jumps_.push_back(jump_id);
		}
	}

	protocols::docking::FaDockingSlideIntoContact slideTogether(jumps_);
	slideTogether.apply(pose);
}

} //namespace ligand_docking
} //namespace protocols
