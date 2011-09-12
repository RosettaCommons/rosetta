// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <basic/Tracer.hh>
#include <core/types.hh>

//Auto Headers
#include <core/pose/util.hh>
#include <iostream>


namespace protocols {
namespace ligand_docking {

static basic::Tracer slide_together_tracer("protocols.ligand_docking.ligand_options.slide_together");

std::string
SlideTogetherCreator::keyname() const
{
	return SlideTogetherCreator::mover_name();
}

protocols::moves::MoverOP
SlideTogetherCreator::create_mover() const {
	return new SlideTogether;
}

std::string
SlideTogetherCreator::mover_name()
{
	return "SlideTogether";
}

SlideTogether::SlideTogether(){}

SlideTogether::SlideTogether(SlideTogether const & that):
	    //utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that ),
		chain_(that.chain_)
{}

SlideTogether::~SlideTogether() {}

protocols::moves::MoverOP SlideTogether::clone() const {
	return new SlideTogether( *this );
}

protocols::moves::MoverOP SlideTogether::fresh_instance() const {
	return new SlideTogether;
}

std::string SlideTogether::get_name() const{
	return "SlideTogether";
}

//@brief parse XML (specifically in the context of the parser/scripting scheme)
void
SlideTogether::parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & /*datamap*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "SlideTogether" ) utility_exit_with_message("This should be impossible");
	if ( ! tag->hasOption("chain") ) utility_exit_with_message("'SlideTogether' mover requires chain tag");

	chain_= tag->getOption<std::string>("chain");//.c_str() ;
}

void
SlideTogether::apply( core::pose::Pose & pose ){
	core::Size jump_id= core::pose::get_jump_id_from_chain(chain_, pose);
	slide_together_tracer<< "chain "<< chain_ << std::endl;
	protocols::docking::FaDockingSlideIntoContact slideTogether(jump_id);
	slideTogether.apply(pose);
}

} //namespace ligand_docking
} //namespace protocols
