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
#include <protocols/ligand_docking/Rotates.hh>
#include <protocols/ligand_docking/RotatesCreator.hh>
#include <protocols/ligand_docking/Rotate.hh>
#include <protocols/ligand_docking/DistributionMap.hh>

#include <protocols/moves/DataMap.hh>
#include <core/pose/util.hh> // includes Pose.hh

// Utility Headers
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>

#include <algorithm>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static basic::Tracer rotates_tracer("protocols.ligand_docking.ligand_options.rotates", basic::t_debug);

std::string
RotatesCreator::keyname() const
{
	return RotatesCreator::mover_name();
}

protocols::moves::MoverOP
RotatesCreator::create_mover() const {
	return new Rotates;
}

std::string
RotatesCreator::mover_name()
{
	return "Rotates";
}

///@brief
Rotates::Rotates(): Mover("Rotates")
{}

Rotates::Rotates(Rotates const & that):
		//utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that )
{}

Rotates::~Rotates() {}

protocols::moves::MoverOP Rotates::clone() const {
	return new Rotates( *this );
}

protocols::moves::MoverOP Rotates::fresh_instance() const {
	return new Rotates;
}

std::string Rotates::get_name() const{
	return "Rotates";
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
Rotates::parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & /*data_map*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose
)
{
	if ( tag->getName() != "Rotates" ){
		utility_exit_with_message("This should be impossible");
	}
	if ( ! tag->hasOption("chain") ) utility_exit_with_message("'Rotates' mover requires 'chain' tag");
	if ( ! tag->hasOption("distribution") ) utility_exit_with_message("'Rotates' mover requires 'distribution' tag");
	if ( ! tag->hasOption("degrees") ) utility_exit_with_message("'Rotates' mover requires 'degrees' tag");
	if ( ! tag->hasOption("cycles") ) utility_exit_with_message("'Rotates' mover requires 'cycles' tag");

	std::string const chain = tag->getOption<std::string>("chain");
	utility::vector1<core::Size> chain_ids = core::pose::get_chain_ids_from_chain(chain, pose);
	std::string const distribution_str= tag->getOption<std::string>("distribution");
	Distribution distribution= get_distribution(distribution_str);
	core::Size const degrees = tag->getOption<core::Size>("degrees");
	core::Size const cycles = tag->getOption<core::Size>("cycles");

	foreach(core::Size chain_id, chain_ids){
		Rotate_info rotate_info;
		rotate_info.chain = chain;
		rotate_info.chain_id= chain_id;
		rotate_info.jump_id= core::pose::get_jump_id_from_chain_id(chain_id, pose);
		rotate_info.distribution= distribution;
		rotate_info.degrees = degrees;
		rotate_info.cycles = cycles;
		rotates_.push_back( new Rotate(rotate_info) );
	}
}

void Rotates::apply(core::pose::Pose & pose){
	foreach(RotateOP rotate, rotates_){
		rotate->apply(pose);
	}
}



} //namespace ligand_docking
} //namespace protocols
