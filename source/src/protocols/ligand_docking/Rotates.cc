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

#include <basic/datacache/DataMap.hh>
#include <core/pose/util.hh> // includes Pose.hh

// Utility Headers
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>

#include <algorithm>

#include <utility/excn/Exceptions.hh>
#include <boost/foreach.hpp>

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static thread_local basic::Tracer rotates_tracer( "protocols.ligand_docking.ligand_options.rotates", basic::t_debug );

std::string
RotatesCreator::keyname() const
{
	return RotatesCreator::mover_name();
}

protocols::moves::MoverOP
RotatesCreator::create_mover() const {
	return protocols::moves::MoverOP( new Rotates );
}

std::string
RotatesCreator::mover_name()
{
	return "Rotates";
}

/// @brief
Rotates::Rotates(): Mover("Rotates")
{}

Rotates::Rotates(Rotates const & that):
		//utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that )
{}

Rotates::~Rotates() {}

protocols::moves::MoverOP Rotates::clone() const {
	return protocols::moves::MoverOP( new Rotates( *this ) );
}

protocols::moves::MoverOP Rotates::fresh_instance() const {
	return protocols::moves::MoverOP( new Rotates );
}

std::string Rotates::get_name() const{
	return "Rotates";
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
Rotates::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & /*data_map*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose
)
{
	if ( tag->getName() != "Rotates" ){
		throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");
	}
	if ( ! tag->hasOption("distribution") ) throw utility::excn::EXCN_RosettaScriptsOption("'Rotates' mover requires 'distribution' tag");
	if ( ! tag->hasOption("degrees") ) throw utility::excn::EXCN_RosettaScriptsOption("'Rotates' mover requires 'degrees' tag");
	if ( ! tag->hasOption("cycles") ) throw utility::excn::EXCN_RosettaScriptsOption("'Rotates' mover requires 'cycles' tag");
	if( tag->hasOption("chain") && tag->hasOption("chains") ) throw utility::excn::EXCN_RosettaScriptsOption("'Rotates' mover cannot have both a 'chain' and a 'chains' tag");
	if( ! (tag->hasOption("chain") || tag->hasOption("chains") ) ) throw utility::excn::EXCN_RosettaScriptsOption("'Rotates' mover requires either a 'chain' or a 'chains' tag");

	utility::vector1<std::string> chain_strs;
	if( tag->hasOption("chain")){
		chain_strs.push_back( tag->getOption<std::string>("chain") );
	}
	else if( tag->hasOption("chains") ){
		std::string const chains_str = tag->getOption<std::string>("chains");
		chain_strs= utility::string_split(chains_str, ',');
	}

	std::string const distribution_str= tag->getOption<std::string>("distribution");
	Distribution distribution= get_distribution(distribution_str);
	core::Size const degrees = tag->getOption<core::Size>("degrees");
	core::Size const cycles = tag->getOption<core::Size>("cycles");

	BOOST_FOREACH(std::string chain, chain_strs){
		utility::vector1<core::Size> chain_ids = core::pose::get_chain_ids_from_chain(chain, pose);
		BOOST_FOREACH(core::Size chain_id, chain_ids){
			Rotate_info rotate_info;
			rotate_info.chain_id = chain_id;
			rotate_info.jump_id = core::pose::get_jump_id_from_chain_id(chain_id, pose);
			rotate_info.distribution= distribution;
			rotate_info.degrees = degrees;
			rotate_info.cycles = cycles;
			rotates_.push_back( protocols::ligand_docking::RotateOP( new Rotate(rotate_info) ) );
		}
	}
}

void Rotates::apply(core::pose::Pose & pose){
	BOOST_FOREACH(RotateOP rotate, rotates_){
		rotate->apply(pose);
	}
}


} //namespace ligand_docking
} //namespace protocols
