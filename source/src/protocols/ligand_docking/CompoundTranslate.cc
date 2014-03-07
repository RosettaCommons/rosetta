// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Gordon Lemmon (glemmon@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/CompoundTranslate.hh>
#include <protocols/ligand_docking/CompoundTranslateCreator.hh>
#include <protocols/ligand_docking/Translate.hh>
// AUTO-REMOVED #include <protocols/ligand_docking/grid_functions.hh>
// AUTO-REMOVED #include <protocols/rigid/RB_geometry.hh>
// AUTO-REMOVED #include <protocols/rigid/RigidBodyMover.hh>

#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <numeric/random/random_permutation.hh>

// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

#include <utility/excn/Exceptions.hh>
#include <boost/foreach.hpp>

//Auto Headers
using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static basic::Tracer translate_tracer("protocols.ligand_docking.CompoundTranslate", basic::t_debug);


std::string
CompoundTranslateCreator::keyname() const
{
	return CompoundTranslateCreator::mover_name();
}

protocols::moves::MoverOP
CompoundTranslateCreator::create_mover() const {
	return new CompoundTranslate;
}

std::string
CompoundTranslateCreator::mover_name()
{
	return "CompoundTranslate";
}

///@brief
CompoundTranslate::CompoundTranslate():
		//utility::pointer::ReferenceCount(),
		Mover("CompoundTranslate")
{}

CompoundTranslate::CompoundTranslate(CompoundTranslate const & that):
		//utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that ),
		translates_(that.translates_),
		randomize_order_(that.randomize_order_),
		allow_overlap_(that.allow_overlap_)
{}

CompoundTranslate::~CompoundTranslate() {}

protocols::moves::MoverOP CompoundTranslate::clone() const {
	return new CompoundTranslate( *this );
}

protocols::moves::MoverOP CompoundTranslate::fresh_instance() const {
	return new CompoundTranslate;
}

std::string CompoundTranslate::get_name() const{
	return "CompoundTranslate";
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
CompoundTranslate::parse_my_tag(
		utility::tag::TagCOP const tag,
		basic::datacache::DataMap & datamap,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
)
{
	if ( tag->getName() != "CompoundTranslate" ){
		throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");
	}
	if ( ! tag->hasOption("randomize_order")){
		throw utility::excn::EXCN_RosettaScriptsOption("CompoundTranslate needs a 'randomize_order' option");
	}
	if ( ! tag->hasOption("allow_overlap")){
		throw utility::excn::EXCN_RosettaScriptsOption("CompoundTranslate needs an 'allow_overlap' option");
	}
	{// parsing randomize_order tag
		std::string allow_overlap_string= tag->getOption<std::string>("randomize_order");
		if(allow_overlap_string == "true" || allow_overlap_string == "True")
			randomize_order_= true;
		else if(allow_overlap_string == "false" || allow_overlap_string == "False")
			randomize_order_= false;
		else throw utility::excn::EXCN_RosettaScriptsOption("'randomize_order' option takes arguments 'true' or 'false'");
	}
	{// parsing allow_overlap tag
		std::string allow_overlap_string= tag->getOption<std::string>("allow_overlap");
		if(allow_overlap_string == "true" || allow_overlap_string == "True")
			allow_overlap_= true;
		else if(allow_overlap_string == "false" || allow_overlap_string == "False")
			allow_overlap_= false;
		else throw utility::excn::EXCN_RosettaScriptsOption("'allow_overlap' option takes arguments 'true' or 'false'");
	}

	BOOST_FOREACH(utility::tag::TagCOP tag, tag->getTags()){
		std::string const name= tag->getName();
		if( name == "Translate"){
			TranslateOP translate = new Translate();
			translate->parse_my_tag(tag, datamap, filters, movers, pose);
			translates_.push_back(translate);
		}
		else if( name == "Translates"){
			if ( ! tag->hasOption("chain") ) throw utility::excn::EXCN_RosettaScriptsOption("'Translates' mover requires chain tag");
			if ( ! tag->hasOption("distribution") ) throw utility::excn::EXCN_RosettaScriptsOption("'Translates' mover requires distribution tag");
			if ( ! tag->hasOption("angstroms") ) throw utility::excn::EXCN_RosettaScriptsOption("'Translates' mover requires angstroms tag");
			if ( ! tag->hasOption("cycles") ) throw utility::excn::EXCN_RosettaScriptsOption("'Translates' mover requires cycles tag");

			std::string const chain = tag->getOption<std::string>("chain");
			utility::vector1<core::Size> chain_ids = core::pose::get_chain_ids_from_chain(chain, pose);

			BOOST_FOREACH(core::Size chain_id, chain_ids){
				Translate_info translate_info;
				translate_info.chain_id= chain_id;
				translate_info.jump_id = core::pose::get_jump_id_from_chain_id(chain_id, pose);
				std::string distribution_str= tag->getOption<std::string>("distribution");
				translate_info.distribution= get_distribution(distribution_str);
				translate_info.angstroms = tag->getOption<core::Real>("angstroms");
				translate_info.cycles = tag->getOption<core::Size>("cycles");
				if(tag->hasOption("force")){
					if(tag->getOption<std::string>("force") == "true")
						translate_info.force= true;
					else if(tag->getOption<std::string>("force") != "false")
						throw utility::excn::EXCN_RosettaScriptsOption("'force' option is true or false");
				}
				translates_.push_back(new Translate(translate_info));
			}
		}
		else{
			throw utility::excn::EXCN_RosettaScriptsOption("CompoundTranslate only takes Translate or Translates child tags");
		}
	}
}

void CompoundTranslate::apply(core::pose::Pose & pose) {
	if(randomize_order_)
		numeric::random::random_permutation(translates_, numeric::random::RG);

	std::set<core::Size> chains_to_translate;

	// TranslateOPs::iterator begin= translates_.begin(); // Unused variable causes warning.
	// TranslateOPs::iterator const end= translates_.end(); // Unused variable causes warning.

	BOOST_FOREACH(TranslateOP translate, translates_){
		core::Size chain_id= translate->get_chain_id(pose);
		chains_to_translate.insert(chain_id);
	}

	if(allow_overlap_){
		BOOST_FOREACH(TranslateOP translate, translates_){
			translate->add_excluded_chains(chains_to_translate.begin(), chains_to_translate.end());
			translate->apply(pose);
		}
	}
	else{ // remove each chain from the exclusion list so that placed chains are in the grid
		BOOST_FOREACH(TranslateOP translate, translates_){
			translate->add_excluded_chains(chains_to_translate.begin(), chains_to_translate.end());
			translate->apply(pose);
			core::Size chain_id= translate->get_chain_id(pose);
			chains_to_translate.erase(chain_id);
		}
	}
}

} //namespace ligand_docking
} //namespace protocols
