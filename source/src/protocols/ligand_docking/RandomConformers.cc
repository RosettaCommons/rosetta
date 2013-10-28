// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocol/ligand_docking/RandomConformers.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com), adapted from the ResfileReader code
/// by Steven Lewis (smlewi@unc.edu) and Andrew Leaver-Fay

// Project Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

// Unit Headers
#include <protocols/ligand_docking/RandomConformers.hh>
#include <protocols/ligand_docking/RandomConformersCreator.hh>

#include <protocols/ligand_docking/RandomConformerMover.hh>
// AUTO-REMOVED #include <protocols/ligand_docking/UnconstrainedTorsionsMover.hh>

#include <core/pose/util.hh>

// Utility Headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

///////////////////////////////////////////////////////////////////////
///@brief

static basic::Tracer random_conformer_tracer("protocols.ligand_docking.ligand_options.RandomConformers", basic::t_debug);

std::string
RandomConformersCreator::keyname() const
{
	return RandomConformersCreator::mover_name();
}

protocols::moves::MoverOP
RandomConformersCreator::create_mover() const {
	return new RandomConformers;
}

std::string
RandomConformersCreator::mover_name()
{
	return "RandomConformers";
}

RandomConformers::RandomConformers():
		//utility::pointer::ReferenceCount(),
		Mover("RandomConformers")
{}

RandomConformers::RandomConformers(RandomConformers const & that):
		//utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that ),
		chain_(that.chain_)
{}

RandomConformers::~RandomConformers() {}

protocols::moves::MoverOP RandomConformers::clone() const {
	return new RandomConformers( *this );
}

protocols::moves::MoverOP RandomConformers::fresh_instance() const {
	return new RandomConformers;
}

std::string RandomConformers::get_name() const{
	return "RandomConformers";
}

//void RandomConformers::set_chain(std::string chain)
//{
//	chain_ = chain;
//}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
RandomConformers::parse_my_tag(
		utility::tag::TagCOP const tag,
		basic::datacache::DataMap & /*datamap*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "RandomConformers" ){
		throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");
	}
	if ( ! tag->hasOption("chain") ) throw utility::excn::EXCN_RosettaScriptsOption("'RandomConformers' mover requires chain tag");

	chain_ = tag->getOption<std::string>("chain");
}

void RandomConformers::apply(core::pose::Pose & pose) {
	core::Size chain_id= core::pose::get_chain_id_from_chain(chain_, pose);
	core::Size i = pose.conformation().chain_begin(chain_id);
	core::Size end = pose.conformation().chain_end(chain_id);

	for (; i != end; ++i) {
		apply_residue(i, pose);
	}
}

void RandomConformers::apply_residue(core::Size const residue_id, core::pose::Pose & pose) {
	using namespace protocols::moves;
	using core::conformation::ResidueOP;
	RandomConformerMoverOP rcm = new RandomConformerMover(residue_id);
	rcm->apply(pose);
/// TODO accomplish the below code within the scripter
//	UnconstrainedTorsionsMoverOP utm =
//			new UnconstrainedTorsionsMover(rcm, ligand_torsion_restraints_);
//	utm->apply(pose);
}

} //namespace ligand_docking
} //namespace protocols
