// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/FavorSymmetricSequence.cc
/// @brief apply constraints to enforce a symmetric sequence
/// @author Sam DeLuca

#include <protocols/simple_moves/FavorSymmetricSequence.hh>
#include <protocols/simple_moves/FavorSymmetricSequenceCreator.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ResidueTypeLinkingConstraint.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

namespace protocols {
namespace simple_moves {

static thread_local basic::Tracer TR( "protocols.simple_moves.FavorSymmetricSequence" );

std::string FavorSymmetricSequenceCreator::keyname() const
{
	return FavorSymmetricSequenceCreator::mover_name();
}

protocols::moves::MoverOP FavorSymmetricSequenceCreator::create_mover() const
{
	return protocols::moves::MoverOP( new FavorSymmetricSequence );
}

std::string FavorSymmetricSequenceCreator::mover_name()
{
	return "FavorSymmetricSequence";
}

FavorSymmetricSequence::FavorSymmetricSequence() : symmetric_units_(0), penalty_(0.0)
{

}

FavorSymmetricSequence::FavorSymmetricSequence(core::Size symmetric_units, core::Real penalty) : symmetric_units_(symmetric_units), penalty_(penalty)
{

}

FavorSymmetricSequence::FavorSymmetricSequence(FavorSymmetricSequence const & src) : Mover(src),
	symmetric_units_(src.symmetric_units_),
	penalty_(src.penalty_)
{

}

protocols::moves::MoverOP FavorSymmetricSequence::clone() const
{
	return protocols::moves::MoverOP( new FavorSymmetricSequence(*this) );
}

void FavorSymmetricSequence::apply(core::pose::Pose & pose)
{
	//Probably a faster way of doing this but it doesn't really matter I don't think
	core::Size residue_count = pose.n_residue();
	for ( core::Size rsd1_index = 1; rsd1_index <= residue_count; ++rsd1_index ) {
		for ( core::Size rsd2_index = rsd1_index; rsd2_index <= residue_count; ++rsd2_index ) {
			if ( rsd1_index % (residue_count/symmetric_units_) == rsd2_index % (residue_count/symmetric_units_) && rsd1_index != rsd2_index ) {
				pose.add_constraint(core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::ResidueTypeLinkingConstraint(pose,rsd1_index,rsd2_index,penalty_) ) ));
				if ( TR.visible(basic::t_debug) ) {
					TR.Debug << "Enforcing Sequence Symmetry between residues: " << rsd1_index << " and " << rsd2_index << std::endl;
				}
			}
		}
	}
}

std::string FavorSymmetricSequence::get_name() const
{
	return "FavorSymmetricSequence";
}

void FavorSymmetricSequence::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	if ( ! tag->hasOption("penalty") ) throw utility::excn::EXCN_RosettaScriptsOption("'FavorSymmetricSequence' mover requires penalty tag");
	if ( ! tag->hasOption("symmetric_units") ) throw utility::excn::EXCN_RosettaScriptsOption("'FavorSymmetricSequence' mover requires symmetric_units tag");

	penalty_ = tag->getOption<core::Real>("penalty");
	symmetric_units_ = tag->getOption<core::Size>("symmetric_units");

	if ( symmetric_units_ < 2 ) {
		throw utility::excn::EXCN_RosettaScriptsOption("'FavorSymmetricSequence' symmetric_units tag should specify at least 2 symmetric units");
	}

	for ( std::map< std::string, utility::pointer::ReferenceCountOP >::const_iterator it = (data)[ "scorefxns" ].begin(); it!=(data)[ "scorefxns" ].end(); ++it ) {
		core::scoring::ScoreFunctionOP scorefxn( data.get_ptr< core::scoring::ScoreFunction >( "scorefxns", it->first ) );
		if ( scorefxn->get_weight( core::scoring::res_type_linking_constraint ) == 0.0 ) {
			scorefxn->set_weight( core::scoring::res_type_linking_constraint, 1.0 );
			TR<<"Setting res_type_linking_constraint weight in scorefxn "<<it->first<<" to "<<1.0<<'\n';
		}

	}
	TR<<std::endl;
}

}
}

