// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/constraints/InverseRotamersRCG.cc
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, april 2010
/// @modified Tom Linsky, tlinsky@uw.edu

//unit headers
#include <protocols/forge/constraints/InverseRotamersRCG.hh>
#include <protocols/forge/constraints/InverseRotamersCstGeneratorCreator.hh>

//protocol headers
#include <protocols/forge/build/Interval.hh>

//project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/id/SequenceMapping.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>

//utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.forge.constraints.InverseRotamersRCG" );

namespace protocols {
namespace forge {
namespace constraints {

protocols::moves::MoverOP
InverseRotamersCstGeneratorCreator::create_mover() const
{
	return protocols::moves::MoverOP( new InverseRotamersRCG() );
}

std::string
InverseRotamersCstGeneratorCreator::keyname() const
{
	return InverseRotamersCstGeneratorCreator::mover_name();
}

std::string
InverseRotamersCstGeneratorCreator::mover_name()
{
	return "InverseRotamersCstGenerator";
}

InverseRotamersRCG::InverseRotamersRCG()
: RemodelConstraintGenerator(),
	constraint_func_( /* NULL */ ),
	func_sd_( 0.4 )
{}

InverseRotamersRCG::InverseRotamersRCG( InverseRotamersRCG const & rval )
: RemodelConstraintGenerator( rval ),
	constraint_func_( rval.constraint_func_ ),
	func_sd_( rval.func_sd_ )
{}

InverseRotamersRCG::InverseRotamersRCG(
	core::Size const lstart,
	core::Size const lstop,
	std::list< core::conformation::ResidueCOP > const & inverse_rotamers )
: RemodelConstraintGenerator(),
	constraint_func_(/* NULL */),
	func_sd_(0.4)
{
	init( lstart, lstop, inverse_rotamers );
}

InverseRotamersRCG::~InverseRotamersRCG(){}

void
InverseRotamersRCG::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	RemodelConstraintGenerator::parse_my_tag( tag, data, filters, movers, pose );
	//nothing here right now
}

std::string
InverseRotamersRCG::get_name() const
{
	return InverseRotamersCstGeneratorCreator::mover_name();
}


protocols::moves::MoverOP
InverseRotamersRCG::fresh_instance() const
{
	return protocols::moves::MoverOP( new InverseRotamersRCG() );
}

protocols::moves::MoverOP
InverseRotamersRCG::clone() const
{
	return protocols::moves::MoverOP( new InverseRotamersRCG( *this ) );
}

void
InverseRotamersRCG::generate_remodel_constraints( core::pose::Pose const & pose )
{
	//tr << "Generating remodel constraints" << std::endl;
	//using namespace core::scoring::constraints;
	//safeguard against bad user input
	if ( inverse_rotamers_.size() == 0 ) {
		std::cerr << "WARNING: InverseRotamersRCG is asked to produce constraints but was not given any inverse rotamers. Something's probably wrong somewhere." << std::endl;
		return;
	}

	//if no constraint func has been set, we'll create a default one
	if ( !constraint_func_ ) {
		constraint_func_ = core::scoring::func::FuncOP( new core::scoring::constraints::BoundFunc( 0, 0.05, func_sd_, "invrot") );
	}
	utility::vector1< core::Size > seqpos;
	for ( core::Size i(1); i <= intervals_.size(); ++i ) {
		//eventually remap intervals according to vlb seqmap
		if ( this->seqmap() ) {
			intervals_[i].left = (*(this->seqmap() ))[ intervals_[i].left ];
			intervals_[i].right = (*(this->seqmap() ))[ intervals_[i].right ];
		}
		for ( core::Size remres( intervals_[i].left ); remres <= intervals_[i].right; ++remres ) {
			seqpos.push_back( remres );
		}
	}
	core::scoring::constraints::ConstraintCOPs csts;
	//tr << "adding the constraint to RCG" << std::endl;
	csts.push_back( protocols::toolbox::match_enzdes_util::constrain_pose_res_to_invrots( inverse_rotamers_, seqpos, pose, constraint_func_ ) );
	//tr << "clearing inverse rotamers" << std::endl;

	//we can probably delete the inverse rotamers now, to save some memory
	this->clear_inverse_rotamers();
	//tr << "done generating remodel constraints!" << std::endl;
	add_constraints( csts );
}

void
InverseRotamersRCG::set_constraint_func(
	core::scoring::func::FuncOP constraint_func ){
	constraint_func_ = constraint_func;
}

void
InverseRotamersRCG::clear_inverse_rotamers()
{
	inverse_rotamers_.clear();
}

void
InverseRotamersRCG::init( core::Size const lstart,
	core::Size const lstop,
	std::list< core::conformation::ResidueCOP > const & inverse_rotamers )
{
	intervals_.clear();
	inverse_rotamers_.clear();
	intervals_.push_back( forge::build::Interval( lstart, lstop ) );
	for ( std::list< core::conformation::ResidueCOP >::const_iterator rot_it( inverse_rotamers.begin() ), rot_end( inverse_rotamers.end() );
			rot_it != rot_end; ++rot_it ) {
		inverse_rotamers_.push_back( *rot_it );
	}
}


} //namespace remodel
} //namespace forge
} //namespace protocols
