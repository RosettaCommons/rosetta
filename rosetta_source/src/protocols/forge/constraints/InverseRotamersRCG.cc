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

#include <protocols/forge/constraints/InverseRotamersRCG.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/id/SequenceMapping.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


static basic::Tracer tr( "protocols.forge.constraints.InverseRotamersRCG" );

namespace protocols{
namespace forge{
namespace constraints{

InverseRotamersRCG::InverseRotamersRCG(
	core::Size lstart,
	core::Size lstop,
	std::list< core::conformation::ResidueCOP > const & inverse_rotamers
	) : constraint_func_(NULL), func_sd_(0.4)
{
	intervals_.clear();
	inverse_rotamers_.clear();
	intervals_.push_back( forge::build::Interval( lstart, lstop ) );
	for( std::list< core::conformation::ResidueCOP >::const_iterator rot_it( inverse_rotamers.begin() ), rot_end( inverse_rotamers.end() );
			 rot_it != rot_end; ++rot_it ){
		inverse_rotamers_.push_back( *rot_it );
	}
}

InverseRotamersRCG::~InverseRotamersRCG(){}

void
InverseRotamersRCG::generate_remodel_constraints(
	core::pose::Pose const & pose )
{
	//using namespace core::scoring::constraints;
	//safeguard against bad user input
	if( inverse_rotamers_.size() == 0 ){
		std::cerr << "WARNING: InverseRotamersRCG is asked to produce constraints but was not given any inverse rotamers. Something's probably wrong somewhere." << std::endl;
		return;
	}

	//if no constraint func has been set, we'll create a default one
	if( !constraint_func_ ){
		constraint_func_ = new core::scoring::constraints::BoundFunc( 0, 0.05, func_sd_, "invrot");
	}
	utility::vector1< core::Size > seqpos;
	for( core::Size i(1); i <= intervals_.size(); ++i ){
		//eventually remap intervals according to vlb seqmap
		if( this->seqmap() ){
			intervals_[i].left = (*(this->seqmap() ))[ intervals_[i].left ];
			intervals_[i].right = (*(this->seqmap() ))[ intervals_[i].right ];
		}
		for( core::Size remres( intervals_[i].left ); remres <= intervals_[i].right; ++remres ){
			seqpos.push_back( remres );
		}
	}
	this->add_constraint( protocols::toolbox::match_enzdes_util::constrain_pose_res_to_invrots( inverse_rotamers_, seqpos, pose, constraint_func_ ) );

	//we can probably delete the inverse rotamers now, to save some memory
	this->clear_inverse_rotamers();
}

void
InverseRotamersRCG::set_constraint_func(
	core::scoring::constraints::FuncOP constraint_func ){
	constraint_func_ = constraint_func;
}

void
InverseRotamersRCG::clear_inverse_rotamers()
{
	inverse_rotamers_.clear();
}

} //namespace remodel
} //namespace forge
} //namespace protocols
