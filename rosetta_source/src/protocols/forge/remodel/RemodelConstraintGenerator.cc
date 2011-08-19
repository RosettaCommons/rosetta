// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/RemodelConstraintGenerator.cc
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, april 2009


// AUTO-REMOVED #include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/id/SequenceMapping.hh>

namespace protocols{
namespace forge{
namespace remodel{


RemodelConstraintGenerator::RemodelConstraintGenerator()
	: seqmap_(NULL), vlb_(NULL)
{}


void
RemodelConstraintGenerator::add_remodel_constraints_to_pose(
	core::pose::Pose & pose )
{

	generate_remodel_constraints(	pose );

	//safeguard against an RCG not generating anything
	if( remodel_csts_.size() == 0 ) return;

	remodel_csts_ = pose.add_constraints( remodel_csts_ );

}

void
RemodelConstraintGenerator::remove_remodel_constraints_from_pose(
	core::pose::Pose & pose
) const
{

	//safeguard against an RCG not generating anything
	if( remodel_csts_.size() == 0 ) return;

	if( ! pose.remove_constraints( remodel_csts_, true ) ){
		utility_exit_with_message("Remodel constraints somehow got lost among the way");
	}
}

void
RemodelConstraintGenerator::add_constraint(	core::scoring::constraints::ConstraintCOP cst )
{
	remodel_csts_.push_back( cst );
}

void
RemodelConstraintGenerator::add_constraints( core::scoring::constraints::ConstraintCOPs csts )
{
	for( core::scoring::constraints::ConstraintCOPs::const_iterator cst_it = csts.begin();
			 cst_it != csts.end(); ++cst_it ){
		add_constraint( (*cst_it) );
	}
}

void
RemodelConstraintGenerator::clear_constraints()
{
	remodel_csts_.clear();
}

RemodelConstraintGenerator::VarLengthBuildCAP
RemodelConstraintGenerator::vlb() const {
	return vlb_; }

void
RemodelConstraintGenerator::set_vlb(
	VarLengthBuildCAP vlb )
{
	vlb_ = vlb;
}


void
RemodelConstraintGenerator::set_seqmap(
	core::id::SequenceMappingCOP seqmap )
{
	seqmap_ = seqmap;
}

core::id::SequenceMappingCOP
RemodelConstraintGenerator:: seqmap() const {
	return seqmap_;
}



} //namespace remodel
} //namespace forge
} //namespace protocols
