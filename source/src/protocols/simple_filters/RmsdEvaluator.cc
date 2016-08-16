// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PoseEvaluator
/// @brief PoseEvaluator
/// @details
///
///
/// @author Oliver Lange


// Unit Headers
#include <protocols/simple_filters/RmsdEvaluator.hh>

// Package Headers

// Project Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>
#include <protocols/jd2/util.hh>
#include <protocols/loops/loops_main.hh>
// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>
#include <core/scoring/rms_util.hh>

// option key includes


#include <protocols/evaluation/util.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

//Auto Headers
static THREAD_LOCAL basic::Tracer tr( "protocols.evalution.RMSD" );
namespace protocols {
namespace simple_filters {
using namespace core;


RmsdEvaluator::RmsdEvaluator( core::pose::PoseCOP pose, Size start, Size end, std::string tag, bool bGDT /*default true*/ )
: evaluation::SingleValuePoseEvaluator< Real > ("rms"+tag ),
	rmsd_pose_( pose ),
	start_( start ),
	end_( end ),
	bGDT_ ( bGDT ),
	tag_ ( tag ),
	report_gdt_components_ ( false )
{
	runtime_assert( start >=1 );
	runtime_assert( end <= pose -> total_residue() );
}

RmsdEvaluator::~RmsdEvaluator() {}

RmsdEvaluator::RmsdEvaluator( core::pose::PoseCOP pose, std::string tag, bool bGDT /* default true */ )
: evaluation::SingleValuePoseEvaluator< Real > ("rms"+tag),
	rmsd_pose_( pose ),
	start_( 1 ),
	end_( pose->total_residue() ),
	bGDT_ ( bGDT ),
	tag_( tag ),
	report_gdt_components_ ( false )
{}


void
RmsdEvaluator::apply( core::pose::Pose& pose, std::string tag, core::io::silent::SilentStruct &pss) const {
	basic::Tracer tr( "protocols.Evaluator.RMSD" );
	tr.Debug << "compute RMSD for " << tag << " for residues " << start_ << "..." << end_ << std::endl;

	core::Real rmsd( apply( pose ) );
	pss.add_energy( "rms"+tag_, rmsd );


	if ( bGDT_ ) { //move into GDT_Evaluator
		tr.Debug << "compute GDT-score for " << tag << std::endl;
		core::Real m_1_1, m_2_2, m_3_3, m_4_3, m_7_4;
		core::Real gdtmm = core::scoring::CA_gdtmm( *rmsd_pose_, pose, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
		pss.add_energy ( "gdtmm"+tag_, gdtmm );
		if ( report_gdt_components_ ) {
			pss.add_energy ( "m11", m_1_1 );
			pss.add_energy ( "m22", m_2_2 );
			pss.add_energy ( "m33", m_3_3 );
			pss.add_energy ( "m43", m_4_3 );
			pss.add_energy ( "m74", m_7_4 );
		}
		tr.Debug << "compute maxsub for " << tag << std::endl;
		int maxsub = core::scoring::CA_maxsub( *rmsd_pose_, pose );
		pss.add_energy( "maxsub", maxsub );
	}
}

Real RmsdEvaluator::apply( core::pose::Pose& pose ) const {
	runtime_assert( rmsd_pose_ != 0 );
	core::Real rmsd;
	if ( start_ == 1 && end_ == rmsd_pose_->total_residue() ) {
		runtime_assert( pose.total_residue() >= end_ );
		rmsd = core::scoring::CA_rmsd( *rmsd_pose_, pose );
	} else {
		runtime_assert( pose.total_residue() >= end_ );
		rmsd = core::scoring::CA_rmsd( *rmsd_pose_, pose, start_, end_ );
	}
	return rmsd;
}


SelectRmsdEvaluator::SelectRmsdEvaluator( core::pose::PoseCOP pose, std::list< Size > const& selection, std::string tag, bool CAonly )
: evaluation::SingleValuePoseEvaluator< Real >( "rms"+tag ),
	rmsd_pose_( pose ),
	selection_( selection ),
	tag_ ( tag ),
	CAonly_( CAonly )
{

}

SelectRmsdEvaluator::SelectRmsdEvaluator( core::pose::PoseCOP pose, utility::vector1< Size> const& selection, std::string tag, bool CAonly )
: evaluation::SingleValuePoseEvaluator< Real >( "rms"+tag ),
	rmsd_pose_( pose ),
	tag_( tag ),
	CAonly_( CAonly )
{
	copy( selection.begin(), selection.end(), std::back_inserter( selection_ ) );
}


SelectRmsdEvaluator::SelectRmsdEvaluator( core::pose::PoseCOP pose, std::string tag, bool CAonly )
: evaluation::SingleValuePoseEvaluator< Real >( "rms"+tag ),
	rmsd_pose_( pose ),
	tag_( tag ),
	CAonly_( CAonly )
{
	if ( pose ) evaluation::find_existing_residues( pose, tag, selection_ );
}

SelectRmsdEvaluator::SelectRmsdEvaluator( core::pose::Pose const& pose, std::string tag, bool CAonly  )
: evaluation::SingleValuePoseEvaluator< Real >( "rms"+tag ),
	rmsd_pose_( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( pose ) ) ) ),
	tag_( tag ),
	CAonly_( CAonly )
{
	evaluation::find_existing_residues( rmsd_pose_, tag, selection_ );
}

Real
SelectRmsdEvaluator::apply( core::pose::Pose& pose ) const {
	core::pose::PoseCOP target_pose = rmsd_pose_;
	if ( !target_pose ) {
		runtime_assert( jd2::jd2_used() );
		target_pose = jd2::get_current_jobs_starting_pose();
	}
	if ( !target_pose ) utility_exit_with_message(" no target pose for rmsd simple_filters "+tag_ );
	core::Real rmsd;
	if ( selection_.size() > 0 ) {
		if ( CAonly_ ) {
			rmsd = core::scoring::CA_rmsd( *target_pose, pose, selection_ );
		} else {
			rmsd = core::scoring::all_atom_rmsd( *target_pose, pose, selection_ );
		}
	} else {
		if ( CAonly_ ) {
			rmsd = core::scoring::CA_rmsd( *target_pose, pose );
		} else {
			rmsd = core::scoring::all_atom_rmsd( *target_pose, pose );
		}
	}
	return rmsd;
}

SelectGdtEvaluator::~SelectGdtEvaluator(){}
SelectGdtEvaluator::SelectGdtEvaluator( core::pose::PoseCOP pose, std::list< Size > const& selection, std::string tag )
: evaluation::SingleValuePoseEvaluator< Real >( "gdtmm"+tag ),
	rmsd_pose_( pose ),
	selection_( selection ),
	tag_ ( tag )
{

}

SelectGdtEvaluator::SelectGdtEvaluator( core::pose::PoseCOP pose, utility::vector1< Size> const& selection, std::string tag )
: evaluation::SingleValuePoseEvaluator< Real >( "gdtmm"+tag ),
	rmsd_pose_( pose ),
	tag_( tag )
{
	copy( selection.begin(), selection.end(), std::back_inserter( selection_ ) );
}

SelectGdtEvaluator::SelectGdtEvaluator( core::pose::PoseCOP pose, std::string tag )
: evaluation::SingleValuePoseEvaluator< Real >( "gdtmm"+tag ),
	rmsd_pose_( pose ),
	tag_( tag )
{
	if ( pose ) evaluation::find_existing_residues( pose, tag, selection_ );
}

SelectGdtEvaluator::SelectGdtEvaluator( core::pose::Pose const& pose, std::string tag )
: evaluation::SingleValuePoseEvaluator< Real >( "gdtmm"+tag ),
	rmsd_pose_( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( pose ) ) ) ),
	tag_( tag )
{
	evaluation::find_existing_residues( rmsd_pose_, tag, selection_ );
}


Real
SelectGdtEvaluator::apply( core::pose::Pose& pose ) const {
	core::pose::PoseCOP target_pose = rmsd_pose_;
	core::Real m_1_1, m_2_2, m_3_3, m_4_3, m_7_4;
	if ( !target_pose ) {
		runtime_assert( jd2::jd2_used() );
		target_pose = jd2::get_current_jobs_starting_pose();
	}
	if ( !target_pose ) utility_exit_with_message(" no target pose for rmsd simple_filters "+tag_ );

	if ( selection_.size() ) {
		return core::scoring::CA_gdtmm( *rmsd_pose_, pose, selection_, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
	}
	return core::scoring::CA_gdtmm( *rmsd_pose_, pose, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
}


SelectMaxsubEvaluator::SelectMaxsubEvaluator( core::pose::PoseCOP pose, std::list< Size > const& selection, std::string tag, core::Real rmsd_threshold )
: evaluation::SingleValuePoseEvaluator< Real >( "maxsub"+tag ),
	rmsd_pose_( pose ),
	selection_( selection ),
	tag_ ( tag ),
	rmsd_threshold_( rmsd_threshold )
{

}

SelectMaxsubEvaluator::SelectMaxsubEvaluator( core::pose::PoseCOP pose, utility::vector1< Size> const& selection, std::string tag, core::Real rmsd_threshold )
: evaluation::SingleValuePoseEvaluator< Real >( "maxsub"+tag ),
	rmsd_pose_( pose ),
	tag_( tag ),
	rmsd_threshold_( rmsd_threshold )
{
	copy( selection.begin(), selection.end(), std::back_inserter( selection_ ) );
}

SelectMaxsubEvaluator::SelectMaxsubEvaluator( core::pose::PoseCOP pose, std::string tag, core::Real rmsd_threshold )
: evaluation::SingleValuePoseEvaluator< Real >( "maxsub"+tag ),
	rmsd_pose_( pose ),
	tag_( tag ),
	rmsd_threshold_( rmsd_threshold )
{
	if ( pose ) evaluation::find_existing_residues( pose, tag, selection_ );
}

SelectMaxsubEvaluator::SelectMaxsubEvaluator( core::pose::Pose const& pose, std::string tag, core::Real rmsd_threshold )
: evaluation::SingleValuePoseEvaluator< Real >( "maxsub"+tag ),
	rmsd_pose_( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( pose ) ) ) ),
	tag_( tag ),
	rmsd_threshold_( rmsd_threshold )
{
	evaluation::find_existing_residues( rmsd_pose_, tag, selection_ );
}


Real
SelectMaxsubEvaluator::apply( core::pose::Pose& pose ) const {
	runtime_assert( rmsd_pose_ != 0 );
	//  core::Real m_1_1, m_2_2, m_3_3, m_4_3, m_7_4;
	core::Real maxsub = core::scoring::CA_maxsub( *rmsd_pose_, pose, selection_, rmsd_threshold_ );
	return maxsub;
}


SymmetricRmsdEvaluator::SymmetricRmsdEvaluator( core::pose::PoseCOP pose, std::string tag )
: evaluation::SingleValuePoseEvaluator< Real > ("symmetric_rms"+tag ),
	rmsd_pose_( pose ) {}

SymmetricRmsdEvaluator::~SymmetricRmsdEvaluator(){}

Real
SymmetricRmsdEvaluator::apply( core::pose::Pose& pose ) const {

	return core::scoring::CA_rmsd_symmetric( *rmsd_pose_, pose );
}

LoopRmsdEvaluator::LoopRmsdEvaluator( core::pose::PoseCOP pose, protocols::loops::Loops loops, std::string tag, bool CAonly, bool superimpose )
: evaluation::SingleValuePoseEvaluator< Real >( "looprms"+tag ),
	rmsd_pose_( pose ),
	loops_( loops ),
	CAonly_( CAonly ),
	superimpose_( superimpose )
{}

LoopRmsdEvaluator::LoopRmsdEvaluator( core::pose::PoseCOP pose, protocols::loops::Loops loops, protocols::loops::Loops core, std::string tag, bool CAonly, bool superimpose )
: evaluation::SingleValuePoseEvaluator< Real >( "looprms"+tag ),
	rmsd_pose_( pose ),
	loops_( loops ),
	core_( core ),
	CAonly_( CAonly ),
	superimpose_( superimpose )
{}


Real
LoopRmsdEvaluator::apply( core::pose::Pose& pose ) const {
	core::pose::PoseCOP target_pose = rmsd_pose_;
	bool temp_pose=false;
	if ( !target_pose ) {
		runtime_assert( jd2::jd2_used() );
		temp_pose=true;
		target_pose = jd2::get_current_jobs_starting_pose();
	}
	if ( !target_pose ) utility_exit_with_message(" no target pose for rmsd simple_filters "+name(0) );
	core::Real rmsd;
	if ( superimpose_ ) {
		rmsd = protocols::loops::loop_rmsd_with_superimpose_core( *target_pose, pose, loops_, core_, CAonly_ /*CA_only*/, false /*bb_only*/ );
	} else {
		rmsd = protocols::loops::loop_rmsd( *target_pose, pose, loops_, CAonly_ /*CA_only*/, false /*bb_only*/ );
	}
	if ( temp_pose ) target_pose=NULL;
	return rmsd;
}

}
}
