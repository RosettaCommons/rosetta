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
/// @author 
//
#include <string>
#include <sstream>
#include <fstream>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <devel/init.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <basic/Tracer.hh>

#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>

#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/hotspot_hashing/PlaceMinimizeSearch.hh>
#include <protocols/hotspot_hashing/SearchPattern.hh>


namespace protocols {
namespace hotspot_hashing {

static basic::Tracer tr( "protocols.hotspot_hashing.PlaceMinimizeSearch" );

PlaceMinimizeSearch::PlaceMinimizeSearch(
  core::pose::Pose const & target_pose,
  core::conformation::ResidueCOP target_residue,
  SearchPatternOP search_pattern,
	protocols::moves::MoverOP relax_mover,
	protocols::filters::FilterOP triage_filter,
	std::string output_tag) :
    target_pose_(target_pose),
    target_residue_(target_residue),
    search_pattern_(search_pattern),
		relax_mover_(relax_mover),
		triage_filter_(triage_filter),
		output_tag_(output_tag),
    current_job_(jd2::JobDistributor::get_instance()->current_job()),
    current_job_outputter_(jd2::JobDistributor::get_instance()->job_outputter())
  {
  }

PlaceMinimizeSearch::PlaceMinimizeSearch(
  core::pose::Pose target_pose,
  core::conformation::ResidueCOP target_residue,
  SearchPatternOP search_pattern,
	protocols::moves::MoverOP relax_mover,
	protocols::filters::FilterOP triage_filter,
	std::string output_tag,
  protocols::jd2::JobOP current_job,
  protocols::jd2::JobOutputterOP current_job_outputter):
    target_pose_(target_pose),
    target_residue_(target_residue),
    search_pattern_(search_pattern),
		relax_mover_(relax_mover),
		triage_filter_(triage_filter),
		output_tag_(output_tag),
    current_job_(current_job),
    current_job_outputter_(current_job_outputter)
  {
  }

void PlaceMinimizeSearch::execute()
{
  utility::vector1<TransformPair> searchpoints = search_pattern_->Searchpoints();
  tr.Info << "Initialized search pattern. Points: " << searchpoints.size() << std::endl;

  core::Size transformindex = 1;

  foreach(TransformPair transform, searchpoints)
  {
    // Clone target pose
    core::pose::Pose pose(target_pose_);

    core::Size residuejumpindex;
    core::Size residueindex;

    placeResidueAtTransform(pose, *target_residue_, transform, residuejumpindex, residueindex);

    // Add additional scoring parameters
    core::pose::add_variant_type_to_pose_residue( pose, "SHOVE_BB", residueindex );

    logPreMinPose(pose, transformindex, residuejumpindex, residueindex);

    if(relax_mover_)
    {
      relax_mover_()->apply( pose );
    }

    bool passed_triage(true);

    if(triage_filter_)
    {
      passed_triage = triage_filter_()->apply(pose);
    }

    if(passed_triage)
    {
			logPostMinPose(pose, transformindex, residuejumpindex, residueindex);
    }

		transformindex += 1;
  }
}

void PlaceMinimizeSearch::placeResidueAtTransform( core::pose::Pose & pose, core::conformation::Residue const & residue, TransformPair transform, core::Size & residuejumpindex, core::Size & residueindex )
{
  // Places residue at last jump & residue number
  placeResidueOnPose(pose, residue);
  residueindex = pose.total_residue();
  residuejumpindex = pose.num_jump();

  core::kinematics::Stub upstreamstub = pose.conformation().upstream_jump_stub(residuejumpindex);
  core::kinematics::Jump residuejump = pose.jump(residuejumpindex);
  core::kinematics::Jump newjump(residuejump);

  // Place residue in "canonical" position
  // Calculate transforms to move target atom to 0,0,0 and rotate into canonical position
  TransformPair residuetransform = residueStubCentroidTransform(pose.residue( residueindex ));

  if (residuetransform.translation.length() != 0)
  {
    newjump.translation_along_axis(upstreamstub, residuetransform.translation, residuetransform.translation.length());
  }

  newjump.rotation_by_matrix(upstreamstub, Vector(), residuetransform.rotation);

  // Apply target transformation
  newjump.rotation_by_matrix(upstreamstub, Vector(), transform.rotation);
  if (transform.translation.length() != 0)
  {
    newjump.translation_along_axis(upstreamstub, transform.translation, transform.translation.length());
  }

  pose.set_jump(residuejumpindex, newjump);

  tr.Debug << "Placed residue at anchor location: " << pose.residue( residueindex ).xyz("CA") << std::endl;
}

TransformPair PlaceMinimizeSearch::residueStubCentroidTransform(core::conformation::Residue const & residue)
{
  // Canonical transform aligns CA atom to <0, 0, 0>
  // CA->SC heavyatom centroid vector along <1,0,0>
  // CA->C vector on the XY plane (<CA->C> * <0,0,1> == 0)
  
  Vector position = -residue.xyz("CA");
  
  Vector xunit = Vector(1, 0, 0);
  Vector yunit = Vector(0, 1, 0);

  Vector cacentroid_vector = residueStubCentroid(residue) + position;
  Matrix cacentroid_rotation = rotation_matrix( cacentroid_vector.cross(xunit), angle_of(cacentroid_vector, xunit));

  Vector cac_vector_prime = cacentroid_rotation * (residue.xyz("C") + position);
  Vector cac_vector_zyprojection = Vector(0, cac_vector_prime.y(), cac_vector_prime.z());
  Matrix cac_rotation = rotation_matrix( cac_vector_zyprojection.cross(yunit), angle_of(cac_vector_zyprojection, yunit));

  return TransformPair(position, cac_rotation * cacentroid_rotation);
}

Vector PlaceMinimizeSearch::residueStubCentroid(core::conformation::Residue const & residue)
{
  Vector centroid;
  centroid = 0;

  if (residue.first_sidechain_atom() > residue.nheavyatoms())
  {
    //TODO Generate pseudocentroid from mainchain atoms
    return centroid;
  }

  for (core::Size i = residue.first_sidechain_atom(); i <= residue.nheavyatoms(); ++i)
  {
    centroid += residue.xyz(i);
  }

  centroid /= (1 + residue.nheavyatoms() - residue.first_sidechain_atom());

  return centroid;
}

void PlaceMinimizeSearch::placeResidueOnPose(core::pose::Pose & pose, core::conformation::Residue const & residue)
{
  pose.append_residue_by_jump( residue, pose.total_residue(), "", residue.atom_name(residue.nbr_atom()), true );

  tr.Debug << "Appended residue on pose. Residue: " << pose.total_residue() << " Jump: " << pose.num_jump() << " Anchor atom: " << residue.atom_name(residue.nbr_atom()) << std::endl;
}

void PlaceMinimizeSearch::logPreMinPose(core::pose::Pose & pose, core::Size transformindex, core::Size residuejumpindex, core::Size residueindex)
{
  // Copy subset of the full pose
  core::pose::Pose residuePose(pose, residueindex, residueindex);
	std::stringstream tag (std::stringstream::in | std::stringstream::out);
	tag << output_tag_ << ".pre." << transformindex;

	current_job_outputter_->other_pose( current_job_, residuePose, tag.str());
}

void PlaceMinimizeSearch::logPostMinPose(core::pose::Pose & pose, core::Size transformindex, core::Size residuejumpindex, core::Size residueindex)
{
  // Copy subset of the full pose
  core::pose::Pose residuePose(pose, residueindex, residueindex);
	std::stringstream tag (std::stringstream::in | std::stringstream::out);
	tag << output_tag_ << ".post." << transformindex;

	current_job_outputter_->other_pose( current_job_, residuePose, tag.str());
}

}
}
