// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file HelixAssemblyJob.cc
///
/// @brief

/// @author Tim jacobs

/*
 * HelixAssemblyJob.cc
 *
 *  Created on: Aug 25, 2011
 *      Author: tjacobs2
 */

#include <apps/pilot/tjacobs/HelixAssemblyJob.hh>

core::Size HelixAssemblyJob::getJob_id() const
{
    return job_id_;
}

core::Size HelixAssemblyJob::getJob_name() const
{
    return job_name_;
}

void HelixAssemblyJob::setJob_id(core::Size job_id_)
{
    this->job_id_ = job_id_;
}

void HelixAssemblyJob::setJob_name(core::Size job_name_)
{
    this->job_name_ = job_name_;
}

core::Size HelixAssemblyJob::getFrag1_end() const
{
    return frag1_end_;
}

core::Size HelixAssemblyJob::getFrag1_start() const
{
    return frag1_start_;
}

core::Size HelixAssemblyJob::getFrag2_end() const
{
    return frag2_end_;
}

core::Size HelixAssemblyJob::getFrag2_start() const
{
    return frag2_start_;
}

core::pose::Pose HelixAssemblyJob::getQuery_pose() const
{
    return query_pose_;
}

void HelixAssemblyJob::setFrag1_end(core::Size frag1_end_)
{
    this->frag1_end_ = frag1_end_;
}

void HelixAssemblyJob::setFrag1_start(core::Size frag1_start_)
{
    this->frag1_start_ = frag1_start_;
}

void HelixAssemblyJob::setFrag2_end(core::Size frag2_end_)
{
    this->frag2_end_ = frag2_end_;
}

void HelixAssemblyJob::setFrag2_start(core::Size frag2_start_)
{
    this->frag2_start_ = frag2_start_;
}

void HelixAssemblyJob::setQuery_pose(core::pose::Pose query_pose_)
{
    this->query_pose_ = query_pose_;
}
