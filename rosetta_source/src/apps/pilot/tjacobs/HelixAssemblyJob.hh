// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file HelixAssemblyJob.hh
///
/// @brief

/// @author Tim jacobs

#ifndef HELIXASSEMBLYJOB_HH_
#define HELIXASSEMBLYJOB_HH_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

class HelixAssemblyJob {
public:
    core::Size getFrag1_end() const;
    core::Size getFrag1_start() const;
    core::Size getFrag2_end() const;
    core::Size getFrag2_start() const;
    core::pose::Pose getQuery_pose() const;
    void setFrag1_end(core::Size frag1_end_);
    void setFrag1_start(core::Size frag1_start_);
    void setFrag2_end(core::Size frag2_end_);
    void setFrag2_start(core::Size frag2_start_);
    void setQuery_pose(core::pose::Pose query_pose_);
    core::Size getJob_id() const;
    core::Size getJob_name() const;
    void setJob_id(core::Size job_id_);
    void setJob_name(core::Size job_name_);



private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & job_id_;
      ar & job_name_;
      ar & query_pdb_;
      ar & frag1_start_;
      ar & frag1_end_;
      ar & frag2_start_;
      ar & frag2_end_;
    }

  core::Size job_id_;
  core::Size job_name_;
  std::string query_pdb_;
  core::Size frag1_start_;
  core::Size frag1_end_;
  core::Size frag2_start_;
  core::Size frag2_end_;
};

#endif /* HELIXASSEMBLYJOB_HH_ */
