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

/// @author Tim Jacobs

#ifndef HELIXASSEMBLYJOB_HH_
#define HELIXASSEMBLYJOB_HH_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

//Unit headers
#include <devel/init.hh>

//core library
#include <core/init.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArrayInitializer.hh>

//utility & numeric
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/vector1.functions.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <numeric/xyzVector.hh>

//JD2
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// option key includes
#include <basic/options/util.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>

// C++ headers
#include <sstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <utility>
#include <algorithm>
#include <iterator>

class HelixAssemblyJob {
public:

  HelixAssemblyJob();
  HelixAssemblyJob(HelixAssemblyJob const & old_job);
  HelixAssemblyJob(core::Size job_id, std::string job_name, core::Size round, std::string query_structure, std::string search_structure,
      core::Size frag1_start, core::Size frag1_end, core::Size frag2_start, core::Size frag2_end);

  core::Size get_job_id() const;
  std::string get_job_name() const;
  core::Size get_round() const;
  core::Size get_frag1_end() const;
  core::Size get_frag1_start() const;
  core::Size get_frag2_end() const;
  core::Size get_frag2_start() const;
  std::string get_query_structure() const;
  std::string get_search_structure() const;
  void set_job_id(core::Size job_id);
  void set_job_name(std::string job_name);
  void set_round(core::Size round);
  void set_frag1_end(core::Size frag1_end);
  void set_frag1_start(core::Size frag1_start);
  void set_frag2_end(core::Size frag2_end);
  void set_frag2_start(core::Size frag2_start);
  void set_query_structure(std::string query_structure);
  void set_search_structure(std::string search_structure);

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
      ar & round_;
      ar & query_structure_;
      ar & search_structure_;
      ar & frag1_start_;
      ar & frag1_end_;
      ar & frag2_start_;
      ar & frag2_end_;
    }

  core::Size job_id_;
  std::string job_name_;
  core::Size round_;
  std::string query_structure_;
  std::string search_structure_;
  core::Size frag1_start_;
  core::Size frag1_end_;
  core::Size frag2_start_;
  core::Size frag2_end_;
};

#endif /* HELIXASSEMBLYJOB_HH_ */
