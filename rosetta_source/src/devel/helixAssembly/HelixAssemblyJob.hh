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

//Unit headers
#include <devel/helixAssembly/HelicalFragment.hh>

//external library
#include <boost/serialization/vector.hpp>

//core library
#include <core/types.hh>

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
  HelixAssemblyJob(core::Size id, std::string name, core::Size round, bool direction_needed, std::string query_structure,
      std::string search_structure, core::Size search_index, core::Size query_frag_1_index,
      core::Size query_frag_2_index, std::vector<HelicalFragment> fragments, bool first_round);

  core::Size get_id() const;
  std::string get_name() const;
  core::Size get_remaining_rounds() const;
  core::Size get_search_index() const;
  std::string get_query_structure() const;
  std::string get_search_structure() const;
  core::Size get_query_frag_1_index() const;
  core::Size get_query_frag_2_index() const;
  std::vector<HelicalFragment> get_fragments() const;
  bool get_direction_needed() const;
  bool get_first_round() const;
  HelicalFragment get_query_frag_1() const;
  HelicalFragment get_query_frag_2() const;
  void set_id(core::Size job_id);
  void set_name(std::string job_name);
  void set_remaining_rounds(core::Size remaining_rounds);
  void set_search_index(core::Size search_index);
  void set_query_structure(std::string query_structure);
  void set_search_structure(std::string search_structure);
  void set_query_frag_1_index(core::Size query_frag_1_index);
  void set_query_frag_2_index(core::Size query_frag_2_index);
  void set_fragments(std::vector<HelicalFragment> fragments);
  void set_direction_needed(bool direction);
  void set_first_round(bool first_round);
  void add_fragment(HelicalFragment new_fragment);


private:
  friend class boost::serialization::access;

  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & id_;
    ar & name_;
    ar & remaining_rounds_;
    ar & query_structure_;
    ar & search_structure_;
    ar & search_index_;
    ar & query_frag_1_index_;
    ar & query_frag_2_index_;
    ar & fragments_;
    ar & direction_needed_;
    ar & first_round_;
  }

  core::Size id_;
  std::string name_;
  core::Size remaining_rounds_;
  std::string query_structure_;
  std::string search_structure_;
  core::Size search_index_;
  core::Size query_frag_1_index_;
  core::Size query_frag_2_index_;
  std::vector<HelicalFragment> fragments_;
  bool direction_needed_; //toggles between true and false to tell HelixAssemblyMover which direction the helix should be built
  bool first_round_;
};

#endif /* HELIXASSEMBLYJOB_HH_ */
