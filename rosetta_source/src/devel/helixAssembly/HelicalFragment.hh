// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file HelicalFragment.hh
///
/// @brief

/// @author Tim jacobs

#ifndef HELICALFRAGMENT_HH_
#define HELICALFRAGMENT_HH_

//Core
#include <core/types.hh>

//External
#include <boost/serialization/access.hpp>

//C++ Headers
#include <string>

class HelicalFragment{

public:

  HelicalFragment();

  HelicalFragment(core::Size start, core::Size end);

  virtual ~HelicalFragment();

  core::Size get_end() const;
  std::string get_pdb_source() const;
  core::Size get_start() const;
  core::Size get_size() const;
  void set_end(core::Size end_);
  void set_pdb_source(std::string pdb_source_);
  void set_start(core::Size start_);

private:
  friend class boost::serialization::access;

  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & start_;
    ar & end_;
    ar & pdb_source_;
  }

  core::Size start_;
  core::Size end_;
  std::string pdb_source_;

};

#endif /* HELICALFRAGMENT_HH_ */
