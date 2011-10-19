// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file NativeAtom.hh
///
/// @brief

/// @author Tim jacobs

#ifndef NATIVEATOM_HH_
#define NATIVEATOM_HH_

//C++ Headers
#include <string>

//External
#include <boost/serialization/access.hpp>

class NativeAtom{

public:

  NativeAtom();

  ///@brief
  NativeAtom(const unsigned short type, const float x, const float y, const float z);

  ~NativeAtom();

  float x();
  float y();
  float z();

  std::string print() const;

private:
  friend class boost::serialization::access;

  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & type_;
    ar & x_;
    ar & y_;
    ar & z_;
  }

  unsigned short type_;
  float x_;
  float y_;
  float z_;
};

#endif /* NATIVEATOM_HH_ */
