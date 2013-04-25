// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file NativeResidue.hh
///
/// @brief

/// @author Tim jacobs

#ifndef NATIVERESIDUE_HH_
#define NATIVERESIDUE_HH_

//Devel
#include <devel/helixAssembly/NativeAtom.hh>

//External
#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

//C++ Headers
#include <string>
#include <vector>

class NativeResidue{

public:

  NativeResidue();

  NativeResidue(const std::string & res_type, const std::vector<NativeAtom> & atoms);

  ~NativeResidue();

  std::string res_type();

  std::string print() const;

private:
  friend class boost::serialization::access;

  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & res_type_;
    ar & atoms_;
  }

  std::string res_type_;
  std::vector<NativeAtom> atoms_;

};

#endif /* NATIVERESIDUE_HH_ */
