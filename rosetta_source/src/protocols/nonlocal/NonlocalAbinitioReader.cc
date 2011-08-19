// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/NonlocalAbinitioReader.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/NonlocalAbinitioReader.hh>

// C/C++ headers
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

// External headers
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>

// Utility headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>

// Package headers
#include <protocols/nonlocal/NLFragment.hh>
#include <protocols/nonlocal/NLFragmentGroup.hh>
#include <protocols/nonlocal/NLGrouping.hh>

namespace protocols {
namespace nonlocal {

void NonlocalAbinitioReader::read(const std::string& filename,
                                  utility::vector1<NLGrouping> *groups) {
  using core::Real;
  using core::Size;
  using std::string;
  using std::vector;

  assert(groups);

  std::ifstream in(filename.c_str());
  if (!in.is_open()) {
    std::cerr << "Unable to open file " << filename << std::endl;
    exit(1);
  }

  NLGrouping grouping;
  NLFragmentGroup group;

  string line;
  while (in.good()) {
    getline(in, line);

    if (line == "{") {
      grouping.clear_groups();
    } else if (line == "}") {
      // Syntactic sugar: the grammar actually requires that the user end a fragment
      // group with an empty line. This is troublesome for the last FragmentGroup in
      // a NonlocalGrouping. The last two lines would then be:
      //
      // ...
      //
      // }
      emit_group(&group, &grouping);

      // Recursively sort the contents of the non-local grouping. This ensures
      // a consistent ordering of elements. Specifically, NLGroupings are sorted
      // according to the starting position of their NLFragments.
      grouping.sort();
      groups->push_back(grouping);
    } else if (line == "") {
      emit_group(&group, &grouping);
    } else {
      vector<string> tokens;
      boost::split(tokens, line, boost::is_any_of(","));
      vector<string>::iterator i;

      // remove leading and trailing spaces from each token
      for (i = tokens.begin(); i != tokens.end(); ++i)
        boost::trim(*i);

      if (tokens.size() != 7) {
        std::cerr << "Invalid line: " << line << std::endl;
        exit(1);
      }

      Size pos = utility::string2int(tokens[0]);
      Real phi = utility::string2float(tokens[1]);
      Real psi = utility::string2float(tokens[2]);
      Real omega = utility::string2float(tokens[3]);
      Real x = utility::string2float(tokens[4]);
      Real y = utility::string2float(tokens[5]);
      Real z = utility::string2float(tokens[6]);

      group.add_entry(NLFragment(pos, phi, psi, omega, x, y, z));
    }
  }
  in.close();
}

void NonlocalAbinitioReader::emit_group(NLFragmentGroup* group,
                                        NLGrouping* grouping) {
  assert(group);
  assert(grouping);

  if (group->num_entries() > 0) {
    group->sort();
    grouping->add_group(*group);
    group->clear_entries();
  }
}

}  // namespace nonlocal
}  // namespace protocols
