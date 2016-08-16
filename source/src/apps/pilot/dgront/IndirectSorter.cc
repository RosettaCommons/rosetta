// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <iostream>
#include <algorithm>
#include <vector>

#include <utility/IndirectSorter.hh>
#include <utility/vector1.hh>

#include <core/types.hh>


class SimpleComparison {
  public:
  bool operator() (double i,double j) { return (i<j);}
};

int main(int , char ** ) {

  using namespace utility;
  using namespace core;

  vector1<Real> data;
  vector1<Size> order;
  for(int i=0;i<10;++i) {
    double r = ((double)(rand() % 1000 ) / 100.0);
    data.push_back(r);
  }

  SimpleComparison comp;
  IndirectSorter< vector1<Real>, SimpleComparison > *sorter = new
    IndirectSorter< vector1<Real>, SimpleComparison > (data,comp);
  sorter->sort(order);
  for(int i=0;i<10;++i) {
    std::cout << i<<" "<<order[i]<<" "<<data[order[i]]<<std::endl;
  }
}
