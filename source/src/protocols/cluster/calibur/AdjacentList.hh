// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file external/calibur/AdjacentList.hh
/// @author SC Li & YK Ng (kalngyk@gmail.com)

#ifndef external_calibur_AdjacentList_HH
#define external_calibur_AdjacentList_HH

#include <vector>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <assert.h>
#include <time.h>
#ifndef __WIN32__
#include <sys/resource.h>
#endif


#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace cluster {
namespace calibur {

enum ADJ_LIST_MODE { MATRIX, LIST , LITE };

#ifndef _LARGE_DECOY_SET_
typedef unsigned short LIST_TYPE;
#else  //number of input decoys large than 65535
	typedef unsigned int LIST_TYPE;
#endif

class AdjacentList
{
public:
	static ADJ_LIST_MODE list_mode_;
	int which_;         // index of the decoy this AdjacentList is for
	int num_neighbors_;      // synchronized with the size of neigh
	std::vector<int> neigh; // keep a record of all the neighbors
	std::vector<double> dist;
	std::vector<LIST_TYPE> reverse_index_; // References index of the decoy within the array
	AdjacentList();
	AdjacentList(int which, int num);
	~AdjacentList();
	void add(int n, double d, bool ifNeigh);
#ifdef _ADD_LITE_MODE_
    void add(int num);
#endif
	void remove(int n);
	double getD(int n) const;
};

typedef utility::pointer::shared_ptr< AdjacentList > AdjacentListOP;

}
}
}

#endif
