// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/kalngyk/InitCluster.cc
/// @author SC Li & YK Ng (kalngyk@gmail.com)


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

//#include <protocols/cluster/calibur/SimpPDB.hh>
//#include <protocols/cluster/calibur/rmsd.hh>
#include <protocols/cluster/calibur/AdjacentList.hh>

namespace protocols {
namespace cluster {
namespace calibur {

#define _OVER_RMSD_ 9999999

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

//============================= Adjacent List =================================
// We have two representation if AdjacentList;
// one is list, the other is matrix (arrays), the latter is used when the input
// size is small
// Each AdjacentList contains all the neighbors of a decoy.
// It is similar to the "neighbors" variable in cluster_info_silent.c.

// for LIST mode
AdjacentList::AdjacentList()
{
	which_ = 0;
	//reverse_index_ = nullptr; // not used
	num_neighbors_ = 0;
}

// for MATRIX mode
AdjacentList::AdjacentList(int which, int size)
{
	which_ = which;
	reverse_index_.resize( size, 0 ); //new LIST_TYPE [size]; // mPDB[i]'s index in neigh/dist
	dist.resize( size );
	neigh.resize( size ); // not 100% certain this is the RIGHT size...
	//memset(reverse_index_, 0, size*sizeof(LIST_TYPE));
	for ( int i=0; i < size; i++ ) dist[i] = _OVER_RMSD_;
	num_neighbors_ = 0;
}

AdjacentList::~AdjacentList()
{
	//if (reverse_index_)  delete [] reverse_index_;
}

double
AdjacentList::getD(int which) const
{
	if ( list_mode_ == LIST
#ifdef _ADD_LITE_MODE_
		   || list_mode_ == LITE
#endif
			)
	return _OVER_RMSD_;
	return dist[which];
}

void
AdjacentList::add(int n, double d, bool ifNeigh)
{
#ifdef _ADD_LITE_MODE_
	if (list_mode_ == LITE) return;
#endif
	if ( list_mode_ != LIST ) {
		dist[n] = d;
	}
	if ( ifNeigh ) {
		if ( list_mode_ == LIST ) {
			neigh.push_back(n);
			dist.push_back(d);
		} else {
			neigh[num_neighbors_] = n;
			reverse_index_[n] = num_neighbors_;
		}
		num_neighbors_++;
	}
}

// Used only in LITE mode
#ifdef _ADD_LITE_MODE_
void
AdjacentList::add(int num)
{
	if (list_mode_ == LITE) num_neighbors_ += num;
}
#endif

void
AdjacentList::remove(int n)
{
#ifdef _ADD_LITE_MODE_
	if (list_mode_ == LITE) return;
#endif
	//if (num_neighbors_ < 1)
	// return;
	if ( list_mode_ == LIST ) {
		int i;
		for ( i=0; i < num_neighbors_; i++ ) { // this is slow
			if ( neigh[i] == n ) break;
		}
		if ( i >= num_neighbors_ ) { // not found! should return error
			return;
		}

		// move the last entry to this one
		neigh[i] = neigh[num_neighbors_-1];
		neigh.pop_back();
		dist[i] = dist[num_neighbors_-1];
		dist.pop_back();
		num_neighbors_--;
	} else {
		// we assume there is at least one element in each AdjacentList: itself
		LIST_TYPE index = reverse_index_[n]; // danger! could be out of bound
		std::cout<<"Removing "<<index<<"-th ("<<n<<") in "<<which_<<"'s LIST"<<std::endl;
		int last = neigh[num_neighbors_-1];
		std::cout<<" swap with " << (num_neighbors_-1) << " with id " << last << std::endl;
		neigh[index] = last;
		reverse_index_[last] = index;
		num_neighbors_--;
		if ( num_neighbors_ < 0 ) {
			std::cout << "OMG" << std::endl;
			exit(0);
		}
	}
}
//============================ Adjacent List ==================================

}
}
}

