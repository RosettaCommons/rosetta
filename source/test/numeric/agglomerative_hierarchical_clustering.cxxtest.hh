// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/agglomerative_hierarchical_clustering.cxxtest.hh
/// @brief
/// @author Domini Gront

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <numeric/ClusteringTreeNode.hh>
#include <numeric/agglomerative_hierarchical_clustering.hh>

#include <core/types.hh>
#include <utility/vector1.hh>

#include <cmath>
#include <iostream>

using namespace numeric;

class agglomerative_hierarchical_clusteringTest : public CxxTest::TestSuite {
public:

	agglomerative_hierarchical_clusteringTest() {};

	// Shared initialization goes here.
	void setUp() {
	}

	void test_int_clustering() {

	    int data[] = {0,1,2,4,5,3,8,11,6,0};
	    core::Size nData = 9;
	    utility::vector1<core::Size> data_in;
	    for(core::Size i=1;i<=nData;i++)
		data_in.push_back( data[i] );
	    utility::vector1<utility::vector1<Real> > dm;
	    for(Size i=1;i<=nData;i++) {
		utility::vector1<Real> row(nData);
		for(Size j=1;j<=nData;j++) {
		    Real val = abs(data[i] - data[j]);
		    row[j] = val;
		}
		dm.push_back(row);
	    }
	    TS_ASSERT_EQUALS( dm.size(), nData)
	    TS_ASSERT_EQUALS( dm[1].size(),nData)
		
		SingleLinkClusterer slc;
		AverageLinkClusterer alc;
		CompleteLinkClusterer clc;
		
	    //---------- try single link
	    utility::vector1<ClusteringTreeNodeOP> clusters = slc.cluster(dm,4);
	    TS_ASSERT_EQUALS( clusters.size(), 4)
	    for(Size i=1;i<=4;i++) {
		utility::vector1<core::Size> data_out;
		get_cluster_data(data_in,clusters[i],data_out);
		/*
		std::cout<<clusters[i]->id()<<" "<<clusters[i]->size()<<" "<<clusters[i]->distance()<<" : ";
		for(Size j=1;j<=data_out.size();j++)
		    std::cout<<data_out[j]<<" ";
		std::cout<<"\n";
		*/
	    }

	    //---------- try complete link
	    clusters = clc.cluster(dm,4);
	    TS_ASSERT_EQUALS( clusters.size(), 4)
	    for(Size i=1;i<=4;i++) {
		utility::vector1<core::Size> data_out;
		get_cluster_data(data_in,clusters[i],data_out);
		/*
		std::cout<<clusters[i]->id()<<" "<<clusters[i]->size()<<" "<<clusters[i]->distance()<<" : ";
		for(Size j=1;j<=data_out.size();j++)
		    std::cout<<data_out[j]<<" ";
		std::cout<<"\n";
		*/
	    }

	    //---------- try average link
	    clusters = alc.cluster(dm,4);
	    TS_ASSERT_EQUALS( clusters.size(), 4)
	    for(Size i=1;i<=4;i++) {
		utility::vector1<core::Size> data_out;
		get_cluster_data(data_in,clusters[i],data_out);
		/*
		std::cout<<clusters[i]->id()<<" "<<clusters[i]->size()<<" "<<clusters[i]->distance()<<" : ";
		for(Size j=1;j<=data_out.size();j++)
		    std::cout<<data_out[j]<<" ";
		std::cout<<"\n";
		*/
	    }
	}

	// Shared finalization goes here.
	void tearDown() {}
};
