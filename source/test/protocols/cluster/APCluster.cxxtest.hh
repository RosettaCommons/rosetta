// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

#include <core/types.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>
#include <numeric/xyzVector.hh>
#include <protocols/cluster/APCluster.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

#include <algorithm>
// AUTO-REMOVED #include <fstream>

//Auto Headers



static basic::Tracer TR("protocols.cluster.APClusterTest.cxxtest");


using namespace protocols::cluster;
using namespace core;
using utility::vector1;
using namespace utility::tools;


class APClusterTest : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_cluster_random_3D_points() {
		vector1< Vector > pts;
		for(Size i = 1; i <= 1000; ++i) {
			pts.push_back( Vector( 10*numeric::random::rg().uniform(), 10*numeric::random::rg().uniform(), 10*numeric::random::rg().uniform() ) );
		}

		// These get us reasonably close to convergence, but also run pretty fast.
		core::Size const maxits = 200;
		core::Size const convits = 10;
		core::Real const lambda = 0.90;

		{ // first test: dense similarities
			APCluster cluster(pts.size());
			vector1< Real > sims;
			for(Size i = 1; i <= pts.size(); ++i) {
				for(Size k = i+1; k <= pts.size(); ++k) {
					Real const negdist2 = - pts[i].distance_squared( pts[k] );
					cluster.set_sim(i, k, negdist2);
					cluster.set_sim(k, i, negdist2);
					sims.push_back(negdist2);
				}
			}
			std::sort( sims.begin(), sims.end() );
			Real const min_sim = sims[1];
			Real const med_sim = sims[ sims.size() / 2 ];
			TR << "Minimum similarity: " << min_sim << std::endl;
			TR << "Median similarity: " << med_sim << std::endl;
			TR << "Maximum similarity: " << sims[ sims.size() ] << std::endl;
			// Small number of clusters
			for(Size k = 1; k <= pts.size(); ++k) {
				cluster.set_sim(k, k, min_sim);
			}
			cluster.cluster(maxits, convits, lambda);
			write_cluster(pts, cluster, "dense_min.kin");
			// Moderate number of clusters
			for(Size k = 1; k <= pts.size(); ++k) {
				cluster.set_sim(k, k, med_sim);
			}
			cluster.cluster(maxits, convits, lambda);
			write_cluster(pts, cluster, "dense_med.kin");
			// Save to binary and read back in
			cluster.save_binary("dense_med.bin");
			APCluster cluster2(0);
			cluster2.load_binary("dense_med.bin");
			write_cluster(pts, cluster2, "dense_med_reconstr.kin");
		}

		{ // second test: sparse similarities
			APCluster cluster(pts.size(), pts.size()/20);
			vector1< Real > sims;
			for(Size i = 1; i <= pts.size(); ++i) {
				for(Size k = i+1; k <= pts.size(); ++k) {
					Real const negdist2 = - pts[i].distance_squared( pts[k] );
					cluster.set_sim(i, k, negdist2);
					cluster.set_sim(k, i, negdist2);
					sims.push_back(negdist2);
				}
			}
			std::sort( sims.begin(), sims.end() );
			Real const min_sim = sims[1];
			Real const med_sim = sims[ sims.size() / 2 ];
			TR << "Minimum similarity: " << min_sim << std::endl;
			TR << "Median similarity: " << med_sim << std::endl;
			TR << "Maximum similarity: " << sims[ sims.size() ] << std::endl;
			// Small number of clusters
			for(Size k = 1; k <= pts.size(); ++k) {
				cluster.set_sim(k, k, min_sim);
			}
			cluster.cluster(maxits, convits, lambda);
			write_cluster(pts, cluster, "sparse_min.kin");
			// Moderate number of clusters
			for(Size k = 1; k <= pts.size(); ++k) {
				cluster.set_sim(k, k, med_sim);
			}
			cluster.cluster(maxits, convits, lambda);
			write_cluster(pts, cluster, "sparse_med.kin");
			// Save to binary and read back in
			cluster.save_binary("sparse_med.bin");
			APCluster cluster2(0);
			cluster2.load_binary("sparse_med.bin");
			write_cluster(pts, cluster2, "sparse_med_reconstr.kin");
		}

		TR << "FIXME:  APClusterTest exercises the code but doesn't assert anything about the results!" << std::endl;
	}

	void write_cluster(vector1< Vector > const & pts, APCluster const & cluster, std::string const & filename)
	{
		vector1< std::string > colors = make_vector1<std::string>("red", "orange", "gold", "yellow", "lime", "green", "sea", "cyan", "sky", "blue", "purple", "magenta", "hotpink");
		std::ofstream out( filename.c_str() );
		out << "@kinemage {" << filename << "}\n";
		out << "@onewidth\n";
		out << "@balllist {points} radius= 0.1 color= gray\n";
		for(Size i = 1; i <= pts.size(); ++i) {
			out << "{" << i << "} " << pts[i].x() << " " << pts[i].y() << " " << pts[i].z() << "\n";
		}
		out << "@group {clusters}\n";
		vector1< Size > exemplars, members;
		cluster.get_all_exemplars(exemplars);
		TR << exemplars.size() << " clusters for " << cluster.num_pts() << " points" << std::endl;
		for(Size ii = 1; ii <= exemplars.size(); ++ii) {
			std::string color = colors[ (ii % colors.size()) + 1 ];
			out << "@subgroup {cluster " << ii << "} dominant\n";
			Size const i = exemplars[ii];
			cluster.get_cluster_for(i, members);
			out << "@vectorlist {cluster " << ii << "} color= " << color << "\n";
			for(Size jj = 1; jj <= members.size(); ++jj) {
				Size const j = members[jj];
				out << "{ cluster" << ii << " exemplar" << i << " }P " << pts[j].x() << " " << pts[j].y() << " " << pts[j].z() << "\n";
				out << "{ cluster" << ii << " exemplar" << i << " }  " << pts[i].x() << " " << pts[i].y() << " " << pts[i].z() << "\n";
			}
			out << "@balllist {cluster " << ii << "} radius= 0.1 color= " << color << "\n";
			for(Size jj = 1; jj <= members.size(); ++jj) {
				Size const j = members[jj];
				out << "{ cluster" << ii << " exemplar" << i << " } " << pts[j].x() << " " << pts[j].y() << " " << pts[j].z() << "\n";
			}
		}
		out.close();
	}
};

