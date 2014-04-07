// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/cmiles/kcluster.cc
/// @author Christopher Miles (cmiles@uw.edu)

// C/C++ headers
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>

// External headers
#include <boost/format.hpp>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/cmiles.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <numeric/random/random.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>

using core::pose::Pose;
using utility::vector1;

class PairComparator {
 public:
  bool operator() (const std::pair<int, int>& a, const std::pair<int, int>& b) const {
    if (a.first < b.first) {
      return true;
    } else if (a.first > b.first) {
      return false;
    } else {  // a.first == b.first
      return a.second < b.second;
    }
  }
};

/// @detail Computes the distance between two structures using gdtmm
double dist(const vector1<Pose>& models, int i, int j) {
  using std::map;
  using std::pair;

  static map<pair<int, int>, double, PairComparator> cache;

  pair<int, int> key(i, j);
  pair<int, int> reverse_key(j, i);

  double d;
  if (cache.find(key) == cache.end()) {
    d = (i == j) ? 0 : 1 - core::scoring::CA_gdtmm(models[i], models[j]);
    cache[key] = d;
    cache[reverse_key] = d;
  }

  return cache[key];
}

/// @detail Selects k cluster centers using a simple rule: the next cluster
/// center to be added is as far as possible from the centers chosen so far
void choose_centroids(const int k, const vector1<Pose>& models, vector1<int>* centroids) {
  using std::cout;
  using std::endl;

  const int num_models = models.size();
  const int initial_centroid = numeric::random::random_range(1, num_models);

  centroids->push_back(initial_centroid);
  cout << "Added cluster center => " << initial_centroid << endl;

  for (int c = 2; c <= k; ++c) {
    int max_index = 0;
    double max_dist = std::numeric_limits<double>::min();

    for (int i = 1; i <= num_models; ++i) {
      double dist_nearest = std::numeric_limits<double>::max();

      for (vector1<int>::const_iterator j = centroids->begin(); j != centroids->end(); ++j) {
        double distance = dist(models, i, *j);
        if (distance < dist_nearest) {
          dist_nearest = distance;
        }
      }

      if (dist_nearest > max_dist) {
        max_dist = dist_nearest;
        max_index = i;
      }
    }

    centroids->push_back(max_index);
    cout << "Added cluster center => " << max_index << endl;
  }
}

/// @detail Assign each model to the nearest cluster center
void assign_models(const vector1<Pose>& models, const vector1<int>& centroids, vector1<int>* assignments) {
  assignments->resize(models.size());

  for (int i = 1; i <= models.size(); ++i) {
    double min_dist = std::numeric_limits<double>::max();
    int min_index = 0;

    for (int j = 1; j <= centroids.size(); ++j) {
      double distance = dist(models, i, centroids[j]);
      if (distance < min_dist) {
        min_dist = distance;
        min_index = j;
      }
    }

    (*assignments)[i] = min_index;
  }
}

/// @brief Generate per-cluster silent files
void write_output_files(const vector1<Pose>& models, const vector1<int>& centroids, const vector1<int>& assignments) {
using core::io::silent::SilentFileData;
  using core::io::silent::SilentStructFactory;
  using core::io::silent::SilentStructOP;
  using std::string;

  vector1<string> filenames;
  for (int i = 1; i <= centroids.size(); ++i) {
    filenames.push_back(str(boost::format("c.%d.out") % i));
  }

  for (int i = 1; i <= models.size(); ++i) {
    SilentStructOP silent = SilentStructFactory::get_instance()->get_silent_struct_out();
    silent->fill_struct(models[i]);

    SilentFileData sfd;
    sfd.write_silent_struct(*silent, filenames[assignments[i]]);
  }
}

/// @detail Writes per-cluster distance distributions to file
void write_distances(const vector1<Pose>& models, const vector1<int>& centroids, const vector1<int>& assignments) {
  using std::endl;
  using std::ofstream;

  vector1<vector1<double> > distances;
  distances.resize(centroids.size());

  for (int i = 1; i <= models.size(); ++i) {
    int cluster_id = assignments[i];
    int cluster_model = centroids[cluster_id];
    distances[cluster_id].push_back(dist(models, i, cluster_model));
  }

  for (int i = 1; i <= distances.size(); ++i) {  // foreach cluster
    ofstream out(str(boost::format("c.%d.dists") % i).c_str());

    for (int j = 1; j <= distances[i].size(); ++j) {  // foreach cluster member
      out << distances[i][j] << endl;
    }

    out.close();
  }
}

int main(int argc, char* argv[]) {
    try {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::chemical::ChemicalManager;
  using core::chemical::ResidueTypeSetCAP;
  using core::import_pose::pose_stream::MetaPoseInputStream;
  devel::init(argc, argv);

  // TODO(cmiles) use -in:file:residue_type_set parameter instead of centroid
  ResidueTypeSetCAP typeset = ChemicalManager::get_instance()->residue_type_set("centroid");
  MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();

  vector1<Pose> models;

  Pose pose;
  while (input.has_another_pose()) {
    input.fill_pose(pose, *typeset);
    models.push_back(pose);
  }

  vector1<int> centroids;
  choose_centroids(option[OptionKeys::cmiles::kcluster::num_clusters](), models, &centroids);

  vector1<int> assignments;
  assign_models(models, centroids, &assignments);
  write_output_files(models, centroids, assignments);
  write_distances(models, centroids, assignments);
    } catch ( utility::excn::EXCN_Base const & e ) {
                std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
     }
  return 0;
}
