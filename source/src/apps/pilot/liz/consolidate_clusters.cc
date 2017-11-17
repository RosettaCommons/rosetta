
#include <protocols/moves/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>
#include <iostream>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <devel/init.hh>
#include <devel/init.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <fstream>
#include <ObjexxFCL/FArray2D.hh>
#include <protocols/toolbox/DecoySetEvaluation.hh>
#include <protocols/toolbox/superimpose.hh>
#include <numeric/util.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>

#include <utility/file/FileName.hh>
#include <utility/excn/Exceptions.hh>

/**
   for each cluster, assign a occupancy and a radius.
   if the cluster has lower than the pre-defined population occupancy,
   then find the closest cluster in rmsd (which is above the population occupancy threshold)
   and then assign it to that cluster.
   at the end, you should have a set of clusters with occupancies and radii,
   then can monitor statistics associated with those sets of clusters.
**/

using namespace core;
using namespace std;
using namespace basic::options;
using namespace ObjexxFCL;

OPT_1GRP_KEY(String, analysis, occupancy_data)
OPT_1GRP_KEY(RealVector, analysis, min_occupancies)
OPT_1GRP_KEY(Real, analysis, transition_threshold)
OPT_1GRP_KEY(String, analysis, output_file)
OPT_1GRP_KEY(String, analysis, output_map)

struct Cluster {
public:
  Cluster( FArray2D<double> & coords, core::Real & transition_threshold, core::Size & occupancy, std::string&  tag ) {
    using namespace core::pose::datacache;
    coords_ = coords;
    radius_ = transition_threshold;
    occupancy_ = occupancy;
    //basic::datacache::CacheableString const dc_str = static_cast< basic::datacache::CacheableString const & >( (pose.data()).get( CacheableDataType::JOBDIST_OUTPUT_TAG ) );
    //tag_ = dc_str.str();
    tag_ = tag;
  }

  //core::pose::Pose pose_;
  FArray2D<double> coords_;
  core::Real radius_;
  core::Size occupancy_;
  std::string tag_;
};

struct SetOfClusters {
public:
  utility::vector1< core::Size > occupancy_;
  utility::vector1< core::Real > radius_;
  protocols::moves::mc_convergence_checks::Pool_RMSD_OP structures_;
  utility::vector1< std::string > tags_;


  SetOfClusters() {
    structures_ = new protocols::moves::mc_convergence_checks::Pool_RMSD();
    structures_->clear();
  }

  void clear() {
    occupancy_.clear();
    radius_.clear();
    structures_->clear();
    tags_.clear();
  }

  Cluster get( core::Size ii ) {
    FArray2D<double> coords;
    structures_->get( ii, coords );
    return Cluster( coords, radius_[ ii ], occupancy_[ ii ], tags_[ ii ] );
  }

  core::Size size() {
    runtime_assert( occupancy_.size() == radius_.size() &&
		    radius_.size() == tags_.size() );
    return radius_.size();
  }

  void add( Cluster cluster ){
    structures_->add( cluster.coords_, cluster.coords_.u2(), cluster.tag_ );
    occupancy_.push_back( cluster.occupancy_ );
    radius_.push_back( cluster.radius_ );
    tags_.push_back( cluster.tag_ );
  }

  void consolidate( Cluster cluster ) {
    // find closest cluster
    std::string closest_cluster;
    core::Real best_rmsd;
    core::Size ndx = structures_->evaluate( cluster.coords_,  closest_cluster, best_rmsd );
    std::cout << "ndx of best-match for " << cluster.tag_ << " is " << closest_cluster << " with rms: " << best_rmsd << std::endl;
    std::cout << "adding cluster center " << cluster.tag_ << " to " <<  tags_[ ndx ] << std::endl;
    //for now, must consolidate with closest cluster
    occupancy_[ ndx ] += cluster.occupancy_;
    bool use_weighted_mean = true;
    bool use_max_dist = false;

    if( use_weighted_mean ) {
      radius_[ ndx ] = ((radius_[ ndx ] * occupancy_[ndx]) + (( cluster.radius_ + best_rmsd) *cluster.occupancy_)) /
	( occupancy_[ndx] + cluster.occupancy_ );
    } else if( use_max_dist ) {
      if( best_rmsd > radius_[ ndx ] ) {
	radius_[ ndx ] = best_rmsd;
      }
    }

    /**
    if( rmsd > radius_[ ndx ] ) {
      occupancy_[ ndx ] += cluster.occupancy_;
      radius_[ ndx ] += (radius_[ ndx ] + cluster.radius_ - rmsd)/2;
    } else {
      structures_.add( cluster.pose_ );
      occupancy_.push_back( cluster.occupancy_ );
      radius_.push_back( cluster.radius_ );
    }
    **/
  }


};

void read_occupancy_data( SetOfClusters & set, std::string & occupancy_data, std::string const & silentfile, core::Real & transition_threshold ){
  using namespace basic::options;
  ifstream odata;
  odata.open( occupancy_data.c_str() );

  io::silent::SilentFileData sfd;
  sfd.read_file( silentfile );
  core::pose::Pose pose;

  while( ! odata.eof() ) {
    core::Size occupancy = -1;
    std::string silent_string = "empty";
    odata >>  silent_string >> occupancy;
    if( silent_string.compare("empty") != 0 ) {
      std::cout << "read in silent_string: " << silent_string << " and occupancy: " << occupancy << std::endl;
      std::cout << "looking for decoy with tag: " << silent_string << std::endl;

      core::io::silent::SilentStructOP ss = sfd[ silent_string ];

      ss->fill_pose( pose );
      FArray2D<double> coords(3, pose.size(), 0.0 );
      protocols::toolbox::fill_CA_coords( pose, coords );
      std::string decoy_tag = ss->decoy_tag();
      core::Real threshold = (core::Real)(transition_threshold);

      set.add( Cluster( coords,  threshold, occupancy, decoy_tag ) );
    }
  }
}


int main( int argc, char *argv[] ) {
    try {

  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace utility::file;


  //  devel::init( argc, argv );

  NEW_OPT( analysis::occupancy_data, "file which  contains the cluster-centers and their occupancies", "empty");
  NEW_OPT( analysis::min_occupancies, "list of minimum occupancies to try", utility::vector1<core::Real> (1,0));
  NEW_OPT( analysis::transition_threshold, "transition threshold to assume for all clusters", 1.5 );
  NEW_OPT( analysis::output_file, "output file for statistics", "output.txt" );
  NEW_OPT( analysis::output_map, "output map for cluster-assignments", "output.map");


  devel::init( argc, argv );

  core::Real transition_threshold = option[ analysis::transition_threshold ]();
  SetOfClusters dataset;
  FileName silentfilename = option[ in::file::silent ]()[ 1 ];
  std::string const silentfilestring = silentfilename.name();
  std::string occupancy_data = option[ analysis::occupancy_data ]();
  read_occupancy_data( dataset, occupancy_data, silentfilestring /* cluster centers */, transition_threshold);

  SetOfClusters below_min_occupancy;
  SetOfClusters above_min_occupancy;

  utility::vector1< core::Real > min_occupancy = option[ analysis::min_occupancies ]();

  ofstream outdata;
  std::string output_file = option[ analysis::output_file ]();
  outdata.open( output_file.c_str() );

  //ofstream omap;
  //std::string output_file = option[ analysis::output_map ]();
  //outdata.open( output_file.c_str() );


  for( int jj = 1; jj <= min_occupancy.size(); jj++ ) {
    below_min_occupancy.clear();
    above_min_occupancy.clear();


    for( int ii = 1; ii <= dataset.size(); ii++ ) {
      if( dataset.get( ii ).occupancy_ < min_occupancy[ 1 ] ) {
	below_min_occupancy.add( dataset.get( ii ) );
      } else {
	above_min_occupancy.add( dataset.get( ii ) );
      }
    }

    for( int ii = 1; ii <= below_min_occupancy.size(); ii++ ) {
      Cluster tmp = below_min_occupancy.get(ii);
      //std::cout << "consolidating ii " << ii << " with occupancy " << tmp.occupancy_ << " radius " << tmp.radius_ << " and tag: " << tmp.tag_   << std::endl;
      above_min_occupancy.consolidate( below_min_occupancy.get( ii ) );
    }


    for( int ii = 1; ii <= above_min_occupancy.size(); ii++ ) {
      Cluster clu =  above_min_occupancy.get( ii );
      std::cout << clu.occupancy_ << " " << clu.radius_ << " " << clu.tag_ << std::endl;
    }


    // need to output:
    //   min radius, median/mean radius, max radius
    //   min occupancy, median/mean occupancy, max occupancy
    //   also output: tags of cluster-centers left after consolidation, their occupancies, and their radii
    core::Real min_radius = 10000000000;
    core::Real max_radius = 0;
    core::Size min_occ = 1000000000000;
    core::Size max_occ = 0;

    for( core::Size kk = 1; kk <= above_min_occupancy.size(); kk++ ) {
      if( above_min_occupancy.occupancy_[ kk ] <= min_occ ) {
	min_occ = above_min_occupancy.occupancy_[ kk ];
      }
      if( above_min_occupancy.occupancy_[ kk ] >= max_occ ) {
	max_occ = above_min_occupancy.occupancy_[ kk ];
      }
      if( above_min_occupancy.radius_[ kk ] <= min_radius ) {
	min_radius = above_min_occupancy.radius_[ kk ];
      }
      if( above_min_occupancy.radius_[ kk ] >= max_radius ) {
	max_radius =  above_min_occupancy.radius_[ kk ];
      }
    }

    outdata << " " << min_occupancy[ jj ]
            << " " << above_min_occupancy.size() //number of clusters left
	    << " " << min_occ
	    << " " << numeric::median(above_min_occupancy.occupancy_)
	    << " " <<  max_occ
	    << " " <<  min_radius
	    << " " <<  numeric::median( above_min_occupancy.radius_ )
	    << " " << max_radius << std::endl;

  }

  outdata.close();

    } catch (utility::excn::Exception const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}


