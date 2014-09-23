// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/brunette/alignmentClustering
///
/// @brief  Divide input alns into clusters based on gdtmm comparison of partial models using ranked alignments.
/// @author TJ Brunettes

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>



//#include <devel/init.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>

#include <basic/Tracer.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/scoring/rms_util.hh>
#include <core/kinematics/FoldTree.hh>

#include <protocols/comparative_modeling/AlignmentClustering.hh>

// AUTO-REMOVED #include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/coord_util.hh>
#include <protocols/comparative_modeling/PartialThreadingMover.hh>

//utilities

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// AUTO-REMOVED #include <ObjexxFCL/format.hh>

#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <algorithm>

#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end



static thread_local basic::Tracer tr( "AlignmentClustering" );

namespace protocols {
namespace comparative_modeling {

using utility::vector1;
using core::Size;
using core::Real;
using std::string;
using std::map;
using std::set;
using std::multimap;
using core::pose::Pose;

using namespace core::sequence;
///////////////////////////////////////////////////////////////////////////////
/// @detail Creates a AlignmentCluster
///////////////////////////////////////////////////////////////////////////////
AlignmentCluster::AlignmentCluster(SequenceAlignment & aln_in){
  alns.push_back(aln_in);
}
///////////////////////////////////////////////////////////////////////////////
/// @detail Deletes an AlignmentCluster object
///////////////////////////////////////////////////////////////////////////////
AlignmentCluster::~AlignmentCluster(){
  alns.clear();
}
///////////////////////////////////////////////////////////////////////////////
/// @detail add alignment
///////////////////////////////////////////////////////////////////////////////
void AlignmentCluster::add_aln(SequenceAlignment & aln_in){
  alns.push_back(aln_in);
}
///////////////////////////////////////////////////////////////////////////////
/// @detail get aln
///////////////////////////////////////////////////////////////////////////////
SequenceAlignment AlignmentCluster::get_aln(Size index){
  return (alns[index]);
}
///////////////////////////////////////////////////////////////////////////////
/// @detail get aln
///////////////////////////////////////////////////////////////////////////////
Real AlignmentCluster::size(){
  return alns.size();
}
///////////////////////////////////////////////////////////////////////////////
/// @detail get cluster center -- assumes lowest energy point
///////////////////////////////////////////////////////////////////////////////
//SequenceAlignment AlignmentCluster::get_clusterCenter(){
SequenceAlignment AlignmentCluster::get_clusterCenter(){
  return(alns[1]);
}
///////////////////////////////////////////////////////////////////////////////
/// @detail outputs cluster to file
///////////////////////////////////////////////////////////////////////////////
void AlignmentCluster::output(std::ostream & alignment_out){
  for ( Size ii = 1; ii <= alns.size(); ++ii ) {
    alns[ii].printGrishinFormat(alignment_out);
    string const aln_id( alns[ii].sequence(2)->id() );
    tr << aln_id << "," ;
  }
  tr << std::endl;
}
///////////////////////////////////////////////////////////////////////////////
/// @detail merges two clusters making sure the alignments are unique
///////////////////////////////////////////////////////////////////////////////
void AlignmentCluster::merge(AlignmentClusterOP cluster_in){
  set<string> alignmentsInCluster;
  set<string>::iterator location;
  for(Size ii=1; ii<= alns.size(); ++ii){
    alignmentsInCluster.insert(alns[ii].sequence(2)->id());
  }
  for(Size jj= 1; jj <= cluster_in->size(); ++jj){
    SequenceAlignment tempAln = cluster_in->get_aln(jj);
    location = alignmentsInCluster.find(tempAln.sequence(2)->id());
    if (location == alignmentsInCluster.end()){
      add_aln(tempAln);
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
/// @detail gets the percent of alignments in the input cluster that are in the current cluster
///////////////////////////////////////////////////////////////////////////////
Real AlignmentCluster::overlap(AlignmentClusterOP cluster_in){
  set<string> alignmentsInCluster;
  set<string>::iterator location;
  for(Size ii=1; ii<= alns.size(); ++ii)
    alignmentsInCluster.insert(alns[ii].sequence(2)->id());
  //loop through and see what percent of cluster_in is in the current cluster;
  Size count_inClust = 0;
  for(Size ii=1; ii<=cluster_in->size(); ++ii){
    location = alignmentsInCluster.find(cluster_in->get_aln(ii).sequence(2)->id());
    if (location != alignmentsInCluster.end())
      count_inClust++;
  }
  return((Real)count_inClust/(Real)cluster_in->size());
}
///////////////////////////////////////////////////////////////////////////////
/// @detail creates an alignment clustering object, this is in charge of
///    clustering
///////////////////////////////////////////////////////////////////////////////
AlignmentClustering::AlignmentClustering(){
  using core::pose::Pose;
  using basic::options::option;

  using utility::file::FileName;
  using protocols::comparative_modeling::gather_coords;
  using ObjexxFCL::string_of;

  using namespace protocols;
  using namespace core::chemical;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace core::sequence;
  using namespace ObjexxFCL::format;
  Real THRESHOLD_FOR_E_VAL = 1e-30;
  Real MIN_GDT = .70;
  Real MAX_GDT = .90;
  Real INCREMENT_GDT = .05;
  Size GOAL_NUMB_CLUSTERS = 5;
  Real MAX_CLUSTER_OVERLAP = .70;
  Size MIN_POSE_SIZE = 5;

  vector1< std::string > align_fns = option[ in::file::alignment ]();
  //vector1< SequenceAlignment > alns;
  map<string,SequenceAlignment> alns;
  map<string,SequenceAlignment>::iterator location_alns;
  for ( Size ii = 1; ii <= align_fns.size(); ++ii ) {
    vector1< SequenceAlignment > tmp_alns = core::sequence::read_aln(
		option[ cm::aln_format ](), align_fns[ii]
		);
    for ( Size jj = 1; jj <= tmp_alns.size(); ++jj ) {
      string aln_id = tmp_alns[jj].sequence(2)->id();
      location_alns = alns.find(aln_id);
      while(location_alns != alns.end()){//if alignment exists change the name by adding 100000 to the name.
	string fixed_aln_id = aln_id + "X";
	tmp_alns[jj].sequence(2)->id(fixed_aln_id);
	aln_id = tmp_alns[jj].sequence(2)->id();
	location_alns = alns.find(aln_id);
      }
      aln_id = (tmp_alns[jj].sequence(2)->id() );
      alns.insert(std::pair<string,SequenceAlignment>(aln_id,tmp_alns[jj]));
    }
  }
  //rank alignments
  vector1<SequenceAlignment> rankedAlignments = generateRankedAlignments(alns,THRESHOLD_FOR_E_VAL);
  vector1<SequenceAlignment> rankedAlignments_valid;
  //generate poses
  map< string, Pose > template_poses = poses_from_cmd_line(
		option[ in::file::template_pdb ]());
  using core::sequence::read_fasta_file;
  string query_sequence (
	read_fasta_file( option[ in::file::fasta ]()[1])[1]->sequence()	);
  Size max_pose_len(0);
  vector1< Pose >   poses;
  vector1< string > aln_ids;
  for ( Size ii = 1; ii <= rankedAlignments.size(); ++ii ) {
    string const aln_id( rankedAlignments[ii].sequence(2)->id() );
    string const template_id( aln_id.substr(0,5) );
    map< string, Pose >::iterator pose_it = template_poses.find( template_id );
    if ( pose_it == template_poses.end() ) {
      string msg( "Error: can't find pose (id = "
		  + template_id + ")"
		  );
      //utility_exit_with_message(msg);
      std::cout << msg << std::endl;
      continue;
    } else {
      core::chemical::ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );
      core::pose::Pose query_pose;
      core::pose::make_pose_from_sequence(
			      query_pose, query_sequence, *(rsd_set)
			      );
      core::pose::Pose template_pose = pose_it->second;
      using namespace protocols::comparative_modeling;
      std::cout << "building incomplete model with " << std::endl
		<< aln_id << std::endl;
      PartialThreadingMover mover( rankedAlignments[ii], template_pose );
      mover.apply( query_pose );
      max_pose_len = std::max( max_pose_len, query_pose.total_residue() );
      if(query_pose.total_residue() >= MIN_POSE_SIZE){
	poses.push_back( query_pose );
	aln_ids.push_back( aln_id );
	rankedAlignments_valid.push_back(rankedAlignments[ii]);
      }
      else
	tr << aln_id << "has only "<<query_pose.total_residue() << "which is below size threshold of " << MIN_POSE_SIZE << std::endl;
      if ( option[ run::debug ]() ) {
	query_pose.dump_pdb( aln_id + ".pdb" );
      }
    }
  } // for rankedAlignmentsalns
  Size total_comparisons_done(0);
  vector1< vector1< Real > > gdtmms(
				    rankedAlignments_valid.size(), vector1< Real >(rankedAlignments_valid.size(), 0.0 )
				    );
  for ( Size ii = 1; ii <= rankedAlignments_valid.size(); ++ii ) {
    gdtmms[ii][ii] = 1.0;
    for ( Size jj = ii + 1; jj <= rankedAlignments_valid.size(); ++jj ) {
      ++total_comparisons_done;
      SequenceAlignment aln( align_poses_naive( poses[ii], poses[jj] ) );
      int n_atoms;
      ObjexxFCL::FArray2D< Real > p1a, p2a;
      protocols::comparative_modeling::gather_coords(poses[ii], poses[jj],
						     aln,
						     n_atoms, p1a, p2a
						     );
      Real const coverage( (Real) n_atoms / (Real) max_pose_len );
      //Real const gdtmm( xyz_gdtmm( p1a, p2a ) );
      using core::scoring::xyz_gdtmm;
      Real const gdtmm( xyz_gdtmm( p1a, p2a ) );
      //Real const gdtmm_adj( coverage * gdtmm );
      //gdtmms[ii][jj] = gdtmm_adj;
      //gdtmms[jj][ii] = gdtmm_adj;
      gdtmms[ii][jj] = gdtmm;
      gdtmms[jj][ii] = gdtmm;
      if ( option[ run::debug ]() ) {
	std::cerr << "coverage = " << coverage << ", n_atoms = " << n_atoms
		  << ", max_pose_len = " << max_pose_len << ", gdtmm = " << gdtmm << std::endl;
	std::cerr << "pose1 seq = " << poses[ii].sequence() << std::endl;
	std::cerr << "pose2 seq = " << poses[jj].sequence() << std::endl;
	std::cerr << aln << std::endl;
	std::cerr << "sim(" << aln_ids[ii] << "," << aln_ids[jj] << ") = " << gdtmm << std::endl;
	std::cerr << "--------------------------------------------------------------------------------"
		  << std::endl;
      }
      if ( total_comparisons_done % 50 == 0 ) {
	std::cout << "." << std::flush;
	if ( total_comparisons_done % 1000 == 0 ) {
	  std::cout << " finished with "
		    << total_comparisons_done << "."	<< std::endl;
	}
      }
    } // jj
  } // ii
  Size number_clusters = 999999;
  vector1<AlignmentClusterOP> cluster_v;
  Real threshold_gdt = MAX_GDT;
  set<Size> mergeSet;
  set<Size>::iterator location_mergeSet;
  multimap<Size,Size> mergeMap;
  multimap<Size,Size>::iterator start_mergeMap,stop_mergeMap;
  while((number_clusters > GOAL_NUMB_CLUSTERS) && (threshold_gdt >= MIN_GDT)){
    cluster_v.clear();
    cluster_v = cluster(gdtmms,rankedAlignments_valid,threshold_gdt);
    //===Merging clusters that have overlap above MAX_CLUSTER_OVERLAP
    mergeSet.clear();
    mergeMap.clear();
    //If cluster has been merged than ignore it.
    for (Size ii = 1; ii <= cluster_v.size(); ++ii){
      location_mergeSet = mergeSet.find(ii);
      if(location_mergeSet == mergeSet.end()){
	for(Size jj=ii+1; jj<=cluster_v.size(); ++jj){
	  location_mergeSet = mergeSet.find(jj);
	  if(location_mergeSet == mergeSet.end()){
	    if (cluster_v[ii]->overlap(cluster_v[jj])>MAX_CLUSTER_OVERLAP){
	      mergeSet.insert(jj);
	      mergeMap.insert(std::pair<int,int>(ii,jj));
	    }
	  }
	}
      }
    }
    start_mergeMap = mergeMap.begin();
    stop_mergeMap = mergeMap.end();
    while(start_mergeMap != stop_mergeMap){
      cluster_v[start_mergeMap->first]->merge(cluster_v[start_mergeMap->second]);
      start_mergeMap++;
    }
    //---end merge
    std::cout << "cluster_v.size"  << cluster_v.size() << "merged clusters" << mergeSet.size() << std::endl;

    number_clusters = cluster_v.size() - mergeSet.size();
    tr << "threshold_gdt" << threshold_gdt << "number_clusters" << number_clusters << std::endl;

    if(number_clusters > GOAL_NUMB_CLUSTERS)
      threshold_gdt = threshold_gdt - INCREMENT_GDT;
  }

  //Output final clusters
  Size numbOutput = 1;
  for (Size ii = 1; ii <= cluster_v.size(); ++ii){
    location_mergeSet = mergeSet.find(ii);
    if(location_mergeSet == mergeSet.end()){
      std::cout << "--------" << std::endl;
      std::stringstream convert_to_string;
      convert_to_string << numbOutput;
      std::string filename = "alignmentCluster_" + convert_to_string.str() + ".filt";
      utility::io::ozstream alignment_out(filename );
      cluster_v[ii]->output(alignment_out);
      numbOutput++;
      alignment_out.close();
    }
  }
  //output results of cluster. For now just output to file.
}

///////////////////////////////////////////////////////////////////////////////
/// @detail Deletes AlignmentClustering object
///////////////////////////////////////////////////////////////////////////////
AlignmentClustering::~AlignmentClustering(){
}
///////////////////////////////////////////////////////////////////////////////
/// @detail  Does the clustering
///////////////////////////////////////////////////////////////////////////////
vector1<AlignmentClusterOP> AlignmentClustering::cluster(vector1< vector1< Real > > & gdtmms, vector1<SequenceAlignment> & rankedAlignments, Real threshold_gdt){
  vector1<AlignmentClusterOP> cluster_v;
  vector1<bool> alignmentInCluster;
  //the first item is automatically in the first cluster.
  for ( Size ii = 1; ii <= rankedAlignments.size(); ++ii ) {
    alignmentInCluster.push_back(false);
  }
  bool allClustered = false;
  Size clusterCenter = 1;
  //This is written so alignments may end up in multiple clusters.
  while(!allClustered){
    AlignmentClusterOP tempCluster( new AlignmentCluster(rankedAlignments[clusterCenter]) );
    alignmentInCluster[clusterCenter] = true;
    for( Size ii = 1; ii <= rankedAlignments.size(); ++ii) {
      if(ii != clusterCenter){
	//string const aln_id(rankedAlignments[ii].sequence(2)->id() );
	//string const aln_center_id(rankedAlignments[clusterCenter].sequence(2)->id());
	//tr << "distance from cluster center" << gdtmms[clusterCenter][ii] << " aln_id:" << aln_id << "cluster_center" << aln_center_id << std::endl;
	if (gdtmms[clusterCenter][ii] >= threshold_gdt){
	  tempCluster->add_aln(rankedAlignments[ii]);
	  alignmentInCluster[ii] = true;
	}
      }
    }
    cluster_v.push_back(tempCluster);
    int nxtClusterCenter = -1;
    for(Size jj = 1; ((jj <= alignmentInCluster.size()) && (nxtClusterCenter==-1))  ; ++jj){
      if(alignmentInCluster[jj] == false){
	nxtClusterCenter = jj;
      }
    }
    if(nxtClusterCenter == -1)
      allClustered = true;
    else
      clusterCenter = nxtClusterCenter;
  }
  return(cluster_v);
}
///////////////////////////////////////////////////////////////////////////////
/// @detail gathers the poses
///////////////////////////////////////////////////////////////////////////////
map< string, Pose > AlignmentClustering::poses_from_cmd_line(utility::vector1< std::string > const & fn_list){

  using utility::file::file_exists;
  using core::import_pose::pose_from_pdb;
  using namespace core::chemical;

  ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );
  map< string, Pose > poses;
  typedef vector1< string >::const_iterator iter;
  for ( iter it = fn_list.begin(), end = fn_list.end(); it != end; ++it ) {
    if ( file_exists(*it) ) {
      Pose pose;
      core::import_pose::pose_from_pdb( pose, *rsd_set, *it );
      string name = utility::file_basename( *it );
      name = name.substr( 0, 5 );
      poses[name] = pose;
    }
  }
  return poses;
}
///////////////////////////////////////////////////////////////////////////////
/// @detail gets the aligment from file and ranks them.  Some models do not have e-vals in ev_map because thy  are too far down in ranking. As long as 10 models are in the hh_map or ev_map the ranking comes from there.
///////////////////////////////////////////////////////////////////////////////
  vector1< SequenceAlignment > AlignmentClustering::generateRankedAlignments(map <string,SequenceAlignment> & alns,Real THRESHOLD_FOR_E_VAL){
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using utility::file::FileName;
  vector1<SequenceAlignment> rankedAlignments;
  multimap<Real,string> rankedTemplate_map; //templates are from ev_map or hh_map
  map<string,Real> template_map;
  vector1<SequenceAlignment>::iterator start_rankedAln,stop_rankedAln;
  map<string,SequenceAlignment>::iterator start_alns, stop_alns;
  //mjo commenting out 'hhsearch_ranking' because it is unused and causes a warning
  //bool hhsearch_ranking = 0;
  //mjo commenting out 'evmap_ranking' because it is unused and causes a warning
  //bool evmap_ranking = 0;
  Real low_e_val = 9999;
  Size minNumbRankedModels = 10;
  //need to get alignments from file and rank depending on type
  if(option[ cm::ev_map ].user()){
    multimap<Real,string>::iterator start_rankedTemplate,stop_rankedTemplate;
    map<string,Real>::iterator start_template, stop_template,location_template;
    vector1< FileName > ev_mapFile( option[ cm::ev_map ]() );
    if (file_exists(ev_mapFile[1])){
      tr << "ev_map being read in" << std::endl;
      std::string const & filename = ev_mapFile[1];
      utility::io::izstream data(filename);
      if(!data){
	utility_exit_with_message(" Warning: can't open file" + filename + "!");
      }
      string line, pdbid;
      Real e_val;
      while (getline(data,line)){
	std::istringstream line_stream(line);
	line_stream >> pdbid >> e_val;
	if(pdbid != "template"){
	  if(pdbid.size()==6){
	    pdbid = pdbid.substr(0,4) + pdbid.substr(5,1);
	  }
	  if(pdbid.size()==7){
	    string temp_string = "X";
	    temp_string[0] = toupper(pdbid[5]);
	    pdbid = pdbid.substr(1,4) + temp_string;
	  }
	  if(e_val < low_e_val)
	    low_e_val = e_val;

	  location_template = template_map.find(pdbid);
	  if (location_template == template_map.end())
	    template_map.insert(std::pair<string,Real>(pdbid,e_val));
	  else{
	    if (location_template->second > e_val)
	      location_template->second = e_val;
	  }
	}
      }
      start_template =template_map.begin();
      stop_template = template_map.end();
      while(start_template!=stop_template){
	rankedTemplate_map.insert(std::pair<Real,string>(start_template->second,start_template->first));
	start_template++;
      }
      if(low_e_val <= THRESHOLD_FOR_E_VAL){
	start_rankedTemplate = rankedTemplate_map.begin();
	stop_rankedTemplate = rankedTemplate_map.end();
	while(start_rankedTemplate != stop_rankedTemplate){
	  start_alns = alns.begin();
	  stop_alns = alns.end();
	  while (start_alns != stop_alns){
	    string aln_id =  start_alns->second.sequence(2)->id();
	    string template_id = aln_id.substr(0,5);
	    if(template_id == start_rankedTemplate->second){
	      rankedAlignments.push_back(start_alns->second);
	    }
	    start_alns++;
	  }
	  start_rankedTemplate++;
	}
      }
    }
  }
  //get hhsearch ranking
  if((option[ cm::hh_map ].user()&& !option[ cm::ev_map ].user()) || (option[ cm::hh_map ].user() && (low_e_val >= THRESHOLD_FOR_E_VAL))){
    rankedTemplate_map.clear();
    rankedAlignments.clear();
    vector1< FileName > hh_mapFile( option[ cm::hh_map ]() );
    if (file_exists(hh_mapFile[1])){
      multimap<Real,string>::reverse_iterator start_rankedTemplate,stop_rankedTemplate;
      tr << "hh_map being read in" << std::endl;
      std::string const & filename = hh_mapFile[1];
      utility::io::izstream data(filename);
      if(!data){
	utility_exit_with_message(" Warning: can't open file" + filename + "!");
      }
      string line, pdbid;
      Real hh_val;
      while (getline(data,line)){
	std::istringstream line_stream(line);
	line_stream >> pdbid >> hh_val;
	rankedTemplate_map.insert(std::pair<Real,string>(hh_val,pdbid));
      }
      start_rankedTemplate =rankedTemplate_map.rbegin();
      stop_rankedTemplate =rankedTemplate_map.rend();
      while(start_rankedTemplate  != stop_rankedTemplate){
	start_alns = alns.begin();
	stop_alns = alns.end();
	while (start_alns != stop_alns){
	  string const aln_id( start_alns->second.sequence(2)->id() );
	  string const template_id( aln_id.substr(0,5) );
	  if(template_id == start_rankedTemplate->second){
	    rankedAlignments.push_back(start_alns->second);
	  }
	  start_alns++;
	}
	start_rankedTemplate++;
      }
    }
  }
  //If the alignments were not ranked in previous two steps.
  if((rankedAlignments.size() < alns.size())&&(rankedAlignments.size()<minNumbRankedModels)){
    rankedAlignments.clear();
    start_alns = alns.begin();
    stop_alns = alns.end();
    while (start_alns != stop_alns){
      rankedAlignments.push_back(start_alns->second);
      start_alns++;
    }
  }
  return rankedAlignments;
  }
} // comparative_modeling
} // protocols
