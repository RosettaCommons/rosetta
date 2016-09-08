// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/brunette/identify_homolog_inaccuracies.cc
///
/// @brief  Reads in an alignment file, ev_map file and native. Generates partial threading models for all templates and outputs contact probabilites.
/// @author TJ Brunette

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <devel/init.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/Tracer.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/scoring/rms_util.hh>
#include <core/kinematics/FoldTree.hh>

//#include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/coord_util.hh>
#include <protocols/comparative_modeling/PartialThreadingMover.hh>

//utilities
#include <utility/io/ozstream.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>

#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <algorithm>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>

using namespace core::sequence;
using utility::vector1;
using core::Size;
using core::Real;
using std::string;
using std::map;
using std::multimap;
using core::pose::Pose;
static THREAD_LOCAL basic::Tracer tr( "contact_map_from_homologs.main" );
///////////////////////////////////////////////////////////////////////////////
/// @detail converts string to upper case
///////////////////////////////////////////////////////////////////////////////
//void string_to_upper(string& str){
//	std::transform(str.begin(),str.end(),str.begin(), &toupper );
//}
///////////////////////////////////////////////////////////////////////////////
/// @detail gets the aligment from file and ranks them.
///////////////////////////////////////////////////////////////////////////////
vector1< SequenceAlignment > generateRankedAlignments(map <string,SequenceAlignment> & alns,Real THRESHOLD_FOR_E_VAL){
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using utility::file::FileName;
	using std::multimap;
  vector1<SequenceAlignment> rankedAlignments;
  multimap<Real,string> rankedTemplate_map; //templates are from ev_map or hh_map
  map<string,Real> template_map;
  vector1<SequenceAlignment>::iterator start_rankedAln,stop_rankedAln;
  map<string,SequenceAlignment>::iterator start_alns, stop_alns;
  Real low_e_val = 9999;
	char c;
  //need to get alignments from file and rank depending on type
  if(option[ cm::ev_map ].user()){
    map<Real,string>::iterator start_rankedTemplate,stop_rankedTemplate;
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
					if(pdbid.size()==6){//takes care of cases with _
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
      map<Real,string>::reverse_iterator start_rankedTemplate,stop_rankedTemplate;
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
	if(rankedAlignments.size() < alns.size()){
    rankedAlignments.clear();
    start_alns = alns.begin();
    stop_alns = alns.end();
    while (start_alns != stop_alns){
			string const aln_id( start_alns->second.sequence(2)->id() );
			string const template_id( aln_id.substr(0,5) );
			rankedAlignments.push_back(start_alns->second);
      start_alns++;
    }
  }
	return rankedAlignments;
}


///////////////////////////////////////////////////////////////////////////////
/// @detail gathers the poses
///////////////////////////////////////////////////////////////////////////////
map< string, Pose > poses_from_cmd_line(utility::vector1< std::string > const & fn_list){

  using utility::file::file_exists;
  using core::import_pose::pose_from_file;
  using namespace core::chemical;

  ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
  map< string, Pose > poses;
  typedef vector1< string >::const_iterator iter;
  for ( iter it = fn_list.begin(), end = fn_list.end(); it != end; ++it ) {
    if ( file_exists(*it) ) {
      Pose pose;
      core::import_pose::pose_from_file( pose, *rsd_set, *it , core::import_pose::PDB_file);
      string name = utility::file_basename( *it );
      name = name.substr( 0, 5 );
      poses[name] = pose;
    }
  }
  return poses;
}
///////////////////////////////////////////////////////////////////////////////
/// @detail calculates contact map without having a mapping.
///////////////////////////////////////////////////////////////////////////////
multimap <Size,Size> calc_contacts(core::pose::Pose pose,SequenceAlignment aln, Real cutoff_dist, Size aa_separation){
	using namespace core::id;
  multimap<Size,Size> contact_map;
	SequenceMapping mapping( aln.sequence_mapping(2,1) );
	for ( Size i = 1; i <= pose.size(); ++i ) {
		for ( Size j = i+1; j <= pose.size(); ++j ) {
			Size const i_native(mapping[i]);
			Size const j_native(mapping[j]);
			if(j_native-i_native >= aa_separation){
				const core::conformation::Residue& resi = pose.residue(i);
				const core::conformation::Residue& resj = pose.residue(j);
				Real distance = resi.xyz("CA").distance( resj.xyz("CA"));
				if (distance < cutoff_dist){
					contact_map.insert(std::pair<Size,Size>(i_native,j_native));
				}
			}
		}
	}
	return(contact_map);
}
///////////////////////////////////////////////////////////////////////////////
/// @detail calculates contact map
///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] ) {
	try {

	using core::Size;
	using core::Real;
	using std::map;
	using std::string;
	using utility::vector1;
	using ObjexxFCL::string_of;
	using core::sequence::read_fasta_file;
	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::sequence;
  using namespace ObjexxFCL::format;
	using protocols::comparative_modeling::gather_coords;
	using namespace protocols::comparative_modeling;
	Real THRESHOLD_FOR_E_VAL = 1e-30;
	Size MIN_POSE_SIZE = 5;
	Real CUTOFF_DISTANCE = 8;
	Size AA_SEPARATION = 12;
	Pose native_pose;
	std::string decoy_aln_id;
	Size decoy_aln_len = 1;
	Real decoy_aln_perc = 1;
	Size maxNumbContacts = 0;
	SequenceAlignment decoy_aln;
	devel::init(argc, argv);
	multimap<Size,Size>::iterator contact_map_start,contact_map_stop;
	core::import_pose::pose_from_file(
					 native_pose,
					 *(rsd_set_from_cmd_line()),
					 option[ in::file::native ]()
					 );
	string query_sequence (
												 read_fasta_file( option[ in::file::fasta ]()[1])[1]->sequence()	);
	//---Get all alignments ------
	vector1< std::string > align_fns = option[ in::file::alignment ]();
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
  map< string, Pose > template_poses = poses_from_cmd_line(option[ in::file::template_pdb ]());


	Size max_pose_len(0);
  vector1< Pose >   poses;
  vector1< string > aln_ids;
	multimap<Size,Size> tmpContact_map;
	vector1 <vector1 < Size > > contactCount_array2D(
																									 query_sequence.size(),
																									 vector1<Size>(query_sequence.size(),0));
	vector1 <vector1 < Size > > topcontactCount_array2D(
																									 query_sequence.size(),
																									 vector1<Size>(query_sequence.size(),0));
	vector1 <vector1 < Real > > contactScaled_array2D(
																									 query_sequence.size(),
																									 vector1<Real>(query_sequence.size(),0));
	for ( Size ii = 1; ii <= rankedAlignments.size(); ++ii ) {
		string aln_id( rankedAlignments[ii].sequence(2)->id() );
		string template_id( aln_id.substr(0,5) );
		map< string, Pose >::iterator pose_it = template_poses.find( template_id );
    if ( pose_it == template_poses.end() ) {
      string msg( "Error: can't find pose (id = "
		  + template_id + ")"
		  );
      //utility_exit_with_message(msg);
      std::cout << msg << std::endl;
      continue;
    } else {
			Pose query_pose;
			core::pose::make_pose_from_sequence(
																					query_pose, query_sequence, *(rsd_set_from_cmd_line())
																					);
			core::pose::Pose template_pose = pose_it->second;


			std::cout << "generating contacts " << std::endl
								<< aln_id << std::endl;
			PartialThreadingMover mover( rankedAlignments[ii], template_pose );
			mover.apply( query_pose );
			max_pose_len = std::max( max_pose_len, query_pose.size() );
			if(query_pose.size() >= MIN_POSE_SIZE){
				tmpContact_map.clear();
				if ( option[ run::debug ]() ) {
					query_pose.dump_pdb( aln_id + ".pdb" );
				}
				SequenceAlignment aln(align_poses_naive(native_pose, query_pose));
				tmpContact_map = calc_contacts(query_pose,aln,CUTOFF_DISTANCE, AA_SEPARATION);
				contact_map_start = tmpContact_map.begin();
				contact_map_stop = tmpContact_map.end();
				while(contact_map_start != contact_map_stop){
					contactCount_array2D[contact_map_start->first][contact_map_start->second]++;
					if(ii==1){
						topcontactCount_array2D[contact_map_start->first][contact_map_start->second]++;}
					if(maxNumbContacts <contactCount_array2D[contact_map_start->first][contact_map_start->second])
						maxNumbContacts = contactCount_array2D[contact_map_start->first][contact_map_start->second];
					contact_map_start++;
				}
			}
		}
	}
	string filename(option[ out::file::o ]());
	utility::io::ozstream output(filename.c_str());

	for(Size jj = 1; jj<=query_sequence.size(); jj++)
		for(Size kk = jj+AA_SEPARATION; kk<=query_sequence.size(); kk++){
			if(contactCount_array2D[jj][kk] != 0){
				contactScaled_array2D[jj][kk] = (Real)contactCount_array2D[jj][kk]/(Real)maxNumbContacts;
			}
			output << jj << " " << kk << " " << 	contactScaled_array2D[jj][kk] << std::endl;
		}
	SequenceAlignment aln(align_poses_naive(native_pose, native_pose));
	multimap<Size,Size> nativeContact_map = calc_contacts(native_pose,aln,CUTOFF_DISTANCE, AA_SEPARATION);
	contact_map_start = nativeContact_map.begin();
	contact_map_stop = nativeContact_map.end();
	while(contact_map_start != contact_map_stop){
		output << contact_map_start->first << " " << contact_map_start->second << std::endl;
		contact_map_start++;
	}
	output.close();
	//---also output the top ranked structure
	filename ="top_" + filename;
	utility::io::ozstream topOutput(filename.c_str());

	for(Size jj = 1; jj<=query_sequence.size(); jj++)
		for(Size kk = jj+AA_SEPARATION; kk<=query_sequence.size(); kk++){
			topOutput << jj << " " << kk << " " << topcontactCount_array2D[jj][kk] << std::endl;
		}
	contact_map_start = nativeContact_map.begin();
	contact_map_stop = nativeContact_map.end();
	while(contact_map_start != contact_map_stop){
		topOutput << contact_map_start->first << " " << contact_map_start->second << std::endl;
		contact_map_start++;
	}
	topOutput.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}


