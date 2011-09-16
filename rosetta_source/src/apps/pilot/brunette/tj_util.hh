#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>

#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/conformation/Residue.functions.hh>
#include <core/conformation/Residue.fwd.hh>

#include <core/id/AtomID_Map.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/optimizeH.hh>
#include <core/pack/pack_missing_sidechains.hh>

#include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/coord_util.hh>
#include <protocols/comparative_modeling/PartialThreadingMover.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <map>

using utility::vector1;
using namespace core::sequence;
using std::map;
using std::string;
using core::pose::Pose;
using core::Size;
using core::Real;

static basic::Tracer tr2("brunette.tj_util");
///@brief inputs the sequence alignments and keeps them ordered the way they
/// were in the alingment.filt file.

vector1<SequenceAlignment> input_alignments(){
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace core::sequence;
  vector1< std::string > align_fns = option[ in::file::alignment ]();
  vector1< SequenceAlignment > alns = core::sequence::read_aln(
				    option[ cm::aln_format ](), align_fns[1]
				  );
  return(alns);
}

/// @brief inputs the sequence alignments and maps them to to either pdbid or alignment file name.
map<string,SequenceAlignment> input_alignmentsMapped(bool mapToPdbid){
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  map<string,SequenceAlignment> alns;
  map<string,SequenceAlignment>::iterator location_alns;
  vector1< std::string > align_fns = option[ in::file::alignment ]();
  for ( Size ii = 1; ii <= align_fns.size(); ++ii ) {
    vector1< SequenceAlignment > tmp_alns = core::sequence::read_aln(
			     option[ cm::aln_format ](), align_fns[ii]
								     );
    for ( Size jj = 1; jj <= tmp_alns.size(); ++jj ) {
			string mapToName;
			if(mapToPdbid == true){
				string aln_id = tmp_alns[jj].sequence(2)->id();
				string pdbid = aln_id.substr(0,aln_id.length()-2);
				mapToName = pdbid;
			}
			else{
				std::stringstream numbConvert;
				string alnBaseName = utility::file_basename(option[ in::file::alignment ]()[ii]);
				numbConvert << alnBaseName << "_" << jj;
				string aln_id = numbConvert.str();
				mapToName = aln_id;
					}
			alns.insert(std::pair<string,SequenceAlignment>(mapToName,tmp_alns[jj]));
    }
  }
  return(alns);
}

//gets the template poses from the cmd line.
std::map< std::string, core::pose::Pose > poses_from_cmd_line(
	 utility::vector1< std::string > const & fn_list) {
  using std::string;
  using core::pose::Pose;
  using utility::file::file_exists;
  using core::import_pose::pose_from_pdb;
  using namespace core::chemical;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
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

///@add side chains onto partial threads
void add_side_chains_partialthread(Pose pose){
  using namespace core::scoring;
  ScoreFunctionOP scorefxn( getScoreFunction() );
  // repack missing sidechains
  core::id::AtomID_Mask missing( true );
  core::pose::initialize_atomid_map( missing, pose );
  tr2.Debug << "repacking residues on pose with ScoreFunction: " << std::endl;
  core::pack::pack_missing_sidechains( pose, missing );
  tr2.Debug << "setting up ideal hydrogen geometry on all residues."
	   << std::endl;
  for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
    core::conformation::ResidueOP iires = pose.residue( ii ).clone();
    core::conformation::idealize_hydrogens( *iires, pose.conformation() );
    pose.replace_residue( ii, *iires, false );
  }
  tr2.Debug << "optimizing hydrogen placement with the packer."
	   << std::endl;
  core::pack::optimize_H_and_notify( pose, missing );
  scorefxn->set_weight( core::scoring::peptide_bond, 1.0 );
  (*scorefxn)(pose);
  scorefxn->show( tr2.Debug, pose );
}

/// @brief creates a partial thread for an input list of alignments
map<string,Pose> generate_partial_threads(map<string,SequenceAlignment> alnData, map< string, Pose > templateData,string query_sequence, bool add_sidechains){
  using namespace protocols::comparative_modeling;
  using namespace core::chemical;
  map<string,Pose> partialThreads;
  map<string,SequenceAlignment>::iterator location_aln;
  map<string,Pose>::iterator location_template;
  location_aln=alnData.begin();
  while(location_aln != alnData.end()){
    string pdbid = location_aln->first;
    core::pose::Pose query_pose;
    core::pose::make_pose_from_sequence(
		    query_pose, query_sequence, *(rsd_set_from_cmd_line()));
    location_template = templateData.find(pdbid);
    PartialThreadingMover mover( location_aln->second, location_template->second );
    mover.apply( query_pose );
    if(add_sidechains)
      add_side_chains_partialthread(query_pose);
    partialThreads.insert(std::pair<string,Pose>(pdbid,query_pose));
    location_aln++;
  }
  return(partialThreads);
}

/// @brief calculates burial 
vector1<bool> calculate_surface_exposure(Pose pose){
  vector1<core::Real> normalized_rsd_sasa;
  core::id::AtomID_Map< core::Real > atom_sasa;
  vector1< core::Real > residue_sasa;
  vector1< bool > surface_exposed;
  core::Real const probe_radius(1.4);
  //1.4 is the radius of water
  core::scoring::calc_per_atom_sasa( pose, atom_sasa, residue_sasa, probe_radius);
  //some constants measured by Yifan. These account for the size difference of the amino acids.  If using a partial thread make sure side chains have been added
  utility::vector1< core::Real > exposed_rsd_sasa(20);
  exposed_rsd_sasa[  1]  = 170; // 1 A
  exposed_rsd_sasa[  2]  = 170; // 2 C
  exposed_rsd_sasa[  3]  = 210; // 3 D
  exposed_rsd_sasa[  4]  = 250; // 4 E
  exposed_rsd_sasa[  5]  = 290; // 5 F
  exposed_rsd_sasa[  6]  = 170; // 6 G
  exposed_rsd_sasa[  7]  = 220; // 7 H
  exposed_rsd_sasa[  8]  = 230; // 8 I
  exposed_rsd_sasa[  9]  = 260; // 9 K
  exposed_rsd_sasa[ 10]  = 230; // 10 L
  exposed_rsd_sasa[ 11]  = 240; // 11 M
  exposed_rsd_sasa[ 12]  = 190; // 12 N
  exposed_rsd_sasa[ 13]  = 220; // 13 P
  exposed_rsd_sasa[ 14]  = 220; // 14 Q
  exposed_rsd_sasa[ 15]  = 260; // 15 R
  exposed_rsd_sasa[ 16]  = 180; // 16 S
  exposed_rsd_sasa[ 17]  = 200; // 17 T
  exposed_rsd_sasa[ 18]  = 200; // 18 V
  exposed_rsd_sasa[ 19]  = 300; // 19 W
  exposed_rsd_sasa[ 20]  = 290; // 20 Y
  for(int ii=1; ii<= pose.total_residue(); ++ii){
    normalized_rsd_sasa.push_back(residue_sasa[ii]/exposed_rsd_sasa[pose.residue(ii).type().aa()]);
    if(normalized_rsd_sasa[ii]>.1)
      surface_exposed.push_back(true);
    else
      surface_exposed.push_back(false);
  }
  //to output bfactors uncomment
  /*core::id::AtomID_Map< Real > bfactors;
  core::pose::initialize_atomid_map( bfactors, pose, 0.0 );
  for(int ii=1; ii<= pose.total_residue(); ++ii){
    for ( Size jj = 1, natoms = pose.residue_type(ii).natoms(); jj <= natoms; ++jj ){
      bfactors[core::id::AtomID(jj,ii)] = normalized_rsd_sasa[ii];
    }
    std::cout << "ii" << ii << "-" << normalized_rsd_sasa[ii] << std::endl;
  }
  string output_fn("b_factor.pdb" );
  utility::io::ozstream output_stream( output_fn );
  core::io::pdb::dump_bfactor_pdb( pose, bfactors, output_stream );
  output_stream.close();
  */
  return(surface_exposed);
}


/// @brief calculates burial 
bool gap_res_re(Size gapStartRes, Size gapEndRes, Size resToStart, Size resToEnd, Size goalRes, Pose pose){
  Size CA_CA_DIST = 3.16; //I measured several beta-sheets and found that the length between CA was exactly 3.16 in 2 cases.(Both cases 5 residues 15.8ang.  In loops the longest I got was 1.37 angstroms.
  distStartToGoal =
  distEndToGoal = 
   
  pose(position gap start) dist to end
    pose(position gap end dist to end
