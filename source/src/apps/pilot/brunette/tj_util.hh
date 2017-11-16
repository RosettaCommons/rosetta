#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/conformation/Residue.functions.hh>
#include <core/conformation/Residue.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>

#include <core/pose/annotated_sequence.hh>

#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

#include <core/pack/optimizeH.hh>
#include <core/pack/pack_missing_sidechains.hh>

#include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/coord_util.hh>
#include <protocols/comparative_modeling/PartialThreadingMover.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>

#include <protocols/star/Extender.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

using utility::vector1;
using namespace core::sequence;
using std::map;
using std::string;
using core::pose::Pose;
using core::Size;
using core::Real;

static basic::Tracer tr2( "brunette.tj_util" );
/// @brief inputs the sequence alignments and keeps them ordered the way they
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
			string aln_id;
			if ( mapToPdbid == true ) {
				aln_id = tmp_alns[jj].sequence(2)->id();
				std::transform(aln_id.begin(), aln_id.end(), aln_id.begin(), toupper);
			} else {
				std::stringstream numbConvert;
				string alnBaseName = utility::file_basename(option[ in::file::alignment ]()[ii]);
				numbConvert << alnBaseName << "_" << jj;
				aln_id = numbConvert.str();
			}
			alns.insert(std::pair<string,SequenceAlignment>(aln_id,tmp_alns[jj]));
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
	using core::import_pose::pose_from_file;
	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

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

//gets the template poses from the cmd line.
std::map< std::string, core::pose::Pose > poses_from_cmd_line_noPDBtag(
	utility::vector1< std::string > const & fn_list) {
	using std::string;
	using core::pose::Pose;
	using utility::file::file_exists;
	using core::import_pose::pose_from_file;
	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
	map< string, Pose > poses;

	typedef vector1< string >::const_iterator iter;
	for ( iter it = fn_list.begin(), end = fn_list.end(); it != end; ++it ) {
		if ( file_exists(*it) ) {
			Pose pose;
			core::import_pose::pose_from_file( pose, *rsd_set, *it , core::import_pose::PDB_file);
			string name = utility::file_basename( *it );
			name = name.substr( 0, 9 );
			std::transform(name.begin(), name.end(), name.begin(), toupper);
			poses[name] = pose;
		}
	}
	return poses;
}

/// @add side chains onto partial threads
void add_side_chains_partialthread(Pose pose){
	using namespace core::scoring;
	ScoreFunctionOP scorefxn( get_score_function() );
	// repack missing sidechains
	core::id::AtomID_Mask missing( true );
	core::pose::initialize_atomid_map( missing, pose );
	tr2.Debug << "repacking residues on pose with ScoreFunction: " << std::endl;
	core::pack::pack_missing_sidechains( pose, missing );
	tr2.Debug << "setting up ideal hydrogen geometry on all residues."
		<< std::endl;
	for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
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
	while ( location_aln != alnData.end() ) {
		string pdbid = location_aln->first;
		core::pose::Pose query_pose;
		core::pose::make_pose_from_sequence(
			query_pose, query_sequence, *(rsd_set_from_cmd_line()));
		location_template = templateData.find(pdbid);
		PartialThreadingMover mover( location_aln->second, location_template->second );
		mover.apply( query_pose );
		if ( add_sidechains ) {
			add_side_chains_partialthread(query_pose);
		}
		partialThreads.insert(std::pair<string,Pose>(pdbid,query_pose));
		location_aln++;
	}
	return(partialThreads);
}

/// @brief get all unaligned loops from set of alignments
map<string, const protocols::loops::Loops> get_unalignedLoopsMapped(map<string,SequenceAlignment> alnDataMapped, Size numResidues){
	using namespace protocols::star;
	using protocols::loops::Loops;
	map<string,const Loops> unalignedLoops;
	map<string, SequenceAlignment>::iterator alnDataMapped_itr;
	for ( alnDataMapped_itr=alnDataMapped.begin(); alnDataMapped_itr != alnDataMapped.end(); alnDataMapped_itr++ ) {
		core::sequence::SequenceAlignmentCOP tmpAlignment = alnDataMapped_itr->second.clone();
		Extender extender(tmpAlignment, numResidues);
		const Loops& unaligned = *(extender.unaligned());
		unalignedLoops.insert(std::pair<string,const Loops>(alnDataMapped_itr->first,unaligned));
	}
	return(unalignedLoops);
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
	for ( int ii=1; ii<= pose.size(); ++ii ) {
		normalized_rsd_sasa.push_back(residue_sasa[ii]/exposed_rsd_sasa[pose.residue(ii).type().aa()]);
		if ( normalized_rsd_sasa[ii]>.1 ) {
			surface_exposed.push_back(true);
		} else {
			surface_exposed.push_back(false);
		}
	}
	//to output bfactors uncomment
	/*core::id::AtomID_Map< Real > bfactors;
	core::pose::initialize_atomid_map( bfactors, pose, 0.0 );
	for(int ii=1; ii<= pose.size(); ++ii){
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

/// @brief superimposes aligned residues
void superimpose_pose_using_aln(Pose & mod_pose, Pose const & ref_pose, SequenceAlignment aln, Size firstAln){
	using namespace core::id;
	using namespace core::sequence;
	using namespace core::scoring;
	Size secondAln;
	if ( firstAln == 1 ) {
		secondAln = 2;
	} else {
		secondAln = 1;
	}
	//was 2,1
	SequenceMapping map = aln.sequence_mapping(firstAln,secondAln);
	core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, mod_pose, core::id::AtomID::BOGUS_ATOM_ID() ); // maps every atomid to bogus atom
	for ( Size ii=1; ii<=mod_pose.size(); ++ii ) {
		if ( map[ii] != 0 ) {
			core::id::AtomID const id1( mod_pose.residue(ii).atom_index("CA"), ii );
			core::id::AtomID const id2( ref_pose.residue(map[ii]).atom_index("CA"), map[ii] );
			atom_map[ id1 ] = id2;
		}
	}
	superimpose_pose( mod_pose, ref_pose, atom_map);
}

/// @brief inputs the first N poses

utility::vector1< core::pose::PoseOP > input_subset_poses(Size n){
	using namespace core::chemical;
	using namespace core::import_pose::pose_stream;
	utility::vector1< core::pose::PoseOP > pose_vector;
	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
	Size count = 0;
	MetaPoseInputStream input = streams_from_cmd_line();
	while ( input.has_another_pose() && count < n ) {
		core::pose::PoseOP input_poseOP;
		input_poseOP = new core::pose::Pose();
		input.fill_pose(*input_poseOP,*rsd_set);
		pose_vector.push_back(input_poseOP);
		count++;
	}
	return(pose_vector);
}

//@brief gets central pose...assumes all poses are the same length
core::pose::PoseOP get_central_pose(utility::vector1< core::pose::PoseOP> pose_vector){
	using namespace core::scoring;
	vector1< vector1 <Real> > gdt_matrix(pose_vector.size(),vector1< Real>(pose_vector.size(),0));
	for ( Size ii=1; ii<pose_vector.size(); ++ii ) {
		for ( Size jj=ii+1; jj<=pose_vector.size(); ++jj ) {
			Real tmp_gdt = CA_gdtmm(*pose_vector[ii],*pose_vector[jj]);
			gdt_matrix[ii][jj] = tmp_gdt;
			gdt_matrix[jj][ii] = tmp_gdt;
		}
	}
	map<Real,Size>  avg_gdt_pose_map;
	for ( Size ii=1; ii<=pose_vector.size(); ++ii ) {
		Real total_gdt = 0;
		for ( Size jj=1; jj<=pose_vector.size(); ++jj ) {
			total_gdt += gdt_matrix[ii][jj];
		}
		avg_gdt_pose_map[total_gdt/(Real)(pose_vector.size()-1)] = ii;
	}
	return(pose_vector[avg_gdt_pose_map.rbegin()->second]);
}

void get_disordered_regions(core::pose::PoseOP centralPose, utility::vector1< core::pose::PoseOP> pose_vector,Real threshold){
	using namespace core::scoring;
	vector1 <Real> per_position_totalDist(centralPose->total_residue(),0);
	vector1 <Real> per_position_avgDist(centralPose->total_residue(),0);
	for ( Size ii=1; ii<=pose_vector.size(); ++ii ) {
		calpha_superimpose_pose(*pose_vector[ii],*centralPose);
		for ( Size jj=1; jj<centralPose->total_residue(); ++jj ) {
			per_position_totalDist[jj] += pose_vector[ii]->residue(jj).xyz("CA").distance(centralPose->residue(jj).xyz("CA"));
		}
	}
	for ( Size ii=1; ii<=centralPose->total_residue(); ++ii ) {
		per_position_avgDist[ii] = per_position_totalDist[ii]/(Real)pose_vector.size();
		std::cout <<ii <<"-" << per_position_avgDist[ii] <<std::endl;
	}
}

