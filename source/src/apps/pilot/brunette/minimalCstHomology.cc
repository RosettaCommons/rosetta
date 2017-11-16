// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file:   apps/pilot/brunette/minimalCstHomology.cc
///
/// @brief:  Generates the minimal coordinate constraints necessary for fixing homology model. Also creates fasta with virtual.
/// @usage: -in::file::alignment alignment.filt -in::file::template_pdb <pdb names from within the alignment> -in:file:coordCsts <pdb names from within the alignment> -in:file:fasta <fasta name>
/// @author TJ Brunette {tjbrunette@gmail.com}


#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <core/types.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <basic/Tracer.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>

#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>
// Unit headers
#include <core/id/SequenceMapping.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/cst_util.hh>


//utilities
#include <utility/file/FileName.hh>
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>


#include <devel/init.hh>
//#include <devel/init.hh>
// option key includes
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>

#include <fstream>
#include <map>
#include <set>
#include <sstream>

#include <devel/cstEnergyBalance/minimalCstRelaxUtil.hh>

static basic::Tracer tr( "minimalCstHomology" );

using utility::vector1;
using core::Size;
using core::Real;
using std::string;
using std::map;
using std::set;
using namespace core::sequence;
using core::pose::Pose;

namespace minimalCstHomology {
basic::options::FileVectorOptionKey coordCstFiles("minimalCstHomology:coordCstFiles");
basic::options::BooleanOptionKey only_res_out("minimalCstHomology:only_res_out");
}

void output_alignments(vector1 <SequenceAlignment> alns, std::ostream & out){
	for ( Size ii = 1; ii <= alns.size(); ++ii ) {
		alns[ii].printGrishinFormat(out);
	}
}

//stores the coordinate constraints mapped to the pdb assuming the function is name pdbid_.coordCsts
map<string,std::set<Size> > input_coordCstsMapped(vector1< std::string >const & fn_list){
	using namespace devel::cstEnergyBalance;
	using utility::file::file_exists;
	map<string,std::set<Size> > coordCsts;
	typedef vector1< string >::const_iterator iter;
	for ( iter it = fn_list.begin(), end = fn_list.end(); it != end; ++it ) {
		if ( file_exists(*it) ) {
			string name = utility::file_basename( *it );
			string pdbid = name.substr( 0, 5 );
			std::set<Size> tmp_coordCsts = input_coordCsts(*it);
			coordCsts.insert(std::pair<string , std::set<Size> >(pdbid,tmp_coordCsts));
		}
	}
	return(coordCsts);
}
//gets the poses from the cmd line.
std::map< std::string, core::pose::Pose >
poses_from_cmd_line(
	utility::vector1< std::string > const & fn_list
) {
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
//removec coordinate constraints near gap
std::set< Size > removeConstraintsNearGap(std::set< Size > caAtomsToConstrainAll,SequenceAlignment aln, Size gapSize){
	using namespace core::id;
	using namespace core::sequence;
	using namespace devel::cstEnergyBalance;
	std::set< Size >::iterator caAtomsToConstrainAll_start,caAtomsToConstrainAll_stop;
	//*identify location of gaps----
	Size EXCLUDE_NRES_TAIL = 4;
	core::id::SequenceMapping mapping_( aln.sequence_mapping( 2, 1 ) );
	Size nres = mapping_.size1();
	std::set< Size> caAtomsToConstrain;
	std::set< Size > unaligned_residues;
	//gets the first and last residues
	Size firstRes,lastRes;
	get_terminal_aln_res(aln,2,firstRes,lastRes);
	for ( Size resi = firstRes+EXCLUDE_NRES_TAIL; resi <= lastRes-EXCLUDE_NRES_TAIL; resi++ ) {
		Size t_resi = mapping_[ resi ];
		// gap checks
		//fpd  First check to see if there is a gap in the alignment
		bool gap_exists =
			t_resi == 0 || // query residue maps to a gap (aln1)
			( resi > 1    && mapping_[ resi - 1 ] != t_resi - 1 ) || // last residue was gapped
			( resi < nres && mapping_[ resi + 1 ] != t_resi + 1 ); // next residue is gapped
		if ( gap_exists ) {
			unaligned_residues.insert( resi );
		}
	}
	//**remove gaps from atom sequence---
	caAtomsToConstrainAll_start = caAtomsToConstrainAll.begin();
	caAtomsToConstrainAll_stop = caAtomsToConstrainAll.end();
	while ( caAtomsToConstrainAll_start != caAtomsToConstrainAll_stop ) {
		bool gap_exists = false;
		//check for residues near a gap or within a gap
		for ( int ii = -1*(int)gapSize; ii<=(int)gapSize; ii++ ) {
			if ( unaligned_residues.find(*caAtomsToConstrainAll_start+ii) !=unaligned_residues.end() ) {
				gap_exists = true;
			}
		}
		//check for residues not in the alignment but first remap. A position of 0 is not in the template
		Size templateResidue_id = aln.sequence_mapping(2,1)[*caAtomsToConstrainAll_start];
		if ( !gap_exists && (templateResidue_id ==0) ) {
			gap_exists = true;
		}
		if ( !gap_exists ) {
			caAtomsToConstrain.insert(*caAtomsToConstrainAll_start);
		}
		caAtomsToConstrainAll_start++;

	}
	return caAtomsToConstrain;
}


//one alignment and coord csts for each protein and r
bool calc_outputCoordCsts(
	SequenceOP fastaSequenceOP,
	map<string, SequenceAlignment> alnData,
	map< string, Pose > poseData,
	map<string, set<Size> > coordCstsData,
	Size nResFromGapExclude,
	Size minNumCoordCsts,
	bool only_res_out )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace devel::cstEnergyBalance;
	map<string,set <Size> >::iterator location_caAtomsToConstrain;
	map<string, Pose>::iterator location_pdb;
	map<string,SequenceAlignment>::iterator location_aln;
	// Check for the presence of all coordinate constraints.
	// (Some may be blank.
	// The blank coordinate constraints confer meaning that the minimal coordinate constraints have been computed.)
	location_aln=alnData.begin();
	while ( location_aln != alnData.end() ) {
		string pdbid = location_aln->second.sequence(2)->id().substr(0,5);
		location_pdb = poseData.find(pdbid);
		if ( location_pdb == poseData.end() ) {
			utility_exit_with_message( "Missing required pdb"+pdbid);
		}
		location_caAtomsToConstrain = coordCstsData.find(pdbid);
		if ( location_caAtomsToConstrain == coordCstsData.end() ) {
			utility_exit_with_message( "Missing required CA constraints"+pdbid);
		}
		string alnName = location_aln->first;
		std::ofstream alnOut( alnName.c_str() );
		vector1 <SequenceAlignment> alns;
		alns.push_back(location_aln->second);
		output_alignments(alns,alnOut);
		alnOut.close();
		//orig name : std::string outCoordCstsName = alnName+".coordCsts";
		std::string outCoordCstsName = location_aln->second.sequence(2)->id() + ".coordCsts";
		std::ofstream coordCstsOut( outCoordCstsName .c_str() );
		std::set< Size > caAtomsToConstrainAll = location_caAtomsToConstrain->second;
		Pose templatePose =location_pdb->second;

		std::set<Size> caAtomsToConstrain;
		caAtomsToConstrain = removeConstraintsNearGap(caAtomsToConstrainAll,alns[1],nResFromGapExclude);
		// If there are <= 3 coordinate constraints then coordinate constraints are not useful.
		if ( caAtomsToConstrain.size()>minNumCoordCsts ) {
			output_coordCsts(caAtomsToConstrain,
				coordCstsOut,
				templatePose,
				alns[1],
				fastaSequenceOP->sequence(),
				only_res_out);
		}
		location_aln++;
	}
	return true;  // There was no return in this non-void function. ~ Labonte
}


int main( int argc, char * argv [] ) {
	try {

		using namespace devel::cstEnergyBalance;
		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::io::silent;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::import_pose::pose_stream;
		using utility::file_basename;
		option.add (minimalCstHomology::coordCstFiles, "coordinate cst files");
		option.add (minimalCstHomology::only_res_out, "which way to make coordinate constraints").def(false);
		devel::init(argc, argv);
		Size N_RES_FROM_GAP_EXCLUDE = 4;
		Size MIN_NUM_COORDCSTS = 3;
		//Get inputs---
		map<string,std::set<Size> > coordCstData = input_coordCstsMapped(option[minimalCstHomology::coordCstFiles]());
		vector1<SequenceAlignment> alnData = input_alignments();
		map< string, SequenceAlignment> alnDataMapped = input_alignmentsMapped(false);
		map< string, Pose > poseData = poses_from_cmd_line(
			option[ in::file::template_pdb ]());
		vector1< SequenceOP > fastaSequenceOP = read_fasta_file(option[ in::file::fasta ]()[1]);
		//calc coordinate constraints
		calc_outputCoordCsts(fastaSequenceOP[1],alnDataMapped,poseData,coordCstData,N_RES_FROM_GAP_EXCLUDE,MIN_NUM_COORDCSTS,option[minimalCstHomology::only_res_out]());

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
