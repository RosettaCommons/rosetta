// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Sean Joseph Wu
/// @author James Thompson

// libRosetta headers

#include <core/types.hh>

#include <devel/init.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/rms_util.hh>
#include <core/pose/PDB_Info.hh>
//#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>

#include <protocols/enzdes/AddorRemoveCsts.hh>

#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <numeric/xyzVector.hh>

#include <core/sequence/util.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <protocols/comparative_modeling/ExtraThreadingMover.hh>
#include <protocols/enzdes/enzdes_util.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>

using core::Size;
using core::Real;
using utility::vector1;

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers
#include <core/pose/util.hh>

int
main( int argc, char * argv [] ) {
	try {
	using namespace core::chemical;
	using namespace basic::options::OptionKeys;
	using namespace basic::options;
	using namespace protocols::enzdes::enzutil;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose;
	using core::id::AtomID;
	using core::id::AtomID_Map;
	using utility::vector1;
	using core::id::SequenceMapping;
	using core::sequence::SequenceAlignment;
	using protocols::comparative_modeling::ExtraThreadingMover;

	devel::init( argc, argv );

	// setup residue types
	ResidueTypeSetCAP rsd_set =
		ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	PoseInputStreamOP pdb_input(
		new PDBPoseInputStream( option[ in::file::s ]() )
	);

	core::pose::Pose comparative_modeling_pose, template_pose;
	pdb_input->fill_pose( comparative_modeling_pose, *rsd_set );
	pdb_input->fill_pose( template_pose, *rsd_set );

	vector1< std::string > align_fns = option[ in::file::alignment ]();
	vector1< SequenceAlignment > alns = core::sequence::read_aln(
		option[ cm::aln_format ](), align_fns.front()
	);
	SequenceAlignment aln = alns.front();

	// read enzyme constraints into template_pose
	protocols::enzdes::AddOrRemoveMatchCsts cstmover;
	cstmover.set_cst_action( protocols::enzdes::ADD_NEW );
	cstmover.set_accept_blocks_missing_header( false );
	cstmover.apply(template_pose);

	//Saves old comparative model length
  int old_comparative_model_length = comparative_modeling_pose.total_residue();

	//Steals Ligand: User designated or residues 357 through 360
	vector1 < Size > residues_to_steal;
	if ( option[ in::target_residues ].user() ) {
		residues_to_steal = option[ in::target_residues ]();
	} else {
		for ( Size ii = 357; ii <= 360; ++ii ) {
			//std::cout << "stealing residue " << ii << std::endl;
  		residues_to_steal.push_back(ii); // peptide from KUMAMO_PQLP
  	}
	}

	ExtraThreadingMover mover( aln, template_pose, residues_to_steal );
	mover.apply(comparative_modeling_pose); //Does a superimpose_pose

	char ligand_chain = core::pose::chr_chains[comparative_modeling_pose.residue(old_comparative_model_length).chain() + 1];

	//Increases pose length in the pdbinfo to accomodate ligand
	comparative_modeling_pose.pdb_info() -> resize_residue_records( comparative_modeling_pose.total_residue() );
  //comparative_modeling_pose.pdb_info() -> resize_atom_records( comparative_modeling_pose );

	//renumbers pdbinfo for the ligand. No longer changes numbers in PDB file...Don't know why
	//Current Fix: Adds on ligand as next number in the chain. If protein ends at 360, then ligand start will be 361
	for(Size i = old_comparative_model_length + 1; i <= comparative_modeling_pose.total_residue(); ++i) {
		//comparative_modeling_pose.pdb_info() -> set_resinfo(i, ligand_chain, i - old_comparative_model_length, ' ');
		comparative_modeling_pose.pdb_info() -> set_resinfo(i, ligand_chain, i, ' ');
		//std::cout << comparative_modeling_pose.pdb_info() -> pose2pdb(i) << std::endl;
		//std::cout << "Renumbering to " << (i - old_comparative_model_length) << " at " << i << std::endl;
	}

	//Creates Remarks object from template pose remarks
	core::pose::Remarks remarks( template_pose.pdb_info()->remarks() );
	std::cout << "currently have " << remarks.size() << " remarks." << std::endl;

	//Create sequence alignment from template(2) to query(1)
	SequenceMapping template_to_query( aln.sequence_mapping(2,1) );

	//Creates comparative remarks in a std::vector by converting from template remarks
	std::vector <core::pose::RemarkInfo> remark_vector;
  for( core::Size i = 0; i < remarks.size(); ++i ){
 		int pdbposA(0), pdbposB(0);
	  std::string chainA(""), chainB(""), resA(""), resB("");
		core::Size cst_block(0), exgeom_id( 0 );
		//protocols::enzdes::enzutil::split_up_remark_line( remarks[i].value, chainA, resA, pdbposA, chainB, resB, pdbposB, cst_block, exgeom_id );
		protocols::toolbox::match_enzdes_util::split_up_remark_line( remarks[i].value, chainA, resA, pdbposA, chainB, resB, pdbposB, cst_block, exgeom_id );

		//std::cout<<"chainA:"<<chainA<<"	resA:"<<resA<<" pdbposA:"<<pdbposA<<"	chainB:" << chainB<< "	resB:"<<resB<<" pdbposB:"<<pdbposB<<std::endl;

    std::string query_chainA, query_chainB, query_resA, query_resB;
		int query_posA, query_posB;
		char model_last_chain = core::pose::chr_chains[comparative_modeling_pose.residue(old_comparative_model_length).chain()];
		//Maps new residue positions for comparative remarks
		//Changes remarks to be consistent with ligand numbering
		if( chainA[0] <= model_last_chain) {
			query_posA = template_to_query[pdbposA];
		} else {
			query_posA = pdbposA + old_comparative_model_length;
			//query_posA = pdbposA;
		}
    if( chainB[0] <= model_last_chain) {
      query_posB = template_to_query[pdbposB];
    } else {
			query_posB = pdbposB + old_comparative_model_length;
    	//query_posB = pdbposB;
    }
    query_resA = resA;
    query_chainA = chainA;
		query_resB = resB;
    query_chainB = chainB;
		core::pose::RemarkInfo cm_remark;
		cm_remark.num = 0;
		//cm_remark.value = protocols::enzdes::enzutil::assemble_remark_line(query_chainA, query_resA,
		cm_remark.value = protocols::toolbox::match_enzdes_util::assemble_remark_line(query_chainA, query_resA,
		query_posA, query_chainB, query_resB, query_posB, cst_block, exgeom_id);
		remark_vector.push_back(cm_remark);
		std::cout << "cm_remark: " << cm_remark << std::endl;
		//Mutates the residues to fit catalytic constraints
		if(query_chainA[0] <= model_last_chain && comparative_modeling_pose.residue(query_posA).name3() != template_pose.residue(pdbposA).name3()) {
			comparative_modeling_pose.replace_residue(query_posA, template_pose.residue(pdbposA), true);
		}
	  if(query_chainB[0] <= model_last_chain && comparative_modeling_pose.residue(query_posB).name3() != template_pose.residue(pdbposB).name3()) {
			comparative_modeling_pose.replace_residue(query_posB, template_pose.residue(pdbposB), true);
		}
  }


//
//Attempt to renumber ligand to 1 B, 2 B, etc...
/*
	vector1< Size > renumber_vector;
	for(Size i = 1; i <= comparative_modeling_pose.total_residue(); ++i) {
		if( i <= old_comparative_model_length ) {
			renumber_vector.push_back(i);
		} else {
			renumber_vector.push_back(i-old_comparative_model_length);
		}
	}

	for(Size i=1; i<= comparative_modeling_pose.total_residue(); ++i ) {
		std::cout<<"testing:"<<renumber_vector[i]<<std::endl;
		((core::conformation::Residue)(comparative_modeling_pose.residue(i))).update_sequence_numbering(renumber_vector);
		std::cout<<comparative_modeling_pose.residue(i)<<std::endl;
  }
*/
//
	//Superimposes again to make the mutated residues in the correct spot.
  vector1< Size > residues;

  if ( option[ in::target_residues ].user() ) {
    // user-specified residues
    residues = option[ in::target_residues ]();
  } else {
    // use all aligned residues
    for ( Size ii = 1; ii <= template_pose.total_residue(); ++ii ) {
      residues.push_back(ii);
    }
  }

	//Place remarks into comparative modeling pose
  core::pose::Remarks query_remarks;
	for (core::Size ii = 0; ii < remark_vector.size(); ++ii) {
		std::cout << remark_vector[ii] << std::endl;
		query_remarks.push_back(remark_vector[ii]);
	}
	comparative_modeling_pose.pdb_info()->remarks(query_remarks);

	//Apply constraints to comparative modeling pose
	cstmover.apply(comparative_modeling_pose);

	vector1< std::string > fns = option[ in::file::s ]();
	std::string output_name( fns[1] + ".super.pdb" );
	if ( option[ out::file::o ].user() ) {
		output_name = option[ out::file::o ]();
	}

	comparative_modeling_pose.dump_pdb( output_name );
	std::cout << "wrote pdb with name " << output_name << std::endl;

	//Testing
	//comparative_modeling_pose.pdb_info() -> show(std::cout);

	 } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
} // int main
