// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/scoring/motif/util.hh>
#include <core/pose/motif/reference_frames.hh>
#include <numeric/xyzTransform.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <utility/io/ozstream.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <devel/init.hh>

//typedef numeric::xyzTransform<double> Xform;

int main(int argc, char *argv[]) {
	try{

	devel::init(argc,argv);

	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::scoring::motif;

	core::pose::Pose pose;
	core::import_pose::pose_from_file(pose,option[in::file::s]()[1], core::import_pose::PDB_file);
	core::scoring::dssp::Dssp(pose).insert_ss_into_pose(pose);

	core::scoring::motif::MotifHashManager & mman(*core::scoring::motif::MotifHashManager::get_instance());


	double score = 0.0;
	for(size_t ir = 1; ir <= pose.size(); ++ir){
		Xform const ibb_stub = core::pose::motif::get_backbone_reference_frame(pose,ir);
		char ss1 = pose.secstruct(ir);
		char aa1 = pose.residue(ir).name1();
		for(size_t jr = ir+1; jr <= pose.size(); ++jr){
			char ss2 = pose.secstruct(jr);
			char aa2 = pose.residue(jr).name1();
			// std::cout << ss1 << ss2 << " " << aa1 << aa2 << std::endl;
			Xform const jbb_stub = core::pose::motif::get_backbone_reference_frame(pose,jr);
			Xform const Xbb = ibb_stub.inverse() * jbb_stub;
			core::scoring::motif::XformScoreCAP xs_bb_fxn1 = mman.get_xform_score_BB_BB(ss1,ss2,aa1,aa2);
			core::scoring::motif::XformScoreCAP xs_bb_fxn2 = mman.get_xform_score_BB_BB(ss2,ss1,aa2,aa1);
			score += xs_bb_fxn1->score_of_bin(Xbb);
			score += xs_bb_fxn2->score_of_bin(Xbb.inverse());
		}
	}
	std::cout << "total motif score (higher is better): " << score << std::endl;

	///////////////// dump matching motifs example //////////////////
	ResPairMotifQuery opt(pose);
	opt.interface_only() = false;
	MotifHits hits;
	MotifHashManager::get_instance()->get_matching_motifs(opt,hits);
	if( hits.size() ){
		std::string outfile = "motif_test.pdb";
		std::cout << "dump " << hits.size() << " (before pruning) hits to " << outfile << std::endl;
		utility::io::ozstream pdbout( outfile );
		pdbout << "MODEL MAIN" << std::endl;
		pose.dump_pdb(pdbout);
		pdbout << "ENDMDL" << std::endl;
		hits.dump_motifs_pdb(pdbout);
		pdbout.close();
	} else {
		std::cout << "no mitif hits? maybe you need to pass -mh:path:motifs?" << std::endl;
	}


    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;

    }
    return 0;
}


