// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

// headers

	#include <protocols/motif_hash/motif_hash_stuff.hh>

	#include <ObjexxFCL/format.hh>
	#include <ObjexxFCL/string.functions.hh>
	#include <basic/Tracer.hh>
	#include <basic/database/open.hh>
	#include <basic/options/keys/in.OptionKeys.gen.hh>
	#include <basic/options/keys/out.OptionKeys.gen.hh>
	#include <basic/options/keys/mh.OptionKeys.gen.hh>
	#include <basic/options/option_macros.hh>
	#include <core/chemical/AtomType.hh>
	#include <core/chemical/ChemicalManager.hh>
	#include <core/conformation/symmetry/util.hh>
	#include <core/conformation/ResidueFactory.hh>
	#include <core/import_pose/import_pose.hh>
	#include <core/io/silent/SilentFileData.hh>
	#include <core/pose/PDBInfo.hh>
	#include <core/pose/Pose.hh>
	#include <core/pose/annotated_sequence.hh>
	#include <core/pose/util.hh>
	#include <core/io/pdb/pose_io.hh>
	#include <core/scoring/Energies.hh>
	#include <core/scoring/EnergyGraph.hh>
	#include <core/scoring/ScoreFunctionFactory.hh>
	#include <core/scoring/ScoreTypeManager.hh>
	#include <core/scoring/dssp/Dssp.hh>
	#include <core/scoring/hbonds/HBondOptions.hh>
	#include <core/scoring/methods/EnergyMethodOptions.hh>
	#include <core/scoring/packing/compute_holes_score.hh>
	#include <core/scoring/rms_util.hh>
	#include <core/scoring/sasa.hh>
	#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
	#include <devel/init.hh>
	#include <numeric/conversions.hh>
	#include <numeric/model_quality/rms.hh>
	#include <numeric/random/random.hh>
	#include <numeric/xyz.functions.hh>
	#include <numeric/xyz.io.hh>
	#include <numeric/xyzVector.hh>
	#include <protocols/sic_dock/RigidScore.hh>
	#include <protocols/sic_dock/SICFast.hh>
	#include <protocols/sic_dock/util.hh>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>
	#include <utility/fixedsizearray1.hh>
	// #include <protocols/sic_dock/designability_score.hh>

	#include <numeric/geometry/hashing/SixDHasher.hh>
	#include <numeric/HomogeneousTransform.hh>

	#include <apps/pilot/will/will_util.ihh>

	#include <protocols/sic_dock/read_biounit.hh>

	#include <boost/foreach.hpp>

	// #include <boost/archive/text_oarchive.hpp>
	// #include <boost/archive/text_iarchive.hpp>
	// #include <boost/archive/binary_oarchive.hpp>
	// #include <boost/archive/binary_iarchive.hpp>

	// // #include <boost/serialization/detail/stack_constructor.hpp>
	// #include <boost/serialization/hash_map.hpp>
	// #include <boost/serialization/hash_collections_save_imp.hpp>
	// #include <boost/serialization/hash_collections_load_imp.hpp>

static thread_local basic::Tracer TR( "motif_hash_util" );



int main(int argc, char *argv[]) {
	try{

	devel::init(argc,argv);

	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace protocols::motif_hash;

	// tmp_test_euler();
	// utility_exit_with_message("EULER");

	if(option[mh::dump_input_pdb].user()) {
		TR << "dumping input_pdbs from  -mh:dump_input_pdb" << endl;
		vector1<string> const & fnames( option[mh::dump_input_pdb]() );
		for(int ifile = 1; ifile <= (int)fnames.size(); ++ifile ){
			Pose pose; vector1<Real> bfactors,occupancy;
			utility::vector1<int> pdbres;
			std::map<int,char> pdbchain;
			int nresmodel1;
			if( !protocols::sic_dock::read_biounit(fnames[ifile],pose,bfactors,occupancy,pdbres,pdbchain,nresmodel1,999999,true) ){
				TR.Error << "FAIL TO READ " << fnames[ifile] << endl;
				continue;
			}
		}
	}
	// if(option[mh::print_motifs_boost].user()){
	// 	TR << "reading motifs in files following -mh:print_motifs_boost" << endl;
	// 	print_motifs_boost(TR);
	// }
	if(option[mh::harvest_motifs].user()){
		TR << "harvesting motifs in files following -mh:harvest_motifs" << endl;
		harvest_motifs();
	}
	if(option[mh::print_motifs].user()){
		TR << "reading motifs in files following -mh:print_motifs" << endl;
		ResPairMotif::print_header(cout);
		print_motifs(cout);
	}
	if(option[mh::merge_motifs].user()){
		TR << "reading motifs in files following -mh:merge_motifs" << endl;
		merge_motifs();
	}
	if(option[mh::dump_motif_pdbs].user()){
		TR << "reading motifs in files following -mh:dump_motif_pdbs" << endl;
		dump_motif_pdbs();
	}
	if(option[mh::harvest_scores].user()){
		TR << "reading motif counts in files following -mh:harvest_scores" << endl;
		harvest_scores();
	}
	if(option[mh::print_scores].user()){
		TR << "reading motifs in files following -mh:print_scores" << endl;
		print_scores();
	}
	// if(option[mh::score_pdbs].user()){
	// 	if(!option[mh::xform_score_data].user()) utility_exit_with_message("specify -mh:score_data <datafile>.xh.bin.gz");
	// 	TR << "reading pdbs to score following -mh:score_data " << option[mh::xform_score_data]() << endl;
	// 	TR << "reading scoring data from file " << option[mh::xform_score_data]() << " following -mh:score_data" << endl;
	// 	score_with_counts();
	// }
	// if(option[mh::sequence_recovery].user()){
	// 	if(!option[mh::input_motifs].user()) utility_exit_with_message("specify -mh:input_motifs <datafile>.xh.bin.gz");
	// 	TR << "reading pdbs to score following -mh:input_motifs " << option[mh::input_motifs]() << endl;
	// 	TR << "reading scoring data from file " << option[mh::input_motifs]() << " following -mh:input_motifs" << endl;
	// 	sequence_recovery();
	// }
	if(option[mh::dump_matching_motifs].user()){
		if(!option[mh::input_motifs].user()) utility_exit_with_message("specify -mh:input_motifs <datafile>.xh.bin.gz");
		TR << "reading pdbs to score following -mh:input_motifs " << option[mh::input_motifs]() << endl;
		TR << "reading scoring data from file " << option[mh::input_motifs]() << " following -mh:input_motifs" << endl;
		dump_matching_motifs();
	}

	// {
	// 	utility::io::ozstream ozs("test.mfh.gz");
	// 	boost::archive::binary_oarchive out_archive(ozs);
	// 	out_archive << mh;
	// }

	// utility::io::ozstream ofs1("test1.txt"); ofs1 << mh << endl; ofs1.close();
	// MotifHash mh2;
	// utility::io::izstream ifs1("test1.txt"); ifs1 >> mh2        ; ifs1.close();
	// utility::io::ozstream ofs2("test2.txt"); ofs2 << mh2 << endl; ofs2.close();

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
    }
    return 0;
}


