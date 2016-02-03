#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/conformation/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
// #include <devel/init.hh>


OPT_KEY( Integer, cdsf_max_res )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( cdsf_max_res ,"" , 12 );
}

void get_hb_info(core::pose::Pose & pose, int start, int stop, int & exhb, int & inhb, int & inhbsc) {
	core::scoring::ScoreFunctionOP shb = new core::scoring::ScoreFunction;
	shb->set_weight(core::scoring::hbond_lr_bb,1.0);
	shb->set_weight(core::scoring::hbond_sr_bb,1.0);
	shb->set_weight(core::scoring::hbond_bb_sc,1.0);
	shb->set_weight(core::scoring::hbond_sc   ,1.0);
	shb->score(pose);
	core::scoring::hbonds::HBondSet hbset;
 	core::scoring::hbonds::fill_hbond_set( pose, false, hbset, false, true, true, true );
 	for(int i = 1; i <= (int)hbset.nhbonds(); ++i) {
 		int dr = hbset.hbond(i).don_res();
		int ar = hbset.hbond(i).acc_res();
		bool din = start <= dr && dr <= stop;
		bool ain = start <= ar && ar <= stop;
		if( !din && !ain ) continue;
		if( din && ain ) {
			if( hbset.hbond(i).don_hatm_is_protein_backbone() && hbset.hbond(i).acc_atm_is_protein_backbone())
				inhb++;
			else
				inhbsc++;
		} else {
				exhb++;
		}
	}
}


int main(int argc, char *argv[]) {

	try {

	register_options();
	devel::init(argc,argv);
	using namespace std;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::io::izstream idone("done.list");
	string s;
	set<string> done;
	while(idone >> s) done.insert(s);

	utility::vector1<string> fnames = option[in::file::s]();
	for(int ifile=1; ifile <= (int)fnames.size(); ++ifile) {
		string fn = utility::file_basename(fnames[ifile]);
		if(done.count(fn)!=0) continue;
		cout << "BEGIN " << fn << endl;
		core::pose::Pose pose;
		core::import_pose::pose_from_file(pose,fnames[ifile], core::import_pose::PDB_file);
		for(int ir=1; ir <= (int)pose.n_residue(); ++ir) {
			if(pose.residue(ir).aa() != core::chemical::aa_cys) continue;
			for(int jr=ir+3; jr < ir+option[cdsf_max_res](); ++jr) {
				if(jr > (int)pose.n_residue()) break;
				if(pose.residue(jr).aa() != core::chemical::aa_cys) continue;
				if(!core::conformation::is_disulfide_bond( pose.conformation(),ir,jr)) continue;
				bool terminus = false;
				for(int i=ir+1; i < jr; ++i) terminus |= pose.residue(i).is_terminus();
				if(terminus) continue;
				core::pose::Pose tmp;
				tmp.append_residue_by_jump(pose.residue(ir),1);
				for(int i=ir+1; i <= jr; ++i) tmp.append_residue_by_bond(pose.residue(i));
				string outfile = option[out::file::o]+"/"+fn+"_"+ObjexxFCL::string_of(ir)+"_"+ObjexxFCL::string_of(jr)+".pdb";


				int exhb=0,inhb=0,inhbsc=0;
				get_hb_info(pose,ir,jr,exhb,inhb,inhbsc);


				using namespace ObjexxFCL::format;
				cout << "HIT " << " " << I(4,ir) << " " << I(4,jr) << " " << I(4,exhb) << " " << I(4,inhb) << " " << I(4,inhbsc) << " " << outfile << endl;
				tmp.dump_pdb(outfile);
			}
		}
		cout << "DONE " << fn << endl;
	}

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

