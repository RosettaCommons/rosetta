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
#include <core/id/AtomID.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
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
OPT_KEY( File, cdsf_match_pdb )
OPT_KEY( Real, cdsf_rms_cut )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( cdsf_max_res ,"" , 12 );
	NEW_OPT( cdsf_match_pdb ,"" , "" );
	NEW_OPT( cdsf_rms_cut ,"" , 0.2 );
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
 	for(int i = 1; i <= hbset.nhbonds(); ++i) {
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


	core::pose::Pose svd;
	core::import_pose::pose_from_file(svd,option[cdsf_match_pdb](), core::import_pose::PDB_file);
	core::pose::remove_lower_terminus_type_from_pose_residue(svd,         1     );
	core::pose::remove_upper_terminus_type_from_pose_residue(svd,svd.size());

	if(svd.size()!=3) utility_exit_with_message("oiarestn");

	utility::vector1<string> fnames = option[in::file::s]();
	for(int ifile=1; ifile <= fnames.size(); ++ifile) {
		string fn = utility::file_basename(fnames[ifile]);
		core::pose::Pose pose;
		core::import_pose::pose_from_file(pose,fnames[ifile], core::import_pose::PDB_file);
		for(int ir=2; ir <= pose.size()-4; ++ir) {
			std::map< core::id::AtomID, core::id::AtomID > m;
			for(int ii=1;ii<=5;++ii) m[core::id::AtomID(ii,ir+0)] = core::id::AtomID(ii,1);
			for(int ii=1;ii<=5;++ii) m[core::id::AtomID(ii,ir+1)] = core::id::AtomID(ii,2);
			for(int ii=1;ii<=5;++ii) m[core::id::AtomID(ii,ir+2)] = core::id::AtomID(ii,3);
			core::Real rms = rms_at_corresponding_atoms(pose,svd,m);
			if(rms < option[cdsf_rms_cut]()) {
				std::cout << rms << " " << fn << " " << ir << std::endl;
				core::pose::Pose tmp(pose);
				tmp.replace_residue(ir+0,svd.residue(1),true);
				tmp.replace_residue(ir+1,svd.residue(2),true);
				tmp.replace_residue(ir+2,svd.residue(3),true);
				tmp.dump_pdb(option[out::file::o]()+"/"+fn+"_"+ObjexxFCL::string_of(ir)+"_"+"svd.pdb");
			}

		}
	}

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

