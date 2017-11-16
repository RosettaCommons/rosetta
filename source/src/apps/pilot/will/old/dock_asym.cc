#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <apps/pilot/will/will_util.ihh>


OPT_1GRP_KEY( Integer      , cxdock, sphere       )
OPT_1GRP_KEY( Real         , cxdock, ang_samp     )
OPT_1GRP_KEY( Real         , cxdock, clash_dis    )
OPT_1GRP_KEY( Real         , cxdock, contact_dis  )
OPT_1GRP_KEY( Real         , cxdock, hb_dis       )
OPT_1GRP_KEY( Real         , cxdock, num_contacts )
OPT_1GRP_KEY( Boolean      , cxdock, dumpfirst    )


void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::s );
	NEW_OPT( cxdock::clash_dis   ,"min acceptable contact dis", 3.5 );
	NEW_OPT( cxdock::contact_dis ,"max acceptable contact dis", 7.0 );
	NEW_OPT( cxdock::num_contacts ,"required no. contacts", 30.0 );
	NEW_OPT( cxdock::sphere       ,"sph points", 1812 );
	NEW_OPT( cxdock::ang_samp     ,"num degrees", 3.0 );
	NEW_OPT( cxdock::dumpfirst    ,"stop on first hit", false );
}


typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
using core::id::AtomID;
using basic::options::option;
using core::pose::Pose;
using core::Real;
using core::scoring::ScoreFunctionOP;
using core::Size;
using numeric::max;
using numeric::min;
using numeric::random::gaussian;
using numeric::random::uniform;
using numeric::rotation_matrix_degrees;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
using std::cerr;
using std::cout;
using std::string;
using utility::io::izstream;
using utility::io::ozstream;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_file;
using Xform;
using core::conformation::ResidueOP;


static basic::Tracer TR( "dock_asym" );
static core::io::silent::SilentFileData sfd;

#include <apps/pilot/will/sicfast.ihh>

struct Hit {
	int islide,irt,cbc,sym;
	Xform s1,s2;
	Hit(int is, int ir, int cb, int sm) : islide(is),irt(ir),cbc(cb),sym(sm) {}
};
bool cmp(Hit i,Hit j) { return i.cbc > j.cbc; }

void dump_points_pdb(vector1<Vec> & p, string fn) {
	std::ofstream o(fn.c_str());
	for ( Size i = 1; i <= p.size(); ++i ) {
		string rn = "VIZ";
		o<<"HETATM"<<I(5,i)<<' '<<" CA "<<' '<<rn<<' '<<"A"<<I(4,i)<<"    "<<F(8,3,p[i].x())<<F(8,3,p[i].y())<<F(8,3,p[i].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	}
	o.close();
}

void dock(Pose const & init1,
	Pose const & init2,
	Pose const & native,
	string fn1,
	string fn2,
	vector1<Vec> const & ssamp
) {
	using namespace basic::options;

	Real cbcth = basic::options::option[basic::options::OptionKeys::cxdock::num_contacts]();

	vector1<Real> asamp;
	for ( Real i = -5.0; i <= 180/option[OptionKeys::cxdock::ang_samp]()+5.0; ++i  ) {
		asamp.push_back( i*option[OptionKeys::cxdock::ang_samp]() ); // just 0-179
	}

	vector1<int> cstres_A,cstres_B;
	cstres_A.push_back( 4);
	cstres_A.push_back( 6);
	cstres_A.push_back( 7);
	cstres_A.push_back( 8);
	cstres_A.push_back( 9);
	cstres_A.push_back(10);
	cstres_A.push_back(11);
	cstres_A.push_back(13);
	cstres_A.push_back(14);
	cstres_A.push_back(15);
	cstres_A.push_back(16);
	cstres_A.push_back(17);
	cstres_A.push_back(18);
	cstres_A.push_back(19);
	cstres_A.push_back(28);
	cstres_A.push_back(31);
	cstres_A.push_back(33);
	cstres_A.push_back(34);
	cstres_A.push_back(35);
	cstres_A.push_back(37);

	cstres_B.push_back( 9);
	cstres_B.push_back(10);
	cstres_B.push_back(11);
	cstres_B.push_back(14);
	cstres_B.push_back(15);
	cstres_B.push_back(16);
	cstres_B.push_back(17);
	cstres_B.push_back(22);
	cstres_B.push_back(23);
	cstres_B.push_back(24);
	cstres_B.push_back(28);
	cstres_B.push_back(35);
	cstres_B.push_back(39);
	cstres_B.push_back(40);

	// set up n-ca-c-o-cb ord arrays
	vector1<Vec> bbi1tmp,cbi1tmp;
	for ( int ir = 1; ir <= init1.size(); ++ir ) {
		if ( !init1.residue(ir).is_protein() ) continue;
		for ( int ia = 1; ia <= 5; ++ia ) bbi1tmp.push_back(init1.xyz(AtomID( min(ia, (init1.residue(ir).has("CB"))?5:4 ) ,ir)));
		if ( std::find(cstres_A.begin(),cstres_A.end(),ir)!=cstres_A.end() ) if ( init1.residue(ir).has("CB") ) cbi1tmp.push_back(init1.xyz(AtomID(5,ir)));
	} vector1<Vec> const bbi1(bbi1tmp), cbi1(cbi1tmp);

	vector1<Vec> bbi2tmp,cbi2tmp;
	for ( int ir = 1; ir <= init2.size(); ++ir ) {
		if ( !init2.residue(ir).is_protein() ) continue;
		for ( int ia = 1; ia <= 5; ++ia ) bbi2tmp.push_back(init2.xyz(AtomID( min(ia, (init2.residue(ir).has("CB"))?5:4 ) ,ir)));
		if ( std::find(cstres_B.begin(),cstres_B.end(),ir)!=cstres_B.end() ) if ( init2.residue(ir).has("CB") ) cbi2tmp.push_back(init2.xyz(AtomID(5,ir)));
	} vector1<Vec> const bbi2(bbi2tmp), cbi2(cbi2tmp);

#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic,1)
#endif
	for ( int islide = 1; islide <= ssamp.size(); ++islide ) {
#ifdef USE_OPENMP
		#pragma omp critical
#endif
		if ( islide%1==0 ) TR << islide << " of " << ssamp.size() << endl;

		// set init1 with silde axis along Z
		Mat const Rslide( rotation_matrix_degrees(
			Vec(0,0,1).cross(ssamp[islide]).normalized(), // axis
			-acos(numeric::max(-1.0,numeric::min(1.0,Vec(0,0,1).dot(ssamp[islide]))))*180.0/numeric::constants::d::pi )); // angle
		vector1<Vec> bb1(bbi1); for ( vector1<Vec>::iterator i = bb1.begin(); i != bb1.end(); ++i ) *i = Rslide*(*i);
		vector1<Vec> cb1(cbi1); for ( vector1<Vec>::iterator i = cb1.begin(); i != cb1.end(); ++i ) *i = Rslide*(*i);

		// loop over all rotations of init2
		for ( int isrot = 1; isrot <= ssamp.size(); ++isrot ) {
			for ( int irot = 1; irot <= asamp.size(); ++irot ) {
				Mat const Rrot = rotation_matrix_degrees(ssamp[isrot],asamp[irot]);
				vector1<Vec> bb2(bbi2); for ( vector1<Vec>::iterator i = bb2.begin(); i != bb2.end(); ++i ) *i = Rrot*(*i);
				vector1<Vec> cb2(cbi2); for ( vector1<Vec>::iterator i = cb2.begin(); i != cb2.end(); ++i ) *i = Rrot*(*i);

				Real cbc;
				Real t = sicfast(bb1,bb2,cb1,cb2,cbc);

				if ( cbc > option[OptionKeys::cxdock::num_contacts]() ) {
					string tag = str(uniform());
#ifdef USE_OPENMP
					#pragma omp critical
#endif
					{
						Pose p(init1),q(init2);
						alignaxis(p,Vec(0,0,1),ssamp[islide],Vec(0,0,0));
						rot_pose(q,ssamp[isrot],asamp[irot]);
						trans_pose(p,Vec(0,0,t));
						p.dump_pdb("hit"+tag+"A.pdb");
						q.dump_pdb("hit"+tag+"B.pdb");

						p.append_residue_by_jump(q.residue(1),1);
						for ( Size i=2; i <= q.size(); ++i ) p.append_residue_by_bond(q.residue(i));
						// p.dump_pdb("test.pdb");
						// native.dump_pdb("nat.pdb");
						// utility_exit_with_message("oraistne");
						Real rms = core::scoring::CA_rmsd(native,p);

						TR << "HIT! " << tag << " " << rms << " " << cbc << " " << islide << " " << isrot << " " << irot << endl;
					}
				}
			}
		}
	}

}


int main(int argc, char *argv[]) {

	try {

		register_options();
		devel::init(argc,argv);
		using namespace basic::options::OptionKeys;


		Size const NSS = basic::options::option[basic::options::OptionKeys::cxdock::sphere]();
		vector1<Vec> ssamp(NSS);
		{
			izstream is;
			basic::database::open(is,"sampling/spheres/sphere_"+str(NSS)+".dat");
			for ( int i = 1; i <= NSS; ++i ) {
				Real x,y,z;
				is >> x >> y >> z;
				ssamp[i] = Vec(x,y,z);
			}
			is.close();
		}

		Pose native;
		core::import_pose::pose_from_file(native,"input/native.pdb", core::import_pose::PDB_file);


		for ( Size ifn = 1; ifn <= option[in::file::s]().size(); ifn+=2 ) {
			string fn1 = option[in::file::s]()[ifn+0];
			string fn2 = option[in::file::s]()[ifn+1];
			Pose pnat1,pnat2;
			TR << "searching " << fn1 << " vs " << fn2 << std::endl;
			{
				core::import_pose::pose_from_file(pnat1,fn1, core::import_pose::PDB_file);
				trans_pose(pnat1,-center_of_geom(pnat1,1,pnat1.size()));
				core::scoring::dssp::Dssp dssp(pnat1);
				dssp.insert_ss_into_pose(pnat1);
			}{
				core::import_pose::pose_from_file(pnat2,fn2, core::import_pose::PDB_file);
				trans_pose(pnat2,-center_of_geom(pnat2,1,pnat2.size()));
				core::scoring::dssp::Dssp dssp(pnat2);
				dssp.insert_ss_into_pose(pnat2);
			}
			dock(pnat1,pnat2,native,fn1,fn2,ssamp);
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

