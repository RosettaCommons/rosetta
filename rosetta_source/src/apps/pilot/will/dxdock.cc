#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/init.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
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

#include <set>

OPT_1GRP_KEY( Integer      , cxdock, sphere       )
//OPT_1GRP_KEY( IntegerVector, cxdock, syms         )
OPT_1GRP_KEY( Real         , cxdock, clash_dis    )
OPT_1GRP_KEY( Real         , cxdock, contact_dis  )
OPT_1GRP_KEY( Real         , cxdock, num_contacts )
OPT_1GRP_KEY( Boolean      , cxdock, dumpfirst )
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::s );
//	NEW_OPT( cxdock::syms	 ,"CX symmitries", utility::vector1< Size >() );
	NEW_OPT( cxdock::clash_dis   ,"min acceptable contact dis", 3.5 );
	NEW_OPT( cxdock::contact_dis ,"max acceptable contact dis", 7.0 );	
	NEW_OPT( cxdock::num_contacts ,"required no. contacts", 20.0 );
	NEW_OPT( cxdock::sphere       ,"sph points", 8192 );
	NEW_OPT( cxdock::dumpfirst    ,"stop on first hit", false );
}


typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Matf;
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
using ObjexxFCL::fmt::F;
using ObjexxFCL::fmt::I;
using ObjexxFCL::string_of;
using std::cerr;
using std::cout;
using std::string;
using utility::io::izstream;
using utility::io::ozstream;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_pdb;
using core::kinematics::Stub;
using core::conformation::ResidueOP;


static basic::Tracer TR("dxdock");
static core::io::silent::SilentFileData sfd;

#include <apps/pilot/will/sicfast.ihh>


int num_hbonds(core::pose::Pose const & pose) {
	core::scoring::hbonds::HBondSet hbset;
	core::scoring::hbonds::fill_hbond_set( pose, false, hbset, false );
	int nhb = 0;
	for(Size i = 1; i <= hbset.nhbonds(); ++i) {
		if( !hbset.hbond(i).don_hatm_is_protein_backbone() ) continue;
		if( !hbset.hbond(i).acc_atm_is_protein_backbone() ) continue;
		nhb++;
	}
	return nhb;
}



struct Hit {
	int iss,irt,cbc,sym;
	core::kinematics::Stub s1,s2;
	Hit(int is, int ir, int cb, int sm) : iss(is),irt(ir),cbc(cb),sym(sm) {}
};
bool cmp(Hit i,Hit j) { return i.cbc > j.cbc; }

void dock(Pose const init, std::string const & fn, vector1<Vec> const & ssamp) {
	using namespace basic::options;
	
	ScoreFunctionOP shb = new core::scoring::ScoreFunction;
	shb->set_weight(core::scoring::hbond_lr_bb,1.0);
	
	vector1<Size> syms;// = basic::options::option[basic::options::OptionKeys::cxdock::syms]();
	syms.push_back(2);
	if(syms.size()==0) utility_exit_with_message("you must specify cbdock::syms");
	vector1<Real> asamp; for(Real i = 0; i < 180; ++i) asamp.push_back(i); // just 0-179
	// set up n-ca-c-o-cb ord arrays
	vector1<Vec> bb0tmp,cb0tmp;
	for(int ir = 1; ir <= init.n_residue(); ++ir) {
		if(!init.residue(ir).is_protein()) continue;
		for(int ia = 1; ia <= ((init.residue(ir).has("CB"))?5:4); ++ia) {
			bb0tmp.push_back(init.xyz(AtomID(ia,ir)));
		}
		if(init.secstruct(ir)=='H') {
			if(init.residue(ir).has("CB")) cb0tmp.push_back(init.xyz(AtomID(5,ir)));
			else													 cb0tmp.push_back(init.xyz(AtomID(4,ir)));
		}
	}
	vector1<Vec> const bb0(bb0tmp);
	vector1<Vec> const cb0(cb0tmp);

	vector1<Matf> Rsym(syms.size());
	for(Size ic = 1; ic <= syms.size(); ++ic) Rsym[ic] = rotation_matrix_degrees(Vec(1,0,0),360.0/(Real)syms[ic]);
	vector1<vector1<Hit> > hits(syms.size());

#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
	for(int iss = 1; iss <= ssamp.size(); ++iss) {
		if(iss%100==0) {
			TR << iss << " of " << ssamp.size() << " nhits:";
			for(Size ic = 1; ic <= hits.size(); ++ic) TR << " " << lss(hits[ic].size(),4);
			TR << endl;
		}
		Vec axs = ssamp[iss];
		for(int irt = 1; irt <= asamp.size(); ++irt) {
			//if(axs.x() < 0.0) continue;
			//Real const u = uniform();
			Real const rot = asamp[irt];
			Matf const R = rotation_matrix_degrees(axs,rot);

			vector1<Vec> bb1 = bb0;
			vector1<Vec> cb1 = cb0;
			for(vector1<Vec>::iterator i = bb1.begin(); i != bb1.end(); ++i) *i = R*(*i);
			for(vector1<Vec>::iterator i = cb1.begin(); i != cb1.end(); ++i) *i = R*(*i);

			for(int ic = 1; ic <= syms.size(); ic++) {
				vector1<Vec> bb2 = bb1;
				vector1<Vec> cb2 = cb1;
				for(vector1<Vec>::iterator i = bb2.begin(); i != bb2.end(); ++i) *i = Rsym[ic]*(*i);
				for(vector1<Vec>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = Rsym[ic]*(*i);
				Real cbc;
				Real t = sicfast(bb1,bb2,cb1,cb2,cbc);
				
				if(cbc >= basic::options::option[basic::options::OptionKeys::cxdock::num_contacts]()) {
// 					Hit h(iss,irt,cbc,syms[ic]);
// 					h.s1.from_four_points(bb1[1],bb1[1],bb1[2],bb1[3]);
// 					h.s2.from_four_points(bb2[1],bb2[1],bb2[2],bb2[3]);
// 					h.s1.v += t*Vec(0,0,1);
// #ifdef USE_OPENMP
// #pragma omp critical
// #endif
// 					hits[ic].push_back(h);

					if( option[OptionKeys::out::file::o].user() || option[OptionKeys::cxdock::dumpfirst] ) {	
						string tag = utility::file_basename(fn)+"_C"+str(syms[ic])+"_"+str(iss)+"_"
						+str(basic::options::option[basic::options::OptionKeys::cxdock::sphere]())+"_"+str(irt)+"_"+str(cbc);
						{
							option[OptionKeys::symmetry::symmetry_definition]("input/sym/C"+str(syms[ic])+"_Z.sym");							
						  Pose p(init);
							Real basehb = shb->score(p);
							int basenhb = num_hbonds(p);
						  rot_pose(p,ssamp[iss],asamp[irt]);
						  trans_pose(p,Vec(0,0,t));
						  Vec cen = Vec(0,t/2.0/tan( numeric::constants::d::pi / (Real)syms[ic] ),t/2.0);
						  trans_pose(p,-cen);
							rot_pose(p,Vec(0,1,0),90.0); // align sym Z
						  //core::pose::symmetry::make_symmetric_pose(p);
							
							Mat R180z = rotation_matrix_degrees(Vec(0.0,0.0,1.0),180.0);
							Mat R180x = rotation_matrix_degrees(Vec(1.0,0.0,0.0),180.0);
							Mat R90y  = rotation_matrix_degrees(Vec(0.0,1.0,0.0), 90.0);
							Mat R1z   = rotation_matrix_degrees(Vec(0.0,0.0,1.0),  2.0);

							// make C2 coords							
							vector1<Vec> c2bb1,c2cb1;
							for(int ir = 1; ir <= init.n_residue(); ++ir) {
								if(!p.residue(ir).is_protein()) continue;
								for(int ia = 1; ia <= ((p.residue(ir).has("CB"))?5:4); ++ia) {
									c2bb1.push_back(          p.xyz(AtomID(ia,ir))  );
									c2bb1.push_back(  R180z * p.xyz(AtomID(ia,ir))  );									
								}
								if(p.secstruct(ir)=='H') {
									Size i = p.residue(ir).has("CB") ? 5 : 4;
									c2cb1.push_back(          p.xyz(AtomID(i,ir))  );
									c2cb1.push_back(  R180z * p.xyz(AtomID(i,ir))  );									
								}
							}
							
							// make partner c2 coords
							vector1<Vec> c2bb2(c2bb1),c2cb2(c2cb1);
							for(vector1<Vec>::iterator i = c2bb2.begin(); i != c2bb2.end(); ++i) *i = R180x*(*i);
							for(vector1<Vec>::iterator i = c2cb2.begin(); i != c2cb2.end(); ++i) *i = R180x*(*i);

							for(Size if4 = 0; if4 < 2; ++if4) {
								for(vector1<Vec>::iterator i = c2bb1.begin(); i != c2bb1.end(); ++i) *i = R1z*(*i);
								for(vector1<Vec>::iterator i = c2cb1.begin(); i != c2cb1.end(); ++i) *i = R1z*(*i);								
								for(vector1<Vec>::iterator i = c2bb2.begin(); i != c2bb2.end(); ++i) *i = R1z*(*i);
								for(vector1<Vec>::iterator i = c2cb2.begin(); i != c2cb2.end(); ++i) *i = R1z*(*i);								
								rot_pose(p,R180x);

								for(Size ir4 = 0; ir4 < 90; ir4+=2) {
									for(vector1<Vec>::iterator i = c2bb2.begin(); i != c2bb2.end(); ++i) *i = R1z*(*i);
									for(vector1<Vec>::iterator i = c2cb2.begin(); i != c2cb2.end(); ++i) *i = R1z*(*i);
									rot_pose(p,R1z);

									Real cbc2;
									Real t = sicfast(c2bb2,c2bb2,c2cb2,c2cb2,cbc2);

									if(cbc2 >= basic::options::option[basic::options::OptionKeys::cxdock::num_contacts]()) {
										Pose tmp(p);
										trans_pose(tmp ,Vec(0,0,t/2));
										option[OptionKeys::symmetry::symmetry_definition]("input/sym/D2.sym");
										core::pose::symmetry::make_symmetric_pose(tmp);
										Pose tmp2(tmp);
										core::pose::symmetry::make_asymmetric_pose(tmp2);
										string tag2 = tag + "_"+str(if4)+"_"+str(ir4);
										Real hbsc = shb->score(tmp2) - basehb*4.0;
										Real nhb = (Real(num_hbonds(tmp2))-4*Real(basenhb))/2.0;
										
										if( nhb > 3 || cbc2/2 > basic::options::option[basic::options::OptionKeys::cxdock::num_contacts]() ) {
											cout << "HIT " << cbc << " " << cbc2/2 << " " << F(7,4,hbsc) << " " << F(3,1,nhb) << " "
											     << utility::file_basename(fn)+" C"+str(syms[ic])+" "+str(iss)+" "
											      +str(basic::options::option[basic::options::OptionKeys::cxdock::sphere]())
											      +" "+str(irt)+" "+str(cbc)+" "+str(if4)+" "+str(ir4) << std::endl;
											tmp.dump_pdb(option[OptionKeys::out::file::o]()+"/"+tag2+".pdb");
										}
//										utility_exit_with_message("oraiestn");								
										if(option[OptionKeys::cxdock::dumpfirst]()) utility_exit_with_message("dumpfirst");
									}
								}
							}
						}
					}
				}
			}
		}
		//break;
	}

	for(int ic = 1; ic <= syms.size(); ic++) {
		std::sort(hits[ic].begin(),hits[ic].end(),cmp);
		for(int i = 1; i <= min((Size)10,hits[ic].size()); ++i) {
			Hit & h(hits[ic][i]);
			cout << "RESULT " << fn << " " << h.sym << " " << h.iss << " " << basic::options::option[basic::options::OptionKeys::cxdock::sphere]()	<< " " << h.irt << " " << h.cbc << endl;
		}
	}

}


int main(int argc, char *argv[]) {
	register_options();
	core::init(argc,argv);
	using namespace basic::options::OptionKeys;


	Size const NSS = basic::options::option[basic::options::OptionKeys::cxdock::sphere]();
	vector1<Vec> ssamp(NSS);
	{
		izstream is;
		basic::database::open(is,"sampling/spheres/sphere_"+str(NSS)+".dat");
		for(int i = 1; i <= NSS; ++i) {
			Real x,y,z;
			is >> x >> y >> z;
			ssamp[i] = Vec(x,y,z);
		}
		is.close();
	}

	for(Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn) {
		string fn = option[in::file::s]()[ifn];
		Pose pnat;
		TR << "searching " << fn << std::endl;
		core::import_pose::pose_from_pdb(pnat,fn);
		trans_pose(pnat,-center_of_geom(pnat,1,pnat.n_residue()));
		core::scoring::dssp::Dssp dssp(pnat);
		dssp.insert_ss_into_pose(pnat);
		//if( pnat.n_residue() > 150 ) continue;
		Size cyscnt=0, nhelix=0;
		for(Size ir = 2; ir <= pnat.n_residue()-1; ++ir) {
			if(pnat.secstruct(ir) == 'H') nhelix++;
			//if(!pnat.residue(ir).is_protein()) goto cont1;
			if(pnat.residue(ir).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
			if(pnat.residue(ir).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
			// if(pnat.residue(ir).name3()=="CYS") { if(++cyscnt > 3) goto cont1; }
		}// goto done1; cont1: TR << "skipping " << fn << std::endl; continue; done1:
		// if( nhelix < 20 ) continue;
		Pose pala(pnat);
		dock(pala,fn,ssamp);
	}
}

