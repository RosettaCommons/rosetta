#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/annotated_sequence.hh>
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
OPT_1GRP_KEY( Real         , cxdock, hb_dis    )
OPT_1GRP_KEY( Real         , cxdock, num_contacts )
OPT_1GRP_KEY( Boolean      , cxdock, dumpfirst )
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::s );
//	NEW_OPT( cxdock::syms	 ,"CX symmitries", utility::vector1< Size >() );
	NEW_OPT( cxdock::clash_dis   ,"min acceptable contact dis", 3.5 );
	NEW_OPT( cxdock::contact_dis ,"max acceptable contact dis", 7.0 );	
	NEW_OPT( cxdock::hb_dis      ,"max acceptable hb dis", 4.7 );	
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
using core::import_pose::pose_from_pdb;
using core::kinematics::Stub;
using core::conformation::ResidueOP;


static basic::Tracer TR("dxdock");
static core::io::silent::SilentFileData sfd;

#include <apps/pilot/will/sicfast.ihh>

void dump_points_pdb(vector1<Vec> & p, string fn, vector1<int> const & d, vector1<int> const & a) {
	// for(Size i = 1; i <= d.size(); ++i) cout << "DON " << d[i] << endl;
	// for(Size i = 1; i <= a.size(); ++i) cout << "ACC " << a[i] << endl;	
	std::ofstream o(fn.c_str());
	for(Size i = 1; i <= p.size(); ++i) {
		string rn = "VIZ";
		if(std::find(d.begin(),d.end(),i)!=d.end()) rn="DON";
		if(std::find(a.begin(),a.end(),i)!=a.end()) rn="ACC";		
		o<<"HETATM"<<I(5,i)<<' '<<" CA "<<' '<<rn<<' '<<"A"<<I(4,i)<<"    "<<F(8,3,p[i].x())<<F(8,3,p[i].y())<<F(8,3,p[i].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	}
	o.close();
}



void get_avail_don_acc(Pose & pose, std::set<Size> & actdon, std::set<Size> & actacc) {
	core::scoring::dssp::Dssp dssp(pose);
	dssp.insert_ss_into_pose(pose);
	ScoreFunctionOP shb = new core::scoring::ScoreFunction;
	shb->set_weight(core::scoring::hbond_lr_bb,1.0);
	shb->set_weight(core::scoring::hbond_sr_bb,1.0);
	shb->score(pose);
	core::scoring::hbonds::HBondSet hbset;
 	core::scoring::hbonds::fill_hbond_set( pose, false, hbset, false, true, true, true );
 	for(Size i = 1; i <= hbset.nhbonds(); ++i) {
		actdon.insert( hbset.hbond(i).don_res() );
		actacc.insert( hbset.hbond(i).acc_res() );		
	}	
}


struct Hit {
	int iss,irt,cbc,sym;
	core::kinematics::Stub s1,s2;
	Hit(int is, int ir, int cb, int sm) : iss(is),irt(ir),cbc(cb),sym(sm) {}
};
bool cmp(Hit i,Hit j) { return i.cbc > j.cbc; }

void dock(Pose & init, std::string const & fn, vector1<Vec> const & ssamp) {
	using namespace basic::options;
	
	std::set<Size> actdon,actacc; {
		get_avail_don_acc(init,actdon,actacc);
	}
	vector1<int> availdon,availacc;
	
	Real cbcth = basic::options::option[basic::options::OptionKeys::cxdock::num_contacts]();
	
	// Pose ala;
	// core::pose::make_pose_from_sequence(ala,"A",init.residue(1).residue_type_set(),false);
	// Pose pala(init);
	// for(Size i = 1; i <= pala.n_residue(); ++i) {
	// 	if(!pala.residue(i).is_protein()) continue;		
	// 	if(pala.residue(i).aa()!=core::chemical::aa_gly) {
	// 		pala.replace_residue(i,ala.residue(1),true);
	// 	}
	// }
	
	
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
			if( ir > 1 && ir < init.n_residue() && ia==1 && actacc.find(ir)==actacc.end() &&
			    (init.secstruct(ir-1)=='E'||init.secstruct(ir)=='E'||init.secstruct(ir+1)=='E') ) availdon.push_back(bb0tmp.size());
			if( ir > 1 && ir < init.n_residue() && ia==4 && actacc.find(ir)==actacc.end() &&
		    	(init.secstruct(ir-1)=='E'||init.secstruct(ir)=='E'||init.secstruct(ir+1)=='E')	) availacc.push_back(bb0tmp.size());
		}
		if(init.secstruct(ir)=='H') {
			if(init.residue(ir).has("CB")) cb0tmp.push_back(init.xyz(AtomID(5,ir)));
			else													 cb0tmp.push_back(init.xyz(AtomID(4,ir)));
		}
	}
	vector1<Vec> const & bb0(bb0tmp);
	vector1<Vec> const & cb0(cb0tmp);
		
	vector1<int> c2availdon(availdon),c2availacc(availacc);
	for(Size i=1,s=availdon.size(); i <= s; ++i) c2availdon.push_back(availdon[i]+bb0.size());
	for(Size i=1,s=availacc.size(); i <= s; ++i) c2availacc.push_back(availacc[i]+bb0.size());	

	vector1<Matf> Rsym(syms.size());
	for(Size ic = 1; ic <= syms.size(); ++ic) Rsym[ic] = rotation_matrix_degrees(Vec(1,0,0),360.0/(Real)syms[ic]);
	vector1<vector1<Hit> > hits(syms.size());

#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
	for(int iss = 1; iss <= ssamp.size(); ++iss) {
		if(iss%500==0) {
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
				Real cbc,nhb;
				Real t = sicfast(bb1,bb2,cb1,cb2,availdon,availacc,cbc,nhb);
				
				if(          cbc>=cbcth*1.0||nhb> 0.0&&cbc>=cbcth*0.9||nhb> 2.0&&cbc>=cbcth*0.8||nhb> 4.0&&cbc>=cbcth*0.5||
				   nhb> 6.0&&cbc>=cbcth*0.3||nhb> 8.0      )//&&cbc>=cbcth*0.5||nhb>10.0&&cbc>=cbcth*0.4||nhb>12.0&&cbc>=cbcth*0.3||
//			   nhb>14.0&&cbc>=cbcth*0.2||nhb>16.0&&cbc>=cbcth*0.1||nhb>18.0&&cbc>=cbcth*0.0  )
				{
					// if( mindis(bb1,bb2,t) < basic::options::option[basic::options::OptionKeys::cxdock::clash_dis]()-0.3	) {
					// 		cerr << "clash fail!" << endl;
					// 		continue;
					// }
					if( option[OptionKeys::out::file::o].user() || option[OptionKeys::cxdock::dumpfirst] ) {	
						string tag = utility::file_basename(fn)+"_C"+str(syms[ic])+"_"+str(iss)+"_"
						+str(basic::options::option[basic::options::OptionKeys::cxdock::sphere]())+"_"+str(irt)+"_"+str(cbc);
						{
							option[OptionKeys::symmetry::symmetry_definition]("input/sym/C"+str(syms[ic])+"_Z.sym");							
						  Pose p(init);
						  rot_pose(p,ssamp[iss],asamp[irt]);
						  trans_pose(p,Vec(0,0,t));
						  Vec cen = Vec(0,t/2.0/tan( numeric::constants::d::pi / (Real)syms[ic] ),t/2.0);
						  trans_pose(p,-cen);
							rot_pose(p,Vec(0,1,0),90.0); // align sym Z

							// for(Size i = 1; i <= bb1.size(); ++i) bb1[i] += Vec(0,0,t);
							// dump_points_pdb(bb1,"test1.pdb",availdon,availacc);
							// dump_points_pdb(bb2,"test2.pdb",availdon,availacc);
							// cout << cbc << " " << nhb << endl;
							// core::pose::symmetry::make_symmetric_pose(p);
							// p.dump_pdb("test.pdb");
							// utility_exit_with_message("arst");
							
							Mat R180z = rotation_matrix_degrees(Vec(0.0,0.0,1.0),180.0);
							Mat R180x = rotation_matrix_degrees(Vec(1.0,0.0,0.0),180.0);
							Mat R90y  = rotation_matrix_degrees(Vec(0.0,1.0,0.0), 90.0);
							Mat R2z   = rotation_matrix_degrees(Vec(0.0,0.0,1.0),  2.0);

							// make C2 coords							
							vector1<Vec> c2bb1,c2cb1;
							for(int ir = 1; ir <= init.n_residue(); ++ir) {
								if(!p.residue(ir).is_protein()) continue;
								for(int ia = 1; ia <= ((p.residue(ir).has("CB"))?5:4); ++ia) {
									c2bb1.push_back(          p.xyz(AtomID(ia,ir))  );
								}
								if(p.secstruct(ir)=='H') {
									Size i = p.residue(ir).has("CB") ? 5 : 4;
									c2cb1.push_back(          p.xyz(AtomID(i,ir))  );
								}
							}
							for(int ir = 1; ir <= init.n_residue(); ++ir) {
								if(!p.residue(ir).is_protein()) continue;
								for(int ia = 1; ia <= ((p.residue(ir).has("CB"))?5:4); ++ia) {
									c2bb1.push_back(  R180z * p.xyz(AtomID(ia,ir))  );			
								}
								if(p.secstruct(ir)=='H') {
									Size i = p.residue(ir).has("CB") ? 5 : 4;
									c2cb1.push_back(  R180z * p.xyz(AtomID(i,ir))  );									
								}
							}
							
							for(Size if4 = 0; if4 < 2; ++if4) {
								for(vector1<Vec>::iterator i = c2bb1.begin(); i != c2bb1.end(); ++i) *i = R180x*(*i);
								for(vector1<Vec>::iterator i = c2cb1.begin(); i != c2cb1.end(); ++i) *i = R180x*(*i);								
								rot_pose(p,                                                               R180x);

								// make partner c2 coords
								vector1<Vec> c2bb2(c2bb1),c2cb2(c2cb1);
								for(vector1<Vec>::iterator i = c2bb2.begin(); i != c2bb2.end(); ++i) *i = R180x*(*i);
								for(vector1<Vec>::iterator i = c2cb2.begin(); i != c2cb2.end(); ++i) *i = R180x*(*i);

								for(Size ir4 = 2; ir4 <= 90; ir4+=2) {
									for(vector1<Vec>::iterator i = c2bb2.begin(); i != c2bb2.end(); ++i) *i = R2z*(*i);
									for(vector1<Vec>::iterator i = c2cb2.begin(); i != c2cb2.end(); ++i) *i = R2z*(*i);

									Real cbc2,nhb2;
									Real t2 = sicfast(c2bb1,c2bb2,c2cb1,c2cb2,c2availdon,c2availacc,cbc2,nhb2);
									cbc2 *= 0.5;
									nhb2 /= 2.0;
									
									if(           cbc2>=cbcth*1.0||nhb2> 0.0&&cbc2>=cbcth*0.9||nhb2> 2.0&&cbc2>=cbcth*0.8||nhb2> 4.0&&cbc2>=cbcth*0.5||
									   nhb2> 6.0&&cbc2>=cbcth*0.3||nhb2> 8.0      )//&&cbc2>=cbcth*0.5||nhb2>10.0&&cbc2>=cbcth*0.4||nhb2>12.0&&cbc2>=cbcth*0.3||
//									   nhb2>14.0&&cbc2>=cbcth*0.2||nhb2>16.0&&cbc2>=cbcth*0.1||nhb2>18.0&&cbc2>=cbcth*0.0  )
									{
										// if( mindis(c2bb1,c2bb2,t2) < basic::options::option[basic::options::OptionKeys::cxdock::clash_dis]()-0.3	) {
										// 	cerr << "clash fail!" << endl;
										// 	continue;
										// }

										Pose tmp(p);
										rot_pose(tmp,Vec(0,0,1),-Real(ir4)/2.0); // half
										trans_pose(tmp ,Vec(0,0,t2/2));
										option[OptionKeys::symmetry::symmetry_definition]("input/sym/D2.sym");
										core::pose::symmetry::make_symmetric_pose(tmp);
										// Pose tmp2(tmp);
										// core::pose::symmetry::make_asymmetric_pose(tmp2);
										string tag2 = tag + "_"+str(if4)+"_"+str(ir4);

										cout << "HIT " << cbc << " " << cbc2 << " " << nhb << " " << nhb2 << " " << tag2 << std::endl;											
										tmp.dump_pdb(option[OptionKeys::out::file::o]()+"/"+tag2+".pdb.gz");

										// for(Size i = 1; i <= c2bb1.size(); ++i) c2bb1[i] += Vec(0,0, t2/2);
										// for(Size i = 1; i <= c2bb1.size(); ++i) c2bb2[i] += Vec(0,0,-t2/2);
										// dump_points_pdb(c2bb1,"test1.pdb",c2availdon,c2availacc);
										// dump_points_pdb(c2bb2,"test2.pdb",c2availdon,c2availacc);
										// p.dump_pdb("test.pdb");
										// utility_exit_with_message("aortn");

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

	try {

	register_options();
	devel::init(argc,argv);
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

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}

