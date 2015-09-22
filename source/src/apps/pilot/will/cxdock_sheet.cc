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
OPT_1GRP_KEY( IntegerVector, cxdock, syms         )
OPT_1GRP_KEY( Real         , cxdock, clash_dis    )
OPT_1GRP_KEY( Real         , cxdock, contact_dis  )
OPT_1GRP_KEY( Real         , cxdock, hb_dis  )
OPT_1GRP_KEY( Real         , cxdock, num_contacts )
OPT_1GRP_KEY( Boolean      , cxdock, dumpfirst )
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::s );
	NEW_OPT( cxdock::syms	 ,"CX symmitries", utility::vector1< Size >() );
	NEW_OPT( cxdock::clash_dis   ,"min acceptable contact dis", 3.5 );
	NEW_OPT( cxdock::contact_dis ,"max acceptable contact dis", 7.0 );
	NEW_OPT( cxdock::hb_dis      ,"max acceptable hb dis", 4.7 );	 // very loose!
	NEW_OPT( cxdock::num_contacts ,"required no. contacts", 30.0 );
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


static THREAD_LOCAL basic::Tracer TR( "pentcb" );
static core::io::silent::SilentFileData sfd;

#include <apps/pilot/will/sicfast.ihh>


void get_avail_don_acc(Pose & pose, std::set<Size> & actdon, std::set<Size> & actacc) {
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

	Real cbcth = basic::options::option[basic::options::OptionKeys::cxdock::num_contacts]();

	std::set<Size> actdon,actacc; {
		get_avail_don_acc(init,actdon,actacc);
	}
	vector1<int> availdon,availacc;

	vector1<Size> syms = basic::options::option[basic::options::OptionKeys::cxdock::syms]();
	if(syms.size()==0) utility_exit_with_message("you must specify cbdock::syms");
	vector1<Real> asamp; for(Real i = 0.0; i < 180.0; i+=1.5) asamp.push_back(i); // just 0-179
	// set up n-ca-c-o-cb ord arrays
	vector1<Vec> bb0tmp,cb0tmp;
	for(int ir = 1; ir <= init.n_residue(); ++ir) {
		if(!init.residue(ir).is_protein()) continue;
		for(int ia = 1; ia <= ((init.residue(ir).has("CB"))?5:4); ++ia) {
			bb0tmp.push_back(init.xyz(AtomID(ia,ir)));
			if( ir > 1 && ir < init.n_residue() && ia==1 && actdon.find(ir)==actdon.end() &&
			    (init.secstruct(ir-1)=='E'&&actdon.find(ir-1)==actdon.end()||
					 init.secstruct(ir  )=='E'&&actdon.find(ir  )==actdon.end()||
					 init.secstruct(ir+1)=='E'&&actdon.find(ir+1)==actdon.end())
			) availdon.push_back(bb0tmp.size());
			if( ir > 1 && ir < init.n_residue() && ia==4 && actacc.find(ir)==actacc.end() &&
		    	(init.secstruct(ir-1)=='E'&&actacc.find(ir-1)==actacc.end()||
				 	 init.secstruct(ir  )=='E'&&actacc.find(ir  )==actacc.end()||
				 	 init.secstruct(ir+1)=='E'&&actacc.find(ir+1)==actacc.end())
			) availacc.push_back(bb0tmp.size());
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

	for(int iss = 1; iss <= ssamp.size(); ++iss) {
		if(iss%1000==0) {
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
				Real cbc, nhb;
				Real t = sicfast(bb1,bb2,cb1,cb2,availdon,availacc,cbc,nhb);

				if(nhb > 3) {
				// if(nhb>0.0&&cbc>=cbcth*1.0||nhb>1.0&&cbc>=cbcth*0.9||nhb>2.0&&cbc>=cbcth*0.8||nhb>3.0&&cbc>=cbcth*0.7||
				//    nhb>4.0&&cbc>=cbcth*0.6||nhb>5.0&&cbc>=cbcth*0.5||nhb>6.0&&cbc>=cbcth*0.4||nhb>7.0&&cbc>=cbcth*0.3||
				//    nhb>8.0&&cbc>=cbcth*0.2||nhb>9.0&&cbc>=cbcth*0.1||nhb>10.0&&cbc>=cbcth*0.0  ){
// 					Hit h(iss,irt,cbc,syms[ic]);
// 					h.s1.from_four_points(bb1[1],bb1[1],bb1[2],bb1[3]);
// 					h.s2.from_four_points(bb2[1],bb2[1],bb2[2],bb2[3]);
// 					h.s1.v += t*Vec(0,0,1);
// #ifdef USE_OPENMP
// #pragma omp critical
// #endif
// 					hits[ic].push_back(h);

					if( option[OptionKeys::out::file::o].user() || option[OptionKeys::cxdock::dumpfirst] ) {
						string tag = utility::file_basename(fn)+"_C"+str(syms[ic])+"_"+str(iss)+"_"+str(basic::options::option[basic::options::OptionKeys::cxdock::sphere]())+"_"+str(irt)+"_"+str(cbc);
						{

							cout << "HIT " << nhb << " " << cbc << " " << tag << std::endl;

							option[OptionKeys::symmetry::symmetry_definition]("input/sym/C"+str(syms[ic])+".sym");
						  Pose p(init);
						  rot_pose(p,ssamp[iss],asamp[irt]);
						  trans_pose(p,Vec(0,0,t));
						  Vec cen = Vec(0,t/2.0/tan( numeric::constants::d::pi / (Real)syms[ic] ),t/2.0);
						  trans_pose(p,-cen);
							rot_pose(p,Vec(0,1,0),90.0); // align sym Z
							p.dump_pdb(option[OptionKeys::out::file::o]+"/"+tag+"_aln_mono.pdb.gz");

							core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
						  ss_out->fill_struct(p,option[OptionKeys::out::file::o]+"/"+tag+"_aln_mono.pdb.gz");
						  ss_out->add_energy("sym",syms[ic]);
						  ss_out->add_energy("num_contact",cbc);
						  ss_out->add_energy("num_res",p.n_residue());
						  sfd.write_silent_struct( *ss_out, option[OptionKeys::out::file::o]() + "/" + option[ OptionKeys::out::file::silent ]() );

						}
						// {
						// 	Pose p(init),q(init);
						// 	core::kinematics::Stub s(init.xyz(AtomID(1,1)),init.xyz(AtomID(2,1)),init.xyz(AtomID(3,1)));
						// 	xform_pose_rev(p,s); xform_pose(p,h.s1);
						// 	xform_pose_rev(q,s); xform_pose(q,h.s2);
						// 	p.dump_pdb(option[OptionKeys::out::file::o]+"/"+tag+"_A.pdb.gz");
						// 	q.dump_pdb(option[OptionKeys::out::file::o]+"/"+tag+"_B.pdb.gz");
						// }
						if(option[OptionKeys::cxdock::dumpfirst]()) utility_exit_with_message("dumpfirst");
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

	Size NTOP = 6;
	for(Size i = 1; i <= ssamp.size(); ++i) {
		vector1<Real> mnv(NTOP,9e9);
		vector1<Size> mni(NTOP,0);
		for(Size j = 1; j <= ssamp.size(); ++j) {
			if(i==j) continue;
			for(Size k = 1; k <= NTOP; k++) {
				Real d = ssamp[i].distance_squared(ssamp[j]);
				if(mnv[k] > d) {
					for(Size l = NTOP; l > k; --l) {
						mnv[l] = mnv[l-1];
						mni[l] = mni[l-1];
					}
					mnv[k] = d;
					mni[k] = j;
					break;
				}
			}
		}
		cerr << i;
		for(Size j = 1; j <= NTOP; ++j) {
			if( j > 1 && ssamp[i].distance(ssamp[mni[j]]) > 1.5*ssamp[i].distance(ssamp[mni[j-1]]) ) {
				cerr << " 0";
			} else {
				cerr << " " << mni[j];
			}
		}
		cerr << endl;
	}

	utility_exit_with_message("oarstne");

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
		} goto done1; cont1: TR << "skipping " << fn << std::endl; continue; done1:
		// if( nhelix < 20 ) continue;
		Pose pala(pnat);
		dock(pala,fn,ssamp);
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

