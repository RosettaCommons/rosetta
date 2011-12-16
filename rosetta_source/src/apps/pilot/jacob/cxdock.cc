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
OPT_1GRP_KEY( Real         , cxdock, ang_samp )
OPT_1GRP_KEY( Boolean      , cxdock, dumpfirst )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::s );
	NEW_OPT( cxdock::syms	 ,"CX symmitries", utility::vector1< Size >() );
	NEW_OPT( cxdock::clash_dis   ,"min acceptable contact dis", 3.6 );
	NEW_OPT( cxdock::contact_dis ,"max acceptable contact dis", 6.0 );	
	NEW_OPT( cxdock::hb_dis       ,"max acceptable hb dis", 6.0 );
	NEW_OPT( cxdock::num_contacts ,"required no. contacts", 20.0 );
	NEW_OPT( cxdock::ang_samp     ,"num degrees", 1.0 );
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


static basic::Tracer TR("cxdock");
static core::io::silent::SilentFileData sfd;

#include <apps/pilot/will/sicfast.ihh>


struct Hit {
	Real cbc;
	int iss,irt,sym;
	Hit(Real _cbc, int _iss, int _irt, int _sym) {
		cbc = _cbc;
		iss = _iss;
		irt = _irt;
		sym = _sym;
	}
};


void dock(Pose const init, std::string const & fn, vector1<Vec> const & ssamp, ObjexxFCL::FArray2D<int> const & snbr) {
	using namespace basic::options;
	
	vector1<Size> syms = basic::options::option[ basic::options::OptionKeys::cxdock::syms ]();
	if(syms.size()==0) utility_exit_with_message("you must specify cbdock::syms");
	vector1<Real> asamp;
	for(Real i = -5.0; i <= 180/option[OptionKeys::cxdock::ang_samp]()+5.0; ++i  ) {
		asamp.push_back( i*option[OptionKeys::cxdock::ang_samp]() ); // just 0-179
	}
	// set up n-ca-c-o-cb ord arrays
	vector1<Vec> bb0tmp,cb0tmp;
	for(int ir = 1; ir <= init.n_residue(); ++ir) {
		if(!init.residue(ir).is_protein()) continue;
		for(int ia = 1; ia <= ((init.residue(ir).has("CB"))?5:4); ++ia) {
			bb0tmp.push_back(init.xyz(AtomID(ia,ir)));
		}
		if( init.secstruct(ir)=='H' || init.secstruct(ir)=='E' ) {
			if(init.residue(ir).has("CB")) cb0tmp.push_back(init.xyz(AtomID(5,ir)));
			else													 cb0tmp.push_back(init.xyz(AtomID(4,ir)));
		}
	}
	vector1<Vec> const bb0(bb0tmp);
	vector1<Vec> const cb0(cb0tmp);

	vector1<ObjexxFCL::FArray2D<float> > grids( syms.size(), ObjexxFCL::FArray2D<float>(asamp.size(),ssamp.size(),0.0f)  );

	vector1<Matf> Rsym(syms.size());
	for(Size ic = 1; ic <= syms.size(); ++ic) Rsym[ic] = rotation_matrix_degrees(Vec(1,0,0),360.0/(Real)syms[ic]);

	for(int iss = 1; iss <= ssamp.size(); ++iss) {
		if(iss%100==0) {
			TR << iss << " of " << ssamp.size();
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

			for(int ic = 1; ic <= syms.size(); ic++) { // loop over syms
				vector1<Vec> bb2 = bb1;
				vector1<Vec> cb2 = cb1;
				for(vector1<Vec>::iterator i = bb2.begin(); i != bb2.end(); ++i) *i = Rsym[ic]*( *i ); // *i is a Vec
				for(vector1<Vec>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = Rsym[ic]*( *i );
				Real cbc;
				Real t = sicfast(bb1,bb2,cb1,cb2,cbc);
				grids[ic](irt,iss) = cbc;
				
				// if(cbc >= basic::options::option[basic::options::OptionKeys::cxdock::num_contacts]()) {
				// 
				// 	if( option[OptionKeys::out::file::o].user() || option[OptionKeys::cxdock::dumpfirst] ) {	
				// 		/// this is the filename
				// 		string tag =  utility::file_basename(fn)	+"_C"+str(syms[ic])		+"_"+str(iss)
				// 			+"_"	+str(basic::options::option[basic::options::OptionKeys::cxdock::sphere]())
				// 			+"_"	+str(irt)		+"_"		+str(cbc);
				// 		{
				// 			option[OptionKeys::symmetry::symmetry_definition]("input/sym/C"+str(syms[ic])+".sym");							
				// 		  Pose p(init);
				// 		  rot_pose(p,ssamp[iss],asamp[irt]);
				// 		  trans_pose(p,Vec(0,0,t));
				// 		  Vec cen = Vec(0,t/2.0/tan( numeric::constants::d::pi / (Real)syms[ic] ),t/2.0);
				// 		  trans_pose(p,-cen);
				// 			rot_pose(p,Vec(0,1,0),90.0); // align sym Z
				// 		  //core::pose::symmetry::make_symmetric_pose(p);
				// 			p.dump_pdb(option[OptionKeys::out::file::o]+"/"+tag+"_aln_mono.pdb.gz");
				// 			
				// 			core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
				// 		  ss_out->fill_struct(p,option[OptionKeys::out::file::o]+"/"+tag+"_aln_mono.pdb.gz");
				// 		  ss_out->add_energy("sym",syms[ic]);
				// 		  ss_out->add_energy("num_contact",cbc);
				// 		  ss_out->add_energy("num_res",p.n_residue());
				// 		  sfd.write_silent_struct( *ss_out, option[OptionKeys::out::file::o]() + "/" + option[ OptionKeys::out::file::silent ]() );
				// 		  							
				// 		}
				// 		// {
				// 		// 	Pose p(init),q(init);
				// 		// 	core::kinematics::Stub s(init.xyz(AtomID(1,1)),init.xyz(AtomID(2,1)),init.xyz(AtomID(3,1)));
				// 		// 	xform_pose_rev(p,s); xform_pose(p,h.s1);
				// 		// 	xform_pose_rev(q,s); xform_pose(q,h.s2);
				// 		// 	p.dump_pdb(option[OptionKeys::out::file::o]+"/"+tag+"_A.pdb.gz");
				// 		// 	q.dump_pdb(option[OptionKeys::out::file::o]+"/"+tag+"_B.pdb.gz");
				// 		// }
				// 		if(option[OptionKeys::cxdock::dumpfirst]()) utility_exit_with_message("dumpfirst");
				// 	}
				// }
			}
		}
		//break;
	}
	


	for(Size ic = 1; ic <= syms.size(); ++ic) {
		Size nlocalmax = 0;
		Size NTOP = 20;
		vector1<Hit> tophits(NTOP,Hit(-9e9,0,0,0));
		
		for(Size iss = 1; iss <= ssamp.size(); ++iss) {
			for(Size irt = 6; irt <= asamp.size()-5; ++irt) {
				Real cbc = grids[ic](irt,iss);
				
				bool is_local_max = true;
				for(Size inb = 1; inb <= 6; inb++) {
					Real ncbc0 = grids[ic]( irt-1, snbr(inb,iss) );
					Real ncbc1 = grids[ic]( irt  , snbr(inb,iss) );
					Real ncbc2 = grids[ic]( irt+1, snbr(inb,iss) );										
					if( ncbc0 > cbc || ncbc1 > cbc || ncbc2 > cbc ) is_local_max = false;
				}
				Real ancbc1 = grids[ic](irt-1,iss);
				Real ancbc2 = grids[ic](irt+1,iss);
				Real ancbc3 = grids[ic](irt-2,iss);
				Real ancbc4 = grids[ic](irt+2,iss);
				if( ancbc1 > cbc || ancbc2 > cbc || ancbc3 > cbc || ancbc4 > cbc) is_local_max = false;
				
				if(is_local_max==false) continue;
				nlocalmax++;
				Hit h(cbc,iss,irt,syms[ic]);
				for(Size k = 1; k <= NTOP; k++) {
					if(tophits[k].cbc < h.cbc) {
						for(Size l = NTOP; l > k; --l) {
							tophits[l] = tophits[l-1];
						}
						tophits[k] = h;
						break;
					}
				}
			}
		}

		TR << "N LOCAL MAX " << syms[ic] << ": " << nlocalmax << " of " << ssamp.size()*asamp.size() << endl;

		for(Size k = 1; k <= NTOP; k++) {
			Hit h = tophits[k];
			TR << "TOP HIT " << h.sym << " " << h.cbc << " " << h.iss << " " << h.irt << " " << ssamp.size() << " " << option[OptionKeys::cxdock::ang_samp]() << endl;
			
			Matf const R = rotation_matrix_degrees(ssamp[h.iss],asamp[h.irt]);
			vector1<Vec> bb1 = bb0;
			vector1<Vec> cb1 = cb0;
			for(vector1<Vec>::iterator i = bb1.begin(); i != bb1.end(); ++i) *i = R*(*i);
			for(vector1<Vec>::iterator i = cb1.begin(); i != cb1.end(); ++i) *i = R*(*i);
			vector1<Vec> bb2 = bb1;
			vector1<Vec> cb2 = cb1;
			for(vector1<Vec>::iterator i = bb2.begin(); i != bb2.end(); ++i) *i = Rsym[ic]*( *i ); // *i is a Vec
			for(vector1<Vec>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = Rsym[ic]*( *i );
			Real tmp;
			Real t = sicfast(bb1,bb2,cb1,cb2,tmp);
			
		  Pose p(init);
		  rot_pose(p,ssamp[h.iss],asamp[h.irt]);
		  trans_pose(p,Vec(0,0,t));
		  Vec cen = Vec(0,t/2.0/tan( numeric::constants::d::pi / (Real)syms[ic] ),t/2.0);
		  trans_pose(p,-cen);
			rot_pose(p,Vec(0,1,0),90.0); // align sym Z
			p.dump_pdb(option[OptionKeys::out::file::o]+"/hit_"+str(k)+".pdb.gz");

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
	ObjexxFCL::FArray2D<int> snbr(6,NSS);
	{
		izstream is;
		basic::database::open(is,"sampling/spheres/sphere_"+str(NSS)+"_neighbors.dat");
		for(int i = 1; i <= NSS; ++i) {
			int tmp;
			is >> tmp >> snbr(1,i) >> snbr(2,i) >> snbr(3,i) >> snbr(4,i) >> snbr(5,i) >> snbr(6,i);
			;
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
		} goto done1; cont1: TR << "skipping " << fn << std::endl; continue; done1:
		// if( nhelix < 20 ) continue;
		Pose pala(pnat);
		dock(pala,fn,ssamp,snbr);
	}
}

