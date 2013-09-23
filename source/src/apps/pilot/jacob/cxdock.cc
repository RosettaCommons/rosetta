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
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/sasa.hh>
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
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <apps/pilot/will/will_util.ihh>

#include <utility/excn/Exceptions.hh>


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
	NEW_OPT( cxdock::clash_dis   ,"min acceptable contact dis", 3.6 );
	NEW_OPT( cxdock::contact_dis ,"max acceptable contact dis", 6.0 );	
	NEW_OPT( cxdock::hb_dis       ,"max acceptable hb dis", 6.0 );
	NEW_OPT( cxdock::num_contacts ,"required no. contacts", 20.0 );
	NEW_OPT( cxdock::sphere       ,"sph points", 8192 );
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
using core::import_pose::pose_from_pdb;
using core::kinematics::Stub;
using core::conformation::ResidueOP;


static basic::Tracer TR("cxdock");
static core::io::silent::SilentFileData sfd;

#include <apps/pilot/will/sicfast.ihh>


class SphereSampler {
	vector1<Vec> all_;
	vector1<int> half_;
	vector1<int> half_w_nbr_;
	vector1<vector1<int> > nbr_;
public:
	SphereSampler(int NSS) {
		all_.resize(NSS);		{
			izstream is;
			basic::database::open(is,"sampling/spheres/sphere_"+str(NSS)+".dat");
			for(int i = 1; i <= NSS; ++i) {
				Real x,y,z;
				is >> x >> y >> z;
				all_[i] = Vec(x,y,z);
				if( all_[i].z() > 0) half_.push_back(i);
				else if( all_[i].z()==0 && all_[i].y() > 0) half_.push_back(i);
				else if( all_[i].z()==0 && all_[i].y()==0 && all_[i].x() > 0) half_.push_back(i);
			}
			is.close();
		}
		nbr_.resize(NSS);		{
			izstream is;
			basic::database::open(is,"sampling/spheres/sphere_"+str(NSS)+"_neighbors.dat");
			if(!is.good()) utility_exit_with_message("sphere neighbors not found!");
			for(int i = 1; i <= NSS; ++i) {
				int tmp;
				nbr_[i].resize(6);
				is >> tmp >> nbr_[i][1] >> nbr_[i][2] >> nbr_[i][3] >> nbr_[i][4] >> nbr_[i][5] >> nbr_[i][6];
				if( ! nbr_[i].back() ) nbr_[i].pop_back();
			}
			is.close();
		}
		// inefficient not to use set...
		half_w_nbr_ = half_;
		for(vector1<int>::const_iterator i = half_.begin(); i != half_.end(); ++i) {
			for(vector1<int>::const_iterator j = nbr_[*i].begin(); j != nbr_[*i].end(); ++j) {
				if( std::find(half_w_nbr_.begin(),half_w_nbr_.end(),*j) == half_w_nbr_.end() ) {
					half_w_nbr_.push_back(*j);
				}
			}
		}
		cout << "SphereSampler " << all_.size() << " " << half_.size() << " " << half_w_nbr_.size() << endl;
		if( half_.size() != all_.size()/2 ) utility_exit_with_message("bad half!!");
	}
	Vec const & operator[](size_t i) const { return all_[i]; }
	size_t size() const { return all_.size(); }
	size_t half_w_nbr_size() const { return half_w_nbr_.size(); }
	vector1<int> const & half_w_nbr() const { return half_w_nbr_; }
	vector1<int>::const_iterator nbr_begin(size_t const & i) const { return nbr_[i].begin(); }
	vector1<int>::const_iterator nbr_end  (size_t const & i) const { return nbr_[i].end()  ; }
	vector1<int>::const_iterator half_begin() const { return half_.begin(); }
	vector1<int>::const_iterator half_end  () const { return half_.end()  ; }
	vector1<int>::const_iterator half_w_nbr_begin() const { return half_w_nbr_.begin(); }
	vector1<int>::const_iterator half_w_nbr_end  () const { return half_w_nbr_.end()  ; }
};

void dump_points_pdb(vector1<Vec> const & p, string fn) {
	std::ofstream o(fn.c_str());
	for(Size i = 1; i <= p.size(); ++i) {
		string rn = "VIZ";
		o<<"HETATM"<<I(5,i)<<' '<<" CA "<<' '<<rn<<' '<<"A"<<I(4,i)<<"    "<<F(8,3,p[i].x())<<F(8,3,p[i].y())<<F(8,3,p[i].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	}
	o.close();
}

void dump_points_pdb_contacts(vector1<Vec> const & p, string fn) {
	std::ofstream o(fn.c_str());
	for(Size i = 1; i <= p.size(); ++i) {
		string aname = " CA ";
		if( CONTACTS.count(i) ) aname = "CONT";
		o<<"HETATM"<<I(5,i)<<' '<<aname<<' '<<"VIZ"<<' '<<"A"<<I(4,i)<<"    "<<F(8,3,p[i].x())<<F(8,3,p[i].y())<<F(8,3,p[i].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	}
	o.close();
}



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

utility::vector1<core::Real> get_ang_samp() {
	using numeric::constants::d::pi;
	vector1<Real> asamp;
	Real requested = basic::options::option[basic::options::OptionKeys::cxdock::sphere]();
	requested = sqrt( 4*pi*(180.0/pi)*(180.0/pi) / requested );
	Real angsamp = 360.0 / round(360.0/requested);
	TR << "using angsamp: " << angsamp << " (closest feasible value to requested "<<requested <<")" << endl;
	for(int i = 0; i < int(360.0/angsamp); ++i  ) {
		asamp.push_back( (Real(i))*angsamp );
	}	
	// for(Size i = 1; i <= asamp.size(); ++i) {
	// 	if( fabs( asamp[i]+angsamp/2.0 - (360.0-asamp[asamp.size()-i+1]) ) > 0.00001 ) utility_exit_with_message("angle samp issue!!!!");
	// }
	return asamp;
}


void dock(Pose const & init, std::string const & fn, SphereSampler const & ssamp) {
	using namespace basic::options;
	
	vector1<Size> syms = basic::options::option[ basic::options::OptionKeys::cxdock::syms ]();
	if(syms.size()==0) utility_exit_with_message("you must specify cbdock::syms");

	vector1<Real> asamp( get_ang_samp() );

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

	vector1<ObjexxFCL::FArray2D<float> > grids( syms.size(), ObjexxFCL::FArray2D<float>(asamp.size(),ssamp.size(),-12345.0f)  );

	vector1<Mat> Rsym(syms.size());
	for(Size ic = 1; ic <= syms.size(); ++ic) Rsym[ic] = rotation_matrix_degrees(Vec(1,0,0),360.0/(Real)syms[ic]);

	int count = 0;
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif	
	//for(vector1<int>::const_iterator issi = ssamp.half_w_nbr_begin(); issi != ssamp.half_w_nbr_end(); ++issi) {
	for(Size issi = 1; issi <= ssamp.half_w_nbr_size(); ++issi) {
		int const iss( ssamp.half_w_nbr()[issi] );
		count++; 
		if(count%100==0) {
#ifdef USE_OPENMP
#pragma omp critical
#endif									
			TR << " ..." << F(4,1, Real(count)/Real(ssamp.half_w_nbr_size())*100.0) << "\%" << endl;
		}
		Vec axs = ssamp[iss];
		for(int irt = 1; irt <= asamp.size(); ++irt) {
			//if(axs.x() < 0.0) continue;
			//Real const u = uniform();
			Real const rot = asamp[irt];
			Mat const R = rotation_matrix_degrees(axs,rot);

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
				
				//cout << I(5,iss) << " " << I(3,irt) << " " << cbc << endl;
				
				grids[ic](irt,iss) = cbc;
				//grids[ic](irt,iss) = -fabs(t);//cbc;				
			}
		}
	}
	TR << endl;



	dump_points_pdb_contacts(bb0,"contacts.pdb");




	
	for(Size ic = 1; ic <= syms.size(); ++ic) {
		Size nlocalmax = 0;
		Size NTOP = 10;
		vector1<Hit> tophits(NTOP,Hit(-9e9,0,0,0));
		
		float mx = -9e9; { 
			for(Size iss = 1; iss <= ssamp.size(); ++iss) 
				for(Size irt = 2; irt <= asamp.size()-1; ++irt) 
					mx = max(mx,grids[ic](irt,iss));
		}
		//cout << "sym " << ic << " " << syms[ic] << " cbc mx: " << mx << endl;
		
		for(vector1<int>::const_iterator issi = ssamp.half_begin(); issi != ssamp.half_end(); ++issi) {
			int const iss(*issi);
			for(Size irt = 1; irt <= asamp.size(); ++irt) {
				float cbc = grids[ic](irt,iss);
				if( cbc < 0.0 ) utility_exit_with_message("BAD CBC!!!");
				if( cbc < 1.0 ) continue;
				
				bool is_local_max = true;
				for(vector1<int>::const_iterator inb = ssamp.nbr_begin(iss); inb != ssamp.nbr_end(iss); inb++) {
					float ncbc = grids[ic]( irt, *inb );
					if( ncbc < 0.0 ) utility_exit_with_message("BAD NCBC!!!");
					if( ncbc > cbc ) is_local_max = false;
				}
				float ancbc1 = grids[ic]( irt==1 ? asamp.size() : irt-1 , iss ); if( ancbc1 < 0.0 ) utility_exit_with_message("BAD ANCBC1!!!");
				float ancbc2 = grids[ic]( irt==asamp.size() ? 1 : irt+1 , iss ); if( ancbc2 < 0.0 ) utility_exit_with_message("BAD ANCBC2!!!");
				if( ancbc1 > cbc || ancbc2 > cbc ) is_local_max = false;
				
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
			Vec a = ssamp[h.iss];
			Real ang = asamp[h.irt];
			// if( a.x() < 0.0 ) {
			// 	a = -a;
			// 	ang = 360.0 - ang;
			// }
			TR << "TOPHIT C" << h.sym << " " << F(11,7,h.cbc) << "  axs: " << F(16,13,a.x())<<","<<F(16,13,a.y())<<","<<F(16,13,a.z()) << "  ang: " << F(20,16,ang) << endl;
			
			Mat const R = rotation_matrix_degrees(ssamp[h.iss],asamp[h.irt]);
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
			
			option[OptionKeys::symmetry::symmetry_definition]("input/sym/C"+str(h.sym)+"_Z.sym");			
		  Pose p(init);
		  rot_pose(p,ssamp[h.iss],asamp[h.irt]);
		  trans_pose(p,Vec(0,0,t));
		  Vec cen = Vec(0,t/2.0/tan( numeric::constants::d::pi / (Real)syms[ic] ),t/2.0);
		  trans_pose(p,-cen);
			rot_pose(p,Vec(0,1,0),90.0); // align sym Z
			core::pose::symmetry::make_symmetric_pose(p);
			p.dump_pdb(option[OptionKeys::out::file::o]+"/hit_"+str(k)+".pdb.gz");

		}


	}



}


void testone(vector1<Vec> const & bb0,
             vector1<Vec> const & cb0,
						 int i, int j, int k
) {
	Vec axs0(-0.7591034465040, 0.1958853313150, 0.6207985941360);
	Real ang0(309.0265486725663777);
	Vec X = Vec(1,0,0).cross(axs0);
	Vec Y = X.cross(axs0);
	int N = 20;
	Real step = 0.5;
	Mat Rsym = rotation_matrix_degrees(Vec(1,0,0),120.0);

	Vec axs = axs0;
	axs = rotation_matrix_degrees(X,Real(i)*step) * axs;
	axs = rotation_matrix_degrees(Y,Real(j)*step) * axs;
	Mat R = rotation_matrix_degrees(axs,ang0+Real(k)*step);

	vector1<Vec> bb1 = bb0;
	vector1<Vec> cb1 = cb0;
	for(vector1<Vec>::iterator it = bb1.begin(); it != bb1.end(); ++it) *it = R*(*it);
	for(vector1<Vec>::iterator it = cb1.begin(); it != cb1.end(); ++it) *it = R*(*it);
	vector1<Vec> bb2 = bb1;
	vector1<Vec> cb2 = cb1;
	for(vector1<Vec>::iterator it = bb2.begin(); it != bb2.end(); ++it) *it = Rsym*(*it);
	for(vector1<Vec>::iterator it = cb2.begin(); it != cb2.end(); ++it) *it = Rsym*(*it);

	Real cbc;
	Real t = sicfast(bb1,bb2,cb1,cb2,cbc);
		
	TR << i << " " << j << " " << k << " " << cbc << " " << mindis(bb1,bb2,t) << endl;
	

	for(vector1<Vec>::iterator it = bb1.begin(); it != bb1.end(); ++it) *it += Vec(0,0,t);
	for(vector1<Vec>::iterator it = cb1.begin(); it != cb1.end(); ++it) *it += Vec(0,0,t);
	dump_points_pdb(bb1,"test"+str(k)+"A.pdb");
	dump_points_pdb(bb2,"test"+str(k)+"B.pdb");	
	//TR << i+N+1 << " " << j+N+1 << " " << k+N+1 << endl;
	
}

void visualize(Pose const & init) {
 
	//Vec axs0(-0.7475852223468591,0.2098162208303284,0.6301535438327303);
	//Real ang0(308.327); 
	// Vec axs0(-0.6719036788257763  ,   0.03696683782422633   ,   0.7397154177665476 );
	// Real ang0(297.32700000);
	Vec axs0(-0.8906229656427180  ,    0.1436574080368423  ,    0.4314548437389231 );
	Real ang0(305.83700000);
	
	Vec X = Vec(1,0,0).cross(axs0);
	Vec Y = X.cross(axs0);
	int N = 30;
	Real step = 0.001;
	
	// set up n-ca-c-o-cb ord arrays
	vector1<Vec> bb0tmp,cb0tmp;
	for(int ir = 1; ir <= init.n_residue(); ++ir) {
		if(!init.residue(ir).is_protein()) continue;
		for(int ia = 1; ia <= ((init.residue(ir).has("CB"))?5:4); ++ia) bb0tmp.push_back(init.xyz(AtomID(ia,ir)));
		if( init.secstruct(ir)=='H' || init.secstruct(ir)=='E' ) {
			if(init.residue(ir).has("CB")) cb0tmp.push_back(init.xyz(AtomID(5,ir)));
			else													 cb0tmp.push_back(init.xyz(AtomID(4,ir)));
		}
	}
	vector1<Vec> const bb0(bb0tmp);
	vector1<Vec> const cb0(cb0tmp);
	
	Mat Rsym = rotation_matrix_degrees(Vec(1,0,0),120.0);
	
	ObjexxFCL::FArray3D<float> grid(2*N+1,2*N+1,2*N+1,0.0f);
	
	// testone(bb0,cb0,12-N-1,33-N-1,11-N-1);
	// testone(bb0,cb0,12-N-1,33-N-1,11-N);
	// 
	// utility_exit_with_message("aorstn");
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
	for(int i = -N; i <= N; i++) {
		cout << i << " "; cout.flush();
		for(int j = -N; j <= N; j++) {
			for(int k = -N; k <= N; k++) {
				Vec axs = axs0;
				axs = rotation_matrix_degrees(X,Real(i)*step) * axs;
				axs = rotation_matrix_degrees(Y,Real(j)*step) * axs;
				Mat R = rotation_matrix_degrees(axs,ang0+Real(k)*step);
	
				vector1<Vec> bb1 = bb0;
				vector1<Vec> cb1 = cb0;
				for(vector1<Vec>::iterator it = bb1.begin(); it != bb1.end(); ++it) *it = R*(*it);
				for(vector1<Vec>::iterator it = cb1.begin(); it != cb1.end(); ++it) *it = R*(*it);
				vector1<Vec> bb2 = bb1;
				vector1<Vec> cb2 = cb1;
				for(vector1<Vec>::iterator it = bb2.begin(); it != bb2.end(); ++it) *it = Rsym*(*it);
				for(vector1<Vec>::iterator it = cb2.begin(); it != cb2.end(); ++it) *it = Rsym*(*it);

				Real cbc;
				Real t = sicfast(bb1,bb2,cb1,cb2,cbc);
				grid(i+N+1,j+N+1,k+N+1) = cbc;
				
			//mndis(i+N+1,j+N+1,k+N+1) = mindis(bb1,bb2,t);
				//TR << i+N+1 << " " << j+N+1 << " " << k+N+1 << endl;
			}
		}
	}
	cout << endl;
	
	// ObjexxFCL::FArray3D<float> gridsmooth(2*N+1,2*N+1,2*N+1,0.0f);
	// for(int k = 2; k < 2*N+1; k++) {
	// for(int j = 2; j < 2*N+1; j++) {
	// for(int i = 2; i < 2*N+1; i++) {
	// 	float wtot=0.0, smval=0.0;
	// 	for(int di = -1; di <= 1; ++di) {
	// 	for(int dj = -1; dj <= 1; ++dj) {
	// 	for(int dk = -1; dk <= 1; ++dk) {						
	// 		float wt = 2.0;
	// 		if( di!=0 || dj!=0 || dk!=0 ) {
	// 			wt = 1.0 / sqrt( di*di +  )
	// 		}
	// 	}
	// 	}
	// 	}
	// 	smval /= wtot;
	// 	gridsmooth(i,j,k) = smval;
	// }
	// }
	// }	

	for(int k = 1; k <= 2*N+1; k++) {
	for(int j = 1; j <= 2*N+1; j++) {
	for(int i = 1; i <= 2*N+1; i++) {
		cerr << grid(i,j,k) << endl;
	}
	}
	}	

	int nlocalmax = 0, nsamp=0;
	float mxcbc = -9e9;
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif	
	for(int i = 3; i < 2*N; i++) {
	for(int j = 3; j < 2*N; j++) {
	for(int k = 3; k < 2*N; k++) {
		nsamp++;
		float cbc = grid(i,j,k);
		bool is_local_max = true;
		if(cbc < 2.0) continue;
		mxcbc = max(mxcbc,cbc);
		float nbavg=0.0,nbmin=9e9,nbmax=0.0;
		for(int di = -1; di <= 1; ++di) {
		for(int dj = -1; dj <= 1; ++dj) {
		for(int dk = -1; dk <= 1; ++dk) {						
			if(di==0&&dj==0&&dk==0) continue;
			float nbcbc = grid(i+di,j+dj,k+dk);
			nbavg += nbcbc;
			nbmin = min(nbmin,nbcbc);
			nbmax = max(nbmax,nbcbc);			
			if( nbcbc > cbc ) is_local_max = false;
		}
		}
		}
		if(!is_local_max) continue;

		nbavg /= 26.0;
		float avg_slope = (cbc-nbavg)/step;
		float min_slope = (cbc-nbmax)/step;
		float max_slope = (cbc-nbmin)/step;		
		if( min_slope < 0.1 ) continue;
		if( avg_slope < 1.0 ) continue;
		
		
		Vec axs = axs0;
		axs = rotation_matrix_degrees(X,Real(i-N-1)*step) * axs;
		axs = rotation_matrix_degrees(Y,Real(j-N-1)*step) * axs;
		Real ang = ang0+Real(k-N-1)*step;

#ifdef USE_OPENMP
#pragma omp critical
#endif						
		TR<<"local max "<<F(7,3,cbc)<<" "<<F(6,3,min_slope)<<" "<<F(6,3,avg_slope)<<" "<<F(6,3,max_slope)<<"     "<<i<<" "<<j<<" "<<k<<" "<<axs<<" "<<F(12,8,ang)<<endl;
		nlocalmax++;
	}
	}
	}
	TR << "max " << mxcbc << " local max " << nlocalmax << " of " << nsamp << endl;
	
}


int main(int argc, char *argv[]) {
	try {

	register_options();
	devel::init(argc,argv);
	using namespace basic::options::OptionKeys;

	Size const NSS = basic::options::option[basic::options::OptionKeys::cxdock::sphere]();
	SphereSampler ssamp(NSS);
	
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
		//visualize(pnat);
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	
	return 0;
}

