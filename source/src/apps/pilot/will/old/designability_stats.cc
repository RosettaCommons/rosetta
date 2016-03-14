#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
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

#include <protocols/sic_dock/designability_score.hh>

typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;

using core::id::AtomID;
using basic::options::option;
using namespace basic::options::OptionKeys;
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
using core::kinematics::Stub;

#define MAX_CYS_RES 0
#define MAX_NRES 1000

OPT_1GRP_KEY( Integer , dstat, min_seq_sep )
OPT_1GRP_KEY( File    , dstat, make_hist   )
OPT_1GRP_KEY( File   , dstat, test_score  )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( dstat::min_seq_sep    ,"",  20 );
	NEW_OPT( dstat::make_hist      ,"",  "" );
	NEW_OPT( dstat::test_score     ,"",  "" );
}

// void
// get_xform_stats(
// 	core::kinematics::Stub const & sir,
// 	core::kinematics::Stub const & sjr,
// 	Real& dx, Real& dy, Real& dz,
// 	Real& ex, Real& ey, Real& ez
// ){
// 	Vec d = sir.global2local(sjr.v);
// 	dx = d.x();
// 	dy = d.y();
// 	dz = d.z();
// 	Mat R = sir.M.transposed()*sjr.M;
// 	// Real ang;
// 	// Vec axis = rotation_axis(R,ang);
// 	// ang = numeric::conversions::degrees(ang);
// 	// std::cout << axis << " " << ang << std::endl;

// 	Real phi,psi,theta;
// 	if( fabs(R.zx())-1.0 < 0.00001 ) {
// 		theta = -asin(R.zx());
// 		Real const ctheta = cos(theta);
// 		psi = atan2(R.zy()/ctheta,R.zz()/ctheta);
// 		phi = atan2(R.yx()/ctheta,R.xx()/ctheta);
// 	} else {
// 		if( R.zx() < 0.0 ) {
// 			theta = numeric::constants::d::pi / 2.0;
// 			psi = atan2(R.xy(),R.xz());
// 		} else {
// 			theta = -numeric::constants::d::pi / 2.0;
// 			psi = atan2(-R.xy(),-R.xz());
// 		}
// 		phi = 0;
// 	}
// 	ex = psi;
// 	ey = theta;
// 	ez = phi;
// }

// struct
// XfoxmScore
// {
// 	char *hh,*he,*hl,*ee,*el,*ll;
// 	XfoxmScore(
// 		std::string datadir
// 	){
// 		hh = new char[16*16*16*24*12*24];
// 		he = new char[16*16*16*24*12*24];
// 		hl = new char[16*16*16*24*12*24];
// 		ee = new char[16*16*16*24*12*24];
// 		el = new char[16*16*16*24*12*24];
// 		ll = new char[16*16*16*24*12*24];
// 		fillarray(hh,datadir+"/hhpb.dat.gz.hist6.dat.gz.bin.gz");
// 		fillarray(he,datadir+"/hepb.dat.gz.hist6.dat.gz.bin.gz");
// 		fillarray(hl,datadir+"/hlpb.dat.gz.hist6.dat.gz.bin.gz");
// 		fillarray(ee,datadir+"/eepb.dat.gz.hist6.dat.gz.bin.gz");
// 		fillarray(el,datadir+"/elpb.dat.gz.hist6.dat.gz.bin.gz");
// 		fillarray(ll,datadir+"/llpb.dat.gz.hist6.dat.gz.bin.gz");
// 		// makebinary(hh,datadir+"/hhpb.dat.gz.hist6.dat.gz"); // uncomment to gen binary files
// 		// makebinary(he,datadir+"/hepb.dat.gz.hist6.dat.gz"); // uncomment to gen binary files
// 		// makebinary(hl,datadir+"/hlpb.dat.gz.hist6.dat.gz"); // uncomment to gen binary files
// 		// makebinary(ee,datadir+"/eepb.dat.gz.hist6.dat.gz"); // uncomment to gen binary files
// 		// makebinary(el,datadir+"/elpb.dat.gz.hist6.dat.gz"); // uncomment to gen binary files
// 		// makebinary(ll,datadir+"/llpb.dat.gz.hist6.dat.gz"); // uncomment to gen binary files
// 		// utility_exit_with_message("made binary files");
// 	}
// 	void
// 	fillarray(
// 		char *a, std::string fname
// 	){
// 		std::cout << "reading " << fname << std::endl;
// 		utility::io::izstream in(fname,std::ios::binary);
// 		if(!in.good()) utility_exit_with_message("bad file");
// 		in.read(a,16*16*16*24*12*24);
// 		in.close();
// 	}
// 	void
// 	makebinary(
// 		char *a, std::string fname
// 	){
// 		std::cout << "reading " << fname << std::endl;
// 		int Nstat = 16*16*16*24*12*24;
// 		{
// 			utility::io::izstream in(fname);
// 			if(!in.good()) utility_exit_with_message("bad file");
// 			for(int i = 0; i < Nstat; ++i){
// 				float f;
// 				in >> f;
// 				int tmp = (int)((log(f)+2.5)*12.0);
// 				a[i] =  (f==0.0f) ? ((char)-127) : ((char)max(-126,min(128,tmp)));
// 				// std::cout << f << " " << (int)a[i] << std::endl;
// 				// if( i > 100) utility_exit_with_message("FOO");
// 			}
// 			in.close();
// 		}
// 		std::cout << "writing: " << fname+".bin.gz" << std::endl;
// 		utility::io::ozstream out(fname+".bin.gz",std::ios::out | std::ios::binary);
// 		out.write(a,Nstat);
// 		out.close();
// 		{ // verify file
// 			std::cout << "verifying: " << fname+".bin.gz" << std::endl;
// 			utility::io::izstream in(fname+".bin.gz",std::ios::binary);
// 			char *tmp = new char[Nstat];
// 			in.read(tmp,Nstat);
// 			for(int i = 0; i < Nstat; ++i) {
// 				if( tmp[i] != a[i] ) utility_exit_with_message("binary file is wrong!");
// 			}
// 			in.close();
// 		}
// 		// std::cout << "writing: " << fname+".txt.gz for testing" << std::endl;
// 		// utility::io::ozstream out2(fname+".txt.gz",std::ios::out);
// 		// for(int i = 0; i < 16*16*16*24*12*24; ++i) out2 << (int)a[i] << std::endl;
// 		// out2.close();
// 	}
// 	float
// 	score(
// 		core::kinematics::Stub const & s1,
// 		core::kinematics::Stub const & s2,
// 		char ss1, char ss2
// 	) const {
// 		using numeric::constants::d::pi_2;
// 		if( s1.global2local(s2.v).length_squared() > 64.0 ) return 0.0;
// 		char *a = hh;
// 		if( ss1=='L' && ss2=='L' ) a = ee;
// 		if( ss1=='L' && ss2=='L' ) a = ll;
// 		if( ss1=='H' && ss2=='E' || ss1=='E' && ss2=='H' ) a = he;
// 		if( ss1=='H' && ss2=='L' || ss1=='L' && ss2=='H' ) a = hl;
// 		if( ss1=='E' && ss2=='L' || ss1=='L' && ss2=='E' ) a = el;
// 		Real dx,dy,dz,ex,ey,ez;
// 		if(ss1=='E' && ss2=='H' || ss1=='L' && ss2=='H' || ss1=='L' && ss2=='E'){
// 			   get_xform_stats(s2,s1,dx,dy,dz,ex,ey,ez); // reverse
// 		} else get_xform_stats(s1,s2,dx,dy,dz,ex,ey,ez);
// 		int idx = dx+8.0; // 0.0-0.999 -> 10
// 		int idy = dy+8.0;
// 		int idz = dz+8.0;
// 		int iex = ex/pi_2*24.0 + 12.0;
// 		int iey = ey/pi_2*24.0 +  6.0;
// 		int iez = ez/pi_2*24.0 + 12.0;
// 		int index = idx + 16*idy + 16*16*idz + 16*16*16*iex + 16*16*16*24*iey + 16*16*16*24*12*iez;
// 		if( 0 > index || index >= 16*16*16*24*12*24 ) utility_exit_with_message("FOO");
// 		// expensive memory lookup
// 		char val = a[index];
// 		//
// 		// std::cout << (int)val << "'" << val << "'"<< std::endl;
// 		return val == -127 ? 0.0 : exp(((float)val)/12.0-2.5);
// 	}
// 	float
// 	score(
// 		core::pose::Pose const & pose,
// 		Size rsd1,
// 		Size rsd2
// 	) const {
// 		if(!pose.residue(rsd1).is_protein()) return -1.0;
// 		if(!pose.residue(rsd2).is_protein()) return -1.0;
// 		if(!pose.residue(rsd1).has("CB")) return -1.0;
// 		if(!pose.residue(rsd2).has("CB")) return -1.0;
// 		Vec CBi = pose.residue(rsd1).xyz("CB");
// 		Vec CAi = pose.residue(rsd1).xyz("CA");
// 		Vec  Ni = pose.residue(rsd1).xyz( "N");
// 		core::kinematics::Stub sir(CBi,CAi,Ni);
// 		Vec CBj = pose.residue(rsd2).xyz("CB");
// 		Vec CAj = pose.residue(rsd2).xyz("CA");
// 		Vec  Nj = pose.residue(rsd2).xyz( "N");
// 		if( CBi.distance_squared(CBj) > 64.0 ) return -1.0;
// 		core::kinematics::Stub sjr(CBj,CAj,Nj);
// 		return score(sir,sjr,pose.secstruct(rsd1),pose.secstruct(rsd2));
// 	}
// 	float
// 	score(
// 		core::pose::Pose & pose,
// 		bool compute_ss = true
// 	) const {
// 		if(compute_ss){
// 			core::scoring::dssp::Dssp dssp(pose);
// 			dssp.insert_ss_into_pose(pose);
// 		}
// 		float tot_score = 0.0;
// 		for(Size ir = 1; ir <= pose.n_residue(); ++ir){
// 			for(Size jr = ir+1; jr <= pose.n_residue(); ++jr){
// 				float s1 = score(pose,ir,jr);
// 				float s2 = score(pose,jr,ir);
// 				// if( s1 < 0.0f ) std::cout << s1 << std::endl;
// 				// if( s2 < 0.0f ) std::cout << s2 << std::endl;
// 				tot_score += (s1<0.0f) ? 0.0f : s1;
// 				tot_score += (s2<0.0f) ? 0.0f : s2;
// 			}
// 		}
// 		return tot_score;
// 	}
// 	float
// 	score(
// 		core::pose::Pose const & pose
// 	) const {
// 		float tot_score = 0.0;
// 		for(Size ir = 1; ir <= pose.n_residue(); ++ir){
// 			for(Size jr = ir+1; jr <= pose.n_residue(); ++jr){
// 				float s1 = score(pose,ir,jr);
// 				float s2 = score(pose,jr,ir);
// 				// if( s1 < 0.0f ) std::cout << s1 << std::endl;
// 				// if( s2 < 0.0f ) std::cout << s2 << std::endl;
// 				tot_score += (s1<0.0f) ? 0.0f : s1;
// 				tot_score += (s2<0.0f) ? 0.0f : s2;
// 			}
// 		}
// 		return tot_score;
// 	}
// };


void
collect_stats(
	core::pose::Pose & pose,
	std::string tag
){
	Pose ala;
 	core::pose::make_pose_from_sequence(ala,"A",core::chemical::FA_STANDARD,false);
	remove_lower_terminus_type_from_pose_residue(ala,1);
	remove_upper_terminus_type_from_pose_residue(ala,1);
	for(Size ir = 1; ir <= pose.n_residue(); ++ir){
		if(!pose.residue(ir).is_protein()) continue;
		if(!pose.residue(ir).has("CB")) continue;
		pose.replace_residue(ir,ala.residue(1),true);

	}
	// pose.dump_pdb("test.pdb");
	for(Size ir = 1; ir <= pose.n_residue(); ++ir){
		if(!pose.residue(ir).is_protein()) continue;
		if(!pose.residue(ir).has("CB")) continue;
		Vec CBi = pose.residue(ir).xyz("CB");
		Vec CAi = pose.residue(ir).xyz("CA");
		Vec  Ni = pose.residue(ir).xyz( "N");
		core::kinematics::Stub sir(CBi,CAi,Ni);
		for(Size jr = ir + option[dstat::min_seq_sep](); jr <= pose.n_residue(); ++jr){
			if(!pose.residue(jr).is_protein()) continue;
			if(!pose.residue(jr).has("CB")) continue;
			Vec CBj = pose.residue(jr).xyz("CB");
			Vec CAj = pose.residue(jr).xyz("CA");
			Vec  Nj = pose.residue(jr).xyz( "N");
			if( CBi.distance_squared(CBj) > 100.0 ) continue;
			core::kinematics::Stub sjr(CBj,CAj,Nj);
			Real dx,dy,dz,ex,ey,ez;
			protocols::sic_dock::get_xform_stats(sir,sjr,dx,dy,dz,ex,ey,ez);
			std::cout<<"RES_RES_XFORM"<<tag<<" "<<ir<<" "<<jr<<" "<<pose.secstruct(ir)<<" "<<pose.secstruct(jr)<<" "<<dx<<" "<<dy<<" "<<dz<<" "<<ex<<" "<<ey<<" "<<ez<<std::endl;
			// utility_exit_with_message("foo");
		}
	}
}

int main(int argc, char *argv[]) {

	try {

	register_options();
	devel::init(argc,argv);
	using namespace basic::options::OptionKeys;
	using numeric::constants::d::pi_2;

	if( option[dstat::test_score].user() ){
		protocols::sic_dock::XfoxmScore xfs(option[dstat::test_score]());
		for(Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn) {
			string fn = option[in::file::s]()[ifn];
			Pose pose;
			core::import_pose::pose_from_file(pose,fn, core::import_pose::PDB_file);
			// std::cout << fn << " " << xfs.score(pose,true) << std::endl;
			Size count = 0;
			Real tot = 0.0;
			std::cout << "computing scoeres" << std::endl;
			for(Real rotx = 0; rotx < 360.0; rotx += 15.0 ){
				std::cout << rotx << " " << count << std::endl;
				for(Real roty = 0; roty < 360.0; roty += 15.0 ){
					for(Real rotz = 0; rotz < 360.0; rotz += 15.0 ){
						using namespace numeric;
						Mat R = x_rotation_matrix_degrees(rotx)*y_rotation_matrix_degrees(roty)*z_rotation_matrix_degrees(rotz);
				 		for(Size ir = 1; ir <= pose.n_residue(); ++ir){
							if(!pose.residue(ir).is_protein()) continue;
							if(pose.residue(ir).aa()==core::chemical::aa_gly) continue;
		 					for(Size jr = ir+1; jr <= pose.n_residue(); ++jr){
								if(!pose.residue(jr).is_protein()) continue;
								if(pose.residue(jr).aa()==core::chemical::aa_gly) continue;
								Vec CBi = pose.residue(ir).xyz(5);
								Vec CAi = pose.residue(ir).xyz(2);
								Vec  Ni = pose.residue(ir).xyz(1);
								core::kinematics::Stub sir(CBi,CAi,Ni);
								Vec CBj = pose.residue(jr).xyz(5);
								Vec CAj = pose.residue(jr).xyz(2);
								Vec  Nj = pose.residue(jr).xyz(1);
								if( CBi.distance_squared(CBj) > 64.0 ) continue;
								core::kinematics::Stub sjr(CBj,CAj,Nj);
								sjr.M = R*sjr.M;
								count++;
								tot += xfs.score(sir,sjr,pose.secstruct(ir),pose.secstruct(jr));
							}
						}
					}
				}
			}
			std::cout << "DONE " << tot << " " << count << std::endl;
		}
	} else if( option[dstat::make_hist].user() ){
		utility::io::izstream in(option[dstat::make_hist]());
		double dx,dy,dz,ex,ey,ez;
		float* hist = new float[16*16*16*24*12*24];
		int count = 0;
		while(in >> dx >> dy >> dz >> ex >> ey >> ez) {

			if( count % 100000 == 0 ) std::cout << "filling hist " << count << std::endl;

			int idx = dx+8.0; // 0.0-0.999 -> 10
			int idy = dy+8.0;
			int idz = dz+8.0;
			int iex = ex/pi_2*24.0 + 12.0;
			int iey = ey/pi_2*24.0 +  6.0;
			int iez = ez/pi_2*24.0 + 12.0;
			// if( iey <  0) utility_exit_with_message("FOO");
			// if( iey > 11) utility_exit_with_message("FOO");
			// if( idx <  0 || idy <  0 || idz <  0 ) continue;
			// if( idx > 15 || idy > 15 || idz > 15 ) continue;
			// if( iex <  0 || iez <  0 ) continue;
			// if( iex > 23 || iez > 23 ) continue;

			if(0) {
				int index = idx + 16*idy + 16*16*idz + 16*16*16*iex + 16*16*16*24*iey + 16*16*16*24*12*iez;
				if( 0 > index || index >= 16*16*16*24*12*24 ) utility_exit_with_message("FOO");
				hist[index] += 1.0;
			} else {
				for( int didx = -1; didx <= 1; ++didx ) {
				for( int didy = -1; didy <= 1; ++didy ) {
				for( int didz = -1; didz <= 1; ++didz ) {
				for( int diex = -1; diex <= 1; ++diex ) {
				for( int diey = -1; diey <= 1; ++diey ) {
				for( int diez = -1; diez <= 1; ++diez ) {

					int i = (idx+didx);
					int j = (idy+didy);
					int k = (idz+didz);
					int l = (iex+diex);
					int m = (iey+diey);
					int n = (iez+diez);

					if( 0 > i || i > 15 ) continue;
					if( 0 > j || j > 15 ) continue;
					if( 0 > k || k > 15 ) continue;
					if( 0 > l || l > 23 ) continue;
					if( 0 > m || m > 11 ) continue;
					if( 0 > n || n > 23 ) continue;

					int index = i + 16*j + 16*16*k + 16*16*16*l + 16*16*16*24*m + 16*16*16*24*12*n;
					if( 0 > index || index >= 16*16*16*24*12*24 ) utility_exit_with_message("bad index");

					double cdx = (double)i - 7.5;
					double cdy = (double)j - 7.5;
					double cdz = (double)k - 7.5;
					double cex = ((double)l - 11.5); // / 24.0 * pi_2;
					double cey = ((double)m -  5.5); // / 24.0 * pi_2;
					double cez = ((double)n - 11.5); // / 24.0 * pi_2;

					double ndelt_i = (cdx-dx);
					double ndelt_j = (cdy-dy);
					double ndelt_k = (cdz-dz);
					double ndelt_l = (cex-ex/pi_2*24.0);
					double ndelt_m = (cey-ey/pi_2*24.0);
					double ndelt_n = (cez-ez/pi_2*24.0);

					// std::cout <<ndelt_i<<" "<<ndelt_j<<" "<<ndelt_k<<" "<<ndelt_l<<" "<<ndelt_m<<" "<<ndelt_n<<std::endl;

					double dist2 = ndelt_i*ndelt_i+ndelt_j*ndelt_j+ndelt_k*ndelt_k+ndelt_l*ndelt_l+ndelt_m*ndelt_m+ndelt_n*ndelt_n;

					// std::cout << sqrt(dist2) << " " << exp(-dist2) << std::endl;

					// if(count > 10) utility_exit_with_message("arst");

					hist[index] += exp(-dist2);

				}}}}}}
			}

			// std::cout << dx << " " << dy << " " << dz << " " << ex << " " << ey << " " << ez << std::endl;
			// utility_exit_with_message("FOO");
			count++;
		}
		in.close();
		std::cout << "read " << count << " points" << std::endl;
		if( count < 100 ) utility_exit_with_message("no data!");
		std::string outfn = utility::file_basename(option[dstat::make_hist]()) + ".hist6.dat.gz";
		std::cout << "writting " << outfn << std::endl;
		std::ostringstream ssout;
		for(int i = 0; i < 16*16*16*24*12*24; ++i){
			if( i % 1000000 == 0 ) std::cout << "writting... " << i << " of " << 16*16*16*24*12*24 << std::endl;
			ssout << hist[i] << std::endl;
		}
		utility::io::ozstream out(outfn);
		out << ssout.str();
		out.close();
	} else {
		// loop over input files, do some checks, call dock
		for(Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn) {
			string fn = option[in::file::s]()[ifn];
			Pose p;
			core::import_pose::pose_from_file(p,fn, core::import_pose::PDB_file);
			core::scoring::dssp::Dssp dssp(p);
			dssp.insert_ss_into_pose(p);
			if( p.n_residue() > MAX_NRES ) continue;
			Size cyscnt=0;
			for(Size ir = 2; ir <= p.n_residue()-1; ++ir) {
				if(p.residue(ir).name3()=="CYS") { if(++cyscnt > MAX_CYS_RES) goto cont1; }
			} goto done1; cont1: std::cerr << "skipping " << fn << std::endl; continue; done1:
			std::cout << "searching " << fn << std::endl;
			collect_stats(p,utility::file_basename(fn));
		}
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
