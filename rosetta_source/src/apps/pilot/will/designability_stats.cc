#define MAX_CYS_RES 0
#define MAX_NRES 1000
#define MINSEQSEP 20

#include <basic/options/option.hh>
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

typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;

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


void
get_xform_stats(
	core::kinematics::Stub const & sir,
	core::kinematics::Stub const & sjr,
	Real& dx, Real& dy, Real& dz,
	Real& ex, Real& ey, Real& ez
){
	Vec d = sir.global2local(sjr.v);
	dx = d.x();
	dy = d.y();
	dz = d.z();
	Mat R = sir.M.transposed()*sjr.M;
	// Real ang;
	// Vec axis = rotation_axis(R,ang);
	// ang = numeric::conversions::degrees(ang);
	// std::cout << axis << " " << ang << std::endl;

	Real phi,psi,theta;
	if( fabs(R.zx())-1.0 < 0.00001 ) {
		theta = -asin(R.zx());
		Real const ctheta = cos(theta);
		psi = atan2(R.zy()/ctheta,R.zz()/ctheta);
		phi = atan2(R.yx()/ctheta,R.xx()/ctheta);
	} else {
		if( R.zx() < 0.0 ) {
			theta = numeric::constants::d::pi / 2.0;
			psi = atan2(R.xy(),R.xz());
		} else {
			theta = -numeric::constants::d::pi / 2.0;
			psi = atan2(-R.xy(),-R.xz());			
		}
		phi = 0;
	}
	ex = psi;
	ey = theta;
	ez = phi;
}

void
collect_stats(core::pose::Pose & pose, std::string tag) {
	Pose ala;
 	core::pose::make_pose_from_sequence(ala,"A","fa_standard",false);
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
		for(Size jr = ir + MINSEQSEP; jr <= pose.n_residue(); ++jr){
			if(!pose.residue(jr).is_protein()) continue;
			if(!pose.residue(jr).has("CB")) continue;
			Vec CBj = pose.residue(jr).xyz("CB");
			Vec CAj = pose.residue(jr).xyz("CA");
			Vec  Nj = pose.residue(jr).xyz( "N");
			if( CBi.distance_squared(CBj) > 100.0 ) continue;
			core::kinematics::Stub sjr(CBj,CAj,Nj);
			Real dx,dy,dz,ex,ey,ez;
			get_xform_stats(sir,sjr,dx,dy,dz,ex,ey,ez);
			std::cout<<"RES_RES_XFORM"<<tag<<" "<<ir<<" "<<jr<<" "<<pose.secstruct(ir)<<" "<<pose.secstruct(jr)<<" "<<dx<<" "<<dy<<" "<<dz<<" "<<ex<<" "<<ey<<" "<<ez<<std::endl;
			// utility_exit_with_message("foo");
		}
	}
}

int main(int argc, char *argv[]) {
	devel::init(argc,argv);
	using namespace basic::options::OptionKeys;

	// loop over input files, do some checks, call dock
	for(Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn) {
		string fn = option[in::file::s]()[ifn];
		Pose pnat;
		core::import_pose::pose_from_pdb(pnat,fn);
		core::scoring::dssp::Dssp dssp(pnat);
		dssp.insert_ss_into_pose(pnat);
		if( pnat.n_residue() > MAX_NRES ) continue;
		Size cyscnt=0;
		for(Size ir = 2; ir <= pnat.n_residue()-1; ++ir) {
			if(pnat.residue(ir).name3()=="CYS") { if(++cyscnt > MAX_CYS_RES) goto cont1; }
		} goto done1; cont1: std::cerr << "skipping " << fn << std::endl; continue; done1:
		std::cout << "searching " << fn << std::endl;
		collect_stats(pnat,utility::file_basename(fn));
	}
}
