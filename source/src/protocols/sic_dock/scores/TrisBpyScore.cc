// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:

#include <protocols/sic_dock/scores/TrisBpyScore.hh>
	#include <protocols/sic_dock/xyzStripeHashPose.hh>
	#include <protocols/sic_dock/util.hh>
	#include <basic/options/keys/sicdock.OptionKeys.gen.hh>
	#include <basic/options/option.hh>
	#include <basic/options/option_macros.hh>
	#include <numeric/constants.hh>
	#include <numeric/xyz.functions.hh>
	#include <numeric/xyz.io.hh>
	#include <ObjexxFCL/format.hh>
	#include <ObjexxFCL/string.functions.hh>
	#include <utility/string_util.hh>
	#include <basic/Tracer.hh>
	#include <boost/foreach.hpp>
	#include <utility/io/ozstream.hh>

namespace protocols {
namespace sic_dock {
namespace scores {


using core::Size;
	using core::Real;
	using numeric::min;
	using core::id::AtomID;
	using std::cout;
	using std::endl;
	using utility::vector1;
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using ObjexxFCL::format::F;
	using ObjexxFCL::format::I;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef Xform Xform;
	typedef numeric::xyzVector<Real> Vec;
	typedef numeric::xyzMatrix<Real> Mat;
	typedef utility::vector1<Vec> Vecs;
	typedef utility::vector1<Real> Reals;
	typedef utility::vector1<Size> Sizes;
	typedef numeric::Xforms Xforms;



	// byp tilt 32.5

static thread_local basic::Tracer TR( "protocols.sic_dock.scores.TrisBpyScore" );

static Vec   const N  ( 2.649,-1.227, 1.777);
	static Vec   const CA ( 3.946,-0.851, 1.187);
	static Vec   const C  ( 4.775, 0.037, 2.144);
	static Vec   const O  ( 4.410, 0.106, 3.359);
	static Vec   const CB ( 3.569,-0.000, 0.000);
	static Vec   const CG ( 3.018,-0.767,-1.287);
	static Vec   const CD ( 1.984,-0.166,-1.966);
	static Vec   const C28( 3.529,-1.952,-1.749);
	static Vec   const N5 ( 0.291,-1.619,-5.187);
	static Vec   const N6 ( 1.455,-0.715,-3.073);
	static Vec   const C21(-0.325,-1.956,-6.336);
	static Vec   const C22( 0.087,-3.028,-7.100);
	static Vec   const C23( 1.133,-3.793,-6.665);
	static Vec   const C24( 1.768,-3.468,-5.482);
	static Vec   const C25( 1.336,-2.377,-4.760);
	static Vec   const C29( 3.015,-2.510,-2.888);
	static Vec   const C30( 1.975,-1.875,-3.536);
	static Vec   const FE ( 0.000, 0.000,-4.157);
	static Vec   const clash_coords[12] = {CB,CD,C28,N5,N6,C21,C22,C23,C24,C25,C29,C30 };
	static Xform const bpybb(CA,N,CA,C);
	static Vec   const CBlocal(~bpybb*CB);
	static Mat   const R3z1(numeric::z_rotation_matrix_degrees(  0.0));
	static Mat   const R3z2(numeric::z_rotation_matrix_degrees(120.0));
	static Mat   const R3z3(numeric::z_rotation_matrix_degrees(240.0));
	static Mat   const R3z[3] = { R3z1, R3z2, R3z3 };
	static Mat   const chiral1(Mat::cols(1,0,0,0, 1,0,0,0, 1));
	static Mat   const chiral2(Mat::cols(1,0,0,0, 1,0,0,0,-1));
	static Mat   const chiral3(Mat::cols(1,0,0,0,-1,0,0,0, 1));
	static Mat   const chiral4(Mat::cols(1,0,0,0,-1,0,0,0,-1));
	static Mat   const chiral[4] = { chiral1, chiral2, chiral3, chiral4 };
	static Vec   const nbr2(C29);
	static Vec   const nbr1(C22);
	// static Vec   const nbr1(2.496000,-1.330833,-2.416500);
	// static Vec   const nbr2(0.715000,-2.706833,-5.921667);
	static Vec   const nbr[6] = { nbr1, nbr2, R3z2*nbr1, R3z2*nbr2, R3z3*nbr1, R3z3*nbr2 };

void printbpy(std::ostream & out, Xform const & x,int resi=1000){
	int ano=0;
	for(int i = 0; i < 3; ++i){
		out << "ATOM  "<<I(5,++ano)<<"  CB  BPY A"<<I(4,resi+i)<<"    "<<F(8,3,(R3z[i]*(x*CB )).x())<<F(8,3,(R3z[i]*(x*CB )).y())<<F(8,3,(R3z[i]*(x*CB )).z())<<endl;
		out << "ATOM  "<<I(5,++ano)<<"  CG  BPY A"<<I(4,resi+i)<<"    "<<F(8,3,(R3z[i]*(x*CG )).x())<<F(8,3,(R3z[i]*(x*CG )).y())<<F(8,3,(R3z[i]*(x*CG )).z())<<endl;
		out << "ATOM  "<<I(5,++ano)<<"  CD  BPY A"<<I(4,resi+i)<<"    "<<F(8,3,(R3z[i]*(x*CD )).x())<<F(8,3,(R3z[i]*(x*CD )).y())<<F(8,3,(R3z[i]*(x*CD )).z())<<endl;
		out << "ATOM  "<<I(5,++ano)<<"  N5  BPY A"<<I(4,resi+i)<<"    "<<F(8,3,(R3z[i]*(x*N5 )).x())<<F(8,3,(R3z[i]*(x*N5 )).y())<<F(8,3,(R3z[i]*(x*N5 )).z())<<endl;
		out << "ATOM  "<<I(5,++ano)<<"  N6  BPY A"<<I(4,resi+i)<<"    "<<F(8,3,(R3z[i]*(x*N6 )).x())<<F(8,3,(R3z[i]*(x*N6 )).y())<<F(8,3,(R3z[i]*(x*N6 )).z())<<endl;
		out << "ATOM  "<<I(5,++ano)<<"  C21 BPY A"<<I(4,resi+i)<<"    "<<F(8,3,(R3z[i]*(x*C21)).x())<<F(8,3,(R3z[i]*(x*C21)).y())<<F(8,3,(R3z[i]*(x*C21)).z())<<endl;
		out << "ATOM  "<<I(5,++ano)<<"  C22 BPY A"<<I(4,resi+i)<<"    "<<F(8,3,(R3z[i]*(x*C22)).x())<<F(8,3,(R3z[i]*(x*C22)).y())<<F(8,3,(R3z[i]*(x*C22)).z())<<endl;
		out << "ATOM  "<<I(5,++ano)<<"  C23 BPY A"<<I(4,resi+i)<<"    "<<F(8,3,(R3z[i]*(x*C23)).x())<<F(8,3,(R3z[i]*(x*C23)).y())<<F(8,3,(R3z[i]*(x*C23)).z())<<endl;
		out << "ATOM  "<<I(5,++ano)<<"  C24 BPY A"<<I(4,resi+i)<<"    "<<F(8,3,(R3z[i]*(x*C24)).x())<<F(8,3,(R3z[i]*(x*C24)).y())<<F(8,3,(R3z[i]*(x*C24)).z())<<endl;
		out << "ATOM  "<<I(5,++ano)<<"  C25 BPY A"<<I(4,resi+i)<<"    "<<F(8,3,(R3z[i]*(x*C25)).x())<<F(8,3,(R3z[i]*(x*C25)).y())<<F(8,3,(R3z[i]*(x*C25)).z())<<endl;
		out << "ATOM  "<<I(5,++ano)<<"  C28 BPY A"<<I(4,resi+i)<<"    "<<F(8,3,(R3z[i]*(x*C28)).x())<<F(8,3,(R3z[i]*(x*C28)).y())<<F(8,3,(R3z[i]*(x*C28)).z())<<endl;
		out << "ATOM  "<<I(5,++ano)<<"  C29 BPY A"<<I(4,resi+i)<<"    "<<F(8,3,(R3z[i]*(x*C29)).x())<<F(8,3,(R3z[i]*(x*C29)).y())<<F(8,3,(R3z[i]*(x*C29)).z())<<endl;
		out << "ATOM  "<<I(5,++ano)<<"  C30 BPY A"<<I(4,resi+i)<<"    "<<F(8,3,(R3z[i]*(x*C30)).x())<<F(8,3,(R3z[i]*(x*C30)).y())<<F(8,3,(R3z[i]*(x*C30)).z())<<endl;
		out << "ATOM  "<<I(5,++ano)<<"  FE  BPY A"<<I(4,resi+i)<<"    "<<F(8,3,(R3z[i]*(x*FE )).x())<<F(8,3,(R3z[i]*(x*FE )).y())<<F(8,3,(R3z[i]*(x*FE )).z())<<endl;
	}
}

using core::Size;
	using core::Real;
	using numeric::min;
	using core::id::AtomID;
	using std::cout;
	using std::endl;
	using utility::vector1;
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef Xform Xform;
	typedef numeric::xyzVector<Real> Vec;
	typedef numeric::xyzMatrix<Real> Mat;
	typedef utility::vector1<Vec> Vecs;
	typedef utility::vector1<Real> Reals;
	typedef utility::vector1<Size> Sizes;
	typedef numeric::Xforms Xforms;
	typedef utility::vector1<RigidScoreCOP> Scores;


TrisBpyScore::TrisBpyScore(
	Pose const & _pose,
	Real const & tolerance,
	Real const & min_bpy_contacts
):
	pose_(_pose),
	tolerance_(tolerance),
	min_bpy_contacts_(min_bpy_contacts)
{
	for(Size ir = 1; ir <= pose_.total_residue(); ++ir){
		Vec  N = pose_.residue(ir).xyz(1);
		Vec CA = pose_.residue(ir).xyz(2);
		Vec  C = pose_.residue(ir).xyz(3);
		bbx_.push_back(Xform(CA,N,CA,C));
		if(pose_.secstruct(ir)!='L') cb_.push_back(pose_.residue(ir).nbr_atom_xyz());
		else                         cb_.push_back(Vec(9e9,9e9,9e9));
	}
	cc_ = new xyzStripeHashPose(pose_,BB,3.2);

}

TrisBpyScore::~TrisBpyScore(){
}

void TrisBpyScore::show(std::ostream & out, int /*width*/) const {
	out << " bpy_cbc bpy_err crl";
}
void TrisBpyScore::show(std::ostream & out, Xforms const & x1s, Xforms const & x2s, int /*width*/) const {
	Real cbc=0,err=0;
	int crl=0;
	Xform tmp;
	score_extra(x1s,x2s,cbc,err,crl,tmp);
	out <<" "<< F(7,3,cbc) << " " << F(7,3,err) << " " << I(3,crl);
}

core::Real
TrisBpyScore::score( Xforms const & x1s, Xforms const & x2s ) const {
	Real cbc=0,err=0;
	int crl=0;
	Xform tmp;
	return score_extra(x1s,x2s,cbc,err,crl,tmp);
}

core::Real
TrisBpyScore::score_extra( Xforms const & x1s, Xforms const & x2s, Real&cbc,Real&err,int&crl, Xform&xbpy ) const {
	if(x1s.size()!=1) utility_exit_with_message("foo");
	if(x2s.size()!=1) utility_exit_with_message("foo");
	BOOST_FOREACH(Xform const & x1,x1s){
		for(core::Size ir = 1; ir <= bbx_.size(); ++ir){
			Xform const & bbx(bbx_[ir]);
			Vec cb = x1*bbx*CBlocal;
			Vec tmp(cb.x(),cb.y(),0);
			Real dxy = tmp.length();
			err = fabs(dxy-3.569);
			if( err > tolerance_ ) continue;
			Vec ca = x1*bbx.t;
			Real ang = numeric::dihedral(cb,Vec(0,0,0),Vec(0,0,1), Vec(1,0,0));
			Xform xbpy0( numeric::z_rotation_matrix_degrees(-ang), Vec(0,0,cb.z()) );
			for(int ichiral=0; ichiral<4; ++ichiral){
				xbpy = xbpy0*chiral[ichiral];
				Xform const xclash = ~x1*xbpy;
				Vec const cg = xbpy*CG;
				Real angerr = fabs(numeric::angle_degrees(ca,cb,cg)-116.7);
				if( angerr > tolerance_*8.8 ) continue;
				err += angerr/8.8;
				bool clash=false;
				for(int ia=0; ia < 11; ++ia){
					if(cc_->clash_not_resid(xclash*R3z[0]*clash_coords[ia],ir)){ clash=true; break; }
					if(cc_->clash          (xclash*R3z[1]*clash_coords[ia]   )){ clash=true; break; }
					if(cc_->clash          (xclash*R3z[2]*clash_coords[ia]   )){ clash=true; break; }
				}
				if(clash) continue;
				cbc = 0;
				Vec const nbrcb[6] = { xclash*nbr[0],xclash*nbr[1],xclash*nbr[2],xclash*nbr[3],xclash*nbr[4],xclash*nbr[5] };
				// for(int jr=   1; jr <= ir-3       ; ++jr) for(int inbr=0; inbr<6; ++inbr) cbc += CBScore_dist_score(cb_[jr].distance_squared(nbrcb[inbr]),6.0,12.0);
				// for(int jr=ir+3; jr <= bbx_.size(); ++jr) for(int inbr=0; inbr<6; ++inbr) cbc += CBScore_dist_score(cb_[jr].distance_squared(nbrcb[inbr]),6.0,12.0);
				for(int jr=4; jr <= (int)bbx_.size()-3; ++jr) for(int inbr=0; inbr<6; ++inbr) cbc += CBScore_dist_score(cb_[jr].distance_squared(nbrcb[inbr]),7.0,15.0);
				if(cbc < min_bpy_contacts_) continue;
				// if(numeric::random::uniform()<01.01){
				// 	#ifdef USE_OPENMP
				// 		#pragma omp critical
				// 	#endif
				// 	{
				// 		utility::io::ozstream out("bpytest.pdb");
				// 		printbpy(out,xbpy);
				// 		out.close();
				// 		Pose tmp(pose_);
				// 		xform_pose(tmp,x1);
				// 		tmp.dump_pdb("test1.pdb");
				// 		utility_exit_with_message("foo");
				// 	}
				// }
				crl=ichiral;
				return 10000.0 + 3.0*cbc;
			}
		}
	}
	return 0.0;
}

void TrisBpyScore::dump_trisbpy( Xforms const & x1s, Xforms const & x2s, std::string const & fname ) const{
	Real cbc=0,err=0;
	int crl=0;
	Xform xbpy;
	score_extra(x1s,x2s,cbc,err,crl,xbpy);
	utility::io::ozstream out(fname);
	printbpy(out,xbpy);
	out.close();
}



} // namespace scores
} // namespace sic_dock
} // namespace protocols
