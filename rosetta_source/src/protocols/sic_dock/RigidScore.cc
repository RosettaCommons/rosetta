// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:

#include <protocols/sic_dock/RigidScore.hh>
#include <protocols/sic_dock/xyzStripeHashPose.hh>
#include <protocols/sic_dock/util.hh>
#include <protocols/sic_dock/loophash_util.hh>

#include <basic/options/keys/sicdock.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/sasa.hh>
#include <core/pose/util.hh>
#include <protocols/loophash/BackboneDB.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace sic_dock {

static basic::Tracer TR( "propocols.sic_dock.RigidScore" );

using core::Size;
using core::Real;
using numeric::min;
using core::id::AtomID;
using std::cout;
using std::endl;
typedef core::Real Real;
typedef core::Size Size;
typedef core::pose::Pose Pose;
typedef core::kinematics::Stub Stub;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
typedef utility::vector1<Vec> Vecs;
typedef utility::vector1<Real> Reals;
typedef utility::vector1<Size> Sizes;
typedef utility::vector1<Stub> Stubs;
typedef utility::vector1<RigidScoreCOP> Scores;

inline double dist_score( double const & sqdist, double const & start, double const & stop ) {
	if( sqdist > stop*stop ) {
		return 0.0;
	} else if( sqdist < start*start ) {
		return 1.0;
	} else {
		double dist = sqrt( sqdist );
		return (stop-dist)/(stop-start);
		//return sqr(1.0	- sqr( (dist - start) / (stop - start) ) );
	}
}

std::ostream & operator<< ( std::ostream & out, Vec3 v ){
	out << v.a << " " << v.b << " " << v.c;
	return out;
}

CBScore::CBScore(
	Pose const & pose1,
	Pose const & pose2,
	Real clash_dis,
	Real contact_dis
):
   	hash_pose1_( pose1_has_most_CBs(pose1,pose2) ),
	clash_dis_(clash_dis),
    contact_dis_(contact_dis),
    weights_( cb_weights_from_pose( hash_pose1_?pose2:pose1 ) ),
    points_(           get_CB_Vecs( hash_pose1_?pose2:pose1 ) ),
    xyzhash_(
    	contact_dis,
    	hash_pose1_?pose1:pose2,
    	cb_weight_map_from_pose( hash_pose1_?pose1:pose2 )
    )
    // ,
    // pose1_(pose1),
    // pose2_(pose2)
{}

core::Real
CBScore::score(
	Stub const & x1,
	Stub const & x2
) const {
	Stub const & xh( hash_pose1_?x1:x2 );
	Stub const & xp( hash_pose1_?x2:x1 );	

	Vecs xpoints(points_);
	for(Vecs::iterator i = xpoints.begin(); i != xpoints.end(); ++i){
		*i = xp.local2global(*i);
	}

	xyzStripeHashPose const & H(xyzhash_);
	Real score = 0.0;
	Reals::const_iterator iwb = weights_.begin();
	for(Vecs::const_iterator i = xpoints.begin(); i != xpoints.end(); ++i,++iwb) {
		Vec v = xh.global2local((*i)) + H.translation();
		if( v.x() < -H.grid_size_ || v.y() < -H.grid_size_ || v.z() < -H.grid_size_ ) continue; // worth it?
		if( v.x() >  H.xmx_       || v.y() >  H.ymx_       || v.z() >  H.zmx_       ) continue; // worth it?
		int const ix  = (v.x()<0.0) ? 0 : (int)(numeric::min(H.xdim_-1,(int)(v.x()/H.grid_size_)));
		int const iy0 = (v.y()<0.0) ? 0 : (int)(v.y()/H.grid_size_);
		int const iz0 = (v.z()<0.0) ? 0 : (int)(v.z()/H.grid_size_);
		int const iyl = numeric::max(0,iy0-1);
		int const izl = numeric::max(0,iz0-1);
		int const iyu = numeric::min((int)H.ydim_,     iy0+2);
		int const izu = numeric::min((int)H.zdim_,(int)iz0+2);
		for(int iy = iyl; iy < iyu; ++iy) {
			for(int iz = izl; iz < izu; ++iz) {
				int const ig = ix+H.xdim_*iy+H.xdim_*H.ydim_*iz;
				assert(ig < H.xdim_*H.ydim_*H.zdim_ && ix < H.xdim_ && iy < H.ydim_ && iz < H.zdim_);
				int const igl = H.grid_stripe_[ig].x;
				int const igu = H.grid_stripe_[ig].y;
				for(int i = igl; i < igu; ++i) {
					numeric::geometry::hashing::xyzStripeHash<double>::float4 const & a2 = H.grid_atoms_[i];
					Real const d2 = (v.x()-a2.x)*(v.x()-a2.x) + (v.y()-a2.y)*(v.y()-a2.y) + (v.z()-a2.z)*(v.z()-a2.z);
					if( d2 <= H.grid_size2_ ) {
						score += dist_score(d2, clash_dis_, contact_dis_ ) * a2.w * (*iwb);
					}
				}
			}
		}
	}

	// // // test
	// utility::vector1<core::Real> w1,w2;
	// Pose tmp1(pose1_),tmp2(pose2_);
	// for(Size i = 1; i <= tmp1.n_residue(); ++i) w1.push_back(numeric::min(1.0,(double)neighbor_count(tmp1,i)/20.0));
	// for(Size i = 1; i <= tmp2.n_residue(); ++i) w2.push_back(numeric::min(1.0,(double)neighbor_count(tmp2,i)/20.0));
	// tmp1.dump_pdb("test0.pdb");
	// tmp2.dump_pdb("test1.pdb");	
	// xform_pose(tmp1,x1);
	// xform_pose(tmp2,x2);
	// tmp1.dump_pdb("test2.pdb");
	// tmp2.dump_pdb("test3.pdb");	
	// Real score2;
	// for(Size ir = 1; ir <= tmp1.n_residue(); ++ir){
	// 	if(!tmp1.residue(ir).has("CB")) continue;
	// 	Vec p1 = tmp1.residue(ir).xyz("CB");
	// 	for(Size jr = 1; jr <= tmp2.n_residue(); ++jr){
	// 		if(!tmp2.residue(jr).has("CB")) continue;
	// 		Vec p2 = tmp2.residue(jr).xyz("CB");
	// 		Real d2 = p2.distance_squared(p1);
	// 		score2 += dist_score(d2, clash_dis_, contact_dis_ ) * w1[ir] * w2[jr];
	// 	}
	// }
	// std::TR << score << " " << score2 << std::endl;
	// utility_exit_with_message("FOO");

	return score;
}



LinkerScore::LinkerScore(
	Pose const & pose1,
	Pose const & pose2,
	Size max_loop_len,
	Size lookup_radius
):
	loopsizes_( range(3,max_loop_len+1) ),
	lookup_radius_( lookup_radius ),
	pose1_(pose1),
	pose2_(pose2)
{
	if( ! basic::options::option[ basic::options::OptionKeys::lh::db_path ].user() ){
		utility_exit_with_message("user must specify -lh:db_path for loop hash database");
	}

	loop_hash_library_ = new protocols::loophash::LoopHashLibrary ( loopsizes_, 1, 0 );
	TR << "loading loophash data" << std::endl;
	if(loopsizes_.size()) loop_hash_library_->load_mergeddb();
	TR << "done loading loophash data" << std::endl;

	get_termini_from_pose(pose1_,uppers1_,lowers1_);
	get_termini_from_pose(pose2_,uppers2_,lowers2_);

	TR << lowers1_.size() << " " << uppers1_.size() << " " << lowers2_.size() << " " << uppers2_.size() << std::endl;
	// utility_exit_with_message("TEST LINK SCORE");
}



core::Real
LinkerScore::score(
	Stub const & x1,
	Stub const & x2
) const {
	using namespace basic::options::OptionKeys;
	using namespace protocols::loophash;
	using namespace numeric::geometry::hashing;
	using namespace core::kinematics;

	Real lkscore = 0.0;
	for(TermInfo::const_iterator c1 = lowers1_.begin(); c1 != lowers1_.end(); ++c1){		
	for(TermInfo::const_iterator n2 = uppers2_.begin(); n2 != uppers2_.end(); ++n2){
		Stub const lower = vec3_to_stub(x1,c1->second);
		Stub const upper = vec3_to_stub(x2,n2->second);
		Size n = count_linkers( lower, upper, loop_hash_library_, loopsizes_, lookup_radius_ );
		Real s = linker_count2score(n);
		lkscore += s;
	}}
	for(TermInfo::const_iterator c2 = lowers2_.begin(); c2 != lowers2_.end(); ++c2){		
	for(TermInfo::const_iterator n1 = uppers1_.begin(); n1 != uppers1_.end(); ++n1){
		Stub const lower = vec3_to_stub(x2,c2->second);
		Stub const upper = vec3_to_stub(x1,n1->second);
		Size n = count_linkers( lower, upper, loop_hash_library_, loopsizes_, lookup_radius_ );
		Real s = linker_count2score(n);
		lkscore += s;
	}}

	if(lkscore > 3.0){
	 	Size ndumped = 0;
		for(TermInfo::const_iterator c1 = lowers1_.begin(); c1 != lowers1_.end(); ++c1){		
		for(TermInfo::const_iterator n2 = uppers2_.begin(); n2 != uppers2_.end(); ++n2){
			Stub const lower = vec3_to_stub(x1,c1->second);
			Stub const upper = vec3_to_stub(x2,n2->second);
			ndumped += dump_loophash_linkers( lower, upper, loop_hash_library_, loopsizes_, lookup_radius_ );
		}}
		for(TermInfo::const_iterator c2 = lowers2_.begin(); c2 != lowers2_.end(); ++c2){
		for(TermInfo::const_iterator n1 = uppers1_.begin(); n1 != uppers1_.end(); ++n1){
			Stub const lower = vec3_to_stub(x2,c2->second);
			Stub const upper = vec3_to_stub(x1,n1->second);
			ndumped += dump_loophash_linkers( lower, upper, loop_hash_library_, loopsizes_, lookup_radius_ );
		}}
		if( ndumped > 0 ){
			Pose tmp1(pose1_),tmp2(pose2_);
			xform_pose(tmp1,x1);
		 	xform_pose(tmp2,x2);
		 	tmp1.dump_pdb("comp1.pdb");
		 	tmp2.dump_pdb("comp2.pdb");	 	
			utility_exit_with_message("TEST LINKER");
		}
	}

	return 10.0*lkscore;
}

// void
// LinkerScore::dump_linkers(
// 	Stub const & x1,
// 	Stub const & x2
// ) const {
// 	using namespace basic::options::OptionKeys;
// 	using namespace protocols::loophash;
// 	using namespace numeric::geometry::hashing;
// 	using namespace core::kinematics;

// 	for(TermInfo::const_iterator c1 = lowers1_.begin(); c1 != lowers1_.end(); ++c1){		
// 	for(TermInfo::const_iterator n2 = uppers2_.begin(); n2 != uppers2_.end(); ++n2){
// 		dump_linkers( vec3_to_stub(x1,*c1), vec3_to_stub(x2,*n2), loop_hash_library_, loopsizes_ ) );
// 	}}
// 	for(TermInfo::const_iterator c2 = lowers2_.begin(); c2 != lowers2_.end(); ++c2){		
// 	for(TermInfo::const_iterator n1 = uppers1_.begin(); n1 != uppers1_.end(); ++n1){
// 		dump_linkers( vec3_to_stub(x2,*c2), vec3_to_stub(x1,*n1), loop_hash_library_, loopsizes_ ) );
// 	}}
// }


JointScore::JointScore(
	Scores scores,
	Reals weights
):
	scores_(scores),
	weights_(weights)
{
	if(scores_.size() != weights_.size()) utility_exit_with_message("bad score/weight");
}

void
JointScore::add_score(
	RigidScoreCOP score,
	Real weight
){
	scores_.push_back(score);
	weights_.push_back(weight);
}

core::Real
JointScore::score(
	Stub const & x1,
	Stub const & x2
) const {
	Real score = 0.0;
	Reals::const_iterator w = weights_.begin();
	for(Scores::const_iterator s = scores_.begin(); s != scores_.end(); ++s,++w){
		score += (*s)->score(x1,x2) * (*w);
	}
	return score;
}



} // namespace sic_dock
} // namespace protocols
