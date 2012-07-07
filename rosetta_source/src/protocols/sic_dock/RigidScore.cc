// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:

#include <protocols/sic_dock/RigidScore.hh>
#include <protocols/sic_dock/xyzStripeHashPoseWithMeta.hh>
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

using platform::Size;
using platform::Real;
using numeric::min;
using core::id::AtomID;
using std::cout;
using std::endl;
typedef platform::Real Real;
typedef platform::Size Size;
typedef core::pose::Pose Pose;
typedef core::kinematics::Stub Stub;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
typedef utility::vector1<Vec> Vecs;
typedef utility::vector1<Real> Reals;
typedef utility::vector1<Size> Sizes;
typedef utility::vector1<Stub> Stubs;
typedef utility::vector1<RigidScoreCOP> Scores;

inline Real dist_score( Real const & sqdist, Real const & start, Real const & stop ) {
	if( sqdist > stop*stop ) {
		return 0.0;
	} else if( sqdist < start*start ) {
		return 1.0;
	} else {
		Real dist = sqrt( sqdist );
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

struct CBScoreVisitor {
	Real const clash_dis,contact_dis;
	Real score;
	CBScoreVisitor(
		Real const & clash_dis_in,
		Real const & contact_dis_in
	): 
		clash_dis(clash_dis_in),
		contact_dis(contact_dis_in),
		score(0.0)
	{}
	inline
	void
	visit(
		numeric::xyzVector<Real> const & /*v*/,
		                   Real  const & vm,
		numeric::xyzVector<Real> const & /*c*/, 
		                   Real  const & cm,
		                   Real  const & d2
	){
		score += dist_score(d2, clash_dis, contact_dis ) * vm * cm;
	}
};

platform::Real
CBScore::score(
	Stub const & x1,
	Stub const & x2
) const {
	Stub const xhp(multstubs(invstub(hash_pose1_?x1:x2),hash_pose1_?x2:x1));
	CBScoreVisitor hash_visitor(clash_dis_,contact_dis_);
	Reals::const_iterator iwb = weights_.begin();
	for(Vecs::const_iterator i = points_.begin(); i != points_.end(); ++i,++iwb) {
		xyzhash_.visit(xhp.local2global(*i),*iwb,hash_visitor);
	}
	return hash_visitor.score;
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



platform::Real
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

platform::Real
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
