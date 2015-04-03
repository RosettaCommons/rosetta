// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <boost/foreach.hpp>

namespace protocols {
namespace sic_dock {

protocols::loophash::LoopHashLibraryOP LinkerScore::loop_hash_library_ = NULL;

static thread_local basic::Tracer TR( "protocols.sic_dock.RigidScore" );

using core::Size;
using core::Real;
using numeric::min;
using core::id::AtomID;
using std::cout;
using std::endl;

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


std::ostream & operator<< ( std::ostream & out, Vec3 v ){
	out << v.a << " " << v.b << " " << v.c;
	return out;
}


void RigidScore::show(std::ostream & out                                      , int width) const { out << ObjexxFCL::format::RJ(width,type()); }
void RigidScore::show(std::ostream & out, Xforms const & x1, Xforms const & x2, int width) const { out << ObjexxFCL::format::F(width,3,score(x1,x2)); }


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
    points_(           get_CB_Vecs_from_pose( hash_pose1_?pose2:pose1 ) ),
    xyzhash_(
    	contact_dis,
    	hash_pose1_?pose1:pose2,
    	cb_weight_map_from_pose( hash_pose1_?pose1:pose2 )
    )
    // ,
    // pose1_(pose1),
    // pose2_(pose2)
{}

CBScore::CBScore(
	Pose const & pose1,
	Pose const & pose2,
	Real clash_dis,
	Real contact_dis,
	core::id::AtomID_Map<core::Real> const & weights1,
	core::id::AtomID_Map<core::Real> const & weights2
):
   	hash_pose1_( pose1_has_most_CBs(pose1,pose2) ),
	clash_dis_(clash_dis),
    contact_dis_(contact_dis),
    weights_( cb_weights_from_map( hash_pose1_?pose2:pose1, hash_pose1_?weights2:weights1 ) ),
    points_( get_CB_Vecs_from_map( hash_pose1_?pose2:pose1, hash_pose1_?weights2:weights1 ) ),
    xyzhash_(
    	contact_dis,
    	hash_pose1_?   pose1:   pose2,
    	hash_pose1_?weights1:weights2
    )
    // ,
    // pose1_(pose1),
    // pose2_(pose2)
{}

struct CBScoreVisitor {
	float const clash_dis,contact_dis;
	float score;
	CBScoreVisitor(
		float const & clash_dis_in,
		float const & contact_dis_in
	):
		clash_dis(clash_dis_in),
		contact_dis(contact_dis_in),
		score(0.0)
	{}
	inline
	void
	visit(
		numeric::xyzVector<float> const & /*v*/,
		                   float  const & vm,
		numeric::xyzVector<float> const & /*c*/,
		                   float  const & cm,
		                   float  const & d2
	){
		score += CBScore_dist_score(d2, 6.0f, contact_dis ) * vm * cm;
	}
};

core::Real
CBScore::score(
	Xforms const & x1s,
	Xforms const & x2s
) const {
	CBScoreVisitor hash_visitor(clash_dis_,contact_dis_);
	BOOST_FOREACH(Xform const & x1,x1s){
		BOOST_FOREACH(Xform const & x2,x2s){
			Xform const xhp(multstubs(invstub(hash_pose1_?x1:x2),hash_pose1_?x2:x1));
			Reals::const_iterator iwb = weights_.begin();
			for(Vecs::const_iterator i = points_.begin(); i != points_.end(); ++i,++iwb) {
				xyzhash_.visit(xhp*(*i),*iwb,hash_visitor);
			}
		}
	}
	return hash_visitor.score;
}

LinkerScore::LinkerScore(
	Pose const & pose1,
	Pose const & pose2,
	Size max_loop_len,
	Size lookup_radius,
	std::string const & outtag
):
	loopsizes_( range(3,max_loop_len+1) ),
	lookup_radius_( lookup_radius ),
	pose1_(pose1),
	pose2_(pose2),
	outtag_(outtag)
{
	if( ! basic::options::option[ basic::options::OptionKeys::lh::db_path ].user() ){
		utility_exit_with_message("user must specify -lh:db_path for loop hash database");
	}

	if(!loop_hash_library_){
		loop_hash_library_ = protocols::loophash::LoopHashLibraryOP( new protocols::loophash::LoopHashLibrary ( loopsizes_, 1, 0 ) );
		TR << "loading loophash data" << std::endl;
		if(loopsizes_.size()) loop_hash_library_->load_mergeddb();
		TR << "done loading loophash data" << std::endl;
	}

	get_termini_from_pose(pose1_,uppers1_,lowers1_);
	get_termini_from_pose(pose2_,uppers2_,lowers2_);

	TR << "LinkerScore termini: " << lowers1_.size() << " " << uppers1_.size() << " " << lowers2_.size() << " " << uppers2_.size() << std::endl;
	// utility_exit_with_message("TEST LINK SCORE");
}


core::Real
LinkerScore::score(
	Xforms const & x1s,
	Xforms const & x2s
) const {
	using namespace basic::options::OptionKeys;
	using namespace protocols::loophash;
	using namespace numeric::geometry::hashing;
	using namespace core::kinematics;

	Real lkscore = 0.0;

	BOOST_FOREACH(Xform const & x1,x1s){
		BOOST_FOREACH(Xform const & x2,x2s){

		for(TermInfo::const_iterator c1 = lowers1_.begin(); c1 != lowers1_.end(); ++c1){
		for(TermInfo::const_iterator n2 = uppers2_.begin(); n2 != uppers2_.end(); ++n2){
			Xform const lower = vec3_to_stub(x1,c1->second);
			Xform const upper = vec3_to_stub(x2,n2->second);
			Size n = count_linkers( lower, upper, loop_hash_library_, loopsizes_, lookup_radius_ );
			Real s = linker_count2score(n);
			lkscore += s;
		}}
		for(TermInfo::const_iterator c2 = lowers2_.begin(); c2 != lowers2_.end(); ++c2){
		for(TermInfo::const_iterator n1 = uppers1_.begin(); n1 != uppers1_.end(); ++n1){
			Xform const lower = vec3_to_stub(x2,c2->second);
			Xform const upper = vec3_to_stub(x1,n1->second);
			Size n = count_linkers( lower, upper, loop_hash_library_, loopsizes_, lookup_radius_ );
			Real s = linker_count2score(n);
			lkscore += s;
		}}

		// if(lkscore > 3.0){
		// 	TR << "LinkerScore test " << lkscore << std::endl;
		//  	Size ndumped = 0;
		// 	for(TermInfo::const_iterator c1 = lowers1_.begin(); c1 != lowers1_.end(); ++c1){
		// 	for(TermInfo::const_iterator n2 = uppers2_.begin(); n2 != uppers2_.end(); ++n2){
		// 		Xform const lower = vec3_to_stub(x1,c1->second);
		// 		Xform const upper = vec3_to_stub(x2,n2->second);
		// 		ndumped += dump_loophash_linkers( lower, upper, loop_hash_library_, loopsizes_, lookup_radius_, outtag_ );
		// 	}}
		// 	for(TermInfo::const_iterator c2 = lowers2_.begin(); c2 != lowers2_.end(); ++c2){
		// 	for(TermInfo::const_iterator n1 = uppers1_.begin(); n1 != uppers1_.end(); ++n1){
		// 		Xform const lower = vec3_to_stub(x2,c2->second);
		// 		Xform const upper = vec3_to_stub(x1,n1->second);
		// 		ndumped += dump_loophash_linkers( lower, upper, loop_hash_library_, loopsizes_, lookup_radius_, outtag_ );
		// 	}}
		// 	if( ndumped > 0 ){
		// 		Pose tmp1(pose1_),tmp2(pose2_);
		// 		xform_pose(tmp1,x1);
		// 	 	xform_pose(tmp2,x2);
		// 	 	tmp1.dump_pdb(outtag_+"comp1.pdb");
		// 	 	tmp2.dump_pdb(outtag_+"comp2.pdb");
		// 	 	std::cerr << "dumped "	 	 << ndumped << std::endl;
		// 	 	throw 12345;
		// 		utility_exit_with_message("TEST LINKER");
		// 	}
		// }
	}
	}
	return 10.0*lkscore;
}

bool
LinkerScore::dump_linkers(
	Xform const & x1,
	Xform const & x2,
	std::string const & out_perfix
) const {
	// TR << "LinkerScore test " << lkscore << std::endl;
 	Size ndumped = 0;
	for(TermInfo::const_iterator c1 = lowers1_.begin(); c1 != lowers1_.end(); ++c1){
	for(TermInfo::const_iterator n2 = uppers2_.begin(); n2 != uppers2_.end(); ++n2){
		Xform const lower = vec3_to_stub(x1,c1->second);
		Xform const upper = vec3_to_stub(x2,n2->second);
		ndumped += dump_loophash_linkers( lower, upper, loop_hash_library_, loopsizes_, lookup_radius_, out_perfix );
	}}
	for(TermInfo::const_iterator c2 = lowers2_.begin(); c2 != lowers2_.end(); ++c2){
	for(TermInfo::const_iterator n1 = uppers1_.begin(); n1 != uppers1_.end(); ++n1){
		Xform const lower = vec3_to_stub(x2,c2->second);
		Xform const upper = vec3_to_stub(x1,n1->second);
		ndumped += dump_loophash_linkers( lower, upper, loop_hash_library_, loopsizes_, lookup_radius_, out_perfix );
	}}
	// if( ndumped > 0 ){
	// 	Pose tmp1(pose1_),tmp2(pose2_);
	// 	xform_pose(tmp1,x1);
	//  	xform_pose(tmp2,x2);
	//  	tmp1.dump_pdb(outtag_+"comp1.pdb");
	//  	tmp2.dump_pdb(outtag_+"comp2.pdb");
	//  	std::cerr << "dumped "	 	 << ndumped << std::endl;
	//  	throw 12345;
	// 	utility_exit_with_message("TEST LINKER");
	// }
	return ndumped > 0;
}


class NoPoseXYX_Func : public core::scoring::func::XYZ_Func {
 public:
 	boost::unordered_map<core::id::AtomID,Vec,AtomIDHashFunction> const & start_coords_;
 	Xform const & x1_;
	Xform const & x2_;
 	Pose const * dummy_pose_;
 	mutable Vec result_hack_awful_1;
 	mutable Vec result_hack_awful_2;
	NoPoseXYX_Func(
		boost::unordered_map<core::id::AtomID,Vec,AtomIDHashFunction> const & start_coords,
		Xform const & x1,
		Xform const & x2
	) :
		start_coords_(start_coords), x1_(x1), x2_(x2), dummy_pose_(NULL)
	{}
	Vec const &	operator()( AtomID const & id ) const {
		if(id.rsd()>1000000){
			result_hack_awful_1 = x2_ * start_coords_.find(id)->second;
			return result_hack_awful_1;
		} else {
			result_hack_awful_2 = x1_ * start_coords_.find(id)->second;
			return result_hack_awful_2;
		}
	}

	virtual
	Residue const &
	residue( Size ) const {
		utility_exit_with_message("no Residue for NoPoseXYX_Func");
		return dummy_pose_->residue(1234567890);
	}

};


// NOTE!!!!
// my convention here is that AtomIDs for pose2 are residue number 1000000 + actual
//
	ConstraintSetScore::ConstraintSetScore(
		Pose const & pose1,
		Pose const & pose2,
		core::scoring::constraints::ConstraintSet const & cstset
	):
		// pose1_(pose1),
		// pose2_(pose2),
		// cstset_(cstset),
		csts_(cstset.get_all_constraints())
	{
		using namespace core::scoring::constraints;
		for(ConstraintCOPs::const_iterator i = csts_.begin(); i != csts_.end(); ++i){
			Constraint const & cst(**i);
			for(Size j = 1; j <= cst.natoms(); ++j){
				AtomID aid = cst.atom(j);
				Vec xyz = aid.rsd() > 1000000 ? pose2.xyz(AtomID(aid.atomno(),aid.rsd()-1000000)) : pose1.xyz(aid);
				start_coords_.insert(std::make_pair(aid,xyz));
			}
		}
	}
	core::Real ConstraintSetScore::score( Xforms const & x1s, Xforms const & x2s ) const {
		using namespace core::scoring::constraints;
		Real s = 0.0;
		BOOST_FOREACH(Xform const & x1,x1s){
			BOOST_FOREACH(Xform const & x2,x2s){
				NoPoseXYX_Func tmp_xyz_func(start_coords_,x1,x2);
				for(ConstraintCOPs::const_iterator i = csts_.begin(); i != csts_.end(); ++i){
					Constraint const & cst(**i);
					// cout << cst.atom(1) << " " << cst.atom(2) << " " << cst.dist(tmp_xyz_func) << " " << tmp_xyz_func(cst.atom(1)) << " " << tmp_xyz_func(cst.atom(2)) << endl;
					s += cst.get_func().func(cst.dist(tmp_xyz_func));
				}
			}
		}
		return s;
	}


// void
// LinkerScore::dump_linkers(
// 	Xform const & x1,
// 	Xform const & x2
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


JointScore::JointScore()
{}

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
	Xforms const & x1,
	Xforms const & x2
) const {
	Real score = 0.0;
	Reals::const_iterator w = weights_.begin();
	for(Scores::const_iterator s = scores_.begin(); s != scores_.end(); ++s,++w){
		if( fabs(*w) == 0 ) continue;
		score += (*s)->score(x1,x2) * (*w);
		// if(score < minscore_hack_ || score > maxscore_hack_ ) break;
	}
	return score;
}


void JointScore::show(std::ostream & out, int width) const {
	for(Scores::const_iterator s = scores_.begin(); s != scores_.end(); ++s){
		(*s)->show(out,width);
		out << " ";
	}
}
void JointScore::show(std::ostream & out, Xforms const & x1, Xforms const & x2, int width) const {
	for(Scores::const_iterator s = scores_.begin(); s != scores_.end(); ++s){
		(*s)->show(out,x1,x2,width);
		out << " ";
	}
}


} // namespace sic_dock
} // namespace protocols
