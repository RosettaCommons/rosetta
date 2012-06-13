// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <protocols/sic_dock/util.hh>
#include <protocols/sic_dock/loophash_util.hh>
#include <protocols/sic_dock/RigidScore.hh>
#include <protocols/sic_dock/SICFast.hh>
#include <core/pose/util.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/HomogeneousTransform.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>

namespace protocols {
namespace sic_dock {

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

Vec3
get_leap_lower_stub(
	core::pose::Pose const & pose,
	Size ir
){
	Vec3 lower;
	lower.a = pose.residue(ir  ).xyz("CA");
	lower.b = pose.residue(ir  ).xyz( "N");
	lower.c = pose.residue(ir-1).xyz( "C");
	return lower;
}

Vec3
get_leap_upper_stub(
	core::pose::Pose const & pose,
	Size ir
){
	Vec3 upper;
	upper.a = pose.residue(ir+1).xyz("CA");
	upper.b = pose.residue(ir+1).xyz( "N");
	upper.c = pose.residue(ir+1).xyz( "H");
	return upper;
}

core::kinematics::Stub vec3_to_stub(Vec3 const & v3){
	return Stub(v3.a,v3.b,v3.c);
}

Stub vec3_to_stub(core::kinematics::Stub const & xform, Vec3 const & v3){
	return Stub( xform.local2global(v3.a), xform.local2global(v3.b), xform.local2global(v3.c) );
}

void
get_termini_from_pose(
	core::pose::Pose const & pose,
	Size ir,
	TermInfo & lowers,
	TermInfo & uppers
){
	if( pose.residue(ir).is_upper_terminus() && 
	    ir > 1 && 
	   !pose.residue(ir  ).is_lower_terminus() && 
	   // !pose.residue(ir-1).is_lower_terminus() &&
	   !pose.residue(ir-1).is_upper_terminus() &&
	    pose.residue(ir-1).is_protein()
	){
		lowers.push_back( std::make_pair(ir,get_leap_lower_stub(pose,ir)) );
		// TR << "add lower " << ir << "/" <<pose.n_residue() <<" "<< get_leap_lower_stub(pose,ir) << endl;
	}
	if( pose.residue(ir).is_lower_terminus() &&
	   !pose.residue(ir).is_upper_terminus() &&
	    ir < pose.n_residue() &&
	   !pose.residue(ir+1).is_lower_terminus() &&		   
	   !pose.residue(ir+1).is_upper_terminus() &&
	    pose.residue(ir+1).is_protein()
	){
		uppers.push_back( std::make_pair(ir,get_leap_upper_stub(pose,ir)) );
		// TR << "add upper " << ir << "/" <<pose.n_residue() <<" "<< get_leap_upper_stub(pose,ir) << endl;
	}
}
void
get_termini_from_pose(
	core::pose::Pose const & pose,
	TermInfo & lowers,
	TermInfo & uppers
){
	for(Size ir = 1; ir <= pose.n_residue(); ++ir){
		get_termini_from_pose(pose,ir,uppers,lowers);
	}
}

numeric::geometry::hashing::Real6
get_leap_6dof(
	Stub const & lower,
	Stub const & upper
){
	core::kinematics::RT leap(lower,upper);
	numeric::HomogeneousTransform< Real > ht( leap.get_rotation(), leap.get_translation() );
	numeric::xyzVector < Real > euler_angles =  ht.euler_angles_rad();
	numeric::geometry::hashing::Real6 rt_6;
	rt_6[1] = leap.get_translation().x();
	rt_6[2] = leap.get_translation().y();
	rt_6[3] = leap.get_translation().z();
	rt_6[4] = euler_angles.x()*180.0/numeric::constants::d::pi;
	rt_6[5] = euler_angles.y()*180.0/numeric::constants::d::pi;
	rt_6[6] = euler_angles.z()*180.0/numeric::constants::d::pi;	
	return rt_6;
}

Size
count_linkers(
	Stub const & lower,
	Stub const & upper,
	protocols::loophash::LoopHashLibraryOP loop_hash_library,
	Sizes const & loopsizes,
	core::Size radius
){
	if( loopsizes.size() == 0 ) return 0;

	Real dist2 = core::kinematics::RT(lower,upper).get_translation().length_squared();
	if( dist2 >= (10.0*Real(loopsizes.back())*Real(loopsizes.back())) ) return 0;

	numeric::geometry::hashing::Real6 rt_6 = get_leap_6dof(lower,upper);

	Size count0=0;
	for(Sizes::const_iterator i = loopsizes.begin(); i != loopsizes.end(); ++i){
		if( dist2 < (10.0*Real(*i)*Real(*i)) ){
			count0 += loop_hash_library->gethash(*i).radial_count(radius,rt_6);
		} // else too far, don't bother lookup
	}
	return (Real)count0;
}

core::Size
dump_loophash_linkers(
	Stub const & lower,
	Stub const & upper,
	// core::pose::Pose const & pose1,
	// core::pose::Pose const & pose2,	
	protocols::loophash::LoopHashLibraryOP loop_hash_library,
	Sizes const & loopsizes,
	Size radius
){
	// TermInfo lowers1,uppers1,lowers2,uppers2;
	// get_termini_from_pose(pose1_,uppers1_,lowers1_);
	// get_termini_from_pose(pose2_,uppers2_,lowers2_);


	protocols::loophash::BackboneDB const & bbdb_ = loop_hash_library->backbone_database();
	protocols::loophash::BackboneSegment backbone_;
	core::Size ndumped = 0;
	for(Sizes::const_iterator ils = loopsizes.begin(); ils != loopsizes.end(); ++ils){
		Size loopsize(*ils);
	
		protocols::loophash::LoopHashMap & hashmap( loop_hash_library->gethash(loopsize) );

		numeric::geometry::hashing::Real6 rt_6 = get_leap_6dof(lower,upper);
		if( rt_6[1]*rt_6[1]+rt_6[2]*rt_6[2]+rt_6[3]*rt_6[3] > 10.0*Real(loopsize)*Real(loopsize) ) continue;
		std::vector < core::Size > leap_index_list;
		for(Sizes::const_iterator i = loopsizes.begin(); i != loopsizes.end(); ++i){
			hashmap.radial_lookup(radius,rt_6,leap_index_list);
		}

		std::vector< protocols::loophash::BackboneSegment > bs_vec_;
		for( std::vector < core::Size >::const_iterator itx = leap_index_list.begin(); itx != leap_index_list.end(); ++itx ){
			core::Size bb_index = *itx;
			protocols::loophash::LeapIndex cp = hashmap.get_peptide( bb_index );
			bbdb_.get_backbone_segment( cp.index, cp.offset , loopsize , backbone_ );
			bs_vec_.push_back( backbone_ );
		}

		Pose tmp;
		{
			core::chemical::ResidueTypeSetCAP rs = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
			for(Size i = 1; i <= loopsize+1; ++i){
				core::conformation::ResidueOP new_rsd( NULL );
				new_rsd = core::conformation::ResidueFactory::create_residue( rs->name_map("ALA") );
				// cout << "apending residue " << new_rsd->name() << std::endl;
				if(1==i) tmp.append_residue_by_jump( *new_rsd, 1 );
				else     tmp.append_residue_by_bond( *new_rsd, true );
				tmp.set_phi  ( tmp.n_residue(), 180.0 );
				tmp.set_psi  ( tmp.n_residue(), 180.0 );
				tmp.set_omega( tmp.n_residue(), 180.0 );
			}
			// tmp.dump_pdb("test.pdb");
		}

		int count = 0;
		for ( std::vector<protocols::loophash::BackboneSegment>::const_iterator i = bs_vec_.begin(), ie = bs_vec_.end(); i != ie; ++i) {
			std::vector<core::Real> phi   = i->phi();
			std::vector<core::Real> psi   = i->psi();
			std::vector<core::Real> omega = i->omega();
			Size seg_length = (*i).length();
			for ( Size i = 0; i < seg_length; i++){
				Size ires = 2+i;  // this is terrible, due to the use of std:vector.  i has to start from 0, but positions offset by 1.
				if (ires > tmp.total_residue() ) break;
				tmp.set_phi  ( ires, phi[i]  );
				tmp.set_psi  ( ires, psi[i]  );
				tmp.set_omega( ires, omega[i]);
			}
			Stub s = vec3_to_stub(get_leap_lower_stub(tmp,2));
			xform_pose_rev(tmp,s);
			xform_pose(tmp,lower);
			// tmp.delete_polymer_residue(tmp.n_residue());
			tmp.dump_pdb("test_lh_"+ObjexxFCL::string_of(loopsize)+"_"+ObjexxFCL::string_of(++count)+".pdb");
			++ndumped;
		}
	}
	return ndumped;
}


inline Stub multstubs(Stub const & a, Stub const & b){
	return Stub( a.M*b.M, a.M*b.v+a.v );
}
inline Stub invstub(Stub const & a){
	Mat const MR = a.M.transposed();
	return Stub( MR, MR * -a.v );
}

Real linker_count2score(Size count){
	return 0.8*sqrt(count)+0.2*Real(count);
}




} // sic_dock
} // protocols
