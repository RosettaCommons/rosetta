// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <core/pose/xyzStripeHashPose.hh>

#include <basic/options/keys/sicdock.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID_Map.hh>

#include <platform/types.hh>

namespace core {
namespace pose {

	core::id::AtomID_Map<platform::Real> make_atom_map( core::pose::Pose const & p, PoseCoordPickMode m ){
		using namespace core::id;
		AtomID_Map<platform::Real> amap;
		core::pose::initialize_atomid_map(amap,p,0.0);
		for(int ir = 1; ir <= (int)p.n_residue(); ++ir) {
			core::conformation::Residue const & r(p.residue(ir));
			using namespace basic::options;
			// if( option[OptionKeys::sicdock::rose::ignore_residues_far_from_origin]()+10.0 < r.nbr_atom_xyz().length() ) continue;

			switch(m){
	    	case PoseCoordPickMode_CA:
				if(r.has("CA")) amap[AtomID(p.residue(ir).atom_index("CA"),ir)] = 1.0;
				break;
			case PoseCoordPickMode_CB:
				if(r.has("CB")) amap[AtomID(p.residue(ir).atom_index("CB"),ir)] = 1.0;
				break;
			case PoseCoordPickMode_CB_else_CA:
				if(     r.has("CB")) amap[AtomID(p.residue(ir).atom_index("CB"),ir)] = 1.0;
				else if(r.has("CA")) amap[AtomID(p.residue(ir).atom_index("CA"),ir)] = 1.0;
				break;
			case PoseCoordPickMode_NBR:
				amap[AtomID(r.nbr_atom(),ir)] = 1.0;
				break;
			case PoseCoordPickMode_BB:
				if(r.has( "N")) amap[AtomID(p.residue(ir).atom_index( "N"),ir)] = 1.0;
				if(r.has("CA")) amap[AtomID(p.residue(ir).atom_index("CA"),ir)] = 1.0;
				if(r.has( "C")) amap[AtomID(p.residue(ir).atom_index( "C"),ir)] = 1.0;
				if(r.has( "O")) amap[AtomID(p.residue(ir).atom_index( "O"),ir)] = 1.0;
				if(r.has("CB")) amap[AtomID(p.residue(ir).atom_index("CB"),ir)] = 1.0;
				break;
			case PoseCoordPickMode_N_CA_C_CB:
				if(r.has( "N")) amap[AtomID(p.residue(ir).atom_index( "N"),ir)] = 1.0;
				if(r.has("CA")) amap[AtomID(p.residue(ir).atom_index("CA"),ir)] = 1.0;
				if(r.has( "C")) amap[AtomID(p.residue(ir).atom_index( "C"),ir)] = 1.0;
				if(r.has("CB")) amap[AtomID(p.residue(ir).atom_index("CB"),ir)] = 1.0;
				break;
			case PoseCoordPickMode_N_CA_C:
				if(r.has( "N")) amap[AtomID(p.residue(ir).atom_index( "N"),ir)] = 1.0;
				if(r.has("CA")) amap[AtomID(p.residue(ir).atom_index("CA"),ir)] = 1.0;
				if(r.has( "C")) amap[AtomID(p.residue(ir).atom_index( "C"),ir)] = 1.0;
				break;
			case PoseCoordPickMode_N_C_O:
				if(r.has( "N")) amap[AtomID(p.residue(ir).atom_index( "N"),ir)] = 1.0;
				if(r.has( "C")) amap[AtomID(p.residue(ir).atom_index( "C"),ir)] = 1.0;
				if(r.has( "O")) amap[AtomID(p.residue(ir).atom_index( "O"),ir)] = 1.0;
				break;
			case PoseCoordPickMode_HVY:
				for(Size ia = 1; ia <= r.nheavyatoms(); ++ia) {
					if(ia==4) continue;
					core::id::AtomID const aid(ia,ir);
					amap[aid] = 1.0;
				}
				break;
			case PoseCoordPickMode_HVY_IF_NP:
				for(Size ia = 1; ia <= r.nheavyatoms(); ++ia) {
					core::id::AtomID const aid(ia,ir);
					if(ia==4) continue;
					if( ia <= 5 || !r.is_polar()) amap[aid] = 1.0;
				}
				break;
			default:
				for(Size ia = 1; ia <= r.natoms(); ++ia) {
					core::id::AtomID const aid(ia,ir);
					amap[aid] = 1.0;
				}
				break;
			}
		}
		return amap;
	}

	xyzStripeHashPose::xyzStripeHashPose(
		platform::Real radius
	):
		numeric::geometry::hashing::xyzStripeHash(radius),
		initialized_(false)
	{}

	xyzStripeHashPose::xyzStripeHashPose(
		core::pose::Pose const & p,
		PoseCoordPickMode m,
		platform::Real radius
	):
		numeric::geometry::hashing::xyzStripeHash(radius),
		initialized_(false)
	{
		add_pose(p,m);
		init_posehash();
	}

	xyzStripeHashPose::xyzStripeHashPose(
		core::pose::Pose const & p,
		utility::vector1<int> const & resmap,
		PoseCoordPickMode m,
		platform::Real radius
	):
		numeric::geometry::hashing::xyzStripeHash(radius),
		initialized_(false)
	{
		id::AtomID_Map<platform::Real> amap = make_atom_map(p,m);
		using core::id::AtomID;
		for(int ir = 1; ir <= (int)p.n_residue(); ++ir) {
			if(resmap[ir] == 0) continue;
			for(int ia = 1; ia <= (int)amap.n_atom(ir); ia++) {
				if(amap[AtomID(ia,ir)] > 0) {
					balls_.resize(balls_.size()+1);
					balls_.back().x() = p.xyz(AtomID(ia,ir)).x();
					balls_.back().y() = p.xyz(AtomID(ia,ir)).y();
					balls_.back().z() = p.xyz(AtomID(ia,ir)).z();
					balls_.back().resid_  = resmap[ir];
					balls_.back().atomno_ = ia;
					balls_.back().radius( p.residue(ir).atom_type(ia).lj_radius() );
				}
			}
		}
		init_posehash();
	}


	xyzStripeHashPose::xyzStripeHashPose(
		core::pose::Pose const & p,
		core::id::AtomID_Map<platform::Real> const & amap,
		platform::Real radius
	):
		numeric::geometry::hashing::xyzStripeHash(radius),
		initialized_(false)
	{
		add_pose(p,amap);
		init_posehash();
	}

	void
	xyzStripeHashPose::add_pose(
		core::pose::Pose const & p,
		core::id::AtomID_Map<platform::Real> const & amap
	){
		if(initialized_) utility_exit_with_message("already initialized!");
		extract_pose_balls(p, balls_, amap);
	}

	void
	xyzStripeHashPose::add_pose(
		core::pose::Pose const & p,
		PoseCoordPickMode m
	){
		if(initialized_) utility_exit_with_message("already initialized!");
		extract_pose_balls(p, balls_, m);
	}

	void 
	xyzStripeHashPose::extract_pose_balls(
			core::pose::Pose const & p,
			utility::vector1<numeric::geometry::hashing::Ball> & balls,
			PoseCoordPickMode m
	){
		extract_pose_balls(p, balls, make_atom_map(p,m));
	}

	void 
	xyzStripeHashPose::extract_pose_balls(
			core::pose::Pose const & p,
			utility::vector1<numeric::geometry::hashing::Ball> & balls,
			core::id::AtomID_Map<platform::Real> const & amap
	){
		using core::id::AtomID;

		for(int ir = 1; ir <= (int)p.n_residue(); ++ir) {
			for(int ia = 1; ia <= (int)amap.n_atom(ir); ia++) {
				if(amap[AtomID(ia,ir)] > 0) {
					balls.resize(balls.size()+1);
					balls.back().x() = p.xyz(AtomID(ia,ir)).x();
					balls.back().y() = p.xyz(AtomID(ia,ir)).y();
					balls.back().z() = p.xyz(AtomID(ia,ir)).z();
					balls.back().resid_  = ir;
					balls.back().atomno_ = ia;
					balls.back().radius( p.residue(ir).atom_type(ia).lj_radius() );
			}
		}
	}
	}

	void
	xyzStripeHashPose::init_posehash(){
		if(initialized_) utility_exit_with_message("can't init twice!");
		init(balls_);
		initialized_ = true;
	}

} // namespace pose
} // namespace core

