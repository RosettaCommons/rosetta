// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <protocols/sic_dock/xyzStripeHashPose.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID_Map.hh>

#include <platform/types.hh>

namespace protocols {
namespace sic_dock {


	xyzStripeHashPose::xyzStripeHashPose(
		double radius
	):
		numeric::geometry::hashing::xyzStripeHash(radius),
		initialized_(false)
	{}

	xyzStripeHashPose::xyzStripeHashPose(
		core::pose::Pose const & p,
		PoseCoordPickMode m,
		double radius
	):
		numeric::geometry::hashing::xyzStripeHash(radius),
		initialized_(false)
	{
		add_pose(p,m);
		init_posehash();
	}

	xyzStripeHashPose::xyzStripeHashPose(
		core::pose::Pose const & p,
		core::id::AtomID_Map<double> const & amap,
		double radius
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
		core::id::AtomID_Map<double> const amap
	){
		if(initialized_) utility_exit_with_message("already initialized!");
		using core::id::AtomID;
		for(int ir = 1; ir <= (int)p.n_residue(); ++ir) {
			for(int ia = 1; ia <= (int)amap.n_atom(ir); ia++) {
				if(amap[AtomID(ia,ir)] > 0) {
					balls_.resize(balls_.size()+1);
					balls_.back().x() = p.xyz(AtomID(ia,ir)).x();
					balls_.back().y() = p.xyz(AtomID(ia,ir)).y();
					balls_.back().z() = p.xyz(AtomID(ia,ir)).z();
					balls_.back().resid_  = ir;
					balls_.back().atomno_ = ia;
					balls_.back().radius( p.residue(ir).atom_type(ia).lj_radius() );
				}
			}
		}
	}

	void
	xyzStripeHashPose::add_pose(
		core::pose::Pose const & p,
		PoseCoordPickMode m
	){
		using namespace core::id;
		AtomID_Map<double> amap;
		core::pose::initialize_atomid_map(amap,p,0.0);
		for(int ir = 1; ir <= (int)p.n_residue(); ++ir) {
			core::conformation::Residue const & r(p.residue(ir));
	    	if(CA==m){
				if(r.has("CA")) amap[AtomID(p.residue(ir).atom_index("CA"),ir)] = 1.0;
			} else if(CB==m) {
				if(r.has("CB")) amap[AtomID(p.residue(ir).atom_index("CB"),ir)] = 1.0;
			} else if(NBR==m) {
				               amap[AtomID(r.nbr_atom()                   ,ir)] = 1.0;
			} else if(BB==m) {
				if(r.has( "N")) amap[AtomID(p.residue(ir).atom_index( "N"),ir)] = 1.0;
				if(r.has("CA")) amap[AtomID(p.residue(ir).atom_index("CA"),ir)] = 1.0;
				if(r.has( "C")) amap[AtomID(p.residue(ir).atom_index( "C"),ir)] = 1.0;
				if(r.has( "O")) amap[AtomID(p.residue(ir).atom_index( "O"),ir)] = 1.0;
				if(r.has("CB")) amap[AtomID(p.residue(ir).atom_index("CB"),ir)] = 1.0;
			} else if(NCAC==m) {
				if(r.has( "N")) amap[AtomID(p.residue(ir).atom_index( "N"),ir)] = 1.0;
				if(r.has("CA")) amap[AtomID(p.residue(ir).atom_index("CA"),ir)] = 1.0;
				if(r.has( "C")) amap[AtomID(p.residue(ir).atom_index( "C"),ir)] = 1.0;
			} else if(NCO==m) {
				if(r.has( "N")) amap[AtomID(p.residue(ir).atom_index( "N"),ir)] = 1.0;
				if(r.has( "C")) amap[AtomID(p.residue(ir).atom_index( "C"),ir)] = 1.0;
				if(r.has( "O")) amap[AtomID(p.residue(ir).atom_index( "O"),ir)] = 1.0;
			} else {
				int natom = (ALL==m) ? r.natoms() : r.nheavyatoms();
				for(int ia = 1; ia <= natom; ++ia) {
					core::id::AtomID const aid(ia,ir);
					amap[aid] = 1.0;
				}
			}
		}
		add_pose(p,amap);
	}

	void
	xyzStripeHashPose::init_posehash(){
		if(initialized_) utility_exit_with_message("can't init twice!");
		init(balls_);
		initialized_ = true;
	}

} // namespace sic_dock
} // namespace protocols

