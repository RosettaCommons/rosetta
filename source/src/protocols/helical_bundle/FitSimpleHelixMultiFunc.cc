// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/FitSimpleHelixMultiFunc.cc
/// @brief  Multifunction class implementation for fitting a simple (straight) helix to the Crick parameters.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

/// Unit headers
#include <protocols/helical_bundle/FitSimpleHelixMultiFunc.hh>

/// Package headers
#include <numeric/crick_equations/HelixParams.hh>

/// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <numeric/model_quality/rms.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/optimization/Multifunc.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/util.tmpl.hh>

/// Utility headers
#include <cmath>
#include <utility>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

static basic::Tracer TR("protocols.helical_bundle.FitSimpleHelixMultiFunc");

namespace protocols {
namespace helical_bundle {

using namespace core;
using namespace core::optimization;

FitSimpleHelixMultiFunc::~FitSimpleHelixMultiFunc() = default;

FitSimpleHelixMultiFunc::FitSimpleHelixMultiFunc(
	core::pose::Pose const &pose,
	std::string atom_name,
	core::Size const first_res_index,
	core::Size const res_per_repeat,
	core::Size const start_index,
	core::Size const end_index,
	core::Size const minimization_mode,
	core::Real const &rms_offset
) :
	pose_( pose ),
	atom_name_(std::move( atom_name )),
	first_res_index_( first_res_index ),
	residues_per_repeat_( res_per_repeat ),
	start_index_(start_index),
	end_index_(end_index),
	minimization_mode_(minimization_mode),
	rms_offset_(rms_offset)
{
	runtime_assert_string_msg( start_index_ > 0, "In FitSimpleHelixMultiFunc constructor function: the starting index is out of range." );
	runtime_assert_string_msg( first_res_index_ >= start_index_ && first_res_index_ <= end_index_, "In FitSimpleHelixMultiFunc constructor function: the first residue index is not in the residue range specified." );
	runtime_assert_string_msg( end_index_ <= pose_.size(), "In FitSimpleHelixMultiFunc constructor function: the ending index is past the end of the pose." );
	runtime_assert_string_msg( start_index_ + 1 < end_index_, "In FitSimpleHelixMultiFunc constructor function: the starting index must be less than the ending index, and they must define at least three residues." );
	runtime_assert_string_msg( minimization_mode_==0 || minimization_mode_==1, "In FitSimpleHelixMultiFunc constructor function: invalid value for minimization_mode_.");
	runtime_assert_string_msg( residues_per_repeat_ > 0, "In FitSimpleHelixMultiFunc constructor function: there must be at least one residue per repeating unit." );

	for ( core::Size ir=first_res_index_; ir<=end_index_; ir += residues_per_repeat_ ) {
		runtime_assert_string_msg(pose.residue(ir).has(atom_name_), "In FitSimpleHelixMultiFunc constructor function: all helix residues must have the named atom (\"" + atom_name_ + "\")." );
	}

}

Real
FitSimpleHelixMultiFunc::operator ()( Multivec const & vars ) const
{
	//The vars Multivec contains the following:
	// vars[1] = r1 (helix radius)
	// vars[2] = omega1 (twist per residue)
	// vars[3] = dz1 (rise per residue)
	// vars[4] = delta_omega1 (omega offset for atoms other than CA)
	// vars[5] = delta_z1 (z offset for atoms other than CA)

	core::Real r1=0.0, omega1=0.0, dz1=0.0, delta_omega1=0.0, delta_z1=0.0;
	vars_to_params( vars, r1, omega1, dz1, delta_omega1, delta_z1 );

	core::Real rms = 0.0;

	core::Real t = -1.0*static_cast<core::Real>(end_index_ - start_index_ + 1) / 2.0;
	t += static_cast<core::Real>(first_res_index_-start_index_);

	//Set the helix conformation and store x,y,z coordinates of appropriate atoms from the pose:
	utility::vector1< Vector > p1_coords, p2_coords;
	for ( core::Size ir = first_res_index_; ir <= end_index_; ir += residues_per_repeat_ ) { //Loop through all residue indices
		p1_coords.push_back( numeric::crick_equations::xyz( r1, omega1, t, dz1, (minimization_mode_==1 ? delta_omega1 : 0.0  ), (minimization_mode_==1 ? delta_z1 : 0.0) ) );
		p2_coords.push_back( pose_.residue(ir).xyz(atom_name_) );
		t += static_cast<core::Real>( residues_per_repeat_ );
	}

	//Add an arbitrary offset to make alignment work better:
	for ( core::Size i=1, imax=p1_coords.size(); i<=imax; ++i ) {
		p1_coords[i].x() += 5;
		p1_coords[i].y() += 6;
		p1_coords[i].z() += 7;
		// p2_coords[i].x() += 5;
		// p2_coords[i].y() += 6;
		// p2_coords[i].z() += 7;
	}

	//This function needs to be super-efficient, since it will be evaluated over and over during a line search.
	//I don't like the fact that the function in model_quality returns rms instead of rms squared.  It means unnecessarily
	//taking a square root and re-squaring.

	if ( minimization_mode_==0 ) {
		//Calculate RMS with alignment.
		rms=numeric::model_quality::calc_rms( p1_coords, p2_coords );
		rms = rms*rms;
	} else { //if minimization_mode_==1
		//Calculate RMS with no alignment.
		rms=0.0;
		for ( core::Size i=1, imax=p1_coords.size(); i<=imax; ++i ) {
			rms += pow( p1_coords[i].x()-p2_coords[i].x() , 2 );
			rms += pow( p1_coords[i].y()-p2_coords[i].y() , 2 );
			rms += pow( p1_coords[i].z()-p2_coords[i].z() , 2 );
		}
	}

	if ( TR.Debug.visible() ) TR.Debug << "r1=" << r1 << " omega1=" << omega1 << " dz1=" << dz1 << " delta_omega1=" << delta_omega1 << " delta_z1=" << delta_z1 << " rms_sq=" << rms << std::endl;

	return rms;
}

void
FitSimpleHelixMultiFunc::dfunc( Multivec const & vars, Multivec & dE_dvars ) const
{
	core::Real r1=0.0, omega1=0.0, dz1=0.0, delta_omega1=0.0, delta_z1=0.0;
	vars_to_params( vars, r1, omega1, dz1, delta_omega1, delta_z1 );

	if ( dE_dvars.size() != vars.size() ) dE_dvars.resize(vars.size());

	//A temporary pose to store the ideal helix generated from the Crick equations:
	core::pose::Pose pose_copy(pose_);

	core::id::AtomID_Map< core::id::AtomID > amap;
	core::pose::initialize_atomid_map(amap, pose_copy, core::id::AtomID::BOGUS_ATOM_ID());

	core::Real t = -1.0*static_cast<core::Real>(end_index_ - start_index_ + 1)/2.0;
	t += static_cast<core::Real>(first_res_index_-start_index_);

	numeric::xyzVector <core::Real> offsetvect;
	offsetvect.x() = 5;
	offsetvect.y() = 6;
	offsetvect.z() = 7;
	for ( core::Size ir = first_res_index_; ir <= end_index_; ir += residues_per_repeat_ ) {
		pose_copy.set_xyz( core::id::NamedAtomID(atom_name_,ir), numeric::crick_equations::xyz( r1, omega1 , t, dz1, (minimization_mode_==1 ? delta_omega1 : 0.0 ), (minimization_mode_==1 ? delta_z1 : 0.0) ) + offsetvect);
		amap[core::id::AtomID(pose_copy.residue(ir).atom_index(atom_name_),ir)] = core::id::AtomID(pose_.residue(ir).atom_index(atom_name_),ir);
		t += static_cast<core::Real>( residues_per_repeat_ );
	}

	//Superimpose the ideal helix on the original helix:
	if ( minimization_mode_==0 ) core::scoring::superimpose_pose( pose_copy, pose_, amap, rms_offset_, true );

	//Accumulators for the derivatives:
	core::Real dE_dr1 = 0.0;
	core::Real dE_domega1 = 0.0;
	core::Real dE_ddz1 = 0.0;
	core::Real dE_ddelta_omega1 = 0.0;
	core::Real dE_ddelta_z1 = 0.0;

	//Reset t to the first index
	t = -1.0*static_cast<core::Real>(end_index_ - start_index_ + 1)/2.0;
	t += static_cast<core::Real>(first_res_index_-start_index_);

	//Calculate the derivatives:
	for ( core::Size ir = first_res_index_; ir <= end_index_; ir += residues_per_repeat_ ) {
		core::Real const x = pose_copy.residue(ir).xyz(atom_name_).x();
		core::Real const y = pose_copy.residue(ir).xyz(atom_name_).y();
		core::Real const z = pose_copy.residue(ir).xyz(atom_name_).z();
		core::Real const X = pose_.residue(ir).xyz(atom_name_).x();
		core::Real const Y = pose_.residue(ir).xyz(atom_name_).y();
		core::Real const Z = pose_.residue(ir).xyz(atom_name_).z();
		core::Real const xdiff = 2*(x-X);
		core::Real const ydiff = 2*(y-Y);
		core::Real const zdiff = 2*(z-Z);
		dE_dr1 += xdiff*numeric::crick_equations::dx_dr1(r1,omega1,t, dz1, (minimization_mode_==1 ? delta_omega1 : 0.0 ), (minimization_mode_==1 ? delta_z1 : 0.0 ));
		dE_dr1 += ydiff*numeric::crick_equations::dy_dr1(r1,omega1,t, dz1, (minimization_mode_==1 ? delta_omega1 : 0.0 ), (minimization_mode_==1 ? delta_z1 : 0.0 ));
		dE_dr1 += zdiff*numeric::crick_equations::dz_dr1(r1,omega1,t, dz1, (minimization_mode_==1 ? delta_omega1 : 0.0 ), (minimization_mode_==1 ? delta_z1 : 0.0 ));
		if ( minimization_mode_==0 ) {
			dE_domega1 += xdiff*numeric::crick_equations::dx_domega1(r1,omega1,t, dz1, 0.0, 0.0);
			dE_domega1 += ydiff*numeric::crick_equations::dy_domega1(r1,omega1,t, dz1, 0.0, 0.0);
			dE_domega1 += zdiff*numeric::crick_equations::dz_domega1(r1,omega1,t, dz1, 0.0, 0.0);
			dE_ddz1 += xdiff*numeric::crick_equations::dx_ddz1(r1,omega1,t, dz1, 0.0, 0.0);
			dE_ddz1 += ydiff*numeric::crick_equations::dy_ddz1(r1,omega1,t, dz1, 0.0, 0.0);
			dE_ddz1 += zdiff*numeric::crick_equations::dz_ddz1(r1,omega1,t, dz1, 0.0, 0.0);
		} else { //if minimization_mode_==1
			dE_ddelta_omega1 += xdiff*numeric::crick_equations::dx_ddelta_omega1(r1,omega1,t, dz1, delta_omega1, delta_z1);
			dE_ddelta_omega1 += ydiff*numeric::crick_equations::dy_ddelta_omega1(r1,omega1,t, dz1, delta_omega1, delta_z1);
			dE_ddelta_omega1 += zdiff*numeric::crick_equations::dz_ddelta_omega1(r1,omega1,t, dz1, delta_omega1, delta_z1);
			dE_ddelta_z1 += xdiff*numeric::crick_equations::dx_ddelta_z1(r1,omega1,t, dz1, delta_omega1, delta_z1);
			dE_ddelta_z1 += ydiff*numeric::crick_equations::dy_ddelta_z1(r1,omega1,t, dz1, delta_omega1, delta_z1);
			dE_ddelta_z1 += zdiff*numeric::crick_equations::dz_ddelta_z1(r1,omega1,t, dz1, delta_omega1, delta_z1);
		}
		t += static_cast<core::Real>( residues_per_repeat_ );
	}

	//Scale the derivatives down:
	//dE_dr1 *= 0.01;
	//dE_domega1 *= 0.01;
	//dE_ddz1 *= 0.01;
	//dE_ddelta_omega1 *=0.01;
	//dE_ddelta_z1 *=0.01;

	if ( TR.Debug.visible() ) TR.Debug << "dE_dr1=" << dE_dr1 << "   dE_domega1=" << dE_domega1 << "    dE_ddz1=" << dE_ddz1 << "    dE_ddelta_omega1=" << dE_ddelta_omega1 << "    dE_ddelta_z1=" << dE_ddelta_z1 << std::endl;
	//pose_copy.dump_pdb("temp2.pdb"); //DELETE ME!

	//Copy the derivatives to the derivative vector:
	params_derivs_to_vars(dE_dvars, dE_dr1, dE_domega1,dE_ddz1, dE_ddelta_omega1, dE_ddelta_z1);

	return;
}

/// @details Useful debugging code that can be re-enabled by changing the boolean
/// variables at the top of this function.
void
FitSimpleHelixMultiFunc::dump( Multivec const &/*vars*/, Multivec const &/*vars2*/ ) const {
	return;
}

/*************************************************************
PRIVATE MEMBER FUNCTIONS
*************************************************************/

void FitSimpleHelixMultiFunc::vars_to_params( Multivec const &vars, core::Real &r1, core::Real &omega1, core::Real &dz1, core::Real &delta_omega1, core::Real &delta_z1 ) const
{
	r1 = vars[1];
	omega1 = vars[2];
	dz1 = vars[3];
	delta_omega1=vars[4];
	delta_z1=vars[5];
	return;
}

/// @brief Convert the Crick parameter derivatives to the derivative Multivec.
///
void FitSimpleHelixMultiFunc::params_derivs_to_vars( Multivec &deriv_vars, core::Real const &dE_dr1, core::Real const &dE_domega1, core::Real const &dE_ddz1, core::Real const &dE_ddelta_omega1, core::Real const &dE_ddelta_z1 ) const
{
	deriv_vars[1]=dE_dr1;
	deriv_vars[2]=dE_domega1;
	deriv_vars[3]=dE_ddz1;
	deriv_vars[4]=dE_ddelta_omega1;
	deriv_vars[5]=dE_ddelta_z1;
	return;
}


} // namespace helical_bundle
} // namespace protocols

