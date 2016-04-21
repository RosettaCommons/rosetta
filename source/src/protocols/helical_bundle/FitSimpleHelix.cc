// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/FitSimpleHelix.cc
/// @brief  Determines the Crick parameters describing a straight helix, given an input pose.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit Headers
#include <protocols/helical_bundle/FitSimpleHelix.hh>
#include <protocols/helical_bundle/FitSimpleHelixCreator.hh>
#include <protocols/helical_bundle/FitSimpleHelixMultiFunc.hh>
#include <numeric/crick_equations/HelixParams.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <core/scoring/rms_util.hh>
#include <core/pose/util.tmpl.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>

using basic::T;
using basic::Error;
using basic::Warning;

//static numeric::random::RandomGenerator RG(8713093);  // <- Magic number, do not change it!

namespace protocols {
namespace helical_bundle {

static THREAD_LOCAL basic::Tracer TR("protocols.helical_bundle.FitSimpleHelix");

std::string
FitSimpleHelixCreator::keyname() const
{
	return FitSimpleHelixCreator::mover_name();
}

protocols::moves::MoverOP
FitSimpleHelixCreator::create_mover() const {
	return protocols::moves::MoverOP( new FitSimpleHelix );
}

std::string
FitSimpleHelixCreator::mover_name()
{
	return "FitSimpleHelix";
}


/// @brief Creator for FitSimpleHelix mover.
FitSimpleHelix::FitSimpleHelix():
	//utility::pointer::ReferenceCount(),
	Mover("FitSimpleHelix"),
	r1_initial_(1.0),
	omega1_initial_(1.0),
	dz1_initial_(1.0),
	start_index_(0),
	end_index_(0),
	min_type_("lbfgs_armijo_nonmonotone"),
	min_tolerance_(0.00000001),
	reference_atom_("CA"),
	reference_residue_(1),
	residues_per_repeat_(1),
	r1_vals_output_(),
	omega1_val_output_(0.0),
	z1_val_output_(0.0),
	delta_omega1_vals_output_(),
	delta_z1_vals_output_(),
	r1_guesses_(),
	delta_omega1_guesses_(),
	delta_z1_guesses_()
{}


/// @brief Destructor for FitSimpleHelix mover.
FitSimpleHelix::~FitSimpleHelix() {}


/// @brief Clone operator to create a pointer to a fresh FitSimpleHelix object that copies this one.
protocols::moves::MoverOP FitSimpleHelix::clone() const {
	return protocols::moves::MoverOP( new FitSimpleHelix( *this ) );
}


/// @brief Fresh_instance operator to create a pointer to a fresh FitSimpleHelix object that does NOT copy this one.
protocols::moves::MoverOP FitSimpleHelix::fresh_instance() const {
	return protocols::moves::MoverOP( new FitSimpleHelix );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
/// @details At the end of this function, two things will have happened:
/// --The pose will be realigned to the ideal helix.
/// --The fitted helical parameters will have been stored in internal variables for the mover, and can be retrieved with a
/// call to get_crick_parameters().
void FitSimpleHelix::apply (core::pose::Pose & pose)
{
	runtime_assert_string_msg( start_index_ > 1 && start_index_ < pose.n_residue() && end_index_ > 1 && end_index_ < pose.n_residue(), "In FitSimpleHelix::apply(): The start and end residues must be within the pose's residue range, and cannot be the end residues." );
	runtime_assert_string_msg( (end_index_ - start_index_ + 1) % residues_per_repeat() == 0, "In FitSimpleHelix::apply(): The range of residues must represent an integer number of repeats." );
	runtime_assert_string_msg( (end_index_ - start_index_ + 1) / residues_per_repeat() > 2, "In FitSimpleHelix::apply(): At least three repeats must be fitted." );

	//Add arbitrary offset to the pose to make the superimposition work better:
	numeric::xyzVector <core::Real> offsetvect;
	offsetvect.x() = 5;
	offsetvect.y() = 6;
	offsetvect.z() = 7;
	for ( core::Size ir=1, irmax=pose.n_residue(); ir<=irmax; ++ir ) {
		for ( core::Size ia=1, iamax=pose.residue(ir).natoms(); ia<=iamax; ++ia ) {
			pose.set_xyz( core::id::AtomID(ia,ir), pose.xyz( core::id::AtomID(ia,ir) ) + offsetvect );
		}
	}

	//Count mainchain atoms in the residues making up the repeating unit:
	core::Size natoms(0);
	core::Size natoms_before_refres(0);
	for ( core::Size i=0, imax=residues_per_repeat(); i<imax; ++i ) {
		if ( i+1 < reference_residue() ) natoms_before_refres += pose.residue(start_index_+i).n_mainchain_atoms();
		natoms += pose.residue(start_index_+i).n_mainchain_atoms();
	}

	r1_vals_output_.resize( natoms, 0.0 );
	delta_omega1_vals_output_.resize( natoms, 0.0 );
	delta_z1_vals_output_.resize( natoms, 0.0 );

	//The indices of the reference residue and the reference atom:
	core::Size const reference_res_index( start_index_ + reference_residue() - 1 );
	core::Size const reference_atom_index( pose.residue( reference_res_index ).atom_index( reference_atom_ ) ); //Index in its residue.
	core::Size const absolute_reference_atom_index( reference_atom_index + natoms_before_refres ); //Index in the vector of all atoms to be fitted.

	{ //Scope 1: fitting the reference atom in the reference residue.

		runtime_assert_string_msg( reference_residue() <= residues_per_repeat(), "In FitSimpleHelix::apply(): The reference residue's index in the repeating unit is greater than the number of residues per repeating unit." );
		runtime_assert_string_msg( reference_res_index >= start_index_ && reference_res_index <= end_index_, "In FitSimpleHelix::apply(): The reference residue does not fall within the range of residues to be fitted." );
		runtime_assert_string_msg( reference_atom_index <= pose.residue( reference_res_index ).n_mainchain_atoms(), "In FitSimpleHelix::apply(): The reference atom is not a mainchain atom in the reference residue." );

		if ( TR.visible() ) TR << "Fitting reference atom " << reference_atom_ << " in reference residue (residue " << reference_res_index << ")." << std::endl;

		FitSimpleHelixMultiFunc multfunc(pose, reference_atom_, reference_res_index, residues_per_repeat(), start_index_, end_index_, 0, rms_offset_); //Make the multifunc that will be used to fit the reference atom.

		Multivec vars(5); //Change the number of DOFs appropriately and initialize.
		vars[1] = r1_initial_;
		vars[2] = omega1_initial_;
		vars[3] = dz1_initial_;
		vars[4] = 0.0;
		vars[5] = 0.0;

		core::optimization::MinimizerOptions minoptions( min_type_, min_tolerance_, false, false, false );
		core::optimization::Minimizer minimizer( multfunc, minoptions );

		//Actually run the minimizer to do the fitting.
		//(The multifunc is set up to trick the minimizer into acting as a nonlinear least-squares fitter for the Crick equations.  Clever, I know.)
		minimizer.run( vars );

		r1_vals_output_[absolute_reference_atom_index] = vars[1];
		omega1_val_output_ = vars[2];
		z1_val_output_ = vars[3];
		delta_omega1_vals_output_[absolute_reference_atom_index] = vars[4];
		delta_z1_vals_output_[absolute_reference_atom_index] = vars[5];

		if ( TR.visible() ) {
			TR << "Reference atom r1 = " << vars[1] << std::endl; //Delete me
			TR << "Reference atom omega1 = " << vars[2] << std::endl; //Delete me
			TR << "Reference atom dz1 = " << vars[3] << std::endl; //Delete me
			TR << "Reference atom delta_omega1 = " << vars[4] << std::endl; //Delete me
			TR << "Reference atom delta_z1 = " << vars[5] << std::endl; //Delete me
		}

		if ( TR.visible() ) TR << "Aligning pose to ideal helix." << std::endl;

		core::pose::Pose pose_copy = pose; //Make a copy of the pose.

		core::id::AtomID_Map< core::id::AtomID > amap;
		core::pose::initialize_atomid_map(amap, pose, core::id::BOGUS_ATOM_ID);

		core::Real t = -1.0*static_cast<core::Real>(end_index_ - start_index_ + 1)/2.0;
		t += static_cast<core::Real>( reference_res_index-start_index_ );

		for ( core::Size ir = reference_res_index; ir <= end_index_; ir += residues_per_repeat() ) {
			pose_copy.set_xyz( core::id::NamedAtomID(reference_atom_,ir), numeric::crick_equations::xyz( vars[1], vars[2] , t, vars[3], 0.0, 0.0 ) + offsetvect );
			amap[core::id::AtomID(pose.residue(ir).atom_index(reference_atom_),ir)] = core::id::AtomID(pose_copy.residue(ir).atom_index(reference_atom_),ir);
			t += static_cast<core::Real>( residues_per_repeat() );
		}

		//Superimpose the original helix on the ideal helix:
		//pose.dump_pdb("TEMP1.pdb"); pose_copy.dump_pdb("TEMP2.pdb"); //DELETE ME
		core::scoring::superimpose_pose( pose, pose_copy, amap, rms_offset_, true );
		//pose.dump_pdb("TEMP3.pdb"); pose_copy.dump_pdb("TEMP4.pdb"); //DELETE ME
		//exit(1); //DELETE ME

	} //End of Scope 1

	//Scope 2 (within this for loop): repeat for all mainchain atoms
	core::Size cur_res(start_index_);
	core::Size cur_atom(0);
	for ( core::Size ia=1; ia<=natoms; ++ia ) {

		++cur_atom; //Increment the current atom.
		if ( cur_atom > pose.residue(cur_res).n_mainchain_atoms() ) { //If we've gone through all of the atoms defining mainchain torsions in the current residue.
			++cur_res; //Increment the current residue
			cur_atom=1; //Start with the first atom of the next residue.
		}
		if ( cur_res==reference_res_index && cur_atom==reference_atom_index ) {
			runtime_assert_string_msg( ia == absolute_reference_atom_index, "Programming error in protocols::helical_bundle::FitSimpleHelix::apply() function.  This shouldn't ever happen.  Please consult a developer." ); //This should be true.
			continue; //Skip the reference atom -- we've already done it.
		}

		std::string const cur_atom_name = pose.residue(cur_res).atom_name( cur_atom );
		if ( TR.visible() ) TR << "Fitting " << cur_atom_name << " atom." << std::endl;

		FitSimpleHelixMultiFunc multfunc(pose, cur_atom_name, cur_res, residues_per_repeat(), start_index_, end_index_, 1, rms_offset_); //Make the multifunc that will be used to fit the current atom.

		Multivec vars(5); //Change the number of DOFs appropriately and initialize.
		vars[1] = (r1_guesses_.size()==0 ? r1_vals_output_[absolute_reference_atom_index] : r1_guesses_[ia]);
		vars[2] = omega1_val_output_;
		vars[3] = z1_val_output_;
		vars[4] = (delta_omega1_guesses_.size()==0 ? (ia > 1 ? delta_omega1_vals_output_[ia-1] : 0.001) : delta_omega1_guesses_[ia] );
		vars[5] = (delta_z1_guesses_.size()==0 ? (ia > 1 ? delta_z1_vals_output_[ia-1] : 0.001) : delta_z1_guesses_[ia] );

		core::optimization::MinimizerOptions minoptions( min_type_, min_tolerance_, false, false, false );
		core::optimization::Minimizer minimizer( multfunc, minoptions );

		//Actually run the minimizer to do the fitting.
		//(The multifunc is set up to trick the minimizer into acting as a nonlinear least-squares fitter for the Crick equations.  Clever, I know.)
		minimizer.run( vars );

		r1_vals_output_[ia] = vars[1];
		delta_omega1_vals_output_[ia] = vars[4];
		delta_z1_vals_output_[ia] = vars[5];

		if ( TR.visible() ) {
			TR << "Atom " << cur_atom_name << " r1 = " << vars[1] << std::endl; //Delete me
			TR << "omega1 = " << vars[2] << std::endl; //Delete me
			TR << "dz1 = " << vars[3] << std::endl; //Delete me
			TR << "Atom " << cur_atom_name << " delta_omega1 = " << vars[4] << std::endl; //Delete me
			TR << "Atom " << cur_atom_name << " delta_z1 = " << vars[5] << std::endl; //Delete me
		}
	}

	//Subtract the arbitrary offset from the pose to shift it back:
	for ( core::Size ir=1, irmax=pose.n_residue(); ir<=irmax; ++ir ) {
		for ( core::Size ia=1, iamax=pose.residue(ir).natoms(); ia<=iamax; ++ia ) {
			pose.set_xyz( core::id::AtomID(ia,ir), pose.xyz( core::id::AtomID(ia,ir) ) - offsetvect );
		}
	}

	if ( TR.visible() ) TR << "Fit complete." << std::endl;

	return;
}

////////////////////////////////////////////////////////////////////////////////


/// @brief Returns the name of this mover ("FitSimpleHelix").
std::string FitSimpleHelix::get_name() const{
	return "FitSimpleHelix";
}

/// @brief Function to retrieve the final values of the helical parameters, post-fit.
/// @details Call this function after the "apply" function, and pass it containers for
/// the data to be retrieved.
void FitSimpleHelix::get_crick_parameters (
	utility::vector1 < core::Real > &r1_out,
	core::Real &omega1_out,
	core::Real &z1_out,
	utility::vector1 < core::Real > &delta_omega1_out,
	utility::vector1 < core::Real > &delta_z1_out
) const {
	r1_out = r1_vals_output_;
	omega1_out = omega1_val_output_;
	z1_out = z1_val_output_;
	delta_omega1_out = delta_omega1_vals_output_;
	delta_z1_out = delta_z1_vals_output_;
	return;
}

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////
/// @brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
///
void
FitSimpleHelix::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const &filters,
	protocols::moves::Movers_map const &movers,
	core::pose::Pose const & pose
) {
	set_initial_guesses(
		tag->getOption<core::Real>("r1_initial", 1.0 ),
		tag->getOption<core::Real>("omega1_initial", 1.0 ),
		tag->getOption<core::Real>("dz1_initial", 1.0 ) );
	
	set_range(
		tag->getOption<core::Size>("start_index", 0 ),
		tag->getOption<core::Size>("end_index", 0 ) );
	
	set_min_type( tag->getOption<std::string>("min_type", "lbfgs_armijo_nonmonotone" ) );
	set_min_tolerance( tag->getOption<core::Real>("min_tolerance", 0.00000001 ) );
	set_reference_atom( tag->getOption<std::string>("reference_atom", "CA" ) );
	set_reference_residue( tag->getOption<core::Size>("reference_residue", 1 ) );
	set_residues_per_repeat( tag->getOption<core::Size>("residues_per_repeat", 1 ) );
}


////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////


} //namespace helical_bundle
} //namespace protocols
