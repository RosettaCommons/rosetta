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

static basic::Tracer TR("protocols.helical_bundle.FitSimpleHelix");

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

///
///@brief Creator for FitSimpleHelix mover.
FitSimpleHelix::FitSimpleHelix():
		//utility::pointer::ReferenceCount(),
		Mover("FitSimpleHelix"),
		r1_initial_(1.0),
		omega1_initial_(1.0),
		dz1_initial_(1.0),
		start_index_(0),
		end_index_(0),
		min_type_("dfpmin"),
		min_tolerance_(0.00000001),
		reference_atom_("CA"),
		r1_vals_output_(),
		omega1_val_output_(0.0),
		z1_val_output_(0.0),
		delta_omega1_vals_output_(),
		delta_z1_vals_output_()
{}

///
///@brief Destructor for FitSimpleHelix mover.
FitSimpleHelix::~FitSimpleHelix() {}

///
///@brief Clone operator to create a pointer to a fresh FitSimpleHelix object that copies this one.
protocols::moves::MoverOP FitSimpleHelix::clone() const {
	return protocols::moves::MoverOP( new FitSimpleHelix( *this ) );
}

///
///@brief Fresh_instance operator to create a pointer to a fresh FitSimpleHelix object that does NOT copy this one.
protocols::moves::MoverOP FitSimpleHelix::fresh_instance() const {
	return protocols::moves::MoverOP( new FitSimpleHelix );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////

///
/// @brief Actually apply the mover to the pose.
/// @details At the end of this function, two things will have happened:
/// --The pose will be realigned to the ideal helix.
/// --The fitted helical parameters will have been stored in internal variables for the mover, and can be retrieved
/// using a to-be-written function.  TODO
void FitSimpleHelix::apply (core::pose::Pose & pose)
{
	runtime_assert_string_msg( start_index_ > 1 && start_index_ < pose.n_residue() && end_index_ > 1 && end_index_ < pose.n_residue(), "In FitSimpleHelix::apply(): The start and end residues must be within the pose's residue range, and cannot be the end residues." );
	runtime_assert_string_msg( start_index_ + 1 < end_index_, "In FitSimpleHelix::apply(): The end residue must be at least two residues past the start residue." );

	r1_vals_output_.resize( pose.residue(start_index_).n_mainchain_atoms(), 0.0 );
	delta_omega1_vals_output_.resize( pose.residue(start_index_).n_mainchain_atoms(), 0.0 );
	delta_z1_vals_output_.resize( pose.residue(start_index_).n_mainchain_atoms(), 0.0 );

	{ //Scope 1: fiting the reference atom.

		core::Size const reference_atom_index = pose.residue(start_index_).atom_index( reference_atom_ );

		runtime_assert_string_msg( reference_atom_index <= pose.residue(start_index_).n_mainchain_atoms(), "In FitSimpleHelix::apply(): The reference atom is not a mainchain atom." );

		if(TR.visible()) TR << "Fitting reference atom " << reference_atom_ << "." << std::endl;

		FitSimpleHelixMultiFunc multfunc(pose, reference_atom_, start_index_, end_index_, 0); //Make the multifunc that will be used to fit the reference atom.

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

		r1_vals_output_[reference_atom_index] = vars[1];
		omega1_val_output_ = vars[2];
		z1_val_output_ = vars[3];
		delta_omega1_vals_output_[reference_atom_index] = vars[4];
		delta_z1_vals_output_[reference_atom_index] = vars[5];

		if(TR.visible()) {
			TR << "Reference atom r1 = " << vars[1] << std::endl; //Delete me
			TR << "Reference atom omega1 = " << vars[2] << std::endl; //Delete me
			TR << "Reference atom dz1 = " << vars[3] << std::endl; //Delete me
			TR << "Reference atom delta_omega1 = " << vars[4] << std::endl; //Delete me
			TR << "Reference atom delta_z1 = " << vars[5] << std::endl; //Delete me
		}

		if(TR.visible()) TR << "Aligning pose to ideal helix." << std::endl;

		core::pose::Pose pose_copy = pose; //Make a copy of the pose.

		core::id::AtomID_Map< core::id::AtomID > amap;
		core::pose::initialize_atomid_map(amap, pose, core::id::BOGUS_ATOM_ID);

		core::Real t = -1.0*static_cast<core::Real>(end_index_ - start_index_ + 1)/2.0;

		for(core::Size ir=start_index_; ir<=end_index_; ++ir) {
			pose_copy.set_xyz( core::id::NamedAtomID(reference_atom_,ir), numeric::crick_equations::xyz( vars[1], vars[2] , t, vars[3], 0.0, 0.0 ) );
			amap[core::id::AtomID(pose.residue(ir).atom_index(reference_atom_),ir)] = core::id::AtomID(pose_copy.residue(ir).atom_index(reference_atom_),ir);
			t+=1.0;
		}

		//Superimpose the original helix on the ideal helix:
		//pose.dump_pdb("TEMP1.pdb"); pose_copy.dump_pdb("TEMP2.pdb"); //DELETE ME
		core::scoring::superimpose_pose( pose, pose_copy, amap );
		//pose.dump_pdb("TEMP3.pdb"); pose_copy.dump_pdb("TEMP4.pdb"); //DELETE ME
		//exit(1); //DELETE ME

	} //End of Scope 1

	//Scope 2 (within this for loop): repeat for all mainchain atoms
	for(core::Size ia=1, iamax=pose.residue(start_index_).n_mainchain_atoms(); ia<=iamax; ++ia) {
		if(pose.residue(start_index_).atom_index(reference_atom_) == ia) continue; //Skip the reference atom -- we've already done it.

		std::string const cur_atom_name = pose.residue(start_index_).atom_name(ia);

		if(TR.visible()) TR << "Fitting " << cur_atom_name << " atom." << std::endl;

		FitSimpleHelixMultiFunc multfunc(pose, cur_atom_name, start_index_, end_index_, 1); //Make the multifunc that will be used to fit the current atom.

		Multivec vars(5); //Change the number of DOFs appropriately and initialize.
		vars[1] = r1_vals_output_[pose.residue(start_index_).atom_index( reference_atom_ )];
		vars[2] = omega1_val_output_;
		vars[3] = z1_val_output_;
		vars[4] = 0.001;
		vars[5] = 0.001;

		core::optimization::MinimizerOptions minoptions( min_type_, min_tolerance_, false, false, false );
		core::optimization::Minimizer minimizer( multfunc, minoptions );

		//Actually run the minimizer to do the fitting.
		//(The multifunc is set up to trick the minimizer into acting as a nonlinear least-squares fitter for the Crick equations.  Clever, I know.)
		minimizer.run( vars );

		r1_vals_output_[ia] = vars[1];
		delta_omega1_vals_output_[ia] = vars[4];
		delta_z1_vals_output_[ia] = vars[5];

		if(TR.visible()) {
			TR << "Atom " << cur_atom_name << " r1 = " << vars[1] << std::endl; //Delete me
			TR << "omega1 = " << vars[2] << std::endl; //Delete me
			TR << "dz1 = " << vars[3] << std::endl; //Delete me
			TR << "Atom " << cur_atom_name << " delta_omega1 = " << vars[4] << std::endl; //Delete me
			TR << "Atom " << cur_atom_name << " delta_z1 = " << vars[5] << std::endl; //Delete me
		}

	}

	if(TR.visible()) TR << "Fit complete." << std::endl;

	return;
}

////////////////////////////////////////////////////////////////////////////////

///
///@brief Returns the name of this mover ("FitSimpleHelix").
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
///@brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
///
/*void
FitSimpleHelix::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const &movers,
		core::pose::Pose const & pose
) {


	return;
}*/



////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////


} //namespace helical_bundle
} //namespace protocols
