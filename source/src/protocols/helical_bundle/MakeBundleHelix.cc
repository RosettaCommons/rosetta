// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/MakeBundleHelix.cc
/// @brief  Builds a single helix as part of a helical bundle.
/// @details The bundle is centred on the origin, with the outer helix axis pointing along the
/// global z-axis.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit Headers
#include <protocols/helical_bundle/MakeBundleHelix.hh>
#include <protocols/helical_bundle/MakeBundleHelixCreator.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <numeric/crick_equations/BundleParams.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <core/id/TorsionID.hh>
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



namespace protocols {
namespace helical_bundle {

static thread_local basic::Tracer TR("protocols.helical_bundle.MakeBundleHelix");

std::string
MakeBundleHelixCreator::keyname() const
{
	return MakeBundleHelixCreator::mover_name();
}

protocols::moves::MoverOP
MakeBundleHelixCreator::create_mover() const {
	return protocols::moves::MoverOP( new MakeBundleHelix );
}

std::string
MakeBundleHelixCreator::mover_name()
{
	return "MakeBundleHelix";
}


/// @brief Constructor for MakeBundleHelix mover.
MakeBundleHelix::MakeBundleHelix():
		Mover("MakeBundleHelix"),
		reset_pose_(true),
		bundle_parameters_( BundleParametersOP( new BundleParameters ) ),
		helix_length_(10),
		residue_name_(),
		tail_residue_name_(""),
		last_apply_failed_(false)
{
	set_minor_helix_params_from_file("alpha_helix"); //By default, set the minor helix parameters to those of an alpha helix (read in from the database).
	residue_name_.push_back("ALA");
}


/// @brief Copy constructor for MakeBundleHelix mover.
MakeBundleHelix::MakeBundleHelix( MakeBundleHelix const &src ):
		protocols::moves::Mover( src ),
		reset_pose_(src.reset_pose_),
		bundle_parameters_(  utility::pointer::dynamic_pointer_cast< BundleParameters >( src.bundle_parameters()->clone() ) ),
		helix_length_(src.helix_length_),
		residue_name_(src.residue_name_),
		tail_residue_name_(src.tail_residue_name_),
		last_apply_failed_(src.last_apply_failed_)
{
}


/// @brief Destructor for MakeBundleHelix mover.
MakeBundleHelix::~MakeBundleHelix() {}


/// @brief Clone operator to create a pointer to a fresh MakeBundleHelix object that copies this one.
protocols::moves::MoverOP MakeBundleHelix::clone() const {
	return protocols::moves::MoverOP( new MakeBundleHelix ( *this ) );
}


/// @brief Fresh_instance operator to create a pointer to a fresh MakeBundleHelix object that does NOT copy this one.
protocols::moves::MoverOP MakeBundleHelix::fresh_instance() const {
	return protocols::moves::MoverOP( new MakeBundleHelix );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void MakeBundleHelix::apply (core::pose::Pose & pose)
{
	if(TR.visible()) TR << "Building a helix in a helical bundle using the Crick equations." << std::endl;

	//Initial checks:
	runtime_assert_string_msg(
		residues_per_repeat() == residue_name_.size(),
		"In protocols::helical_bundle::MakeBundleHelix::apply(): The number of residues per repeat does not match the size of the list of residue types."
	);

	//Should tail residues be added:
	bool const tail_residues_exist = ( tail_residue_name()!="" );

	//Create the pose object:
	core::pose::Pose helixpose;

	//Build the pose:
	if(TR.Debug.visible()) TR.Debug << "Doing initial build." << std::endl;
	protocols::cyclic_peptide::PeptideStubMover stubmover;
	stubmover.set_reset_mode( true );
	stubmover.reset_mover_data();
	if(tail_residues_exist) stubmover.add_residue ("Append", tail_residue_name(), 0, true, "", 1, 0, "");
	core::Size repeat_index( repeating_unit_offset() ); //Index in the repeating unit making up the minor helix
	for(core::Size i=1, imax=helix_length(); i<=imax; ++i) {
		++repeat_index;
		if(repeat_index > residues_per_repeat()) repeat_index=1;
		stubmover.add_residue ("Append", residue_name(repeat_index), 0, (i==1 && !tail_residues_exist ? true : false), "", 1, 0, "");
	}
	//stubmover.add_residue ("Append", residue_name(), 0, (tail_residues_exist ? false : true), "", helix_length(), 0, "");
	if(tail_residues_exist) stubmover.add_residue ("Append", tail_residue_name(), 0, false, "", 1, 0, "");
	stubmover.apply(helixpose);

	//If there are tail residues, set their backbone dihedral values to something reasonable:
	if(tail_residues_exist) {
		if(TR.Debug.visible()) TR.Debug << "Setting tail dihedrals." << std::endl;
		set_tail_dihedrals(helixpose);
	}

	//Variables for the start residue and end residue of the helix in helixpose:
	core::Size const helix_start = (tail_residues_exist ? 2 : 1);
	core::Size const helix_end = helix_length() + (tail_residues_exist ? 1 : 0 );

	//Generate a vector of vectors of x,y,z coordinates of mainchain atoms.
	//Note that this vector extends one extra residue in each direction.
	utility::vector1 < utility::vector1 < numeric::xyzVector< core::Real > > > atom_positions;
	if(TR.Debug.visible()) TR.Debug << "Generating atom positions." << std::endl;
	bool failed=false;
	generate_atom_positions(atom_positions, helixpose, helix_start, helix_end, r0(),
		omega0(), delta_omega0(), delta_t(), z1_offset(), z0_offset(), invert_helix(), r1_vect(), omega1(), z1(),
		delta_omega1_vect(), delta_omega1_all(), delta_z1_vect(), residues_per_repeat(), atoms_per_residue(), repeating_unit_offset(),
		failed );

	set_last_apply_failed(failed);
	if(failed) {
		if(TR.visible()) TR << "Mover failed.  The Crick parameters do not generate a sensible helix.  Returning input pose." << std::endl;
		return; //At this point, the input pose has not been modified.
	}

	// Make a temporary pose copy and place the atoms using the Crick equations.
	if(TR.Debug.visible()) TR.Debug << "Placing atoms." << std::endl;
	core::pose::Pose helixpose_copy = helixpose;
	place_atom_positions(helixpose_copy, atom_positions, helix_start, helix_end);

	//Copy the bond length values from the pose in which we set the x,y,z coordinates of mainchain atoms based on the Crick equations
	//to the ideal pose.
	if(allow_bondlengths()) {
		if(TR.Debug.visible()) TR.Debug << "Copying bond length values." << std::endl;
		copy_helix_bondlengths(helixpose, helixpose_copy, helix_start, helix_end);
	}

	//Copy the bond angle values from the pose in which we set the x,y,z coordinates of mainchain atoms based on the Crick equations
	//to the ideal pose.
	if(allow_bondangles()) {
		if(TR.Debug.visible()) TR.Debug << "Copying bond angle values." << std::endl;
		copy_helix_bondangles(helixpose, helixpose_copy, helix_start, helix_end);
	}

	//Copy the dihedral values from the pose in which we set the x,y,z coordinates of mainchain atoms based on the Crick equations
	//to the ideal pose.
	if(allow_dihedrals()) {
		if(TR.Debug.visible()) TR.Debug << "Copying dihedral values." << std::endl;
		copy_helix_dihedrals(helixpose, helixpose_copy, helix_start, helix_end);
	}

	//Align the ideal pose to the pose in which we set the x,y,z coordinates of mainchain atoms.
	if(TR.Debug.visible()) TR.Debug << "Aligning to ideal helix." << std::endl;
	align_mainchain_atoms(helixpose, helixpose_copy, helix_start, helix_end);

	//Either reset the pose and replace it with the helix pose, or append the helix pose to the current pose, depending on the reset mode:
	if(reset_pose()) {
		if(TR.Debug.visible()) TR.Debug << "Clearing pose and adding helix to pose." << std::endl;
		pose.clear();
		pose=helixpose;
	} else {
		if(TR.Debug.visible()) TR.Debug << "Appending helix to pose." << std::endl;
		pose.append_pose_by_jump(helixpose, 1);
	}

	if(TR.Debug.visible()) TR.Debug << "Adding Crick parameter data to Conformation oject." << std::endl;
	add_parameter_info_to_pose( pose );

	if(TR.Debug.visible()) TR.Debug << "Finished apply function." << std::endl;

	return;
}

////////////////////////////////////////////////////////////////////////////////


/// @brief Returns the name of this mover ("MakeBundleHelix").
std::string MakeBundleHelix::get_name() const{
	return "MakeBundleHelix";
}

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////
/// @brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
///
/*void
MakeBundleHelix::parse_my_tag(
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

/// @brief If there are tail residues, set their backbone dihedral angles
/// to something reasonable (a helical conformation).
void MakeBundleHelix::set_tail_dihedrals( core::pose::Pose &helixpose ) const
{
		if(helixpose.residue(1).type().is_alpha_aa()) {
			helixpose.set_psi( 1, -47.0 );
			helixpose.set_omega(1, 180.0);
			helixpose.set_phi(helixpose.n_residue(), -57.8);
		} else if (helixpose.residue(1).type().is_beta_aa() ) {
			helixpose.set_theta( 1, 55.6 );
			helixpose.set_psi( 1, -126.4 );
			helixpose.set_omega(1, 180.0);
			helixpose.set_phi(helixpose.n_residue(), -135.5);
			helixpose.set_theta( helixpose.n_residue(), 55.6 );
		}
		if(helixpose.residue(helixpose.n_residue()-1).type().is_alpha_aa() || helixpose.residue(helixpose.n_residue()-1).type().is_beta_aa()) {
			helixpose.set_omega(helixpose.n_residue()-1, 180.0);
		}
		return;
}

/// @brief Add Crick parameter information to the Conformation object within the pose.
/// @details This function updates the bundle_parameters_ object's links to residues within the pose,
/// and then copies the owning pointer into the pose's Conformation object.  The bundle_parameters_
/// object will be added to a new ParameterSet object in the Conformation object.
void MakeBundleHelix::add_parameter_info_to_pose( core::pose::Pose &pose )
{
	bool const has_tails( tail_residue_name()!="" );
	core::Size const first_res( pose.n_residue() - helix_length() - ( has_tails ? 2 : 0 ) + 1 );
	core::Size const last_res( pose.n_residue() - ( has_tails ? 1 : 0 ));

	BundleParametersOP output_parameters( utility::pointer::dynamic_pointer_cast<BundleParameters>( bundle_parameters_->clone() ) );

	output_parameters->reset_residue_list();

	for(core::Size ir=first_res; ir<=last_res; ++ir) { //Loop through all of the helix residues.
		output_parameters->add_residue( pose.conformation().residue_op( ir ) ); //Add owning pointers to the residue objects.
	}

	//Create a new ParametersSet in the Conformation object:
	pose.conformation().create_new_parameters_set();

	//Add the output_parameters ParametersSet to the new ParametersSet
	pose.conformation().parameters_set( pose.conformation().n_parameters_sets() )->add_parameters( utility::pointer::dynamic_pointer_cast<Parameters>(output_parameters) );

	return;
}

} //namespace helical_bundle
} //namespace protocols
