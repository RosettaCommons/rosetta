// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/PerturbBundleHelix.cc
/// @brief  Reads the Crick parameters for a piece of a pose from the input pose and sets the
/// mainchain torsions accordingly.  The parameters are presumed to have been perturbed by
/// another mover.  This mover is intended to be called by the PerturbBundle mover, which handles
/// the perturbation of the Crick parameters.
/// @details The bundle is centred on the origin, with the outer helix axis pointing along the
/// global z-axis.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit Headers
#include <protocols/helical_bundle/PerturbBundleHelix.hh>
#include <protocols/helical_bundle/PerturbBundleHelixCreator.hh>
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

//static numeric::random::RandomGenerator RG(741701);  // <- Magic number, do not change it!

namespace protocols {
namespace helical_bundle {

static THREAD_LOCAL basic::Tracer TR("protocols.helical_bundle.PerturbBundleHelix");

std::string
PerturbBundleHelixCreator::keyname() const
{
	return PerturbBundleHelixCreator::mover_name();
}

protocols::moves::MoverOP
PerturbBundleHelixCreator::create_mover() const {
	return protocols::moves::MoverOP( new PerturbBundleHelix );
}

std::string
PerturbBundleHelixCreator::mover_name()
{
	return "PerturbBundleHelix";
}


/// @brief Creator for PerturbBundleHelix mover.
PerturbBundleHelix::PerturbBundleHelix():
	Mover("PerturbBundleHelix"),
	parameters_set_index_(0),
	parameters_index_(0),
	last_apply_failed_(false)
{
}


/// @brief Destructor for PerturbBundleHelix mover.
PerturbBundleHelix::~PerturbBundleHelix() {}


/// @brief Clone operator to create a pointer to a fresh PerturbBundleHelix object that copies this one.
protocols::moves::MoverOP PerturbBundleHelix::clone() const {
	return protocols::moves::MoverOP( new PerturbBundleHelix( *this ) );
}


/// @brief Fresh_instance operator to create a pointer to a fresh PerturbBundleHelix object that does NOT copy this one.
protocols::moves::MoverOP PerturbBundleHelix::fresh_instance() const {
	return protocols::moves::MoverOP( new PerturbBundleHelix );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void PerturbBundleHelix::apply (core::pose::Pose & pose)
{
	if ( TR.visible() ) TR << "Perturbing a helix in a helical bundle using the Crick equations." << std::endl;

	runtime_assert_string_msg( parameters_set_index()>0  && parameters_set_index()<=pose.conformation().n_parameters_sets(),
		"In protocols::helical_bundle::PerturbBundleHelix::apply() : The index of the ParametersSet object is not an index that exists in the pose." );
	runtime_assert_string_msg( parameters_index()>0  && parameters_index()<=pose.conformation().parameters_set( parameters_set_index() )->n_parameters(),
		"In protocols::helical_bundle::PerturbBundleHelix::apply() : The index of the Parameters object is not an index that exists in the ParametersSet object." );

	bool failed(false);

	BundleParametersSetOP paramset = utility::pointer::dynamic_pointer_cast<BundleParametersSet>( pose.conformation().parameters_set( parameters_set_index() ) );
	runtime_assert_string_msg( paramset, "In protocols::helical_bundle::PerturbBundleHelix::apply() : The ParametersSet object is not a BundleParametersSet." );
	BundleParametersOP params = utility::pointer::dynamic_pointer_cast<BundleParameters>( paramset->parameters( parameters_index() ) );
	runtime_assert_string_msg( params, "In protocols::helical_bundle::PerturbBundleHelix::apply() : The Parameters object is not a BundleParameters object." );

	utility::vector1 < utility::vector1 < numeric::xyzVector< core::Real > > > atom_positions;

	//Generate the atom positions (function defined in util.cc).
	generate_atom_positions(
		atom_positions, //output
		pose, //for reference only
		params->first_residue()->seqpos(),
		params->last_residue()->seqpos(),
		params->r0(),
		params->omega0(),
		params->delta_omega0(),
		params->delta_t(),
		params->z1_offset(),
		params->z0_offset(),
		params->invert_helix(),
		params->r1_vect(),
		params->omega1(),
		params->z1(),
		params->delta_omega1_vect(),
		params->delta_omega1_all(),
		params->delta_z1_vect(),
		params->residues_per_repeat(),
		params->atoms_per_residue(),
		params->repeating_unit_offset(),
		failed
	);

	set_last_apply_failed(failed);
	if ( failed ) {
		if ( TR.visible() ) TR << "Mover failed.  The Crick parameters do not generate a sensible helix.  Returning input pose." << std::endl;
		return; //At this point, the input pose has not been modified.
	}

	// Make a temporary pose copy and place the atoms using the Crick equations.
	if ( TR.Debug.visible() ) TR.Debug << "Placing atoms." << std::endl;
	core::pose::Pose pose_copy(pose);
	place_atom_positions(pose_copy, atom_positions, params->first_residue()->seqpos(), params->last_residue()->seqpos());

	//Copy the bond length values from the pose in which we set the x,y,z coordinates of mainchain atoms based on the Crick equations
	//to the ideal pose.
	if ( params->allow_bondlengths() ) {
		if ( TR.Debug.visible() ) TR.Debug << "Copying bond length values." << std::endl;
		copy_helix_bondlengths(pose, pose_copy, params->first_residue()->seqpos(), params->last_residue()->seqpos());
	}

	//Copy the bond angle values from the pose in which we set the x,y,z coordinates of mainchain atoms based on the Crick equations
	//to the ideal pose.
	if ( params->allow_bondangles() ) {
		if ( TR.Debug.visible() ) TR.Debug << "Copying bond angle values." << std::endl;
		copy_helix_bondangles(pose, pose_copy, params->first_residue()->seqpos(), params->last_residue()->seqpos());
	}

	//Copy the dihedral values from the pose in which we set the x,y,z coordinates of mainchain atoms based on the Crick equations
	//to the ideal pose.
	if ( params->allow_dihedrals() ) {
		if ( TR.Debug.visible() ) TR.Debug << "Copying dihedral values." << std::endl;
		copy_helix_dihedrals(pose, pose_copy, params->first_residue()->seqpos(), params->last_residue()->seqpos());
	}

	//TODO align perturbed helix to ideal helix position.
	align_mainchain_atoms_of_residue_range( pose, pose_copy, params->first_residue()->seqpos(), params->last_residue()->seqpos());

	if ( TR.Debug.visible() ) TR.Debug << "Finished apply function." << std::endl;

	return;
}

////////////////////////////////////////////////////////////////////////////////


/// @brief Returns the name of this mover ("PerturbBundleHelix").
std::string PerturbBundleHelix::get_name() const{
	return "PerturbBundleHelix";
}

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////
/// @brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
///
/*void
PerturbBundleHelix::parse_my_tag(
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
