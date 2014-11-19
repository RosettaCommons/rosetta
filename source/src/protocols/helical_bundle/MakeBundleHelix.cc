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

//static numeric::random::RandomGenerator RG(741701);  // <- Magic number, do not change it!

namespace protocols {
namespace helical_bundle {

static basic::Tracer TR("protocols.helical_bundle.MakeBundleHelix");

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

///
///@brief Creator for MakeBundleHelix mover.
MakeBundleHelix::MakeBundleHelix():
		Mover("MakeBundleHelix"),
		reset_pose_(true),
		bundle_parameters_( BundleParametersOP( new BundleParameters ) ),
		helix_length_(10),
		residue_name_("ALA"),
		tail_residue_name_(""),
		last_apply_failed_(false)
{
	set_minor_helix_params_from_file("alpha_helix"); //By default, set the minor helix parameters to those of an alpha helix (read in from the database).
}

///
///@brief Destructor for MakeBundleHelix mover.
MakeBundleHelix::~MakeBundleHelix() {}

///
///@brief Clone operator to create a pointer to a fresh MakeBundleHelix object that copies this one.
protocols::moves::MoverOP MakeBundleHelix::clone() const {
	return protocols::moves::MoverOP( new MakeBundleHelix( *this ) );
}

///
///@brief Fresh_instance operator to create a pointer to a fresh MakeBundleHelix object that does NOT copy this one.
protocols::moves::MoverOP MakeBundleHelix::fresh_instance() const {
	return protocols::moves::MoverOP( new MakeBundleHelix );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////

///
/// @brief Actually apply the mover to the pose.
void MakeBundleHelix::apply (core::pose::Pose & pose)
{
	if(TR.visible()) TR << "Building a helix in a helical bundle using the Crick equations." << std::endl;

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
	stubmover.add_residue ("Append", residue_name(), 0, (tail_residues_exist ? false : true), "", helix_length(), 0, "");
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
	generate_atom_positions(atom_positions, helixpose, helix_start, helix_end, failed );
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

///
///@brief Returns the name of this mover ("MakeBundleHelix").
std::string MakeBundleHelix::get_name() const{
	return "MakeBundleHelix";
}

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////
///@brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
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

/// @brief Generate the x,y,z coordinates of the mainchain atoms using the Crick equations.
/// @details Coordinates will be returned as a vector of vectors of xyzVectors.  The outer
/// index will refer to residue number, and the inner index will refer to atom number.
/// Returns failed=true if coordinates could not be generated, false otherwise.
void MakeBundleHelix::generate_atom_positions(
	utility::vector1 < utility::vector1 < numeric::xyzVector< core::Real > > > &outvector,
	core::pose::Pose const &helixpose,
	core::Size const helix_start,
	core::Size const helix_end,
	bool &failed
) const {

	outvector.clear();
	failed=false;

	core::Real t = -1.0 * static_cast<core::Real>(helix_length() + 2) / 2.0 + delta_t();

	for(core::Size ir2=helix_start-1; ir2<=helix_end+1; ++ir2) { //Loop through all residues in the helix, padded by one on either side
		core::Size ir = ir2;
		//Repeat start and end residues to pad by one:
		if(ir2==helix_start-1) ir=helix_start;
		if(ir2==helix_end+1) ir=helix_end;
		utility::vector1 < numeric::xyzVector <core::Real> > innervector;
		for(core::Size ia=1, iamax=helixpose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia) { //Loop through all mainchain atoms in the current helix residue
			innervector.push_back( numeric::crick_equations::XYZ_BUNDLE(
				t, r0(), omega0(),
				delta_omega0() + ( invert_helix() ? numeric::constants::d::pi : 0 ),
				r1(ia), omega1()-omega0(), z1(), delta_omega1(ia)+delta_omega1_all(), delta_z1(ia), failed )
			);
			if(invert_helix()) {
				innervector[ia].x( -1.0*innervector[ia].x() );
				//innervector[ia].y( -1.0*innervector[ia].y() );
				innervector[ia].z( -1.0*innervector[ia].z() );
			}
			if(failed) break;
		}
		if(failed) break;
		outvector.push_back(innervector);
		t+=1.0;
	}

	if(failed) outvector.clear();

	return;
}

/// @brief Place the helix mainchain atoms based on the Crick equations.
///
void MakeBundleHelix::place_atom_positions(
	core::pose::Pose &pose,
	utility::vector1 < utility::vector1 < numeric::xyzVector < core::Real >  > > const &atom_positions,
	core::Size const helix_start,
	core::Size const helix_end
) const {

	//Index in the outer vector for the current residue.
	core::Size index=2;

	for(core::Size ir=helix_start; ir<=helix_end; ++ir) {
		for(core::Size ia=1, iamax=atom_positions[index].size(); ia<=iamax; ++ia) {
			pose.set_xyz( core::id::AtomID( ia, ir ), atom_positions[index][ia] );
		}
		++index;
	}

	pose.update_residue_neighbors();

	return;
}

/// @brief Copy backbone bond length values from one pose, where helix mainchain atom coordinates have been
/// set with the Crick equations, to another with ideal geometry.
void MakeBundleHelix::copy_helix_bondlengths(
	core::pose::Pose &pose,
	core::pose::Pose const &ref_pose,
	core::Size const helix_start,
	core::Size const helix_end
) const {
	for(core::Size ir=helix_start; ir<=helix_end; ++ir) {
		for(core::Size ia=1, iamax=ref_pose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia) {
			if(ia==1 && ir==helix_start) continue; //Skip the first atom.
			core::id::AtomID const thisatom( ia, ir );
			core::id::AtomID const prevatom( (ia==1 ? iamax : ia - 1), (ia==1 ? ir-1 : ir) ); //The previous atom is ia-1 in this residue, unless this atom is the first atom, in which case the previous atom is iamax in the previous residue.
			pose.conformation().set_bond_length( thisatom, prevatom, ref_pose.xyz(thisatom).distance( ref_pose.xyz(prevatom) ) );
		}
	}

	//TODO properly handle the first and last residues using the extra residue xyz coordinates that were generated!

	pose.update_residue_neighbors();

	return;
}

/// @brief Copy backbone bond angle values from one pose, where helix mainchain atom coordinates have been
/// set with the Crick equations, to another with ideal geometry.
void MakeBundleHelix::copy_helix_bondangles(
	core::pose::Pose &pose,
	core::pose::Pose const &ref_pose,
	core::Size const helix_start,
	core::Size const helix_end
) const {
	for(core::Size ir=helix_start; ir<=helix_end; ++ir) {
		for(core::Size ia=1, iamax=ref_pose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia) {
			if(ia==1 && ir==helix_start) continue; //Skip the first atom.
			if(ia==iamax && ir==helix_end) continue; //Skip the last atom.
			core::id::AtomID const thisatom( ia, ir );
			core::id::AtomID const prevatom( (ia==1 ? iamax : ia - 1), (ia==1 ? ir-1 : ir) ); //The previous atom is ia-1 in this residue, unless this atom is the first atom in this residue, in which case the previous atom is iamax in the previous residue.
			core::id::AtomID const nextatom( (ia==iamax ? 1 : ia + 1), (ia==iamax ? ir+1 : ir) ); //The next atom is ia+1 in this residue, unless this atom is the last atom in this residue, in which case the next atom is the first atom in the next residue.
			pose.conformation().set_bond_angle( prevatom, thisatom, nextatom, numeric::angle_radians<core::Real>( ref_pose.xyz(prevatom), ref_pose.xyz(thisatom), ref_pose.xyz(nextatom) ) );
		}
	}

	//TODO properly handle the first and last residues using the extra residue xyz coordinates that were generated!

	pose.update_residue_neighbors();

	return;
}

/// @brief Copy backbone dihedral values from one pose, where helix mainchain atom coordinates have been
/// set with the Crick equations, to another with ideal geometry.
void MakeBundleHelix::copy_helix_dihedrals(
	core::pose::Pose &pose,
	core::pose::Pose const &ref_pose,
	core::Size const helix_start,
	core::Size const helix_end
) const {
	for(core::Size ir=helix_start; ir<=helix_end; ++ir) {
		for(core::Size itors=1, itorsmax=ref_pose.residue(ir).mainchain_torsions().size(); itors<=itorsmax; ++itors) {
			pose.conformation().set_torsion( core::id::TorsionID( ir, core::id::BB, itors ), ref_pose.conformation().torsion( core::id::TorsionID( ir, core::id::BB, itors ) ) );
		}
	}

	//TODO properly handle the first and last residues using the extra residue xyz coordinates that were generated!

	pose.update_residue_neighbors();

	return;
}

/// @brief Align mainchain atoms of pose to ref_pose mainchain atoms.
///
void MakeBundleHelix::align_mainchain_atoms(
	core::pose::Pose &pose,
	core::pose::Pose const &ref_pose,
	core::Size const helix_start,
	core::Size const helix_end
) const {
		core::id::AtomID_Map< core::id::AtomID > amap;
		core::pose::initialize_atomid_map(amap, pose, core::id::BOGUS_ATOM_ID);

		for(core::Size ir=helix_start; ir<=helix_end; ++ir) {
			for(core::Size ia=1, iamax=ref_pose.residue(ir).n_mainchain_atoms(); ia<=iamax; ++ia) {
				amap[core::id::AtomID(ia,ir)] = core::id::AtomID(ia,ir);
			}
		}

		core::scoring::superimpose_pose( pose, ref_pose, amap );

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
