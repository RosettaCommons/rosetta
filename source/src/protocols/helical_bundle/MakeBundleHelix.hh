// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/MakeBundleHelix.hh
/// @brief  Headers for MakeBundleHelix.cc.  Builds a single helix as part of a helical bundle.
/// @details The bundle is centred on the origin, with the outer helix axis pointing along the
/// global z-axis.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_helical_bundle_MakeBundleHelix_hh
#define INCLUDED_protocols_helical_bundle_MakeBundleHelix_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/helical_bundle/MakeBundleHelix.fwd.hh>
#include <protocols/helical_bundle/util.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.hh>
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/ParametersSet.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/ContingentFilter.fwd.hh>
#include <protocols/filters/ContingentFilter.hh>

// Project Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/constants.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>

#include <set>

#include <core/grid/CartGrid.fwd.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace helical_bundle {

class MakeBundleHelix : public protocols::moves::Mover
{
public:
		//Typedefs:
		typedef core::conformation::parametric::Parameters Parameters;
		typedef core::conformation::parametric::ParametersOP ParametersOP;
		typedef core::conformation::parametric::ParametersSet ParametersSet;
		typedef core::conformation::parametric::ParametersSetOP ParametersSetOP;

		typedef protocols::helical_bundle::parameters::BundleParameters BundleParameters;
		typedef protocols::helical_bundle::parameters::BundleParametersOP BundleParametersOP;
		typedef protocols::helical_bundle::parameters::BundleParametersCOP BundleParametersCOP;
		typedef protocols::helical_bundle::parameters::BundleParametersSet BundleParametersSet;
		typedef protocols::helical_bundle::parameters::BundleParametersSetOP BundleParametersSetOP;
		typedef protocols::helical_bundle::parameters::BundleParametersSetCOP BundleParametersSetCOP;

public:
	MakeBundleHelix();
	MakeBundleHelix( MakeBundleHelix const &src );
	virtual ~MakeBundleHelix();

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;


	/// @brief Actually apply the mover to the pose.
	virtual void apply(core::pose::Pose & pose);

	virtual std::string get_name() const;

	/*virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const &
	);*/

	/// @brief Set whether the input pose should be reset prior to building a helix.
	///
	void set_reset_pose( bool const reset_in) { reset_pose_ = reset_in; }

	/// @brief Returns a bool indicating whether the input pose will be reset prior to
	/// building a helix.
	bool reset_pose() const { return reset_pose_; }

	/// @brief Set whether the helix direction should be inverted.
	///
	void set_invert_helix(bool const invert_in ) { bundle_parameters_->set_invert_helix(invert_in); return; }

	/// @brief Return whether the helix direction will be inverted.
	///
	bool invert_helix() const { return bundle_parameters_->invert_helix(); }

	/// @brief Set the length of the helix, in residues.
	///
	void set_helix_length( core::Size const helix_length_in ) {
		runtime_assert_string_msg(helix_length_in >=2, "In protocols::helical_bundle::MakeBundleHelix::set_helix_length(): A helix must contain at least two residues.");
		helix_length_ = helix_length_in;
	}

	/// @brief Returns a bool indicating whether the input pose will be reset prior to
	/// building a helix.
	core::Size helix_length() const { return helix_length_; }

	/// @brief Set the residue type (full name, not 3-letter code) that will make up the helix.
	///
	void set_residue_name(std::string const &name) { residue_name_=name; }

	/// @brief Get the name of the residue type (full name, not 3-letter code) that will make up the helix.
	///
	std::string residue_name() const { return residue_name_; }

	/// @brief Set the residue type (full name, not 3-letter code) that will cap the helix.
	/// @details Set this to "" for no cap.
	void set_tail_residue_name(std::string const &name) { tail_residue_name_=name; }

	/// @brief Get the name of the residue type (full name, not 3-letter code) that will cap the helix.
	/// @details Returns "" for no cap.
	std::string tail_residue_name() const { return tail_residue_name_; }

	/// @brief Set the r0 value.
	///
	void set_r0(core::Real const &val) {
		runtime_assert_string_msg( val>=0.0,  "In protocols::helical_bundle::MakeBundleHelix::set_r0(): the value of r0 must be greater than or equal to 0." );
		bundle_parameters_->set_r0(val);
		return;
	}

	/// @brief Get the r0 value.
	///
	core::Real r0() const { return bundle_parameters_->r0(); }

	/// @brief Set the omega0 value.
	///
	void set_omega0(core::Real const &val) { bundle_parameters_->set_omega0(val); return; }

	/// @brief Get the omega0 value.
	///
	core::Real omega0() const { return bundle_parameters_->omega0(); }

	/// @brief Set the delta_omega0 value.
	///
	void set_delta_omega0( core::Real const &delta_omega0_in ) { bundle_parameters_->set_delta_omega0(delta_omega0_in); return; }

	/// @brief Get the delta_omega0 value.
	///
	core::Real delta_omega0() const { return bundle_parameters_->delta_omega0(); }

	/// @brief Get the r1 value for a given mainchain atom.
	///
	core::Real r1( core::Size const index ) const { return bundle_parameters_->r1(index); }

	/// @brief Const-access to the whole r1 vector.
	utility::vector1 < core::Real > const & r1_vect() const { return bundle_parameters_->r1_vect(); }

	/// @brief Set the omega1 value.
	///
	void set_omega1(core::Real const &val) { bundle_parameters_->set_omega1(val); return; }

	/// @brief Get the omega1 value.
	///
	core::Real omega1() const { return bundle_parameters_->omega1(); }

	/// @brief Set the z1 value.
	///
	void set_z1(core::Real const &val) { bundle_parameters_->set_z1(val); return; }

	/// @brief Get the z1 value.
	///
	core::Real z1() const { return bundle_parameters_->z1(); }

	/// @brief Get the delta_omega1 value for a given mainchain atom.
	///
	core::Real delta_omega1( core::Size const index ) const { return bundle_parameters_->delta_omega1(index); }

	/// @brief Const-access to the whole delta_omega1 vector.
	utility::vector1 < core::Real > const & delta_omega1_vect() const { return bundle_parameters_->delta_omega1_vect(); }

	/// @brief Get the delta_z1 value for a given mainchain atom.
	///
	core::Real delta_z1( core::Size const index ) const { return bundle_parameters_->delta_z1(index); }

	/// @brief Const-access to the whole delta_z1 vector.
	utility::vector1 < core::Real > const & delta_z1_vect() const { return bundle_parameters_->delta_z1_vect(); }

	/// @brief Get the delta_t value.
	core::Real delta_t() const { return bundle_parameters_->delta_t(); }

	/// @brief Set the major helix parameters
	void set_major_helix_params (
		core::Real const &r0_in,
		core::Real const &omega0_in,
		core::Real const &delta_omega0_in
	) {
		runtime_assert_string_msg( r0_in>=0.0,  "In protocols::helical_bundle::MakeBundleHelix::set_major_helix_params(): the value of r0 must be greater than or equal to 0." );
		bundle_parameters_->set_r0(r0_in);
		bundle_parameters_->set_omega0(omega0_in);
		bundle_parameters_->set_delta_omega0(delta_omega0_in);
		return;
	}

	/// @brief Set the minor helix parameters
	///
	void set_minor_helix_params (
		utility::vector1 < core::Real > const &r1_in,
		core::Real const &omega1_in,
		core::Real const &z1_in,
		utility::vector1 < core::Real > const &delta_omega1_in,
		utility::vector1 < core::Real > const &delta_z1_in
	) {
		bundle_parameters_->set_r1(r1_in);
		bundle_parameters_->set_omega1(omega1_in);
		bundle_parameters_->set_z1(z1_in);
		bundle_parameters_->set_delta_omega1(delta_omega1_in);
		bundle_parameters_->set_delta_z1(delta_z1_in);
		return;
	}

	/// @brief Set delta_t
	///
	void set_delta_t ( core::Real const &delta_t_in) { bundle_parameters_->set_delta_t(delta_t_in); return; }

	/// @brief Set global omega1 offset
	/// @detail The overall omega1 offset is the sum of this global value and the per-atom
	/// delta_omega1 values.
	void set_delta_omega1_all ( core::Real const &delta_omega1_all_in ) { bundle_parameters_->set_delta_omega1_all(delta_omega1_all_in); return; }

	/// @brief Get global omega1 offset
	/// @detail The overall omega1 offset is the sum of this global value and the per-atom
	/// delta_omega1 values.
	core::Real delta_omega1_all () const { return bundle_parameters_->delta_omega1_all(); }

	/// @brief Set the minor helix parameters by reading them in from a file.
	///
	void set_minor_helix_params_from_file ( std::string const &filename )
	{
		utility::vector1 <core::Real> r1;
		core::Real omega1(0.0);
		core::Real z1(0.0);
		utility::vector1 <core::Real> delta_omega1;
		utility::vector1 <core::Real> delta_z1;
		read_minor_helix_params ( filename, r1, omega1, z1, delta_omega1, delta_z1 );
		bundle_parameters_->set_r1(r1);
		bundle_parameters_->set_omega1(omega1);
		bundle_parameters_->set_z1(z1);
		bundle_parameters_->set_delta_omega1(delta_omega1);
		bundle_parameters_->set_delta_z1(delta_z1);
		return;
	}

	/// @brief Returns "true" if the last call to the apply function failed, false otherwise.
	///
	bool last_apply_failed() const { return last_apply_failed_; }


	/// @brief Sets whether this mover is allowed to set mainchain dihedrals.
	///
	void set_allow_dihedrals( bool const val ) { bundle_parameters_->set_allow_dihedrals(val); return; }

	/// @brief Sets whether this mover is allowed to set mainchain bond angles.
	///
	void set_allow_bondangles( bool const val ) { bundle_parameters_->set_allow_bondangles(val); return; }

	/// @brief Sets whether this mover is allowed to set mainchain bond lengths.
	///
	void set_allow_bondlengths( bool const val ) { bundle_parameters_->set_allow_bondlengths(val); return; }


	/// @brief Returns "true" if and only if this mover is allowed to set mainchain dihedrals.
	///
	bool allow_dihedrals() const { return bundle_parameters_->allow_dihedrals(); }

	/// @brief Returns "true" if and only if this mover is allowed to set mainchain bond angles.
	///
	bool allow_bondangles() const { return bundle_parameters_->allow_bondangles(); }

	/// @brief Returns "true" if and only if this mover is allowed to set mainchain bond lengths.
	///
	bool allow_bondlengths() const { return bundle_parameters_->allow_bondlengths(); }

	/// @brief Const-access to bundle parameters.
	///
	BundleParametersCOP bundle_parameters() const { return bundle_parameters_; }

	/// @brief Set a new set of bundle parameters.
	///
	void set_bundle_parameters( BundleParametersOP newparams ) {  bundle_parameters_ = newparams; return; }

private:
////////////////////////////////////////////////////////////////////////////////
//          PRIVATE DATA                                                      //
////////////////////////////////////////////////////////////////////////////////

	/// @brief Should the input pose be reset?
	/// @details Default true.
	bool reset_pose_;

	/// @brief An owning pointer for the BundleParameters object.
	/// @details The BundleParameters object holds all of the Crick parameters for creating
	/// a helical bundle.  This gets passed to the pose that is created or modified, and
	/// stored in the Conformation object, so that other movers can modify or perturb the
	/// object created (or just read its parameters).
	BundleParametersOP bundle_parameters_;

	/// @brief Length of the helix, in residues.  Defaults to 10.
	///
	core::Size helix_length_;

	/// @brief Name (full-length, not 3-letter code) of the residue type that the helix will be made of.
	/// @details Defaults to "ALA".
	std::string residue_name_;

	/// @brief Name (full-length, not 3-letter code) of the residue type that will cap the helix.
	/// @details If blank, there will be no cap.  Defaults to blank.
	std::string tail_residue_name_;

	/// @brief Did the last apply fail?
	/// @details Initialized to "false"; "true" if the last apply failed.
	bool last_apply_failed_;

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

	/// @brief If there are tail residues, set their backbone dihedral angles
	/// to something reasonable (a helical conformation).
	void set_tail_dihedrals( core::pose::Pose &helixpose ) const;

	/// @brief Set whether the last call to the apply() function failed or not.
	/// @details Should only be called by the apply() function.
	void set_last_apply_failed( bool const val ) { last_apply_failed_=val; return; }

	/// @brief Add Crick parameter information to the Conformation object within the pose.
	/// @details This function updates the bundle_parameters_ object's links to residues within the pose,
	/// and then copies the owning pointer into the pose's Conformation object.
	void add_parameter_info_to_pose( core::pose::Pose &pose );

};

} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_MakeBundleHelix_hh
