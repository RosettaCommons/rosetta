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
	MakeBundleHelix();
	virtual ~MakeBundleHelix();

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	///
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
	void set_invert_helix(bool const invert_in ) { invert_=invert_in; }

	/// @brief Return whether the helix direction will be inverted.
	///
	bool invert_helix() const { return invert_; }

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
		r0_=val;
		return;
	}

	/// @brief Get the r0 value.
	///
	core::Real r0() const { return r0_; }

	/// @brief Set the omega0 value.
	///
	void set_omega0(core::Real const &val) { omega0_=val; return; }

	/// @brief Get the omega0 value.
	///
	core::Real omega0() const { return omega0_; }

	/// @brief Set the delta_omega0 value.
	///
	void set_delta_omega0( core::Real const &delta_omega0_in ) { delta_omega0_=delta_omega0_in; return; }

	/// @brief Get the delta_omega0 value.
	///
	core::Real delta_omega0() const { return delta_omega0_; }

	/// @brief Get the r1 value for a given mainchain atom.
	///
	core::Real r1( core::Size const index ) const { return r1_[index]; }

	/// @brief Set the omega1 value.
	///
	void set_omega1(core::Real const &val) { omega1_=val; return; }

	/// @brief Get the omega1 value.
	///
	core::Real omega1() const { return omega1_; }

	/// @brief Set the z1 value.
	///
	void set_z1(core::Real const &val) { z1_=val; return; }

	/// @brief Get the z1 value.
	///
	core::Real z1() const { return z1_; }

	/// @brief Get the delga_omega1 value for a given mainchain atom.
	///
	core::Real delta_omega1( core::Size const index ) const { return delta_omega1_[index]; }

	/// @brief Get the delta_z1 value for a given mainchain atom.
	///
	core::Real delta_z1( core::Size const index ) const { return delta_z1_[index]; }

	/// @brief Get the delta_t value.
	core::Real delta_t() const { return delta_t_; }

	/// @brief Set the major helix parameters
	void set_major_helix_params (
		core::Real const &r0_in,
		core::Real const &omega0_in,
		core::Real const &delta_omega0_in
	) {
		r0_ = r0_in;
		omega0_ = omega0_in;
		delta_omega0_ = delta_omega0_in;
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
		r1_ = r1_in;
		omega1_ = omega1_in;
		z1_ = z1_in;
		delta_omega1_ = delta_omega1_in;
		delta_z1_ = delta_z1_in;
		return;
	}

	/// @brief Set delta_t_
	///
	void set_delta_t ( core::Real const &delta_t_in) { delta_t_=delta_t_in; return; }

	/// @brief Set global omega1 offset
	/// @detail The overall omega1 offset is the sum of this global value and the per-atom
	/// delta_omega1 values.
	void set_delta_omega1_all ( core::Real const &delta_omega1_all_in ) { delta_omega1_all_=delta_omega1_all_in; return; }

	/// @brief Get global omega1 offset
	/// @detail The overall omega1 offset is the sum of this global value and the per-atom
	/// delta_omega1 values.
	core::Real delta_omega1_all () const { return delta_omega1_all_; }

	/// @brief Set the minor helix parameters by reading them in from a file.
	///
	void set_minor_helix_params_from_file ( std::string const &filename )
	{
		read_minor_helix_params ( filename, r1_, omega1_, z1_, delta_omega1_, delta_z1_ );
		return;
	}

	/// @brief Returns "true" if the last call to the apply function failed, false otherwise.
	///
	bool last_apply_failed() const { return last_apply_failed_; }


	/// @brief Sets whether this mover is allowed to set mainchain dihedrals.
	///
	void set_allow_dihedrals( bool const val ) { allow_dihedrals_=val; return; }

	/// @brief Sets whether this mover is allowed to set mainchain bond angles.
	///
	void set_allow_bondangles( bool const val ) { allow_bondangles_=val; return; }

	/// @brief Sets whether this mover is allowed to set mainchain bond lengths.
	///
	void set_allow_bondlengths( bool const val ) { allow_bondlengths_=val; return; }


	/// @brief Returns "true" if and only if this mover is allowed to set mainchain dihedrals.
	///
	bool allow_dihedrals() const { return allow_dihedrals_; }

	/// @brief Returns "true" if and only if this mover is allowed to set mainchain bond angles.
	///
	bool allow_bondangles() const { return allow_bondangles_; }

	/// @brief Returns "true" if and only if this mover is allowed to set mainchain bond lengths.
	///
	bool allow_bondlengths() const { return allow_bondlengths_; }


private:
////////////////////////////////////////////////////////////////////////////////
//          PRIVATE DATA                                                      //
////////////////////////////////////////////////////////////////////////////////

	/// @brief Should the input pose be reset?
	/// @details Default true.
	bool reset_pose_;

	/// @brief Should the helix direction be inverted?
	/// @details Default false.
	bool invert_;

	/// @brief Length of the helix, in residues.  Defaults to 10.
	///
	core::Size helix_length_;

	/// @brief Name (full-length, not 3-letter code) of the residue type that the helix will be made of.
	/// @details Defaults to "ALA".
	std::string residue_name_;

	/// @brief Name (full-length, not 3-letter code) of the residue type that will cap the helix.
	/// @details If blank, there will be no cap.  Defaults to blank.
	std::string tail_residue_name_;

	/// @brief The r0 value for the major helix.
	///
	core::Real r0_;

	/// @brief The omega0 value for the major helix.
	///
	core::Real omega0_;

	/// @brief The delta_omega0 value for the major helix.
	///
	core::Real delta_omega0_;

	/// @brief The vector of r1 values, one for each mainchain atom in the residue type.
	///
	utility::vector1 < core::Real > r1_;

	/// @brief The omega1 value, constant for all mainchain atoms.
	///
	core::Real omega1_;

	/// @brief The omega1 offset, constant for all mainchain atoms.
	/// @details Though the effective omega1 value ends up being delta_omega1_[atom] + delta_omega1_all_, it is
	/// convenient to keep these separate so that the former can be loaded from a file and the latter can be set
	/// by the user.
	core::Real delta_omega1_all_;

	/// @brief The z1 (rise per residue) value, constant for all mainchain atoms.
	///
	core::Real z1_;

	/// @brief The vector of delta_omega1 values, one for each mainchain atom in the residue type.
	///
	utility::vector1 < core::Real > delta_omega1_;

	/// @brief The vector of delta_z1 values, one for each mainchain atom in the residue type.
	///
	utility::vector1 < core::Real > delta_z1_;

	/// @brief An offset factor for the value of t, the residue index.
	///
	core::Real delta_t_;

	/// @brief Should the mover set dihedral values?
	/// @details True by default.
	bool allow_dihedrals_;

	/// @brief Should the mover set bond angle values?
	/// @details False by default.  If true, the generated backbone is non-ideal but might
	/// conform to a perfect helix of helices more precisely.
	bool allow_bondangles_;

	/// @brief Should the mover set bond length values?
	/// @details False by default.  If true, the generated backbone is non-ideal but might
	/// conform to a perfect helix of helices more precisely.
	bool allow_bondlengths_;

	/// @brief Did the last apply fail?
	/// @details Initialized to "false"; "true" if the last apply failed.
	bool last_apply_failed_;

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

	/// @brief If there are tail residues, set their backbone dihedral angles
	/// to something reasonable (a helical conformation).
	void set_tail_dihedrals( core::pose::Pose &helixpose ) const;

	/// @brief Generate the x,y,z coordinates of the mainchain atoms using the Crick equations.
	/// @details Coordinates will be returned as a vector of vectors of xyzVectors.  The outer
	/// index will refer to residue number, and the inner index will refer to atom number.
	/// Returns failed=true if coordinates could not be generated, false otherwise.
	void generate_atom_positions(
		utility::vector1 < utility::vector1 < numeric::xyzVector< core::Real > > > &outvector,
		core::pose::Pose const &helixpose,
		core::Size const helix_start,
		core::Size const helix_end,
		bool &failed
	) const;

	/// @brief Place the helix mainchain atoms based on the Crick equations.
	///
	void place_atom_positions(
		core::pose::Pose &pose,
		utility::vector1 < utility::vector1 < numeric::xyzVector < core::Real >  > > const &atom_positions,
		core::Size const helix_start,
		core::Size const helix_end
	) const;

	/// @brief Copy backbone bond length values from one pose, where helix mainchain atom coordinates have been
	/// set with the Crick equations, to another with ideal geometry.
	void copy_helix_bondlengths(
		core::pose::Pose &pose,
		core::pose::Pose const &ref_pose,
		core::Size const helix_start,
		core::Size const helix_end
	) const;

	/// @brief Copy backbone bond angle values from one pose, where helix mainchain atom coordinates have been
	/// set with the Crick equations, to another with ideal geometry.
	void copy_helix_bondangles(
		core::pose::Pose &pose,
		core::pose::Pose const &ref_pose,
		core::Size const helix_start,
		core::Size const helix_end
	) const;

	/// @brief Copy backbone dihedral values from one pose, where helix mainchain atom coordinates have been
	/// set with the Crick equations, to another with ideal geometry.
	void copy_helix_dihedrals(
		core::pose::Pose &pose,
		core::pose::Pose const &ref_pose,
		core::Size const helix_start,
		core::Size const helix_end
	) const;

	/// @brief Align mainchain atoms of pose to ref_pose mainchain atoms.
	///
	void align_mainchain_atoms(
		core::pose::Pose &pose,
		core::pose::Pose const &ref_pose,
		core::Size const helix_start,
		core::Size const helix_end
	) const;

	/// @brief Set whether the last call to the apply() function failed or not.
	/// @details Should only be called by the apply() function.
	void set_last_apply_failed( bool const val ) { last_apply_failed_=val; return; }

};

} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_MakeBundleHelix_hh
