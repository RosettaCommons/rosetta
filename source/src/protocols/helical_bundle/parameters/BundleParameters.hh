// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/parameters/BundleParameters.hh
/// @brief  Prototypes and method declarations for the BundleParameters class, a class for holding parameters for helical bundle backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_protocols_helical_bundle_parameters_BundleParameters_hh
#define INCLUDED_protocols_helical_bundle_parameters_BundleParameters_hh

// Unit headers
#include <protocols/helical_bundle/parameters/BundleParameters.fwd.hh>

// Package headers
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Numeric headers

// C++ headers


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace helical_bundle {
namespace parameters {

/// @brief  Parameters class, used to store sets of parameters for parametric backbone generation.
///
class BundleParameters : public core::conformation::parametric::Parameters
{

public: //Typedefs:

	typedef core::conformation::parametric::Parameters Parameters;
	typedef core::conformation::parametric::ParametersOP ParametersOP;
	typedef core::conformation::parametric::ParametersSet ParametersSet;
	typedef core::conformation::parametric::ParametersSetOP ParametersSetOP;

public:

	/// @brief constructors
	///
	BundleParameters();

	BundleParameters( BundleParameters const & src );

	~BundleParameters();

	/// @brief Copy this residue( allocate actual memory for it )
	///
	ParametersOP clone() const;

	///////////////////
	// Setters:
	///////////////////


	/// @brief Sets the major helix radius.
	/// @details In Angstroms.
	void set_r0 ( core::Real const &val) {
		runtime_assert_string_msg(val>=0.0, "In core::conformation::parametric::Parameters::set_r0(): The major helix radius must be greater than or equal to 0.");
		r0_ = val;
		return;
	}

	/// @brief Sets the major helix twist per residue.
	/// @details In radians.
	void set_omega0 ( core::Real const &val) {
		omega0_ = val;
		return;
	}

	/// @brief Sets the major helix angular offset about the major helix (bundle) axis.
	/// @details In radians.
	void set_delta_omega0 ( core::Real const &val) {
		delta_omega0_ = val;
		return;
	}

	/// @brief Sets the number of residues per repeating unit in the minor helix.
	///
	void set_residues_per_repeat( core::Size const val ) {
		runtime_assert_string_msg( val > 0, "In protocols::helical_bundle::parameters::BundleParameters::set_residues_per_repeat(): The number of residues per repeat must be greater than 0." );
		residues_per_repeat_ = val;
		return;
	}

	/// @brief Sets the residue offset in the repeating unit in the minor helix.
	/// @details An offset of 0 means that the first residue in the helix is the first
	/// residue of the repeating unit.  An offset of 1 means that the first residue in
	/// the helix is the second residue in the repeating unit, etc.
	void set_repeating_unit_offset( core::Size const val ) {
		repeating_unit_offset_ = val;
		return;
	}

	/// @brief Set the number of residues in the repeating unit in the helix.
	/// @details Note that there is no check that the size of the atoms_per_residue_ vector matches residues_per_repeat_.
	void set_atoms_per_residue( utility::vector1 <core::Size> const &atoms_per_residue_in ) {
		atoms_per_residue_ = atoms_per_residue_in;
		return;
	}

	/// @brief Resets the minor helix radius list.
	///
	void reset_r1 () {
		r1_.clear();
		return;
	}

	/// @brief Sets the minor helix radius list to a pre-defined list.
	/// @details Values in Angstroms.
	void set_r1( utility::vector1 < core::Real > const &vals) {
		r1_.clear();
		r1_ = vals;
		return;
	}

	/// @brief Adds a value to the minor helix radius list.
	/// @details In Angstroms
	void add_r1_value( core::Real const &val ) {
		runtime_assert_string_msg(val>=0.0, "In core::conformation::parametric::Parameters::add_r1_value(): The minor helix radius must be greater than or equal to 0.");
		r1_.push_back(val);
		return;
	}

	/// @brief Sets the minor helix twist per residue.
	/// @details In radians.
	void set_omega1 ( core::Real const &val) {
		omega1_ = val;
		return;
	}

	/// @brief Sets the minor helix angular offset around the minor helix axis, defined for all atoms (to "roll" the minor helix).
	/// @details In radians.
	void set_delta_omega1_all ( core::Real const &val) {
		delta_omega1_all_ = val;
		return;
	}

	/// @brief Sets the minor helix rise per residue, defined for all atoms.
	/// @details In Angstroms.
	void set_z1 ( core::Real const &val) {
		z1_ = val;
		return;
	}

	/// @brief Resets the list of minor helix angular offsets around the minor helix axis, defined on a per-atom basis.
	///
	void reset_delta_omega1 () {
		delta_omega1_.clear();
		return;
	}

	/// @brief Sets the list of minor helix angular offsets around the minor helix axis, defined on a per-atom basis, to a pre-defined list.
	/// @details Values in radians.
	void set_delta_omega1( utility::vector1 < core::Real > const &vals) {
		delta_omega1_.clear();
		delta_omega1_ = vals;
		return;
	}

	/// @brief Adds a value to the list of minor helix angular offsets around the minor helix axis, defined on a per-atom basis.
	/// @details In radians.
	void add_delta_omega1_value( core::Real const &val ) {
		delta_omega1_.push_back(val);
		return;
	}

	/// @brief Resets the list of minor helix offsets along the minor helix axis, defined on a per-atom basis.
	///
	void reset_delta_z1 () {
		delta_z1_.clear();
		return;
	}

	/// @brief Sets the list of minor helix offsets along the minor helix axis, defined on a per-atom basis, to a pre-defined list.
	/// @details Values in Angstroms.
	void set_delta_z1( utility::vector1 < core::Real > const &vals) {
		delta_z1_.clear();
		delta_z1_ = vals;
		return;
	}

	/// @brief Adds a value to the list of minor helix offsets along the minor helix axis, defined on a per-atom basis.
	/// @details In radians.
	void add_delta_z1_value( core::Real const &val ) {
		delta_z1_.push_back(val);
		return;
	}

	/// @brief Sets the minor helix offset along the minor helix axis, defined on a per-helix basis, in Angstroms.
	///
	void set_z1_offset( core::Real const &val ) { z1_offset_=val; return; }

	/// @brief Sets the minor helix offset along the major helix axis, defined on a per-helix basis, in Angstroms.
	///
	void set_z0_offset( core::Real const &val ) { z0_offset_=val; return; }

	/// @brief Sets the bundle lateral squash factor, defined on a per-helix basis.
	/// @details A value of 0.5 means that the bundle's minor axis, in the y direction, will be half as long as the bundle's
	/// major axis, in the x direction.  (The helices curl around the z-axis).
	void set_epsilon( core::Real const &val ) { epsilon_=val; return; }

	/// @brief Set whether the helix direction be reversed.
	///
	void set_invert_helix( bool const val) { invert_helix_=val; return; }

	/// @brief Set the residue offset, to "slide" the helix along its path.
	/// @details In residues (i.e. a value of 0.5 is a half-residue offset).
	void set_delta_t( core::Real const &val ) { delta_t_=val; return; }

	/// @brief Set whether the generator should be able to set mainchain dihedral values.
	///
	void set_allow_dihedrals( bool const val) { allow_dihedrals_=val; return; }

	/// @brief Set whether the generator should be able to set mainchain bond angle values.
	///
	void set_allow_bondangles( bool const val) { allow_bondangles_=val; return; }

	/// @brief Set whether the generator should be able to set mainchain bond length values.
	///
	void set_allow_bondlengths( bool const val) { allow_bondlengths_=val; return; }

	///////////////////
	// Getters:
	///////////////////

	/// @brief Returns the major helix radius.
	/// @details In Angstroms.
	inline core::Real r0() const { return r0_; }

	/// @brief Returns the major helix twist per residue.
	/// @details In radians.
	inline core::Real omega0() const { return omega0_; }

	/// @brief Returns the major helix angular offset about the major helix (bundle) axis.
	/// @details In radians.
	inline core::Real delta_omega0() const { return delta_omega0_; }

	/// @brief Returns the number of residues per repeating unit in the minor helix.
	///
	inline core::Size residues_per_repeat() const { return residues_per_repeat_; }

	/// @brief Returns the residue offset in the repeating unit in the minor helix.
	/// @details An offset of 0 means that the first residue in the helix is the first
	/// residue of the repeating unit.  An offset of 1 means that the first residue in
	/// the helix is the second residue in the repeating unit, etc.
	inline core::Size repeating_unit_offset() const { return repeating_unit_offset_; }

	/// @brief Return the number of atoms in one of the residues in the repeating unit of the minor helix.
	///
	inline core::Size atoms_per_residue( core::Size const index_in_repeating_unit ) const {
		runtime_assert_string_msg( index_in_repeating_unit <= atoms_per_residue_.size() && index_in_repeating_unit > 0,
			"Error in protocols::helical_bundle::parameters::BundleParameters::atoms_per_residue():  The index provided to atoms_per_residue() is out of range." );
		return atoms_per_residue_[index_in_repeating_unit];
	}

	/// @brief Return the atoms_per_residue vector (const-access).
	///
	inline utility::vector1 < core::Size > const & atoms_per_residue() const { return atoms_per_residue_; }

	/// @brief Returns the minor helix r1 value for the atom with the index value, in Angstroms.
	/// @details The index value is checked to see whether it is in range ONLY in debug-mode compiliation.
	inline core::Real r1( core::Size const index ) const {
		debug_assert( index<=r1_.size() && index>0 );
		return r1_[index];
	}

	/// @brief Const-access to the whole r1 vector.
	utility::vector1 < core::Real > const & r1_vect() const { return r1_; }

	/// @brief Returns the minor helix turn per residue, defined for all atoms.
	/// @details In radians.
	inline core::Real omega1() const { return omega1_; }

	/// @brief Returns the minor helix angular offset around the minor helix axis, defined for all atoms (to "roll" the minor helix).
	/// @details In radians.
	inline core::Real delta_omega1_all() const { return delta_omega1_all_; }

	/// @brief Returns the minor helix rise per residue, defined for all atoms.
	/// @details In Angstroms.
	inline core::Real z1() const { return z1_; }

	/// @brief Returns the minor helix angular offset around the minor helix axis, defined on a per-atom basis, in radians, for the atom with the index value.
	/// @details The index value is checked to see whether it is in range ONLY in debug-mode compiliation.
	inline core::Real delta_omega1( core::Size const index ) const {
		debug_assert( index<=delta_omega1_.size() && index>0 );
		return delta_omega1_[index];
	}

	/// @brief Const-access to the whole delta_omega1 vector.
	utility::vector1 < core::Real > const & delta_omega1_vect() const { return delta_omega1_; }

	/// @brief Returns the minor helix offset along the minor helix axis, defined on a per-atom basis, in Angstroms, for the atom with the index value.
	/// @details The index value is checked to see whether it is in range ONLY in debug-mode compiliation.
	inline core::Real delta_z1( core::Size const index ) const {
		debug_assert( index<=delta_z1_.size() && index>0 );
		return delta_z1_[index];
	}

	/// @brief Returns the minor helix offset along the minor helix axis, defined on a per-helix basis, in Angstroms.
	///
	inline core::Real z1_offset() const { return z1_offset_; }

	/// @brief Returns the minor helix offset along the major helix axis, defined on a per-helix basis, in Angstroms.
	///
	inline core::Real z0_offset() const { return z0_offset_; }

	/// @brief Returns the bundle lateral squash factor, defined on a per-helix basis.
	/// @details A value of 0.5 means that the bundle's minor axis, in the y direction, will be half as long as the bundle's
	/// major axis, in the x direction.  (The helices curl around the z-axis).
	inline core::Real epsilon() const { return epsilon_; }

	/// @brief Const-access to the whole delta_z1 vector.
	utility::vector1 < core::Real > const & delta_z1_vect() const { return delta_z1_; }

	/// @brief Should the helix direction be reversed?
	///
	inline bool invert_helix() const { return invert_helix_; }

	/// @brief Residue offset, to "slide" the helix along its path.
	/// @details In residues (i.e. a value of 0.5 is a half-residue offset).
	inline core::Real delta_t() const { return delta_t_; }

	/// @brief Should the generator set mainchain dihedral values?
	///
	inline bool allow_dihedrals() const { return allow_dihedrals_; }

	/// @brief Should the generator set mainchain bond angle values?
	///
	inline bool allow_bondangles() const { return allow_bondangles_; }

	/// @brief Should the generator set mainchain bond length values?
	///
	inline bool allow_bondlengths() const { return allow_bondlengths_; }

	///////////////////
	//Output:
	///////////////////

	/// @brief Get a summary of this ParametersSet object, for output to remark lines of a PDB file.
	///
	virtual void get_pdb_remark(std::stringstream &remark) const;

private:

	/********************************************************************************
	PRIVATE DATA
	*********************************************************************************/

	///////////////////
	//Major helix data:
	///////////////////

	/// @brief Major helix radius.
	/// @details In Angstroms.
	core::Real r0_;

	/// @brief Major helix twist per residue.
	/// @details In radians.
	core::Real omega0_;

	/// @brief Major helix angular offset about the major helix (bundle) axis.
	/// @details In radians.
	core::Real delta_omega0_;

	///////////////////
	//Minor helix data:
	///////////////////

	/// @brief The number of residues in the repeating unit in the helix.
	/// @details Defaults to 1.
	core::Size residues_per_repeat_;

	/// @brief The reside offset in the repeating unit.
	/// @details Defaults to 0.
	core::Size repeating_unit_offset_;

	/// @brief The number of atoms in each residue in the repeating unit in the helix.
	///
	utility::vector1 < core::Size > atoms_per_residue_;

	/// @brief Minor helix radius, defined on a per-atom basis.
	/// @details In Angstroms.
	utility::vector1 < core::Real > r1_;

	/// @brief Minor helix turn per residue, defined for all atoms.
	/// @details In radians.
	core::Real omega1_;

	/// @brief Minor helix angular offset around the minor helix axis, defined for all atoms (to "roll" the minor helix).
	/// @details In radians.
	core::Real delta_omega1_all_;

	/// @brief Minor helix rise per residue, defined for all atoms.
	/// @details In Angstroms.
	core::Real z1_;

	/// @brief Minor helix angular offset around the minor helix axis, defined on a per-atom basis.
	/// @details In radians.
	utility::vector1 < core::Real > delta_omega1_;

	/// @brief Minor helix offset along the minor helix axis, defined on a per-atom basis.
	/// @details In Angstroms.
	utility::vector1 < core::Real > delta_z1_;

	/// @brief Minor helix offset along the minor helix axis, defined on a per-helix basis.
	/// @details In Angstroms.
	core::Real z1_offset_;

	/// @brief Major helix offset along the major helix axis, defined on a per-helix basis.
	/// @details In Angstroms.
	core::Real z0_offset_;

	/// @brief Bundle lateral squash factor.
	/// @details Dimensionless.
	core::Real epsilon_;

	///////////////////
	//Other data:
	///////////////////

	/// @brief Should the helix direction be reversed?
	///
	bool invert_helix_;

	/// @brief Residue offset, to "slide" the helix along its path.
	/// @details In residues (i.e. a value of 0.5 is a half-residue offset).
	core::Real delta_t_;

	/// @brief Should the generator set mainchain dihedral values?
	///
	bool allow_dihedrals_;

	/// @brief Should the generator set mainchain bond angle values?
	///
	bool allow_bondangles_;

	/// @brief Should the generator set mainchain bond length values?
	///
	bool allow_bondlengths_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class BundleParameters

} // namespace parametric
} // namespace conformation
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_helical_bundle_parameters_BundleParameters )
#endif // SERIALIZATION


#endif
