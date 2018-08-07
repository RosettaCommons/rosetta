// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/BundleParametrizationCalculator.hh
/// @brief  Prototypes and method declarations for the BundleParametrizationCalculator class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_protocols_conformation_parametric_BundleParametrizationCalculator_hh
#define INCLUDED_protocols_conformation_parametric_BundleParametrizationCalculator_hh


// Unit headers
#include <protocols/helical_bundle/BundleParametrizationCalculator.fwd.hh>
#include <core/conformation/parametric/ParametrizationCalculator.hh>
#include <core/conformation/parametric/Parameter.fwd.hh>

// Package headers
#include <protocols/helical_bundle/parameters/BundleParametersSet.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.fwd.hh>
#include <core/conformation/parametric/RealValuedParameter.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Numeric headers

// C++ headers
#include <tuple>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace helical_bundle {

/// @brief Enum for relevant parameters.
/// @details If you add to this list, update the parameter_name_from_enum() function, the short_parameter_description_from_enum() function, the parameter_units_from_enum() function, the
/// parameter_type_from_enum() function, the parameter_properties_from_enum() function, and the parameter_description_from_enum()
/// function.
enum BPC_Parameters {
	//Parameters to be sampled.  LIST THESE FIRST, ALWAYS:
	BPC_r0 = 1, //Keep first -- radius from z-axis (Angstroms)
	BPC_omega0, //Twist per residue of the superhelix (radians)
	BPC_delta_omega0, //Rotation of helix about the z-axis (radians)
	BPC_delta_omega1, //Rotation of helix about its own axis (radians)
	BPC_delta_t, //Offset along the polypeptide backbone (residues)
	BPC_z0_offset, //Offset along the Z-axis (Angstroms)
	BPC_z1_offset, //Offset along the superhelix path (Angstroms)
	BPC_epsilon, //Lateral squash
	BPC_last_parameter_to_be_sampled = BPC_epsilon, //Keep this with previous

	//Invariant parameters (minor helix parameters):
	BPC_residues_per_repeat,
	BPC_repeating_unit_offset,
	BPC_atoms_per_residue,
	BPC_r1_peratom, //Minor helix radius -- vector for each atom, in Angstroms.
	BPC_omega1, //Twist per residue of the minor helix (single value, radians).
	BPC_z1, //rise per residue (single value, Angstroms)
	BPC_delta_omega1_peratom, //Rotational offset about helix axis -- vector for each atom, in radians.
	BPC_delta_z1_peratom,  //Offset along helix axis -- vector for each atom, in radians.

	//Other data (invariant):
	BPC_invert_helix, //Flip the helix top-to-bottom? Boolean.
	BPC_set_dihedrals, //Allow dihedral values to be set.
	BPC_set_bondangles, //Allow bond angles to deviate from ideality.
	BPC_set_bondlengths, //Allow bond lengths to deviate from ideality.

	BPC_unknown_parameter, //Keep second-to-last
	BPC_end_of_list = BPC_unknown_parameter //Keep last
};

/// @brief  BundleParametrizationCalculator class, used for parametric helical bundle backbone generation.
///
class BundleParametrizationCalculator : public core::conformation::parametric::ParametrizationCalculator
{
public:

	typedef core::conformation::parametric::ParameterizationCalculatorProperties ParameterizationCalculatorProperties;

public:

	/// @brief constructors
	///
	BundleParametrizationCalculator( bool const use_degrees=false );

	/// @brief Params constructor.
	BundleParametrizationCalculator( bool const use_degrees, protocols::helical_bundle::parameters::BundleParametersCOP params_in );

	/// @brief Copy constructor.
	/// @details Deep-copies the stored parameters.
	BundleParametrizationCalculator( BundleParametrizationCalculator const & src );

	/// @brief Destructor.
	~BundleParametrizationCalculator() override;

	/// @brief Copy this object ( allocate actual memory for it )
	core::conformation::parametric::ParametrizationCalculatorOP clone() const override;

public: //Functions

	/// @brief Read a Crick params file and initialize this calculator.
	/// @details Triggers a read from disk!
	void init_from_file( std::string const &filename );

	/// @brief Copy the parameter values for parameters that have not been set from the global parameters.
	/// @details This function should be called before build_helix().
	void copy_unset_params_from_globals( BundleParametrizationCalculatorCOP global_calculator );

	/// @brief Copy the parameter values for parameters that copy values from previous helices, from the previous helices.
	/// @details This function should be called before build_helix().
	/// @note This function assumes that the input pose consists *only* of previous helices, and stores N ParametersSet objects
	/// for N helices, each with 1 Parameters object contained.
	/// @returns Returns true for failure, false for success.
	bool copy_params_from_previous_helices_makebundlehelix_style( core::pose::Pose const & prev_helices_pose );

	/// @brief Copy the parameter values for parameters that copy values from previous helices, from the previous helices.
	/// @details This function should be called before build_helix().
	/// @note This function assumes that the input pose is the whole shebang passed to the PerturbBundle mover.  Unlike
	/// copy_params_from_previous_helices_makebundlehelix_style(), this functon does *not* modify the
	/// BundleParametrizationCalculator.  Instead it copies the parameter value for the given parameter from the
	/// indicated previous helix to the indicated parameter object.
	/// @returns Returns true for failure, false for success.
	static bool copy_params_from_previous_helices_perturbbundle_style( parameters::BundleParametersSetCOP paramset, core::Size const copy_from_helix_index, BPC_Parameters const parameter_type, core::conformation::parametric::RealValuedParameterOP parameter_to_copy_to, protocols::helical_bundle::parameters::BundleParametersCOP current_helix_params );

	/// @brief Attempts to build a helix based on the current Crick parameters.
	/// @details Returns true for FAILURE and false for SUCCESS.  The object "pose" is replaced with a
	/// new helix if this operation succeeds.
	bool build_helix( core::pose::Pose &pose ) const;

	/// @brief Attempts to build a helix based on the current Crick parameters.
	/// @details Returns true for FAILURE and false for SUCCESS.  The object "pose" is replaced with a
	/// new helix if this operation succeeds.
	/// @note This version takes parameters for the start and end of the helix.
	bool build_helix( core::pose::Pose &pose, core::Size const helix_start, core::Size const helix_end ) const;

	/// @brief Get the parameter type from the BPC_Parameters enum.
	/// @details Returns PT_invalid_type if invalid.
	static core::conformation::parametric::ParameterType parameter_type_from_enum( BPC_Parameters param_enum );

	/// @brief Given a BPC_Parameters enum, get a short lay-language description (used for annotating output) of the parameter.
	/// @details Returns "INVALID!!!" if invalid.
	static std::string const & parameter_description_from_enum( BPC_Parameters param_enum );

	/// @brief Get the short parameter name from the BPC_Parameters enum.
	static std::string const & short_parameter_description_from_enum( BPC_Parameters param_type );

	/// @brief Get the parameter units from the BPC_Parameters enum.
	static std::string const & parameter_units_from_enum( BPC_Parameters param_type );

	/// @brief Given a BPC_Parameters enum, get whether this parameter is one that the user can set, whether it can be sampled, whether it can be perturbed, and whether it's global
	/// for a ParametersSet (helical bundle).
	static ParameterizationCalculatorProperties const & parameter_properties_from_enum( BPC_Parameters param_enum );

	/// @brief Get the parameter name from the BPC_Parameters enum.
	static std::string const & parameter_name_from_enum( BPC_Parameters param_type );

	/// @brief Get the BPC_Parameters enum from the parameter name.
	/// @details Returns BPC_unknown_parameter if this can't be parsed.
	static BPC_Parameters parameter_enum_from_name( std::string const &name );

	/// @brief Set whether this calculator uses degrees or radians for the user-set angles.
	/// @details Defaults to radians (false).
	void set_use_degrees( bool const setting );

	/// @brief Get whether this calculator expects input in degrees or radians for user-set angles.
	inline bool use_degrees() const { return use_degrees_; }

	/// @brief Set the perturbation type for all parameters.
	void set_perturbation_type_globally( std::string const & perturbation_type_in );

private: //Private functions

	/// @brief Update whether certain parameters expect their inputs in degrees.
	void set_use_degrees_for_parameters();

private:

	/********************************************************************************
	PRIVATE DATA
	*********************************************************************************/

	bool use_degrees_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class BundleParametrizationCalculator

} // namespace helical_bundle
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_helical_bundle_BundleParametrizationCalculator )
#endif // SERIALIZATION


#endif
