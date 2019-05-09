// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/membrane/ImplicitLipidInfo.hh
/// @brief An Implicit Lipid Membrane Model
/// @version menv-franklin2018
///
/// @detail This class defines physical and chemical properties of the membrane environment. It
/// is mainly used by the membrane energy function.
///  1. Parameters of the hydration function (e.g. thickness, rate, of transition, pore size)
///  2. Lipid composition details
///  3. Hydration function smoothing parameters
///  4. Structure-based lipid accessibiity information
///
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_ImplicitLipidInfo_hh
#define INCLUDED_core_conformation_membrane_ImplicitLipidInfo_hh

#include <core/conformation/membrane/ImplicitLipidInfo.fwd.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

// Package headers
#include <core/conformation/membrane/SpanningTopology.fwd.hh>
#include <core/conformation/membrane/AqueousPoreParameters.fwd.hh>
#include <core/types.hh>

// Numeric Headers
#include <numeric/MathMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/cubic_polynomial.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace conformation {
namespace membrane {

///@brief Definition of an implicit membrane with parameters for different lipid compositions
class ImplicitLipidInfo : public utility::pointer::ReferenceCount {

	typedef utility::vector1< numeric::CubicPolynomial > piecewise_poly;

public:

	/// @brief Custom constructor for ImplicitLipidInfo
	ImplicitLipidInfo(
		std::string lipid_composition_name,
		core::Real temperature
	);

	ImplicitLipidInfo( ImplicitLipidInfo const & src );
	ImplicitLipidInfo & operator=( ImplicitLipidInfo const & src );

	virtual ~ImplicitLipidInfo();

	ImplicitLipidInfoOP
	clone() const;

	/// @brief Generate a string representation of information represented by ths Lipid MembraneInfo
	virtual void show( std::ostream & output ) const;

public: // Getters for information about the lipid composition

	/// @brief Number of carbons and degrees of saturation in the chains
	std::string chain_type() const;

	/// @brief Chemical name of the headgroup
	std::string headgroup_type() const;

	/// @brief Abbreviated name of the lipid composition
	std::string lipid_composition_name() const;

	/// @brief Full name of the lipid composiiton
	std::string lipid_composition_name_long() const;

	/// @brief Number of degrees of saturation in the overall lipid
	core::Real degrees_of_saturation() const;

	/// @brief Temperature at which the Db parameter was measured/calculated (in celcius)
	core::Real temperature() const;

public: // Getters for information about the membrane pore

	/// @brief Does this protein have a water-filled pore?
	bool has_pore() const;
	void has_pore( bool const is_there_a_pore );

	/// @brief Is the protein alpha helical or beta barrel
	bool is_helical() const;
	void is_helical( bool const is_helical );

	/// @brief Make a set of parameters for the case with no pore
	void make_no_pore_parameters();

	/// @brief Membrane Aqueous pore center - h parameter
	core::Real pore_center_x( core::Real const zcoord ) const;

	/// @brief Membrane aqueous pore center - k parameter
	core::Real pore_center_y( core::Real const zcoord ) const;

	/// @brief Membrane aqueous pore - major radius
	core::Real pore_major_radius( core::Real const zcoord ) const;

	/// @brief Membrane aqueous pore - minor radius
	core::Real pore_minor_radius( core::Real const zcoord ) const;

	/// @brief Membrane aqueous pore - rotation matrix
	numeric::MathMatrix< core::Real > pore_rotation( core::Real const zcoord ) const;

	/// @brief Set membrane aqueous pore parameters
	void set_aqueous_pore_parameters( AqueousPoreParametersOP aqueous_pore );

public: // Lipid Accessibility Information

	/// @brief Per-atom accessibility to membrane lipids (0 = not exposed, 50 = exposed)
	core::Real
	per_atom_lipid_accessibility( core::Size rsd, core::Size atom ) const;

	/// @brief Per-atom accessibility to membrane lipids (0 = not exposed, 50 = exposed)
	void
	set_per_atom_lipid_accessibility(
		utility::vector1< utility::vector1< core::Real > > v);

public: // Chemical information about this membrane

	/// @brief Water thickness of the membrane
	/// @detail Thickness of the membrane, defined by the pair of z coordinates
	/// with a water density of 50% (SAXS Db value)
	core::Real water_thickness() const;

	/// @brief Change in water density from membrane core to water bulk water
	/// @detail Steepness defined by the number of waters lost per increase in z
	/// from the membrane center (s = b)
	core::Real water_steepness() const;

	/// @brief Access function to change the value of the transition steepness
	/// @detail Only change the steepness for the alpha v. beta case
	void water_steepness( core::Real const v ) { change_in_water_density_ = v; }

	/// @brief Pseudo-thickness parameter
	/// @details A parameter in the expotential membrane definition
	/// t = -(1/b) * ln(1/A)
	core::Real water_pseudo_thickness() const;

public: // Functions to define spatially-dependent hydration

	/// UNDER CONSTRUCTION - Derivative functions ///

	/// @brief Calcuclate the hydration of an atom based on its location in a
	/// lipid-specific implicit bilayer
	core::Real
	f_depth( core::Real const z ) const;

	/// @brief Ccalcculate the gradient of f_depth
	core::Real
	f_depth_gradient( core::Real const z ) const;

	/// @brief Calculate the hydration of an atom based on its location relative to
	/// an aqueous pore or cavity
	core::Real
	f_cavity( numeric::xyzVector< core::Real > const & p ) const;

	/// @brief Calculate the derivative of f_cavity (without any r(x,y,z) dependence)
	core::Real
	f_cavity_gradient( core::Real const r ) const;

	/// @brief Calculate the location of an atom relative to the pore structure
	core::Real
	g_radius( numeric::xyzVector< core::Real > const & p ) const;

	/// @breif Calcuclate the gradient of f_cavity
	core::Real
	g_radius_gradient( numeric::xyzVector< core::Real > const & p ) const;

	/// @brief Overall hydration given the atomic depth and cavity structure
	core::Real
	f_hydration( numeric::xyzVector< core::Real > const & p ) const;

	/// @brief Gradient
	core::Real
	f_hydration_gradient( numeric::xyzVector< core::Real > const & p ) const;

public: // Helper functions - making public for testing

	/// @brief Helper function to initialize the lipid composiiton data
	void
	initialize_implicit_lipid_parameters();

private:

	// Disallow use of the default constructor
	ImplicitLipidInfo();

private:

	// Parameters for the logistic the implicit membrane definition
	core::Real water_thickness_;
	core::Real change_in_water_density_;
	core::Real transformed_water_thickness_;

	// Parameters about this membrane composition
	std::string chain_type_;
	std::string headgroup_type_;
	std::string lipid_composition_name_;
	std::string lipid_composition_name_long_;
	core::Real degrees_of_saturation_;
	core::Real temperature_;

	// Does this protein have a pore?
	bool has_pore_;

	// Is this protein alpha helical or beta barre?
	bool is_helical_;

	// aqueous pore parameters
	AqueousPoreParametersOP pore_params_;
	core::Real pore_transition_steepness_;

	// lipid accessibility
	utility::vector1< utility::vector1< core::Real > > per_atom_lipid_accessibility_;

#ifdef    SERIALIZATION
	friend class cereal::access;
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // membrane
} // conformation
} // core



#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_membrane_ImplicitLipidInfo )
#endif // SERIALIZATION


#endif //INCLUDED_core_conformation_membrane_ImplicitLipidInfo_hh





