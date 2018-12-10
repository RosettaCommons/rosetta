// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/mainchain_potential/GenerateMainchainPotentialOptions.hh
/// @brief Options container for the generator for mainchain potentials.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_mainchain_potential_GenerateMainchainPotentialOptions_hh
#define INCLUDED_protocols_mainchain_potential_GenerateMainchainPotentialOptions_hh

#include <protocols/mainchain_potential/GenerateMainchainPotentialOptions.fwd.hh>
#include <protocols/mainchain_potential/GenerateMainchainPotentialTests.fwd.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace mainchain_potential {

/// @brief Options container for the generator for mainchain potentials.
class GenerateMainchainPotentialOptions : public utility::pointer::ReferenceCount {

	friend class ::GenerateMainchainPotentialTests; //Needed to allow the unit tests for this class access internal data members.

public:

	/// @brief Default constructor.
	/// @details "True" indicates that we should initialize from the global options system.  "False" by default,
	/// meaning that all options are set to default values.
	GenerateMainchainPotentialOptions( bool const initialize_from_globals = false );

	/// @brief Copy constructor.  Set to default.
	GenerateMainchainPotentialOptions(GenerateMainchainPotentialOptions const & /*src*/) = default;

	/// @brief Destructor.
	virtual ~GenerateMainchainPotentialOptions();

	/// @brief Clone function: create a copy of this object, and return an owning pointer to the copy.
	GenerateMainchainPotentialOptionsOP clone() const;

public: //Static functions

	/// @brief Static function used to register releavnt options.
	static void register_options();

public: //Setters

	/// @brief Set the name of the residue type for which we'll be generating a mainchain potential.
	void set_residue_name( std::string const &name_in );

	/// @brief Set whether we're generating a pre-proline potential.
	/// @details If true, then we generate a potential for a position before a proline, sarcosine, or other N-substituted building block (e.g. a peptoid).  If false,
	/// then the next residue's nitrogen is unsubstituted.
	/// @note Determines the patch applied to the C-terminus.
	void set_make_pre_proline_potential( bool const setting );

	/// @brief Set the vector of number of points in each dimension for the mainchain potential that we're generating.
	/// @details This takes a vector of ints, not Sizes, because that's what comes from the global options system.
	void set_dimensions( utility::vector1< int > const & dimensions_in );

	/// @brief Set which mainchain torsions are covered.
	/// @details If this is an empty vector, it indicates that all mainchain torsions are covered.  This takes a vector of
	/// ints, not Sizes, because that's what comes from the global options sytem.
	/// @note Pass this an empty vector to set the default (all mainchain torsions for the residue type.)
	void set_mainchain_torsions_covered( utility::vector1< int > const & mainchain_torsions_in );

	/// @brief Set name of the scorefunction that we'll use.  An empty string indicates that the default should be used.
	void set_scorefxn_filename( std::string const & filename_in );

	/// @brief Set the Boltzmann temperature.
	void set_kbt( core::Real const & kbt_in );

	/// @brief Set whether we're minimizing each rotamer.
	void set_do_minimization( bool const setting );

	/// @brief Set the minimization type.
	void set_minimization_type( std::string const & type_in );

	/// @brief Set the minimization threshold.
	void set_minimization_threshold( core::Real const & threshold_in );

	/// @brief Set the mainchain potential filename to write.
	void set_output_filename( std::string const & filename_in );

	/// @brief Set whether we are writing mainchain potentials for individual scoreterms.
	void set_write_potentials_for_individual_scoreterms( bool const setting );

	/// @brief Set whether the output should be made symmetric.
	void set_symmetrize_output( bool const setting );

public: //Getters

	/// @brief Get the name of the residue type for which we'll be generating a mainchain potential.
	inline std::string const & residue_name() const { return residue_name_; }

	/// @brief Get whether we're generating a pre-proline potential.
	/// @details If true, then we generate a potential for a position before a proline, sarcosine, or other N-substituted building block (e.g. a peptoid).  If false,
	/// then the next residue's nitrogen is unsubstituted.
	/// @note Determines the patch applied to the C-terminus.
	inline bool make_pre_proline_potential() const { return make_pre_proline_potential_; }

	/// @brief Access the vector of the number of points in each dimension for a given residue type.  (Const-access)/
	/// @details The length of this vector must match the number of mainchain dihedral DoFs for which we're
	/// generating a mainchain potential.
	inline utility::vector1< core::Size > const & dimensions() const { return dimensions_; }

	/// @brief Get which mainchain torsions are covered.
	/// @details If this is an empty vector, it indicates that all mainchain torsions are covered.
	inline utility::vector1< core::Size > const & mainchain_torsions_covered() const { return mainchain_torsions_covered_; }

	/// @brief Name of the scorefunction that we'll use.  An empty string indicates that the default should be used.
	inline std::string const & get_scorefxn_filename() const { return scorefxn_filename_; }

	/// @brief Get the Boltzmann temperature.
	inline core::Real const & kbt() const { return kbt_; }

	/// @brief Get whether we're minimizing each rotamer.
	inline bool do_minimization() const { return do_minimization_; }

	/// @brief Get the minimization type.
	std::string const & minimization_type() const { return minimization_type_; }

	/// @brief Get the minimization threshold.
	inline core::Real const & minimization_threshold() const { return minimization_threshold_; }

	/// @brief Get the mainchain potential filename to write.
	inline std::string const & output_filename() const { return output_filename_; }

	/// @brief Are we writing mainchain potentials for individual scoreterms?
	/// @details Default false.
	inline bool write_potentials_for_individual_scoreterms() const { return write_potentials_for_individual_scoreterms_; }

	/// @brief Should the output be made symmetric?
	inline bool symmetrize_output() const { return symmetrize_output_; }

private: //Methods

	/// @brief Initialize this object from the global options system.
	void do_initialization_from_globals();

private: //Data

	/// @brief Name of the residue type for which we'll be generating a mainchain potential.
	std::string residue_name_;

	/// @brief Name of the scorefunction that we'll use.
	/// @details If an empty string (no user-specified scorefunction), then the GenerateMainchainPotential class has a default
	/// that it uses.
	std::string scorefxn_filename_;

	/// @brief If true, then we generate a potential for a position before a proline, sarcosine, or other N-substituted building block (e.g. a peptoid).  If false,
	/// then the next residue's nitrogen is unsubstituted.
	/// @details Determines the patch applied to the C-terminus.
	bool make_pre_proline_potential_;

	/// @brief The number of points in each dimension for a given residue type.
	/// @details The length of this vector must match the number of mainchain dihedral DoFs for which we're
	/// generating a mainchain potential.
	utility::vector1< core::Size > dimensions_;

	/// @brief Which mainchain torsions are covered?
	/// @details If this is an empty vector, it indicates that all mainchain torsions are covered.
	utility::vector1< core::Size > mainchain_torsions_covered_;

	/// @brief The Boltzmann temperature for average energy calculation.
	/// @details Defaults to 0.63 kcal/mol.  I'm not sure whether this is right for the mm_std scorefunction, though.
	core::Real kbt_;

	/// @brief Should each rotamer be minimized?  Default true.
	bool do_minimization_;

	/// @brief Minimization flavour to use to approximate the Hessian matrix.  Defaults to dfpmin (which is really lbfgs).
	std::string minimization_type_;

	/// @brief The minimization threshold at which minimization trajectories are aborted.  Defaults to 1.0E-7.
	core::Real minimization_threshold_;

	/// @brief The mainchain potential filename to write.
	std::string output_filename_;

	/// @brief Are we writing mainchain potentials for individual scoreterms?
	/// @details Default false.
	bool write_potentials_for_individual_scoreterms_;

	/// @brief Should the output be made symmetric?  False by default.
	bool symmetrize_output_;

};


} //protocols
} //mainchain_potential



#endif //INCLUDED_protocols_mainchain_potential_GenerateMainchainPotentialOptions_hh





