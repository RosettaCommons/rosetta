// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/mhc_epitope_energy/MHCEpitopeConstraint.hh
/// @brief Headers for a constraint for MHC epitope scores
/// @details Follows analogous file for Vikram K. Mulligan's NetChargeEnergy
/// The corresponding energy term for this constraint is the MHCEpitopeEnergy
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

#ifndef INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopeConstraint_hh
#define INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopeConstraint_hh

#include <core/scoring/aa_composition_energy/SequenceConstraint.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeConstraint.fwd.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergySetup.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/XYZ_Func.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <string>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION



namespace core {
namespace scoring {
namespace mhc_epitope_energy {

class MHCEpitopeConstraint : public core::scoring::aa_composition_energy::SequenceConstraint {

public: //Constructor, destructor, copy, clone:

	/// @brief Constructor
	MHCEpitopeConstraint();

	/// @brief Copy constructor
	MHCEpitopeConstraint( MHCEpitopeConstraint const & src );

	/// @brief Destructor
	~MHCEpitopeConstraint();

	/// @brief Clone operator
	virtual constraints::ConstraintOP clone() const;

	virtual
	bool operator == ( constraints::Constraint const & /*other*/ ) const;

	virtual
	bool
	same_type_as_me( constraints::Constraint const & other ) const;


public: //Functions that actually do stuff:

	/// @brief Set the selector to be used by this constraint.
	/// @details Clones the input.
	void set_selector( select::residue_selector::ResidueSelectorCOP selector_in );

	/// @brief Const access to the selector.
	/// @details Returns a null pointer if no selector has been specified.
	select::residue_selector::ResidueSelectorCOP selector() const;

	/// @brief Stores the name of the selector used to mask the constraint.
	void set_selector_name( std::string const & selector_name ) { cst_selector_name_ = selector_name; }

	/// @brief Returns the name of the selector used to mask the constraint.
	std::string get_selector_name() { return cst_selector_name_; }

	/// @brief Const access to the MHCEpitopeEnergySetup object.
	MHCEpitopeEnergySetupCOP
	mhc_epitope_energy_setup() const;

	/// @brief Initialize the MHCEpitopeEnergySetup object from a file.
	void initialize_from_file( std::string const &filename );

	/// @brief Initialize the MHCEpitopeEnergySetup object from the contents of a file.
	/// @details Allows external code to initialize a constriant object without having the
	/// object read directly from disk.
	void initialize_from_file_contents( std::string const &filecontents );

	/// @brief Set the cst_weight_
	void set_cst_weight( core::Real cst_weight );

	/// @brief Get the cst_weight_
	core::Real get_cst_weight() const;

	/// @brief Print info on the constraint
	void show_def (std::ostream &TO, pose::Pose const &pose) const;

private:
	// Member variables

	/// @brief Owning pointer to a ResidueSelector.
	/// @details Optional; will serve as a mask for this MHCEpitopeConstraint if provided.
	select::residue_selector::ResidueSelectorOP selector_;

	/// @brief MHCEpitopeEnergySetup object stored by this object.
	/// @details Created on object construction
	MHCEpitopeEnergySetupOP mhc_epitope_setup_;

	/// @brief A weight score for a particular constraint.
	/// @details An additional weighting factor, beyond that set in the scorefunction, to reset the weight of a single constraint without changing the other predictors.
	/// Useful for increasing/decreasing the importance of one predictor over another.
	core::Real cst_weight_;

	/// @brief The name of the residue selector used to restrict the constraint to part of the pose.
	/// @details Intended only for reporting purposes.  Defaults to "full pose" to indicate no mask.
	std::string cst_selector_name_ = "full pose";

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // mhc_epitope_energy
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_mhc_epitope_energy_MHCEpitopeConstraint )
#endif // SERIALIZATION

#endif
