// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/energy_methods/DEEREnergy.hh
/// @brief  Score term for data obtained with double electron-electron resonance (DEER)
/// @details This method is initiated by assigning a weight to deer_decay and providing an input
///       text file via the option -epr_deer:input_files data.txt.
///
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_energy_methods_DEEREnergy_hh
#define INCLUDED_core_energy_methods_DEEREnergy_hh

// Unit headers
#include <core/scoring/epr_deer/DEERDataCache.fwd.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.fwd.hh>
#include <core/energy_methods/DEEREnergy.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzVector.fwd.hh>

// C++ headers

#include <core/scoring/methods/EnergyMethodOptions.fwd.hh> // AUTO IWYU For EnergyMethodOptions
#include <map> // AUTO IWYU For map

namespace core {
namespace energy_methods {

class DEEREnergy : public scoring::methods::ContextDependentLRTwoBodyEnergy
{

public: // Methods

	/// @brief  Constructor
	/// @param  energymethodcreator: Creator object
	DEEREnergy(
		scoring::methods::EnergyMethodCreatorOP energymethodcreator
	);

	/// @brief  Constructor for obtaining command-line options
	DEEREnergy(
		scoring::methods::EnergyMethodOptions const &
	);

	/// @brief  Default Constructor
	DEEREnergy();

	/// @brief   Copy constructor
	DEEREnergy(
		DEEREnergy const & other
	);

	/// @brief Copy function
	/// @return Pointer to new object that is a copy of this object
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/// @brief  Destructor
	~DEEREnergy() override;

	///////////////////////////
	/// INHERITED FUNCTIONS

	/// @brief Inherited function specifying where the data is getting stored
	/// @return Descriptor of this energy method
	scoring::methods::LongRangeEnergyType
	long_range_type() const override;


	/// @brief Set up scoring process
	/// @param pose: Pose being scored
	void
	setup_for_scoring(
		pose::Pose & pose,
		scoring::ScoreFunction const &
	) const override;

	/// @brief Add energy/score for a given pair of residues to map.
	/// @param rsd1: First residue
	/// @param rsd2: Second residue
	/// @param pose: Pose being scored
	/// @param emap: Where we are adding the score
	/// @detail Because some sets of data may be between more than two
	/// @detail residues, this function looks at a map kept in
	/// @detail DEERDataCache specifically maintained for this purpose
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		scoring::ScoreFunction const &,
		scoring::EnergyMap & emap
	) const override;

	/// @brief  A declared but unused function
	void
	setup_for_minimizing(
		pose::Pose &,
		scoring::ScoreFunction const &,
		kinematics::MinimizerMapBase const &
	) const override;

	/// @brief  Returns "false" by default
	/// @return Returns false (DEER cannot evaluate anything intra-residue)
	bool
	defines_intrares_energy(
		scoring::EnergyMap const &
	) const override;

	/// @brief Fetches if two residues are "connected" (i.e. there is a score)
	/// @param pose: Pose being scored
	/// @param rsd1: First residue
	/// @param rsd2: Second residue
	/// @return Whether the pair of residues are linked by data
	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size rsd1,
		Size rsd2
	) const override;

	/// @brief  Calculates the derivatives for minimization. Done numerically
	/// @param pose: Pose being scored/modified
	void
	setup_for_derivatives(
		pose::Pose & pose,
		scoring::ScoreFunction const &
	) const override;

	/// @brief  A declared but unused function
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		scoring::ScoreFunction const &,
		scoring::EnergyMap & emap
	) const override;

	/// @brief  A declared but unused function
	void
	finalize_total_energy(
		pose::Pose &,
		scoring::ScoreFunction const &,
		scoring::EnergyMap &
	) const override;

	/// @brief  Apply atom derivative to atom
	/// @param id: Atom ID in pose
	/// @param pose: Pose being modified
	/// @param emap: Energy map with scores used by derivatives
	/// @param F1 vector to modify for atom
	/// @param F2 vector to modify for atom
	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &,
		scoring::ScoreFunction const & /*sfxn*/,
		scoring::EnergyMap const & emap,
		numeric::xyzVector< Real > & F1,
		numeric::xyzVector< Real > & F2
	) const override;

	/// @brief  A declared but unused function
	void
	indicate_required_context_graphs(
		utility::vector1< bool > &
	) const override;

	/// @brief  Version 1 (as of 3 January 2021)
	/// @return Version 1
	Size
	version() const override;

	////////////////////////////////////////////
	/// NON-VIRTUAL FUNCTIONS

	/// @brief  Returns the DEERDataCache in the pose (const)
	/// @param pose: Pose with DEERDataCache stashed in it
	/// @return Pointer to const DEERDataCache object
	/// @detail This will crash if the DEERDataCache object is absent!
	scoring::epr_deer::DEERDataCacheCOP const
	get_const_deer_data(
		pose::Pose const & pose
	) const;

	/// @brief  Returns the DEERDataCache in the pose (non-const)
	/// @param pose: Pose with DEERDataCache stashed in it
	/// @return Pointer to non-const DEERDataCache object
	/// @detail This will crash if the DEERDataCache object is absent!
	scoring::epr_deer::DEERDataCacheOP
	get_deer_data(
		pose::Pose & pose
	) const;

	/// @brief This checks if multiple sets of spin labels is being used
	/// @param pose: Pose being scored
	/// @param mod: Modifier to X-axis (for calculating derivatives)
	/// @return Map of edges/data indices to scores for those edges/indices
	/// @detail This is used when spin labels are obtained for specific
	/// @detail proteins a priori, e.g., by triangulation/multilateration,
	/// @detail and the solution is underdetermined, so multiple solutions
	/// @detail are required.
	std::map< Size, Real >
	iter_over_labels(
		pose::Pose & pose,
		int const & mod = 0
	) const;

	/// @brief Calculates the score for a specific data set
	/// @param pose: Pose being scored
	/// @param sl: Spin label being used to obtain score
	/// @param edge_id: Index of data being scored (in DEERDataCache)
	/// @param mod: Modifier to X-axis (zero when scored)
	/// @return Score corresponding to the data
	Real
	get_score(
		pose::Pose & pose,
		scoring::epr_deer::EPRSpinLabel & sl,
		Size const & edge_id,
		int const & mod = 0
	) const;

	/// @brief Initializes energy method
	/// @param pose: Pose being evaluated/used for scoring
	void
	initialize_energy_method(
		pose::Pose & pose
	) const;

};


} // namespace energy_methods
} // namespace core

#endif // INCLUDED_core_energy_methods_DEEREnergy_hh
