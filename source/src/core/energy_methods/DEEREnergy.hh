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
#include <core/scoring/epr_deer/DEERData.hh>
#include <core/scoring/epr_deer/DEERDataCache.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/energy_methods/DEEREnergy.fwd.hh>
#include <core/energy_methods/DEEREnergyCreator.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/DenseEnergyContainer.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// C++ headers
#include <iosfwd>
#include <string>

namespace core {
namespace energy_methods {

class DEEREnergy : public scoring::methods::ContextDependentLRTwoBodyEnergy
{

public: // Methods

	/// @brief  Constructor
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
	scoring::methods::EnergyMethodOP
	clone() const override;

	/// @brief  Destructor
	~DEEREnergy() override;

	/// @brief  Returns the parsed input file (const)
	scoring::epr_deer::DEERDataCacheCOP const
	get_const_deer_data(
		pose::Pose const & pose
	) const;

	/// @brief  Returns the parsed input file (non-const)
	scoring::epr_deer::DEERDataCacheOP
	get_deer_data(
		pose::Pose & pose
	) const;

	/// @brief Calculates the score for a specific data set
	/// @detail Called by setup_for_scoring() and setup_for_derivatives.
	///     EPRSpinLabel is used to calculate where the spin label
	///     is/can be. Edge ID corresponds to the specific set of
	///     data being interrogated/scored. Modifier is used for derivative
	///     calculation since the determination of derivatives is not
	///     analytically solvable in many cases
	Real
	get_score(
		pose::Pose & pose,
		scoring::epr_deer::EPRSpinLabel & sl,
		Size const & edge_id,
		Real const & modifier = 0.0
	) const;

	/// @brief Overrides virtual fxn for finding scores for a given pose
	void
	setup_for_scoring(
		pose::Pose & pose,
		scoring::ScoreFunction const &
	) const override;

	/// @brief Initialize energy method
	void
	initialize_energy_method(
		pose::Pose & pose
	) const;

	/// @brief Derived function for specifying where the data is getting stored
	scoring::methods::LongRangeEnergyType
	long_range_type() const override;

	/// @brief  A declared but unused function
	void
	setup_for_minimizing(
		pose::Pose &,
		scoring::ScoreFunction const &,
		kinematics::MinimizerMapBase const &
	) const override;

	/// @brief  Returns "false" by default
	bool
	defines_intrares_energy(
		scoring::EnergyMap const &
	) const override;

	/// @brief Fetches if two residues are "connected" (i.e. there is a score)
	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size rsd1,
		Size rsd2
	) const override;

	/// @brief Get energy/score for a given pair of residues.
	/// @detail Because some sets of data may be between more than two residues,
	///     this function looks at a map kept in DEERDataCache specifically
	///     maintained for this purpose
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		scoring::ScoreFunction const &,
		scoring::EnergyMap & emap
	) const override;

	/// @brief  Calculates the derivatives for minimization. Done numerically
	void
	setup_for_derivatives(
		pose::Pose & pose,
		scoring::ScoreFunction const &
	) const override;

	/// @brief Returns the res-specific F1 and F2 vectors saved in DEERDataCache
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

	/// @brief  A declared but unused function
	void
	indicate_required_context_graphs(
		utility::vector1< bool > &
	) const override;

	/// @brief  Version 1 (as of 30 January 2017)
	Size
	version() const override;

};


} // namespace energy_methods
} // namespace core

#endif // INCLUDED_core_energy_methods_DEEREnergy_hh
