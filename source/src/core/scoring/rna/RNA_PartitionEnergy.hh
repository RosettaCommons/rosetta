// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RNA_PartitionEnergy.hh
/// @brief Compute RNA score for secondary structure partition function. Requires Vienna's RNAfold.
/// @author Ramya Rangan


#ifndef INCLUDED_core_scoring_rna_RNA_PartitionEnergy_hh
#define INCLUDED_core_scoring_rna_RNA_PartitionEnergy_hh


// Package headers
#include <core/scoring/annealing/ResidueArrayAnnealableEnergy.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace scoring {
namespace rna {


class RNA_PartitionEnergy : public methods::WholeStructureEnergy, public core::scoring::annealing::ResidueArrayAnnealableEnergy {
public:
	typedef methods::WholeStructureEnergy  parent1;
	typedef annealing::ResidueArrayAnnealableEnergy parent2;

public:
	/// @brief Default Constructor
	RNA_PartitionEnergy();

	/// @brief Copy Constructor
	RNA_PartitionEnergy( RNA_PartitionEnergy const & src );

	/// clone
	methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const override;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & // context_graphs_required
	) const override {}

	/// @brief Calculate the total energy given a vector of const owning pointers to residues.
	/// @details Called directly by the ResidueArrayAnnealingEvaluator during packer runs.
	core::Real calculate_energy( utility::vector1< core::conformation::ResidueCOP > const & resvect,
		core::Size const substitution_position = 0
	) const override;

	/// @brief Cache data from the pose in this EnergyMethod in anticipation of scoring.
	void set_up_residuearrayannealableenergy_for_packing (
		pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSets const &rotamersets,
		scoring::ScoreFunction const & //sfxn
	) override;

private:


	mutable
		std::map< std::string, Real > partition_cache_;
	/// The following are cached from the Pose for packer energy calculations
	utility::vector1< Size > global_mapping_;
	utility::vector1< Size > res_list_;
	std::string native_sequence_;

	void
	get_score_from_sequences(
		std::string const & native_sequence,
		std::string const & current_sequence,
		Real & score
	) const;

	void
	get_partition_score_from_cache(
		std::string const & global_sequence,
		Real & partition_score
	) const;

	void exec_cmd(
		std::string const & cmd,
		std::string & cmd_result
	) const;

	void
	get_partition_score(
		std::string const & global_sequence,
		Real & partition_score
	) const;

	core::Size version() const override;

};


} //rna
} //scoring
} //core

#endif
