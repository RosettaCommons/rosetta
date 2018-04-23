// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphOptions.hh
/// @brief Options container for the BuriedUnsatPenaltyGraph.  Initialized by the BuriedUnsatPenalty energy method.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


#ifndef INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_graph_BuriedUnsatPenaltyGraphOptions_hh
#define INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_graph_BuriedUnsatPenaltyGraphOptions_hh

// Relevant headers
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphOptions.fwd.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace buried_unsat_penalty {
namespace graph {

///@brief Options container for the BuriedUnsatPenaltyGraph.  Initialized by the BuriedUnsatPenalty energy method.
class BuriedUnsatPenaltyGraphOptions : public utility::pointer::ReferenceCount {

public:

	/// @brief Default constructor -- explicitly deleted.
	BuriedUnsatPenaltyGraphOptions() = delete;

	/// @brief Options constructor.
	BuriedUnsatPenaltyGraphOptions(
		core::Real const angle_exponent,
		core::Real const angle_shift_factor,
		core::Real const dist_exponent,
		core::Real const dist_midpoint,
		core::Real const burial_threshold,
		core::Real const hbond_energy_threshold
	);

	/// @brief Copy constructor.
	BuriedUnsatPenaltyGraphOptions(BuriedUnsatPenaltyGraphOptions const & /*src*/) = default;

	/// @brief Default destructor.
	virtual ~BuriedUnsatPenaltyGraphOptions();

	/// @brief Clone method: return a copy of the original object, by owning pointer.
	BuriedUnsatPenaltyGraphOptionsOP
	clone() const;

	/// @brief Print information about this object.
	void show( std::ostream & out ) const;

public: //Getters

	/// @brief Get the angle exponent, for determining burial by the method of sidechain neighbor cones.
	inline core::Real angle_exponent() const { return angle_exponent_; }

	/// @brief Get the angle shift factor, for determining burial by the method of sidechain neighbor cones.
	inline core::Real angle_shift_factor() const { return angle_shift_factor_; }

	/// @brief Get the distance exponent, for determining burial by the method of sidechain neighbor cones.
	inline core::Real dist_exponent() const { return dist_exponent_; }

	/// @brief Get the distance midpoint, for determining burial by the method of sidechain neighbor cones.
	inline core::Real dist_midpoint() const { return dist_midpoint_; }

	/// @brief Get the minimum number of cones that a point must fall within, for determining burial by the
	/// method of sidechain neighbor cones.
	inline core::Real burial_threshold() const { return burial_threshold_; }

	/// @brief Get the maximum hydrogen bond energy allowed, above which a hydrogen bond is not counted.
	inline core::Real hbond_energy_threshold() const { return hbond_energy_threshold_; }

private:

	/// @brief The angle exponent, for determining burial by the method of sidechain neighbor cones.
	core::Real angle_exponent_;

	/// @brief The angle shift factor, for determining burial by the method of sidechain neighbor cones.
	core::Real angle_shift_factor_;

	/// @brief The distance exponent, for determining burial by the method of sidechain neighbor cones.
	core::Real dist_exponent_;

	/// @brief The distance midpoint, for determining burial by the method of sidechain neighbor cones.
	core::Real dist_midpoint_;

	/// @brief The minimum number of cones that a point must fall within, for determining burial by the
	/// method of sidechain neighbor cones.
	core::Real burial_threshold_;

	/// @brief The maximum hydrogen bond energy allowed, above which a hydrogen bond is not counted.
	core::Real hbond_energy_threshold_;

};


} //core
} //pack
} //guidance_scoreterms
} //buried_unsat_penalty
} //graph



#endif //INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_graph_BuriedUnsatPenaltyGraphOptions_hh
