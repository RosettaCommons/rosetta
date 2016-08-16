// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/optimize_weights/OptEData.fwd.hh
/// @brief  Forward declarations for OptE data structures
/// @author Jim Havranek


#ifndef INCLUDED_protocols_optimize_weights_OptEData_fwd_hh
#define INCLUDED_protocols_optimize_weights_OptEData_fwd_hh

// utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace optimize_weights {

class PNatAAOptERotamerData;
class PNatRotOptERotamerData;
class SingleStructureData;
class ConstraintedOptimizationWeightFunc;

class OptEPositionData;
class PNatAAOptEPositionData;
class PSSMOptEPositionData;
class PNatRotOptEPositionData;
class PNatStructureOptEData;
class DDGMutationOptEData;

class OptEData;

//// Individual pieces of data held by position containers
typedef utility::pointer::shared_ptr< PNatAAOptERotamerData > PNatAAOptERotamerDataOP;
typedef utility::pointer::shared_ptr< PNatRotOptERotamerData > PNatRotOptERotamerDataOP;
typedef utility::pointer::shared_ptr< SingleStructureData > SingleStructureDataOP;
typedef utility::pointer::shared_ptr< ConstraintedOptimizationWeightFunc > ConstraintedOptimizationWeightFuncOP;
typedef utility::pointer::shared_ptr< SingleStructureData const > SingleStructureDataCOP;

typedef utility::pointer::shared_ptr< OptEPositionData > OptEPositionDataOP;

/// Position containers
typedef utility::pointer::shared_ptr< PNatAAOptEPositionData > PNatAAOptEPositionDataOP;
typedef utility::pointer::shared_ptr< PSSMOptEPositionData > PSSMOptEPositionDataOP;
typedef utility::pointer::shared_ptr< PNatRotOptEPositionData > PNatRotOptEPositionDataOP;
typedef utility::pointer::shared_ptr< PNatStructureOptEData > PNatStructureOptEDataOP;
typedef utility::pointer::shared_ptr< DDGMutationOptEData > DDGMutationOptEDataOP;

typedef utility::pointer::shared_ptr< OptEData > OptEDataOP;

typedef utility::vector1< PNatAAOptERotamerDataOP > PNatAAOptERotamerDataOPs;
typedef utility::vector1< PNatRotOptERotamerDataOP > PNatRotOptERotamerDataOPs;
typedef utility::vector1< SingleStructureDataOP > SingleStructureDataOPs;

typedef utility::vector1< OptEPositionDataOP > OptEPositionDataOPs;

enum OptEPositionDataType {
	prob_native_amino_acid = 1,
	pssm_data,
	prob_native_rotamer,
	prob_native_structure,
	prob_native_ligand_pose,
	dG_binding_correlation,
	ddG_mutation_correlation,
	constrained_optimization_weight_func,
	prob_native_amino_acid_with_unfolded_energy,
	ddG_mutation_correlation_with_unfolded_energy,
	ddG_bind_correlation,
	n_optE_data_types = ddG_bind_correlation // keep this guy last
};

} // namespace optimize_weights
} // namespace protocols


#endif // INCLUDED_protocols_optimize_weights_OptEData_FWD_HH
