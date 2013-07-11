// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/datacache/CacheableDataType.hh
/// @brief  enum indexing the data types stored in a Pose's internal DataCache
#ifndef INCLUDED_core_pose_datacache_CacheableDataType_hh
#define INCLUDED_core_pose_datacache_CacheableDataType_hh

#include <utility/vector1.hh>

#include <string>
#include <map>

namespace core {
namespace pose {
namespace datacache {

// hold the enum within a descriptive namespace to avoid name collisions
class CacheableDataType {
public:

/// @brief enum indexing the data types stored in a Pose's internal DataCache
enum Enum {
	BASE_PARTNER = 1, // for indexing into vector1
	CEN_LIST_INFO,
	SIGMOID_WEIGHTED_CEN_LIST,
	SIGMOID_WEIGHTED_D_CEN_LIST,
	RG_MINDATA,
	RG_LOCAL_MINDATA,
	MEMBRANE_TOPOLOGY,
	MEMBRANE_EMBED,
	MEMBRANE_POTENTIAL, //pba
	INTERFACE_INFO, //<protocols/scoring/InterfaceInfoOP
	RB_JUMP,
	SITE_CST, //<protocols/scoring/SiteConstraintOP
	DOCK_ENS_CONF1,
	DOCK_ENS_CONF2,
	SS_INFO,
	SS_KILLHAIRPINS_INFO,
	RNA_SCORING_INFO,
	RNA_SECSTRUCT_INFO,
	JOBDIST_OUTPUT_TAG, //< a CacheableStringOP
	WATER_PACKING_INFO,
	SCORE_MAP, //<DiagnosticDataOP
	FILTER_STAGE2_BEGINNING,
	FILTER_STAGE2_QUARTER,
	FILTER_STAGE2_HALF,
	FILTER_STAGE2_END,
	ARBITRARY_FLOAT_DATA,
	POSE_BEFORE_CAVITIES_ADDED,
	STM_STORED_TASKS,
	STRING_MAP, // string-based annotations about a Pose
	SCORE_LINE_STRINGS, // score entries composed of a pair< string, string >
	HOLES_POSE_INFO,
	SEQUENCE_PROFILE,
	TEMPLATE_HYBRIDIZATION_HISTORY,  // during template hybridization, the source id for each residue
	DAB_SASA_POSE_INFO,
	DAB_SEV_POSE_INFO,
	CHEMICAL_SHIFT_ANISOTROPY_DATA, //NMR Chemical Shift Anisotropy
	RESIDUAL_DIPOLAR_COUPLING_DATA,//NMR Residual Dipolar Coupling data
	RESIDUAL_DIPOLAR_COUPLING_DATA_ROHL,//NMR Residual Dipolar Coupling data - for rdc_rohl
	RESIDUAL_DIPOLAR_COUPLING_SEGMENTS_DATA,
	PSEUDOCONTACT_SHIFT_DATA,//NMR Psuedocontact Shift (PCS) data SHOULD DESEAPPEAR END 2010
	PSEUDOCONTACT_SHIFT_MULTI_DATA,//NMR Psuedocontact Shift (PCS) data, multi paramgnetic center
	CUSTOM_PAIR_DIST_SCORE_INFO,
	GEN_BORN_POSE_INFO,
	FACTS_POSE_INFO,
	MEMBRANE_FAEMBED, //pba high reslution membrane embedding info
	POSITION_CONSERVED_RESIDUES,
	LK_BALL_POSE_INFO,
	STRUCTURAL_CONSERVATION,
	SURFACE_PARAMS,
	PB_LIFETIME_CACHE, // Poisson-boltzmann energy state dependent data (see PoissonBoltzmannEnergy)
	FULL_MODEL_INFO, // protocols/swa/monte_carlo/SubToFullInfo.cc

	// *** IMPORTANT ***
	// The 'num_cacheable_data_types' below must be the last enum, and must
	// always be set equal to the (last-1) enum.  If you append a new enum
	// to the list, remember to change the value below!
  num_cacheable_data_types = FULL_MODEL_INFO
};

static std::string get_name( CacheableDataType::Enum datatype);

private:
	static void initialize_name_map();

	static bool initialized_;

	static utility::vector1< std::string > enum2name_;
	static std::map< std::string, CacheableDataType::Enum > name2enum_;

}; // class CacheableDataType

} // namespace datacache
} // namespace pose
} // namespace core

#endif /* INCLUDED_core_pose_datacache_CacheableDataType_HH */
