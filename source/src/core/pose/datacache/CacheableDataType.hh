// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
		// The enum starts at 1 for indexing into vector1

		// General items
		JOBDIST_OUTPUT_TAG = 1, // a basic::datacache::CacheableString; Set by the jobdistributor
		ARBITRARY_FLOAT_DATA, // a basic::datacache::CacheableStringFloatMap (basically a std::map<string,Real>)
		ARBITRARY_STRING_DATA, // a basic::datacache::CacheableStringMap (basically a std::map<string,string>)
		STRING_MAP, // a basic::datacache::CacheableStringMap; string-based annotations about a Pose  (What's the distinction to ARBITRARY_STRING_DATA?)
		SCORE_LINE_STRINGS, // a basic::datacache::CacheableStringMap; score entries composed of a pair< string, string > (What's the distinction to ARBITRARY_STRING_DATA?)
		SCORE_MAP, // basic::datacache::DiagnosticData; (basically a std::map<string,Real> for scores) // ??? Need more info on usage (What's the distinction to ARBITRARY_FLOAT_DATA?)

		// General pose-associated data
		STM_STORED_TASKS, // a protocols::toolbox::task_operations::STMStoredTask; used by StoreTaskMover and related
		STORED_RESIDUE_SUBSET, // a core::select::residue_selector::CachedResidueSubset -- For storing residue subsets
		CONSTRAINT_GENERATOR, // a protocols::constraint_generator::ConstraintGenerator -- For constraint generator data

		// Pose annotation
		POSE_BEFORE_CAVITIES_ADDED, // a core::pose::datacache::CacheablePoseRawPtr; used by AddCavitiesMover
		TEMPLATE_HYBRIDIZATION_HISTORY,  // a protocols::hybridization:TemplateHistory; during template hybridization, the source id for each residue
		NCS_RESIDUE_MAPPING, // a protocols::symmetry::NCSResMapping; used for Non-Crystollographic Symmetry
		FULL_MODEL_INFO, // a core::pose::full_model_info::FullModelInfo; [map residues/chains to full model, for stepwise buildup]
		VDW_REP_SCREEN_INFO, // a protocols::stepwise::modeler::rns::checker::VDW_CachedRepScreenInfo; for stepwise modeling -- grid of peripheral regions that are sterically disallowed.
		CDR_CLUSTER_INFO, // a protocols::antibody::clusters::BasicCDRClusterSet ( Cacheable Antibody CDR Cluster Information)
		NATIVE_ANTIBODY_SEQ, // a protocols::antibody::design::NativeAntibodySeq -- For keeping track of the near-native sequence during antibody design.
		WRITEABLE_DATA, // a basic::datacache::WriteableCacheableMap Used by the Environment/Broker and stored in SilentStructs
		POSITION_CONSERVED_RESIDUES, // a core::pose::datacache::PositionConservedResiduesStore; used by core::pose::is_position_conserved_residue() (currently unused)
		INTERFACE_DDG_MUTATION, // a InterfaceDDGMutationTask; the mutation to be considered; set by the InterfaceDDGJobInputter (in the interface_ddg_bind app)

		// Energy-function information cache (including experimental data)
		BASE_PARTNER, // a core::scoring::dna::BasePartner; used by DNA energy terms
		CEN_LIST_INFO, // a core::scoring::CenListInfo; used by EnvPairPotential
		SIGMOID_WEIGHTED_CEN_LIST, // a core::scoring::SigmoidWeightedCenListReal; used by SmoothEnvPairPotential
		SIGMOID_WEIGHTED_D_CEN_LIST, // a core::scoring::SigmoidWeightedCenListReal; used by SmoothEnvPairPotential
		RG_MINDATA, // a core::scoring::methods::RG_MinData; used by RG_Energy_Fast
		RG_LOCAL_MINDATA, // a core::scoring::methods::RG_Local_MinData; used by RG_LocalEnergy
		MEMBRANE_TOPOLOGY, // a core::scoring::MembraneTopology; used for (pre-Framework) membrane poses
		MEMBRANE_EMBED, // a core::scoring::MembraneEmbed; used for (pre-Framework) membrane poses
		INTERFACE_INFO, // a protocols/scoring/InterfaceInfo; used for InterchainPotential
		SS_INFO, // a core::scoring::SS_Info; used by SecondaryStructurePotential
		SS_KILLHAIRPINS_INFO, // a a core::scoring::SS_Killhairpins_Info; used by SecondaryStructurePotential
		RNA_SCORING_INFO, // a core::scoring::rna::RNA_ScoringInfo; used by RNA scoring
		RNA_SECSTRUCT_INFO, // a protocols::farna::secstruct::RNA_SecStructLegacyInfo
		WATER_PACKING_INFO, // a core::pack::rotamer_set::WaterPackingInfo; used for making water rotamers
		HOLES_POSE_INFO, // a core::id::CacheableAtomID_MapVector; used by HolesEnergy
		DAB_SASA_POSE_INFO, // a core::id::CacheableAtomID_MapVector; used by SurfEnergy
		DAB_SEV_POSE_INFO, // a core::id::CacheableAtomID_MapVector; used by SurfEnergy
		CHEMICAL_SHIFT_ANISOTROPY_DATA, // a core::scoring::ChemicalShiftAnisotropy; for NMR Chemical Shift Anisotropy
		RESIDUAL_DIPOLAR_COUPLING_DATA, // a core::scoring::ResidualDipolarCoupling; for NMR Residual Dipolar Coupling data
		RESIDUAL_DIPOLAR_COUPLING_DATA_ROHL, // a core::scoring::ResidualDipolarCoupling_Rohl; for NMR Residual Dipolar Coupling data - for rdc_rohl
		RESIDUAL_DIPOLAR_COUPLING_SEGMENTS_DATA, // a core::scoring::ResidualDipolarCouplingRigidSegments
		FIBER_DIFFRACTION_CEN_SCATTERING, // a core::scoring::fiber_diffraction::CentroidScatter; for Fiber Diffraction scattering form factors in centroid mode
		FIBER_DIFFRACTION_FA_SCATTERING, // a core::scoring::fiber_diffraction::FAScatter; for Fiber Diffraction scattering form factors in all-atom mode
		PSEUDOCONTACT_SHIFT_DATA, // a protocols::scoring::methods::pcs::PCS_data; for NMR Psuedocontact Shift (PCS) data SHOULD DESEAPPEAR END 2010
		TS1_PSEUDOCONTACT_SHIFT_DATA, // a protocols::scoring::methods::pcsTs1::PCS_data_Ts1
		TS2_PSEUDOCONTACT_SHIFT_DATA, // a protocols::scoring::methods::pcsTs2::PCS_data_Ts2
		TS3_PSEUDOCONTACT_SHIFT_DATA, // a protocols::scoring::methods::pcsTs3::PCS_data_Ts3
		TS4_PSEUDOCONTACT_SHIFT_DATA, // a protocols::scoring::methods::pcsTs4::PCS_data_Ts4
		PSEUDOCONTACT_SHIFT_MULTI_DATA, // a protocols::scoring::methods::pcs2::PcsDataCenterManager; NMR Psuedocontact Shift (PCS) data, multi paramgnetic center
		GEN_BORN_POSE_INFO, // a core::scoring::GenBornPoseInfo; used by GenBornEnergy
		MULTIPOLE_POSE_INFO, // a core::scoring::MultipoleElecPoseInfo; used by MultipoleElecEnergy
		VDWTINKER_POSE_INFO, // a core::scoring::VdWTinkerPoseInfo; used by VdWTinkerPotential
		FACTS_POSE_INFO, // a core::scoring::FACTSPoseInfo; used by FACTSEnergy
		MEMBRANE_FAEMBED, // a core::scoring::Membrane_FAEmbed; pba high reslution membrane embedding info
		LK_BALL_POSE_INFO, // a core::scoring::lkball::LKB_ResiduesInfo (aka LKB_PoseInfo); used by LK_BallEnergy
		PB_LIFETIME_CACHE, // a core::scoring::methods::PBLifetimeCache; Poisson-boltzmann energy state dependent data (see PoissonBoltzmannEnergy)
		FAELEC_CONTEXT_DATA, // a core::scoring::elec::FAElecContextData; used by FA_GrpElecEnergy
		POSE_SEQUENCE, // a core::scoring::electron_density_atomwise::PoseSeqeunce

		// Old, unused terms
		//  MEMBRANE_POTENTIAL,
		//  RB_JUMP,
		//  SITE_CST,
		//  DOCK_ENS_CONF1,
		//  DOCK_ENS_CONF2,
		//  FILTER_STAGE2_BEGINNING,
		//  FILTER_STAGE2_QUARTER,
		//  FILTER_STAGE2_HALF,
		//  FILTER_STAGE2_END,
		//  VALL_LOOKBACK_DATA,
		//  SEQUENCE_PROFILE,
		//  CUSTOM_PAIR_DIST_SCORE_INFO,
		//  STRUCTURAL_CONSERVATION,
		//  SURFACE_PARAMS,
		//  MULTIPLE_POSE_INFO, // pose/multiple_pose_info/MultiplePoseInfo.cc [info on sister poses, for stepwise buildup]
		//  FLOATING_POINT_CLOUD_INFO,

		// *** IMPORTANT ***  // The 'num_cacheable_data_types' below must be the last enum, and must
		// always be set equal to the (last-2) enum. The 'dummy_cacheable_data_type'
		// below must be the (last-1) enum, so that 'num_cacheable_data_types' does
		// NOT require a manual update everytime a new cachaeable_data_type is added.
		dummy_cacheable_data_type,
		num_cacheable_data_types = dummy_cacheable_data_type - 1
	};

	static std::string get_name( CacheableDataType::Enum datatype);

private:
	static void initialize_name_map();

	static utility::vector1< std::string > & enum2name_();
	static std::map< std::string, CacheableDataType::Enum > &name2enum_();

}; // class CacheableDataType

} // namespace datacache
} // namespace pose
} // namespace core

#endif /* INCLUDED_core_pose_datacache_CacheableDataType_HH */
