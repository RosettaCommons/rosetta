// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pockets/DarcParticleSwarmMinimizer.hh
/// @brief
/// @author Karen R. Khar


#ifndef INCLUDED_protocols_pockets_DarcParticleSwarmMinimizer_hh
#define INCLUDED_protocols_pockets_DarcParticleSwarmMinimizer_hh

#include <core/optimization/ParticleSwarmMinimizer.hh>
#include <protocols/pockets/DarcParticleSwarmMinimizer.fwd.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/types.hh>
#include <core/conformation/Residue.hh>
#include <utility/vector1.hh>
#include <basic/gpu/GPU.hh>

namespace protocols {
namespace pockets {

class DarcParticleSwarmMinimizer : public core::optimization::ParticleSwarmMinimizer {
public:

  using core::optimization::ParticleSwarmMinimizer::run;

  DarcParticleSwarmMinimizer( FingerprintManager & fpmanager_in, NonPlaidFingerprint & nfp_in, PlaidFingerprint & pfp_in,
		core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight,
		core::optimization::Multivec p_min, core::optimization::Multivec p_max) :
		core::optimization::ParticleSwarmMinimizer(p_min, p_max),
		fpmanager_( fpmanager_in ),
		nfp_( nfp_in ),
		pfp_( pfp_in ),
		missing_pt_(missing_point_weight),
		steric_(steric_weight),
		extra_pt_(extra_point_weight){};

  ~DarcParticleSwarmMinimizer() {};

  void score_all_particles(core::optimization::Multifunc & f_fitness, core::optimization::ParticleOPs & particles);

  private:

  void fill_atom_arrays_( core::Size particle_inx, core::conformation::ResidueCOP ligand_rsd, std::vector<basic::gpu::float4> & atoms, std::vector<basic::gpu::float4> & atom_maxmin_phipsi, std::vector<basic::gpu::float4> & ligand_maxmin_phipsi );
  void fill_atom_arrays_multiple_origin_( core::Size particle_inx, core::conformation::ResidueCOP ligand_rsd, std::vector<basic::gpu::float4> & atoms, std::vector<basic::gpu::float4> & atom_maxmin_phipsi, std::vector<basic::gpu::float4> & ligand_maxmin_phipsi, numeric::xyzVector<core::Real> const & origin_offset );
  void fill_atom_arrays_for_electrostatics_( core::Size particle_inx, core::pose::Pose ligand_pose_for_elec_calc, std::vector<basic::gpu::float4> & atoms_coors_and_charge );
  core::Real DarcPSO_fp_compare_( core::Size particle_inx, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, std::vector<basic::gpu::float4> & atoms, std::vector<basic::gpu::float4> & atom_maxmin_phipsi );
  core::Real DarcPSO_fp_compare_multiple_origin_( core::Size particle_inx, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, std::vector<basic::gpu::float4> & atoms, std::vector<basic::gpu::float4> & atom_maxmin_phipsi );
	core::Real DarcPSO_elsts_score_( core::Size particle_inx, core::Size dim_x, core::Size dim_y, core::Size dim_z, core::Real mid_x, core::Real mid_y, core::Real mid_z, core::Real spacing, std::vector < std::vector < std::vector <core::Real> > > espGrid, std::vector < std::vector < std::vector <ElectrostaticpotentialGrid::PtType> > > typGrid, std::vector<basic::gpu::float4> & atom_coors_charge );

	FingerprintManager & fpmanager_;
  NonPlaidFingerprint & nfp_;
  PlaidFingerprint & pfp_;
  core::Real missing_pt_;
  core::Real steric_;
  core::Real extra_pt_;
  core::Size ligand_natoms_;

}; // DarcParticleSwarmMinimizer

} // namespace pockets
} // namespace protocols

#endif // INCLUDED_protocols_pockets_DarcParticleSwarmMinimizer_HH

