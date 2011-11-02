// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RMS_Energy.cc
/// @brief  RMS Energy function. Used to optimize the RMSD between two structures.
/// @author James Thompson


// Unit headers
#include <core/scoring/methods/MembraneLipo.hh>
#include <core/scoring/methods/MembraneLipoCreator.hh>

// Package headers
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/Atom.hh>

#include <core/scoring/MembranePotential.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
//#include <core/scoring/rms_util.hh>

//#include <core/io/pdb/pose_io.hh>

//symmetry
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

//#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>

// Auto-header: duplicate removed #include <core/conformation/Residue.hh>

//#include <basic/prof.hh>
//#include <utility/exit.hh>

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the MembraneLipo class,
/// never an instance already in use
methods::EnergyMethodOP
MembraneLipoCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new MembraneLipo;
}

ScoreTypes
MembraneLipoCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( Mlipo );
	return sts;
}


/// c-tor
MembraneLipo::MembraneLipo() :
	parent( new MembraneLipoCreator ),
	potential_( ScoringManager::get_instance()->get_MembranePotential() )
{}


/// clone
EnergyMethodOP
MembraneLipo::clone() const
{
	return new MembraneLipo();
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/// @brief Calculate the RMS difference between native_pose_ (provided by
/// the option -in::file::native and the given Pose. The actual energy calculation
/// is the difference between the RMSD and the target RMSD. Target RMSD is specified
/// the option -score::rms_target.

///

void
MembraneLipo::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	// compute interpolated number of neighbors at various distance cutoffs
	pose.update_residue_neighbors();
	potential_.compute_centroid_environment( pose );
//	potential_.compute_membrane_embedding( pose );

}

void
MembraneLipo::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {

	//using namespace core::pose::datacache::CacheableDataType;

//	if(potential_.Menv_penalties()) //bw quick hack before putting them as individual scoring terms....
//	{
//		Real tm_projection(0);
//		Real hbond_pen(0);
//		Real termini_pen(0);
//		potential_.tm_projection_penalty(pose,tm_projection);
//		potential_.hbond_penalty(pose,hbond_pen);
//		potential_.termini_penalty(pose,termini_pen);
//		emap[ Menv_hbond ]=hbond_pen;
//		emap[ Menv_termini ]=termini_pen;

	Real lipo(0);

	MembraneTopology const & topology( *( static_cast< MembraneTopology const * >( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY )())));
	if(topology.LipoDefined()) {
		Real cen10Buried(0);
		Real cen10Exposed(0);
		Real cen10Buried_norm(0);
		Real cen10Exposed_norm(0);
		CenListInfo const & cenlist( *( static_cast< CenListInfo const * >( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::CEN_LIST_INFO )() ))); //cenlist_from_pose( pose ));
		for(Size i=1;i<=pose.total_residue();++i) {
			Size rsdSeq(i);
			if ( core::pose::symmetry::is_symmetric( pose ) ) {
				using namespace core::conformation::symmetry;
				SymmetricConformation const & symm_conf (
														 dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
				SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
				if (!symm_info->bb_is_independent(pose.residue(i).seqpos())) {
					rsdSeq = symm_info->bb_follows(pose.residue(i).seqpos());
				}
				if (symm_info->is_virtual(i)) {
					rsdSeq = 0;
				}
			}

			if (rsdSeq ==0 ) continue;
			if ( pose.residue(rsdSeq).aa() == core::chemical::aa_vrt ) continue;
			if(topology.allow_scoring(rsdSeq)) {
				Real B(topology.LipidBurial(rsdSeq));
				Real E(topology.LipidExposure(rsdSeq));
				if(B!=0) {
					cen10Buried+=B*cenlist.fcen10(i);
					cen10Buried_norm+=1;
				}
				if(E!=0) {
					cen10Exposed+=E*cenlist.fcen10(i);
					cen10Exposed_norm+=1;
				}

				//std::cout << i << " " << cenlist.fcen10(i) << " " << topology.tmh_inserted() << "\n";
			}
		}
		Real B_mean(0);
		Real E_mean(0);
		if(cen10Exposed_norm!=0)
			E_mean=cen10Exposed/cen10Exposed_norm;
		if(cen10Buried_norm!=0)
			B_mean=cen10Buried/cen10Buried_norm;
		lipo=(E_mean-B_mean)*topology.tmh_inserted();
	}
	emap[ Mlipo ]=lipo;
	potential_.finalize( pose );
//

//		std::cout << "Menv_penalties (tm_projection+hbond_pen+termini_pen+10) " << tm_projection << " " << hbond_pen << " " << termini_pen << std::endl;
//	}
//	potential_.finalize( pose );
//	totals[ rms ]  = std::abs( rms_target_ - rmsd );

	// PROF_STOP( basic::RMS );
}
/*
MembraneTopology const &
MembraneLipo::MembraneTopology_from_pose( pose::Pose & pose ) const
	{
		return *( static_cast< MembraneTopology const * >( pose.data().get_const_ptr( basic::MEMBRANE_TOPOLOGY )() ));

	}
core::Size
MembraneLipo::version() const
{
	return 1; // Initial versioning
}

	*/
core::Size
MembraneLipo::version() const
{
	return 1; // Initial versioning
}

} // namespace methods
} // namespace scoring
} // namespace core
