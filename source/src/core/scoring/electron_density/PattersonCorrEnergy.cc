// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/PattersonCorrEnergy.cc
/// @brief  Scoring a structure's fit to electron density
/// @author Frank DiMaio


// Unit headers
#include <core/scoring/electron_density/PattersonCorrEnergy.hh>
#include <core/scoring/electron_density/PattersonCorrEnergyCreator.hh>
#include <basic/options/option.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>

#include <core/scoring/Energies.hh>

#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/OneToAllEnergyContainer.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/statistics/functions.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>


// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <basic/options/keys/patterson.OptionKeys.gen.hh>

// Utility headers


#include <basic/Tracer.hh>

#include <core/chemical/AtomType.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <utility/vector1.hh>

#ifdef WIN32
	#define _USE_MATH_DEFINES
	#include <math.h>
#endif

using basic::T;
using basic::Error;
using basic::Warning;

// C++

namespace core {
namespace scoring {
namespace electron_density {


/// @details This must return a fresh instance of the PattersonCorrEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
PattersonCorrEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new PattersonCorrEnergy );
}

ScoreTypes
PattersonCorrEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( patterson_cc );
	return sts;
}

inline core::Real SQ( core::Real N ) { return N*N; }

using namespace core::scoring::methods;
static thread_local basic::Tracer TR( "core.scoring.electron_density.PattersonCorrEnergy" );

methods::LongRangeEnergyType
PattersonCorrEnergy::long_range_type() const { return patterson_corr_energy; }

/// c-tor
PattersonCorrEnergy::PattersonCorrEnergy() : parent( methods::EnergyMethodCreatorOP( new PattersonCorrEnergyCreator ) ) {
	map_loaded = core::scoring::electron_density::getDensityMap().isMapLoaded();  // loads map
	scoreRepacks = basic::options::option[ basic::options::OptionKeys::patterson::use_on_repack ]();
}


/// clone
EnergyMethodOP PattersonCorrEnergy::clone() const {
	return EnergyMethodOP( new PattersonCorrEnergy( *this ) );
}

/////////////////////////////////////////////////////////////////////////////

bool PattersonCorrEnergy::defines_residue_pair_energy(
	pose::Pose const & pose,
	Size res1,
	Size res2
) const {
	return ( pose.residue( res1 ).aa() == core::chemical::aa_vrt || pose.residue( res2 ).aa() == core::chemical::aa_vrt );
}


void
PattersonCorrEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & /* sf */) const {
	if (!pose.is_fullatom()) return;

	core::scoring::electron_density::getDensityMap().matchPoseToPatterson( pose, true );
}

void
PattersonCorrEnergy::finalize_after_derivatives( pose::Pose &, ScoreFunction const &  ) const {
}

void
PattersonCorrEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const &
) const {
	using namespace methods;

	// after packing the pose is rescored
	isRepacking = false;

	// Do we have a map?
	if (!map_loaded) {
		utility_exit_with_message("Density scoring function called but no map loaded.");
	}

	// make sure the root of the FoldTree is a virtual atom and is followed by a jump
	kinematics::Edge const &root_edge ( *pose.fold_tree().begin() );
	int virt_res_idx = root_edge.start();
	conformation::Residue const &root_res( pose.residue( virt_res_idx ) );

	// create LR energy container
	LongRangeEnergyType const & lr_type( long_range_type() );
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		OneToAllEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::OneToAllEnergyContainer > ( lrc ) );
		// make sure size or root did not change
		if ( dec->size() != pose.total_residue() || dec->fixed() != virt_res_idx ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		TR << "Creating new one-to-all energy container (" << pose.total_residue() << ")" << std::endl;
		LREnergyContainerOP new_dec( new OneToAllEnergyContainer( virt_res_idx, pose.total_residue(),  elec_dens_window ) );
		energies.set_long_range_container( lr_type, new_dec );
	}

	pose_is_proper = true;
	if (root_res.aa() != core::chemical::aa_vrt || root_edge.label() < 0) {
		pose_is_proper = false;  // we may be able to recover from this some time but for now just exit
		//utility_exit_with_message("Fold tree is not set properly for density scoring!");
		TR.Warning << "Fold tree is not set properly for patterson function scoring." << std::endl;
		pcc_structure = 0.0;
		return;
	}

	// allocate space for per-AA stats
	int nres = pose.total_residue();
	core::scoring::electron_density::getDensityMap().set_nres( nres );

	// do the actual matching here; split scores among individual residues
	pcc_structure = core::scoring::electron_density::getDensityMap().matchPoseToPatterson( pose, false );

	// # of NON-VRT residues
	nreses = 0;
	for (int i=1; i<=(int)pose.total_residue(); ++i)
		if (pose.residue(i).aa() != core::chemical::aa_vrt)
			nreses++;

	TR.Debug << "PattersonCorrEnergy::setup_for_scoring() returns PCC = " << pcc_structure << std::endl;
}


void
PattersonCorrEnergy::finalize_total_energy(
	pose::Pose const & ,
	ScoreFunction const &,
	EnergyMap &
) const {
	return;
}


void
PattersonCorrEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & /* pose */,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	if (!pose_is_proper) return;

	using namespace numeric::statistics;

	if (rsd1.aa() != core::chemical::aa_vrt && rsd2.aa() != core::chemical::aa_vrt) return;
	if (rsd1.aa() == core::chemical::aa_vrt && rsd2.aa() == core::chemical::aa_vrt) return;

	conformation::Residue const & rsd = (rsd1.aa() == core::chemical::aa_vrt) ? rsd2 : rsd1;

 	if ( isRepacking && scoreRepacks) {
 		// update patterson map quickly only moving sidechain of 'rsd'
 		core::Real cc = core::scoring::electron_density::getDensityMap().rematchResToPatterson( rsd );
 		emap[ patterson_cc ] += -cc;
 		TR.Debug << "Rescore residue " << rsd.seqpos() << " :  patterson CC = " << cc << std::endl;
 	} else {
		// just return score from setup_for_scoring
		emap[ patterson_cc ] += -pcc_structure;
	}

	return;
}

void
PattersonCorrEnergy::update_residue_for_packing(
	pose::Pose &pose,
	Size resid
) const {
	if (!pose_is_proper) return;

	core::scoring::electron_density::getDensityMap().updateCachedDensity( pose.residue(resid) );
	//std::cout << "UPDATE residue " << resid << " " << std::endl;
}


void
PattersonCorrEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & ,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	if (!pose_is_proper) return;
	if (!pose.is_fullatom()) return;

	using namespace numeric::statistics;

	int resid = id.rsd();
	int atmid = id.atomno();

	// if (hydrogen) return
	if ( pose.residue(resid).aa() != core::chemical::aa_vrt && !pose.residue(resid).atom_type(atmid).is_heavyatom() ) return;

	numeric::xyzVector<core::Real> X = pose.xyz(id);
	numeric::xyzVector< core::Real > dCCdx;
	numeric::xyzMatrix< core::Real > R = numeric::xyzMatrix<core::Real>::rows(1,0,0, 0,1,0, 0,0,1);

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		; // TO DO
	} else {
		// ASYMMETRIC CASE

		core::scoring::electron_density::getDensityMap().dCCdx_pat( atmid, resid, X, pose, dCCdx );
		numeric::xyzVector< core::Real > dEdx = -(core::Real)nreses*dCCdx;

		numeric::xyzVector<core::Real> atom_x = X;
		numeric::xyzVector<core::Real> const f2( dEdx );
		numeric::xyzVector<core::Real> atom_y = -f2 + atom_x;
		Vector const f1( atom_x.cross( atom_y ) );

		F1 += weights[ patterson_cc ] * f1;
		F2 += weights[ patterson_cc ] * f2;
	}
}
core::Size
PattersonCorrEnergy::version() const
{
	return 1; // Initial versioning
}


}
}
}
