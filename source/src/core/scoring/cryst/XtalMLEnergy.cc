// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/cryst/XtalMLEnergy.cc
/// @brief  ML target
/// @author Frank DiMaio


#include <core/scoring/cryst/XtalMLEnergy.hh>
#include <core/scoring/cryst/PhenixInterface.hh>

// Unit headers
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/methods/EnergyMethod.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>

#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/pdb/file_data.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/chemical/ChemicalManager.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>


// Utility headers
#include <basic/prof.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

// C++
#include <iomanip>
#include <string>

namespace core {
namespace scoring {
namespace cryst {

// tracer
static basic::Tracer TR("core.scoring.cryst.XtalMLEnergy");

ScoreTypes XtalMLEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( xtal_ml );
	sts.push_back( xtal_rwork );
	sts.push_back( xtal_rfree );
	return sts;
}

methods::EnergyMethodOP XtalMLEnergyCreator::create_energy_method( methods::EnergyMethodOptions const &) const {
	return new XtalMLEnergy();
}

XtalMLEnergy::XtalMLEnergy() :
	parent( new XtalMLEnergyCreator ) {
	ml=0;
	dml_dx.clear();
}

/// clone
methods::EnergyMethodOP XtalMLEnergy::clone() const {
	return new XtalMLEnergy( *this );
}

void XtalMLEnergy::setup_for_scoring( pose::Pose & , ScoreFunction const & ) const {

}

void XtalMLEnergy::finalize_total_energy(pose::Pose & pose, ScoreFunction const &sf, EnergyMap & emap) const {
	using namespace core::scoring;
	if (sf.get_weight( core::scoring::xtal_ml ) > 1e-8 || sf.get_weight( core::scoring::xtal_rwork ) > 1e-8 || sf.get_weight( core::scoring::xtal_rfree) > 1e-8) {
		emap[ xtal_ml ] += getPhenixInterface().getScore( pose );
		emap[ xtal_rwork ] += getPhenixInterface().getR( );
		emap[ xtal_rfree ] += getPhenixInterface().getRfree( );
	}
}

void
XtalMLEnergy::setup_for_minimizing( pose::Pose & pose, ScoreFunction const & sf, kinematics::MinimizerMapBase const & ) const {
	using namespace core::chemical;
	getPhenixInterface().getScore( pose );
	getPhenixInterface().updateSolventMask();
}


void XtalMLEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sf) const {
	using namespace core::chemical;

	if (sf.get_weight( core::scoring::xtal_rwork ) > 1e-8)
		utility_exit_with_message( "xtal_rwork  is not minimizable!" );
	if (sf.get_weight( core::scoring::xtal_rfree ) > 1e-8)
		utility_exit_with_message( "xtal_rfree is not minimizable!" );

	//fpd update mask
	getPhenixInterface().updateSolventMask( pose );
	if (sf.get_weight( core::scoring::xtal_ml ) > 1e-8 )
		ml = getPhenixInterface().getScoreAndDerivs( pose, dml_dx );
	else
		ml = 0;
}

void XtalMLEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & ,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	int resid = id.rsd();
	int atmid = id.atomno();
	core::conformation::Residue const &rsd_i = pose.residue(resid);

	if (rsd_i.aa() == core::chemical::aa_vrt) return;

	if (dml_dx.size() == 0) return; //

	// look up derivative
	numeric::xyzVector< core::Real > dCCdx = dml_dx[resid][atmid];
	numeric::xyzVector<core::Real> X = pose.xyz(id);

	// NUMERIC DERIV CHECK
	// const core::Real NUM_DERIV_H_CEN=0.1;
	// core::pose::Pose pose_copy = pose;
	// numeric::xyzVector<core::Real> dCCdxn;
	//
	// pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0]+NUM_DERIV_H_CEN,X[1],X[2] ) );
	// core::Real CC_px = getPhenixInterface().getScore( pose_copy );
	// pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0]-NUM_DERIV_H_CEN,X[1],X[2] ) );
	// core::Real CC_mx = getPhenixInterface().getScore( pose_copy );
	// pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1]+NUM_DERIV_H_CEN,X[2] ) );
	// core::Real CC_py = getPhenixInterface().getScore( pose_copy );
	// pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1]-NUM_DERIV_H_CEN,X[2] ) );
	// core::Real CC_my = getPhenixInterface().getScore( pose_copy );
	// pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1],X[2]+NUM_DERIV_H_CEN ) );
	// core::Real CC_pz = getPhenixInterface().getScore( pose_copy );
	// pose_copy.set_xyz(id, numeric::xyzVector<core::Real>( X[0],X[1],X[2]-NUM_DERIV_H_CEN ) );
	// core::Real CC_mz = getPhenixInterface().getScore( pose_copy );
	//
	// dCCdxn[0] = (CC_px-CC_mx)/(2*NUM_DERIV_H_CEN);
	// dCCdxn[1] = (CC_py-CC_my)/(2*NUM_DERIV_H_CEN);
	// dCCdxn[2] = (CC_pz-CC_mz)/(2*NUM_DERIV_H_CEN);
	//
	// TR << dCCdx[0] << "," << dCCdx[1] << "," << dCCdx[2] << "  ,  " << dCCdxn[0] << "," << dCCdxn[1] << "," << dCCdxn[2] << std::endl;

	numeric::xyzVector<core::Real> atom_x = X;
	numeric::xyzVector<core::Real> const f2( dCCdx );
	numeric::xyzVector<core::Real> atom_y = -f2 + atom_x;
	Vector const f1( atom_x.cross( atom_y ) );

	F1 += weights[ xtal_ml ] * f1;
	F2 += weights[ xtal_ml ] * f2;
}


core::Size
XtalMLEnergy::version() const
{
	return 1; // Initial versioning
}


} // cryst
} // scoring
} // core

