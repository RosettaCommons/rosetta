// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Set up morphing with electron density map
/// @author Yifan Song

#include <protocols/electron_density/VoxelSpacingRefinementMover.hh>
#include <protocols/electron_density/VoxelSpacingRefinementMoverCreator.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/electron_density/xray_scattering.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/scoring/cryst/util.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/kinematics/MoveMap.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/cryst.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>


#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.electron_density.VoxelSpacingRefinementMover" );

namespace protocols {
namespace electron_density {

using core::scoring::electron_density::poseCoords;
using core::scoring::electron_density::poseCoord;

VoxelSpacingMultifunc::VoxelSpacingMultifunc( core::pose::Pose const &pose)
{
	// get rhoC from pose
	ObjexxFCL::FArray3D< double > rhoMask;

	// convert to poseCoords
	poseCoords litePose;
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		core::conformation::Residue const & rsd_i ( pose.residue(i) );
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
		core::Size natoms = rsd_i.nheavyatoms();
		for ( core::Size j = 1; j <= natoms; ++j ) {
			core::chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );
			poseCoord coord_j;
			coord_j.x_ = rsd_i.xyz( j );
			coord_j.B_ = pose.pdb_info()->temperature( i, j );
			coord_j.elt_ = atom_type_set[ rsd_i.atom_type_index( j ) ].element();
			litePose.push_back( coord_j );
		}
	}

	// get rhoC
	core::scoring::electron_density::getDensityMap().calcRhoC( litePose, 0.0, rhoC_, rhoMask, -1, 1e6 );
}


core::Real
VoxelSpacingMultifunc::operator ()( core::optimization::Multivec const & vars ) const {
	numeric::xyzVector< core::Real > ori ( vars[1],vars[2],vars[3] );
	numeric::xyzVector< core::Real > apix_scale;

	if ( vars.size() == 4 ) {
		apix_scale[0] = apix_scale[1] = apix_scale[2] = vars[4];
	} else {
		apix_scale[0] = vars[4];
		apix_scale[1] = vars[5];
		apix_scale[2] = vars[6];
	}

	// get CC between rhoC and map
	core::Real sumC=0, sumC2=0, sumO=0, sumO2=0, sumCO=0;
	core::Real N = rhoC_.u1()*rhoC_.u2()*rhoC_.u3();

	for ( int k=1; k<=rhoC_.u3(); ++k ) {
		for ( int j=1; j<=rhoC_.u2(); ++j ) {
			for ( int i=1; i<=rhoC_.u1(); ++i ) {
				core::Real rhoC_ijk = rhoC_(i,j,k);

				numeric::xyzVector< core::Real > mapidx(
					ori[0] + 1./apix_scale[0]*(i-1) +1,
					ori[1] + 1./apix_scale[1]*(j-1) +1,
					ori[2] + 1./apix_scale[2]*(k-1) +1
				);

				core::Real rhoO_ijk = core::scoring::electron_density::getDensityMap().get( mapidx );

				sumC  += rhoC_ijk;
				sumC2 += rhoC_ijk*rhoC_ijk;
				sumO  += rhoO_ijk;
				sumO2 += rhoO_ijk*rhoO_ijk;
				sumCO += rhoC_ijk*rhoO_ijk;
			}
		}
	}

	core::Real varC = (sumC2 - sumC*sumC / N );
	core::Real varO = (sumO2 - sumO*sumO / N ) ;
	runtime_assert ( varC != 0 && varO != 0 );
	core::Real CC = (sumCO - sumC*sumO/N) / sqrt( varC * varO );

	return -100*CC;
}

void
VoxelSpacingMultifunc::dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const {
	// analytic
	dE_dvars.clear( );
	dE_dvars.resize( vars.size(), 0.0 );
	numeric::xyzVector< core::Real > ori ( vars[1],vars[2],vars[3] );
	numeric::xyzVector< core::Real > apix_scale;

	if ( vars.size() == 4 ) {
		apix_scale[0] = apix_scale[1] = apix_scale[2] = vars[4];
	} else {
		apix_scale[0] = vars[4];
		apix_scale[1] = vars[5];
		apix_scale[2] = vars[6];
	}


	// get CC between rhoC and map
	core::Real sumC=0, sumC2=0, sumO=0, sumO2=0, sumCO=0;
	core::optimization::Multivec dsumCO(6,0), dsumO(6,0), dsumO2(6,0);
	core::Real N = rhoC_.u1()*rhoC_.u2()*rhoC_.u3();

	for ( int k=1; k<=rhoC_.u3(); ++k ) {
		for ( int j=1; j<=rhoC_.u2(); ++j ) {
			for ( int i=1; i<=rhoC_.u1(); ++i ) {
				core::Real rhoC_ijk = rhoC_(i,j,k);

				numeric::xyzVector< core::Real > offset(
					(i-1)/apix_scale[0] +1,
					(j-1)/apix_scale[1] +1,
					(k-1)/apix_scale[2] +1
				);

				numeric::xyzVector< core::Real > mapidx(
					ori[0] + offset[0],
					ori[1] + offset[1],
					ori[2] + offset[2]
				);

				core::Real rhoO_ijk = core::scoring::electron_density::getDensityMap().get( mapidx );

				sumC  += rhoC_ijk;
				sumC2 += rhoC_ijk*rhoC_ijk;
				sumO  += rhoO_ijk;
				sumO2 += rhoO_ijk*rhoO_ijk;
				sumCO += rhoC_ijk*rhoO_ijk;

				numeric::xyzVector< core::Real > drhoO_ijk = core::scoring::electron_density::getDensityMap().grad( mapidx );
				dsumCO[1] += rhoC_ijk * drhoO_ijk[0];
				dsumCO[2] += rhoC_ijk * drhoO_ijk[1];
				dsumCO[3] += rhoC_ijk * drhoO_ijk[2];
				dsumCO[4] += rhoC_ijk * drhoO_ijk[0] * (-(i-1) / (apix_scale[0]*apix_scale[0]));
				dsumCO[5] += rhoC_ijk * drhoO_ijk[1] * (-(j-1) / (apix_scale[1]*apix_scale[1]));
				dsumCO[6] += rhoC_ijk * drhoO_ijk[2] * (-(k-1) / (apix_scale[2]*apix_scale[2]));

				dsumO[1] += drhoO_ijk[0];
				dsumO[2] += drhoO_ijk[1];
				dsumO[3] += drhoO_ijk[2];
				dsumO[4] += drhoO_ijk[0] * (-(i-1) / (apix_scale[0]*apix_scale[0]));
				dsumO[5] += drhoO_ijk[1] * (-(j-1) / (apix_scale[1]*apix_scale[1]));
				dsumO[6] += drhoO_ijk[2] * (-(k-1) / (apix_scale[2]*apix_scale[2]));

				dsumO2[1] += 2 * rhoO_ijk * drhoO_ijk[0];
				dsumO2[2] += 2 * rhoO_ijk * drhoO_ijk[1];
				dsumO2[3] += 2 * rhoO_ijk * drhoO_ijk[2];
				dsumO2[4] += 2 * rhoO_ijk * drhoO_ijk[0] * (-(i-1) / (apix_scale[0]*apix_scale[0]));
				dsumO2[5] += 2 * rhoO_ijk * drhoO_ijk[1] * (-(j-1) / (apix_scale[1]*apix_scale[1]));
				dsumO2[6] += 2 * rhoO_ijk * drhoO_ijk[2] * (-(k-1) / (apix_scale[2]*apix_scale[2]));
			}
		}
	}

	core::Real varC = sqrt(sumC2 - sumC*sumC / N );
	core::Real varO = sqrt(sumO2 - sumO*sumO / N ) ;
	runtime_assert ( varC != 0 && varO != 0 );
	core::Real CC = (sumCO - sumC*sumO/N) / ( varC * varO );

	for (core::Size i=1; i<=6; ++i) {
		core::Size i_eff = std::min( i, vars.size() );
		core::Real dvarO_i = (dsumO2[i] - 2*sumO*dsumO[i]/N) / (2*varO);
		dE_dvars[i_eff] += -100.0 * ( (dsumCO[i] - sumC*dsumO[i]/N) / ( varC * varO ) - CC*dvarO_i / varO );
	}

	std::cerr << "CC(" << vars[1];
	for (core::Size i=2; i<=vars.size(); ++i) {
		std::cerr << "," << vars[i];
	}
	std::cerr << ") = " << CC << std::endl;
}

void
VoxelSpacingMultifunc::dump( core::optimization::Multivec const &, core::optimization::Multivec const & ) const {
}

void
VoxelSpacingMultifunc::getMapSpacingAndOrigin( core::optimization::Multivec & vars, bool aniso ) {
	// construct the starting point for refinement
	vars.clear();

	// origin
	vars.push_back(0);
	vars.push_back(0);
	vars.push_back(0);

	// scale
	vars.push_back(1);
	if ( aniso ) {
		vars.push_back(1);
		vars.push_back(1);
	}
}

void
VoxelSpacingMultifunc::foldInChanges( core::pose::Pose &pose, core::optimization::Multivec & vars ) {
	// read origin & apix_scale parameters, apply to pose
	numeric::xyzVector< core::Real > apix_scale(1.,1.,1.);
	if ( vars.size() == 4 ) {
		apix_scale[0] = apix_scale[1] = apix_scale[2] = vars[4];
	} else {
		apix_scale[0] = vars[4];
		apix_scale[1] = vars[5];
		apix_scale[2] = vars[6];
	}

	numeric::xyzVector<int> grid = core::scoring::electron_density::getDensityMap().getGrid();
	numeric::xyzMatrix<core::Real> f2c = core::scoring::electron_density::getDensityMap().get_f2c();
	numeric::xyzVector<core::Real> origin = core::scoring::electron_density::getDensityMap().getOrigin(), origin_shift;

	origin_shift = f2c * numeric::xyzVector<core::Real>(
		(vars[1] + (origin[0]+1)*(1.-1./apix_scale[0]))/grid[0],
		(vars[2] + (origin[1]+1)*(1.-1./apix_scale[1]))/grid[1],
		(vars[3] + (origin[2]+1)*(1.-1./apix_scale[2]))/grid[2]
	);

	pose.apply_transform_Rx_plus_v( numeric::xyzMatrix<core::Real>::identity(), origin_shift );

	numeric::xyzVector< core::Real > new_apix =
		core::scoring::electron_density::getDensityMap().get_voxel_spacing();

	new_apix[0] *= apix_scale[0];
	new_apix[1] *= apix_scale[1];
	new_apix[2] *= apix_scale[2];

	core::scoring::electron_density::getDensityMap().set_voxel_spacing( new_apix );
}

//////////////

VoxelSpacingRefinementMover::VoxelSpacingRefinementMover() : moves::Mover() {
	init();
}

void VoxelSpacingRefinementMover::init() {
	aniso_ = false;
	minimizer_ = "lbfgs_armijo";
	max_iter_ = 100;
}

void VoxelSpacingRefinementMover::apply(core::pose::Pose & pose) {
	core::optimization::MinimizerOptions options( minimizer_, 1e-2, true, false, false );
	options.max_iter(max_iter_);

	// set up optimizer
	VoxelSpacingMultifunc f_voxel( pose );
	core::optimization::Minimizer minimizer( f_voxel, options );

	// ... and run
	core::optimization::Multivec y;
	f_voxel.getMapSpacingAndOrigin( y, aniso_ );
	minimizer.run( y );

	f_voxel.foldInChanges( pose, y );

	if ( mapout_.length() > 0 ) {
		core::scoring::electron_density::getDensityMap().writeMRC(mapout_);
	}
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
VoxelSpacingRefinementMover::parse_my_tag(TagCOP const tag, basic::datacache::DataMap &, Filters_map const &,
	moves::Movers_map const &, Pose const &)
{
	if ( tag->hasOption("aniso") ) {
		aniso_ = tag->getOption<bool>("aniso");
	}

	if ( tag->hasOption("mapout") ) {
		mapout_ = tag->getOption<std::string>("mapout");
	}

	if ( tag->hasOption("max_iter") ) {
		max_iter_ = tag->getOption<core::Size>("max_iter");
	}

}

protocols::moves::MoverOP
VoxelSpacingRefinementMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new VoxelSpacingRefinementMover );
}

std::string
VoxelSpacingRefinementMoverCreator::keyname() const
{
	return VoxelSpacingRefinementMoverCreator::mover_name();
}

std::string
VoxelSpacingRefinementMoverCreator::mover_name()
{
	return "VoxelSpacingRefinement";
}

}
}
