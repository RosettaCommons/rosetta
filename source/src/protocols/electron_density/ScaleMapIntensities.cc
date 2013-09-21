// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file

#include <protocols/electron_density/ScaleMapIntensities.hh>
#include <protocols/electron_density/ScaleMapIntensitiesCreator.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
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


#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>

static basic::Tracer TR( "protocols.electron_density.ScaleMapIntensities" );

namespace protocols {
namespace electron_density {

using core::scoring::electron_density::poseCoords;
using core::scoring::electron_density::poseCoord;


ScaleMapIntensities::ScaleMapIntensities() : moves::Mover() {
	init();
}

void ScaleMapIntensities::init() {
	res_low_=1000.0;
	res_high_= res_fade_= 0;
	nresbins_=50;
	scale_by_fsc_=false;
	asymm_only_=false;
	ignore_bs_=false;
	outmap_name_="";
}

void ScaleMapIntensities::apply(core::pose::Pose & pose) {
	utility::vector1< core::Real > mapI, modelI;

	// pose->poseCoords
	poseCoords litePose;
	core::conformation::symmetry::SymmetryInfoOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}
	for (int i=1; i<=pose.total_residue(); ++i) {
		if (asymm_only_ && symm_info && !symm_info->bb_is_independent( i ) ) continue;
		core::conformation::Residue const & rsd_i ( pose.residue(i) );
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;

		core::Size natoms = rsd_i.nheavyatoms();
		for (int j=1; j<=natoms; ++j) {
			core::conformation::Atom const &atom_j( rsd_i.atom(j) );
			core::chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );

			poseCoord coord_j;
			coord_j.x_ = rsd_i.xyz( j );
			if (ignore_bs_) {
				coord_j.B_ = 0.0;
			} else {
				coord_j.B_ = pose.pdb_info()->temperature( i, j );
			}
			coord_j.elt_ = atom_type_set[ rsd_i.atom_type_index( j ) ].element();

			litePose.push_back( coord_j );
		}
	}


	utility::vector1< core::Real > resobins = core::scoring::electron_density::getDensityMap().getResolutionBins(nresbins_, 1.0/res_low_, 1.0/res_high_);
	core::Real resobins_step = resobins[2]-resobins[1];

	// TODO: make masking optional
	core::scoring::electron_density::getDensityMap().getIntensities( litePose, nresbins_, 1.0/res_low_, 1.0/res_high_, modelI );
	mapI = core::scoring::electron_density::getDensityMap().getIntensitiesMasked( litePose, nresbins_, 1.0/res_low_, 1.0/res_high_);

	utility::vector1< core::Real > fade(nresbins_,1.0);
	if (scale_by_fsc_) {
		fade = core::scoring::electron_density::getDensityMap().getFSCMasked( litePose, nresbins_, 1.0/res_low_, 1.0/res_high_ );
	} else if (res_fade_ > 0) {
		// fade over ~10 buckets
		core::Real sigma = (1/resobins_step);
		core::Real inv_fade = (1.0/res_fade_);
		for (Size i=1; i<=nresbins_; ++i) {
			fade[i] = 1/(1+exp(sigma*(resobins[i]-inv_fade)) );
		}
	}

	utility::vector1< core::Real > rescale_factor(nresbins_,0.0);
	for (Size i=1; i<=nresbins_; ++i) {
		if (mapI[i] > 0 && fade[i]*modelI[i] >= 0 )
			rescale_factor[i] = sqrt(fade[i]*modelI[i] / mapI[i]);
	}

	TR << "SCALING MAP:" << std::endl;
	TR << "resbin   fade   model   map   rescale" << std::endl;
	for (Size i=1; i<=nresbins_; ++i)
		TR << resobins[i] << "  " << fade[i] << "  " << modelI[i] << " " << mapI[i] << " " << rescale_factor[i] << std::endl;

	core::scoring::electron_density::getDensityMap().scaleIntensities( rescale_factor, 1.0/res_low_, 1.0/res_high_ );
	if (outmap_name_.length()>0)
		core::scoring::electron_density::getDensityMap().writeMRC( outmap_name_ );
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
ScaleMapIntensities::parse_my_tag(
								   TagPtr const tag,
								   moves::DataMap & datamap,
								   Filters_map const & filters,
								   moves::Movers_map const & movers,
								   Pose const & pose
							   )
{
	if ( tag->hasOption("res_low") ) {
		res_low_ = tag->getOption<core::Real>("res_low");
	}
	if ( tag->hasOption("res_high") ) {
		res_high_ = tag->getOption<core::Real>("res_high");
	}
	if ( tag->hasOption("res_fade") ) {
		res_fade_ = tag->getOption<core::Real>("res_fade");
	}
	if ( tag->hasOption("nresbins") ) {
		nresbins_ = tag->getOption<core::Size>("nresbins");
	}
	if ( tag->hasOption("scale_by_fsc") ) {
		scale_by_fsc_ = tag->getOption<bool>("scale_by_fsc");
	}
	if ( tag->hasOption("asymm_only") ) {
		asymm_only_ = tag->getOption<bool>("asymm_only");
	}
	if ( tag->hasOption("ignore_bs") ) {
		ignore_bs_ = tag->getOption<bool>("ignore_bs");
	}
	if ( tag->hasOption("outmap") ) {
		outmap_name_ = tag->getOption<std::string>("outmap");
	}
}

protocols::moves::MoverOP
ScaleMapIntensitiesCreator::create_mover() const {
	return new ScaleMapIntensities;
}

std::string
ScaleMapIntensitiesCreator::keyname() const
{
	return ScaleMapIntensitiesCreator::mover_name();
}

std::string
ScaleMapIntensitiesCreator::mover_name()
{
	return "ScaleMapIntensities";
}

}
}
