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
	asymm_only_=false;
	ignore_bs_=false;
	mask_=true;
	outmap_name_="";
	bin_squared_=true;
	b_sharpen_=0.0;
}

void ScaleMapIntensities::apply(core::pose::Pose & pose) {
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


	utility::vector1< core::Real > resobins;
	utility::vector1< core::Size > counts;
	utility::vector1< core::Real > rescale_factor(nresbins_,0.0);
	utility::vector1< core::Real > fade(nresbins_,1.0);
	utility::vector1< core::Real > mapI, modelI(nresbins_,0.0);

	core::scoring::electron_density::getDensityMap().getResolutionBins(nresbins_, 1.0/res_low_, 1.0/res_high_, resobins, counts, bin_squared_);

	if (b_sharpen_ == 0) {
		////
		//// use model
		////
		core::scoring::electron_density::getDensityMap().getIntensities( litePose, nresbins_, 1.0/res_low_, 1.0/res_high_, modelI, bin_squared_ );

		if (mask_) {
			mapI = core::scoring::electron_density::getDensityMap().getIntensitiesMasked( litePose, nresbins_, 1.0/res_low_, 1.0/res_high_, bin_squared_);
		} else {
			mapI = core::scoring::electron_density::getDensityMap().getIntensities( nresbins_, 1.0/res_low_, 1.0/res_high_, bin_squared_);
		}

		if (res_fade_ > 0 && res_high_ > res_fade_) {
			core::Real inv_fade_s = 0.5/(res_fade_) + 0.5/(res_high_);
			core::Real fade_width = 1/(res_high_) - 1/(res_fade_);
			core::Real sigma = 12.0/fade_width;
			for (Size i=1; i<=nresbins_; ++i) {
				fade[i] = 1/(1+exp(sigma*(resobins[i]-inv_fade_s)) );
			}
		}

		for (Size i=1; i<=nresbins_; ++i) {
			if (mapI[i] > 0 && fade[i]*modelI[i] >= 0 )
				rescale_factor[i] = fade[i]*sqrt(modelI[i] / mapI[i]);
		}
	} else {
		////
		//// bfactor sharpen
		////
 		mapI = core::scoring::electron_density::getDensityMap().getIntensities( nresbins_, 1.0/res_low_, 1.0/res_high_, bin_squared_);

		if (res_fade_ > 0 && res_high_ < res_fade_) {
			core::Real inv_fade_s = 0.5/(res_fade_) + 0.5/(res_high_);
			core::Real fade_width = 1/(res_high_) - 1/(res_fade_);
			core::Real sigma = 12.0/fade_width;
			for (Size i=1; i<=nresbins_; ++i) {
				fade[i] = 1/(1+exp(sigma*(resobins[i]-inv_fade_s)) );
			}
		}
		for (Size i=1; i<=nresbins_; ++i) {
			rescale_factor[i] = fade[i]*std::sqrt( exp(-b_sharpen_*resobins[i]*resobins[i]) );
		}
	}

	TR << "SCALING MAP:" << std::endl;
	TR << "resbin   fade   model   map   rescale" << std::endl;
	for (Size i=1; i<=nresbins_; ++i)
		TR << resobins[i] << "  " << fade[i] << "  " << modelI[i] << " " << mapI[i] << " " << rescale_factor[i] << std::endl;

	core::scoring::electron_density::getDensityMap().scaleIntensities( rescale_factor, 1.0/res_low_, 1.0/res_high_, bin_squared_ );
	if (outmap_name_.length()>0)
		core::scoring::electron_density::getDensityMap().writeMRC( outmap_name_ );
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
ScaleMapIntensities::parse_my_tag(
								   TagCOP const tag,
								   basic::datacache::DataMap & datamap,
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
	if ( tag->hasOption("b_sharpen") ) {
		b_sharpen_ = tag->getOption<core::Real>("b_sharpen");
	}
	if ( tag->hasOption("res_fade") ) {
		res_fade_ = tag->getOption<core::Real>("res_fade");
	}
	if ( tag->hasOption("mask") ) {
		mask_ = tag->getOption<bool>("mask");
	}
	if ( tag->hasOption("nresbins") ) {
		nresbins_ = tag->getOption<core::Size>("nresbins");
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
	if ( tag->hasOption("bin_squared") ) {
		bin_squared_ = tag->getOption<bool>("bin_squared");
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
