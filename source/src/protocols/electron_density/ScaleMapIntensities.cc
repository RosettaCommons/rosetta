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

static THREAD_LOCAL basic::Tracer TR( "protocols.electron_density.ScaleMapIntensities" );

namespace protocols {
namespace electron_density {

using core::scoring::electron_density::poseCoords;
using core::scoring::electron_density::poseCoord;


ScaleMapIntensities::ScaleMapIntensities() : moves::Mover() {
	init();
}

void ScaleMapIntensities::init() {
	res_low_=10000.0;
	res_high_= 0;
	fade_width_ = 0.1;
	nresbins_=50;
	asymm_only_=false;
	ignore_bs_=false;
	mask_=true;
	mask_output_=false;
	outmap_name_="";
	bin_squared_=true;
	b_sharpen_=0.0;
	truncate_only_=false;
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

	if ( b_sharpen_ == 0 ) {
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			if ( asymm_only_ && symm_info && !symm_info->bb_is_independent( i ) ) continue;
			core::conformation::Residue const & rsd_i ( pose.residue(i) );
			if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;

			core::Size natoms = rsd_i.nheavyatoms();
			for ( core::Size j = 1; j <= natoms; ++j ) {
				core::chemical::AtomTypeSet const & atom_type_set( rsd_i.atom_type_set() );

				poseCoord coord_j;
				coord_j.x_ = rsd_i.xyz( j );
				if ( ignore_bs_ ) {
					coord_j.B_ = 0.0;
				} else {
					coord_j.B_ = pose.pdb_info()->temperature( i, j );
				}
				coord_j.elt_ = atom_type_set[ rsd_i.atom_type_index( j ) ].element();

				litePose.push_back( coord_j );
			}
		}
	}

	utility::vector1< core::Real > resobins;
	utility::vector1< core::Size > counts;
	utility::vector1< core::Real > rescale_factor(nresbins_,0.0);
	utility::vector1< core::Real > fade(nresbins_,1.0);
	utility::vector1< core::Real > mapI, modelI(nresbins_,0.0);

	core::scoring::electron_density::getDensityMap().getResolutionBins(nresbins_, 0.0, 0.0, resobins, counts, bin_squared_);

	if ( b_sharpen_ == 0 ) {
		////
		//// use model
		////
		ObjexxFCL::FArray3D< double > rhoC, rhoMask;
		ObjexxFCL::FArray3D< std::complex<double> > Frho;
		core::scoring::electron_density::getDensityMap().calcRhoC( litePose, res_high_, rhoC, rhoMask );  // truncate mask at highres limit
		numeric::fourier::fft3(rhoC, Frho);

		if ( outmap_name_.length()>0 ) {
			core::scoring::electron_density::ElectronDensity tmp = core::scoring::electron_density::getDensityMap();
			tmp.set_data( rhoC );
			tmp.writeMRC( outmap_name_+"_rho_calc.mrc" );
			tmp.set_data( rhoMask );
			tmp.writeMRC( outmap_name_+"_mask_calc.mrc" );
		}
		core::scoring::electron_density::getDensityMap().getIntensities( Frho, nresbins_, 0.0, 0.0, modelI, bin_squared_ );

		// reuse rhoC
		if ( mask_ ) {
			for ( int i=0; i<rhoC.u1()*rhoC.u2()*rhoC.u3(); ++i ) {
				rhoC[i] = rhoMask[i] * core::scoring::electron_density::getDensityMap().get_data()[i];
			}
			if ( mask_output_ ) core::scoring::electron_density::getDensityMap().set_data( rhoC );
		} else {
			for ( int i=0; i<rhoC.u1()*rhoC.u2()*rhoC.u3(); ++i ) rhoC[i] = core::scoring::electron_density::getDensityMap().get_data()[i];
		}
		numeric::fourier::fft3(rhoC, Frho);
		core::scoring::electron_density::getDensityMap().getIntensities( Frho, nresbins_, 0.0, 0.0, mapI, bin_squared_ );

		for ( Size i=1; i<=nresbins_; ++i ) {
			if ( mapI[i] > 0 ) {
				rescale_factor[i] = sqrt(modelI[i] / mapI[i]);
			}
		}
	} else if ( !truncate_only_ ) {
		////
		//// bfactor sharpen
		////
		ObjexxFCL::FArray3D< std::complex<double> > Frho;
		numeric::fourier::fft3(core::scoring::electron_density::getDensityMap().get_data(), Frho);
		core::scoring::electron_density::getDensityMap().getIntensities( Frho, nresbins_, 0.0, 0.0, mapI, bin_squared_);

		for ( Size i=1; i<=nresbins_; ++i ) {
			rescale_factor[i] = std::sqrt( exp(-b_sharpen_*resobins[i]*resobins[i]/4.0) );
		}
	}

	if ( b_sharpen_ == 0 ) {
		TR << "SCALING MAP:" << std::endl;
		TR << "resbin   model   map   rescale" << std::endl;
		for ( Size i=1; i<=nresbins_; ++i ) {
			TR << resobins[i] << "  " << modelI[i] << " " << mapI[i] << " " << rescale_factor[i] << std::endl;
		}
	} else {
		TR << "SCALING MAP:" << std::endl;
		TR << "resbin    map   rescale" << std::endl;
		for ( Size i=1; i<=nresbins_; ++i ) {
			TR << resobins[i] << "  " << mapI[i] << " " << rescale_factor[i] << std::endl;
		}
	}

	core::scoring::electron_density::getDensityMap().scaleIntensities( rescale_factor, 0.0, 0.0, bin_squared_ );
	core::scoring::electron_density::getDensityMap().reciprocalSpaceFilter( res_low_, res_high_, fade_width_);

	if ( outmap_name_.length()>0 ) {
		core::scoring::electron_density::getDensityMap().writeMRC( outmap_name_+"_scaled.mrc" );
	}
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
ScaleMapIntensities::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &)
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
	if ( tag->hasOption("truncate_only") ) {
		truncate_only_ = tag->getOption<core::Real>("truncate_only");
	}
	if ( tag->hasOption("fade_width") ) {
		fade_width_ = tag->getOption<core::Real>("fade_width");
	}
	if ( tag->hasOption("mask") ) {
		mask_ = tag->getOption<bool>("mask");
	}
	if ( tag->hasOption("mask_output") ) {
		mask_output_ = tag->getOption<bool>("mask_output");
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
	return protocols::moves::MoverOP( new ScaleMapIntensities );
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
