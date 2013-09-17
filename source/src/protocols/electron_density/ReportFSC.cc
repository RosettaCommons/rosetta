// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file

#include <protocols/electron_density/ReportFSC.hh>
#include <protocols/electron_density/ReportFSCCreator.hh>

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
#include <basic/options/keys/edensity.OptionKeys.gen.hh>


#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>

static basic::Tracer TR( "protocols.electron_density.ReportFSC" );

namespace protocols {
namespace electron_density {

using core::scoring::electron_density::poseCoords;
using core::scoring::electron_density::poseCoord;


ReportFSC::ReportFSC() : moves::Mover() {
	init();
}

void ReportFSC::init() {
	res_low_=1000.0;
	res_high_=0.0;
	nresbins_=50;
	asymm_only_=false;
	testmap_ = NULL;
}

void ReportFSC::apply(core::pose::Pose & pose) {
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
			coord_j.B_ = pose.pdb_info()->temperature( i, j );
			coord_j.elt_ = atom_type_set[ rsd_i.atom_type_index( j ) ].element();

			litePose.push_back( coord_j );
		}
	}

	utility::vector1< core::Real > modelmap1FSC(nresbins_,1.0), modelmap2FSC(nresbins_,1.0);
	core::Real fsc1=0.0, fsc2=0.0;

	// train map
	modelmap1FSC = core::scoring::electron_density::getDensityMap().getFSC( litePose, nresbins_, 1.0/res_low_, 1.0/res_high_ );

	for (Size i=1; i<=modelmap1FSC.size(); ++i) fsc1+=modelmap1FSC[i];
	fsc1 /= modelmap1FSC.size();

	if (testmap_ && testmap_->isMapLoaded()) {
		modelmap2FSC = testmap_->getFSC( litePose, nresbins_, 1.0/res_low_, 1.0/res_high_ );
		for (Size i=1; i<=modelmap2FSC.size(); ++i) fsc2+=modelmap2FSC[i];
		fsc2 /= modelmap2FSC.size();
	}

	// tag
	core::pose::RemarkInfo remark;
	std::ostringstream oss;
	oss << "FSC(" << res_low_ << ":" << res_high_ << ") = " << fsc1;
	if (testmap_ && testmap_->isMapLoaded())
		oss << " / " << fsc2;

	TR << oss.str() << std::endl;
	remark.num = 1;	remark.value = oss.str();
	pose.pdb_info()->remarks().push_back( remark );
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
ReportFSC::parse_my_tag(
								   TagPtr const tag,
								   moves::DataMap & datamap,
								   Filters_map const & filters,
								   moves::Movers_map const & movers,
								   Pose const & pose
							   )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( tag->hasOption("res_low") ) {
		res_low_ = tag->getOption<core::Real>("res_low");
	}
	if ( tag->hasOption("res_high") ) {
		res_high_ = tag->getOption<core::Real>("res_high");
	}
	if ( tag->hasOption("nresbins") ) {
		nresbins_ = tag->getOption<core::Size>("nresbins");
	}
	if ( tag->hasOption("asymm_only") ) {
		asymm_only_ = tag->getOption<bool>("asymm_only");
	}
	if ( tag->hasOption("testmap") ) {
		std::string mapfile = tag->getOption<std::string>("testmap");
		core::Real mapreso = option[ edensity::mapreso ]();
		core::Real mapsampling = option[ edensity::grid_spacing ]();
		std::cerr << "Loading alternate density map " << mapfile << std::endl;
		testmap_ = new core::scoring::electron_density::ElectronDensity();
		testmap_->readMRCandResize( mapfile , mapreso , mapsampling );
	}
}

protocols::moves::MoverOP
ReportFSCCreator::create_mover() const {
	return new ReportFSC;
}

std::string
ReportFSCCreator::keyname() const
{
	return ReportFSCCreator::mover_name();
}

std::string
ReportFSCCreator::mover_name()
{
	return "ReportFSC";
}

}
}
