// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/HRF_MSLabelingEnergy.cc
/// @author Melanie Aprahamian (aprahamian.4@osu.edu)

#include <core/scoring/methods/HRF_MSLabelingEnergy.hh>
#include <core/scoring/methods/HRF_MSLabelingEnergyCreator.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/pose/Pose.hh>
#include <basic/prof.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/func/FadeFunc.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>

#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/vector1.hh>
#include <numeric/NumericTraits.hh>

namespace core {
namespace scoring {
namespace methods {

/// @details This must return a fresh instance of the HRF_MSLabelingEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
HRF_MSLabelingEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new HRF_MSLabelingEnergy );
}

ScoreTypes
HRF_MSLabelingEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back(hrf_ms_labeling);
	return sts;
}

void
HRF_MSLabelingEnergy::setup_for_scoring(
	pose::Pose & pose, ScoreFunction const &
) const {
	pose.update_residue_neighbors();

}

// constructor
HRF_MSLabelingEnergy::HRF_MSLabelingEnergy() :
	methods::ContextDependentOneBodyEnergy( methods::EnergyMethodCreatorOP( new HRF_MSLabelingEnergyCreator ) )
{
	init_from_file();
	dist_midpoint_ = basic::options::option[ basic::options::OptionKeys::score::ms_dist_midpoint ]();
	dist_exponent_ = basic::options::option[ basic::options::OptionKeys::score::ms_dist_exponent ]();
	slope_ = basic::options::option[ basic::options::OptionKeys::score::ms_fit_slope ]();
	intercept_ = basic::options::option[ basic::options::OptionKeys::score::ms_fit_intercept ]();
	fade_outer_ = basic::options::option[ basic::options::OptionKeys::score::ms_fade_outer ]();
	fade_dist_ = basic::options::option[ basic::options::OptionKeys::score::ms_fade_dist ]();
}

/// clone
EnergyMethodOP
HRF_MSLabelingEnergy::clone() const {
	return EnergyMethodOP( new HRF_MSLabelingEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
HRF_MSLabelingEnergy::residue_energy(
	core::conformation::Residue const & residue,
	core::pose::Pose const & pose,
	EnergyMap & emap
) const {

	core::Size const number_residues (pose.total_residue());
	core::Real neighbor_count (0.0); // calculated using a logistic function, so non-integer
	core::Real lnPF_from_file (0.0);
	core::Real distance (0.0);

	for ( core::Size j=1; j <= prot_factor_.size(); j++ ) {
		if ( prot_factor_[j].first == residue.seqpos() ) {
			lnPF_from_file = prot_factor_[j].second;
			for ( core::Size res_count_neighbor = 1; res_count_neighbor <= number_residues; res_count_neighbor++ ) {
				if ( residue.seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
					distance = residue.xyz("CEN").distance(pose.residue(res_count_neighbor).xyz("CEN"));
					neighbor_count += 1.0/(1.0 + std::exp(dist_exponent_*(distance-dist_midpoint_)));
				}
			}
			func::FuncOP hrf_ms_labeling_score_ (func::FuncOP(new func::FadeFunc(-fade_outer_, fade_outer_, fade_dist_, -1.0)));
			core::Real const number_of_neighbors_difference (std::abs(neighbor_count - (slope_*lnPF_from_file + intercept_)));
			core::Real const burial_score = hrf_ms_labeling_score_->func(number_of_neighbors_difference);
			emap[hrf_ms_labeling] += burial_score;
		}
	}
}

void
HRF_MSLabelingEnergy::indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
{
	context_graphs_required[ twelve_A_neighbor_graph ] = true;
}


core::Size
HRF_MSLabelingEnergy::version() const
{
	return 1;
}

void HRF_MSLabelingEnergy::init_from_file() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string const & hrf_ms_labeling_fn( option[ in::file::hrf_ms_labeling ]());

	utility::io::izstream input(hrf_ms_labeling_fn);
	if ( !input.good() ) {
		std::string const msg( "Error opening file: " + hrf_ms_labeling_fn );
		utility_exit_with_message( msg );
	}

	std::string line;
	getline(input,line);

	while ( getline(input,line) ) {
		std::istringstream ss(line);
		core::Size resi;
		core::Real lnPF;

		ss >> resi >> lnPF;

		prot_factor_.push_back( std::pair< core::Size, core::Real >(resi, lnPF) );
	}

}

} // methods
} // scoring
} // core
