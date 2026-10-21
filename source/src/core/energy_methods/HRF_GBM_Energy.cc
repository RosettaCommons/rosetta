// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/energy_methods/HRF_GBM_Energy.cc
/// @brief  Score term assesses neighbor count agreement between decoy and expected value generated from external lightGBM model. Model uses experimental HRF data as input feature. Detailed in manuscript published in 2025.
/// @author Elijah Day (ehday@ucla.edu)

#include <core/energy_methods/HRF_GBM_Energy.hh>
#include <core/energy_methods/HRF_GBM_EnergyCreator.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>


#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <numeric/NumericTraits.hh>

#include <basic/options/keys/OptionKeys.hh> // AUTO IWYU For OptionKeys,

namespace core {
namespace energy_methods {


/// @details Returns fresh instance of the HRF_GBM_Energy class.
core::scoring::methods::EnergyMethodOP
HRF_GBM_EnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< HRF_GBM_Energy >( options );
}

core::scoring::ScoreTypes
HRF_GBM_EnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back(hrf_gbm);
	return sts;
}

void
HRF_GBM_Energy::setup_for_scoring(
	pose::Pose & pose, core::scoring::ScoreFunction const &
) const {
	pose.update_residue_neighbors();

}

// constructor
HRF_GBM_Energy::HRF_GBM_Energy( core::scoring::methods::EnergyMethodOptions const & options ) :
	parent( utility::pointer::make_shared< HRF_GBM_EnergyCreator >() ),
	hrf_gbm_input_file_( options.hrf_gbm_input() )
{
	init_from_file();
}

/// clone
core::scoring::methods::EnergyMethodOP
HRF_GBM_Energy::clone() const {
	return utility::pointer::make_shared< HRF_GBM_Energy >( *this );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
HRF_GBM_Energy::residue_energy(
	core::conformation::Residue const & residue,
	core::pose::Pose const & pose,
	core::scoring::EnergyMap & emap
) const {

	core::Size const number_residues (pose.total_residue());
	core::Real pose_neighbor_count (0.0);
	core::Real nc_from_file (0.0);
	core::Real distance (0.0);
	core::Real angle (0.0);


	for ( core::Size j=1; j <= input_nc_.size(); j++ ) {
		if ( input_nc_[j].first == residue.seqpos() ) {
			nc_from_file = input_nc_[j].second; // Predicted neighbor count values were obtained from the external lightGBM model described in the Day 2025 manuscript. Briefly, the model takes in sequence derived features (PSSM profile, predicted SSE, etc.) and experimental mass spec data to generated an expected neighbor count for the residue. Those neighbor counts are passed as an input for the hrf_gbm score term.
			std::string target_atom ("CB");
			if ( residue.type().name1() == 'G' ) {
				target_atom = "1HA";
			}
			for ( core::Size res_count_neighbor = 1; res_count_neighbor <= number_residues; res_count_neighbor++ ) {
				if ( residue.seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
					std::string neighbor_atom ("CB");
					if ( pose.residue(res_count_neighbor).type().name1() == 'G' ) {
						neighbor_atom = "1HA";
					}
					distance = residue.xyz("CA").distance(pose.residue(res_count_neighbor).xyz(neighbor_atom));
					numeric::xyzVector<core::Real> norm_CA_target_vector = (residue.xyz(target_atom) - residue.xyz("CA")).normalize()/(residue.xyz(target_atom).distance(residue.xyz("CA")));
					numeric::xyzVector<core::Real> norm_neighbor_vector = (pose.residue(res_count_neighbor).xyz(neighbor_atom)-residue.xyz("CA"))/(pose.residue(res_count_neighbor).xyz(neighbor_atom).distance(residue.xyz("CA")));
					angle = std::acos(norm_CA_target_vector.dot(norm_neighbor_vector));
					pose_neighbor_count += 1.0/(1.0 + std::exp(1.0*(distance-9.0)))*1.0/(1.0 + std::exp(numeric::NumericTraits< float >::pi()*2.0*(angle-numeric::NumericTraits< float >::pi()/2.0)));
				}
			}
			emap[ core::scoring::hrf_gbm] += -1.0/(1.0 + std::exp(2.9*(std::abs(pose_neighbor_count - nc_from_file)-2.1))); //2.9 and 2.1 were optimized in the Day 2025 manuscript
		}
	}
}

void
HRF_GBM_Energy::indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
{
	context_graphs_required[ core::scoring::twelve_A_neighbor_graph ] = false;
}


core::Size
HRF_GBM_Energy::version() const
{
	return 1;
}

void HRF_GBM_Energy::init_from_file() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::io::izstream input(hrf_gbm_input_file_);
	if ( !input.good() ) {
		std::string const msg( "Error opening file: " + hrf_gbm_input_file_ );
		utility_exit_with_message( msg );
	}

	std::string line;

	while ( getline(input,line) ) {
		if ( line.substr(0,1) == "#" ) continue;
		std::istringstream ss(line);
		core::Size resi;
		core::Real nc;

		ss >> resi >> nc;

		runtime_assert_string_msg( !(ss.fail() || ss.bad()), "Error in HRF_GBM_Energy::init_from_file():  Could not parse line \"" + line + "\"." );
		input_nc_.push_back( std::pair< core::Size, core::Real >(resi, nc) );
	}

}

} // scoring
} // core
