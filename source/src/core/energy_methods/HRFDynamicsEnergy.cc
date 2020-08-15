// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/energy_methods/HRFDynamicsEnergy.cc
/// @brief  Score term calculated based upon neighbor count agreement with hydroxyl radical footprinting labeling data. In a manuscript to be published in 2020.
/// @author Sarah Biehn (biehn.4@osu.edu)

#include <core/energy_methods/HRFDynamicsEnergy.hh>
#include <core/energy_methods/HRFDynamicsEnergyCreator.hh>
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
#include <core/scoring/methods/EnergyMethodOptions.hh>

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

/// @details This must return a fresh instance of the HRFDynamicsEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
HRFDynamicsEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< HRFDynamicsEnergy >( options );
}

ScoreTypes
HRFDynamicsEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back(hrf_dynamics);
	return sts;
}

void
HRFDynamicsEnergy::setup_for_scoring(
	pose::Pose & pose, ScoreFunction const &
) const {
	pose.update_residue_neighbors();

}

// constructor
HRFDynamicsEnergy::HRFDynamicsEnergy( methods::EnergyMethodOptions const & options ) :
	parent( utility::pointer::make_shared< HRFDynamicsEnergyCreator >() ),
	hrf_dynamics_input_file_( options.hrf_dynamics_input() )
{
	init_from_file();
}

/// clone
EnergyMethodOP
HRFDynamicsEnergy::clone() const {
	return utility::pointer::make_shared< HRFDynamicsEnergy >( *this );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
HRFDynamicsEnergy::residue_energy(
	core::conformation::Residue const & residue,
	core::pose::Pose const & pose,
	EnergyMap & emap
) const {

	core::Size const number_residues (pose.total_residue());
	core::Real neighbor_count (0.0);
	core::Real pf_from_file (0.0);
	core::Real nc_from_file (0.0);
	core::Real distance (0.0);
	core::Real angle (0.0);


	for ( core::Size j=1; j <= input_nc_.size(); j++ ) {
		if ( input_nc_[j].first == residue.seqpos() ) {
			pf_from_file = input_nc_[j].second;
			nc_from_file = (pf_from_file*0.88)+4.42; //0.88 is intercept; 4.42 is slope. These values were obtained from lnPF and neighbor count relationship as described in Biehn 2020 manuscript. LnPF was plotted against the averaged neighbor count from 30 mover models generated with nma_relax.xml. The resulting linear fit of the data, ie the prediction equation, is used here to calculated the predicted neighbor count based on the values in the protection factor (mass spec data) file
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
					neighbor_count += 1.0/(1.0 + std::exp(1.0*(distance-9.0)))*1.0/(1.0 + std::exp(numeric::NumericTraits< float >::pi()*2.0*(angle-numeric::NumericTraits< float >::pi()/2.0)));
				}
			}
			emap[hrf_dynamics] += -1.0/(1.0 + std::exp(2.0*(std::abs(neighbor_count - nc_from_file)-3.5))); //3.5 was also optimized in Biehn 2020 manuscript
		}
	}
}

void
HRFDynamicsEnergy::indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
{
	context_graphs_required[ twelve_A_neighbor_graph ] = false;
}


core::Size
HRFDynamicsEnergy::version() const
{
	return 1;
}

void HRFDynamicsEnergy::init_from_file() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::io::izstream input(hrf_dynamics_input_file_);
	if ( !input.good() ) {
		std::string const msg( "Error opening file: " + hrf_dynamics_input_file_ );
		utility_exit_with_message( msg );
	}

	std::string line;

	while ( getline(input,line) ) {
		if ( line.substr(0,1) == "#" ) continue;
		std::istringstream ss(line);
		core::Size resi;
		core::Real nc;

		ss >> resi >> nc;

		runtime_assert_string_msg( !(ss.fail() || ss.bad()), "Error in HRFDynamicsEnergy::init_from_file():  Could not parse line \"" + line + "\"." );
		input_nc_.push_back( std::pair< core::Size, core::Real >(resi, nc) );
	}

}

} // methods
} // scoring
} // core
