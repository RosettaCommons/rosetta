// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/BurialEnergy.cc
/// @author James Thompson

#include <core/scoring/methods/BurialEnergy.hh>
#include <core/scoring/methods/BurialEnergyCreator.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

#include <core/pose/Pose.hh>

#include <basic/prof.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

/// @details This must return a fresh instance of the BurialEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
BurialEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new BurialEnergy );
}

ScoreTypes
BurialEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back(burial);
	return sts;
}

void
BurialEnergy::setup_for_scoring(
	pose::Pose & pose, ScoreFunction const &
) const {
	pose.update_residue_neighbors();
}

/// clone
EnergyMethodOP
BurialEnergy::clone() const {
	return EnergyMethodOP( new BurialEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
BurialEnergy::residue_energy(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose,
	EnergyMap & emap
) const {
	static const core::Real dist_sq_cutoff(100);
	static const core::Size NN_cutoff(16);

	TwelveANeighborGraph const & graph ( pose.energies().twelveA_neighbor_graph() );

	core::conformation::Atom const & atom_i = rsd.atom(rsd.nbr_atom());

	// iterate across neighbors within 12 angstroms, count number <10A
	Size countN(0);
	for ( graph::Graph::EdgeListConstIter
			ir  = graph.get_node( rsd.seqpos() )->const_edge_list_begin(),
			ire = graph.get_node( rsd.seqpos() )->const_edge_list_end();
			ir != ire; ++ir ) {
		Size const j( (*ir)->get_other_ind( rsd.seqpos() ) );
		core::conformation::Residue const & rsd_j( pose.residue(j) );

		core::conformation::Atom const & atom_j = rsd_j.atom(rsd_j.nbr_atom());
		Real sqdist = atom_i.xyz().distance_squared( atom_j.xyz() );
		if ( sqdist < dist_sq_cutoff ) {
			countN++;
		}
	}
	//std::cout << "countN(" << ii << ") = " << countN << std::endl;

	Real const & burial_prediction( pred_burial_[rsd.seqpos()] );
	Real score = burial_prediction;
	if ( countN < NN_cutoff ) {
		emap[ burial ] += -1 * score;
	} else {
		emap[ burial ] += score;
	}
}

core::Size
BurialEnergy::version() const
{
	return 1;
}

void
BurialEnergy::indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
{
	context_graphs_required[ twelve_A_neighbor_graph ] = true;
}

void BurialEnergy::init_from_file() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string const & burial_fn( option[ in::file::burial ]()[1] );
	utility::io::izstream input(burial_fn);
	if ( !input.good() ) {
		std::string const msg( "Error opening file: " + burial_fn );
		utility_exit_with_message( msg );
	}

	using std::string;
	string line;
	getline(input,line); // header
	while ( getline(input,line) ) {
		std::istringstream ss(line);
		core::Size resi;
		char aa;
		ss >> resi >> aa;
		string last_token;
		while ( ss.good() ) ss >> last_token;
		std::istringstream dbl_reader( last_token.substr(2,4) );
		Real burial(0);
		dbl_reader >> burial;
		if ( last_token.substr(0,1) == "0" ) burial *= -1;
		pred_burial_.push_back( burial );
		//std::cout << "line = " << line << std::endl;
		//std::cout << "burial = " << burial << std::endl;
	}
}

} // methods
} // scoring
} // core
