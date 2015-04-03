// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SequenceDependentRefEnergy.hh
/// @brief  Reference energy method implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/scoring/methods/SequenceDependentRefEnergy.hh>
#include <core/scoring/methods/SequenceDependentRefEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>

// C++ Headers
#include <string>
#include <vector>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the SequenceDependentRefEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
SequenceDependentRefEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new SequenceDependentRefEnergy );
}

ScoreTypes
SequenceDependentRefEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( seqdep_ref );
	return sts;
}


SequenceDependentRefEnergy::SequenceDependentRefEnergy() :
	parent( methods::EnergyMethodCreatorOP( new SequenceDependentRefEnergyCreator ) )
{
	read_energy_weight_table();
}

SequenceDependentRefEnergy::SequenceDependentRefEnergy( utility::vector1< utility::vector1< Real > > const & aa_seq_weights_in ):
	parent( methods::EnergyMethodCreatorOP( new SequenceDependentRefEnergyCreator ) )
{
	aa_seq_weights_.clear();
	for ( utility::vector1<utility::vector1< Real > >::const_iterator it = aa_seq_weights_in.begin(); it != aa_seq_weights_in.end(); ++it ) {
			aa_seq_weights_.push_back(*it);
	}
}

SequenceDependentRefEnergy::~SequenceDependentRefEnergy() {}


void SequenceDependentRefEnergy::read_energy_weight_table() {

	using namespace basic::options;

	std::cout << "JL checking for SequenceDependentRefEnergy weights" << std::endl;

	if ( option[ OptionKeys::score::seqdep_refene_fname ].user() ) {
		std::string const in_fname( option[ OptionKeys::score::seqdep_refene_fname ] );
		std::cout << "JL reading SequenceDependentRefEnergy weights from " << in_fname << std::endl;
		utility::io::izstream in_stream( in_fname );
		if (!in_stream.good()) {
			utility_exit_with_message( "[ERROR] Error opening SequenceDependentRefEnergy file" );
		}
		std::string line;
		core::Size seqpos(0);
		while( getline( in_stream, line) ) {
			++seqpos;
			std::cout << "JL got " << seqpos << " line " << line << std::endl;
			utility::vector1< std::string > const tokens ( utility::split( line ) );
			utility::vector1< Real > energies;
			for ( utility::vector1< std::string >::const_iterator it = tokens.begin(); it != tokens.end(); ++it ) {
				energies.push_back( atof(it->c_str()) );
			}
			aa_seq_weights_.push_back(energies);
		}

		if ( option[ OptionKeys::score::secondary_seqdep_refene_fname ].user() ) {
			std::string const in_fname( option[ OptionKeys::score::secondary_seqdep_refene_fname ] );
			std::cout << "JL reading SECONDARY SequenceDependentRefEnergy weights from " << in_fname << std::endl;
			utility::io::izstream in_stream( in_fname );
			if (!in_stream.good()) {
				utility_exit_with_message( "[ERROR] Error opening SECONDARY SequenceDependentRefEnergy file" );
			}
			std::string line;
			core::Size seqpos(0);
			while( getline( in_stream, line) ) {
				++seqpos;
				std::cout << "JL got " << seqpos << " line " << line << std::endl;
				utility::vector1< std::string > const tokens ( utility::split( line ) );
				utility::vector1< Real > energies;
				for ( utility::vector1< std::string >::const_iterator it = tokens.begin(); it != tokens.end(); ++it ) {
					energies.push_back( atof(it->c_str()) );
				}
				core::Size aa(0);
				for ( utility::vector1< Real >::const_iterator it = energies.begin(); it != energies.end(); ++it ) {
					++aa;
					(aa_seq_weights_[seqpos])[ aa ] += *it;
				}
			}
		}

	}

}


EnergyMethodOP
SequenceDependentRefEnergy::clone() const
{
	return EnergyMethodOP( new SequenceDependentRefEnergy( aa_seq_weights_ ) );
}


void
SequenceDependentRefEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const
{
	using namespace chemical;

	if ( aa_seq_weights_.empty() ) return;

	AA const & aa( rsd.aa() );
	Size const & seqpos( rsd.seqpos() );
	if ( seqpos > aa_seq_weights_.size() ) return;

	Real const ene = (aa_seq_weights_[seqpos])[ aa ];
	emap[ seqdep_ref ] += ene;
	//	std::cout << "JL at seqpos " << seqpos << " for aa " << aa << " using seqdep_ref " <<  ene << std::endl;

	return;

}


Real
SequenceDependentRefEnergy::eval_dof_derivative(
	id::DOF_ID const &,
	id::TorsionID const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap const &
) const
{
	return 0.0;
}

/// @brief SequenceDependentRefEnergy is context independent; indicates that no
/// context graphs are required
void
SequenceDependentRefEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}
core::Size
SequenceDependentRefEnergy::version() const
{
	return 1; // Initial versioning
}
} // methods
} // scoring
} // core

