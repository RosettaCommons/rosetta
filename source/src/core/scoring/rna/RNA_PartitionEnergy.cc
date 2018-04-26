// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RNA_PartitionEnergy.cc
/// @brief  Compute RNA score for secondary structure partition function. Requires Vienna's RNAfold.
/// @author Ramya Rangan

// Basic headers
#include <basic/Tracer.hh>

// Unit headers
#include <core/scoring/rna/RNA_PartitionEnergy.hh>
#include <core/scoring/rna/RNA_PartitionEnergyCreator.hh>

// Package headers

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>

// Utility headers
#include <core/scoring/EnergyMap.hh>


// C++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>

#ifdef WIN32
#define popen  _popen
#define pclose _pclose
#endif

static basic::Tracer TR( "core.scoring.rna.RNA_PartitionEnergy" );

namespace core {
namespace scoring {
namespace rna {


/// @details This must return a fresh instance of the RNA_PartitionEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_PartitionEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new RNA_PartitionEnergy );
}

ScoreTypes
RNA_PartitionEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_partition );
	return sts;
}


/// c-tor // TODO initialize partition cache here
RNA_PartitionEnergy::RNA_PartitionEnergy() :
	parent1( methods::EnergyMethodCreatorOP( new RNA_PartitionEnergyCreator ) ),
	parent2()
{}

/// copy
RNA_PartitionEnergy::RNA_PartitionEnergy(
	RNA_PartitionEnergy const & src) :
	parent1( methods::EnergyMethodCreatorOP( new RNA_PartitionEnergyCreator ) ),
	parent2( src ),
	partition_cache_( src.partition_cache_ ),
	global_mapping_( src.global_mapping_ ),
	res_list_( src.global_mapping_ ),
	native_sequence_( src.native_sequence_ )
{}

/// clone
methods::EnergyMethodOP
RNA_PartitionEnergy::clone() const
{
	return methods::EnergyMethodOP( new RNA_PartitionEnergy( *this ) );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
RNA_PartitionEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	using namespace core::pose::full_model_info;

	std::string native_sequence = const_full_model_info( pose ).global_sequence();
	if ( native_sequence.compare("") == 0 ) {
		TR.Warning << "Using partition score but no global sequence provided." << std::endl;
	}
	std::string current_sequence = get_current_global_sequence(pose);
	Real score;
	get_score_from_sequences(native_sequence, current_sequence, score);
	totals[ rna_partition ] = score;

}

core::Real
RNA_PartitionEnergy::calculate_energy( utility::vector1< core::conformation::ResidueCOP > const & resvect,
	core::Size const //substitution_position = 0
) const {
	using namespace core::pose::full_model_info;

	std::string current_sequence = get_current_global_sequence( resvect, global_mapping_, res_list_, native_sequence_);
	Real score;
	get_score_from_sequences(native_sequence_, current_sequence, score);
	return score;
}

void
RNA_PartitionEnergy::set_up_residuearrayannealableenergy_for_packing (
	pose::Pose & pose,
	core::pack::rotamer_set::RotamerSets const &/*rotamersets*/,
	scoring::ScoreFunction const & //sfxn
) {
	using namespace core::pose::full_model_info;


	global_mapping_ = const_full_model_info( pose ).global_mapping();
	res_list_ = const_full_model_info( pose ).res_list();
	native_sequence_ = const_full_model_info( pose ).global_sequence();

}

void
RNA_PartitionEnergy::get_score_from_sequences(
	std::string const & native_sequence,
	std::string const & current_sequence,
	Real & score ) const {
	Real native_score;
	get_partition_score_from_cache( native_sequence, native_score );
	get_partition_score_from_cache( current_sequence, score );
	score = score - native_score;
}

void
RNA_PartitionEnergy::get_partition_score_from_cache(
	std::string const & sequence,
	Real & partition_score
) const {
	auto it = partition_cache_.find(sequence);
	if ( it == partition_cache_.end() ) {
		get_partition_score(sequence, partition_score);
		partition_cache_[ sequence ] = partition_score;
	} else {
		partition_score = it->second;
	}
}

void
RNA_PartitionEnergy::exec_cmd(
	std::string const & cmd,
	std::string & cmd_result
) const {
	std::array<char, 128> buffer;
	std::shared_ptr<FILE> pipe( popen( cmd.c_str(), "r" ), pclose );
	if ( !pipe ) {
		TR.Error << "Partition function shell command failed to run!" << std::endl;
	}
	while ( !feof(pipe.get()) ) {
		if ( fgets(buffer.data(), 128, pipe.get()) != nullptr ) {
			cmd_result += buffer.data();
		}
	}
}


void
RNA_PartitionEnergy::get_partition_score(
	std::string const & global_sequence,
	Real & partition_score
) const {
	// Run Vienna
	std::string command_line("echo " + global_sequence);
	command_line = command_line + " | RNAfold -p";

	// Parse Vienna output

	// get line with FE ensemble from commandline output
	command_line = command_line + " | tail -n 3 | head -n 1";
	// retrieve FE ensemble in kcal/mol from Vienna's output line
	command_line = command_line + " | sed 's/[][]//g' | awk '{ print $2 }'";
	try {
		std::string cmd_result;
		exec_cmd(command_line , cmd_result);
		partition_score = stod( cmd_result );
	} catch ( std::invalid_argument const & ) {
		TR << "Command did not produce a real number:" << std::endl;
		TR << command_line << std::endl;
		try { // Try again once more.
			std::string cmd_result;
			exec_cmd(command_line , cmd_result);
			partition_score = std::stod( cmd_result );
		} catch ( std::invalid_argument const & ) {
			TR << "Command did not produce a real number again." << std::endl;
			partition_score = 0;
		}
	}
	partition_score = -partition_score;
}


core::Size
RNA_PartitionEnergy::version() const
{
	return 1; // Initial versioning
}


} //rna
} //scoring
} //core
