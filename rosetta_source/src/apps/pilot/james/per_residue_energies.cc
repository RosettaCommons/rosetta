// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file per_residue_energies.cc
/// @brief simple application for printing out energies from a Pose.
/// @author James Thompson

#include <devel/init.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/options/option.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <utility/vector1.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <string>

int
main( int argc, char* argv [] ) {
	devel::init( argc, argv );

	using core::Size;
	using core::Real;
	using core::pose::Pose;
	using utility::vector1;
	using namespace basic;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace ObjexxFCL;

	MetaPoseInputStream input = streams_from_cmd_line();
	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
		option[ in::file::residue_type_set ]()
	);

	// set up ScoreFunction, make certain that hydrogen bonding energies
	// are kept in the EnergyGraph.
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
	core::scoring::methods::EnergyMethodOptionsOP emopts(
		new core::scoring::methods::EnergyMethodOptions( scorefxn->energy_method_options() )
	);
	emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
	scorefxn->set_energy_method_options( *emopts );

	SilentFileData sfd;
	core::pose::Pose current_pose;
	while ( input.has_another_pose() ) {
		input.fill_pose( current_pose, *rsd_set );
		(*scorefxn)(current_pose);
		EnergyMap weights( current_pose.energies().weights() );

		for ( Size jj = 1; jj <= current_pose.total_residue(); ++jj ) {
			EnergyMap rsd_energies(
				weights * current_pose.energies().residue_total_energies(jj)
			);

			SilentStructOP ss( new ScoreFileSilentStruct );
			ss->decoy_tag( "residue_" + string_of(jj) );
			ss->add_string_value( "pose_id", core::pose::tag_from_pose(current_pose) );
			Real total(0);
			for ( int ii = 1; ii <= n_score_types; ++ii ) {
				if ( weights[ ScoreType(ii) ] != 0.0 ) {
					Real const value( rsd_energies[ ScoreType(ii) ] );
					std::string const scorename( name_from_score_type( ScoreType(ii) ) );
					total += value;
					ss->add_energy( scorename, value );
				}
			} // for n_score_types
			ss->add_energy( "score", total );
			sfd.write_silent_struct( *ss, option[ out::file::silent ]() );
		} // for current_pose.total_residue()
	} // while ( input.has_another_pose() )

	return 0;
}
