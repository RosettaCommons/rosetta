// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file  apps/pilot/rmoretti/residue_energy_breakdown.cc
/// @brief Output score terms for residue self energies and residue-residue interaction as a silent scorefile
/// @detailed Based on apps/pilot/james/per_residue_energies.cc by James Thompson
/// @author Rocco Moretti (rmoretti@uw.edu)

#include <devel/init.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/options/option.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>

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

	try {

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

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	core::scoring::methods::EnergyMethodOptionsOP emopts(
		new core::scoring::methods::EnergyMethodOptions( scorefxn->energy_method_options() )
	);
	emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
	scorefxn->set_energy_method_options( *emopts );
	// Post-loading commandline scorefunction alterations.
	core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *scorefxn );
	// TODO: Do we need other scorefunction modifications here?

	SilentFileData sfd;
	core::pose::Pose current_pose;
	ScoreTypes const scoretypes( scorefxn->get_nonzero_weighted_scoretypes() );
	EnergyMap weights( scorefxn->weights() );

	while ( input.has_another_pose() ) {
		input.fill_pose( current_pose, *rsd_set );

		// Load commandline pose modifications
		core::scoring::constraints::add_constraints_from_cmdline_to_pose( current_pose );
		// TODO: Others that should go here? Disulfides, rotamer bonuses, etc.?

		(*scorefxn)(current_pose);
		Energies const & pose_energies( current_pose.energies() );
		EnergyGraph const & egraph( pose_energies.energy_graph() );
		std::string tag( core::pose::tag_from_pose(current_pose) );

		//One body energies
		for ( Size ii = 1; ii <= current_pose.total_residue(); ++ii ) {
			EnergyNode const * node(egraph.get_energy_node(ii));
			runtime_assert( node != 0 );
			EnergyMap unwt_residue1b;
			scorefxn->eval_ci_1b(current_pose.residue(ii), current_pose, unwt_residue1b);
			scorefxn->eval_cd_1b(current_pose.residue(ii), current_pose, unwt_residue1b);
			scorefxn->eval_intrares_energy(current_pose.residue(ii), current_pose, unwt_residue1b);
			EnergyMap residue1b( unwt_residue1b * weights );

			SilentStructOP ss( new ScoreFileSilentStruct );
			ss->decoy_tag( tag + "_" + string_of(ii) + "_onebody" );
			ss->add_string_value( "pose_id", tag );
			ss->add_string_value( "resi1", string_of(ii) );
			ss->add_string_value( "restype1", current_pose.residue_type(ii).name3() );
			ss->add_string_value( "resi2", "--" );
			ss->add_string_value( "restype2", "onebody" );
			for ( ScoreTypes::const_iterator iter( scoretypes.begin() ); iter != scoretypes.end(); ++iter ) {
				ss->add_energy(  name_from_score_type( *iter ), residue1b[ *iter ] );
			} // for non-zero score types
			ss->add_energy( "total", residue1b.sum() );
			sfd.write_silent_struct( *ss, option[ out::file::silent ]() );
		}

		// Two body energies
		for ( Size ii = 1; ii <= current_pose.total_residue(); ++ii ) {
			for ( Size jj = ii + 1; jj <= current_pose.total_residue(); ++jj ) {
				EnergyMap pair_energies;
				bool output( false );

				// Short Range
				EnergyEdge const * edge(egraph.find_energy_edge(ii, jj));
				if ( edge != 0 ) {
					EnergyMap unwt_pair_energies( edge->fill_energy_map() );
					pair_energies += unwt_pair_energies * weights;
					output = true;
				}

				// Long Range
				for( Size lr = 1; lr <= core::scoring::methods::n_long_range_types; lr++){
					LREnergyContainerCOP lrec = pose_energies.long_range_container( core::scoring::methods::LongRangeEnergyType( lr ) );
					if( !lrec || lrec->empty()) { continue; }
					for ( ResidueNeighborConstIteratorOP rni( lrec->const_upper_neighbor_iterator_begin( ii ) ),
							end( lrec->const_upper_neighbor_iterator_end( ii ) ); *rni != *end; ++(*rni) ) {
						if( rni->upper_neighbor_id() != jj ) { continue; }
						rni->accumulate_energy( pair_energies );
						output = true;
					}
				}

				//Output
				if( !output ) { continue; }

				SilentStructOP ss( new ScoreFileSilentStruct );
				ss->decoy_tag( tag + "_" + string_of(ii) + "_" + string_of(jj) );
				ss->add_string_value( "pose_id", tag );
				ss->add_string_value( "resi1", string_of(ii) );
				ss->add_string_value( "restype1", current_pose.residue_type(ii).name3() );
				ss->add_string_value( "resi2", string_of(jj) );
				ss->add_string_value( "restype2", current_pose.residue_type(jj).name3() );

				bool nonzero(false);
				for ( ScoreTypes::const_iterator iter( scoretypes.begin() ); iter != scoretypes.end(); ++iter ) {
					ss->add_energy(  name_from_score_type( *iter ), pair_energies[ *iter ] );
					if ( pair_energies[ *iter ] <= -0.001 || 0.001 <= pair_energies[ *iter ] ) {
						nonzero = true;
					}
				} // for non-zero score types
				ss->add_energy( "total", pair_energies.sum() );
				if ( pair_energies.sum() <= -0.001 || 0.001 <= pair_energies.sum() ) {
					nonzero = true;
				}
				if (nonzero) { // Ignore pairs which would print as zero for all components
					sfd.write_silent_struct( *ss, option[ out::file::silent ]() );
				}
			} // for inner residue loop
		} // for Two body current_pose.total_residue()

	} // while ( input.has_another_pose() )
	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
