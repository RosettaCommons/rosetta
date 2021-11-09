// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/public/analysis/atom_energy_breakdown.cc
/// @brief Output score terms for atom self energies and atom-atom interactions as a silent scorefile
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <devel/init.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/id/AtomID.hh>

#include <core/select/residue_selector/OrResidueSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

#include <core/scoring/electron_density/util.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethod.hh>

#include <core/scoring/util.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <protocols/electron_density/SetupForDensityScoringMover.hh>

#include <utility/vector1.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

#include <string>

OPT_KEY( StringVector, residue_selection )

std::string
pdb_id(core::pose::PDBInfo const & pdbinfo, core::Size ii ) {
	std::string pdbid( std::to_string( pdbinfo.number(ii) ) + pdbinfo.chain(ii) );
	if ( pdbinfo.icode(ii) != ' ' ) {
		pdbid += ':';
		pdbid += pdbinfo.icode(ii);
	}
	return pdbid;
}

int
main( int argc, char* argv [] ) {

	try {

		NEW_OPT( residue_selection, "Which residues to look at when outputting. Accepts pose numbered, pdb numbered and ranges.", utility::vector1<std::string>{} );

		devel::init( argc, argv );

		using namespace core::io::silent;
		using namespace core::scoring;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		core::scoring::methods::EnergyMethodOptionsOP emopts( new core::scoring::methods::EnergyMethodOptions( scorefxn->energy_method_options() ) );
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		scorefxn->set_energy_method_options( *emopts );
		// Post-loading commandline scorefunction alterations.
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *scorefxn );
		if ( option[ edensity::mapfile ].user() ) {
			core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn );
		}
		// TODO: Do we need other scorefunction modifications here?

		core::scoring::ScoreTypes supported_scoretypes; // Which types (concievably) have atomistic support.
		core::scoring::ScoreTypes active_scoretypes; // intersection between non_zero_scoretypes & supported_scoretypes

		// Prune to the scoretypes which are relevant
		for ( auto iter = scorefxn->all_energies_begin(); iter != scorefxn->all_energies_end(); ++iter ) {
			if ( (*iter)->has_atomistic_energies() || (*iter)->has_atomistic_pairwise_energies() ) {
				supported_scoretypes.append( (*iter)->score_types() );
			}
		}
		for ( ScoreType st: scorefxn->get_nonzero_weighted_scoretypes() ) {
			if ( supported_scoretypes.contains(st) ) {
				active_scoretypes.push_back(st);
			}
		}

		///////////////////////////////////////////
		// Pose processing

		core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
		core::chemical::ResidueTypeSetCOP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
			option[ in::file::residue_type_set ]()
		);

		core::io::silent::SilentFileOptions opts; // initialized from the command line
		core::io::silent::SilentFileData sfd(opts);

		while ( input.has_another_pose() ) {
			core::pose::Pose pose;
			input.fill_pose( pose, *rsd_set );

			// Load commandline pose modifications
			core::scoring::constraints::add_constraints_from_cmdline_to_pose( pose );
			if ( option[ edensity::mapfile ].user() ) {
				protocols::electron_density::SetupForDensityScoringMover().apply( pose );
			}

			if ( ! pose.pdb_info() ) {
				// Make a default PDB info if one doesn't already exist - this simplifies logic later.
				core::pose::PDBInfoOP new_pdb_info( new core::pose::PDBInfo(pose) );
				pose.pdb_info( new_pdb_info );
			}

			// TODO: Others that should go here? Disulfides, rotamer bonuses, etc.?

			std::string tag( core::pose::tag_from_pose(pose) );

			core::select::residue_selector::ResidueSubsetOP subset = nullptr;
			// Handle residue selection, if any.
			if ( option[ residue_selection ].user() ) {
				core::select::residue_selector::OrResidueSelector selector;
				for (  std::string const & entry: option[ residue_selection ] ) {
					if ( entry.size() == 1 && std::isalpha( entry[0] ) ) {
						// Special case for chain designations.
						selector.add_residue_selector( utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( entry[0] ) );
					} else {
						selector.add_residue_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( entry ) );
					}
				}
				subset = utility::pointer::make_shared< core::select::residue_selector::ResidueSubset >( selector.apply(pose) );
			}

			// Do a base scoring to make sure that the pose is properly set up for scoring.
			// We don't actually use the results of the scoring here --
			// we just want to make sure that all the associated data is properly set up.
			// (This is also needed for the two-body energy check.)
			(*scorefxn)(pose);

			utility::vector1< core::Real > totals( active_scoretypes.size(), 0.0 );

			// Onebody
			for ( auto const & entry: core::scoring::get_single_atom_energies(pose, *scorefxn, active_scoretypes, subset ) ) {
				core::Size seqpos = entry.first.rsd();
				core::Size atm = entry.first.atomno();
				utility::vector1< core::Real > const & onebody = entry.second;

				core::conformation::Residue const & res( pose.residue(seqpos) );

				std::string atom_name = utility::strip( res.atom_name(atm) );

				SilentStructOP ss( new ScoreFileSilentStruct(opts) );
				ss->decoy_tag( tag + "_" + std::to_string(seqpos) + "-" + atom_name + "_onebody" );
				ss->add_string_value( "pose_id", tag );
				ss->add_string_value( "resi1", std::to_string(seqpos) );
				ss->add_string_value( "pdbid1", pdb_id( *pose.pdb_info(), seqpos ) );
				ss->add_string_value( "restype1", res.name3() );
				ss->add_string_value( "atom1", std::to_string(atm) );
				ss->add_string_value( "name1", res.atom_name(atm) );
				ss->add_string_value( "resi2", "--" );
				ss->add_string_value( "pdbid2", "--" );
				ss->add_string_value( "restype2", "onebody" );
				ss->add_string_value( "atom2", "--" );
				ss->add_string_value( "name2", "--" );
				core::Real total = 0.0;
				for ( core::Size ii(1); ii <= active_scoretypes.size(); ++ii ) {
					total += onebody[ ii ];
					totals[ ii ] += onebody[ ii ];
					ss->add_energy(  name_from_score_type( active_scoretypes[ii] ), onebody[ ii ] );
				}
				ss->add_energy( "total", total );
				sfd.write_silent_struct( *ss, option[ out::file::silent ]() );
			}

			// Twobody
			for ( auto const & entry: core::scoring::get_pairwise_atom_energies(pose, *scorefxn, active_scoretypes, subset ) ) {
				core::Size seqpos1 = entry.first.first.rsd();
				core::Size atm1 = entry.first.first.atomno();
				core::Size seqpos2 = entry.first.second.rsd();
				core::Size atm2 = entry.first.second.atomno();
				utility::vector1< core::Real > const & twobody = entry.second;

				core::conformation::Residue const & res1( pose.residue(seqpos1) );
				std::string atom_name1 = utility::strip( res1.atom_name(atm1) );
				core::conformation::Residue const & res2( pose.residue(seqpos2) );
				std::string atom_name2 = utility::strip( res2.atom_name(atm2) );

				SilentStructOP ss( new ScoreFileSilentStruct(opts) );
				ss->decoy_tag( tag + "_" + std::to_string(seqpos1) + "-" + atom_name1 + "_" + std::to_string(seqpos2) + "-" + atom_name2 );
				ss->add_string_value( "pose_id", tag );
				ss->add_string_value( "resi1", std::to_string(seqpos1) );
				ss->add_string_value( "pdbid1", pdb_id( *pose.pdb_info(), seqpos1 ) );
				ss->add_string_value( "restype1", res1.name3() );
				ss->add_string_value( "atom1", std::to_string(atm1) );
				ss->add_string_value( "name1", res1.atom_name(atm1) );
				ss->add_string_value( "resi2", std::to_string(seqpos2) );
				ss->add_string_value( "pdbid2", pdb_id( *pose.pdb_info(), seqpos2 ) );
				ss->add_string_value( "restype2", res2.name3() );
				ss->add_string_value( "atom2", std::to_string(atm2) );
				ss->add_string_value( "name2", res2.atom_name(atm2) );
				core::Real total = 0.0;
				for ( core::Size ii(1); ii <= active_scoretypes.size(); ++ii ) {
					total += twobody[ ii ];
					totals[ ii ] += twobody[ ii ];
					ss->add_energy(  name_from_score_type( active_scoretypes[ii] ), twobody[ ii ] );
				}
				ss->add_energy( "total", total );
				sfd.write_silent_struct( *ss, option[ out::file::silent ]() );
			}

			// Add a line with grand totals (primarily for debugging purposes)
			SilentStructOP ss( new ScoreFileSilentStruct(opts) );
			ss->decoy_tag( tag + "_totals" );
			ss->add_string_value( "pose_id", tag );
			ss->add_string_value( "resi1", "--" );
			ss->add_string_value( "pdbid1", "--" );
			ss->add_string_value( "restype1", "--" );
			ss->add_string_value( "atom1", "--" );
			ss->add_string_value( "name1", "--" );
			ss->add_string_value( "resi2", "--" );
			ss->add_string_value( "pdbid2", "--" );
			ss->add_string_value( "restype2", "--" );
			ss->add_string_value( "atom2", "--" );
			ss->add_string_value( "name2", "--" );
			core::Real total = 0.0;
			for ( core::Size ii(1); ii <= active_scoretypes.size(); ++ii ) {
				total += totals[ ii ];
				ss->add_energy(  name_from_score_type( active_scoretypes[ii] ), totals[ ii ] );
			}
			ss->add_energy( "total", total );
			sfd.write_silent_struct( *ss, option[ out::file::silent ]() );
		} // while ( input.has_another_pose() )
		return 0;

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

}
