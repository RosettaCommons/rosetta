// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author James Thompson

#include <devel/init.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/options/option.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/chemical/AA.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <core/scoring/Energies.hh>
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
#include <numeric/random/random.hh>

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <core/types.hh>
#include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/sasa.hh>

#include <utility/vector1.hh>

#include <string>
#include <iostream>

//Auto Headers
#include <core/pose/util.hh>


core::Real
calculate_backbone_sasa(core::id::AtomID_Map< core::Real > & atom_sasa, core::pose::Pose & current_pose, core::Size & resn) {

	core::Real backbone_sasa = 0.0;

	core::id::AtomID const id_n( current_pose.residue_type(resn).atom_index("N"),resn);
	core::id::AtomID const id_ca( current_pose.residue_type(resn).atom_index("CA"),resn);
	core::id::AtomID const id_c( current_pose.residue_type(resn).atom_index("C"),resn);
	core::id::AtomID const id_o( current_pose.residue_type(resn).atom_index("O"),resn);

	backbone_sasa += atom_sasa[id_n];
	backbone_sasa += atom_sasa[id_ca];
	backbone_sasa += atom_sasa[id_c];
	backbone_sasa += atom_sasa[id_o];

	if (core::chemical::name_from_aa(current_pose.aa(resn)) == "GLY") {
		core::id::AtomID const id_1ha( current_pose.residue_type(resn).atom_index("1HA"),resn);
		core::id::AtomID const id_2ha( current_pose.residue_type(resn).atom_index("2HA"),resn);
		backbone_sasa += atom_sasa[id_1ha];
		backbone_sasa += atom_sasa[id_2ha];
	} else {
		core::id::AtomID const id_ha( current_pose.residue_type(resn).atom_index("HA"),resn);
		backbone_sasa += atom_sasa[id_ha];
	}

	return backbone_sasa;
}

int
main( int argc, char* argv [] ) {

	try {

	// options, random initialization
	devel::init( argc, argv );

	using core::Size;
	using core::Real;
	using core::pose::Pose;
	using namespace core::chemical;
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

	utility::vector1< Size > exposed_residues;
	if ( option[ score::sidechain_exposed ].user() ) {
		exposed_residues = option[ score::sidechain_exposed ]();
	}

	utility::vector1< Size > buried_residues;
	if ( option[ score::sidechain_buried ].user() ) {
		buried_residues = option[ score::sidechain_buried ]();
	}

	runtime_assert( ( exposed_residues.size() + buried_residues.size() ) != 0 );

	SilentFileData sfd;
	core::pose::Pose current_pose;
	while ( input.has_another_pose() ) {
		input.fill_pose( current_pose, *rsd_set );
		(*scorefxn)(current_pose);
		EnergyMap weights( current_pose.energies().weights() );

		core::scoring::dssp::Dssp dssp( current_pose );
		dssp.insert_ss_into_pose( current_pose );
		std::string secstruct(current_pose.secstruct());

		core::id::AtomID_Map< bool > atom_subset;
		core::pose::initialize_atomid_map( atom_subset, current_pose, true );

		core::id::AtomID_Map< Real > atom_sasa;
		utility::vector1< Real > rsd_sasa;

		core::Real probe_radius = 1.5;
		bool use_big_polar_H = false;

		core::scoring::calc_per_atom_sasa( current_pose, atom_sasa, rsd_sasa, probe_radius, use_big_polar_H, atom_subset );


		core::Real sasa_score(0.0);
		SilentStructOP ss( new ScoreFileSilentStruct );

		for ( Size ti = 1; ti <= exposed_residues.size(); ++ti ) {

			Size tr = exposed_residues[ti];

			runtime_assert( tr <= current_pose.size() );

			core::Real sc_sasa = rsd_sasa[tr] - calculate_backbone_sasa(atom_sasa,current_pose,tr);

			if ( sc_sasa == 0.0 ) {
				sasa_score += 1;
			}

		}

		for ( Size ti = 1; ti <= buried_residues.size(); ++ti ) {

			Size tr = buried_residues[ti];

			runtime_assert( tr <= current_pose.size() );

			core::Real sc_sasa = rsd_sasa[tr] - calculate_backbone_sasa(atom_sasa,current_pose,tr);

			if ( sc_sasa >= 0.0 ) {
				sasa_score += 1;
			}

		}

		ss->add_energy("sasa_score", sasa_score );
		ss->set_decoy_tag( core::pose::tag_from_pose(current_pose) );

		sfd.write_silent_struct( *ss, option[ out::file::scorefile ]() );

	} // while ( input.has_another_pose() )

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
