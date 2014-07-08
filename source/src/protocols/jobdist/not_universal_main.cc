// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jobdist/not_universal_main.cc
/// @brief  Simple main method for applying a Mover to a set of
/// input Poses.
/// @author James Thompson

#include <core/types.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <basic/Tracer.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/viewer/viewers.hh>

#include <string>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>



namespace protocols {
namespace jobdist {

static basic::Tracer tr("protocols.jobdist.not_universal_main");

bool pose_matches_user_tag(
	core::pose::Pose & pose,
	utility::vector1< std::string > const & user_tags
) {
	using std::string;
	using utility::vector1;

	bool valid_tag( false );
	std::string this_id( "" );
	get_comment( pose, "user_tag", this_id );
	if ( this_id == "" ) get_comment( pose, "user_ta", this_id );
	if ( this_id == "" ) {
		tr.Error << "can't find user_tag in pose!" << std::endl;
	}

	typedef vector1< string >::const_iterator iter;
	for ( iter it = user_tags.begin(), end = user_tags.end();
				it != end && !valid_tag && this_id != ""; ++it
	) {
		if ( this_id.find( *it ) != std::string::npos ) {
			break;
			valid_tag = true;
			tr.Debug << "processing Pose with user_tag " << this_id << std::endl;
		}
	} // for user_tags

	if ( !valid_tag ) {
		tr.Warning << "skipping Pose with user_tag " << this_id << std::endl;
	}

	return valid_tag;
}

/////////////////////////////////////////////////////////////////////////////
/// @brief main function that applies a mover to a set of input poses and
/// only does silent-file output.  Similar in idea to universal_main, but
/// without all of the code duplication for silent-files/PDB files.
int not_universal_main(
	protocols::moves::Mover & mover
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;

	using core::Size;
	using core::Real;
	using std::string;
	using utility::vector1;

	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
		option[ in::file::residue_type_set ]()
	);

	core::pose::Pose native_pose, current_pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb(
			native_pose, *rsd_set, option[ in::file::native ]()
		);
	}

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	if ( option[ constraints::cst_weight ].user() ) {
		Real const weight( option[ constraints::cst_weight ]() );
		scorefxn->set_weight( core::scoring::angle_constraint,      weight );
		scorefxn->set_weight( core::scoring::dihedral_constraint,   weight );
		scorefxn->set_weight( core::scoring::atom_pair_constraint,  weight );
		scorefxn->set_weight( core::scoring::dunbrack_constraint,   weight );
		scorefxn->set_weight( core::scoring::coordinate_constraint, weight );
	}
	scorefxn->show( tr.Debug );
	tr.flush_all_channels();

	vector1< string > user_tags;
	if ( option[ in::file::user_tags ].user() ) {
		user_tags = option[ in::file::user_tags ]();
	}

	MetaPoseInputStream input = streams_from_cmd_line();

	Size const n_repeats( option[ run::repeat ]() );

	SilentFileData sfd_out;
	while( input.has_another_pose() ) {
		input.fill_pose( current_pose, *rsd_set );

		if ( option[ in::file::user_tags ].user() ) {
			if ( ! pose_matches_user_tag( current_pose, user_tags ) ) {
				continue;
			}
		}

		for ( Size ii = 1; ii <= n_repeats; ++ii ) {
			// only keep input scores if specified by the user.
			if ( !option[ in::file::keep_input_scores ]() ) {
				clearPoseExtraScores( current_pose );
			}

			protocols::viewer::add_conformation_viewer(
				current_pose.conformation(), "start_pose"
			);

			mover.apply( current_pose );

			// rescore the Pose if necessary, on by default.
			if ( option[ in::file::rescore ]() ) {
				(*scorefxn)(current_pose);
			}

			// output
			if ( !option[ out::nooutput ]() ) {
				core::io::silent::SilentStructOP ss
					= SilentStructFactory::get_instance()->get_silent_struct_out();
				ss->fill_struct( current_pose );

				// add model quality information if a native was given.
				if ( option[ in::file::native ].user() ) {
					core::Real CA_rmsd
						= core::scoring::CA_rmsd( native_pose, current_pose );
					core::Real CA_maxsub
						= core::scoring::CA_maxsub( native_pose, current_pose );
					core::Real CA_gdtmm
						= core::scoring::CA_gdtmm( native_pose, current_pose );
					ss->add_energy( "rmsd", CA_rmsd );
					ss->add_energy( "maxsub", CA_maxsub );
					ss->add_energy( "gdtmm", CA_gdtmm );
				}

				sfd_out.write_silent_struct( *ss, option[ out::file::silent ]() );
			} // if ( !option[ out::nooutput ]() )
		}
		tr.flush();
	} // while( input.has_another_pose() )
	tr.flush();

	return 0;
} // not_universal_main

} // jobdist
} // protocols
