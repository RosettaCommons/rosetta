// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief quick demo program for cm_main
/// @author James Thompson

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/comparative_modeling/StealSideChainsMover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <protocols/jobdist/Jobs.hh>

#include <utility/excn/Exceptions.hh>


int
main( int argc, char * argv [] ) {
	try {

	// options, random initialization
	devel::init( argc, argv );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::sequence;
	ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set(
		option[ in::file::residue_type_set ]()
	);

	core::pose::Pose query_pose, template_pose;
	core::import_pose::pose_from_file(
		template_pose,
		*rsd_set,
		option[ in::file::template_pdb ]()[1]
	);

	SequenceAlignment aln;
	if ( option[ in::file::alignment ].user() ) {
		std::string align_fn  = option[ in::file::alignment ]()[1];
		aln.read_from_file( align_fn );
	}
	std::cout << aln << std::endl;

	protocols::moves::MoverOP mover(
		new protocols::comparative_modeling::StealSideChainsMover(
			template_pose, aln.sequence_mapping(1,2)
		)
	);
	protocols::jobdist::not_universal_main( *mover );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
