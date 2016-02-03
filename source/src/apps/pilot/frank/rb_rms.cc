// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <devel/init.hh>
#include <core/types.hh>

#include <protocols/jobdist/standard_mains.hh>
#include <protocols/loops/Loops.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/jobdist/Jobs.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <core/scoring/rms_util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/chemical/ChemicalManager.hh>

#include <utility/excn/Exceptions.hh>


#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/RBSegmentRelax.OptionKeys.gen.hh>

#include <list>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/rbsegment_relax.hh>


int main(int argc, char **argv) {
    try {
    using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;

	using core::Size;
	using std::string;

	devel::init( argc,argv );

	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
		option[ in::file::residue_type_set ]()
	);

	core::pose::Pose native_pose, current_pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_file( native_pose, *rsd_set, option[ in::file::native ]() , core::import_pose::PDB_file);
	}

	// rbsegfile
	utility::vector1< protocols::RBSegment::RBSegment > rbsegs;
	protocols::loops::Loops loops;
	std::string filename( option[ OptionKeys::RBSegmentRelax::rb_file ]().name() );
	protocols::RBSegment::read_RBSegment_file( rbsegs, loops, filename );
 	std::list< core::Size > core_reses;
	for ( core::Size i=1; i <= rbsegs.size(); ++i )
		for ( core::Size j=1, j_end=rbsegs[i].nContinuousSegments(); j<=j_end; ++j)
			for ( core::Size k=rbsegs[i][j].start(), k_end=rbsegs[i][j].end() ; k<=k_end; ++k)
				core_reses.push_back( k );


	MetaPoseInputStream input = streams_from_cmd_line();

	SilentFileData sfd_out;
	while( input.has_another_pose() ) {
		input.fill_pose( current_pose, *rsd_set );


		if ( option[ in::file::native ].user() ) {
			core::Real CA_rmsd = core::scoring::CA_rmsd( native_pose, current_pose );
			core::Real CA_core_rmsd = core::scoring::CA_rmsd( native_pose, current_pose , core_reses );

			std::cout << core::pose::extract_tag_from_pose( current_pose ) << "  " << CA_core_rmsd << "  " << CA_rmsd << std::endl;
		}

	} // while( input.has_another_pose() )
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
	return 0;
}
