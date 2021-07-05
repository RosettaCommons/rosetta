// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
//
/// @file   FlexPepDock.cc
//
/// @brief Main file for applying flexible peptide docking protocol
/// @author Barak Raveh
/// @date August 05, 2008

//#define GL_GRAPHICS

// Mini-Rosetta headers
#include <protocols/flexpep_docking/FlexPepDockingProtocol.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/flexPepDocking.OptionKeys.gen.hh>
#include <basic/citation_manager/CitationManager.hh>
#include <basic/citation_manager/CitationCollection.hh>
#include <core/pose/Pose.hh>

#include <protocols/jd2/JobDistributor.hh>

#include <protocols/moves/Mover.hh>



#include <devel/init.hh>

// C++ headers
//#include <cstdlib>

#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <utility/excn/Exceptions.hh>


#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif



//typedef utility::pointer::owning_ptr< BaseJobDistributor< BasicJobOP > > BaseJobDistributorOP
static basic::Tracer TR( "pilot_apps.FlexPepDock" );


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols;
		using namespace protocols::jobdist;
		using namespace protocols::moves;
		using basic::options::option;
		using namespace basic::options::OptionKeys;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// setup
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		devel::init(argc, argv);

		////////// Register with Citation Manager //////////
		{
			basic::citation_manager::CitationManager * cm ( basic::citation_manager::CitationManager::get_instance() );
			basic::citation_manager::CitationCollectionOP collection(
				utility::pointer::make_shared< basic::citation_manager::CitationCollection >(
				"FlexPepDocking", basic::citation_manager::CitedModuleType::Application
				)
			);
			collection->add_citation( cm->get_citation_by_doi("10.1002/prot.22716") );
			if ( option[basic::options::OptionKeys::flexPepDocking::lowres_abinitio]() ) {
				collection->add_citation( cm->get_citation_by_doi("10.1371/journal.pone.0018934") );
			}
			cm->add_citation( collection );
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// end of setup
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		MoverOP fpDock( utility::pointer::make_shared< flexpep_docking::FlexPepDockingProtocol >(1,true, true) );

		// read native pose: (TODO: look how this should be handled in Job Distributor 2)
		// protocols::jd2::set_native_in_mover(*fpDock);
		if ( option[ in::file::native ].user() ) {
			core::pose::PoseOP native_pose( utility::pointer::make_shared< core::pose::Pose >() );
			core::chemical::ResidueTypeSetCAP rsd_set;
			if ( option[ in::file::centroid_input ].user() ) {
				core::import_pose::centroid_pose_from_pdb( *native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
			} else {
				core::import_pose::pose_from_file( *native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
			}
			fpDock->set_native_pose( native_pose );
		}


		// run:
		protocols::jd2::JobDistributor::get_instance()->go(fpDock);
		// protocols::jobdist::main_plain_mover( *fpDock);
		//protocols::jobdist::universal_main(*fpDock);

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
