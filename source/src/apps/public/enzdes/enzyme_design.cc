// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
///
/// @brief
/// @author Florian Richter (floric@u.washington.edu)

#include <protocols/enzdes/EnzdesFixBBProtocol.hh>
#include <protocols/enzdes/EnzdesFlexBBProtocol.hh>
//#include <devel/enzdes/EnzdesRemodelProtocol.hh>
#include <protocols/enzdes/enzdes_util.hh>
#include <protocols/enzdes/EnzFilters.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/CacheableString.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh> //function for reading cstfiles
//#include <devel/enzdes/DesignSilentStruct.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/JobDistributors.hh>

#include <utility/string_util.hh>


#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

#include <core/io/silent/SilentFileData.hh>
#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <utility/excn/Exceptions.hh>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif


static thread_local basic::Tracer tr( "apps.public.enzdes.enzyme_design" );

using namespace core;

/////////////////////////////////////////

////main function
int
main( int argc, char * argv [])
{
	try {
	protocols::enzdes::EnzdesBaseProtocol::register_options();
	protocols::enzdes::EnzdesFixBBProtocol::register_options();
	protocols::enzdes::EnzdesFlexBBProtocol::register_options();
	basic::options::option.add_relevant( basic::options::OptionKeys::enzdes::flexbb_protocol );
	basic::options::option.add_relevant( basic::options::OptionKeys::enzdes::process_ligrot_separately );
	basic::options::option.add_relevant( basic::options::OptionKeys::out::file::o );

	devel::init(argc, argv);

	time_t overall_start_time = time(NULL);
	utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();

	//#ifndef USEMPI
	// Reduce read contention between processes by randomizing the order in which structures are processed
	// Do not randomize, though, if job distribution is controlled by MPI
	//numeric::random::random_permutation( input_jobs, numeric::random::rg() );
	//#endif

	protocols::jobdist::PlainPdbJobDistributor jobdist( input_jobs );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if( option[ out::nooutput ]() ){
		jobdist.disable_output();
		jobdist.enable_ignorefinished();
	}


	protocols::enzdes::EnzdesBaseProtocolOP enzdes_protocol;
	protocols::enzdes::EnzdesScorefileFilterOP enz_scofile;

	if( option[ OptionKeys::enzdes::flexbb_protocol ] ) enzdes_protocol = protocols::enzdes::EnzdesBaseProtocolOP( new protocols::enzdes::EnzdesFlexBBProtocol() );
	//else if( option[ OptionKeys::enzdes::remodel_protocol ] ) enzdes_protocol = new devel::enzdes::EnzdesRemodelProtocol();
	else enzdes_protocol = protocols::enzdes::EnzdesBaseProtocolOP( new protocols::enzdes::EnzdesFixBBProtocol() );

	std::string scorefile_name("");

	if( option[ OptionKeys::out::file::o ].user() ){
		scorefile_name = option[ OptionKeys::out::file::o ]();
		enz_scofile = protocols::enzdes::EnzdesScorefileFilterOP( new protocols::enzdes::EnzdesScorefileFilter() );
		//enz_scofile->set_cstio( enzdes_protocol->cst_io() );
		if( option[ OptionKeys::out::overwrite ].user() ){
			if( utility::file::file_exists( scorefile_name ) ) utility::file::file_delete( scorefile_name );
		}
	}
	core::io::silent::SilentFileDataOP scorefile( new core::io::silent::SilentFileData() );

	if( option[OptionKeys::enzdes::cstfile].user() ){
		option[OptionKeys::run::preserve_header ].value(true);
		//enzdes_protocol->cst_io()->read_enzyme_cstfile(basic::options::option[basic::options::OptionKeys::enzdes::cstfile]);
	}


	protocols::jobdist::BasicJobOP curr_job, prev_job;
	int curr_nstruct=0, num_structures_processed = 0;
	core::Size prevstruct(0);
	core::pose::PoseOP input_pose; // starts NULL, coords *never* modified!
	jobdist.startup();

	while( jobdist.next_job(curr_job, curr_nstruct) ) {
		time_t pdb_start_time = time(NULL);
		tr << "Starting " << curr_job->output_tag(curr_nstruct) << " ..." << std::endl;

		if( curr_nstruct == 1 ) prevstruct = 0;

		// we read each PDB just once to save on disk I/O
		if( curr_job.get() != prev_job.get() || input_pose.get() == NULL ) {
			input_pose = core::pose::PoseOP( new core::pose::Pose() );
			//core::import_pose::pose_from_pdb( *input_pose, curr_job->input_tag() );
			protocols::enzdes::enzutil::read_pose_from_pdb(  *input_pose, curr_job->input_tag() );
		}

		// Make a modifiable copy of the pose read from disk
		core::pose::PoseOP the_pose( new core::pose::Pose( *input_pose ) );
		std::string outext = "_DE";
		if( basic::options::option[basic::options::OptionKeys::out::suffix].user() ){
			outext = basic::options::option[basic::options::OptionKeys::out::suffix];
		}
		std::string prefix = "";
    if( basic::options::option[basic::options::OptionKeys::out::prefix].user() ){
      prefix = basic::options::option[basic::options::OptionKeys::out::prefix];
    }

		utility::vector1< core::pose::PoseOP > poses_to_process;

		if( option[ OptionKeys::enzdes::process_ligrot_separately ].user() ){

			(*enzdes_protocol->get_scorefxn() )( *the_pose );
			enzdes_protocol->generate_explicit_ligand_rotamer_poses( *the_pose, poses_to_process, enzdes_protocol->get_scorefxn() );

			if( poses_to_process.size() == 0 ) poses_to_process.push_back( the_pose );
			else tr << "For " << curr_job->input_tag() << ", " << poses_to_process.size() << " explicit ligrot poses will be processed." << std::endl;
		}
		else poses_to_process.push_back( the_pose );

		for( core::Size pose_count = 1; pose_count <= poses_to_process.size(); ++pose_count){

			utility::file::FileName out_name( curr_job->input_tag() );
			//std::string outtag = out_name.base() + "_" + utility::to_string( curr_nstruct ) + "_" + outext + utility::to_string( pose_count + prevstruct );
			std::string outtag = prefix + out_name.base() + "_" + outext + "_" + utility::to_string( pose_count + prevstruct );

			if( utility::file::file_exists( outtag+".pdb" ) && ! option[ OptionKeys::out::overwrite ].user() ){
				tr << "File " << outtag+".pdb" << " already exists, skipping structure. Use option -out::overwrite if you want to overwrite existing files." << std::endl;
				continue;
			}

			using namespace basic::datacache;
			(poses_to_process[ pose_count ])->data().set(core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG, DataCache_CacheableData::DataOP( new basic::datacache::CacheableString( outtag ) ) );
			enzdes_protocol->apply( *(poses_to_process[ pose_count ]) );

				//in case we're only interested in scoring and there is an output
				//file, dont write a pdb (can be overriden by nstruct)
				if( !( ( option[OptionKeys::enzdes::enz_score] ) &&
						( option[ OptionKeys::out::file::o ].user() ) &&
						( ! option[ OptionKeys::out::nstruct ].user() ) ) ){
					jobdist.dump_pose_and_map( outtag,  *(poses_to_process[ pose_count ]) );
				}

				//write the score file, eventually maybe even a silent file
				if( basic::options::option[ basic::options::OptionKeys::out::file::o ].user() ){

					protocols::toolbox::match_enzdes_util::EnzConstraintIOCOP cstio( protocols::enzdes::enzutil::get_enzcst_io( *(poses_to_process[ pose_count ] ) ) );
					enz_scofile->set_cstio( cstio );
					core::io::silent::SilentStructOP ss( new core::io::silent::ScoreFileSilentStruct(
								*(poses_to_process[ pose_count ]), outtag ) );
					ss->precision( 2 );
					ss->scoreline_prefix( "" );

					enz_scofile->examine_pose( *(poses_to_process[ pose_count ]) );

					if( option[OptionKeys::enzdes::final_repack_without_ligand] ){
						//pose after ligand repack get's written as model 2
						if( option[OptionKeys::enzdes::dump_final_repack_without_ligand_pdb] ){
							jobdist.dump_pose_and_map( outtag+"_nlrepack", *(enz_scofile->rnl_pose() ) );
						}
					} //if( option[OptionKeys::enzdes::final_repack_without_ligand] )

					ss->silent_energies( enz_scofile->silent_Es() );
					scorefile->write_silent_struct( *ss, scorefile_name, true );
				}
		} // iterator over poses to process

		prevstruct += poses_to_process.size();

		prev_job = curr_job; // pointer assignment, not a copy op
		num_structures_processed += 1;
		time_t pdb_end_time = time(NULL);
		//std::cerr << curr_job->output_tag(curr_nstruct) << " done." << std::endl;
		tr << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;
	} // loop over jobs and nstructs
	jobdist.shutdown();

	time_t overall_end_time = time(NULL);
	tr << "Finished all " << num_structures_processed << " structures in " << (long)(overall_end_time - overall_start_time) << " seconds." << std::endl;
	if ( num_structures_processed == 0 )
		basic::Warning() << "No structures processed.  Existing output files may have been skipped, did you mean to delete them or to use the -overwrite flag?" << std::endl;
	return 0;
	}
	catch (utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

