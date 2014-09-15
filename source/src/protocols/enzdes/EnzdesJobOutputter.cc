// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzdes/EnzdesJobOutputter.cc
/// @brief implementation for EnzdesJobOutputter
/// @author Florian Richter (floric@u.washington.edu), september 2010



//unit headers
#include <protocols/enzdes/EnzdesJobOutputter.hh>
#include <protocols/enzdes/EnzdesJobOutputterCreator.hh>

//package headers
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>

//project headers
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>

#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/enzdes/EnzFilters.hh>

//utility headers
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace enzdes {



EnzdesJobOutputter::EnzdesJobOutputter()
: parent(),
  scorefile_writer_(NULL),
  enz_scofile_(NULL),
	silent_output_(false),
	silent_job_outputter_(NULL)
{

	if(  basic::options::option[ basic::options::OptionKeys::out::file::silent_struct_type ].user() ){
		silent_output_ = true;
		silent_job_outputter_ = new protocols::jd2::SilentFileJobOutputter;
		silent_job_outputter_->set_write_separate_scorefile( false );
		silent_job_outputter_->set_silent_file_name( utility::file::FileName( basic::options::option[ basic::options::OptionKeys::out::file::o ]) );
	}

  if( !this->write_scorefile() ) return;
  scorefile_writer_ = new core::io::silent::SilentFileData();
  enz_scofile_ = new protocols::enzdes::EnzdesScorefileFilter();
  if( basic::options::option[ basic::options::OptionKeys::out::overwrite ].user() ){
    if( utility::file::file_exists( scorefile_name().name() ) ) utility::file::file_delete( scorefile_name().name() );
  }
}

EnzdesJobOutputter::~EnzdesJobOutputter(){}

void
EnzdesJobOutputter::final_pose( protocols::jd2::JobOP job, core::pose::Pose const & pose, std::string const & tag )
{
	if( silent_output_ ){
		silent_job_outputter_->final_pose(job, pose, tag);
		this->scorefile( job, pose );
	}
	else parent::final_pose(job, pose, tag);
}

bool
EnzdesJobOutputter::job_has_completed( protocols::jd2::JobCOP job )
{
	if( silent_output_ ) return silent_job_outputter_->job_has_completed( job );
	return parent::job_has_completed( job );
}

void
EnzdesJobOutputter::scorefile(
  protocols::jd2::JobCOP job,
  core::pose::Pose const & pose,
    std::string, //tag,
    std::string /*suffix_tag*/,
    std::string //scorefile
)
{
  if( !this->write_scorefile() ) return;

  protocols::toolbox::match_enzdes_util::EnzdesCacheableObserverCOP enz_obs( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose ) );
  if( enz_obs ){
    if( enz_obs->cst_cache() ) enz_scofile_->set_cstio( enz_obs->cst_cache()->enzcst_io() );
  }

  enz_scofile_->examine_pose( pose );
  utility::vector1< core::io::silent::SilentEnergy > silent_Es( enz_scofile_->silent_Es() );

  //also add eventual filter reports to the scorefile
  for( protocols::jd2::Job::StringRealPairs::const_iterator it(job->output_string_real_pairs_begin()), it_end(job->output_string_real_pairs_end());
       it != it_end; ++it){
    silent_Es.push_back( core::io::silent::SilentEnergy(  it->first, it->second, 1.0, std::max( 10, (int) (it->first).length() + 3 )  ) );
  }

  std::string outstruct_name( this->output_name( job ) );
  core::io::silent::SilentStructOP ss( new core::io::silent::ScoreFileSilentStruct( pose, outstruct_name ) );
  ss->precision( 2 );
  ss->scoreline_prefix( "" );
  ss->silent_energies( silent_Es );

  scorefile_writer_->write_silent_struct( *ss, scorefile_name().name(), true );

  if( basic::options::option[ basic::options::OptionKeys::enzdes::final_repack_without_ligand] &&
      basic::options::option[ basic::options::OptionKeys::enzdes::dump_final_repack_without_ligand_pdb ] ){
    utility::io::ozstream out( output_name(job)+"_nlrepack.pdb" );
    if ( !out.good() ) utility_exit_with_message( "Unable to open file: " + output_name(job)+"_nlrepack.pdb" + "\n" );
    this->dump_pose( job, *(enz_scofile_->rnl_pose()) , out);
    out.close();
    enz_scofile_->clear_rnl_pose();
  }
}

//CREATOR SECTION
std::string
EnzdesJobOutputterCreator::keyname() const
{
        return "EnzdesJobOutputter";
}

protocols::jd2::JobOutputterOP
EnzdesJobOutputterCreator::create_JobOutputter() const {
        return new EnzdesJobOutputter;
}

}//enzdes
}//protocols




