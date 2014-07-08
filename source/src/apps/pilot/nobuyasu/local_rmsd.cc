// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief
/// @author Nobuyasu Koga

// Unit headers

// Project headers
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/rms_util.hh>

// AUTO-REMOVED #include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

// type headers
#include <core/types.hh>

// Utility Headers
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

// C++ header
#include <fstream>
#include <map>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/excn/Exceptions.hh>


using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

typedef core::Size Size;
typedef std::string String;

static basic::Tracer TR("localrmsd");

namespace localrmsd
{
	//StringOptionKey blueprint( "foldptn:blueprint" );
	StringOptionKey output_name( "localrmsd:output_name" );
}

/// local mover for testing purposes
class LocalRmsd : public protocols::moves::Mover {
public:
	LocalRmsd():Mover("LocalRmsd"), window_( 7 ), skip_(3), initialize_( false )
	{

		// set native pose
		core::pose::Pose native_pose;
		if ( option[ in::file::native ].user() ) {
			core::import_pose::pose_from_pdb( native_pose, option[ in::file::native ]() );
			set_native_pose( new core::pose::Pose(native_pose) );
		}

		// read scorefxn
		//scorefxn_ = core::scoring::get_score_function();
		//		TR << "score will be calculated by remodel_cen.wts" << std::endl;
		scorefxn_ = core::scoring::get_score_function();

		// set output
		std::ostringstream filename;
		if( option[ localrmsd::output_name ].active() ){
			filename <<  option[ localrmsd::output_name ]() << ".localrmsd.dat";
		}else{
			filename <<  "localrmsd.dat";
		}
		output_.open( filename.str().c_str(), std::ios::out );

	}

	virtual ~LocalRmsd(){};

	virtual
	void
	apply( core::pose::Pose & pose ){

		if( ! initialize_ ){

			output_ << "# tag score rmsd  ";
			Size ntrial( pose.total_residue()/skip_ );
			for( Size i=1; i<=ntrial; i++ ){

				Size begin = ( i - 1 )*skip_ + 1;
				Size end = begin + window_ - 1;
				if( end > pose.total_residue() ){
					break;
				}
				output_ << begin << '-' << end;
			}
			output_ << std::endl;
			initialize_ = true;

		}

		using namespace protocols::jd2;
		JobOP job_me( JobDistributor::get_instance()->current_job() );
		String me( JobDistributor::get_instance()->job_outputter()->output_name( job_me ) );

		using namespace core::scoring;
		(*scorefxn_)(pose);
		Real totE = pose.energies().total_energies()[ scoring::total_score ];

		// set native pose
		core::pose::Pose const & native_pose( *get_native_pose() );

		//
		Real overall_rmsd ( CA_rmsd( pose, native_pose ) );

		//
		Size ntrial( pose.total_residue()/skip_ );
		utility::vector1< Real > local_rmsd;

		for( Size i=1; i<=ntrial; i++ ){

			Size begin = ( i - 1 )*skip_ + 1;
			Size end = begin + window_ - 1;
			if( end > pose.total_residue() ){
				break;
			}
			Real rmsd = CA_rmsd( pose, native_pose, begin, end );
			local_rmsd.push_back( rmsd );

		}

		output_ << me << ' ' << totE << ' ' << overall_rmsd << ' ';
		for( utility::vector1< Real>::iterator it=local_rmsd.begin(), ite=local_rmsd.end(); it != ite; ++it ){
			output_ << *it << ' ';
		}
		output_ << std::endl;

		/*
		job_me->add_string("somestuff");
		job_me->add_string_string_pair("alhpa", "btea");
		job_me->add_string_real_pair("theanswer", 42.0);
		JobDistributor::get_instance()->job_outputter()->other_pose( job_me, pose, "lookit");
		JobDistributor::get_instance()->job_outputter()->file( job_me, "hey here's some junk from job " + me + "\n" );
		*/


		return;
	}

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return new LocalRmsd;
	}


private:

	Size window_;
	Size skip_;
	bool initialize_;

	core::scoring::ScoreFunctionOP scorefxn_;
	std::ofstream output_;

};

typedef utility::pointer::owning_ptr< LocalRmsd > LocalRmsdOP;


int
main( int argc, char * argv [] )
{
	try{
	//
	option.add( localrmsd::output_name, "output name" );

	// init
  devel::init(argc, argv);

	// mover
	protocols::moves::MoverOP protocol;
	protocol = new LocalRmsd();

	// run
	protocols::jd2::JobDistributor::get_instance()->go( protocol );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
