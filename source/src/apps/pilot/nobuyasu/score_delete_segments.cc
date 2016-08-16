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
/// @author Nobuyasu Koga

// Package headers
// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>
#include <utility/string_util.hh>
#include <boost/lexical_cast.hpp>

// C++ header
#include <fstream>
#include <map>
#include <utility/excn/Exceptions.hh>


using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

typedef std::string String;


static THREAD_LOCAL basic::Tracer TR( "pick_bab" );

class ThisApplication  {
public:
	ThisApplication(){};
	static void register_options();
};

// OPT_KEY( File, blue )
// OPT_KEY( Integer, abego_level )
OPT_KEY( String, segments )
OPT_KEY( Boolean, ignore_error )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( segments, "deleted segments, specified as like 20-43,55-79", "" );
	NEW_OPT( ignore_error, "ignore error reading silent file, and continue the calculation.", false );
}


/// local mover for testing purposes
class DeleteSegments : public protocols::moves::Mover {
public: // typedef


public: // constructor/deconstructor


	DeleteSegments()
	{
		using utility::string_split;

		// set delted segments
		utility::vector1< String > parts;
		if( option[ segments ]() != "" ) {
			parts = string_split( option[ segments ](), ',' );
		} else {
			utility_exit_with_message( "[ERROR] No segments to be deleted. " );
		}
		if( parts.size() < 1 ) {
			utility_exit_with_message( "[ERROR] Segments are not defined properly: ',' is required to separete segments. " );
		}

		utility::vector1< String > residue_pair;
		for( Size ii=1; ii<=parts.size(); ii++ ) {
			residue_pair = string_split( parts[ ii ], '-' );
			if( residue_pair.size() != 2 ) {
				utility_exit_with_message( "[ERROR] Segment is not defined properly, '-' is required between begin and end positions. " );
			}
			TR << "deleted segment : " << residue_pair[ 1 ] << " " << residue_pair[ 2 ]  << std::endl;
			Size begin = boost::lexical_cast<Size>( residue_pair[ 1 ] );
			Size end = boost::lexical_cast<Size>( residue_pair[ 2 ] );
			runtime_assert( begin < end );
			segments_region_.insert( std::pair< Size, Size > ( begin, end ) );
		}

		// set scorefxn
		scorefxn_ = core::scoring::get_score_function();

		if ( option[ in::file::native ].user() ) {
			core::import_pose::pose_from_file( native_, option[ basic::options::OptionKeys::in::file::native ] , core::import_pose::PDB_file);
		}
	}
	virtual ~DeleteSegments(){}

	virtual std::string get_name() const { return ""; }

	virtual
	protocols::moves::MoverOP
	fresh_instance() const
	{
		return new DeleteSegments;
	}


public: // apply


	virtual
	void
	apply( core::pose::Pose & pose )
	{
		using namespace protocols::jd2;
		using core::util::switch_to_residue_type_set;
		JobOP job_me( JobDistributor::get_instance()->current_job() );
		String me( JobDistributor::get_instance()->job_outputter()->output_name( job_me ) );

		if( ! pose.is_fullatom() ) {
			core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
		}
		Pose copy_pose( pose );

		Size removed_residues( 0 );
		std::map< Size, Size >::const_iterator ite = segments_region_.end();
		for ( std::map< Size, Size >::const_iterator it = segments_region_.begin(); it != ite; ++it ) {
			Size begin = it->first - removed_residues;
			Size end = it->second - removed_residues;

			if( begin < 1 || end > pose.total_residue() ) {
        if( option[ ignore_error ]() ) {
					TR << "Segments to be deleted are out of range of pose. There might be a junk in a readin silent files. " << me << std::endl;
					TR << pose.total_residue() << " " << begin << " " << end << std::endl;
					continue;
				} else {
					utility_exit_with_message(" Segment to be deleted is out of range of pose ");
				}
			}

			pose.conformation().delete_residue_range_slow( begin, end );
			removed_residues += end - begin + 1;
		}

		( *scorefxn_ )( pose );


		if( option[ in::file::native ].user() ) {
			if( copy_pose.total_residue() != native_.total_residue() ) {
				TR << copy_pose.total_residue() << " " << native_.total_residue() << " " << me << std::endl;
			}
			Real rms = core::scoring::CA_rmsd( copy_pose, native_ );
			job_me->add_string_real_pair( "rms", rms );
		}
	}

private: // data

	core::scoring::ScoreFunctionOP scorefxn_;
	std::map< Size, Size > segments_region_;
	Pose native_;

};

typedef utility::pointer::owning_ptr< DeleteSegments > DeleteSegmentsOP;


int
main( int argc, char * argv [] )
{
	try{
	ThisApplication::register_options();

	// init
	devel::init(argc, argv);

	// mover
	protocols::moves::MoverOP protocol;
	protocol = new DeleteSegments();

	// run
	protocols::jd2::JobDistributor::get_instance()->go( protocol );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
