// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   sandbox
/// @brief  apps/pilot/yab/ufv.cc
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// project headers
#include <devel/init.hh>
#include <devel/init.hh>
#include <core/chemical/ChemicalManager.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/ufv.OptionKeys.gen.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/Tracer.hh>
#include <protocols/forge/build/ConnectRight.hh>
#include <protocols/forge/build/RelativeConnectRight.hh>
#include <protocols/forge/build/RelativeSequencePosition.hh>
#include <protocols/forge/build/SegmentInsert.hh>
#include <protocols/forge/build/SegmentRebuild.hh>
#include <protocols/forge/components/BDR.hh>
#include <protocols/forge/remodel/RemodelMover.hh>
#include <protocols/forge/remodel/RemodelGlobalFrame.hh>
#include <protocols/forge/remodel/RemodelGlobalFrame.fwd.hh>
#include <protocols/simple_filters/PoseMetricEvaluator.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/viewer/viewers.hh>

#include <protocols/moves/PyMolMover.hh>

// utility headers
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/options/keys/OptionKey.hh>

// boost headers
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

// C++ headers
#include <string>
#include <vector>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>



// using
using utility::options::OptionKey;


// typedefs
typedef std::string String;
typedef std::vector< OptionKey const * > KeyVec;


// static
static basic::Tracer TR( "apps.pilot.possu.helical_param" );

void fill_required_options( KeyVec & keys ) {
	using namespace basic::options::OptionKeys;

	keys.push_back( &in::path::database );

	keys.push_back( &out::nstruct );

	keys.push_back( &run::max_retry_job );

	keys.push_back( &pose_metrics::neighbor_by_distance_cutoff );
}


void fill_optional_options( KeyVec & keys ) {
	using namespace basic::options::OptionKeys;


	keys.push_back( &in::file::s );
	keys.push_back( &in::file::silent );
	keys.push_back( &in::file::vall );
	keys.push_back( &remodel::blueprint);

	keys.push_back( &out::pdb_gz );
	keys.push_back( &out::pdb );
	keys.push_back( &out::silent_gz );

	keys.push_back( &out::file::o );
	keys.push_back( &out::file::silent );
	keys.push_back( &out::file::silent_struct_type );

	// either the following three or 'ufv_loops' is required; conditional
	// options don't appear easily encodeable, so just mark these as optional
	// and check the conditions explicitly
	/*
	keys.push_back( &Remodel::default_centroid_aa_as );
	keys.push_back( &Remodel::keep_junction_torsions );
	keys.push_back( &Remodel::use_fullmer );
	keys.push_back( &Remodel::centroid_loop_mover );
	keys.push_back( &Remodel::no_neighborhood_design );
	keys.push_back( &Remodel::dr_cycles );
	keys.push_back( &Remodel::centroid_sfx );
	keys.push_back( &Remodel::centroid_sfx_patch );
	keys.push_back( &Remodel::fullatom_sfx );
	keys.push_back( &Remodel::fullatom_sfx_patch );

	keys.push_back( &Remodel::insert::insert_pdb );
	keys.push_back( &Remodel::insert::attached_pdb );
	keys.push_back( &Remodel::insert::connection_scheme );
	*/

}


void register_options( KeyVec & keys ) {
	using basic::options::option;

	for ( KeyVec::const_iterator i = keys.begin(), ie = keys.end(); i != ie; ++i ) {
		option.add_relevant( **i );
	}
}


bool check_required_options( KeyVec & keys ) {
	using basic::options::option;

	bool flags_ok = true;
	for ( KeyVec::const_iterator i = keys.begin(), ie = keys.end(); i != ie; ++i ) {
		flags_ok &= option[ **i ].specified_report();
	}

	return flags_ok;
}


bool check_option_conflicts() {
	using namespace basic::options::OptionKeys;
	using core::pose::annotated_to_oneletter_sequence;
	using basic::options::option;

	bool flags_ok = true;

	//if ( !option[ ufv::ufv_loops ].user() ) {
	//
	//  if ( !option[ ufv::ss ].user() ) {
	//   flags_ok &= false;
	//   TR.Fatal << "no ufv:ufv_loops option specified, so ufv:ss required!" << std::endl;
	//  }
	//
	//  if ( !option[ ufv::left ].user() ) {
	//   flags_ok &= false;
	//   TR.Fatal << "no ufv:ufv_loops option specified, so ufv:left required!" << std::endl;
	//  }
	//
	//  if ( !option[ ufv::right ].user() ) {
	//   flags_ok &= false;
	//   TR.Fatal << "no ufv:ufv_loops option specified, so ufv:right required!" << std::endl;
	//  }
	//
	//  if ( option[ ufv::aa_during_build ].user() ) {
	//   String aa = core::pose::annotated_to_oneletter_sequence( option[ ufv::aa_during_build ] );
	//   flags_ok &= aa.length() == option[ ufv::ss ].value().length();
	//
	//   if ( !flags_ok ) {
	//    TR.Fatal << "length of ufv:ss string and one letter version of ufv::aa_during_build must be equal!" << std::endl;
	//   }
	//  }
	//
	//  if ( option[ ufv::aa_during_design_refine ].user() ) {
	//   String aa = core::pose::annotated_to_oneletter_sequence( option[ ufv::aa_during_design_refine ] );
	//   flags_ok &= aa.length() == option[ ufv::ss ].value().length();
	//
	//   if ( !flags_ok ) {
	//    TR.Fatal << "length of ufv:ss string and one letter version of ufv::aa_during_design_refine must be equal!" << std::endl;
	//   }
	//  }
	// }

	if ( option[ run::max_retry_job ] < 0 ) {
		flags_ok = false;
		TR.Fatal << "run:max_retry_job must be positive!" << std::endl;
	}

	return flags_ok;
}


/// @brief load ufv loops from file
/// @return number of loops read
core::Size load_loops_from_file(
	protocols::forge::components::BDR & bdr,
	utility::file::FileName const & filename
)
{
	using namespace basic::options::OptionKeys;
	using core::Size;
	using basic::options::option;
	using protocols::forge::build::Interval;
	using protocols::forge::build::SegmentRebuild;
	using std::istringstream;
	typedef std::string String;

	utility::io::izstream in( filename );

	Size left, right;
	String ss, aa_during_build, aa_during_design_refine;
	bool use_sequence_biased_fragments = false;

	Size count = 0;
	String line;
	utility::vector1< String > entries;
	while ( getline( in, line ) ) {
		// split by whitespace " \n\t"
		boost::tokenizer< boost::char_separator< char > > tokens( line, boost::char_separator< char >( " \n\t" ) );
		entries.assign( tokens.begin(), tokens.end() );

		// skip empty and comment lines
		if ( entries.size() == 0 || ( entries.size() > 0 && entries.begin()->at( 0 ) == '#' ) ) {
			continue;
		}

		if ( entries.size() >= 3 ) {
			left = boost::lexical_cast< Size >( entries[ 1 ] );
			right = boost::lexical_cast< Size >( entries[ 2 ] );
			ss = entries[ 3 ];

			aa_during_build = entries.size() >= 4 ? entries[ 4 ] : String();
			aa_during_design_refine = entries.size() >= 5 ? entries[ 5 ] : String();

			if ( aa_during_build.length() == 1 && aa_during_build.at( 0 ) == '-' ) {
				aa_during_build = String();
			}

			use_sequence_biased_fragments |= !aa_during_build.empty();

			bdr.add_instruction(
				new SegmentRebuild(
				Interval( left, right ),
				ss, aa_during_build,
				core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ),
				option[ ufv::keep_junction_torsions ]
				),
				aa_during_design_refine
			);

			++count;

			TR << "added loop: [" << left << ", " << right << "]  " << ss << "  "
				<< aa_during_build << "  " << aa_during_design_refine << std::endl;

		} else {
			// error message
			TR.Error << "ERROR: too few columns in line of ufv_loops file: " << line << std::endl;
			std::exit( 1 );
		}
	}

	in.close();

	bdr.use_sequence_bias( use_sequence_biased_fragments );

	return count;
}




void * graphics_main( void * ) {
	using namespace basic::options::OptionKeys;
	using namespace protocols::forge::remodel;
	using core::Size;
	using basic::options::option;
	using protocols::simple_filters::PoseMetricEvaluator;
	using protocols::forge::build::Interval;
	using protocols::forge::build::SegmentRebuild;
	using protocols::forge::remodel::RemodelMover;
	using protocols::forge::remodel::RemodelMover_OP;
	using protocols::forge::components::BDR_OP;
	using protocols::jd2::JobDistributor;

	// option transmute; just for safety
	String aa_during_build;
	String aa_during_design_refine;

	if ( option[ ufv::aa_during_build ].user() ) {
		aa_during_build = option[ ufv::aa_during_build ];
	}

	if ( option[ ufv::aa_during_design_refine ].user() ) {
		aa_during_design_refine = option[ ufv::aa_during_design_refine ];
	}

	RemodelGlobalFrame_OP RGF = new RemodelGlobalFrame(option[remodel::repeat_structure]);

	core::pose::Pose input_pose;
	core::import_pose::pose_from_pdb(input_pose, option[in::file::s]().vector()[0]);

	RGF->align_segment(input_pose);
	std::cout << "dumping..." << std::endl;
	input_pose.dump_pdb("aligned.pdb");

	// run job
	JobDistributor::get_instance()->go( RGF );




	return 0;
}


int main( int argc, char * argv [] ) {
	using namespace basic::options::OptionKeys;
	using basic::options::option;

	// track options
	//KeyVec required_options;
	//fill_required_options( required_options );

	KeyVec optional_options;
	fill_optional_options( optional_options );

	// register options so help file appears correctly
	//register_options( required_options );
	register_options( optional_options );

	// initialize rosetta
	devel::init( argc, argv );
	devel::init(argc, argv);


	protocols::viewer::viewer_main( graphics_main );

	return 0;
}
