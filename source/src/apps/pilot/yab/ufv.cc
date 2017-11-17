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
#include <core/chemical/ChemicalManager.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/ufv.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/Tracer.hh>
#include <protocols/forge/build/ConnectRight.hh>
#include <protocols/forge/build/RelativeConnectRight.hh>
#include <protocols/forge/build/RelativeSequencePosition.hh>
#include <protocols/forge/build/SegmentInsert.hh>
#include <protocols/forge/build/SegmentRebuild.hh>
#include <protocols/forge/components/BDR.hh>
#include <protocols/evaluation/PoseMetricEvaluator.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/viewer/viewers.hh>

// utility headers
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/excn/Exceptions.hh>

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
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/jumping/StrandPairing.hh>


// using
using utility::options::OptionKey;


// typedefs
typedef std::string String;
typedef std::vector< OptionKey const * > KeyVec;


// static
static basic::Tracer TR( "apps.pilot.yab.ufv" );


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

	keys.push_back( &out::pdb_gz );
	keys.push_back( &out::pdb );
	keys.push_back( &out::silent_gz );

	keys.push_back( &out::file::o );
	keys.push_back( &out::file::silent );
	keys.push_back( &out::file::silent_struct_type );

	// either the following three or 'ufv_loops' is required; conditional
	// options don't appear easily encodeable, so just mark these as optional
	// and check the conditions explicitly
	keys.push_back( &ufv::ufv_loops );
	keys.push_back( &ufv::left );
	keys.push_back( &ufv::right );
	keys.push_back( &ufv::ss );

	keys.push_back( &ufv::aa_during_build );
	keys.push_back( &ufv::aa_during_design_refine );
	keys.push_back( &ufv::keep_junction_torsions );
	keys.push_back( &ufv::use_fullmer );
	keys.push_back( &ufv::centroid_loop_mover );
	keys.push_back( &ufv::no_neighborhood_design );
	keys.push_back( &ufv::dr_cycles );
	keys.push_back( &ufv::centroid_sfx );
	keys.push_back( &ufv::centroid_sfx_patch );
	keys.push_back( &ufv::fullatom_sfx );
	keys.push_back( &ufv::fullatom_sfx_patch );

	keys.push_back( &ufv::insert::insert_pdb );
	keys.push_back( &ufv::insert::attached_pdb );
	keys.push_back( &ufv::insert::connection_scheme );

	keys.push_back( &packing::resfile );
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

	if ( !option[ ufv::ufv_loops ].user() ) {

		if ( !option[ ufv::ss ].user() ) {
			flags_ok &= false;
			TR.Fatal << "no ufv:ufv_loops option specified, so ufv:ss required!" << std::endl;
		}

		if ( !option[ ufv::left ].user() ) {
			flags_ok &= false;
			TR.Fatal << "no ufv:ufv_loops option specified, so ufv:left required!" << std::endl;
		}

		if ( !option[ ufv::right ].user() ) {
			flags_ok &= false;
			TR.Fatal << "no ufv:ufv_loops option specified, so ufv:right required!" << std::endl;
		}

		if ( option[ ufv::aa_during_build ].user() ) {
			String aa = core::pose::annotated_to_oneletter_sequence( option[ ufv::aa_during_build ] );
			flags_ok &= aa.length() == option[ ufv::ss ].value().length();

			if ( !flags_ok ) {
				TR.Fatal << "length of ufv:ss string and one letter version of ufv::aa_during_build must be equal!" << std::endl;
			}
		}

		if ( option[ ufv::aa_during_design_refine ].user() ) {
			String aa = core::pose::annotated_to_oneletter_sequence( option[ ufv::aa_during_design_refine ] );
			flags_ok &= aa.length() == option[ ufv::ss ].value().length();

			if ( !flags_ok ) {
				TR.Fatal << "length of ufv:ss string and one letter version of ufv::aa_during_design_refine must be equal!" << std::endl;
			}
		}
	}

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
			TR.Error << "too few columns in line of ufv_loops file: " << line << std::endl;
			std::exit( 1 );
		}
	}

	in.close();

	bdr.use_sequence_bias( use_sequence_biased_fragments );

	return count;
}


void setup_segment_insert( protocols::forge::components::BDR & bdr ) {
	using namespace basic::options::OptionKeys;
	using basic::options::option;
	using core::pose::Pose;
	using core::Size;
	using protocols::forge::build::ConnectRight;
	using protocols::forge::build::CountFromLeft;
	using protocols::forge::build::CountFromLeftOP;
	using protocols::forge::build::Interval;
	using protocols::forge::build::RelativeConnectRight;
	using protocols::forge::build::RelativeConnectRightOP;
	using protocols::forge::build::SegmentInsert;
	using protocols::forge::build::SegmentInsertOP;
	using protocols::forge::build::SegmentInsertConnectionScheme::C;
	using protocols::forge::build::SegmentInsertConnectionScheme::N;
	using protocols::forge::build::SegmentInsertConnectionScheme::RANDOM_SIDE;
	using core::scoring::dssp::Dssp;

	using core::import_pose::pose_from_file;

	// option transmute; just for safety
	String aa_during_build;
	String aa_during_design_refine;
	String connection_scheme_str;

	if ( option[ ufv::aa_during_build ].user() ) {
		aa_during_build = option[ ufv::aa_during_build ];
	}

	if ( option[ ufv::aa_during_design_refine ].user() ) {
		aa_during_design_refine = option[ ufv::aa_during_design_refine ];
	}

	if ( option[ ufv::insert::connection_scheme ].user() ) {
		connection_scheme_str = option[ ufv::insert::connection_scheme ];
		boost::to_upper( connection_scheme_str ); // force uppercase letters
	}

	// Read pdbs.
	// If 'attached' exists, it should already be in proper rigid body
	// orientation wrt to 'insert'.
	bool const use_attached = option[ ufv::insert::attached_pdb ].user();

	Pose insert;
	core::import_pose::pose_from_file( insert, option[ ufv::insert::insert_pdb ] , core::import_pose::PDB_file);
	Dssp dssp_i( insert );
	dssp_i.insert_ss_into_pose( insert );

	RelativeConnectRightOP rcr;

	if ( use_attached ) {
		Pose attached;
		core::import_pose::pose_from_file( attached, option[ ufv::insert::attached_pdb ] , core::import_pose::PDB_file);
		Dssp dssp_a( attached );
		dssp_a.insert_ss_into_pose( attached );

		// Squish together insert <-> attached so we can extract the jump.
		// We attach to the middle of each pdb.  Optimally we should find the
		// closest CAs between the structures, but we're just testing so it's
		// fine for now.
		Pose insert_plus_attached = insert;

		Size const insert_jump_pos = insert.size() % 2 == 0 ? insert.size() / 2 : insert.size() / 2 + 1;
		Size const attached_jump_pos = attached.size() % 2 == 0 ? attached.size() / 2 : attached.size() / 2 + 1;
		Size const shifted_attached_jump_pos = attached_jump_pos + insert.size();

		ConnectRight cr( insert_jump_pos, attached_jump_pos, attached );
		cr.modify( insert_plus_attached );

		// setup RelativeConnectRight
		CountFromLeftOP cfl = new CountFromLeft();
		cfl->p = insert_jump_pos;
		cfl->left_skip = option[ ufv::ss ].value().find( SegmentInsert::insertion_char() );

		rcr = new RelativeConnectRight( cfl, attached_jump_pos, attached );
		rcr->extract_rt( insert_plus_attached, insert_jump_pos, shifted_attached_jump_pos );
		rcr->use_rt( true );
	}

	// setup SegmentInsert
	protocols::forge::build::SegmentInsertConnectionScheme::Enum connection_scheme = RANDOM_SIDE;
	if ( !connection_scheme_str.empty() ) {
		if ( connection_scheme_str == "N2C" ) {
			connection_scheme = N;
		} else if ( connection_scheme_str == "C2N" ) {
			connection_scheme = C;
		}
	}

	SegmentInsertOP si = new SegmentInsert(
		Interval( option[ ufv::left ], option[ ufv::right ] ),
		option[ ufv::ss ], aa_during_build,
		insert,
		false, // force false for now
		connection_scheme
	);

	// add to BDR
	bdr.add_instruction( si, aa_during_design_refine );
	if ( use_attached ) {
		bdr.add_instruction( rcr );
		bdr.create_directed_dependency( si, rcr );
	}
}


void * graphics_main( void * ) {
	using namespace basic::options::OptionKeys;
	using core::Size;
	using basic::options::option;
	using protocols::simple_filters::PoseMetricEvaluator;
	using protocols::forge::build::Interval;
	using protocols::forge::build::SegmentRebuild;
	using protocols::forge::components::BDR;
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

	// init BDR
	BDR_OP bdr = new BDR;
	if ( option[ ufv::ufv_loops ].user() ) {

		load_loops_from_file( *bdr, option[ ufv::ufv_loops ] );

	} else {

		if ( option[ ufv::insert::insert_pdb ].user() ) { // SegmentInsert

			setup_segment_insert( *bdr );

		} else { // SegmentRebuild

			bdr->add_instruction(
				new SegmentRebuild(
				Interval( option[ ufv::left ], option[ ufv::right ] ),
				option[ ufv::ss ], aa_during_build,
				core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ),
				option[ ufv::keep_junction_torsions ]
				),
				aa_during_design_refine
			);

		}

		bdr->use_sequence_bias( !aa_during_build.empty() );
	}

	bdr->use_fullmer( option[ ufv::use_fullmer ].value() );
	bdr->centroid_loop_mover_str( option[ ufv::centroid_loop_mover ] );
	bdr->redesign_loop_neighborhood( !option[ ufv::no_neighborhood_design ].value() );
	bdr->dr_cycles( option[ ufv::dr_cycles ] );

	String centroid_sfx = option[ ufv::centroid_sfx ].user() ? option[ ufv::centroid_sfx ].value() : "remodel_cen";
	String centroid_sfx_patch = option[ ufv::centroid_sfx_patch ].user() ? option[ ufv::centroid_sfx_patch ].value() : String();
	// APL new default -- Talaris2013
	String fullatom_sfx = option[ ufv::fullatom_sfx ].user() ? option[ ufv::fullatom_sfx ].value() : TALARIS_2013;
	String fullatom_sfx_patch = option[ ufv::fullatom_sfx_patch ].user() ? option[ ufv::fullatom_sfx_patch ].value() : String();

	bdr->centroid_scorefunction( core::scoring::ScoreFunctionFactory::create_score_function( centroid_sfx, centroid_sfx_patch ) );
	bdr->fullatom_scorefunction( core::scoring::ScoreFunctionFactory::create_score_function( fullatom_sfx, fullatom_sfx_patch ) );

	// resfile only if requested
	if ( option[ packing::resfile ].user() ) {
		bdr->resfile( option[ packing::resfile ].value().at( 1 ) );
	}

#if defined GL_GRAPHICS
	core::pose::Pose pose;
	core::import_pose::pose_from_file( pose, option[ in::file::s ].value().at( 1 ) , core::import_pose::PDB_file);
	protocols::viewer::add_conformation_viewer( pose.conformation(), "ufv" );
	bdr->apply( pose );
#else
	// setup evaluators (currently only works with silent file output?)
	JobDistributor::get_instance()->job_outputter()->add_evaluation(
		new PoseMetricEvaluator< Size >( BDR::loops_buns_polar_calc_name(), "special_region_bur_unsat_polars" )
	);

	JobDistributor::get_instance()->job_outputter()->add_evaluation(
		new PoseMetricEvaluator< Size >( BDR::neighborhood_buns_polar_calc_name(), "special_region_bur_unsat_polars" )
	);

	// run job
	JobDistributor::get_instance()->go( bdr );
#endif

	return 0;
}


int main( int argc, char * argv [] ) {
	try {
		using namespace basic::options::OptionKeys;
		using basic::options::option;

		// track options
		KeyVec required_options;
		fill_required_options( required_options );

		KeyVec optional_options;
		fill_optional_options( optional_options );

		// register options so help file appears correctly
		register_options( required_options );
		register_options( optional_options );

		// initialize rosetta
		devel::init( argc, argv );

		// check required options are specified
		if ( !check_required_options( required_options ) ) {
			return 1;
		}

		// check option conflicts
		if ( !check_option_conflicts() ) {
			return 1;
		}

		protocols::viewer::viewer_main( graphics_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
