// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/public/remodel.cc
/// @author Yih-En Andrew Ban (yab@u.washington.edu)
/// @author Possu Huang

// project headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/ufv.OptionKeys.gen.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// utility headers
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/options/keys/OptionKey.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/forge/build/ConnectRight.hh>
#include <protocols/forge/build/RelativeConnectRight.hh>
#include <protocols/forge/build/RelativeSequencePosition.hh>
#include <protocols/forge/build/SegmentInsert.hh>
#include <protocols/forge/build/SegmentRebuild.hh>
#include <protocols/forge/components/BDR.hh>
#include <protocols/forge/remodel/RemodelMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/simple_filters/PoseMetricEvaluator.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/PyMolMover.hh>

#include <devel/init.hh>

// boost headers
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

// C++ headers
#include <string>
#include <vector>

#include <utility/excn/Exceptions.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;

typedef std::string String;
typedef std::vector< OptionKey const * > KeyVec;

static THREAD_LOCAL basic::Tracer TR( "apps.public.remodel" );


std::string usage_string;


/// @brief
/// the usage prompt that gets printed when the user doesn't enter any arguments or uses the -h flag since the application
/// specific help Rosetta bring up with the -help flag is, in fact, not helpful.
///
void init_usage_prompt( std::string exe ) {

	// place the prompt up here so that it gets updated easily; global this way, but that's ok
	std::stringstream usage_stream;
	usage_stream
		<< "Usage: " << exe
		<< "\n"
		<< "\n\t-database path/to/rosetta/database"
		<< "\n\t-s pdb | -silent silent_file                        read in input PDB file or silent file"
		<< "\n\t-pdb_gz | -silent_gz                                compress output PDB or silent files"
		<< "\n"
		<< "\n\t-remodel:blueprint blueprint_file                   blueprint file which explains how to run remodel"
		<< "\n"
		<< "\n\t-remodel::num_trajectory                            number of build trajectories to run (default: 10)"
		<< "\n\t-remodel::dr_cycles <int>                           number of design/refine cycles to run (default: 3)"
		<< "\n"
		<< "\n\t-remodel::quick_and_dirty                           only do fragment sampling; bypass refinement of final structures which is slow. useful in early stages of design when one wants to sample different loop lengths to find appropriate setup (default: false)"
		<< "\n"
		<< "\n\t-remodel::bypass_fragments                          skip creation of fragments for remodelling; do refinement only. no extensions or deletions are honored in the blueprint (default: false)"
		<< "\n\t-remodel::use_blueprint_sequence                    find fragments which have the same amino acid sequence and secondary structure as what is specified in the second column of the blueprint file (default: false)"
		<< "\n\t-remodel::use_same_length_fragments                 harvest fragments that match the length of the segment being rebuilt (default: true)"
		<< "\n"
		<< "\n\t-remodel::build_disulf                              use Remodel to find residue pairs - between the \"build\" and \"landing\" residues - which would make good disulfide bonds (default: false)"
		<< "\n\t-remodel::match_rt_limit <float>                    the score cutoff to use for determining how closely a potential disulfide must match observed disulfide distributions (default: 0.4)"
		<< "\n\t-remodel::disulf_landing_range <range>              the range within which Remodel attempts to find disulfides for positions in the \"build\" region (default: none)"
		<< "\n"
		<< "\n\t-remodel::use_pose_relax                            add fast relax to the refinement stage (instead of the default minimization step), but use constraints in a similar way (default: false)"
		<< "\n\t-remodel::run_confirmation                          use kinematic loop closure algorithm for build confirmation (default: false)"
		<< "\n\t-remodel::swap_refine_confirm_protocols             swap protocols used for refinement and confirmation test; i.e. use kinematic loop closure instead of CCD closure (default: false)"
		<< "\n\t-remodel::repeat_structure                          build identical repeats this many times (default: 1)"
		<< "\n"
		<< "\n\t-symmetry::symmetry_definition                      text file describing symmetry setup (default: none)"
		<< "\n"
		<< "\n\t-enzdes::cstfile                                    enzyme design constraints file"
		<< "\n\t-remodel::cstfilter                                 threshold to put on the atom_pair_constraint score type filter during the centroid build phase refinement (default: 10)"
		<< "\n"
		<< "\n\t-remodel::domainFusion::insert_segment_from_pdb     segment PDB file to be inserted into the input structure"
		<< "\n"
		<< "\n\t-remodel::checkpoint                                turns on checkpointing, for use in preemptive scheduling environments. writes out the best pdbs collected after each design step. (default: false)"
		<< "\n\t-remodel::use_clusters                              specifies whether to perform clustering during structure aggregation (default: false)"
		<< "\n\t-remodel::save_top                                  the number of final lowest scoring pdbs to keep (default: 5)"
		<< "\n"
		<< "\n\t-remodel::generic_aa <letter>                       residue type to use as placeholder during centroid phase (default: V)"
		<< "\n\t-remodel::cen_sfxn                                  score function to be used for centroid phase building (default: remodel_cen)"
		<< "\n\t-remodel::cen_minimize                              centroid minimization after fragment building (default: false)"
		<< "\n"
		<< "\n\t-run::chain <letter>                                chain id of the chain to remodel, if a multichain structure is given (default: -)"
		<< "\n"
		<< "\n\t[-overwrite]"
		<< "\n\t[-ignore_unrecognized_res]"
		<< "\n\t[-mute core.io core.conformation core.pack core.scoring]"

		<< "\n\nPlease see the Rosetta Remodel documentation in the Rosetta Manual for more information, including a list of all available options."
		<< "\nhttp://www.rosettacommons.org/manual_guide"

		<< "\n\n";
	usage_string = usage_stream.str();

}


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
	using protocols::forge::build::BuildInstructionOP;
	using protocols::forge::build::Interval;
	using protocols::forge::build::SegmentRebuild;
	using std::istringstream;

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
				BuildInstructionOP( new SegmentRebuild(
				Interval( left, right ),
				ss, aa_during_build,
				core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ),
				option[ ufv::keep_junction_torsions ]
				) ),
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

		Size const insert_jump_pos = insert.n_residue() % 2 == 0 ? insert.n_residue() / 2 : insert.n_residue() / 2 + 1;
		Size const attached_jump_pos = attached.n_residue() % 2 == 0 ? attached.n_residue() / 2 : attached.n_residue() / 2 + 1;
		Size const shifted_attached_jump_pos = attached_jump_pos + insert.n_residue();

		ConnectRight cr( insert_jump_pos, attached_jump_pos, attached );
		cr.modify( insert_plus_attached );

		// setup RelativeConnectRight
		CountFromLeftOP cfl( new CountFromLeft() );
		cfl->p = insert_jump_pos;
		cfl->left_skip = option[ ufv::ss ].value().find( SegmentInsert::insertion_char() );

		rcr = RelativeConnectRightOP( new RelativeConnectRight( cfl, attached_jump_pos, attached ) );
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

	SegmentInsertOP si( new SegmentInsert(
		Interval( option[ ufv::left ], option[ ufv::right ] ),
		option[ ufv::ss ], aa_during_build,
		insert,
		false, // force false for now
		connection_scheme
		) );

	// add to BDR
	bdr.add_instruction( si, aa_during_design_refine );
	if ( use_attached ) {
		bdr.add_instruction( rcr );
		bdr.create_directed_dependency( si, rcr );
	}
}


void* graphics_main( void* ) {

	using namespace basic::options::OptionKeys;
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
	/*String aa_during_build;
	String aa_during_design_refine;

	if ( option[ ufv::aa_during_build ].user() ) {
	aa_during_build = option[ ufv::aa_during_build ];
	}

	if ( option[ ufv::aa_during_design_refine ].user() ) {
	aa_during_design_refine = option[ ufv::aa_during_design_refine ];
	}
	*/
	// init BDR
	RemodelMover_OP rmdl( new RemodelMover );

	// run job
	JobDistributor::get_instance()->go( rmdl );

	protocols::viewer::clear_conformation_viewers();
	exit( 0 ); // ensures graceful exit of graphics.
}


int main( int argc, char * argv [] ) {
	try {
		// check to see if no flags or the -h flag were specified.
		if ( argc == 1 || ( argc > 1 && strcmp(argv[ 1 ], "-h") == 0 ) ) {
			init_usage_prompt( argv[0] );
			std::cout << usage_string;
			exit(0);
		}

		KeyVec optional_options;
		fill_optional_options( optional_options );

		// register options so help file appears correctly
		register_options( optional_options );

		// initialize rosetta
		devel::init( argc, argv );

		// viewer_main() just calls graphics_main with the parameter NULL and returns 0
		protocols::viewer::viewer_main( graphics_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}


