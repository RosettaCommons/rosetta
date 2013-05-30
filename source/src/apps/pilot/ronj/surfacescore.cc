// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file surfaceE_protocols.cc
/// @brief A suite of protocols that do various things using the new surface energy term.
/// @author Ron Jacak


// Unit headers
#include <devel/init.hh>

//project Headers
#include <core/io/pdb/pose_io.hh>
#include <basic/options/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/TenANeighborGraph.hh>

#include <protocols/jobdist/standard_mains.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/toolbox/pose_metric_calculators/SurfaceCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <basic/prof.hh>
#include <utility/file/file_sys_util.hh>
#include <ObjexxFCL/format.hh>

// Numeric Headers

// ObjexxFCL Headers

// C++ headers
#include <sstream>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>



static basic::Tracer TR("surface_app");

using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::pack::task::operation;
using namespace core::pack::task;
using namespace ObjexxFCL::fmt;

// application specific options
namespace surface_app {
	BooleanOptionKey const surfacescore( "surface_app::surfacescore" );
	BooleanOptionKey const surfacedesign( "surface_app::surfacedesign" );
	BooleanOptionKey const surfacescan( "surface_app::surfacescan" );
	BooleanOptionKey const no_repack_before_scoring( "surface_app::no_repack_before_scoring" );
	BooleanOptionKey const no_repack_before_design( "surface_app::no_repack_before_design" );
	BooleanOptionKey const output_by_res_surface_score( "surface_app::output_by_res_surface_score" );
	BooleanOptionKey const repack_mutant_position_only( "surface_app::repack_mutant_position_only" );
	BooleanOptionKey const use_reweighted_score12_with_surfaceE_scorefunction( "surface_app::use_reweighted_score12_with_surfaceE_scorefunction" );
	BooleanOptionKey const use_reweighted_score12_score_function( "surface_app::use_reweighted_score12_score_function" );
	StringOptionKey const native( "surface_app::native" );
}


std::string usage_string;

void init_usage_prompt( std::string exe ) {

	// place the prompt up here so that it gets updated easily; global this way, but that's ok
	std::stringstream usage_stream;
	usage_stream << "No files given: Use either -file:s or -file:l to designate a single pdb or a list of pdbs.\n\n"
			<< "Usage: " << exe
			<< "\n\t-surfacescore | -surfacedesign | -surfacescan"
			<< "\n\t-database path/to/minidb"
			<< "\n\t-s pdb|-l pdbs"
			<< "\n\t-ignore_unrecognized_res"
			<< "\n\t-use_reweighted_score12_with_surfaceE_scorefunction - Use a score12 scorefunction that has the surface term and reference energies, with everything reweighted."
			<< "\n\t-use_reweighted_score12_score_function - Use a reweighted score12 energy function - does not include the surface term."
			<< "\n\t[-output_by_res_surface_score] - Include a printout of the surface score for each residue in the final output."
			<< "\n\t[-ex1 [-ex2]]"
			<< "\n\t[-mute core.io core.conformation]"

			<< "\n\n\t\"-surfacescore\" specific options"
			<< "\n\t[-no_repack_before_scoring]"
			<< "\n\t[-native native.pdb]"

			<< "\n\n\t\"-surfacedesign\" specific options"
			<< "\n\t[-no_repack_before_design]"
			<< "\n\t[-resfile resfile.designall]"
			<< "\n\t[-linmem_ig 20]"
			<< "\n\t[-ndruns 5]"

			<< "\n\n\t\"-surfacescan\" specific options"
			<< "\n\t[-repack_mutant_position_only] - When scanning for mutations that improve the surface score, don't repack the protein."

			<< "\n\n";
	usage_string = usage_stream.str();

}

///
/// @begin print_energies
///
/// @brief
/// Helper method for the main function. Takes in a pose, the scorefunction, surface energies and weights and prints everything
/// out in a pretty format.
///
void print_energies( pose::Pose & pose, scoring::ScoreFunctionOP scorefxn ) {

	scoring::EnergyMap const & wts( scorefxn->weights() );

	std::cout << "res  ";
	// for each energy term, print the weighted energy
	for ( int jj = 1; jj <= scoring::n_score_types; ++jj ) {
		Real const weight = wts[ scoring::ScoreType(jj) ];

		switch( scoring::ScoreType( jj ) ) {
		case scoring::fa_atr:
		case scoring::fa_rep:
		case scoring::fa_sol:
		case scoring::fa_intra_rep:
		case scoring::pro_close:
		case scoring::fa_pair:
		case scoring::hbond_sr_bb:
		case scoring::hbond_lr_bb:
		case scoring::hbond_bb_sc:
		case scoring::hbond_sc:
		case scoring::rama:
		case scoring::omega:
		case scoring::fa_dun:
		case scoring::p_aa_pp:
		case scoring::ref:
		case scoring::unfolded:
		case scoring::surface:
			if ( weight != 0.0 ) {
				std::cout << scoring::ScoreType(jj) << "  ";
			}
			break;
		default:
			break;
		}
	}
	std::cout << std::endl;


	// n_score_types comes in with ScoreType.hh
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {

		// print the string res:
		std::cout << I(2,ii) << "  ";
		scoring::EnergyMap const & unweighted_scores( pose.energies().residue_total_energies( ii ) );

		// for each energy term, print the weighted energy
		for ( int jj = 1; jj <= scoring::n_score_types; ++jj ) {
			Real const weight = wts[ scoring::ScoreType(jj) ];

			switch( scoring::ScoreType( jj ) ) {
			case scoring::fa_atr:
			case scoring::fa_rep:
			case scoring::fa_sol:
			case scoring::fa_intra_rep:
			case scoring::pro_close:
			case scoring::fa_pair:
			case scoring::hbond_sr_bb:
			case scoring::hbond_lr_bb:
			case scoring::hbond_bb_sc:
			case scoring::hbond_sc:
			case scoring::rama:
			case scoring::omega:
			case scoring::fa_dun:
			case scoring::p_aa_pp:
			case scoring::ref:
			case scoring::unfolded:
			case scoring::surface:
				if ( weight != 0.0 ) {
					Real const val = unweighted_scores[ scoring::ScoreType(jj) ];
					std::cout << ObjexxFCL::fmt::F(8,3, weight * val ) << "  ";
				}
				break;
			default:
				break;
			}
		}
		std::cout << std::endl;
	}

}


///@brief helper method which uses the tenA nb graph in the pose object to fill a vector with nb counts
void fill_num_neighbors( pose::Pose & pose, utility::vector1< int > & num_nbs ) {

	// now that it's scored, we can get the number of neighbors each residue has (in the wt structure anyway)
	num_nbs.resize( pose.n_residue(), 0 );
	for ( Size ii=1; ii <= pose.n_residue(); ++ii ) {
		scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
		num_nbs[ ii ] = tenA_neighbor_graph.get_node( ii )->num_neighbors_counting_self();
	}

	return;
}


///@brief helper method which pretty prints the surface scores
void print_surface_score_by_res( pose::Pose & pose, utility::vector1< int > & num_nbs,
		basic::MetricValue< utility::vector1< core::Real > > & pose_residue_surface, Real surface_score_weight ) {

	if ( surface_score_weight == 0.0 ) { std::cout << "Surface energy term not in use." << std::endl; return; }

	Size div2 = pose.n_residue() / 2;
	if ( pose.n_residue() % 2 ) { div2++; }

	for ( Size jj=1; jj < div2; ++jj ) {
		std::cout << "res: " << I(2, jj) << X(3) << "nbs: " << I( 2, num_nbs[jj]) << X(3) << "aa: " << pose.residue(jj).name1() << X(3);
		if ( fabs( pose_residue_surface.value()[jj] - 0.0 ) < 0.001 ) {
			std::cout << "        --- ";
		} else {
			std::cout << F( 12,3, surface_score_weight * pose_residue_surface.value()[jj] );
		}
		if ( jj == div2 && pose.n_residue() % 2 == 1 ) { break; }
		else {
			std::cout << X(5) << "res: " << I(2, jj+div2) << X(3) << "nbs: " << I( 2, num_nbs[jj+div2]) << X(3) << "aa: " << pose.residue(jj+div2).name1() << X(3);
			if ( fabs( pose_residue_surface.value()[jj+div2] - 0.0 ) < 0.001 ) {
				std::cout << "        --- ";
			} else {
				std::cout << F( 12,3, surface_score_weight * pose_residue_surface.value()[jj+div2] );
			}
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


///@brief Takes the input Pose and runs a fast repack protocol on it.
void repack_pose( pose::Pose & pose, scoring::ScoreFunctionOP scorefxn ) {

	pack::task::PackerTaskOP repack_task = pack::task::TaskFactory::create_packer_task( pose );
	//repack_task->set_bump_check( true );
	repack_task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );

	protocols::simple_moves::PackRotamersMoverOP repack_protocol( new protocols::simple_moves::PackRotamersMover( scorefxn, repack_task, 3 /*ndruns value hardcoded*/ ) );
	repack_protocol->apply( pose );

	return;
}


//@brief main method for the surface protocols
int
main( int argc, char* argv[] ) {

	try {


	//
	// add application specific options to options system
	//
	option.add( surface_app::surfacescore, "Score proteins with the surface score and score12." );
	option.add( surface_app::surfacedesign, "Redesign proteins using the score12 score function together with the surface score." );
	option.add( surface_app::surfacescan, "Run a protocol which tries to find mutations that most improve the surface score." );
	option.add( surface_app::no_repack_before_scoring, "Do not repack the pose to remove any steric clashes before scoring." );
	option.add( surface_app::native, "A native structure to use for calculating ddG against." );
	option.add( surface_app::no_repack_before_design, "Do not repack the native pose; clashes in the native pose will not be removed by repacking." );
	option.add( surface_app::output_by_res_surface_score, "Include a printout of the surface score for each residue in the final output." );
	option.add( surface_app::repack_mutant_position_only, "When scanning for mutations that improve the surface score, don't repack the protein." );
	option.add( surface_app::use_reweighted_score12_with_surfaceE_scorefunction, "Use a score12 scorefunction that has the surface term and reference energies, with everything reweighted." );
	option.add( surface_app::use_reweighted_score12_score_function, "Use a reweighted score12 energy function." );


	//
	// options, random initialization
	//
	devel::init( argc, argv );

	//
	// concatenate -s and -l flags together to get total list of PDB files
	// The advantage of parsing -s and -l separately is that users can specify a list and a single structure on the
	// command line.
	//
	using utility::file::file_exists;
	using utility::file::FileName;

	std::vector< FileName > pdb_file_names;
	if ( option[ in::file::s ].active() )
		pdb_file_names = option[ in::file::s ]().vector(); // make a copy (-s)

	std::vector< FileName > list_file_names;
	if ( option[ in::file::l ].active() ) {
		list_file_names = option[ in::file::l ]().vector(); // make a copy (-l)

		//ronj for each file input with the -l switch...
		for ( std::vector< FileName >::iterator i = list_file_names.begin(), i_end = list_file_names.end(); i != i_end; ++i ) {
			std::string listfilename( i->name() );
			std::ifstream data( listfilename.c_str() );
			//ronj try to open that particular file first...
			if ( !data.good() ) {
				utility_exit_with_message( "Unable to open file: " + listfilename + '\n' );
			}
			std::string line;
			//ronj then read all the lines in that file until there are no more
			while( getline(data, line) ) {
				pdb_file_names.push_back( FileName(line) );
			}
			data.close();
		}
	}

	if ( pdb_file_names.size() == 0 ) {
		init_usage_prompt( argv[0] );
		utility_exit_with_message_status( usage_string, 1 );
	}

	//
	// create a surface metric calculator
	//
	pose::metrics::PoseMetricCalculatorOP surface_calculator = new protocols::toolbox::pose_metric_calculators::SurfaceCalculator;
	pose::metrics::CalculatorFactory::Instance().register_calculator( "surface", surface_calculator );

	//
	// create a custom score function as well as a score12 scorefunction
	//
	TR << "Creating score functions." << std::endl;

	scoring::ScoreFunctionOP scorefxn;
	if ( option[ surface_app::use_reweighted_score12_with_surfaceE_scorefunction ] ) {

		TR << "Using reweighted score12 with surfaceE scorefunction (optE 458)." << std::endl;
		scorefxn = scoring::getScoreFunction();

		utility::vector1< Real > refEs_;
		refEs_.resize( chemical::num_canonical_aas, 0.0 );
		Real const rpp_refs[20] = { 0.366095, 0.683858, -0.219389, -0.162503, -0.00261107, 0.508836, 0.427877, -0.437483, -0.12092, -0.0895638,
				0.186069, -0.208825, -0.651193, -0.00180204, -0.0482886, 0.121436, -0.19738, -0.519352, 0.479798, -0.114659 };

		for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
			refEs_[ii] = rpp_refs[ ii-1 ];
		}
		scorefxn->set_method_weights( core::scoring::ref, refEs_ );

		scorefxn->set_weight( core::scoring::fa_atr, 0.496342 );
		scorefxn->set_weight( core::scoring::fa_rep, 0.792435 );
		scorefxn->set_weight( core::scoring::fa_sol, 0.382786 );
		scorefxn->set_weight( core::scoring::fa_intra_rep, 0.062265 );
		scorefxn->set_weight( core::scoring::pro_close, 0.0130755 );
		scorefxn->set_weight( core::scoring::fa_pair, 0.206752 );
		scorefxn->set_weight( core::scoring::hbond_sr_bb, 0.585 );
		scorefxn->set_weight( core::scoring::hbond_lr_bb, 1.17 );
		scorefxn->set_weight( core::scoring::hbond_bb_sc, 0.984159 );
		scorefxn->set_weight( core::scoring::hbond_sc, 0.14829 );
		scorefxn->set_weight( core::scoring::rama, 0.0251727 );
		scorefxn->set_weight( core::scoring::fa_dun, 0.0485666 );
		scorefxn->set_weight( core::scoring::p_aa_pp, 0.469872 );

		scorefxn->set_weight( core::scoring::dslf_ss_dst, 1 );
		scorefxn->set_weight( core::scoring::dslf_cs_ang, 1 );
		scorefxn->set_weight( core::scoring::dslf_ss_dih, 1 );
		scorefxn->set_weight( core::scoring::dslf_ca_dih, 1 );
		scorefxn->set_weight( core::scoring::omega, 0.5 );
		scorefxn->set_weight( core::scoring::ref, 1.0 );
		scorefxn->set_weight( core::scoring::surface, 1.4 );

	} else if ( option[ surface_app::use_reweighted_score12_score_function ] ) {

		TR << "Using reweighted score12 scorefunction." << std::endl;
		scorefxn = scoring::getScoreFunction();

		utility::vector1< Real > refEs_;
		refEs_.resize( chemical::num_canonical_aas, 0.0 );
		Real const rpp_refs[20] = { 1.44647, 3.20455, 1.28033, 1.50335, 2.39638, 1.47398, 2.35098, 1.24274, 1.64901, 1.11337,
				1.68414, 0.895638, 1.68502, 1.43777, 1.28772, 1.26944, 1.03373, 0.903269, 4.10334, 1.8192 };
		for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
			refEs_[ii] = rpp_refs[ ii-1 ];
		}
		scorefxn->set_method_weights( core::scoring::ref, refEs_ );

		scorefxn->set_weight( core::scoring::fa_atr, 0.6 );
		scorefxn->set_weight( core::scoring::fa_rep, 0.676039 );
		scorefxn->set_weight( core::scoring::fa_sol, 0.473486 );
		scorefxn->set_weight( core::scoring::fa_intra_rep, 0.0669174 );
		scorefxn->set_weight( core::scoring::pro_close, 0.0430282 );
		scorefxn->set_weight( core::scoring::fa_pair, 0.155594 );
		scorefxn->set_weight( core::scoring::hbond_sr_bb, 1.34807 );
		scorefxn->set_weight( core::scoring::hbond_lr_bb, 1.79619 );
		scorefxn->set_weight( core::scoring::hbond_bb_sc, 0.855615 );
		scorefxn->set_weight( core::scoring::hbond_sc, 0.34485 );
		scorefxn->set_weight( core::scoring::dslf_ss_dst, 1 );
		scorefxn->set_weight( core::scoring::dslf_cs_ang, 1 );
		scorefxn->set_weight( core::scoring::dslf_ss_dih, 1 );
		scorefxn->set_weight( core::scoring::dslf_ca_dih, 1 );
		scorefxn->set_weight( core::scoring::rama, 0.0407331 );
		scorefxn->set_weight( core::scoring::omega, 0.5 );
		scorefxn->set_weight( core::scoring::fa_dun, 0.0647276 );
		scorefxn->set_weight( core::scoring::p_aa_pp, 0.319602 );
		scorefxn->set_weight( core::scoring::ref, 1.0 );

	} else {

		TR << "Using standard score12 scorefunction." << std::endl;
		scorefxn = scoring::getScoreFunction();

	}

	// used by the surface calculators
	utility::vector1< int > num_neighbors_;

	// if we're in score mode, we might need a native pose handle
	pose::Pose native_pose;
	basic::MetricValue< Real > native_pose_surface_total;
	basic::MetricValue< utility::vector1< core::Real > > native_pose_residue_surface;

	if ( option[ surface_app::native ].user() ) {
		core::import_pose::pose_from_pdb( native_pose, option[ surface_app::native ] );
		repack_pose( native_pose, scorefxn );
		native_pose.metric( "surface", "residue_surface", native_pose_residue_surface );
		Energy native_score = (*scorefxn)( native_pose );
	}

	// iterate through all the structures
	for(std::vector< FileName >::iterator input_pdb_filename = pdb_file_names.begin(), last_pdb = pdb_file_names.end();
		input_pdb_filename != last_pdb; ++input_pdb_filename) {

		// check to make sure the file exist; if not, move on to the next one
		if ( !file_exists( *input_pdb_filename ) ) {
			std::cerr << "Input pdb " << *input_pdb_filename << " not found, skipping" << std::endl;
			continue;
		}

		pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, *input_pdb_filename );

		// score the wt pose, so that the metrics can be computed - inefficient but scoring is really fast so...
		TR << "Scoring pose so that metrics can be calculated." << std::endl;
		(*scorefxn)( pose );

		basic::MetricValue< utility::vector1< core::Real > > wt_pose_residue_surface;
		pose.metric( "surface", "residue_surface", wt_pose_residue_surface );

		// now we get to the part where we have to decide what we're actually doing
		// check what "mode" the user ran: was -surfacescore, -surfacedesign, or -surfacescan specified.

		///
		/// surface design
		///
		if ( option[ surface_app::surfacedesign ] ) {

			clock_t starttime = clock();
			basic::prof_reset();

			if ( !( option[ surface_app::no_repack_before_design ] ) ) {
				repack_pose( pose, scorefxn );
			}

			Energy wt_score = (*scorefxn)( pose );

			// print out the scorefxn information...
			std::cout << "wt score:" << std::endl;
			scorefxn->show( std::cout, pose );
			std::cout << std::endl;

			fill_num_neighbors( pose, num_neighbors_ );

			//print_energies( pose, scorefxn );

			// set up a vector of strings for the total surface and a vector of vectors for the per res surfaces
			Energy design_total_energy;
			Real design_surface;
			basic::MetricValue< utility::vector1< core::Real > > residue_surfaceEs;
			scoring::EnergyMap design_pd_scores;

			// create a custom packer task, for use in creation of pack rotamers mover
			// design positions, as specified by the resfile and command line
			pack::task::PackerTaskOP designtask( pack::task::TaskFactory::create_packer_task( pose ));
			designtask->set_bump_check( true );
			designtask->initialize_from_command_line();
			parse_refile(pose, *designtask);

			protocols::simple_moves::PackRotamersMoverOP design_protocol( new protocols::simple_moves::PackRotamersMover( scorefxn, designtask, (Size)basic::options::option[packing::ndruns].value() ) );
			design_protocol->apply( pose );

			Energy design_score = (*scorefxn)( pose );

			// print out the scorefxn information...
			std::cout << "redesign score:" << std::endl;
			scorefxn->show( std::cout, pose );
			std::cout << std::endl;

			design_total_energy = design_score;
			design_pd_scores = pose.energies().total_energies();

			// output designed structures with '.surface' in the name - the ndrun number (e.g. 2) will also be in the name
			std::stringstream out;
			std::string filename = utility::file::file_basename( *input_pdb_filename );

			core::Real surface_score_weight( scorefxn->get_weight( core::scoring::surface ) );
			if ( surface_score_weight ) {
				out << ".surface." << 'w' << surface_score_weight << '.';
			}
			filename = filename + ".designed" + out.str() + ".pdb";
			pose.dump_scored_pdb( filename, *scorefxn );

			basic::MetricValue< utility::vector1< core::Real > > design_surface_by_res;
			pose.metric( "surface", "residue_surface", design_surface_by_res );
			residue_surfaceEs = design_surface_by_res;

			print_energies( pose, scorefxn );

			TR << "Done with packing. Printing summary..." << std::endl;

			if ( option[ surface_app::output_by_res_surface_score ] ) {
				std::cout << A(35, "redesign surface score by res:") << std::endl;
				print_surface_score_by_res( pose, num_neighbors_, residue_surfaceEs, scorefxn->get_weight( scoring::surface ) );
			}

			clock_t stoptime = clock();
			TR << "Whole run took " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " seconds" << std::endl;


		///
		/// surface scan
		///
		} else if ( option[ surface_app::surfacescan ] ) {

			// we need to repack the entire Pose so that when we go into the loop that tries mutations, we don't get big effects
			// from the other score terms.
			repack_pose( pose, scorefxn );

			Energy wt_repacked_score = (*scorefxn)( pose );

			fill_num_neighbors( pose, num_neighbors_ );

			// and then follow with lots of surface information.
			std::cout << "wild type score summary:" << std::endl;
			scorefxn->show( std::cout, pose );
			std::cout << std::endl;
			if ( option[ surface_app::output_by_res_surface_score ] ) {
				std::cout << A( "surface score by res:" ) << std::endl;
				print_surface_score_by_res( pose, num_neighbors_, wt_pose_residue_surface, scorefxn->get_weight( scoring::surface ) );
			}

			// now we are ready to go into the scan loop.  the optE protocol does this
			std::cout << A( "mutation ddG's:" ) << std::endl;
			std::cout << A( "scorefxn energy" ) << X(3) << A( "scorefxn ddG" ) << X(3) << A( "surface ddG" ) << X(3) << A( "mutation" ) << std::endl;


			for ( Size resid = 1; resid <= pose.n_residue(); ++resid ) {
				for ( Size aa_enum_index = 1; aa_enum_index < chemical::num_canonical_aas; ++aa_enum_index ) {

					if ( pose.residue( resid ).aa() == chemical::AA( aa_enum_index ) ) { continue; }

					pose::Pose pose_copy = pose;

					// need to create a calculator here that I can use to identify neighbors of the mutated residue
					pose::metrics::PoseMetricCalculatorOP mutant_nb_calculator = new toolbox::pose_metric_calculators::NeighborsByDistanceCalculator( resid );
					pose::metrics::CalculatorFactory::Instance().register_calculator( "mutant_nb_calculator", mutant_nb_calculator );

					// the restrict operation class (which in the end is just a TaskOperation) takes a calculator during construction. I've already
					// created that calculator above.  This operation will disable repacking and design at all positions except those in the neighborhood
					// of the mutated position.
					TaskOperationCOP nb_op = new toolbox::task_operations::RestrictToNeighborhoodOperation( "mutant_nb_calculator" );

					// extra operations we want to also include
					// the restrict residue to repacking ops are used to make sure that only repacking and not design is done to the residues in the neighborhood
					InitializeFromCommandlineOP init_op = new InitializeFromCommandline();
					IncludeCurrentOP ic_op = new IncludeCurrent();
					RestrictResidueToRepackingOP mutant_repack_op = new RestrictResidueToRepacking();

					TaskFactoryOP tf = new TaskFactory();

					tf->push_back( init_op );
					tf->push_back( ic_op );
					tf->push_back( nb_op );

					for ( Size ii = 1; ii <= pose_copy.n_residue(); ++ii ) {
						if ( ii == resid ) {
							// do design on this position
							utility::vector1< bool > keep_canonical_aas( chemical::num_canonical_aas, false );
							keep_canonical_aas[ aa_enum_index ] = true;
							RestrictAbsentCanonicalAASOP restrict_op = new RestrictAbsentCanonicalAAS( ii, keep_canonical_aas );
							tf->push_back( restrict_op );
						} else {
							// make this position repackable only; because of the commutativity of packer task ops, only the residues that are in the neighborhood
							// of the mutant will be allowed to repack. the restrict to neighborhood op will disallow packing at all positions not near the mutant.
							mutant_repack_op->include_residue( ii );
						}
					}

					tf->push_back( mutant_repack_op );

					pack::task::PackerTaskOP scan_task = tf->create_task_and_apply_taskoperations( pose_copy );
					scan_task->num_to_be_packed();
					//std::cout << *scan_task << std::endl;  // generates a TON of output

					protocols::simple_moves::PackRotamersMoverOP mutant_repack( new protocols::simple_moves::PackRotamersMover( scorefxn, scan_task, 2 /*ndruns*/) );
					mutant_repack->apply( pose_copy );

					Energy mutant_repacked_score = (*scorefxn)( pose_copy );

					// print out scorefunction information, if desired
					//scorefxn->show( std::cout, pose_copy );
					//std::cout << std::endl;

					Real ddG_mutation = mutant_repacked_score - wt_repacked_score;

					// don't bother printing out destabilizing mutations even if they really improve the surface score
					if ( ddG_mutation > 0.0 ) {
						pose::metrics::CalculatorFactory::Instance().clear_calculators();
						continue;
					}
					Real mutant_surfaceE = pose_copy.energies().total_energies()[ scoring::surface ];
					Real wt_surfaceE = pose.energies().total_energies()[ scoring::surface ];
					Real surface_ddG_mutation = mutant_surfaceE - wt_surfaceE;

					std::stringstream out;
					out << pose.residue( resid ).name1() << resid << oneletter_code_from_aa( chemical::AA( aa_enum_index ) ) ;
					std::string mutation_string = out.str();

					if ( scorefxn->get_weight( scoring::surface ) != 0.0 ) {
						std::cout << F( 9,2,mutant_repacked_score ) << X(3) << F( 9,2,ddG_mutation ) << X(3) << F( 9,2,surface_ddG_mutation ) << X(3) << mutation_string << std::endl;
					} else {
						std::cout << F( 9,2,mutant_repacked_score ) << X(3) << F( 9,2,ddG_mutation ) << X(3) << "---" << X(3) << mutation_string << std::endl;
					}

					// this needs to get recreated each time around
					pose::metrics::CalculatorFactory::Instance().clear_calculators();


				} // over aas
			} // over resid

		///
		/// surface score
		///
		} else /* default behaviour is to score */ {

			if ( !( option[ surface_app::no_repack_before_scoring ] ) ) {
				repack_pose( pose, scorefxn );
			}

			Energy score = (*scorefxn)( pose );
			fill_num_neighbors( pose, num_neighbors_ );

			// print out the scorefxn information...
			scorefxn->show( std::cout, pose );
			std::cout << std::endl;
			if ( option[ surface_app::output_by_res_surface_score ] ) {
				std::cout << A( "surface score by res:" ) << std::endl;
				print_surface_score_by_res( pose, num_neighbors_, wt_pose_residue_surface, scorefxn->get_weight( scoring::surface ) );
			}

		}


	} // now do it all over for another structure


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}
