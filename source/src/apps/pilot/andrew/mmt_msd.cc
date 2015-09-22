// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/andrew/apl_msd.cc
/// @brief  Massively multithreaded multistate design protocol
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef USEMPI
#include <mpi.h>
#endif

// MMT MSD Headers
#include <devel/mmt_msd/MMTReceiver.hh>
#include <devel/mmt_msd/MMTDriver.hh>

/// Pack Daemon headers
//#include <protocols/pack_daemon/EntityCorrespondence.hh>
//#include <protocols/pack_daemon/DynamicAggregateFunction.hh>
//#include <protocols/pack_daemon/MultistateAggregateFunction.hh>
//#include <protocols/pack_daemon/MultistateFitnessFunction.hh>
//#include <protocols/pack_daemon/PackDaemon.hh>
//#include <protocols/pack_daemon/util.hh>

/// Core headers
#include <devel/init.hh>
//#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/MetricValue.hh>

//#include <core/chemical/ResidueType.hh>
//#include <core/conformation/Residue.hh>
//#include <core/pose/Pose.hh>
//#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>
//
//#include <core/scoring/hbonds/HBondOptions.hh>
//#include <core/scoring/methods/EnergyMethodOptions.hh>

////#include <core/pack/task/PackerTask.hh>
////#include <core/pack/interaction_graph/SurfacePotential.hh>
//
//#include <core/import_pose/import_pose.hh>

#include <basic/Tracer.hh>

// Protocols headers
//#include <protocols/genetic_algorithm/Entity.hh>
//#include <protocols/genetic_algorithm/EntityRandomizer.hh>
//#include <protocols/genetic_algorithm/Mutate1Randomizer.hh>
//#include <protocols/genetic_algorithm/GeneticAlgorithm.hh>
//#include <protocols/multistate_design/util.hh>

//#include <protocols/toolbox/pose_metric_calculators/HPatchCalculator.hh>
//#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
//#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>

// Devel headers
//#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>

// Utility headers
#include <utility/mpi_util.hh>
//#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

// Numeric headers
//#include <numeric/numeric.functions.hh>
//#include <numeric/random/random.hh>

// option key includes
#include <basic/options/keys/ms.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>


#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


OPT_1GRP_KEY( String, msd, entity_resfile )
OPT_1GRP_KEY( String, msd, fitness_file )
OPT_1GRP_KEY( StringVector, msd, seed_sequences )
//OPT_1GRP_KEY( Integer, msd, double_lazy_ig_mem_limit )
OPT_1GRP_KEY( Boolean, msd, dont_score_bbhbonds )
OPT_1GRP_KEY( Boolean, msd, exclude_background_energies )
OPT_1GRP_KEY( String, msd, seed_sequence_from_input_pdb )
OPT_1GRP_KEY( String, msd, seed_sequence_using_correspondence_file )
OPT_1GRP_KEY( Integer, msd, n_worker_threads_per_process )

/*  Option( 'double_lazy_ig_mem_limit', 'Integer',
desc="The amount of memory, in MB, that each double-lazy interaction graph should be allowed \
to allocate toward rotamer pair energies.  Using this flag will not trigger the \
use of the double-lazy interaction graph, and this flag is not read in the PackerTask's \
initialize_from_command_line routine.  For use in multistate design",
default='0',
),*/


using basic::t_info;
using basic::t_debug;
static THREAD_LOCAL basic::Tracer TR( "apps.public.design.mpi_msd", t_info );

/*class SimpleDGBindAggregateFunction : public protocols::pack_daemon::MultistateAggregateFunction
{
public:
typedef protocols::pack_daemon::MultistateAggregateFunction parent;

public:
SimpleDGBindAggregateFunction() : parent() {}
virtual ~SimpleDGBindAggregateFunction() {}

virtual core::Real evaluate( utility::vector1< core::Real > const & vec, Entity const &  ) {
if ( vec.size() != 2 ) {
utility_exit_with_message( "SimpleDGBindAggregateFunction expects exactly 2 states" );
}
return vec[ 1 ] - vec[ 2 ];
}

virtual StateIndices select_relevant_states( StateEnergies const & ) {
StateIndices two_vec;
two_vec.push_back( 1 );
two_vec.push_back( 2 );
return two_vec;
}
};*/


//protocols::genetic_algorithm::EntityElements
//entity_elements_from_1letterstring(
// std::string const & input
//)
//{
// protocols::genetic_algorithm::EntityElements elements( input.size() );
// for ( core::Size ii = 0, count = 1; ii < input.size(); ++ii, ++count ) {
//  std::ostringstream output;
//  output << "AA:" << count << ":" << input[ ii ];
//  elements[ count ] = protocols::genetic_algorithm::EntityElementFactory::get_instance()->element_from_string( output.str() );
// }
// return elements;
//}

//std::string
//read_native_sequence_for_entity_elements( core::Size n_designed_positions )
//{
// using namespace basic::options;
// using namespace basic::options::OptionKeys;
// using namespace core;
// using namespace protocols::pack_daemon;
//
// if ( ! option[ msd::seed_sequence_using_correspondence_file ].user() ) {
//  utility_exit_with_message( "Must provide correspondence file to read the native sequence" );
// }
//
// std::string pdb_name = option[ msd::seed_sequence_from_input_pdb ];
// std::string correspondence_file_name = option[ msd::seed_sequence_using_correspondence_file ];
//
// /// Read in the pdb
// pose::Pose pose;
// import_pose::pose_from_pdb( pose, pdb_name );
//
// utility::io::izstream correspondence_file( correspondence_file_name );
// if ( ! correspondence_file ) {
//  utility_exit_with_message( "Could not open correspondence file named: " + correspondence_file_name );
// }
//
// EntityCorrespondenceOP ec = new EntityCorrespondence;
// ec->set_pose( new pose::Pose( pose ));
// ec->set_num_entities( n_designed_positions );
// ec->initialize_from_correspondence_file( correspondence_file );
//
//
// std::map< Size, chemical::AA > aa_for_design_position;
// for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
//  Size ii_entity = ec->entity_for_residue( ii );
//  if ( ii_entity  == 0 ) continue;
//  if ( aa_for_design_position.find( ii_entity ) != aa_for_design_position.end() ) {
//   utility_exit_with_message( "Repeat entity element for native pdb: " + pdb_name + " with correspondence file " +
//    correspondence_file_name + ".  Entity correspondence file should only include each residue once");
//  }
//  aa_for_design_position[ ii_entity ] = pose.residue_type( ii ).aa();
// }
// std::string aa_string( n_designed_positions, 'X' );
// for ( Size ii = 1; ii <= n_designed_positions; ++ii ) {
//  if ( aa_for_design_position.find( ii ) == aa_for_design_position.end() ) {
//   utility_exit_with_message( "Did not find residue assigned to correspond to entity element " +
//    utility::to_string( ii ) + " while reading correspondence file " + correspondence_file_name );
//  }
//  aa_string[ ii-1 ] = oneletter_code_from_aa( aa_for_design_position[ ii ] );
// }
// return aa_string;
//
//
//}

///////////////////////////////////////////////////////////////////////////////
///////////////////////   HPatchNPDCalculator   ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//class HPatchNPDCalculator : public protocols::pack_daemon::NPDPropCalculator
//{
//public:
//
// virtual
// core::Real
// calculate( core::pose::Pose const & p ) {
//  return core::pack::interaction_graph::SurfacePotential::get_instance()->compute_pose_hpatch_score( p );
// }
//
//};
//
//class HPatchNPDCalculatorCreator : public protocols::pack_daemon::NPDPropCalculatorCreator
//{
// virtual
// std::string
// calculator_name() const {return "hpatch"; }
//
// virtual
// protocols::pack_daemon::NPDPropCalculatorOP
// new_calculator() const { return new HPatchNPDCalculator; }
//};

///////////////////////////////////////////////////////////////////////////////
///////////////////////   HPatchByChainNPDCalculator   ////////////////////////
///////////////////////////////////////////////////////////////////////////////


//class HPatchByChainNPDCalculator : public protocols::pack_daemon::NPDPropCalculator
//{
//public:
// virtual
// void
// setup(
//  core::pose::Pose const & pose,
//  core::pack::task::PackerTask const & task
// ){
//  chains_ = pose.split_by_chain();
//  task_ = task.clone();
//  Size last_chain( 1 ), first_residue_for_chain( 1 );
//  resid_2_chain_and_resid_.resize( pose.total_residue() );
//  for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
//   if ( pose.residue( ii ).chain() != last_chain ) {
//    last_chain = pose.residue( ii ).chain();
//    first_residue_for_chain = ii;
//   }
//   resid_2_chain_and_resid_[ ii ] = std::make_pair( pose.residue( ii ).chain(), ii + 1 - first_residue_for_chain );
//  }
// }
//
// virtual
// core::Real
// calculate( core::pose::Pose const & p ) {
//  // MJO COMMENTING OUT BECAUSE IT IS UNUSED:
//  // Size chain_offset = 0;
//  for ( core::Size ii = 1; ii <= p.total_residue(); ++ii ) {
//   if ( ! task_->being_packed( ii ) ) continue;
//   chains_[ resid_2_chain_and_resid_[ ii ].first ]->replace_residue(
//    resid_2_chain_and_resid_[ ii ].second, p.residue( ii ), false );
//  }
//
//  core::Real hpatch_sum = 0;
//  for ( Size ii = 1; ii <= chains_.size(); ++ii ) {
//   hpatch_sum += core::pack::interaction_graph::SurfacePotential::get_instance()->compute_pose_hpatch_score( *chains_[ ii ] );
//  }
//  return hpatch_sum;
// }
//private:
// utility::vector1< core::pose::PoseOP > chains_;
// core::pack::task::PackerTaskOP task_;
// utility::vector1< std::pair< core::Size, core::Size > > resid_2_chain_and_resid_;
//};
//
//class HPatchByChainNPDCalculatorCreator : public protocols::pack_daemon::NPDPropCalculatorCreator
//{
// virtual
// std::string
// calculator_name() const {return "hpatch_by_chain"; }
//
// virtual
// protocols::pack_daemon::NPDPropCalculatorOP
// new_calculator() const { return new HPatchByChainNPDCalculator; }
//};


///////////////////////////////////////////////////////////////////////////////
///////////////////////      BunsCalculator     ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//class NBuriedUnsatsCalcultor : public protocols::pack_daemon::NPDPropCalculator
//{
//public:
// virtual
// void
// setup(
//  core::pose::Pose const & pose,
//  core::pack::task::PackerTask const & task
// ){
//  pose_ = new core::pose::Pose( pose );
//  task_ = task.clone();
//  sfxn_ = core::scoring::get_score_function();
//
//  //register calculators
//  if ( ! core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "sasa" ) ) {
//   core::pose::metrics::PoseMetricCalculatorOP sasa_calculator = new devel::vardist_solaccess::VarSolDistSasaCalculator;
//   core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );
//  }
//
//  if ( ! core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "num_hbonds" ) ) {
//   core::pose::metrics::PoseMetricCalculatorOP num_hbonds_calculator = new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator();
//   core::pose::metrics::CalculatorFactory::Instance().register_calculator( "num_hbonds", num_hbonds_calculator );
//  }
//
//  if ( ! core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "unsat" ) ) {
//   core::pose::metrics::PoseMetricCalculatorOP unsat_calculator = new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("sasa", "num_hbonds");
//   core::pose::metrics::CalculatorFactory::Instance().register_calculator( "unsat", unsat_calculator );
//  }
//
// }
//
// virtual
// core::Real
// calculate( core::pose::Pose const & p ) {
//  for ( core::Size ii = 1; ii <= pose_->total_residue(); ++ii ) {
//   if ( ! task_->being_packed( ii ) ) continue;
//   pose_->replace_residue( ii, p.residue( ii ), false );
//  }
//  (*sfxn_)(*pose_);
//  basic::MetricValue< core::Size > nburied_unsats;
//  pose_->metric( "unsat", "all_bur_unsat_polars", nburied_unsats );
//  return nburied_unsats.value();
// }
//
//private:
// core::pose::PoseOP             pose_;
// core::pack::task::PackerTaskOP task_;
// core::scoring::ScoreFunctionOP sfxn_;
//};
//
//class NBuriedUnsatsCalcultorCreator : public protocols::pack_daemon::NPDPropCalculatorCreator
//{
// virtual
// std::string
// calculator_name() const {return "nbunsats"; }
//
// virtual
// protocols::pack_daemon::NPDPropCalculatorOP
// new_calculator() const { return new NBuriedUnsatsCalcultor; }
//};


///////////////////////////////////////////////////////////////////////////////
///////////////////////           main          ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////


int main( int argc, char ** argv )
{
	using namespace utility;
	using namespace protocols::pack_daemon;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::genetic_algorithm;

	NEW_OPT( msd::entity_resfile, "Resfile for the entity elements which are shared between the multiple states", "" );
	NEW_OPT( msd::fitness_file,   "Fitness function file specifying the multiple states and the objective function to optimize", "" );
	NEW_OPT( msd::dont_score_bbhbonds, "Disable the default activation of the decompose_bb_hb_into_pair_energies flag for hbonds", false );
	NEW_OPT( msd::seed_sequences, "Seed the GA's population with the given input sequences", "" );
	NEW_OPT( msd::seed_sequence_from_input_pdb, "Seed the GA's population with the given input sequences using the sequence already present in a given pdb file; requires the use of the msd::seed_sequence_using_correspondence_file flag ", "" );
	NEW_OPT( msd::seed_sequence_using_correspondence_file, "The name of the correspondence file to guide the seeding of the GA's population with the sequence from a particular pdb", "" );
	NEW_OPT( msd::n_worker_threads_per_process, "The number of worker threads that each MPI process should use", 1 );

	try {


		devel::init( argc, argv );

		//std::string secondary_resfile( "2wo2_secondary.resfile" );
		//std::string bound_pdb( "2wo2.pdb" );
		//std::string unbound_pdb( "2wo2_sep.pdb" );
		//std::string entity_correspondence_file( "entity_map.txt" );
		//std::string daf_filename = "objective_function.daf";

		if ( ! option[ msd::entity_resfile ].user() ) {
			utility_exit_with_message("The entity resfile must be specified for the mpi_msd application with the -msd::entity_resfile <filename> flag" );
		}

		if ( ! option[ msd::fitness_file ].user() ) {
			utility_exit_with_message("The fitness-function file must be specified for the mpi_msd application with the -msd::fitness_file <filename> flag" );
		}

		std::string entity_resfile( option[ msd::entity_resfile ] );
		std::string daf_filename( option[ msd::fitness_file ] );

		if ( mpi_rank() == 0 ) {
			devel::mmt_msd::MMTDriver mmt_driver;
			mmt_driver.set_nworkers( utility::mpi_nprocs() - 1 ); // num workers is one fewer than the number of mpi processes
			mmt_driver.set_ngenerations( option[ ms::generations ] );
			mmt_driver.set_pop_size( option[ ms::pop_size ] );
			mmt_driver.set_frac_by_recomb( option[ ms::fraction_by_recombination ] );
			mmt_driver.set_n_results_to_output( option[ ms::numresults ] );
			mmt_driver.set_daf_fname( daf_filename );
			mmt_driver.set_entity_resfile_fname( entity_resfile );

			// main optimization block
			mmt_driver.setup();
			mmt_driver.run();
			mmt_driver.write_optimal_solutions_to_disk();

		} else {
			devel::mmt_msd::MMTReceiver mmt_receiver;
			mmt_receiver.set_max_capacity( option[ msd::n_worker_threads_per_process ] );
			mmt_receiver.initial_handshake();
			mmt_receiver.main_optimization_loop();
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

#ifdef USEMPI
	MPI_Finalize();
#endif

	return 0;
}


