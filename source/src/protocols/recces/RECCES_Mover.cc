// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/RECCES_Mover.cc
/// @brief Fast Monte Carlo, using MC_Sampler & simulated tempering.
/// @detailed Mover used by: recces_turner (Fang-Chieh Chou), thermal_sampler (K. Kappel + others)
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/recces/RECCES_Mover.hh>
#include <protocols/recces/options/RECCES_Options.hh>
#include <protocols/recces/params/RECCES_Parameters.hh>
#include <protocols/recces/sampler/MC_Comb.hh>
#include <protocols/recces/sampler/MC_Loop.hh>
#include <protocols/recces/sampler/rna/MC_RNA_OneJump.hh>
#include <protocols/recces/sampler/util.hh>
#include <protocols/recces/stdev_util.hh>
#include <protocols/recces/util.hh>
#include <core/id/TorsionID.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/moves/SimulatedTempering.hh>
#include <basic/Tracer.hh>
#include <utility>

static basic::Tracer TR( "protocols.recces.RECCES_Mover" );

///////////////////////////////////////////////////////////////////////////////////////////////////
/// TODO: The post-processing of the energies from this Mover are carried out with scripts in
///        tools/recces/ but its very easy to get tripped up with settings. A big help would be
///        output to disk a JSON file (e.g., "run_info.json") with:
///
///         1. The sequence simulated (e.g., ccc_gg), which is otherwise guessed by the python script.
///         2. score_types & weights (in output order)
///         3. Jumps & torsions that are moved (TorsionID's would be acceptable) and the angle_ranges (or rmsd_cutoff for jump)
///         4. For runs that move base pair jump(s), output of xyz.txt
///         5. Histogram information: histogram_min, histogram_max, histogram_bin_size
///

using namespace core;
using namespace utility;
using namespace utility::tools;
using namespace protocols::recces;
using namespace protocols::toolbox;

namespace protocols {
namespace recces {

//Constructor
RECCES_Mover::RECCES_Mover( options::RECCES_OptionsCOP  options ):
	options_(std::move( options ))
{
	initialize();
}

//Destructor
RECCES_Mover::~RECCES_Mover() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::apply( core::pose::Pose & pose )
{
	runtime_assert( params_ );
	runtime_assert( options_ );
	sampler_ = sampler::initialize_sampler( pose, *options_, *params_ );

	run_sampler( pose );

	save_data_to_disk();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::run_sampler( pose::Pose & pose )
{
	using namespace core::pose;
	using namespace core::scoring;
	using namespace recces::sampler;

	clock_t const time_start( clock() );

	runtime_assert( options_ );

	// Are we doing MC_loop?
	MC_LoopOP loop_sampler = ( sampler_->type() == toolbox::MC_LOOP ) ? MC_LoopOP( utility::pointer::dynamic_pointer_cast< MC_Loop >( sampler_) ) : nullptr;

	if ( options_->show_more_pose_scores() ) scorefxn_->show( pose ); // matching legacy output
	if ( loop_sampler == nullptr ) sampler_->apply( pose ); // matching legacy output

	if ( options_->out_torsions() ) prepare_output_torsion_ids();

	// Simulated Tempering setup
	if ( scorefxn_ == nullptr ) utility_exit_with_message( "Must set scorefunction." );
	if ( options_->setup_base_pair_constraints() ) runtime_assert( scorefxn_->has_nonzero_weight( base_pair_constraint ) );
	moves::SimulatedTempering tempering( pose, scorefxn_, options_->temperatures(), weights_ );
	set_sampler_gaussian_stdev( tempering.temperature(), pose );

	// Setup for data saving for output
	Size curr_counts( 1 );
	utility::vector1< float > scores;
	utility::vector1< ScoreType > const score_types = options_->blank_score_terms() ?  vector1< ScoreType>() : get_scoretypes();
	update_scores( scores, pose, scorefxn_, score_types );

	// Useful counters and variables during the loop
	Size n_accept_total( 0 ), n_t_jumps_accept( 0 );
	Size const t_jump_interval( 10 );
	Size const n_t_jumps( options_->n_cycle() / t_jump_interval );
	Size temp_id( tempering.temp_id() );

	// Min-score pose. And stored_pose, in use in loop mode to restore pose. Not quick.
	Pose min_pose( pose ), stored_pose( pose );
	Real min_score( 99999 );

	TR << "Start the main sampling loop." << std::endl;
	if ( options_->dump_pdb() ) pose.dump_pdb( output_pdb_name( "init" ) );

	Size curr_dump = 1;

	// Main sampling cycle
	Size pct = 0;
	std::map< std::string, Size > num_accepts, num_tries;
	for ( Size n = 1; n <= options_->n_cycle(); ++n ) {

		if ( options_->n_cycle() > 100 && n % ( options_->n_cycle()/100 ) == 0 ) TR << ++pct << "% complete." << std::endl;

		//////////////////////////////
		// Do a move
		//////////////////////////////
		++( *sampler_ ); // prepare DOFS for next conformation
		sampler_->apply( pose );

		bool const did_move( sampler_->found_move() || options_->accept_no_op_moves() );
		if ( loop_sampler != nullptr ) num_tries[ loop_sampler->rotamer()->get_name() ]++;

		//////////////////////////////
		// Check metropolis criterion
		//////////////////////////////
		if ( ( tempering.check_boltzmann( pose ) && did_move ) || n == options_->n_cycle() ) {
			// accept new pose.
			save_history( curr_counts, scores, temp_id ); // save hist & data so far, weighted by curr_counts
			update_scores( scores, pose, scorefxn_, score_types );
			if ( !options_->skip_last_accept() ) increment_accepts( n_accept_total,  num_accepts, loop_sampler );

			if ( n == options_->n_cycle() ) break;

			// critical: tell sampler that new DOFs have been accepted.
			sampler_->update();
			if ( options_->skip_last_accept() ) increment_accepts( n_accept_total,  num_accepts, loop_sampler );
			if ( options_->out_torsions() ) output_torsions( pose, curr_counts );
			dump_stuff( n, pose, scores, min_pose, min_score, curr_dump );
			curr_counts = 1;

			if ( loop_sampler != nullptr ) stored_pose = pose; // eventually deprecate this
		} else {
			// reject new pose
			++curr_counts; // pose stayed the same.
			if ( loop_sampler != nullptr && sampler_->found_move() ) pose = stored_pose; // let's eventually replace this with sampler->restore( pose )
		}

		//////////////////////////
		// Check temperature jump.
		//////////////////////////
		if ( n % t_jump_interval == 0 && tempering.t_jump() ) {
			++n_t_jumps_accept;
			save_history( curr_counts, scores, temp_id ); // before jumping temperature, save history at this pose so far.
			curr_counts = 1;
			set_sampler_gaussian_stdev( tempering.temperature(), pose);
			temp_id = tempering.temp_id();
		}

		more_dump_stuff( pose, n, tempering.temperature() ); // legacy stuff. do we need so much dumping?
	}

	final_dump_stuff( pose, min_pose );

	TR << "n_cycles: " << options_->n_cycle() << std::endl;
	TR << "Accept rate: " << double( n_accept_total ) / options_->n_cycle() << std::endl;
	if ( num_accepts.size() > 0 ) {
		for ( auto const & it : num_accepts ) {
			std::string const & sampler_name( it.first );
			TR << " " <<  sampler_name << " accept rate: " << double( num_accepts[ sampler_name ] )/double( num_tries[ sampler_name ] ) << std::endl;
		}
	}
	TR << "T_jump accept rate: " << double( n_t_jumps_accept ) / n_t_jumps << std::endl;

	Real const time_in_test = static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC;
	TR << "Time in sampler: " <<  time_in_test << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::initialize() {
	runtime_assert( options_ );

	Size const num_temperatures = options_->temperatures().size();
	runtime_assert( num_temperatures != 0 );

	weights_.clear();
	utility::vector1< Real > const & orig_weights( options_->st_weights() );
	if ( ( num_temperatures != orig_weights.size() ) && ( orig_weights.empty() || orig_weights[1] != 0 ) ) weights_.push_back( 0 );
	weights_.append( orig_weights );
	runtime_assert( num_temperatures == weights_.size() );

	data_ = utility::vector1<utility::vector1< float > >( num_temperatures, utility::vector1< float >() );

	// need to fix these magic numbers:
	Real const min( options_->histogram_min() ), max( options_->histogram_max() ), spacing( options_->histogram_spacing() );
	hist_list_ = vector1< Histogram >( num_temperatures, Histogram( min, max, spacing) );

}

//////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::set_sampler_gaussian_stdev( Real const & temperature, pose::Pose const & pose )
{
	runtime_assert( params_ );
	runtime_assert( options_ );
	recces::set_sampler_gaussian_stdev( sampler_, temperature, pose, *options_, *params_ );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// @details curr_counts is a running count of number of times the pose has been the same.
void
RECCES_Mover::save_history(
	Size const & curr_counts,
	vector1< float > const & scores,
	Size const & temp_id )
{
	runtime_assert( options_ );

	if ( options_->save_scores() ) fill_data( data_[ temp_id ], curr_counts, scores );
	hist_list_[ temp_id ].add( scores[ 1 ], curr_counts );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::increment_accepts( Size & n_accept_total,
	std::map< std::string, Size > & num_accepts,
	sampler::MC_LoopCOP loop_sampler ) const
{
	++n_accept_total;
	if ( loop_sampler != nullptr ) num_accepts[ loop_sampler->rotamer()->get_name() ]++;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
std::string
RECCES_Mover::output_pdb_name( std::string const & tag ) const {
	runtime_assert( options_ );

	if ( options_->prefix_each_output_pdb() ) {
		return ( options_->out_prefix() + "_" + tag + ".pdb" );
	}
	return ( tag + ".pdb" );
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::dump_stuff(
	Size const & n,
	pose::Pose const & pose,
	utility::vector1< Real > const & scores,
	pose::Pose & min_pose,
	Real & min_score,
	Size & curr_dump ) const
{
	runtime_assert( options_ );

	if ( options_->dump_pdb() && scores[ 1 ] < min_score ) {
		min_score = scores[ 1 ];
		min_pose = pose;
	}
	if ( options_->n_dump() != 0 &&
			n * (options_->n_dump() + 1) / double(options_->n_cycle()) >= curr_dump ) {
		std::ostringstream oss;
		oss << "intermediate" << '_' << curr_dump << ".pdb";
		pose.dump_pdb(oss.str());
		++curr_dump;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::more_dump_stuff( pose::Pose const & pose,
	Size const & n,
	Real const & current_temperature ) const
{
	runtime_assert( options_ );

	using namespace core::io::silent;
	Size const & dump_interval( options_->dump_freq() );
	if ( (n % dump_interval) == 0 && options_->dump_pdb() ) {
		std::stringstream str;
		str << options_->out_prefix() << "_" << n << ".pdb";
		pose.dump_pdb( str.str() );
	}
	if ( (n % dump_interval) == 0 && options_->dump_silent() ) {
		std::stringstream tag;
		tag << options_->out_prefix() << "_" << n;
		SilentFileOptions opts;
		BinarySilentStruct silent( opts, pose, tag.str() );
		silent.add_energy( "temp", current_temperature );
		SilentFileData silent_file_data( opts );
		silent_file_data.write_silent_struct( silent, options_->silent_file(), false /*write score only*/ );
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::final_dump_stuff( pose::Pose & pose,
	pose::Pose & min_pose ) const
{
	runtime_assert( options_ );

	if ( options_->dump_pdb() ) {
		pose.dump_pdb( output_pdb_name( "end" ) );
		min_pose.dump_pdb( output_pdb_name( "min" ) );
		if ( options_->show_more_pose_scores() ) {
			std::cout << "Score for the end pose" << std::endl;
			scorefxn_->show( pose );
		}
		if ( options_->show_more_pose_scores() ) {
			std::cout << "Score for the min pose" << std::endl;
		}
		scorefxn_->show( min_pose );
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::save_data_to_disk() const
{
	runtime_assert( options_ );

	// this is ridiculous -- need to fix external code to just look at score type names written to a file.
	utility::vector1< core::scoring::ScoreType > const score_types = options_->blank_score_terms() ? vector1< core::scoring::ScoreType>() : get_scoretypes();

	for ( Size i = 1; i <= options_->temperatures().size(); ++i ) {
		if ( options_->save_scores() ) {
			std::ostringstream oss;
			oss << options_->out_prefix() << '_' << std::fixed << std::setprecision( 2 )
				<< options_->temperatures()[ i ] << ".bin.gz";
			Size const data_dim2( data_dim( score_types ) );
			Size const data_dim1( data_[ i ].size() / data_dim2 );
			vector2disk_in2d( oss.str(), data_dim1, data_dim2, data_[ i ] );
		}
		std::ostringstream oss;
		oss << options_->out_prefix() << '_' << std::fixed << std::setprecision( 2 )
			<< options_->temperatures()[ i ] << ".hist.gz";
		utility::vector1< Size > const & hist( hist_list_[ i ].get_hist() );
		utility::vector1< Real > const & scores( hist_list_[ i ].get_scores() );
		vector2disk_in1d( oss.str(), hist );

		std::ostringstream oss1;
		oss1 << options_->out_prefix() << "_hist_scores.gz";
		vector2disk_in1d( oss1.str(), scores );
	}

	if ( bb_tor_out_.is_open() ) bb_tor_out_.close();
	if ( chi_tor_out_.is_open() ) chi_tor_out_.close();

	// wow super verbose legacy -- these will be huge text files.
	if ( options_->output_simple_text_files() ) {
		// Some simple text data files
		std::stringstream str, str2, str3;
		str << options_->out_prefix() << "_data.txt";
		str2 << options_->out_prefix() << "_hist.txt";
		str3 << options_->out_prefix() << "_scores.txt";
		std::ofstream datafile ( (str.str()).c_str() );
		for ( Size i = 1; i <= data_[1].size(); i+=2 ) {
			datafile << data_[1][i] << "   " << data_[1][i+1] << "\n";
		}
		datafile.close();
		std::ofstream histfile ( (str2.str()).c_str() );
		utility::vector1<Size> const & hist( hist_list_[1].get_hist() );
		for ( Size i = 1; i<=hist.size(); i++ ) {
			histfile<< hist[i] << "\n";
		}
		histfile.close();
		utility::vector1<Real> const & escores( hist_list_[1].get_scores() );
		std::ofstream scorefile ( (str3.str()).c_str() );
		for ( Size i = 1; i<=escores.size(); i++ ) {
			scorefile<< escores[i] << "\n";
		}
		scorefile.close();
	}

}

/////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::prepare_output_torsion_ids()
{
	using namespace core::id;
	using namespace core::chemical::rna;

	runtime_assert( options_ );
	//Make a list of all the backbone torsion IDs that are sampled (This is only complete if all the residues are consecutive)
	// REPLACE this with sampler.find(), looking over all torsions!
	utility::vector1< Size > const & residues( options_->sample_residues() );
	bb_torsion_ids_.clear();
	bb_torsion_ids_.emplace_back( residues[1]-1, id::BB, EPSILON );
	bb_torsion_ids_.emplace_back( residues[1]-1, id::BB, ZETA );
	for ( Size const residue : residues ) {
		bb_torsion_ids_.emplace_back( residue, id::BB, ALPHA );
		bb_torsion_ids_.emplace_back( residue, id::BB, BETA );
		bb_torsion_ids_.emplace_back( residue, id::BB, GAMMA );
		bb_torsion_ids_.emplace_back( residue, id::BB, EPSILON );
		bb_torsion_ids_.emplace_back( residue, id::BB, ZETA );
	}
	bb_torsion_ids_.emplace_back( residues.back()+1, id::BB, ALPHA );
	bb_torsion_ids_.emplace_back( residues.back()+1, id::BB, BETA );
	bb_torsion_ids_.emplace_back( residues.back()+1, id::BB, GAMMA );

	chi_torsion_ids_.clear();
	for ( Size const residue : residues ) {
		chi_torsion_ids_.emplace_back( residue , id::CHI, 1 );
	}

	std::stringstream name, name2;
	name << options_->out_prefix() << "_bb_torsions.txt";
	name2 << options_->out_prefix() << "_chi_torsions.txt";
	bb_tor_out_.open( (name.str()).c_str() );
	chi_tor_out_.open( (name2.str()).c_str() );
}

/////////////////////////////////////////////////////////////////////////////
void
RECCES_Mover::output_torsions( core::pose::Pose const & pose, Size const & curr_counts ) const
{
	utility::vector1< Real > const bb_torsions = get_torsions(bb_torsion_ids_, pose);
	utility::vector1< Real > const chi_torsions = get_torsions(chi_torsion_ids_, pose);
	bb_tor_out_ << curr_counts << "  " << bb_torsions <<"\n";
	chi_tor_out_ << curr_counts << "  " << chi_torsions <<"\n";
}

} //recces
} //protocols
