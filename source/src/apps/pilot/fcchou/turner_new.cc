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

// libRosetta headers
#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

#include <protocols/viewer/viewers.hh>

#include <protocols/moves/SimulatedTempering.hh>
#include <protocols/stepwise/modeler/rna/helix/RNA_HelixAssembler.hh>
#include <protocols/stepwise/sampler/rna/RNA_MC_Suite.hh>
#include <protocols/stepwise/sampler/rna/RNA_MC_MultiSuite.hh>

#include <utility/io/ozstream.hh>

// C++ headers
#include <string>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Exception handling
#include <utility/excn/Exceptions.hh>

using namespace core;
using namespace core::pose;
using namespace protocols;
using namespace protocols::stepwise;
using namespace basic::options::OptionKeys;
using namespace basic::options;

OPT_KEY( String, seq1 )
OPT_KEY( String, seq2 )
OPT_KEY( Integer, n_cycle )
OPT_KEY( Real, a_form_range )
OPT_KEY( RealVector, temps )
OPT_KEY( RealVector, st_weights )
OPT_KEY( String, out_prefix )
OPT_KEY( Boolean, save_score_terms )
OPT_KEY( Boolean, dump_pdb )
//////////////////////////////////////////////////////////////////////////////
// Histogram class for accumulating samples
class Histogram {
public:
	Histogram( Real const min, Real const max, Real const spacing ):
		min_( min ),
		max_( max ),
		spacing_( spacing )
	{
		runtime_assert( max > min );
		n_elem_ = static_cast<Size>( ( max - min ) / spacing ) + 1;
		for ( Size i = 1; i <= n_elem_; ++i ) hist_.push_back( 0 );
	}

	void add( float const value, Size const n_items ) {
		Size bin_index;
		if ( value <= min_ ) {
			bin_index = 1;
		} else if ( value >= max_ ) {
			bin_index = n_elem_;
		} else {
			bin_index = static_cast<Size>( ( value - min_ ) / spacing_ ) + 1;
		}
		hist_[bin_index] += n_items;
	}

	void clear() {
		for ( Size i = 0; i <= n_elem_; ++i ) hist_[i] = 0;
	}

	utility::vector1<Real> get_scores() const {
		utility::vector1<Real> scores;
		for ( Size i = 1; i <= n_elem_; ++i )
			scores.push_back( min_ + spacing_ * ( i - 0.5 ) );
		return scores;
	}

	utility::vector1<Size> get_hist() const { return hist_;	}

private:
	Real const min_, max_, spacing_;
	Size n_elem_;
	utility::vector1<Size> hist_;
};
//////////////////////////////////////////////////////////////////////////////
// score types to be recorded
utility::vector1<scoring::ScoreType> const & get_scoretypes() {
	using namespace scoring;
	static utility::vector1<ScoreType> scoretypes;
	if ( !scoretypes.empty() ) return scoretypes;

	// List of score types to be cached
	scoretypes.push_back( fa_atr );
	scoretypes.push_back( fa_rep );
	scoretypes.push_back( fa_intra_rep );
	scoretypes.push_back( fa_stack );
	scoretypes.push_back( rna_torsion );
	scoretypes.push_back( hbond_sc );
	scoretypes.push_back( lk_nonpolar );
	scoretypes.push_back( geom_sol_fast );
	scoretypes.push_back( stack_elec );
	scoretypes.push_back( fa_elec_rna_phos_phos );
	return scoretypes;
}
//////////////////////////////////////////////////////////////////////////////
void update_scores(
	utility::vector1<float> & scores,
	Pose & pose,
	scoring::ScoreFunctionOP const scorefxn
) {
	using namespace scoring;
	scores.clear();
	scores.push_back( ( *scorefxn )( pose ) );
	utility::vector1<ScoreType> const & score_types( get_scoretypes() );
	for ( Size i = 1; i<= score_types.size(); ++i )
			scores.push_back( scorefxn->score_by_scoretype(
						pose, score_types[i], false /*weighted*/ ) );
}
//////////////////////////////////////////////////////////////////////////////
void fill_data(
	utility::vector1<float> & data,
	Size const count,
	utility::vector1<float> & scores
) {
	using namespace scoring;
	data.push_back( count );
	data.insert( data.end(), scores.begin(), scores.end() );
}
//////////////////////////////////////////////////////////////////////////////
// Simple heuristic for gaussian stdev
Real gaussian_stdev( Real const n_rsd, Real const temp, bool const is_bp ) {
	// Negative temp is infinite
	if ( temp < 0 ) return -1;
	if ( is_bp ) return 5 * temp / n_rsd;
	return 6 * pow( temp / n_rsd, 0.75 );
}
//////////////////////////////////////////////////////////////////////////////
template<typename T>
void vector2disk_in1d(
	std::string const & out_filename,
	utility::vector1<T> const & out_vector
) {
		utility::io::ozstream out (
				out_filename.c_str(), std::ios::out | std::ios::binary );
		out.write( (const char*) &out_vector[1], sizeof(T) * out_vector.size() );
		out.close();
}
//////////////////////////////////////////////////////////////////////////////
template<typename T>
void vector2disk_in2d(
	std::string const & out_filename,
	Size const dim1,
	Size const dim2,
	utility::vector1<T> const & out_vector
) {
		utility::io::ozstream out (
				out_filename.c_str(), std::ios::out | std::ios::binary );
		runtime_assert( dim1 * dim2 == out_vector.size() );
		out.write( (const char*) &dim1, sizeof(Size) );
		out.write( (const char*) &dim2, sizeof(Size) );
		out.write( (const char*) &out_vector[1], sizeof(T) * out_vector.size() );
		out.close();
}
//////////////////////////////////////////////////////////////////////////////

void set_gaussian_stdev(
	sampler::rna::RNA_MC_MultiSuite & sampler,
	moves::SimulatedTempering const & tempering,
	utility::vector1<Size> const & bp_rsd,
	utility::vector1<Size> const & dangling_rsd
) {
	Size const n_rsd( bp_rsd.size() + dangling_rsd.size() );
	Real const temp( tempering.temperature() );
	Real const bp_stdev( gaussian_stdev( n_rsd, temp, true ) );
	Real const dangling_stdev( gaussian_stdev( n_rsd, temp, false ) );
	for ( Size i = 1; i <= bp_rsd.size(); ++ i )
			sampler.set_gaussian_stdev( bp_stdev, bp_rsd[i] );
	for ( Size i = 1; i <= dangling_rsd.size(); ++ i )
			sampler.set_gaussian_stdev( dangling_stdev, dangling_rsd[i] );
	std::cout << "std:" << bp_stdev << std::endl;
}
//////////////////////////////////////////////////////////////////////////////
PoseOP pose_setup(
	std::string const & seq1,
	std::string const & seq2,
	Size const len1
) {
	protocols::stepwise::modeler::rna::helix::RNA_HelixAssembler assembler;
	assembler.use_phenix_geo( true );
	PoseOP pose( assembler.build_init_pose( seq1, seq2 ) );
	add_variant_type_to_pose_residue( *pose, chemical::VIRTUAL_PHOSPHATE, 1 );
	if ( seq1 != "" && seq2 != "" )	add_variant_type_to_pose_residue(
			*pose, chemical::VIRTUAL_PHOSPHATE, len1 + 1 );
	return pose;
}
////////////////////////////////////////////////////////////////////////////////
Size data_dim() {
	using namespace scoring;
	utility::vector1<ScoreType> const & score_types( get_scoretypes() );
	return score_types.size() + 2;
}
//////////////////////////////////////////////////////////////////////////////
void
MC_run () {
	using namespace protocols::stepwise::sampler::rna;
	using namespace protocols::moves;
	using namespace scoring;

	clock_t const time_start( clock() );

	utility::vector1<Real> const & temps_( option[ temps ]() );
	runtime_assert( temps_.size() != 0 );

	utility::vector1<Real> weights_;
	utility::vector1<Real> const & orig_weights( option[ st_weights ]() );
	if ( temps_.size() != orig_weights.size() )
			weights_.push_back( 0 );
	weights_.insert( weights_.end(), orig_weights.begin(),
			orig_weights.end() );
	runtime_assert( temps_.size() == weights_.size() );

	Size const n_cycle_( option[n_cycle]() );
	std::string const & seq1_( option[seq1]() );
	std::string const & seq2_( option[seq2]() );
	Size const len1( get_sequence_len( seq1_ ) );
	Size const len2( get_sequence_len( seq2_ ) );

	// Score function setup
	ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user() ) {
		scorefxn = get_score_function();
	} else {
		scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	}

	// Pose setup
	Pose pose( *pose_setup( seq1_, seq2_, len1 ) );

	// Figure out bp and dangling residues
	utility::vector1<Size> bp_rsd, dangling_rsd;
	Size const n_bp( std::min( len1, len2 ) );
	Size const total_len( len1 + len2 );
	for ( Size i = 1; i <= total_len; ++i ) {
		if ( i > n_bp && i <= total_len - n_bp ) {
			dangling_rsd.push_back( i );
		} else {
			bp_rsd.push_back( i );
		}
	}

	// Sampler setup
	RNA_MC_MultiSuite sampler;
	for ( Size i = 1; i <= total_len; ++i ) {
		bool const sample_near_a_form( bp_rsd.has_value( i ) );
		if ( i == 1 || ( i > len1 && i != total_len ) ) {
			RNA_MC_SuiteOP suite_sampler( new RNA_MC_Suite( i ) );
			suite_sampler->set_sample_bb( i != 1 );
			suite_sampler->set_sample_lower_nucleoside( true );
			suite_sampler->set_sample_upper_nucleoside( false );
			suite_sampler->set_sample_near_a_form( sample_near_a_form );
			suite_sampler->set_a_form_range( option[a_form_range]() );
			sampler.add_external_loop_rotamer( suite_sampler );
		} else {
			RNA_MC_SuiteOP suite_sampler( new RNA_MC_Suite( i - 1 ) );
			suite_sampler->set_sample_bb( len1 == total_len || i != total_len );
			suite_sampler->set_sample_lower_nucleoside( false );
			suite_sampler->set_sample_upper_nucleoside( true );
			suite_sampler->set_sample_near_a_form( sample_near_a_form );
			suite_sampler->set_a_form_range( option[a_form_range]() );
			sampler.add_external_loop_rotamer( suite_sampler );
		}
	}
	sampler.init();
	sampler.apply( pose );

	// Simulated Tempering setup
	SimulatedTempering tempering(	pose,	scorefxn,	temps_,	weights_ );
	// tempering.set_rep_cutoff( 100 );
	set_gaussian_stdev( sampler, tempering, bp_rsd, dangling_rsd );

	// Setup for data saving for output
	Size curr_counts( 1 );
	utility::vector1<float> scores;
	update_scores( scores, pose, scorefxn );
	utility::vector1<float> const null_arr_;
	utility::vector1<utility::vector1<float> > data(
			temps_.size(), null_arr_ );

	Real const min( -100.05 ), max( 800.05 ), spacing( 0.1 );
	Histogram null_hist( min, max, spacing);
	utility::vector1<Histogram> hist_list( temps_.size(), null_hist );

	// Useful coounters and variables during the loop
	Size n_accept_total( 0 ), n_t_jumps_accept( 0 );
	Size const t_jump_interval( 10 );
	Size const n_t_jumps( n_cycle_ / t_jump_interval );
	Size temp_id( tempering.temp_id() );
	bool const is_save_scores( option[save_score_terms]() );

	// Min-score pose
	Pose min_pose = pose;
	Real min_score( 99999 );

	std::cout << "Start the main sampling loop." << std::endl;
	if ( option[dump_pdb]() ) pose.dump_pdb( "init.pdb" );

	// Main sampling cycle
	for ( Size n = 1; n <= n_cycle_; ++n ) {
		++sampler;
		sampler.apply( pose );
		if ( tempering.boltzmann( pose ) || n == n_cycle_ ) {
			if ( is_save_scores ) fill_data( data[temp_id], curr_counts, scores );
			hist_list[temp_id].add( scores[1], curr_counts );
			update_scores( scores, pose, scorefxn );
			if ( n == n_cycle_ ) break;
			sampler.update();
			++n_accept_total;
			curr_counts = 1;
			if ( option[dump_pdb]() && scores[1] < min_score ) {
				min_score = scores[1];
				min_pose = pose;
			}
		} else {
			++curr_counts;
		}

		if ( n % t_jump_interval == 0 && tempering.t_jump() ) {
			++n_t_jumps_accept;
			if ( is_save_scores )	fill_data( data[temp_id], curr_counts, scores );
			hist_list[temp_id].add( scores[1], curr_counts );
			curr_counts = 1;
			set_gaussian_stdev( sampler, tempering, bp_rsd, dangling_rsd );
			temp_id = tempering.temp_id();
		}
	}
	if ( option[dump_pdb]() ) {
		pose.dump_pdb( "end.pdb" );
		min_pose.dump_pdb( "min.pdb" );
		scorefxn->show( min_pose );
	}

	// Output simple statistics and the data
	//pose.dump_pdb("final.pdb");

	std::cout << "n_cycles: " << n_cycle_ << std::endl;
	std::cout << "Accept rate: " << double( n_accept_total ) / n_cycle_
			<< std::endl;
	std::cout << "T_jump accept rate: " << double( n_t_jumps_accept ) / n_t_jumps
			<< std::endl;
	Real const time_in_test = static_cast<Real>( clock() - time_start )
			/ CLOCKS_PER_SEC;
	std::cout << "Time in sampler: " <<  time_in_test << std::endl;

	for (Size i = 1; i <= temps_.size(); ++i) {
		if ( is_save_scores ) {
			std::ostringstream oss;
			oss << option[out_prefix]() << '_' << std::fixed << std::setprecision(2)
					<< temps_[i] << ".bin.gz";
			Size const data_dim2( data_dim() );
			Size const data_dim1( data[i].size() / data_dim2 );
			vector2disk_in2d( oss.str(), data_dim1, data_dim2, data[i] );
		}
		std::ostringstream oss;
		oss << option[out_prefix]() << '_' << std::fixed << std::setprecision(2)
				<< temps_[i] << ".hist.gz";
		utility::vector1<Size> const & hist( hist_list[i].get_hist() );
		utility::vector1<Real> const & scores( hist_list[i].get_scores() );
		vector2disk_in1d( oss.str(), hist );

		std::ostringstream oss1;
		oss1 << option[out_prefix]() << "_hist_scores.gz";
		vector2disk_in1d( oss1.str(), scores );
	}
}
//////////////////////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	MC_run();
        protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}
//////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	using namespace core;
	utility::vector1< Real > null_real_vector;
	NEW_OPT( seq1, "sequence 1 to model, 3' to 5' ", "" );
	NEW_OPT( seq2, "sequence 2 to model, 3' to 5' ", "" );
	NEW_OPT( n_cycle, "cycle number for random sampling", 0 );
	NEW_OPT( a_form_range, "Range of sampling near A-form for duplexes.", 60.0);
	NEW_OPT( temps, "Simulated tempering temperatures", null_real_vector );
	NEW_OPT( st_weights, "Simulated tempering weights", null_real_vector );
	NEW_OPT( out_prefix, "prefix for the out file", "turner" );
	NEW_OPT( save_score_terms,
			"Save scores and individual score terms"
			" of all sampled conformers", false );
	NEW_OPT( dump_pdb, "Dump pdb files", false );

	try {
		devel::init ( argc, argv );
		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }
  return 0;
}

