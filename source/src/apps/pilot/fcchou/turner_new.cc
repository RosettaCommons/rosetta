// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

// libRosetta headers
#include <core/init/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

#include <protocols/viewer/viewers.hh>

#include <protocols/moves/SimulatedTempering.hh>
#include <protocols/stepwise/sampling/rna/helix/RNA_HelixAssembler.hh>
#include <protocols/rotamer_sampler/rna/RNA_McSuite.hh>
#include <protocols/rotamer_sampler/rna/RNA_McMultiSuite.hh>

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
using namespace basic::options::OptionKeys;
using namespace basic::options;

OPT_KEY( String, seq1 )
OPT_KEY( String, seq2 )
OPT_KEY( Integer, n_cycle )
OPT_KEY( RealVector, temps )
OPT_KEY( RealVector, st_weights )
OPT_KEY( String, out_prefix )
OPT_KEY( Boolean, save_terms )
//////////////////////////////////////////////////////////////////////////////
utility::vector1<scoring::ScoreType> const & get_scoretypes() {
	using namespace scoring;
	static utility::vector1<ScoreType> scoretypes;
	if ( !scoretypes.empty() ) return scoretypes;
	if ( option[save_terms]() ) {
		scoretypes.push_back( fa_atr );
		scoretypes.push_back( fa_rep );
		scoretypes.push_back( hbond_sc );
		scoretypes.push_back( rna_torsion );
		scoretypes.push_back( fa_stack );
		scoretypes.push_back( geom_sol_fast );
		scoretypes.push_back( lk_nonpolar );
		scoretypes.push_back( stack_elec );
	}
	return scoretypes;
}
//////////////////////////////////////////////////////////////////////////////
Size data_dim() {
	using namespace scoring;
	utility::vector1<ScoreType> const & score_types( get_scoretypes() );
	return score_types.size() + 2;
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
						pose, score_types[i] ) );
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
void set_gaussian_stdev(
	rotamer_sampler::rna::RNA_McMultiSuite & sampler,
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
}
//////////////////////////////////////////////////////////////////////////////
PoseOP pose_setup(
	std::string const & seq1,
	std::string const & seq2,
	Size const len1
) {
	protocols::stepwise::sampling::rna::helix::RNA_HelixAssembler assembler;
	assembler.use_phenix_geo( true );
	PoseOP pose( assembler.build_init_pose( seq1, seq2 ) );
	add_variant_type_to_pose_residue( *pose, "VIRTUAL_PHOSPHATE", 1 );
	if ( seq1 != "" && seq2 != "" )	add_variant_type_to_pose_residue(
			*pose, "VIRTUAL_PHOSPHATE", len1 + 1 );
	return pose;
}
//////////////////////////////////////////////////////////////////////////////
void
MC_run () {
	using namespace protocols::rotamer_sampler::rna;
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
		scorefxn = getScoreFunction();
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
	RNA_McMultiSuite sampler;
	for ( Size i = 1; i <= total_len; ++i ) {
		bool const sample_near_a_form( bp_rsd.has_value( i ) );
		if  ( i == 1 || ( i > len1 && i != total_len ) ) {
			RNA_McSuiteOP suite_sampler( new RNA_McSuite( i ) );
			suite_sampler->set_sample_bb( i != 1 );
			suite_sampler->set_sample_lower_nucleoside( true );
			suite_sampler->set_sample_upper_nucleoside( false );
			suite_sampler->set_sample_near_a_form( sample_near_a_form );
			sampler.add_rotamer( suite_sampler );
		} else {
			RNA_McSuiteOP suite_sampler( new RNA_McSuite( i - 1 ) );
			suite_sampler->set_sample_bb( len1 == total_len || i != total_len );
			suite_sampler->set_sample_lower_nucleoside( false );
			suite_sampler->set_sample_upper_nucleoside( true );
			suite_sampler->set_sample_near_a_form( sample_near_a_form );
			sampler.add_rotamer( suite_sampler );
		}
	}
	sampler.init();
	sampler.apply( pose );

	// Simulated Tempering setup
	SimulatedTempering tempering(	pose,	scorefxn,	temps_,	weights_ );
	tempering.set_rep_cutoff( 100 );
	set_gaussian_stdev( sampler, tempering, bp_rsd, dangling_rsd );

	// Setup for data saving for output
	Size curr_counts( 1 );
	utility::vector1<float> scores;
	update_scores( scores, pose, scorefxn );
	utility::vector1<float> const null_arr_;
	utility::vector1<utility::vector1<float> > data(
			temps_.size(), null_arr_ );

	// Useful coounters and variables during the loop
	Size n_accept_total( 0 ), n_t_jumps_accept( 0 );
	Size const t_jump_interval( 10 );
	Size const n_t_jumps( n_cycle_ / t_jump_interval );
	Size temp_id( tempering.temp_id() );

	std::cout << "Start the main sampling loop." << std::endl;
	//pose.dump_pdb( "init.pdb" );

	// Main sampling cycle
	for ( Size n = 1; n <= n_cycle_; ++n ) {
		++sampler;
		sampler.apply( pose );
		if ( tempering.boltzmann( pose ) ) {
			sampler.update();
			++n_accept_total;
			fill_data( data[temp_id], curr_counts, scores );
			curr_counts = 1;
			update_scores( scores, pose, scorefxn );
		} else {
			++curr_counts;
		}

		if ( n % t_jump_interval == 0 && tempering.t_jump() ) {
			++n_t_jumps_accept;
			fill_data( data[temp_id], curr_counts, scores );
			curr_counts = 1;
			set_gaussian_stdev( sampler, tempering, bp_rsd, dangling_rsd );
			temp_id = tempering.temp_id();
		}
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
		std::ostringstream oss;
		oss << option[out_prefix]() << '_' << std::fixed << std::setprecision(2)
				<< temps_[i] << ".bin.gz";
		utility::io::ozstream out (
				oss.str().c_str(), std::ios::out | std::ios::binary );
		Size const data_dim1( data_dim() );
		Size const data_dim2( data[i].size() / data_dim1 );
		runtime_assert( data_dim1 * data_dim2 == data[i].size() );
		out.write( (const char*) &data_dim2, sizeof(Size) );
		out.write( (const char*) &data_dim1, sizeof(Size) );
		out.write( (const char*) &data[i][1], sizeof(float) * data[i].size() );
		out.close();
	}
}
//////////////////////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	MC_run();
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
	NEW_OPT( n_cycle, "cycle number for Random sampling", 0 );
	NEW_OPT( temps, "Simulated tempering temperatures", null_real_vector );
	NEW_OPT( st_weights, "Simulated tempering weights", null_real_vector );
	NEW_OPT( out_prefix, "prefix for the out file", "turner" );
	NEW_OPT( save_terms, "Save individual score terms", true );

	try {
		core::init::init ( argc, argv );
		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }
  return 0;
}

