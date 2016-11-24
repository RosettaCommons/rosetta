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
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/data/RNA_ChemicalMappingEnergy.hh>
#include <core/scoring/rna/data/RNA_DataInfo.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/database/open.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/util.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/align/util.hh>
#include <protocols/farna/setup/RNA_DeNovoParameters.hh>
#include <protocols/farna/util.hh>
#include <core/io/rna/RNA_DataReader.hh>
#include <core/pose/PDBInfo.hh>

#include <protocols/stepwise/sampler/rna/RNA_MC_Suite.hh>
#include <protocols/stepwise/sampler/rna/RNA_MC_MultiSuite.hh>
#include <protocols/moves/SimulatedTempering.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/stepwise/sampler/rna/RNA_MC_KIC_Sampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_KIC_Sampler.hh>
#include <core/id/TorsionID.hh>
#include <protocols/stepwise/sampler/MC_OneTorsion.hh>
#include <utility/io/ozstream.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <numeric/random/random.hh>

#include <utility/excn/Exceptions.hh>


// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

using namespace core::pose;
using namespace basic::options;


using namespace core;
using namespace protocols;
using namespace protocols::stepwise;
using namespace protocols::moves;
using namespace basic::options::OptionKeys;
using utility::vector1;

OPT_KEY( String, seq1 )
OPT_KEY( String, seq2 )
OPT_KEY( Integer, n_cycle )
OPT_KEY( RealVector, temps )
OPT_KEY( RealVector, st_weights )
OPT_KEY( String, out_prefix )
OPT_KEY( Boolean, save_terms )
OPT_KEY( Boolean, save_scores )
OPT_KEY( Boolean, dump_pdb )
OPT_KEY( Boolean, dump_silent )
OPT_KEY( IntegerVector, sample_residues )
OPT_KEY( IntegerVector, free_residues )
OPT_KEY( Real, angle_range_chi )
OPT_KEY( Real, angle_range_free_chi )
OPT_KEY( Real, angle_range_bb )
OPT_KEY( Real, angle_range_free_bb )
OPT_KEY( Real, chi_stdev )
OPT_KEY( Real, bb_stdev )
OPT_KEY( Real, standard_bb_stdev )
OPT_KEY( Boolean, out_torsions )
OPT_KEY( Integer, dump_freq )
//OPT_KEY( RealVector, sample_residues)

///// First a bunch of stuff copied from Fang's turner_new code ///////
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
	// Don't do this for now
	//if ( option[save_terms]() ) {
	//	scoretypes.push_back( fa_atr );
	//	scoretypes.push_back( fa_rep );
	//	scoretypes.push_back( fa_intra_rep );
	//	scoretypes.push_back( fa_stack );
	//	scoretypes.push_back( rna_torsion );
	//	scoretypes.push_back( hbond_sc );
	//	scoretypes.push_back( lk_nonpolar );
	//	scoretypes.push_back( geom_sol_fast );
	//	scoretypes.push_back( stack_elec );
	//	scoretypes.push_back( fa_elec_rna_phos_phos );
	//}
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
utility::vector1<core::Real> get_torsions(
	utility::vector1<core::id::TorsionID> & torsion_ids,
	const Pose & pose
) {
	utility::vector1<core::Real> curr_torsions;
	for (Size i = 1; i <= torsion_ids.size(); ++i) {
		curr_torsions.push_back( pose.torsion( torsion_ids[i] ) );
	}
	return curr_torsions;
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
		// What if the vector is empty because of a short trajectory? Let's say
		// we just skip.
		if ( out_vector.size() != 0 ) out.write( (const char*) &out_vector[1], sizeof(T) * out_vector.size() );
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
		// What if the vector is empty because of a short trajectory? Let's say
		// we just skip.
		if ( out_vector.size() != 0 ) out.write( (const char*) &out_vector[1], sizeof(T) * out_vector.size() );
		out.close();
}
//////////////////////////////////////////////////////////////////////////////
void set_gaussian_stdevs( 
	utility::vector1<protocols::stepwise::sampler::rna::RNA_MC_KIC_SamplerOP> & internal_bb_sampler,
	utility::vector1<protocols::stepwise::sampler::MC_OneTorsionOP> & chi_sampler,
	sampler::rna::RNA_MC_MultiSuite & standard_bb_sampler,
	moves::SimulatedTempering const & tempering,
	Size const & total_rsd,
	Size const & sampled_rsd,
	utility::vector1<bool> is_free
) {
	Real const temp( tempering.temperature() );
	Real internal_bb_stdev( 0.1 * pow( temp, 0.25 ) + 0.1);
	Real free_chi_stdev( 55 * pow( temp, 0.5 ) + 50 );
	Real chi_stdev( 5 * pow( temp , 0.5) + 15 );
	Real standard_bb_stdev( 8 * pow( temp, 0.5 ) / (2 * total_rsd + sampled_rsd) );
	if ( temp < 0 ) {
		internal_bb_stdev = 0.5 ;
		free_chi_stdev = -1 ;
		chi_stdev = -1 ;
		standard_bb_stdev = -1 ;
	} 
	for ( Size i = 1; i <= internal_bb_sampler.size(); ++i ) {
		internal_bb_sampler[i]->set_gaussian_stdev( internal_bb_stdev );
	}
	for ( Size i = 1; i <= chi_sampler.size(); ++i ) {
		if ( is_free[i] )
			chi_sampler[i]->set_gaussian_stdev( free_chi_stdev );
		else chi_sampler[i]->set_gaussian_stdev( chi_stdev );
	}
	standard_bb_sampler.set_gaussian_stdev( standard_bb_stdev );
}
/////////////////////////////////////////////////////////////////////////////////
Size data_dim() {
	using namespace scoring;
	utility::vector1<ScoreType> const & score_types( get_scoretypes() );
	return score_types.size() + 2;
}
///////////////////////////////////////////////////////////////////////////////
void
thermal_sampler()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
    using namespace core::chemical::rna;
	using namespace core::scoring::rna::data;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose::full_model_info;
	using namespace protocols::stepwise::modeler;

    
    using namespace protocols::stepwise::setup;
	using namespace protocols::stepwise::sampler::rna;
	using namespace protocols::moves;
	using namespace core::id;
	using namespace protocols::stepwise::sampler;

	clock_t const time_start( clock() );

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD /*RNA*/ );

	// input stream
	PoseInputStreamOP input;
	if ( option[ in::file::silent ].user() ) {
		if ( option[ in::file::tags ].user() ) {
			input = PoseInputStreamOP( new SilentFilePoseInputStream(
				option[ in::file::silent ](),
				option[ in::file::tags ]()
			) );
		} else {
			input = PoseInputStreamOP( new SilentFilePoseInputStream( option[ in::file::silent ]() ) );
		}
	} else {
		input = PoseInputStreamOP( new PDBPoseInputStream( option[ in::file::s ]() ) );
	}

	// score function setup
	core::scoring::ScoreFunctionOP scorefxn;
	if ( basic::options::option[ basic::options::OptionKeys::score::weights ].user() ) {
		scorefxn = core::scoring::get_score_function();
	} else {
		scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::RNA_HIRES_WTS );
	}

	// Silent file output setup
	std::string const silent_file = option[ out::file::silent  ]();
	SilentFileData silent_file_data;

	FullModelInfoOP my_model;
	utility::vector1< pose::PoseOP > other_poses;

	// if trying to compute stem RMSD
	RNA_ChemicalMappingEnergyOP rna_chemical_mapping_energy;
	pose::Pose pose,start_pose;

	input->fill_pose( pose, *rsd_set );

	if ( !option[ in::file::silent ].user() ) cleanup( pose );

	if ( !full_model_info_defined( pose ) || option[ in::file::fasta ].user() ){
			fill_full_model_info_from_command_line( pose, other_poses ); // only does something if -in:file:fasta specified.
	}


	// do it
	//(*scorefxn)( pose );
	scorefxn->show( pose );


	Size const total_len( pose.total_residue() );
	//bool const sample_near_a_form( false );


	//Setting up the internal move sampler

//////////THIS IS THE CODE THAT was working before, but causing malloc error upon function exit/////////////
//	pose::Pose ref_pose_;
//	pose::PoseOP ref_pose;
//	ref_pose_ = pose;
//	ref_pose = PoseOP (&ref_pose_);
/////////DONE//////////////
	///// This fixes the error -- but have not tested extensively to make sure it is really ok /////

	// Obtains base pairs and constraints from pose
	utility::vector1< std::pair< Size, Size > > pairings;
	protocols::farna::get_base_pairing_list( pose, pairings );
    protocols::farna::setup_base_pair_constraints( pose, pairings );

	pose::PoseOP ref_pose( new pose::Pose( pose ) );
	

	utility::vector1<int> residues( option[sample_residues]() );
	Size const total_sampled( residues.size() );
	utility::vector1<int> free_rsd( option[free_residues]() );
	utility::vector1<bool> is_free;
	std::cout << "Sample residues: " << option[sample_residues]() <<std::endl;

	//Assign "free" residues (get bigger gaussian stdev)
	for ( Size i = 1; i <= residues.size(); ++i ) {
		if ( std::find( free_rsd.begin(), free_rsd.end(), residues[i]) != free_rsd.end() ) {
			is_free.push_back( true );
		} else is_free.push_back( false );
	}
	
	//Set up the internal move samplers
	utility::vector1<RNA_MC_KIC_SamplerOP> sampler;
	for ( Size i = 1; i<= residues.size(); ++i ) {
		RNA_MC_KIC_SamplerOP suite_sampler( new RNA_MC_KIC_Sampler( ref_pose, residues[i]-1, residues[i]) );
		suite_sampler->init();
		if ( is_free[i] || ( i > 1 && is_free[i-1] ) ) {
			suite_sampler->set_angle_range_from_init_torsions( option[angle_range_free_bb]() );
		} else {
			suite_sampler->set_angle_range_from_init_torsions( option[angle_range_bb]() );
		}
		sampler.push_back( suite_sampler );
	}

	//Set up the chi samplers
	utility::vector1<MC_OneTorsionOP> chi_sampler;
	core::Real init_torsion; 
	utility::vector1<TorsionID> chi_torsion_ids;
	for ( Size i = 1; i<= residues.size(); ++i ) {
		TorsionID chi_ID (TorsionID( residues[i] , id::CHI, 1));
		chi_torsion_ids.push_back( chi_ID );
		init_torsion = ref_pose->torsion( chi_ID );
		MC_OneTorsionOP chi_torsion( new MC_OneTorsion( chi_ID, ref_pose->torsion( chi_ID)));
		chi_torsion->init();
        
		if ( is_free[i] ) {
            chi_torsion->set_angle_range( init_torsion - option[angle_range_free_chi]() , init_torsion + option[angle_range_free_chi]() );
        }
        else {
			chi_torsion->set_angle_range( init_torsion - option[angle_range_chi]() , init_torsion + option[angle_range_chi]() );
        }
        
		chi_sampler.push_back( chi_torsion );
	}

	//Set up the regular torsion samplers, necessary for free energy comparisons
	//Don't sample chi torsions here because there is a separate chi sampler
	RNA_MC_MultiSuite standard_sampler;
	for ( Size i = 1; i <= residues.size(); ++i ) {
		if ( i == 1 || ( i > 1 && residues[i] != residues[i-1]+1 ) ) { 
			//create samplers for [i]-1 and [i]
			RNA_MC_SuiteOP suite_sampler_1( new RNA_MC_Suite( residues[i] - 1 ) );
			suite_sampler_1->set_init_from_pose( pose );
			suite_sampler_1->set_sample_bb( true );
			suite_sampler_1->set_sample_lower_nucleoside( false );
			suite_sampler_1->set_sample_upper_nucleoside( false );
			suite_sampler_1->set_pucker_flip_rate( 0 );
			suite_sampler_1->set_sample_near_a_form( false );
			if ( is_free[i] ) {
				suite_sampler_1->set_angle_range_from_init( option[angle_range_free_bb]() );
			} else {
				suite_sampler_1->set_angle_range_from_init( option[angle_range_bb]() );
			}
			standard_sampler.add_external_loop_rotamer( suite_sampler_1 );
		} 
		RNA_MC_SuiteOP suite_sampler( new RNA_MC_Suite( residues[i] ) );
		suite_sampler->set_init_from_pose( pose );
		suite_sampler->set_sample_bb( true );
		suite_sampler->set_sample_lower_nucleoside( false );
		suite_sampler->set_sample_upper_nucleoside( false );
		suite_sampler->set_pucker_flip_rate( 0 );
		suite_sampler->set_sample_near_a_form( false );
        
		if ( is_free[i] ) {
			suite_sampler->set_angle_range_from_init( option[angle_range_free_bb]() );
		} else {
			suite_sampler->set_angle_range_from_init( option[angle_range_bb]() );
		}
		standard_sampler.add_external_loop_rotamer( suite_sampler );
	}
    
    
	standard_sampler.init();

	//Make a list of all the backbone torsion IDs that are sampled (This is only complete if all the residues are consecutive)
	utility::vector1<TorsionID> bb_torsion_ids;
	bb_torsion_ids.push_back( TorsionID( residues[1]-1, id::BB, EPSILON ) );
	bb_torsion_ids.push_back( TorsionID( residues[1]-1, id::BB, ZETA ) );
	for ( Size i = 1; i <= residues.size(); ++i ) {
		bb_torsion_ids.push_back( TorsionID( residues[i], id::BB, ALPHA) );	
		bb_torsion_ids.push_back( TorsionID( residues[i], id::BB, BETA) );	
		bb_torsion_ids.push_back( TorsionID( residues[i], id::BB, GAMMA) );	
		bb_torsion_ids.push_back( TorsionID( residues[i], id::BB, EPSILON) );	
		bb_torsion_ids.push_back( TorsionID( residues[i], id::BB, ZETA) );	
	}
	bb_torsion_ids.push_back( TorsionID( residues.back()+1, id::BB, ALPHA ) );
	bb_torsion_ids.push_back( TorsionID( residues.back()+1, id::BB, BETA ) );
	bb_torsion_ids.push_back( TorsionID( residues.back()+1, id::BB, GAMMA ) );

	utility::vector1<core::Real> bb_torsions;
	utility::vector1<core::Real> chi_torsions;
	
	
	bool const is_save_scores( true );

        utility::vector1<Real> const & temps_( option[ temps ]() );
        runtime_assert( temps_.size() != 0 );

        utility::vector1<Real> weights_;
        utility::vector1<Real> const & orig_weights( option[ st_weights ]() );
        if ( temps_.size() != orig_weights.size() )
                        weights_.push_back( 0 );
        weights_.insert( weights_.end(), orig_weights.begin(),
                        orig_weights.end() );
        runtime_assert( temps_.size() == weights_.size() );

	Size curr_counts( 1 );
	utility::vector1<float> scores;
	update_scores( scores, pose, scorefxn );
	utility::vector1<float> const null_arr_;
	utility::vector1<utility::vector1<float> > data(
	                temps_.size(), null_arr_ );

	Real const min( -100.05 ), max( 100.05 ), spacing( 0.1 );
	Histogram null_hist( min, max, spacing);
	utility::vector1<Histogram> hist_list( temps_.size(), null_hist );

	// Min-score pose
	Pose min_pose = pose;
	Real min_score( 99999 );

	SimulatedTempering tempering(   pose,   scorefxn,       temps_, weights_ );
	MonteCarlo mc(pose, *scorefxn, 1. );
	Vector center_vector = Vector( 0.0 );
	Size n_accept_total( 0 ), n_accept_backbone( 0 );
	Size n_accept_chi( 0 ), n_accept_standard( 0 ), n_t_jumps_accept( 0 );
	Size const t_jump_interval( 10 );
	int const dump_interval( option[dump_freq]() );
	bool did_move;
	Size const n_cycle_( option[n_cycle]() );
	int index;
	Size temp_id( 1 ); 

	std::stringstream name, name2;
	name << option[out_prefix]() << "_bb_torsions.txt";
	name2 << option[out_prefix]() << "_chi_torsions.txt";
	std::ofstream bb_tor_file ( (name.str()).c_str() );
	std::ofstream chi_tor_file ( (name2.str()).c_str() );


	std::stringstream initstr;
	initstr << option[out_prefix]() << "_init.pdb";
	pose.dump_pdb( initstr.str() );

	set_gaussian_stdevs( sampler, chi_sampler, standard_sampler, 
		tempering, total_len, total_sampled, is_free );

	Pose stored_pose_ = pose;
	//protocols::viewer::add_conformation_viewer ( pose.conformation(), "current", 700, 700, false, false , center_vector );
	// Main sampling cycle
    Size pct = 0;
	for ( Size n = 1; n <= n_cycle_; ++n ) {
        
        if ( n % ( n_cycle_ / 100 ) == 0 ) {
            ++pct;
            std::cout << pct << "% complete." << std::endl;
        }
        
		did_move = true;
		index = numeric::random::rg().random_range(1,sampler.size());
		if ( (n % 10) == 0 ) {
			// We need to update the torsions in case they were updated by the other samplers
			standard_sampler.set_angle( pose );
			++standard_sampler;
			standard_sampler.apply( pose );
		}
		else if ( (n % 2) == 0 ) {
			sampler[index]->next( pose ); //This function also updates the stored torsions
			sampler[index]->apply( pose );
			if ( !(sampler[index]->check_moved() ) ) {
				did_move = false;
			}
		} else {
			++(*chi_sampler[index]);
			chi_sampler[index]->apply( pose );
		}
		if ( ( tempering.boltzmann( pose ) && did_move ) || n == n_cycle_ ) {
			stored_pose_ = pose;
			if ( is_save_scores ) fill_data( data[temp_id], curr_counts, scores );
			hist_list[temp_id].add( scores[1], curr_counts );
			update_scores( scores, pose, scorefxn );
			if ( n == n_cycle_ ) break;
			if ( (n % 10) == 0 ) {
				standard_sampler.update();
				++n_accept_standard;
			}
			else if ( (n % 2) == 0 ) {
				sampler[index]->update( pose ); // don't think this is necessary
				++n_accept_backbone;
			} else {
				chi_sampler[index]->update();
				++n_accept_chi;
			}
			if ( option[out_torsions]() ) {
				bb_torsions = get_torsions(bb_torsion_ids, pose);
				chi_torsions = get_torsions(chi_torsion_ids, pose);
				bb_tor_file << curr_counts<< "  "<< bb_torsions <<"\n";
				chi_tor_file << curr_counts<< "  "<< chi_torsions <<"\n";
			}
			++n_accept_total;
			curr_counts = 1;
			if ( scores[1] < min_score ) {
				min_score = scores[1];
				min_pose = pose;
			}
		} else {
			++curr_counts;
			if ( did_move ) pose = stored_pose_;
		}

		if ( n % t_jump_interval == 0 && tempering.t_jump() ) {
			++n_t_jumps_accept;
			fill_data( data[temp_id], curr_counts, scores );
                        hist_list[temp_id].add( scores[1], curr_counts );
                        curr_counts = 1;
			set_gaussian_stdevs( sampler, chi_sampler, standard_sampler, 
				tempering, total_len, total_sampled, is_free);
                        temp_id = tempering.temp_id();
                }
		if ( (n % dump_interval) == 0 && option[dump_pdb]() ) {
			std::stringstream str;
			str << option[out_prefix]() << "_" << n << ".pdb";
			pose.dump_pdb( str.str() );
		}
		if ( (n % dump_interval) == 0 && option[dump_silent]() ) {
			std::stringstream tag;
			tag << option[out_prefix]() << "_" << n;
			BinarySilentStruct silent( pose, tag.str() );
			silent_file_data.write_silent_struct( silent, silent_file, false /*write score only*/ );
		}

	}


	std::cout << "Score for the end pose" << std::endl;
	scorefxn->show( pose );
	std::stringstream endstr, minstr;
	endstr << option[out_prefix]() << "_end.pdb";
	minstr << option[out_prefix]() << "_min.pdb";
	pose.dump_pdb( endstr.str() );
	min_pose.dump_pdb( minstr.str() );
	std::cout << "Score for the min pose" << std::endl;
	scorefxn->show( min_pose );

	std::cout << "n_cycles: " << n_cycle_ << std::endl;
	std::cout << "Total accept rate: " << double( n_accept_total ) / n_cycle_
		<< std::endl;
	std::cout << "Backbone accept rate: " << double( n_accept_backbone ) / ((n_cycle_ / 2) - (n_cycle_ / 10)) << std::endl;
	std::cout << "Chi accept rate: " << double( n_accept_chi ) / (n_cycle_ / 2) << std::endl;
	std::cout << "Standard accept rate: " << double( n_accept_standard ) / (n_cycle_ / 10) << std::endl;
	std::cout << "Temp jump accept rate: " << double( n_t_jumps_accept ) / (n_cycle_ / t_jump_interval) << std::endl;

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
	bb_tor_file.close();
	chi_tor_file.close();


// Some simple text data files
	std::stringstream str, str2, str3;
	str << option[out_prefix]() << "_data.txt";
	str2 << option[out_prefix]() << "_hist.txt";
	str3 << option[out_prefix]() << "_scores.txt";
	std::ofstream datafile ( (str.str()).c_str() );
	for ( Size i = 1; i <= data[1].size(); i+=2 ) {
		datafile << data[1][i] << "   " << data[1][i+1] << "\n";}
	datafile.close();
	std::ofstream histfile ( (str2.str()).c_str() );
	utility::vector1<Size> const & hist( hist_list[1].get_hist() );
	for ( Size i = 1; i<=hist.size(); i++){
		histfile<< hist[i] << "\n"; }
	histfile.close();
	utility::vector1<Real> const & escores( hist_list[1].get_scores() );
	std::ofstream scorefile ( (str3.str()).c_str() );
	for ( Size i = 1; i<=escores.size(); i++){
		scorefile<< escores[i] << "\n"; }
	scorefile.close();


        Real const time_in_test = static_cast<Real>( clock() - time_start )
                        / CLOCKS_PER_SEC;
        std::cout << "Time in sampler: " <<  time_in_test << std::endl;



//	// tag
//	std::string tag = tag_from_pose( pose );
//	BinarySilentStruct s( pose, tag );
//
//
//	std::cout << "Outputting " << tag << " to silent file: " << silent_file << std::endl;
//	silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );

}



///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	thermal_sampler();

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
    try {
        using namespace basic::options;

        std::cout << std::endl << "Basic usage:  " << argv[0] << "  -s <pdb file> " << std::endl;
        std::cout              << "              " << argv[0] << "  -in:file:silent <silent file> " << std::endl;
        std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

				utility::vector1< int > null_int_vector;
				utility::vector1< core::Real > null_real_vector;
				utility::vector1< Size > blank_size_vector;
				utility::vector1< std::string > blank_string_vector;
				option.add_relevant( score::weights );
				option.add_relevant( in::file::s );
				option.add_relevant( in::file::silent );
				option.add_relevant( in::file::tags );
				option.add_relevant( in::file::fasta );
				option.add_relevant( in::file::input_res );
				option.add_relevant( full_model::cutpoint_open );
				option.add_relevant( score::weights );
				NEW_OPT( out_prefix, "prefix for the out file", "thermal" );
				NEW_OPT( sample_residues, "residues to sample", null_int_vector );
				NEW_OPT( free_residues, "residues that are 'free,' affects stdev of chi sampler", null_int_vector );
				NEW_OPT( n_cycle, "cycle number for Random sampling", 0 );
				NEW_OPT( angle_range_bb, "range bb torsions are allowed to move", 20 );
				NEW_OPT( angle_range_free_bb, "range free bb torsions are allowed to move", 180 );
				NEW_OPT( angle_range_chi, "range chi torsions are allowed to move", 20 );
				NEW_OPT( angle_range_free_chi, "range chi torsions are allowed to move", 180 );
				NEW_OPT( chi_stdev, "standard deviation for chi sampler", 20 );
				NEW_OPT( bb_stdev, "standard deviation for backbone sampler", 1 ); 
				NEW_OPT( standard_bb_stdev, "standard deviation for standard backbone sampler", 1 ); 
				NEW_OPT( dump_pdb, "Dump pdb files", false );
				NEW_OPT( dump_silent, "Dump structures to a silent file", false );
				NEW_OPT( out_torsions, "Print out torsion angles", false );
				NEW_OPT( temps, "Simulated tempering temperatures", null_real_vector );
				NEW_OPT( st_weights, "Simulated tempering weights", null_real_vector );
				NEW_OPT( dump_freq, "Frequency to dump pdb or silent files", 500 );

        ////////////////////////////////////////////////////////////////////////////
        // setup
        ////////////////////////////////////////////////////////////////////////////
        devel::init(argc, argv);
				option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_BASE" );
				option[ OptionKeys::chemical::patch_selectors ].push_back( "TERMINAL_PHOSPHATE" );
				option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_RNA_RESIDUE" );

        ////////////////////////////////////////////////////////////////////////////
        // end of setup
        ////////////////////////////////////////////////////////////////////////////
        protocols::viewer::viewer_main( my_main );
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cout << "caught exception " << e.msg() << std::endl;
				return -1;
    }
}

