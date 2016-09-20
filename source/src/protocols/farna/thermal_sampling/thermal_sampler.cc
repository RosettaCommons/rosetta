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
#include <core/init/init.hh>
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
#include <core/io/rna/RNA_DataReader.hh>
#include <core/pose/PDBInfo.hh>

#include <protocols/stepwise/sampler/rna/RNA_MC_Suite.hh>
#include <protocols/stepwise/sampler/rna/RNA_MC_MultiSuite.hh>
#include <protocols/moves/SimulatedTempering.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/stepwise/sampler/rna/RNA_MC_KIC_Sampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_KIC_Sampler.hh>
#include <protocols/farna/thermal_sampling/util.hh>
#include <protocols/farna/thermal_sampling/thermal_sampler.hh>

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
using namespace basic::options::OptionKeys::rna::farna::thermal_sampling;
using namespace protocols::farna::thermal_sampling;
using utility::vector1;

namespace protocols {
namespace farna {
namespace thermal_sampling {


//////////////////////////////////////////////////////////////////////////////
utility::vector1<core::Real> get_torsions(
	utility::vector1<core::id::TorsionID> & torsion_ids,
	const Pose & pose
) {
	utility::vector1<core::Real> curr_torsions;
	for ( Size i = 1; i <= torsion_ids.size(); ++i ) {
		curr_torsions.push_back( pose.torsion( torsion_ids[i] ) );
	}
	return curr_torsions;
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
		if ( is_free[i] ) {
			chi_sampler[i]->set_gaussian_stdev( free_chi_stdev );
		} else chi_sampler[i]->set_gaussian_stdev( chi_stdev );
	}
	standard_bb_sampler.set_gaussian_stdev( standard_bb_stdev );
}

///////////////////////////////////////////////////////////////////////////////
void
thermal_sampler( Pose & pose )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring::rna::data;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose::full_model_info;
	using namespace protocols::stepwise::modeler;

	using namespace protocols::stepwise::sampler::rna;
	using namespace protocols::moves;
	using namespace core::id;
	using namespace protocols::stepwise::sampler;
	using namespace protocols::stepwise::setup;

	clock_t const time_start( clock() );

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
	pose::Pose start_pose;

	if ( !option[ in::file::silent ].user() ) cleanup( pose );

	if ( !full_model_info_defined( pose ) || option[ in::file::fasta ].user() ) {
		fill_full_model_info_from_command_line( pose, other_poses ); // only does something if -in:file:fasta specified.
	}

	// do it
	scorefxn->show( pose );

	Size const total_len( pose.total_residue() );
	bool const sample_near_a_form( false );

	//Setting up the internal move sampler
	pose::PoseOP ref_pose(new pose::Pose(pose));

	utility::vector1<int> residues( option[sample_residues]() );
	Size const total_sampled( residues.size() );
	utility::vector1<int> free_rsd( option[free_residues]() );
	utility::vector1<bool> is_free;
	std::sort( residues.begin(), residues.end() );
	std::cout << "Sample residues: " << residues <<std::endl;

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
		suite_sampler->set_angle_range_from_init_torsions( option[angle_range_bb]() );
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
		chi_torsion->set_angle_range( init_torsion - option[angle_range_chi]() , init_torsion + option[angle_range_chi]() );
		chi_sampler.push_back( chi_torsion );
	}

	//Set up the regular torsion samplers, necessary for free energy comparisons
	//Don't sample chi torsions here because there is a separate chi sampler
	RNA_MC_MultiSuite standard_sampler;
	for ( Size i = 1; i <= residues.size(); ++i ) {
		if ( i == 1 || residues[i] != residues[i-1]+1 ) {
			//create samplers for [i]-1 and [i]
			RNA_MC_SuiteOP suite_sampler_1( new RNA_MC_Suite( residues[i] - 1 ) );
			suite_sampler_1->set_init_from_pose( pose );
			suite_sampler_1->set_sample_bb( true );
			suite_sampler_1->set_sample_lower_nucleoside( false );
			suite_sampler_1->set_sample_upper_nucleoside( false );
			suite_sampler_1->set_pucker_flip_rate( 0 );
			suite_sampler_1->set_sample_near_a_form( sample_near_a_form );
			suite_sampler_1->set_angle_range_from_init( option[angle_range_bb]() );
			standard_sampler.add_external_loop_rotamer( suite_sampler_1 );
		}
		RNA_MC_SuiteOP suite_sampler( new RNA_MC_Suite( residues[i] ) );
		suite_sampler->set_init_from_pose( pose );
		suite_sampler->set_sample_bb( true );
		suite_sampler->set_sample_lower_nucleoside( false );
		suite_sampler->set_sample_upper_nucleoside( false );
		suite_sampler->set_pucker_flip_rate( 0 );
		suite_sampler->set_sample_near_a_form( sample_near_a_form );
		suite_sampler->set_angle_range_from_init( option[angle_range_bb]() );
		standard_sampler.add_external_loop_rotamer( suite_sampler );
	}
	standard_sampler.init();

	//Make a list of all the backbone torsion IDs that are sampled (This is only complete if all the residues are consecutive)
	utility::vector1<TorsionID> bb_torsion_ids;
	/*
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
	*/

	// OK, how do we do this better? We need EZ from -1 of any in set, and ABG from + 1
	// and eliminate dupes
	// Trivial insight: explicit EZ addition only needed if residues[i-1] != residues[i]-1 isn't
	// BTW we sort residues to avoid any problems there.

	for ( Size i = 1; i <= residues.size(); ++i ) {
		if ( i == 1 || residues[ i - 1 ] != residues[ i ] - 1 ) {
			bb_torsion_ids.push_back( TorsionID( residues[1]-1, id::BB, EPSILON ) );
			bb_torsion_ids.push_back( TorsionID( residues[1]-1, id::BB, ZETA ) );
		}
		bb_torsion_ids.push_back( TorsionID( residues[i], id::BB, ALPHA) );
		bb_torsion_ids.push_back( TorsionID( residues[i], id::BB, BETA) );
		bb_torsion_ids.push_back( TorsionID( residues[i], id::BB, GAMMA) );
		bb_torsion_ids.push_back( TorsionID( residues[i], id::BB, EPSILON) );
		bb_torsion_ids.push_back( TorsionID( residues[i], id::BB, ZETA) );
		if ( i == residues.size() || residues[ i + 1 ] != residues[ i ] + 1 ) {
			bb_torsion_ids.push_back( TorsionID( residues.back()+1, id::BB, ALPHA ) );
			bb_torsion_ids.push_back( TorsionID( residues.back()+1, id::BB, BETA ) );
			bb_torsion_ids.push_back( TorsionID( residues.back()+1, id::BB, GAMMA ) );
		}
	}

	bool const is_save_scores( true );

	// AMW: this is where MC_run starts, effectively, in recces_turner
	// Set up temperatures and ST weights
	utility::vector1<Real> const & temps_( option[ temps ]() );
	runtime_assert( temps_.size() != 0 );

	utility::vector1<Real> weights_;
	utility::vector1<Real> const & orig_weights( option[ st_weights ]() );
	if ( temps_.size() != orig_weights.size() ) {
		weights_.push_back( 0 );
	}
	weights_.insert( weights_.end(), orig_weights.begin(),
		orig_weights.end() );
	runtime_assert( temps_.size() == weights_.size() );

	Size curr_counts( 1 );
	utility::vector1<float> scores;
	update_scores( scores, pose, scorefxn );
	utility::vector1<float> const null_arr_;
	utility::vector1<utility::vector1<float> > data( temps_.size(), null_arr_ );

	Real const min( -100.05 ), max( 100.05 ), spacing( 0.1 );
	Histogram null_hist( min, max, spacing);
	utility::vector1<Histogram> hist_list( temps_.size(), null_hist );

	// Min-score pose
	Pose min_pose = pose;
	Real min_score( 99999 );

	SimulatedTempering tempering( pose, scorefxn, temps_, weights_ );
	MonteCarlo mc(pose, *scorefxn, 1. );
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
	std::ofstream bb_tor_file( name.str() );
	std::ofstream chi_tor_file( name2.str() );


	std::stringstream initstr;
	initstr << option[out_prefix]() << "_init.pdb";
	pose.dump_pdb( initstr.str() );

	set_gaussian_stdevs( sampler, chi_sampler, standard_sampler,
		tempering, total_len, total_sampled, is_free );

	Pose stored_pose_ = pose;
	// Vector center_vector = Vector( 0.0 );
	// protocols::viewer::add_conformation_viewer ( pose.conformation(), "current", 700, 700, false, false , center_vector );

	// Main sampling cycle
	for ( Size n = 1; n <= n_cycle_; ++n ) {
		did_move = true;
		// Pick from among samplers
		index = numeric::random::rg().random_range(1,sampler.size());
		// Sampler updates
		if ( (n % 10) == 0 ) {
			// We need to update the torsions in case they were updated by the other samplers
			standard_sampler.set_angle( pose );
			++standard_sampler;
			standard_sampler.apply( pose );
		} else if ( (n % 2) == 0 ) {
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
			} else if ( (n % 2) == 0 ) {
				sampler[index]->update( pose ); // don't think this is necessary
				++n_accept_backbone;
			} else {
				chi_sampler[index]->update();
				++n_accept_chi;
			}
			if ( option[out_torsions]() ) {
				utility::vector1<core::Real> bb_torsions = get_torsions(bb_torsion_ids, pose);
				utility::vector1<core::Real> chi_torsions = get_torsions(chi_torsion_ids, pose);
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

	for ( Size i = 1; i <= temps_.size(); ++i ) {
		if ( is_save_scores ) {
			std::ostringstream oss;
			oss << option[out_prefix]() << '_' << std::fixed << std::setprecision(2)
				<< temps_[i] << ".bin.gz";
			Size const data_dim2( scorefxn->get_nonzero_weighted_scoretypes().size() + 2 );
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
		datafile << data[1][i] << "   " << data[1][i+1] << "\n";
	}
	datafile.close();
	std::ofstream histfile ( (str2.str()).c_str() );
	utility::vector1<Size> const & hist( hist_list[1].get_hist() );
	for ( Size i = 1; i<=hist.size(); i++ ) {
		histfile<< hist[i] << "\n";
	}
	histfile.close();
	utility::vector1<Real> const & escores( hist_list[1].get_scores() );
	std::ofstream scorefile ( (str3.str()).c_str() );
	for ( Size i = 1; i<=escores.size(); i++ ) {
		scorefile<< escores[i] << "\n";
	}
	scorefile.close();


	Real const time_in_test = static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC;
	std::cout << "Time in sampler: " <<  time_in_test << std::endl;
}

}
}
}

