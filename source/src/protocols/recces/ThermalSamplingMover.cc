// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/ThermalSamplingMover.cc
/// @brief Use a simulated tempering simulation to refine a pose
/// @author Andy Watkins (amw579@nyu.edu)

// Unit headers
#include <protocols/recces/ThermalSamplingMover.hh>
#include <protocols/recces/ThermalSamplingMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/data/RNA_ChemicalMappingEnergy.hh>
#include <core/scoring/rna/data/RNA_DataInfo.hh>
#include <core/kinematics/MoveMap.hh>
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
#include <core/id/TorsionID.hh>
#include <core/pose/PDBInfo.hh>

#include <protocols/recces/sampler/rna/MC_RNA_Suite.hh>
#include <protocols/recces/sampler/rna/MC_RNA_MultiSuite.hh>
#include <protocols/moves/SimulatedTempering.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/recces/sampler/rna/MC_RNA_KIC_Sampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_KIC_Sampler.hh>
#include <protocols/recces/util.hh>
#include <protocols/recces/thermal_sampler.hh>

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/recces.OptionKeys.gen.hh>
#include <basic/options/keys/recces.OptionKeys.gen.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.recces.ThermalSamplingMover" );

namespace protocols {
namespace recces {

using namespace core::chemical::rna;
using namespace core;
using namespace core::pose;

/// @brief at least in the helix-strand case, will successfully distinguish a bad bond length
/// from actual breaks
Size find_likely_first_chain_ending( Pose const & pose ) {
	for ( Size ii = 1; ii <= pose.size() - 1; ++ii ) {
		// Where is the bond length from ii to ii+1 too long?
		if ( pose.residue( ii ).xyz( "O3'" ).distance_squared( pose.residue( ii + 1 ) .xyz( "P" ) ) > 4 ) {
			return ii;
		}
	}
	return pose.size();
}

/////////////////////
/// Constructors ///
/////////////////////

/// @brief Default constructor
ThermalSamplingMover::ThermalSamplingMover():
	protocols::moves::Mover( ThermalSamplingMover::mover_name() ),
	residues_( basic::options::option[ basic::options::OptionKeys::recces::sample_residues ]() ),
	free_rsd_( basic::options::option[ basic::options::OptionKeys::recces::free_residues ]() ),
	loop_rsd_( basic::options::option[ basic::options::OptionKeys::recces::loop_residues ]() ),
	n_cycle_( basic::options::option[ basic::options::OptionKeys::recces::n_cycle ] ),
	dump_pdb_( basic::options::option[ basic::options::OptionKeys::recces::dump_pdb ]() ),
	dump_silent_( basic::options::option[ basic::options::OptionKeys::recces::dump_silent ]() ),
	angle_range_chi_( basic::options::option[ basic::options::OptionKeys::recces::angle_range_chi ]() ),
	angle_range_bb_( basic::options::option[ basic::options::OptionKeys::recces::angle_range_bb ]() ),
	temps_( basic::options::option[ basic::options::OptionKeys::recces::temps ]() ),
	st_weights_( basic::options::option[ basic::options::OptionKeys::recces::st_weights ]() )
{
	total_sampled_ = residues_.size();
	// Insert zero weights at the beginning if there is a mismatch.
	while ( temps_.size() != st_weights_.size() ) {
		utility::vector1<Real> const orig_weights( st_weights_ );
		st_weights_.clear();
		st_weights_.push_back( 0 );
		st_weights_.insert( st_weights_.end(), orig_weights.begin(), orig_weights.end() );
	}
	runtime_assert( temps_.size() == st_weights_.size() );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
ThermalSamplingMover::ThermalSamplingMover( ThermalSamplingMover const & src ):
	protocols::moves::Mover( src )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
ThermalSamplingMover::~ThermalSamplingMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
ThermalSamplingMover::apply( core::pose::Pose & pose ) {
	runtime_assert( temps_.size() != 0 );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::recces;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose::full_model_info;
	using namespace protocols::stepwise::modeler;

	using namespace protocols::recces::sampler::rna;
	using namespace protocols::moves;
	using namespace core::id;
	using namespace protocols::recces::sampler;
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
	std::string const silent_file = option[ out::file::silent ]();
	SilentFileData silent_file_data;

	// if trying to compute stem RMSD
	pose::Pose start_pose;

	// do it
	scorefxn->show( pose );

	Size const total_len( pose.size() );
	bool const sample_near_a_form( false );

	//Setting up the internal move sampler
	pose::PoseOP ref_pose(new pose::Pose(pose));

	utility::vector1<bool> is_free;
	std::sort( residues_.begin(), residues_.end() );
	TR << "Sample residues: " << residues_ <<std::endl;

	//Assign "free" residues (get bigger gaussian stdev)
	for ( Size i = 1; i <= residues_.size(); ++i ) {
		if ( std::find( free_rsd_.begin(), free_rsd_.end(), residues_[i] ) != free_rsd_.end() ) {
			is_free.push_back( true );
		} else is_free.push_back( false );
	}

	//Assign "loop" residues (get bigger bb torsion range)
	utility::vector1<bool> is_loop;
	for ( Size i = 1; i <= residues_.size(); ++i ) {
		if ( std::find( loop_rsd_.begin(), loop_rsd_.end(), residues_[i] ) != loop_rsd_.end() ) {
			is_loop.push_back( true );
		} else is_loop.push_back( false );
	}

	utility::vector1<MC_RNA_KIC_SamplerOP> sampler;
	utility::vector1<MC_OneTorsionOP> chi_sampler;
	utility::vector1<core::id::TorsionID> chi_torsion_ids;
	if ( !recces_turner_mode_ ) {
		//Set up the internal move samplers
		for ( Size i = 1; i <= residues_.size(); ++i ) {
			Size const seqpos = residues_[ i ];
			// AMW TODO: make it so the MC_KIC_Sampler can handle terminal residue?
			if ( seqpos == pose.size() ) continue; // we need one torsion from i+1
			if ( seqpos == 1 ) continue; // no seqpos - 1

			MC_RNA_KIC_SamplerOP suite_sampler( new MC_RNA_KIC_Sampler( ref_pose, seqpos-1, seqpos ) );
			suite_sampler->init();
			suite_sampler->set_angle_range_from_init_torsions( is_loop[ i ] ? angle_range_loop_bb_ : angle_range_bb_ );
			sampler.push_back( suite_sampler );
		}

		//Set up the chi samplers
		core::Real init_torsion;
		for ( Size const seqpos : residues_ ) {
			// skip if movemap exists but dof off
			if ( mm_ && !mm_->get_chi( seqpos ) ) continue;

			core::id::TorsionID chi_ID (core::id::TorsionID( seqpos, id::CHI, 1));
			chi_torsion_ids.push_back( chi_ID );
			init_torsion = ref_pose->torsion( chi_ID );
			MC_OneTorsionOP chi_torsion( new MC_OneTorsion( chi_ID, ref_pose->torsion( chi_ID ) ) );
			chi_torsion->init();
			chi_torsion->set_angle_range( init_torsion - angle_range_chi_ , init_torsion + angle_range_chi_ );
			chi_sampler.push_back( chi_torsion );
		}
	}

	// AMW TODO: look at recces turner for clever use of upper and lower nucleoside

	//Set up the regular torsion samplers, necessary for free energy comparisons
	//Don't sample chi torsions here because there is a separate chi sampler
	MC_RNA_MultiSuite standard_sampler;

	// Preserve recces_turner mode until we can figure out how best to unify it.
	if ( recces_turner_mode_ ) {
		Size const len1 = find_likely_first_chain_ending( pose );
		for ( Size i = 1; i <= pose.size(); ++i ) {
			//bool const sample_near_a_form( bp_rsd.has_value( i ) );
			// bp_rsd == ! free_rsd in the recces_turner case, so...
			bool const sample_near_a_form( ! free_rsd_.has_value( i ) );
			if ( i == 1 || ( i > len1 && i != pose.size() ) ) {
				MC_RNA_SuiteOP suite_sampler( new MC_RNA_Suite( i ) );
				suite_sampler->set_sample_bb( i != 1 );
				suite_sampler->set_sample_lower_nucleoside( true );
				suite_sampler->set_sample_upper_nucleoside( false );
				suite_sampler->set_sample_near_a_form( sample_near_a_form );
				// Since we're trying to preserve classic recces. add another flag if we want (a_form_range)
				suite_sampler->set_a_form_range( 60.0 );
				// What was called "sampler" in recces_turner is "standard_sampler" here
				standard_sampler.add_external_loop_rotamer( suite_sampler );
			} else {
				MC_RNA_SuiteOP suite_sampler( new MC_RNA_Suite( i - 1 ) );
				suite_sampler->set_sample_bb( len1 == pose.size() || i != pose.size() );
				suite_sampler->set_sample_lower_nucleoside( false );
				suite_sampler->set_sample_upper_nucleoside( true );
				suite_sampler->set_sample_near_a_form( sample_near_a_form );
				// Since we're trying to preserve classic recces. add another flag if we want (a_form_range)
				suite_sampler->set_a_form_range( 60.0 );
				// What was called "sampler" in recces_turner is "standard_sampler" here
				standard_sampler.add_external_loop_rotamer( suite_sampler );
			}
		}
	} else {
		for ( Size i = 1; i <= residues_.size(); ++i ) {
			if ( residues_[i] != 1 && ( i == 1 || residues_[i] != residues_[i-1]+1 ) ) {
				//create samplers for [i]-1 and [i]

				// Can't just have a continue condition because we might want to set up the second sampler
				TR << "Going to attempt to sample " << residues_[i] - 1 << " with a MC_RNA_Suite and pose size is " << pose.size() << std::endl;
				MC_RNA_SuiteOP suite_sampler_1( new MC_RNA_Suite( residues_[i] - 1 ) );
				suite_sampler_1->set_init_from_pose( pose );
				suite_sampler_1->set_sample_bb( true );
				suite_sampler_1->set_sample_lower_nucleoside( false );
				suite_sampler_1->set_sample_upper_nucleoside( false );
				suite_sampler_1->set_pucker_flip_rate( 0 );
				suite_sampler_1->set_sample_near_a_form( sample_near_a_form );
				// We don't consider the possibility that this is a loop residue that deserves aggressive sampling, because
				// if it were, it would be getting added as a sampler in the clause below!
				//suite_sampler_1->set_angle_range_from_init( is_loop[ i ] ? angle_range_loop_bb_ : angle_range_bb_ );
				//suite_sampler_1->set_angle_range_from_init( angle_range_bb_ );
				standard_sampler.add_external_loop_rotamer( suite_sampler_1 );
			}

			// can't do this!
			if ( residues_[i] == pose.size() ) continue;

			TR << "Going to attempt to sample " << residues_[i] << " with a MC_RNA_Suite\tand pose size is " << pose.size() << std::endl;
			MC_RNA_SuiteOP suite_sampler( new MC_RNA_Suite( residues_[i] ) );
			suite_sampler->set_init_from_pose( pose );
			suite_sampler->set_sample_bb( true );
			suite_sampler->set_sample_lower_nucleoside( false );
			suite_sampler->set_sample_upper_nucleoside( false );
			suite_sampler->set_pucker_flip_rate( 0 );
			suite_sampler->set_sample_near_a_form( sample_near_a_form );
			//suite_sampler->set_angle_range_from_init( is_loop[ i ] ? angle_range_loop_bb_ : angle_range_bb_ );
			standard_sampler.add_external_loop_rotamer( suite_sampler );
		}
	}
	standard_sampler.init();

	//Make a list of all the backbone torsion IDs that are sampled
	// Since AMW changes this will always be complete
	// AMW: notably, this is NOT modified based on the above MoveMap conditionals
	// but, that's okay for now. If you're supplying a movemap, you are probably not using
	// the data-dumping functionality anyway (minimizer replacement mode)
	utility::vector1<core::id::TorsionID> bb_torsion_ids;

	for ( Size i = 1; i <= residues_.size(); ++i ) {
		if ( i == 1 || residues_[ i - 1 ] != residues_[ i ] - 1 ) {
			bb_torsion_ids.push_back( core::id::TorsionID( residues_[1]-1, id::BB, EPSILON ) );
			bb_torsion_ids.push_back( core::id::TorsionID( residues_[1]-1, id::BB, ZETA ) );
		}
		bb_torsion_ids.push_back( core::id::TorsionID( residues_[i], id::BB, ALPHA) );
		bb_torsion_ids.push_back( core::id::TorsionID( residues_[i], id::BB, BETA) );
		bb_torsion_ids.push_back( core::id::TorsionID( residues_[i], id::BB, GAMMA) );
		bb_torsion_ids.push_back( core::id::TorsionID( residues_[i], id::BB, EPSILON) );
		bb_torsion_ids.push_back( core::id::TorsionID( residues_[i], id::BB, ZETA) );
		if ( i == residues_.size() || residues_[ i + 1 ] != residues_[ i ] + 1 ) {
			bb_torsion_ids.push_back( core::id::TorsionID( residues_.back()+1, id::BB, ALPHA ) );
			bb_torsion_ids.push_back( core::id::TorsionID( residues_.back()+1, id::BB, BETA ) );
			bb_torsion_ids.push_back( core::id::TorsionID( residues_.back()+1, id::BB, GAMMA ) );
		}
	}

	bool const is_save_scores( true );

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

	SimulatedTempering tempering( pose, scorefxn, temps_, st_weights_ );
	MonteCarlo mc(pose, *scorefxn, 1. );
	Size n_accept_total( 0 ), n_accept_backbone( 0 );
	Size n_accept_chi( 0 ), n_accept_standard( 0 ), n_t_jumps_accept( 0 );
	Size const t_jump_interval( 10 );
	int const dump_interval( option[dump_freq]() );
	bool did_move;
	int index;
	Size temp_id( 1 );

	std::stringstream name, name2;
	name << option[out_prefix]() << "_bb_torsions.txt";
	name2 << option[out_prefix]() << "_chi_torsions.txt";
	std::ofstream bb_tor_file( name.str() );
	std::ofstream chi_tor_file( name2.str() );

	if ( dumping_app_ ) {
		std::stringstream initstr;
		initstr << option[out_prefix]() << "_init.pdb";
		pose.dump_pdb( initstr.str() );
	}

	set_gaussian_stdevs( sampler, chi_sampler, standard_sampler, tempering, total_len, total_sampled_, is_free );

	Pose stored_pose_ = pose;
	// Vector center_vector = Vector( 0.0 );
	// protocols::viewer::add_conformation_viewer ( pose.conformation(), "current", 700, 700, false, false , center_vector );

	//debug_assert( sampler.size() );

	// Main sampling cycle
	// AMW: Had to change the structure of the main loop so that all conditional
	// references to "index" were together
	Size pct = 0;
	for ( Size n = 1; n <= n_cycle_; ++n ) {

		if ( n % ( n_cycle_ / 100 ) == 0 ) {
			++pct;
			std::cout << pct << "% complete." << std::endl;
		}

		did_move = true;
		// Sampler updates
		bool boltz = false;

		if ( recces_turner_mode_ ) {
			standard_sampler.set_angle( pose );
			++standard_sampler;
			standard_sampler.apply( pose );
			if ( ( tempering.check_boltzmann( pose ) && did_move ) || n == n_cycle_ ) {
				stored_pose_ = pose;
				if ( is_save_scores ) fill_data( data[temp_id], curr_counts, scores );
				hist_list[temp_id].add( scores[1], curr_counts );
				update_scores( scores, pose, scorefxn );
				if ( n == n_cycle_ ) break;
				boltz = true;
			}
		} else {
			if ( (n % 10) == 0 ) {
				// We need to update the torsions in case they were updated by the other samplers
				standard_sampler.set_angle( pose );
				++standard_sampler;
				standard_sampler.apply( pose );
				if ( ( tempering.check_boltzmann( pose ) && did_move ) || n == n_cycle_ ) {
					stored_pose_ = pose;
					if ( is_save_scores ) fill_data( data[temp_id], curr_counts, scores );
					hist_list[temp_id].add( scores[1], curr_counts );
					update_scores( scores, pose, scorefxn );
					if ( n == n_cycle_ ) break;
					boltz = true;
				}
			} else if ( (n % 2) == 0 && sampler.size() > 0 ) {
				// Pick from among samplers
				index = numeric::random::rg().random_range(1,sampler.size());
				sampler[index]->next( pose ); //This function also updates the stored torsions
				sampler[index]->apply( pose );
				if ( !(sampler[index]->check_moved() ) ) {
					did_move = false;
				}
				if ( ( tempering.check_boltzmann( pose ) && did_move ) || n == n_cycle_ ) {
					stored_pose_ = pose;
					if ( is_save_scores ) fill_data( data[temp_id], curr_counts, scores );
					hist_list[temp_id].add( scores[1], curr_counts );
					update_scores( scores, pose, scorefxn );
					if ( n == n_cycle_ ) break;
					standard_sampler.update();
					++n_accept_standard;
					boltz = true;
				}
			} else if ( chi_sampler.size() > 0 ) {
				// Pick from among samplers
				index = numeric::random::rg().random_range(1,chi_sampler.size());
				++(*chi_sampler[index]);
				chi_sampler[index]->apply( pose );
				if ( ( tempering.check_boltzmann( pose ) && did_move ) || n == n_cycle_ ) {
					stored_pose_ = pose;
					if ( is_save_scores ) fill_data( data[temp_id], curr_counts, scores );
					hist_list[temp_id].add( scores[1], curr_counts );
					update_scores( scores, pose, scorefxn );
					if ( n == n_cycle_ ) break;
					sampler[index]->update( pose ); // don't think this is necessary
					++n_accept_backbone;
					chi_sampler[index]->update();
					++n_accept_chi;
					boltz = true;
				}
			}
		}

		if ( ( boltz && did_move ) || n == n_cycle_ ) {
			if ( option[out_torsions]() && !recces_turner_mode_ ) {
				utility::vector1<core::Real> bb_torsions = get_torsions(bb_torsion_ids, pose);
				utility::vector1<core::Real> chi_torsions = get_torsions(chi_torsion_ids, pose);
				bb_tor_file << curr_counts<< "\t"<< bb_torsions <<"\n";
				chi_tor_file << curr_counts<< "\t"<< chi_torsions <<"\n";
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
			set_gaussian_stdevs( sampler, chi_sampler, standard_sampler, tempering, total_len, total_sampled_, is_free);
			temp_id = tempering.temp_id();
		}
		if ( (n % dump_interval) == 0 && dump_pdb_ ) {
			std::stringstream str;
			str << option[out_prefix]() << "_" << n << ".pdb";
			pose.dump_pdb( str.str() );
		}
		if ( (n % dump_interval) == 0 && dump_silent_ ) {
			std::stringstream tag;
			tag << option[out_prefix]() << "_" << n;
			BinarySilentStruct silent( pose, tag.str() );
			silent_file_data.write_silent_struct( silent, silent_file, false /*write score only*/ );
		}
	}

	TR << "Score for the end pose" << std::endl;
	scorefxn->show( pose );
	if ( dumping_app_ ) {
		std::stringstream endstr, minstr;
		endstr << option[out_prefix]() << "_end.pdb";
		minstr << option[out_prefix]() << "_min.pdb";
		pose.dump_pdb( endstr.str() );
		min_pose.dump_pdb( minstr.str() );
	}
	TR << "Score for the min pose" << std::endl;
	scorefxn->show( min_pose );

	TR << "n_cycles: " << n_cycle_ << std::endl;
	TR << "Total accept rate: " << double( n_accept_total ) / n_cycle_ << std::endl;
	TR << "Backbone accept rate: " << double( n_accept_backbone ) / ((n_cycle_ / 2) - (n_cycle_ / 10)) << std::endl;
	TR << "Chi accept rate: " << double( n_accept_chi ) / (n_cycle_ / 2) << std::endl;
	TR << "Standard accept rate: " << double( n_accept_standard ) / (n_cycle_ / 10) << std::endl;
	TR << "Temp jump accept rate: " << double( n_t_jumps_accept ) / (n_cycle_ / t_jump_interval) << std::endl;

	// If we're being called from an app that wants to dump data after sampling
	// like recces_turner or thermal_sampler
	if ( dumping_app_ ) {
		for ( Size i = 1; i <= temps_.size(); ++i ) {
			if ( is_save_scores ) {
				std::ostringstream oss;
				oss << option[out_prefix]() << '_' << std::fixed << std::setprecision(2)
					<< temps_[i] << ".bin.gz";
				Size const data_dim2( get_scoretypes().size() + 2 );
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
		// In rare circumstances, I guess data[1] has odd parity -- saw a bounds issue
		for ( Size i = 1; i <= data[1].size() - 1; i+=2 ) {
			datafile << data[1][i] << "\t " << data[1][i+1] << "\n";
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
	}

	Real const time_in_test = static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC;
	std::cout << "Time in sampler: " << time_in_test << std::endl;


	// If we're using this as a mover in a larger protocol, the most reasonable
	// choice for the final pose is the min score pose
	pose = min_pose;
}

/// @brief Sets thermal sampling level given a MoveMap
void
ThermalSamplingMover::set_residue_sampling_from_pose_and_movemap(
	core::pose::Pose const & pose,
	core::kinematics::MoveMap const & mm
) {
	// Stupid starting point: if delta and epsilon are free you're a free_rsd
	// otherwise you're just a rsd
	residues_.clear();
	free_rsd_.clear();

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( mm.get( core::id::TorsionID( ii, id::BB, BETA ) ) && mm.get( core::id::TorsionID( ii, id::BB, GAMMA ) ) ) {
			residues_.push_back( ii );
			if ( mm.get( core::id::TorsionID( ii, id::BB, DELTA ) ) && mm.get( core::id::TorsionID( ii, id::BB, EPSILON ) ) ) {
				free_rsd_.push_back( ii );
			}
		}
	}
	TR << "Movemap is " << mm << std::endl;
	TR << "Residues is " << residues_ << std::endl;
	TR << "Free residues is " << free_rsd_ << std::endl;
	// total_sampled_ needs to be updated
	total_sampled_ = residues_.size();

	mm_ = std::make_shared< core::kinematics::MoveMap >( mm );
}


////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
ThermalSamplingMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
ThermalSamplingMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new ThermalSamplingMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
ThermalSamplingMover::clone() const
{
	return protocols::moves::MoverOP( new ThermalSamplingMover( *this ) );
}

// XRW TEMP std::string
// XRW TEMP ThermalSamplingMover::get_name() const
// XRW TEMP {
// XRW TEMP  return ThermalSamplingMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ThermalSamplingMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "ThermalSamplingMover";
// XRW TEMP }

void
ThermalSamplingMover::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, ThermalSamplingMover const & mover )
{
	mover.show(os);
	return os;
}

/////////////// Creator ///////////////

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ThermalSamplingMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new ThermalSamplingMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ThermalSamplingMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return ThermalSamplingMover::mover_name();
// XRW TEMP }

std::string ThermalSamplingMover::get_name() const {
	return mover_name();
}

std::string ThermalSamplingMover::mover_name() {
	return "ThermalSamplingMover";
}

void ThermalSamplingMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	//No attributes!
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string ThermalSamplingMoverCreator::keyname() const {
	return ThermalSamplingMover::mover_name();
}

protocols::moves::MoverOP
ThermalSamplingMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ThermalSamplingMover );
}

void ThermalSamplingMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ThermalSamplingMover::provide_xml_schema( xsd );
}



} //recces
} //protocols

