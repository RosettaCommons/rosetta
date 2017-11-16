// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/denovo/recces/ThermalMinimizer.cc
/// @brief Use a simulated tempering simulation to refine a pose
/// @author Andy Watkins (amw579@nyu.edu)

// Unit headers
#include <protocols/recces/scratch/ThermalMinimizer.hh>
#include <protocols/recces/scratch/ThermalMinimizerCreator.hh>

// Core headers
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/rms_util.hh>
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
#include <protocols/recces/scratch/thermal_sampler.hh>

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/recces.OptionKeys.gen.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.recces.scratch.ThermalMinimizer" );

using namespace protocols::stepwise::sampler::rna;
using namespace protocols::recces::sampler;
using namespace protocols::recces::sampler::rna;

using namespace core;

namespace protocols {
namespace recces {
namespace scratch {


void set_gaussian_stdevs(
	utility::vector1<MC_RNA_KIC_SamplerOP> & internal_bb_sampler,
	utility::vector1<MC_OneTorsionOP> & chi_sampler,
	Real const temp,
	Size const //total_sampled
) {
	Real internal_bb_stdev( 0.1 * pow( temp, 0.25 ) + 0.1);
	Real chi_stdev( 5 * pow( temp , 0.5) + 15 );
	if ( temp < 0 ) {
		internal_bb_stdev = 0.5 ;
		chi_stdev = -1 ;
	}
	for ( auto & sampler : internal_bb_sampler ) {
		sampler->set_gaussian_stdev( internal_bb_stdev );
	}
	for ( auto & sampler : chi_sampler ) {
		if ( sampler->get_torsion_id().type() == core::id::BB ) {
			sampler->set_gaussian_stdev( internal_bb_stdev );
		} else {
			sampler->set_gaussian_stdev( chi_stdev );
		}
	}
}

void set_gaussian_stdevs(
	utility::vector1<MC_RNA_KIC_SamplerOP> & internal_bb_sampler,
	utility::vector1<MC_OneTorsionOP> & chi_sampler,
	MC_RNA_MultiSuite & standard_bb_sampler,
	Real const temp,
	Size const total_sampled
) {
	// Shrinking variation by 2/3, expected need for 3/2 as many cycles vaguely.
	//Real internal_bb_stdev( 0.1 * pow( temp, 0.25 ) + 0.1);
	//Real chi_stdev( 5 * pow( temp , 0.5) + 15 );
	//Real standard_bb_stdev( 8 * pow( temp, 0.5 ) / ( 3 * total_sampled ) );
	Real internal_bb_stdev( 0.067 * pow( temp, 0.25 ) + 0.067);
	Real chi_stdev( 3.33 * pow( temp , 0.5) + 10 );
	Real standard_bb_stdev( 16 * pow( temp, 0.5 ) / ( 9 * total_sampled ) );
	if ( temp < 0 ) {
		internal_bb_stdev = 0.5 ;
		chi_stdev = -1 ;
		standard_bb_stdev = -1 ;
	}
	for ( Size i = 1; i <= internal_bb_sampler.size(); ++i ) {
		internal_bb_sampler[i]->set_gaussian_stdev( internal_bb_stdev );
	}
	for ( Size i = 1; i <= chi_sampler.size(); ++i ) {
		chi_sampler[i]->set_gaussian_stdev( chi_stdev );
	}
	standard_bb_sampler.set_gaussian_stdev( standard_bb_stdev );
}


/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ThermalMinimizer::ThermalMinimizer():
	protocols::moves::Mover( ThermalMinimizer::mover_name() ),
	n_cycle_( basic::options::option[ basic::options::OptionKeys::recces::n_cycle ] ),
	angle_range_chi_( basic::options::option[ basic::options::OptionKeys::recces::thermal_sampling::angle_range_chi ]() ),
	angle_range_bb_( basic::options::option[ basic::options::OptionKeys::recces::thermal_sampling::angle_range_bb ]() ),
	kic_sampling_( true ),
	output_min_pose_( true )
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
ThermalMinimizer::ThermalMinimizer( ThermalMinimizer const & src ):
	protocols::moves::Mover( src )
{
}

// ii - 1 === moving, ii === chainbreak
bool
mm_compatible_with_kic( core::kinematics::MoveMapOP mm, Size const ii ) {

	// At the moment, let's just assume KIC is always wrong.
	//return false;

	using namespace core::chemical::rna;
	using namespace core::id;

	if ( !mm->get( TorsionID( ii - 1, BB, EPSILON ) ) ) return false;
	if ( !mm->get( TorsionID( ii - 1, BB, ZETA ) ) ) return false;
	if ( !mm->get( TorsionID( ii, BB, ALPHA ) ) ) return false;
	if ( !mm->get( TorsionID( ii, BB, BETA ) ) ) return false;
	if ( !mm->get( TorsionID( ii, BB, GAMMA ) ) ) return false;
	if ( !mm->get( TorsionID( ii, BB, EPSILON ) ) ) return false;
	if ( !mm->get( TorsionID( ii, BB, ZETA ) ) ) return false;
	if ( !mm->get( TorsionID( ii + 1, BB, ALPHA ) ) ) return false;
	if ( !mm->get( TorsionID( ii + 1, BB, BETA ) ) ) return false;
	if ( !mm->get( TorsionID( ii + 1, BB, GAMMA ) ) ) return false;

	return true;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
ThermalMinimizer::~ThermalMinimizer(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
ThermalMinimizer::apply( core::pose::Pose & pose ) {

	if ( ! mm_ ) { TR.Warning << "No movemap provided! Returning." << std::endl; }
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

	using namespace protocols::stepwise::sampler::rna;
	using namespace protocols::moves;
	using namespace core::id;
	using namespace protocols::stepwise::sampler;
	using namespace protocols::stepwise::setup;

	utility::vector1< Size > residues_sampled;
	Size total_sampled = 0;
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		// Check for a free, sample-able dof
		if ( mm_->get( core::id::TorsionID( ii, id::BB, 1 ) )
				|| mm_->get( core::id::TorsionID( ii, id::BB, 2 ) )
				|| mm_->get( core::id::TorsionID( ii, id::BB, 3 ) )
				|| mm_->get( core::id::TorsionID( ii, id::BB, 4 ) )
				|| mm_->get( core::id::TorsionID( ii, id::BB, 5 ) )
				|| mm_->get( core::id::TorsionID( ii, id::BB, 6 ) )
				|| mm_->get( core::id::TorsionID( ii, id::CHI, 1 ) )
				|| mm_->get( core::id::TorsionID( ii, id::CHI, 4 ) )
				) {
			++total_sampled;
			residues_sampled.push_back( ii );
		}
	}

	if ( total_sampled == 0 ) return;

	clock_t const time_start( clock() );

	// do it
	TR << "Score for the start pose: " << (*score_fxn_)( pose ) << std::endl;
	score_fxn_->show( TR.Debug, pose );
	Real const start_score = pose.energies().total_energy();

	auto const vec = figure_out_moving_cutpoints_closed_from_moving_res( pose, residues_sampled );
	TR << "Residues sampled " << residues_sampled << std::endl;
	//TR << "moving cutpoint res  " << vec << std::endl;

	//Size const total_len( pose.size() );
	//bool const sample_near_a_form( false );

	//Setting up the internal move sampler
	pose::PoseOP ref_pose(new pose::Pose(pose));

	MC_RNA_MultiSuite standard_sampler;
	utility::vector1<MC_RNA_KIC_SamplerOP> kic_sampler;
	utility::vector1<MC_OneTorsionOP> tors_sampler;
	//utility::vector1<MC_RNA_SugarOP> sugar_sampler;

	std::set< Size > sampled_by_kic;

	if ( kic_sampling_ ) {
		//Set up the internal move samplers
		for ( Size const seqpos : vec ) {
			if ( seqpos < 2 || seqpos > pose.size() - 1 ) {
				TR.Warning << "Oddly, it seems " << seqpos << " was marked as KIC-compatible " << std::endl;
				continue;
			}
			//}Size ii = 2; ii <= pose.size() - 1; ++ii ) {
			// use the movemap, luke
			//TR << "Do I work with " << seqpos << "? " << std::endl;
			//TR << pose.residue_type( seqpos ) << std::endl;
			//mm_->show();
			if ( !mm_compatible_with_kic( mm_, seqpos ) ) continue;

			MC_RNA_KIC_SamplerOP suite_sampler( new MC_RNA_KIC_Sampler( ref_pose, seqpos-1, seqpos, false /*change_ft*/ ) );
			suite_sampler->init();
			suite_sampler->set_angle_range_from_init_torsions( angle_range_bb_ );
			//TR << "Going to attempt to sample " << seqpos << " with a MC_RNA_KIC\tand pose size is " << pose.size() << std::endl;
			kic_sampler.push_back( suite_sampler );
			// Keep track of what to skip for suites
			sampled_by_kic.insert(seqpos-1);
			sampled_by_kic.insert(seqpos);
		}
	}

	//Set up the torsion samplers
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		auto tid = core::id::TorsionID( ii, id::CHI, 1 );
		if ( mm_->get( tid ) ) {
			Real init_torsion = ref_pose->torsion( tid );
			MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
			tors_mc->init();
			tors_mc->set_angle_range( init_torsion - angle_range_chi_ , init_torsion + angle_range_chi_ );
			tors_sampler.push_back( tors_mc );
		}
		tid = core::id::TorsionID( ii, id::CHI, 4 );
		if ( mm_->get( tid ) ) {
			Real init_torsion = ref_pose->torsion( tid );
			MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
			tors_mc->init();
			tors_mc->set_angle_range( init_torsion - angle_range_chi_ , init_torsion + angle_range_chi_ );
			tors_sampler.push_back( tors_mc );
		}

		// Don't set up backbone OneTorsion samplers for which KIC could do the job.
		if ( sampled_by_kic.find( ii ) != sampled_by_kic.end() ) continue;

		continue;
		/*
		if ( !pose.residue_type(ii).has_variant_type( core::chemical::VIRTUAL_PHOSPHATE ) ) {
		// skip if movemap exists but dof off
		tid = core::id::TorsionID( ii, id::BB, 1 );
		if ( mm_->get( tid ) ) {
		Real init_torsion = ref_pose->torsion( tid );
		MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
		tors_mc->init();
		tors_mc->set_angle_range( init_torsion - angle_range_bb_ , init_torsion + angle_range_bb_ );
		tors_sampler.push_back( tors_mc );
		}
		tid = core::id::TorsionID( ii, id::BB, 2 );
		if ( mm_->get( tid ) ) {
		Real init_torsion = ref_pose->torsion( tid );
		MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
		tors_mc->init();
		tors_mc->set_angle_range( init_torsion - angle_range_bb_ , init_torsion + angle_range_bb_ );
		tors_sampler.push_back( tors_mc );
		}
		tid = core::id::TorsionID( ii, id::BB, 3 );
		if ( mm_->get( tid ) ) {
		Real init_torsion = ref_pose->torsion( tid );
		MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
		tors_mc->init();
		tors_mc->set_angle_range( init_torsion - angle_range_bb_ , init_torsion + angle_range_bb_ );
		tors_sampler.push_back( tors_mc );
		}
		}
		tid = core::id::TorsionID( ii, id::BB, 4 );
		if ( mm_->get( tid ) ) {
		Real init_torsion = ref_pose->torsion( tid );
		MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
		tors_mc->init();
		tors_mc->set_angle_range( init_torsion - angle_range_bb_ , init_torsion + angle_range_bb_ );
		tors_sampler.push_back( tors_mc );
		}
		tid = core::id::TorsionID( ii, id::BB, 5 );
		if ( mm_->get( tid ) ) {
		Real init_torsion = ref_pose->torsion( tid );
		MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
		tors_mc->init();
		tors_mc->set_angle_range( init_torsion - angle_range_bb_ , init_torsion + angle_range_bb_ );
		tors_sampler.push_back( tors_mc );
		}
		tid = core::id::TorsionID( ii, id::BB, 6 );
		if ( mm_->get( tid ) ) {
		Real init_torsion = ref_pose->torsion( tid );
		MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
		tors_mc->init();
		tors_mc->set_angle_range( init_torsion - angle_range_bb_ , init_torsion + angle_range_bb_ );
		tors_sampler.push_back( tors_mc );
		}
		*/
	}


	for ( Size i = 1; i <= residues_sampled.size(); ++i ) {
		if ( sampled_by_kic.find( residues_sampled[i]-1 ) == sampled_by_kic.end() &&residues_sampled[i] != 1 && ( i == 1 || residues_sampled[i] != residues_sampled[i-1]+1 ) ) {
			//create samplers for [i]-1 and [i]

			// Can't just have a continue condition because we might want to set up the second sampler
			//TR << "Going to attempt to sample " << residues_sampled[i] - 1 << " with a MC_RNA_Suite and pose size is " << pose.size() << std::endl;
			MC_RNA_SuiteOP suite_sampler_1( new MC_RNA_Suite( residues_sampled[i] - 1 ) );
			suite_sampler_1->set_init_from_pose( pose );
			suite_sampler_1->set_sample_bb( true );
			suite_sampler_1->set_sample_lower_nucleoside( false );
			suite_sampler_1->set_sample_upper_nucleoside( false );
			suite_sampler_1->set_pucker_flip_rate( 0 );
			suite_sampler_1->set_sample_near_a_form( false );
			suite_sampler_1->set_angle_range_from_init( 60 );
			standard_sampler.add_rotamer( suite_sampler_1 );
		}

		// can't do this!
		if ( residues_sampled[i] == pose.size() ) continue;
		if ( sampled_by_kic.find( residues_sampled[i] ) != sampled_by_kic.end() ) continue;

		//TR << "Going to attempt to sample " << residues_sampled[i] << " with a MC_RNA_Suite\tand pose size is " << pose.size() << std::endl;
		MC_RNA_SuiteOP suite_sampler( new MC_RNA_Suite( residues_sampled[i] ) );
		suite_sampler->set_init_from_pose( pose );
		suite_sampler->set_sample_bb( true );
		suite_sampler->set_sample_lower_nucleoside( false );
		suite_sampler->set_sample_upper_nucleoside( false );
		suite_sampler->set_pucker_flip_rate( 0 );
		suite_sampler->set_sample_near_a_form( false );
		suite_sampler->set_angle_range_from_init( 60 );
		standard_sampler.add_rotamer( suite_sampler );
	}
	//standard_sampler.init();

	//for ( Size i = 1; i <= residues_sampled.size(); ++i ) {
	//MC_RNA_SugarOP sampler( new MC_RNA_Sugar( residues_sampled[i], 1, 0 /*ANY_PUCKER*/ ) );
	// sugar_sampler.push_back( sampler );
	// }

	// Min-score pose
	Pose min_pose = pose;
	Real min_score( 99999 );

	auto st_temps = utility::tools::make_vector1< Real >( temp_ );
	auto st_weights = utility::tools::make_vector1< Real >( 1 );

	SimulatedTempering tempering( pose, score_fxn_, st_temps, st_weights );
	MonteCarlo mc(pose, *score_fxn_, 1.0 );
	Size n_accept_total( 0 ), n_accept_backbone( 0 ), n_accept_onetors( 0 ), n_accept_kic( 0 );
	Size /*n_total( 0 ), */n_backbone( 0 ), n_onetors( 0 ), n_kic( 0 );
	bool did_move;
	int index;
	//Size temp_id( 1 );

	set_gaussian_stdevs( kic_sampler, tors_sampler, standard_sampler, temp_, total_sampled );

	Pose stored_pose_ = pose;

	// Main sampling cycle
	for ( Size n = 1; n <= n_cycle_; ++n ) {
		did_move = true;
		bool boltz = false;
		// Sampler updates
		/*if ( (n % 100) == 0 ) {
		// Pick from among samplers
		index = numeric::random::rg().random_range(1,sugar_sampler.size());
		++sugar_sampler[index];
		sugar_sampler[index]->apply( pose );
		//if ( !(sugar_sampler[index]->check_moved() ) ) {
		// did_move = false;
		//}
		++n_sugar;
		if ( ( tempering.check_boltzmann( pose ) && did_move ) || n == n_cycle_ ) {
		stored_pose_ = pose;
		if ( n == n_cycle_ ) break;
		sugar_sampler[index]->update( pose ); // don't think this is necessary
		++n_accept_sugar;
		boltz = true;
		}
		} else if ( (n % 10) == 0 ) {
		// Pick from among samplers
		standard_sampler.set_angle( pose );
		++standard_sampler;
		standard_sampler.apply( pose );
		++n_backbone;
		if ( ( tempering.check_boltzmann( pose ) && did_move ) || n == n_cycle_ ) {
		stored_pose_ = pose;
		if ( n == n_cycle_ ) break;
		boltz = true;
		++n_accept_backbone;
		}
		} else */if ( (n % 2) == 0 && kic_sampler.size() > 0 ) {
			// Pick from among samplers
			index = numeric::random::rg().random_range(1,kic_sampler.size());
			// Let's check on the cutpoint status of the involved residues
			//TR << "About to sample cut between " << kic_sampler[index]->moving_suite() << " and " << kic_sampler[index]->chainbreak_suite() << std::endl;
			//TR << pose.fold_tree() << std::endl;
			kic_sampler[index]->next( pose ); //This function also updates the stored torsions
			kic_sampler[index]->apply( pose );
			if ( !(kic_sampler[index]->check_moved() ) ) {
				did_move = false;
			}
			++n_kic;
			if ( ( tempering.check_boltzmann( pose ) && did_move ) || n == n_cycle_ ) {
				stored_pose_ = pose;
				if ( n == n_cycle_ ) break;
				kic_sampler[index]->update(); // don't think this is necessary
				++n_accept_kic;
				boltz = true;
			}
		} else if ( tors_sampler.size() > 0 ) {
			// Pick from among samplers
			index = numeric::random::rg().random_range(1,tors_sampler.size());
			++(*tors_sampler[index]);
			tors_sampler[index]->apply( pose );
			++n_onetors;
			if ( ( tempering.check_boltzmann( pose ) && did_move ) || n == n_cycle_ ) {
				stored_pose_ = pose;
				if ( n == n_cycle_ ) break;
				tors_sampler[index]->update();
				// this logic/naming does not make sense -- here, n_accept_chi includes bb's.
				++n_accept_onetors;
				boltz = true;
			}
		}

		if ( ( boltz && did_move ) || n == n_cycle_ ) {
			++n_accept_total;
			if ( ( *score_fxn_ )( pose )< min_score ) {
				min_score = pose.energies().total_energy();
				min_pose = pose;
			}
		} else { // boltz = false: roll back the move.
			if ( did_move ) pose = stored_pose_;
		}
	}
	TR.Debug << "Score for the end pose" << std::endl;
	score_fxn_->show( TR.Debug, pose );
	TR.Debug << "Score for the min pose" << std::endl;
	score_fxn_->show( TR.Debug, min_pose );

	Real const end_score = pose.energies().total_energy();
	Real const minned_score = min_pose.energies().total_energy();

	TR << " start score " << start_score << " end " << end_score << " min " << minned_score << std::endl;

	TR << "n_cycles: " << n_cycle_ << std::endl;
	TR << "Total accept rate:       " << double( n_accept_total ) / n_cycle_ << std::endl;
	TR << "Suite accept rate:       " << double( n_accept_backbone ) / n_backbone << std::endl;
	TR << "One-torsion accept rate: " << double( n_accept_onetors ) / n_onetors << std::endl;
	TR << "KIC accept rate:         " << double( n_accept_kic ) / n_kic << std::endl;

	Real const time_in_test = static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC;
	TR << "Time in sampler: " <<  time_in_test << std::endl;

	if ( output_min_pose_ ) pose = min_pose;
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
ThermalMinimizer::parse_my_tag(
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
ThermalMinimizer::fresh_instance() const
{
	return protocols::moves::MoverOP( new ThermalMinimizer );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
ThermalMinimizer::clone() const
{
	return protocols::moves::MoverOP( new ThermalMinimizer( *this ) );
}

void
ThermalMinimizer::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, ThermalMinimizer const & mover )
{
	mover.show(os);
	return os;
}

/////////////// Creator ///////////////

std::string ThermalMinimizer::get_name() const {
	return mover_name();
}

std::string ThermalMinimizer::mover_name() {
	return "ThermalMinimizer";
}

void ThermalMinimizer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	//Doesn't parse anything!
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Use a simulated tempering simulation to refine a pose",
		attlist );
}

std::string ThermalMinimizerCreator::keyname() const {
	return ThermalMinimizer::mover_name();
}

protocols::moves::MoverOP
ThermalMinimizerCreator::create_mover() const {
	return protocols::moves::MoverOP( new ThermalMinimizer );
}

void ThermalMinimizerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ThermalMinimizer::provide_xml_schema( xsd );
}



} //scratch
} //recces
} //protocols

