// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farna/thermal_sampling/ThermalMinimizer.cc
/// @brief Use a simulated tempering simulation to refine a pose
/// @author Andy Watkins (amw579@nyu.edu)

// Unit headers
#include <protocols/stepwise/modeler/ThermalMinimizer.hh>
#include <protocols/stepwise/modeler/ThermalMinimizerCreator.hh>

// Core headers
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
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

#include <protocols/stepwise/sampler/rna/RNA_MC_Suite.hh>
#include <protocols/stepwise/sampler/rna/RNA_MC_MultiSuite.hh>
#include <protocols/moves/SimulatedTempering.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/stepwise/sampler/rna/RNA_MC_KIC_Sampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_KIC_Sampler.hh>
#include <protocols/farna/thermal_sampling/util.hh>
#include <protocols/farna/thermal_sampling/thermal_sampler.hh>

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.ThermalMinimizer" );

namespace protocols {
namespace stepwise {
namespace modeler {


void set_gaussian_stdevs(
	utility::vector1<protocols::stepwise::sampler::rna::RNA_MC_KIC_SamplerOP> & internal_bb_sampler,
	utility::vector1<protocols::stepwise::sampler::MC_OneTorsionOP> & chi_sampler,
	Real const temp,
	Size const //total_sampled
) {
	Real internal_bb_stdev( 0.1 * pow( temp, 0.25 ) + 0.1);
	Real chi_stdev( 5 * pow( temp , 0.5) + 15 );
	if ( temp < 0 ) {
		internal_bb_stdev = 0.5 ;
		chi_stdev = -1 ;
	}
	for ( Size i = 1; i <= internal_bb_sampler.size(); ++i ) {
		internal_bb_sampler[i]->set_gaussian_stdev( internal_bb_stdev );
	}
	for ( Size i = 1; i <= chi_sampler.size(); ++i ) {
		chi_sampler[i]->set_gaussian_stdev( chi_stdev );
	}
}



/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ThermalMinimizer::ThermalMinimizer():
	protocols::moves::Mover( ThermalMinimizer::class_name() ),
	n_cycle_( basic::options::option[ basic::options::OptionKeys::rna::farna::thermal_sampling::n_cycle ] ),
	angle_range_chi_( basic::options::option[ basic::options::OptionKeys::rna::farna::thermal_sampling::angle_range_chi ]() ),
	angle_range_bb_( basic::options::option[ basic::options::OptionKeys::rna::farna::thermal_sampling::angle_range_bb ]() )
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

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::rna::farna::thermal_sampling;
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

	Size total_sampled = 0;
	if ( mm_ ) {
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
			}
		}
	}

	if ( total_sampled == 0 ) return;

	clock_t const time_start( clock() );

	// do it
	score_fxn_->show( pose );

	//Size const total_len( pose.size() );
	//bool const sample_near_a_form( false );

	//Setting up the internal move sampler
	pose::PoseOP ref_pose(new pose::Pose(pose));

	utility::vector1<RNA_MC_KIC_SamplerOP> sampler;
	utility::vector1<MC_OneTorsionOP> tors_sampler;
	//Set up the internal move samplers
	for ( Size ii = 2; ii <= pose.size() - 1; ++ii ) {
		// use the movemap, luke
		if ( !mm_ || !mm_compatible_with_kic( mm_, ii ) ) continue;

		RNA_MC_KIC_SamplerOP suite_sampler( new RNA_MC_KIC_Sampler( ref_pose, ii-1, ii ) );
		suite_sampler->init();
		suite_sampler->set_angle_range_from_init_torsions( angle_range_bb_ );
		sampler.push_back( suite_sampler );
	}

	//Set up the torsion samplers
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		// Don't set up OneTorsion samplers for which KIC could do the job.
		if ( !mm_ || mm_compatible_with_kic( mm_, ii ) ) continue;

		// skip if movemap exists but dof off
		auto tid = core::id::TorsionID( ii, id::BB, 1 );
		if ( mm_ && mm_->get( tid ) ) {
			Real init_torsion = ref_pose->torsion( tid );
			MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
			tors_mc->init();
			tors_mc->set_angle_range( init_torsion - angle_range_bb_ , init_torsion + angle_range_bb_ );
			tors_sampler.push_back( tors_mc );
		}
		tid = core::id::TorsionID( ii, id::BB, 2 );
		if ( mm_ && mm_->get( tid ) ) {
			Real init_torsion = ref_pose->torsion( tid );
			MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
			tors_mc->init();
			tors_mc->set_angle_range( init_torsion - angle_range_bb_ , init_torsion + angle_range_bb_ );
			tors_sampler.push_back( tors_mc );
		}
		tid = core::id::TorsionID( ii, id::BB, 3 );
		if ( mm_ && mm_->get( tid ) ) {
			Real init_torsion = ref_pose->torsion( tid );
			MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
			tors_mc->init();
			tors_mc->set_angle_range( init_torsion - angle_range_bb_ , init_torsion + angle_range_bb_ );
			tors_sampler.push_back( tors_mc );
		}
		tid = core::id::TorsionID( ii, id::BB, 4 );
		if ( mm_ && mm_->get( tid ) ) {
			Real init_torsion = ref_pose->torsion( tid );
			MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
			tors_mc->init();
			tors_mc->set_angle_range( init_torsion - angle_range_bb_ , init_torsion + angle_range_bb_ );
			tors_sampler.push_back( tors_mc );
		}
		tid = core::id::TorsionID( ii, id::BB, 5 );
		if ( mm_ && mm_->get( tid ) ) {
			Real init_torsion = ref_pose->torsion( tid );
			MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
			tors_mc->init();
			tors_mc->set_angle_range( init_torsion - angle_range_bb_ , init_torsion + angle_range_bb_ );
			tors_sampler.push_back( tors_mc );
		}
		tid = core::id::TorsionID( ii, id::BB, 6 );
		if ( mm_ && mm_->get( tid ) ) {
			Real init_torsion = ref_pose->torsion( tid );
			MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
			tors_mc->init();
			tors_mc->set_angle_range( init_torsion - angle_range_bb_ , init_torsion + angle_range_bb_ );
			tors_sampler.push_back( tors_mc );
		}
		tid = core::id::TorsionID( ii, id::CHI, 1 );
		if ( mm_ && mm_->get( tid ) ) {
			Real init_torsion = ref_pose->torsion( tid );
			MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
			tors_mc->init();
			tors_mc->set_angle_range( init_torsion - angle_range_chi_ , init_torsion + angle_range_chi_ );
			tors_sampler.push_back( tors_mc );
		}
		tid = core::id::TorsionID( ii, id::CHI, 4 );
		if ( mm_ && mm_->get( tid ) ) {
			Real init_torsion = ref_pose->torsion( tid );
			MC_OneTorsionOP tors_mc( new MC_OneTorsion( tid, init_torsion ) );
			tors_mc->init();
			tors_mc->set_angle_range( init_torsion - angle_range_chi_ , init_torsion + angle_range_chi_ );
			tors_sampler.push_back( tors_mc );
		}
	}

	// Min-score pose
	Pose min_pose = pose;
	Real min_score( 99999 );

	auto st_temps = utility::tools::make_vector1< Real >( temp_ );
	auto st_weights = utility::tools::make_vector1< Real >( 0 );

	SimulatedTempering tempering( pose, score_fxn_, st_temps, st_weights );
	MonteCarlo mc(pose, *score_fxn_, 1.0 );
	Size n_accept_total( 0 ), n_accept_backbone( 0 );
	Size n_accept_chi( 0 );
	bool did_move;
	bool boltz = false;
	int index;
	//Size temp_id( 1 );

	set_gaussian_stdevs( sampler, tors_sampler, temp_, total_sampled );

	Pose stored_pose_ = pose;

	// Main sampling cycle
	for ( Size n = 1; n <= n_cycle_; ++n ) {
		did_move = true;
		// Sampler updates
		if ( (n % 2) == 0 && sampler.size() > 0 ) {
			// Pick from among samplers
			index = numeric::random::rg().random_range(1,sampler.size());
			sampler[index]->next( pose ); //This function also updates the stored torsions
			sampler[index]->apply( pose );
			if ( !(sampler[index]->check_moved() ) ) {
				did_move = false;
			}
			if ( ( tempering.boltzmann( pose ) && did_move ) || n == n_cycle_ ) {
				stored_pose_ = pose;
				if ( n == n_cycle_ ) break;
				sampler[index]->update( pose ); // don't think this is necessary
				++n_accept_backbone;
				boltz = true;
			}
		} else if ( tors_sampler.size() > 0 ) {
			// Pick from among samplers
			index = numeric::random::rg().random_range(1,tors_sampler.size());
			++(*tors_sampler[index]);
			tors_sampler[index]->apply( pose );
			if ( ( tempering.boltzmann( pose ) && did_move ) || n == n_cycle_ ) {
				stored_pose_ = pose;
				if ( n == n_cycle_ ) break;
				tors_sampler[index]->update();
				++n_accept_chi;
				boltz = true;
			}
		}

		if ( boltz ) {
			++n_accept_total;
			if ( ( *score_fxn_ )( pose )< min_score ) {
				min_score = pose.energies().total_energy();
				min_pose = pose;
			}
		} else {
			if ( did_move ) pose = stored_pose_;
		}
	}

	TR << "Score for the end pose" << std::endl;
	score_fxn_->show( pose );
	TR << "Score for the min pose" << std::endl;
	score_fxn_->show( min_pose );

	TR << "n_cycles: " << n_cycle_ << std::endl;
	TR << "Total accept rate: " << double( n_accept_total ) / n_cycle_ << std::endl;
	TR << "Backbone accept rate: " << double( n_accept_backbone ) / ((n_cycle_ / 2) - (n_cycle_ / 10)) << std::endl;
	TR << "Chi accept rate: " << double( n_accept_chi ) / (n_cycle_ / 2) << std::endl;

	Real const time_in_test = static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC;
	std::cout << "Time in sampler: " <<  time_in_test << std::endl;

	pose = min_pose;
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

std::string
ThermalMinimizer::get_name() const
{
	return ThermalMinimizer::class_name();
}

std::string
ThermalMinimizer::class_name()
{
	return "ThermalMinimizer";
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

protocols::moves::MoverOP
ThermalMinimizerCreator::create_mover() const
{
	return protocols::moves::MoverOP( new ThermalMinimizer );
}

std::string
ThermalMinimizerCreator::keyname() const
{
	return ThermalMinimizer::class_name();
}


} //modeler
} //stepwise
} //protocols

