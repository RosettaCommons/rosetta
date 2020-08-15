// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/normalmode/NormalModeRelaxMover.cc
/// @brief   Normal Mode Perturbation + Relax
/// @details
/// @author  Hahnbeom Park

#include <protocols/normalmode/NormalMode.hh>
#include <protocols/normalmode/NormalModeRelaxMover.hh>
#include <protocols/normalmode/NormalModeRelaxMoverCreator.hh>

//Option
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
//
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>

// Pose stuffs
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

//Relax, optimization
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/util.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

// Silent struct
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>

// Rosetta scripts
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility>
#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

//can't figure out how to feed a random number generator to the shuffle func
#include <random>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

static basic::Tracer TR( "protocols.normalmode.NormalModeRelaxMover" );

namespace protocols {
namespace normalmode {

using namespace ObjexxFCL::format;
using namespace core;

// Creator stuffs first


///////////////////////////

NormalModeRelaxMover::NormalModeRelaxMover()
{
	set_default();
}

NormalModeRelaxMover::NormalModeRelaxMover(
	core::scoring::ScoreFunctionCOP sfxn,
	bool const cartesian ):
	nsample_( 10 )
{
	set_default();
	cartesian_ = cartesian;

	mm_ = utility::pointer::make_shared< core::kinematics::MoveMap >();
	mm_->set_bb( true );
	mm_->set_chi( true );
	mm_->set_jump( true );

	NM_ = NormalMode( "CA", 10.0 );
	if ( !cartesian_ ) NM_.torsion( true );

	// Scorefunction
	sfxn_ = sfxn->clone();
}

// advanced
NormalModeRelaxMover::NormalModeRelaxMover(
	core::scoring::ScoreFunctionCOP sfxn,
	bool cartesian,
	core::kinematics::MoveMapCOP mm,
	std::string const relaxmode,
	core::Real const distcut
):
	nsample_( 10 )
{
	set_default();
	relaxmode_ = relaxmode;
	cartesian_ = cartesian;

	NM_ = NormalMode( "CA", distcut );
	if ( !cartesian_ ) NM_.torsion( true );

	if ( mm ) {
		mm_ = mm->clone();
	} else {
		mm_ = nullptr;
	}

	// Scorefunction
	sfxn_ = sfxn->clone();
}

NormalModeRelaxMover::~NormalModeRelaxMover()= default;

void
NormalModeRelaxMover::set_default()
{
	relaxmode_ = "relax";
	cartesian_ = true;
	minoption_ = utility::pointer::make_shared< optimization::MinimizerOptions >( "lbfgs_armijo_nonmonotone", 0.02, true, false, false );
	//minoption_->

	sfxn_ = scoring::ScoreFunctionFactory::create_score_function( "talaris2013_cart" );
	sfxn_cen_ = scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" );
}


void
NormalModeRelaxMover::apply( pose::Pose &pose )
{

	core::kinematics::MoveMapOP mm;
	if ( mm_ ) {
		mm = mm_;
	} else if ( mmf_ ) {
		mm = mmf_->create_movemap_from_pose( pose );
	} else {
		mm = utility::pointer::make_shared< core::kinematics::MoveMap >();
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_jump( true );
	}
	initialize_NM_for_movemap( pose, *mm );

	if ( centroid_ ) {
		protocols::moves::MoverOP tocen( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );
		tocen->apply( pose );
	}

	core::Size nstruct = mix_modes_ ? nsample_ : nmodes_*2;

	utility::vector1< core::Size > modes( nmodes_ );
	for ( core::Size i = 1; i <= nmodes_; ++i ) modes[i] = i;

	// setup modes first
	utility::vector1< utility::vector1< Real > > modescales( nstruct );
	for ( core::Size i_comb = 1; i_comb <= nstruct; ++i_comb ) {
		utility::vector1< core::Real > scalev( nmodes_, 0.0 );
		if ( mix_modes_ ) {
			for ( core::Size i = 1; i <= nmodes_; ++i ) scalev[i] = 2.0*numeric::random::rg().uniform() - 1.0;
		} else { // regular; -1 or +1
			core::Size i_mode = (i_comb+1)/2;
			if ( i_comb%2 == 0 ) {
				scalev[i_mode] = 1.0;
			} else {
				scalev[i_mode] = -1.0;
			}
		}
		modescales[i_comb] = scalev;
	}

	// apply
	for ( core::Size i_comb = 1; i_comb <= modescales.size(); ++i_comb ) {
		set_mode( modes, modescales[i_comb] );

		TR << "Normal mode for mode no " << i_comb;
		for ( core::Size i = 1; i <= nmodes_; ++i ) TR << " " << F(8,5,mode_scale_[i]);
		TR << std::endl;

		pose::Pose pose_tmp( pose );
		set_extrapolate_scale( pertscale_ );
		apply_on_pose( pose_tmp, mm );

		Real score;
		if ( centroid_ ) {
			runtime_assert( pose_tmp.is_centroid() );
			score = sfxn_cen_->score( pose_tmp );
		} else {
			runtime_assert( !pose_tmp.is_centroid() );
			score = sfxn_->score( pose_tmp );
		}
		TR << "score: " << F(8,3,score) << std::endl;

		output_pose_OPs_.push_back( utility::pointer::make_shared< pose::Pose >( pose_tmp ) );
	}


	if ( dump_silent_ ) {
		core::io::silent::SilentFileOptions opts;
		core::io::silent::SilentFileData sfd( opts );
		for ( core::Size i_pose = 1; i_pose <= output_pose_OPs_.size(); ++i_pose ) {
			core::io::silent::SilentStructOP ss =
				core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary", opts);
			std::stringstream tag;
			tag << "mode_" << i_pose;
			ss->fill_struct( *output_pose_OPs_[i_pose], tag.str() );
			sfd.write_silent_struct( *ss, outsilent_ );
		}
	}

	//if the user wants a random selection, shuffle it. Otherwise sort by descending score.
	//Using descending order because additional poses will be retrieved off the back of the
	//vector in get_additional_output and the best pose will be in the back (if sorted).
	if ( randomselect_ ) {
		//This is a lazy use of auto because I had a mildly hard time finding out
		//what type this should be. Forgive me.
		auto eng = std::default_random_engine {};
		std::shuffle(output_pose_OPs_.begin(), output_pose_OPs_.end(), eng);
	} else {
		//sorts pose OPs for descending score
		if ( centroid_ ) {
			//assertions for centroid-ness were done above
			std::sort(output_pose_OPs_.begin(), output_pose_OPs_.end(),
				[&](core::pose::PoseOP a, core::pose::PoseOP b) -> bool {
				return sfxn_cen_->score( *a ) > sfxn_cen_->score( *b );
				} );
		} else {
			std::sort(output_pose_OPs_.begin(), output_pose_OPs_.end(),
				[&](core::pose::PoseOP a, core::pose::PoseOP b) -> bool {
				return sfxn_->score( *a ) > sfxn_->score( *b );
				} );
		}
	}

	//this should produce the approprate result wrt randomselect_
	pose = *output_pose_OPs_.back();
	TR << "final pose with score " << F(8,3,(centroid_ ? sfxn_cen_->score(pose) : sfxn_->score(pose))) << std::endl;
	output_pose_OPs_.pop_back();
}

// Set single mode
void
NormalModeRelaxMover::set_mode( core::Size const i_mode )
{
	mode_using_.resize( 0 );
	mode_scale_.resize( 0 );
	mode_using_.push_back( i_mode );
	mode_scale_.push_back( 1.0 );
}

// Set combined mode
void
NormalModeRelaxMover::set_mode( utility::vector1< core::Size > const mode_using,
	utility::vector1< Real > const mode_scales )
{
	TR << "printing mode_scales" << std::endl << mode_scales << std::endl;
	mode_using_ = mode_using;
	mode_scale_ = mode_scales;

	// Normalize mode_scales
	Real scalesum( 0.0 );
	for ( core::Size i = 1; i <= mode_scales.size(); ++i ) {
		// Assert if values are valid
		debug_assert( mode_using[i] <= NM().nmode() );

		scalesum += mode_scales[i]*mode_scales[i];
	}

	for ( core::Size i = 1; i <= mode_scale_.size(); ++i ) { mode_scale_[i] = mode_scale_[i]/std::sqrt(scalesum); }
}

// Set randomized mode
void
NormalModeRelaxMover::set_random_mode(std::string const select_option,
	Real const importance_portion )
{
	mode_using_.resize( 0 );
	mode_scale_.resize( 0 );

	// Normalize mode_scales
	Real scalesum( 0.0 );
	for ( core::Size i = 1; i <= nmodes_; ++i ) {
		// Multiply scale by its importance
		Real scale1 = numeric::random::rg().uniform();
		if ( select_option == "probabilistic" ) {
			scale1 *= importance_portion*NM().get_importance( i ) + (1.0 - importance_portion);
		}

		mode_using_.push_back( i );
		mode_scale_.push_back( scale1 );
		scalesum += mode_scale_[i];
	}

	TR << "Mode setup: ";
	for ( core::Size i = 1; i <= mode_scale_.size(); ++i ) {
		mode_scale_[i] = mode_scale_[i]/scalesum;
		TR << " " << std::setw(8) << mode_scale_[i];
	}
	TR << std::endl;
}

void
NormalModeRelaxMover::initialize_NM_for_movemap(core::pose::Pose const & pose,
	core::kinematics::MoveMap const & movemap )
{
	// For Torsional NormalMode
	NM_.clear_torsions_using();
	for ( core::Size ires = 1; ires <= pose.size(); ++ires ) {
		if ( movemap.get_bb( ires ) ) {
			NM_.set_torsions_using( ires );
		}
	}
}

void
NormalModeRelaxMover::set_harmonic_constants( Real const k_uniform )
{
	NM_.set_harmonic_constants( k_uniform );
}

void
NormalModeRelaxMover::set_harmonic_constants( Real const k_connected,
	Real const k_segment,
	Real const k_long )
{
	NM_.set_harmonic_constants( k_connected, k_segment, k_long );
}

void
NormalModeRelaxMover::apply_on_pose( pose::Pose &pose, core::kinematics::MoveMapOP mm )
{
	pose::Pose const pose_init( pose );

	// define local minimization incase it's not differentiable
	core::scoring::ScoreFunctionOP sfxn_min = sfxn_->clone();
	core::scoring::ScoreFunctionOP sfxn_cen_min
		= scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" );

	sfxn_min->set_weight( core::scoring::coordinate_constraint, 1.0 );
	sfxn_cen_min->set_weight( core::scoring::coordinate_constraint, 1.0 );

	// NormalMode setup
	// Don't solve again NormalMode until "refresh_normalmode" is called
	if ( refresh_normalmode_ ) {
		TR << "Solving Normal Mode for given pose... " << std::endl;
		clock_t starttime = clock();
		NM_.solve( pose );
		clock_t endtime = clock();
		// important for bigger system to trace time limiting step
		TR << "NM solved in " << core::Real(endtime - starttime)/CLOCKS_PER_SEC << " sec." << std::endl;
		//refresh_normalmode_ = false;
	}

	// update only if torsional NM
	pose::Pose expose;
	if ( !cartesian_ ) expose = extrapolate_mode_on_pose( pose );
	utility::vector1< Vector > excrd = extrapolate_mode_on_crd( pose );
	gen_coord_constraint( pose, excrd );

	// Relax setup
	if ( relaxmode_ == "relax" ) {
		// This cannot run in centroid level
		runtime_assert( !pose.is_centroid() );

		protocols::relax::FastRelax relax_prot( sfxn_min, 1 );
		relax_prot.set_movemap( mm );
		relax_prot.min_type("lbfgs_armijo_nonmonotone");

		relax_prot.cartesian( cartesian_minimize_ );
		relax_prot.apply( pose );

	} else if ( relaxmode_ == "min" ) {
		//optimization::CartesianMinimizer minimizer;

		if ( !cartesian_ ) pose = expose; // start from perturbed

		core::scoring::ScoreFunctionOP sfxn_loc;
		// Pick proper scorefunction
		if ( pose.is_centroid() ) {
			sfxn_loc = sfxn_cen_min;
		} else {
			sfxn_loc = sfxn_min;
		}

		// Minimize
		if ( cartesian_minimize_ ) {
			optimization::CartesianMinimizer minimizer;
			minimizer.run( pose, *mm, *sfxn_loc, *minoption_ );

		} else {
			optimization::AtomTreeMinimizer minimizer;
			minimizer.run( pose, *mm, *sfxn_loc, *minoption_ );
		}

	} else if ( relaxmode_ == "extrapolate" ) {
		pose = expose;
	} else {
		TR.Warning << "unknown relaxmode " << relaxmode_ << ", doing nothing!" << std::endl;
	}

	// Measure RMSD to target
	Real rmsd = get_RMSD( excrd, pose );

	// Measure RMSD to init pose
	Real rmsd_init = core::scoring::CA_rmsd( pose_init, pose );

	TR << "Rmsd: initial vs relaxed: " << rmsd_init << std::endl;
	TR << "Rmsd: extrapolated vs relaxed: " << rmsd << std::endl;
}

///@detail normal mode relaxed poses are remembered as PoseOPs in a score-sorted fashion. Each time this is called,
///the best is removed from memory and returned until no more remain. Assumes randomselect_/order of
///output_pose_OPs_ isn't changed between calls because this would not be reflected in how addnl poses are retrived
///(sorting/shuffling is done only once)

core::pose::PoseOP
NormalModeRelaxMover::get_additional_output()
{
	core::pose::PoseOP out_pose_OP;

	if ( output_pose_OPs_.size() == 0 ) {
		return nullptr;
	}

	out_pose_OP = output_pose_OPs_.back();
	//pop_back normally calls destructors but we're storing PoseOPs so we should be alright
	output_pose_OPs_.pop_back();

	TR << "returning additional pose with score " << F(8,3,centroid_ ? sfxn_cen_->score(*out_pose_OP) : sfxn_->score(*out_pose_OP)) << std::endl;

	return out_pose_OP;

}

// Gen Coordinate Constraint
void
NormalModeRelaxMover::gen_coord_constraint( pose::Pose &pose,
	utility::vector1< Vector > const &excrd ) const
{
	// make sure there is no other constraint
	//pose.remove_constraints();

	core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, cst_sdev() ) );
	for ( core::Size i_atm = 1; i_atm <= NM().natm(); ++i_atm ) {
		if ( cartesian_ ) {
			pose.add_constraint(
				scoring::constraints::ConstraintCOP( utility::pointer::make_shared< core::scoring::constraints::CoordinateConstraint >
				( NM().get_atomID()[i_atm], NM().get_atomID()[i_atm], excrd[ i_atm ], fx ) ) );
		} else {
			core::Size resno( NM().get_atomID()[i_atm].rsd() );
			core::Size atmno( NM().get_atomID()[i_atm].atomno() );

			pose.add_constraint(
				scoring::constraints::ConstraintCOP( utility::pointer::make_shared< core::scoring::constraints::CoordinateConstraint >
				( NM().get_atomID()[i_atm], NM().get_atomID()[i_atm],
				pose.residue(resno).xyz(atmno), fx ) ) );
		}
	}
}

// Report rmsd to excrd
Real
NormalModeRelaxMover::get_RMSD( utility::vector1< Vector > const excrd,
	pose::Pose const &pose ) const
{
	Real rmsd( 0.0 );
	for ( core::Size ica = 1; ica <= excrd.size(); ++ica ) {
		core::Size resno( NM().get_atomID()[ica].rsd() );
		core::Size atmno( NM().get_atomID()[ica].atomno() );
		Vector const dxyz( excrd[ ica ] - pose.residue(resno).xyz(atmno) );
		rmsd += (dxyz[0]*dxyz[0]+dxyz[1]+dxyz[1]+dxyz[2]+dxyz[2]);
	}
	return std::sqrt(rmsd/excrd.size());
}


utility::vector1< Vector >
NormalModeRelaxMover::extrapolate_mode_on_crd( pose::Pose const &pose ) const
{
	utility::vector1< Vector > excrd( NM().natm(), Vector( 0.0 ) );
	utility::vector1< Vector > movevec( NM().natm(), Vector( 0.0 ) );
	//excrd.resize( NM().natm() );
	//movevec.resize( NM().natm() );

	// First, build moving vector by doing fusions of eigenvectors used
	for ( core::Size i_mode = 1; i_mode <= mode_using_.size(); ++i_mode ) {
		Real const &scale_i( mode_scale_[i_mode] );
		core::Size const &modeno( mode_using_[i_mode] );

		utility::vector1< Vector > const &eigvec( NM().get_eigvec_cart( modeno ) );

		for ( core::Size i_atm = 1; i_atm <= NM().natm(); ++i_atm ) {
			movevec[i_atm][0] += scale_i*eigvec[i_atm][0];
			movevec[i_atm][1] += scale_i*eigvec[i_atm][1];
			movevec[i_atm][2] += scale_i*eigvec[i_atm][2];
		}
	}

	// Then, measure rmsd for desired "normalized" moving vector
	Real rmsd( 0.0 );
	for ( core::Size ica = 1; ica <= movevec.size(); ++ica ) {
		Vector const &dxyz( movevec[ica] );
		rmsd += (dxyz[0]*dxyz[0]+dxyz[1]*dxyz[1]+dxyz[2]*dxyz[2]);
	}
	rmsd = std::sqrt(rmsd/movevec.size());

	// Multiply scale to excrd making desired movement as input
	scale_dynamic_ = moving_distance_/rmsd;

	TR << "CNM, Dynamic scale for set as " << scale_dynamic_ << ", given moving distance " << moving_distance_ << std::endl;

	for ( core::Size ica = 1; ica <= movevec.size(); ++ica ) {
		Vector const &dxyz( movevec[ica] );
		core::Size resno( NM().get_atomID()[ica].rsd() );
		core::Size atmno( NM().get_atomID()[ica].atomno() );
		Vector const &CAcrd( pose.residue( resno ).xyz( atmno ) );

		excrd[ ica ] =  CAcrd + direction_ * scale_dynamic_ * dxyz;
	}

	return excrd;
}

pose::Pose
NormalModeRelaxMover::extrapolate_mode_on_pose( pose::Pose const &pose ) const
{
	// Reset dtor
	dtor_.resize( NM().ntor() );
	for ( core::Size i_tor = 1; i_tor <= NM().ntor(); ++i_tor ) dtor_[i_tor] = 0.0;

	utility::vector1< id::TorsionID > const torIDs = NM().get_torID();

	// Get torsion perturbation from given mode setup
	for ( core::Size i_mode = 1; i_mode <= mode_using_.size(); ++i_mode ) {
		Real const &scale_i( mode_scale_[i_mode] );
		core::Size const &modeno( mode_using_[i_mode] );
		utility::vector1< Real > const &eigvec( NM().get_eigvec_tor( modeno ) );

		for ( core::Size i_tor = 1; i_tor <= NM().ntor(); ++i_tor ) {
			dtor_[i_tor] += direction_ * scale_i * eigvec[i_tor];
		}
	}

	// Find maxval
	Real maxval( 0.0 );
	for ( core::Size i_tor = 1; i_tor <= dtor_.size(); ++i_tor ) {
		if ( std::abs(dtor_[i_tor]) > maxval ) maxval = std::abs(dtor_[i_tor]);
	}

	// Trial scale - the scale making max torsion change as 5.0 degree
	Real trial_scale = 5.0/maxval;
	pose::Pose pose_trial( pose );
	for ( core::Size i_tor = 1; i_tor <= NM().ntor(); ++i_tor ) {
		Real tor = pose.torsion( torIDs[i_tor] ) + trial_scale*dtor_[i_tor];
		TR.Debug << "itor/rsd/dtor " << i_tor << " " << torIDs[i_tor].rsd();
		TR.Debug << " " << trial_scale*dtor_[i_tor] << std::endl;
		pose_trial.set_torsion( torIDs[i_tor], tor );
	}
	Real rmsd_trial = core::scoring::CA_rmsd( pose, pose_trial );

	// Adjust the scale of angles so as to bring RMSD to desired moving_distance_
	scale_dynamic_ = 5.0 * moving_distance_/rmsd_trial;
	Real const max_scale( 99.0 );
	if ( scale_dynamic_ > max_scale ) scale_dynamic_ = max_scale;
	if ( scale_dynamic_ < -max_scale ) scale_dynamic_ = -max_scale;

	// Apply to expose
	pose::Pose expose( pose );
	for ( core::Size i_tor = 1; i_tor <= NM().ntor(); ++i_tor ) {
		Real tor = pose.torsion( torIDs[i_tor] ) + scale_dynamic_/maxval*dtor_[i_tor];
		expose.set_torsion( torIDs[i_tor], tor );
	}

	core::scoring::CA_rmsd( pose, expose );

	TR << "TNM, Dynamic scale for set as " << scale_dynamic_ << ", given moving distance " << moving_distance_ << std::endl;
	return expose;
}

// RosettaScripts stuffs
// currently, its simply returning randomized structure
void NormalModeRelaxMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & data
) {
	cartesian_ = tag->getOption< bool >( "cartesian", true );
	centroid_ = tag->getOption< bool >( "centroid", false );
	nmodes_       = tag->getOption< core::Size >( "nmodes", 5 );
	mix_modes_    = tag->getOption< bool >( "mix_modes", false );
	pertscale_    = tag->getOption< Real >( "pertscale", 1.0 );
	randomselect_ = tag->getOption< bool >( "randomselect", false );
	relaxmode_    = tag->getOption< std::string >( "relaxmode", "min" );
	selection_kT_ = tag->getOption< Real >( "selection_kT", 1e6 ); // currently just placeholder
	cartesian_minimize_ = tag->getOption< bool >( "cartesian_minimize", false );
	nsample_     = tag->getOption< core::Size >( "nsample", nmodes_*2 );

	bool weighted_k( false );
	Real k_hbond( 1.0 ), k_long( 1.0 ), k_short( 1.0 );
	if ( tag->hasOption( "k_hbond" ) ) {
		k_hbond = tag->getOption< Real >( "k_hbond" );
		weighted_k = true;
	}
	if ( tag->hasOption( "k_short" ) ) {
		k_short = tag->getOption< Real >( "k_short" );
		weighted_k = true;
	}
	if ( tag->hasOption( "k_long" ) ) {
		k_long = tag->getOption< Real >( "k_long" );
		weighted_k = true;
	}

	if ( tag->hasOption( "outsilent" ) ) {
		dump_silent_ = true;
		outsilent_ = tag->getOption< std::string >( "outsilent" );
	}

	if ( tag->hasOption( "scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption< std::string >( "scorefxn" ) );
		if ( centroid_ ) {
			sfxn_cen_ = data.get_ptr<core::scoring::ScoreFunction>( "scorefxns", scorefxn_name );
		} else {
			sfxn_     = data.get_ptr<core::scoring::ScoreFunction>( "scorefxns", scorefxn_name );
		}
	}

	NM_ = NormalMode( "CA", 10.0 );
	if ( !cartesian_ ) NM_.torsion( true );
	if ( weighted_k ) set_harmonic_constants( k_short, k_hbond, k_long );

	core::select::movemap::MoveMapFactoryOP mmf( new core::select::movemap::MoveMapFactory );
	bool const chi( tag->getOption< bool >( "chi", true ) ), bb( tag->getOption< bool >( "bb", true ) );
	mmf->all_chi( chi );
	mmf->all_bb( bb );
	mmf_ = protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data, false, mmf );
}

std::string NormalModeRelaxMover::get_name() const {
	return mover_name();
}

std::string NormalModeRelaxMover::mover_name() {
	return "NormalModeRelax";
}

void NormalModeRelaxMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( [] (std::string const& name) {
		return mover_name() + "_subelement_" + name + "Type";
	});
	rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy(xsd, subelements);


	XMLSchemaRestriction relaxmode_enum;
	relaxmode_enum.name("relaxmode_name");
	relaxmode_enum.base_type(xs_string);
	relaxmode_enum.add_restriction( xsr_enumeration, "min" );
	relaxmode_enum.add_restriction( xsr_enumeration, "relax" );
	relaxmode_enum.add_restriction( xsr_enumeration, "extrapolate" );
	xsd.add_top_level_element(relaxmode_enum);

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"cartesian", xsct_rosetta_bool,
		"Use Cartesian normal model, which is recommended over Torsional normal mode",
		"true");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"centroid", xsct_rosetta_bool,
		"Use centroid description in the Mover. relaxmode = \"relax\" is prevented in this case",
		"false");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"nmodes", xsct_non_negative_integer,
		"How many normal modes to try",
		"5");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"mix_modes", xsct_rosetta_bool,
		"Whether use mixture of modes",
		"false");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"pertscale", xsct_real,
		"By how much to perturb the backbone in Angstrom",
		"1.0");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"randomselect", xsct_rosetta_bool,
		"Do random selection among sampled structures instead of picking best scoring one",
		"false");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"relaxmode", "relaxmode_name",
		"The way of relaxing structure. \"relax\" calls FastRelax, \"min\" Minimizer",
		"min");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"selection_kT", xsct_real,
		"XRW TO DO",
		"1e6");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"cartesian_minimize", xsct_rosetta_bool,
		"Use cartesian minimization for relax or minimization",
		"false");
	attlist + XMLSchemaAttribute(
		"nsample", xsct_non_negative_integer,
		"How many structures to try. This only works when mix_modes = true; otherwise will try 2*nmodes");
	attlist + XMLSchemaAttribute(
		"k_bond", xsct_real,
		"XRW TO DO");
	attlist + XMLSchemaAttribute(
		"k_short", xsct_real,
		"XRW TO DO");
	attlist + XMLSchemaAttribute(
		"k_long", xsct_real,
		"XRW TO DO");
	attlist + XMLSchemaAttribute(
		"outsilent", xs_string,
		"Specifies name of output silent file to store all the sampled structures. "
		"By default the Mover won't output any silent");
	attlist + XMLSchemaAttribute(
		"scorefxn", xs_string,
		"Scorefxn to relax poses and select among multiple solutions. For centroid mode, stage4");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"chi", xsct_rosetta_bool,
		"XRW TO DO",
		"true");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"bb", xsct_rosetta_bool,
		"XRW TO DO",
		"true");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(),
		"NormalModeRelax relaxes a structures to the coordinate directions derived from Anisotropic "
		"Network Model (ANM). The way how it works is: a) generate a extrapolated coordinates, "
		"b) put coordinate restraints to that coordinates, and c) run minimization or FastRelax. "
		""
		"Current implementation tries multiple normal and returns `nsamples`. Use a `MultiplePoseMover` to handle all the output samples."
		"If a `MultiplePoseMover` is not used, only the top scoring pose is returned (unless `randomselect` is True).",
		attlist, subelements );
}

std::string NormalModeRelaxMoverCreator::keyname() const {
	return NormalModeRelaxMover::mover_name();
}

protocols::moves::MoverOP
NormalModeRelaxMoverCreator::create_mover() const {
	return utility::pointer::make_shared< NormalModeRelaxMover >();
}

void NormalModeRelaxMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	NormalModeRelaxMover::provide_xml_schema( xsd );
}


} // namespace normalmode
} // namespace protocols
