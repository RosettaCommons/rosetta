// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/normalmode/NormalModeRelaxMover.cc
/// @brief   Normal Mode Perturbation + Relax
/// @detailed
/// @author  Hahnbeom Park

#include <protocols/normalmode/NormalMode.hh>
#include <protocols/normalmode/NormalModeRelaxMover.hh>

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
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

//Relax, optimization
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/util.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>

//Temporary
//#include <core/io/pdb/pose_io.hh>
//#include <sys/time.h>

static basic::Tracer TR("protocols.normalmode.NormalModeRelaxMover");

namespace protocols{
namespace normalmode{

static numeric::random::RandomGenerator RG( 151234 ); //Magic number??

void
NormalModeRelaxMover::set_harmonic_constants( Real const k_uniform )
{
  NM_.set_harmonic_constants( k_uniform );
}

// Set modes from NM
void
NormalModeRelaxMover::set_mode( utility::vector1< Size > const mode_using,
				utility::vector1< Real > const mode_scales )
{
  mode_using_ = mode_using;
	mode_scale_ = mode_scales;

  // Normalize mode_scales
  Real scalesum( 0.0 );
  for( Size i = 1; i <= mode_scales.size(); ++i ){
    // Assert if values are valid
    assert( mode_scales[i] > 0.0 );
    assert( mode_using[i] <= NM().nmode() );

    scalesum += mode_scales[i];
  }

  for( Size i = 1; i <= mode_scale_.size(); ++i ){ mode_scale_[i] = mode_scale_[i]/scalesum; }
}

void
NormalModeRelaxMover::set_movemap( core::pose::Pose const & pose,
																	 core::kinematics::MoveMapCOP movemap )
{
	mm_ = movemap->clone();

	// For Torsional NormalMode
	NM_.clear_torsions_using();
	for( Size ires = 1; ires <= pose.total_residue(); ++ires ){
		if( mm_->get_bb( ires ) )
			NM_.set_torsions_using( ires );
	}
}

void
NormalModeRelaxMover::set_mode( Size const i_mode )
{
  mode_using_.resize( 0 );
	mode_scale_.resize( 0 );
	mode_using_.push_back( i_mode );
	mode_scale_.push_back( 1.0 );
}

// Set modes randomly
void
NormalModeRelaxMover::set_random_mode( Size const nmode,
																			 std::string const select_option,
																			 Real const importance_portion )
{
  mode_using_.resize( 0 );
	mode_scale_.resize( 0 );

  // Normalize mode_scales
  Real scalesum( 0.0 );
  for( Size i = 1; i <= nmode; ++i ){
		// Multiply scale by its importance
		Real scale1 = RG.uniform();
		if( select_option.compare("probabilistic") ==  0 ){
			scale1 *= importance_portion*NM().get_importance( i ) + (1.0 - importance_portion);
		}

		mode_using_.push_back( i );
		mode_scale_.push_back( scale1 );
    scalesum += mode_scale_[i];
  }

	TR << "Mode setup: ";
  for( Size i = 1; i <= mode_scale_.size(); ++i ){
		mode_scale_[i] = mode_scale_[i]/scalesum;
		TR << " " << std::setw(8) << mode_scale_[i];
	}
	TR << std::endl;
}

void
NormalModeRelaxMover::set_harmonic_constants( Real const k_connected,
					      Real const k_segment,
					      Real const k_long )
{
  NM_.set_harmonic_constants( k_connected, k_segment, k_long );
}

///////////
//Cartesian
///////////
// Simple constructor
CartesianNormalModeMover::CartesianNormalModeMover(
		 core::pose::Pose const &, //pose
		 core::scoring::ScoreFunctionCOP sfxn,
		 std::string const relaxmode )
{
  NM_ = NormalMode( "CA", 10.0 );
  moving_distance_ = 1.0; // Move by 1.0 Angstrom
  refresh_normalmode_ = true;
	direction_ = 1.0;
	cst_sdev_= 1.0;
	relaxmode_ = relaxmode;
	cartesian_minimize_ = true;

	//set_movemap( pose, mm );
	set_default_minoption( );

	// Scorefunction
  sfxn_ = sfxn->clone();
	sfxn_cen_ = scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" );

	sfxn_->set_weight( core::scoring::coordinate_constraint, 1.0 );
	sfxn_cen_->set_weight( core::scoring::coordinate_constraint, 1.0 );

  // Default: Use mode 1 with 100% contribution
  utility::vector1< Size > mode1( 1, 1 );
  utility::vector1< Real > scale1( 1, 1.0 );
  set_mode( mode1, scale1 );
}

// Advanced constructor
CartesianNormalModeMover::CartesianNormalModeMover(
     core::pose::Pose const & pose,
		 core::scoring::ScoreFunctionCOP sfxn,
		 core::kinematics::MoveMapCOP mm,
		 std::string const mode,
		 Real const distcut,
		 std::string const relaxmode )
{
  NM_ = NormalMode( mode, distcut );
  moving_distance_ = 1.0; // Move by 1.0 Angstrom
  refresh_normalmode_ = true;
	direction_ = 1.0;
	cst_sdev_= 1.0;
	relaxmode_ = relaxmode;
	cartesian_minimize_ = true;

	set_movemap( pose, mm );
	set_default_minoption( );

	// Scorefunction
  sfxn_ = sfxn->clone();
	sfxn_cen_ = scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" );

	sfxn_->set_weight( core::scoring::coordinate_constraint, 1.0 );
	sfxn_cen_->set_weight( core::scoring::coordinate_constraint, 1.0 );

  // Default: Use mode 1 with 100% contribution
  utility::vector1< Size > mode1( 1, 1 );
  utility::vector1< Real > scale1( 1, 1.0 );
  set_mode( mode1, scale1 );
}

void
CartesianNormalModeMover::set_default_minoption(){
	minoption_ =
		new optimization::MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.02, true, false, false );
}

CartesianNormalModeMover::~CartesianNormalModeMover(){}

void
CartesianNormalModeMover::apply( pose::Pose &pose )
{
	pose::Pose const pose_init( pose );

  // NormalMode setup
  // Don't solve again NormalMode until "refresh_normalmode" is called
  if( refresh_normalmode_ ){
		TR << "Solving Normal Mode for given pose... " << std::endl;
    NM_.solve( pose );
    refresh_normalmode_ = false;
  }

	// Set extrapolation coordinate/ Coordinate constraint
	utility::vector1< Vector > excrd = extrapolate_mode( pose );
	gen_coord_constraint( pose, excrd );

	// Relax setup
	if( relaxmode_.compare("relax") == 0 ){
		// This cannot run in centroid level
		runtime_assert( !pose.is_centroid() );

		protocols::relax::FastRelax relax_prot( sfxn_, 1 );
		relax_prot.set_movemap( mm_ );
		relax_prot.min_type("lbfgs_armijo_nonmonotone");

		relax_prot.cartesian( true ); // force to use cartesian

		relax_prot.apply( pose );

	} else if ( relaxmode_.compare("min") == 0 ){
		optimization::CartesianMinimizer minimizer;

		core::scoring::ScoreFunctionOP sfxn_loc;
		// Pick proper scorefunction
		if( pose.is_centroid() ){
			sfxn_loc = sfxn_cen_->clone();
		} else {
			sfxn_loc = sfxn_->clone();
		}

		// Minimize
		if( cartesian_minimize_ ){
			optimization::CartesianMinimizer minimizer;
			minimizer.run( pose, *mm_, *sfxn_loc, *minoption_ );

		} else {
			optimization::AtomTreeMinimizer minimizer;
			minimizer.run( pose, *mm_, *sfxn_loc, *minoption_ );

		}
	}

	// Measure RMSD to target
	Real rmsd = get_RMSD( excrd, pose );

	// Measure RMSD to init pose
	Real rmsd_init = core::scoring::CA_rmsd( pose_init, pose );

  TR << "Rmsd: initial vs relaxed: " << rmsd_init << std::endl;
  TR << "Rmsd: extrapolated vs relaxed: " << rmsd << std::endl;
}

utility::vector1< Vector >
CartesianNormalModeMover::extrapolate_mode( pose::Pose const &pose )
{
  utility::vector1< Vector > excrd;
  utility::vector1< Vector > movevec;
  excrd.resize( NM().natm() );
  movevec.resize( NM().natm() );

  // First, build moving vector by doing fusions of eigenvectors used
  for( Size i_mode = 1; i_mode <= mode_using_.size(); ++i_mode ){
    Real const &scale_i( mode_scale_[i_mode] );
    Size const &modeno( mode_using_[i_mode] );

    utility::vector1< Vector > const &eigvec( NM().get_eigvec_cart( modeno ) );

    for( Size i_atm = 1; i_atm <= NM().natm(); ++i_atm ){
      movevec[i_atm] += scale_i*eigvec[i_atm];
    }
  }

  // Then, measure rmsd for desired "normalized" moving vector
  Real rmsd( 0.0 );
  for( Size ica = 1; ica <= movevec.size(); ++ica ){
    Vector const &dxyz( movevec[ica] );
    rmsd += (dxyz[0]*dxyz[0]+dxyz[1]*dxyz[1]+dxyz[2]*dxyz[2]);
  }
  rmsd = std::sqrt(rmsd/movevec.size());

  // Multiply scale to excrd making desired movement as input
  scale_dynamic_ = moving_distance_/rmsd;

	TR << "CNM, Dynamic scale for set as " << scale_dynamic_ << ", given moving distance " << moving_distance_ << std::endl;

  for( Size ica = 1; ica <= movevec.size(); ++ica ){
    Vector const &dxyz( movevec[ica] );
    Size resno( NM().get_atomID()[ica].rsd() );
    Size atmno( NM().get_atomID()[ica].atomno() );
    Vector const &CAcrd( pose.residue( resno ).xyz( atmno ) );

    excrd[ ica ] =  CAcrd + direction_ * scale_dynamic_ * dxyz;
  }

  return excrd;
}

// Gen Coordinate Constraint
void
CartesianNormalModeMover::gen_coord_constraint( pose::Pose &pose,
				      utility::vector1< Vector > const &excrd )
{
	// make sure there is no other constraint
	//pose.remove_constraints();

  for( Size i_atm = 1; i_atm <= NM().natm(); ++i_atm ){
    pose.add_constraint(
		new core::scoring::constraints::CoordinateConstraint
	  ( NM().get_atomID()[i_atm], NM().get_atomID()[i_atm], excrd[ i_atm ],
			new core::scoring::constraints::HarmonicFunc( 0.0, cst_sdev() ) ) );
  }
}

// Report rmsd to excrd
Real
CartesianNormalModeMover::get_RMSD( utility::vector1< Vector > const excrd,
				pose::Pose const &pose ){

  Real rmsd( 0.0 );
  for( Size ica = 1; ica <= excrd.size(); ++ica ){
    Size resno( NM().get_atomID()[ica].rsd() );
    Size atmno( NM().get_atomID()[ica].atomno() );
    Vector const dxyz( excrd[ ica ] - pose.residue(resno).xyz(atmno) );
    rmsd += (dxyz[0]*dxyz[0]+dxyz[1]+dxyz[1]+dxyz[2]+dxyz[2]);
  }
  return std::sqrt(rmsd/excrd.size());
}

///////////
//Torsion
///////////
// Simple constructor
TorsionNormalModeMover::TorsionNormalModeMover(
		 core::pose::Pose const &, //pose
		 core::scoring::ScoreFunctionCOP sfxn,
		 std::string const relaxmode )

{
  NM_ = NormalMode( "CA", 10.0 );
	NM_.torsion( true );

  moving_distance_ = 1.0; // Move by 1.0 Angstrom
  refresh_normalmode_ = true;
	direction_ = 1.0;
	cst_sdev_ = 1.0; // in Angstrom
	relaxmode_ = relaxmode;
	cartesian_minimize_ = false;

	//set_movemap( pose, mm );
	set_default_minoption( );

	// Scorefunction
  sfxn_ = sfxn->clone();
	sfxn_cen_ = scoring::ScoreFunctionFactory::create_score_function( "score4_smooth" );

	sfxn_->set_weight( core::scoring::coordinate_constraint, 1.0 );
	sfxn_cen_->set_weight( core::scoring::coordinate_constraint, 1.0 );

  // Default: Use mode 1 with 100% contribution
  utility::vector1< Size > mode1( 1, 1 );
  utility::vector1< Real > scale1( 1, 1.0 );
  set_mode( mode1, scale1 );
}

// Advanced constructor
TorsionNormalModeMover::TorsionNormalModeMover(
     core::pose::Pose const & pose,
		 core::scoring::ScoreFunctionCOP sfxn,
		 core::kinematics::MoveMapCOP mm,
		 std::string const mode,
		 Real const distcut,
		 std::string const relaxmode )

{
  NM_ = NormalMode( mode, distcut );
	NM_.torsion( true );

  moving_distance_ = 1.0; // Move by 1.0 Angstrom
  refresh_normalmode_ = true;
	direction_ = 1.0;
	cst_sdev_ = 1.0; // in Angstrom
	relaxmode_ = relaxmode;
	cartesian_minimize_ = false;

	set_movemap( pose, mm );
	set_default_minoption( );

	// Scorefunction
  sfxn_ = sfxn->clone();
	sfxn_cen_ = scoring::ScoreFunctionFactory::create_score_function( "score4_smooth" );

	sfxn_->set_weight( core::scoring::coordinate_constraint, 1.0 );
	sfxn_cen_->set_weight( core::scoring::coordinate_constraint, 1.0 );

  // Default: Use mode 1 with 100% contribution
  utility::vector1< Size > mode1( 1, 1 );
  utility::vector1< Real > scale1( 1, 1.0 );
  set_mode( mode1, scale1 );
}

TorsionNormalModeMover::~TorsionNormalModeMover(){}

void
TorsionNormalModeMover::set_default_minoption(){
	//optimization::MinimizerOptionsCOP minoption =
	//	new ( "dfpmin", 0.02, true, false, false );
	//set_minoption( minoption );
	minoption_ = new optimization::MinimizerOptions( "dfpmin", 0.02, true, false, false );
}

void
TorsionNormalModeMover::apply( pose::Pose &pose )
{
	pose::Pose const pose_init( pose );
	protocols::moves::MoverOP tocen
		= new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID );
	protocols::moves::MoverOP tofa
		= new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD );

  // NormalMode setup
  // Don't solve again NormalMode until "refresh_normalmode" is called
  if( refresh_normalmode_ ){
		TR << "Solving Normal Mode for given pose... " << std::endl;
    NM_.solve( pose );
    refresh_normalmode_ = false;
  }

	// Get weight of perturbation given mode setup
	//if( mode_changed_ )	set_dynamic_scale( pose );
	//mode_changed_ = false; // Store information not to update dynamic scale again

	pose::Pose expose = extrapolate_mode( pose );
	utility::vector1< Vector > dummy;

	// Relax setup
	if( relaxmode_.compare("relax") == 0 ){
		// This cannot run in centroid level
		runtime_assert( !pose.is_centroid() );

		protocols::relax::FastRelax relax_prot( sfxn_, 1 );
		relax_prot.set_movemap( mm_ );
		relax_prot.min_type("lbfgs_armijo_nonmonotone");

		gen_coord_constraint( pose, dummy );

		relax_prot.apply( pose );

	} else if ( relaxmode_.compare("min") == 0 ){
		pose = expose; // Start from perturbed
		gen_coord_constraint( pose, dummy );

		core::scoring::ScoreFunctionOP sfxn_loc;
		// Pick proper scorefunction
		if( pose.is_centroid() ){
			sfxn_loc = sfxn_cen_->clone();
		} else {
			sfxn_loc = sfxn_->clone();
		}

		// Minimize
		if( cartesian_minimize_ ){
			optimization::CartesianMinimizer minimizer;
			minimizer.run( pose, *mm_, *sfxn_loc, *minoption_ );

		} else {
			optimization::AtomTreeMinimizer minimizer;
			minimizer.run( pose, *mm_, *sfxn_loc, *minoption_ );

		}

	} else if ( relaxmode_.compare("extrapolate") == 0 ){
		pose = expose;

	}

	// Measure RMSD to target
	Real rmsd = core::scoring::CA_rmsd( pose, expose );
	// Measure RMSD to init pose
	Real rmsd_init = core::scoring::CA_rmsd( pose_init, pose );

  TR << "Rmsd: initial vs relaxed: " << rmsd_init << std::endl;
  TR << "Rmsd: extrapolated vs relaxed: " << rmsd << std::endl;

}

// Gen Coordinate Constraint
void
TorsionNormalModeMover::gen_coord_constraint( pose::Pose &pose,
					utility::vector1< Vector > const & )
{
	// make sure there is no other constraint
	//pose.remove_constraints();

  for( Size i_atm = 1; i_atm <= NM().natm(); ++i_atm ){
    Size resno( NM().get_atomID()[i_atm].rsd() );
    Size atmno( NM().get_atomID()[i_atm].atomno() );

		// Just put current position
    pose.add_constraint(
		new core::scoring::constraints::CoordinateConstraint
	  ( NM().get_atomID()[i_atm], NM().get_atomID()[i_atm],
			pose.residue(resno).xyz(atmno),
			new core::scoring::constraints::HarmonicFunc( 0.0, cst_sdev() ) ) );
  }
}


pose::Pose
TorsionNormalModeMover::extrapolate_mode( pose::Pose const &pose )
{
	// Reset dtor
	dtor_.resize( NM().ntor() );
	for( Size i_tor = 1; i_tor <= NM().ntor(); ++i_tor ) dtor_[i_tor] = 0.0;

	utility::vector1< id::TorsionID > const torIDs = NM().get_torID();

	// Get torsion perturbation from given mode setup
	for( Size i_mode = 1; i_mode <= mode_using_.size(); ++i_mode ){
		Real const &scale_i( mode_scale_[i_mode] );
		Size const &modeno( mode_using_[i_mode] );
    utility::vector1< Real > const &eigvec( NM().get_eigvec_tor( modeno ) );

		for( Size i_tor = 1; i_tor <= NM().ntor(); ++i_tor ){
			dtor_[i_tor] += direction_ * scale_i * eigvec[i_tor];
		}
	}

	// Find maxval
	Real maxval( 0.0 );
	for( Size i_tor = 1; i_tor <= dtor_.size(); ++i_tor ){
		if( std::abs(dtor_[i_tor]) > maxval ) maxval = std::abs(dtor_[i_tor]);
	}

	// Trial scale - the scale making max torsion change as 5.0 degree
	Real trial_scale = 5.0/maxval;
	pose::Pose pose_trial( pose );
	for( Size i_tor = 1; i_tor <= NM().ntor(); ++i_tor ){
		Real tor = pose.torsion( torIDs[i_tor] ) + trial_scale*dtor_[i_tor];
		TR.Debug << "itor/rsd/dtor " << i_tor << " " << torIDs[i_tor].rsd();
		TR.Debug << " " << trial_scale*dtor_[i_tor] << std::endl;
		pose_trial.set_torsion( torIDs[i_tor], tor );
	}
	Real rmsd_trial = core::scoring::CA_rmsd( pose, pose_trial );

	// Adjust the scale of angles so as to bring RMSD to desired moving_distance_
	scale_dynamic_ = 5.0 * moving_distance_/rmsd_trial;
	Real const max_scale( 99.0 );
	if( scale_dynamic_ > max_scale ) scale_dynamic_ = max_scale;
	if( scale_dynamic_ < -max_scale ) scale_dynamic_ = -max_scale;

	// Apply to expose
	pose::Pose expose( pose );
	for( Size i_tor = 1; i_tor <= NM().ntor(); ++i_tor ){
		Real tor = pose.torsion( torIDs[i_tor] ) + scale_dynamic_/maxval*dtor_[i_tor];
		expose.set_torsion( torIDs[i_tor], tor );
	}

	Real rmsd_est = core::scoring::CA_rmsd( pose, expose );

	TR << "TNM, Dynamic scale for set as " << scale_dynamic_ << ", given moving distance " << moving_distance_ << std::endl;
	return expose;
}

} // namespace normalmode
} // namespace protocols
