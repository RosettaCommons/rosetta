// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file DockingInitialPerturbation.cc
/// @brief initial position functions
/// @detailed
///		This contains the functions that create initial positions for docking
/// @author Ingemar Andre

#include <protocols/simple_moves/symmetry/SymDockingInitialPerturbation.hh>

// Rosetta Headers
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/kinematics/Jump.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

// for symmetry
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <string>

//Utility Headers
#include <numeric/trig.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.simple_moves_symmetry.SymDockingInitialPerturbation");
static core::Size trans ( 1 ), rot ( 2 );

using namespace core;
using namespace conformation::symmetry;

namespace protocols {
namespace simple_moves{
namespace symmetry {

static numeric::random::RandomGenerator RG(4227034);

// Symmetric version of initial perturbation on one of the partners
// the type of perturbation is defined in the options
// some of the options are randomize1 or randomize2 (they are the same),
// dock_pert
//------------------------------------------------------------------------------
//
//     there are several ways to perturb the structure before beginning
//     the search; they are controlled through command-line flags
//
//     at the end, partners are slid into contact and scored
//
	// default constructor
	SymDockingInitialPerturbation::SymDockingInitialPerturbation() : protocols::moves::Mover()
	{
		slide_ = false;
		protocols::moves::Mover::type( "SymmDockingInitialPerturbation" );
	}

	// constructor with arguments
	SymDockingInitialPerturbation::SymDockingInitialPerturbation(
		bool const slide_in
	) : protocols::moves::Mover(),
			slide_(slide_in)
	{
		protocols::moves::Mover::type( "SymmDockingInitialPerturbation" );
	}

SymDockingInitialPerturbation::~SymDockingInitialPerturbation(){}

////////////////////////////////////////////////////////////////////////////////
/// @begin initial_perturbation
///
/// @brief   Make starting perturbations for rigid body moves
///
/////////////////////////////////////////////////////////////////////////////////
void SymDockingInitialPerturbation::apply( core::pose::Pose & pose )
{
	using namespace moves;
	using namespace basic::options;
	//////////////////////////////////
	//The options are -symmetry::initialize_rigid_body_dofs
	//								-symmetry::perturb_rigid_body_dofs <trans> <rot>
	//						for dock_pert, also need to get the normal perturbation
	// In the future this is specified in the symmetry definition file
	//////////////////////////////////
	assert( core::pose::symmetry::is_symmetric( pose ));
	SymmetricConformation & symm_conf (
        dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

	std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );
	SymSlideInfo const & slide_info( symm_conf.Symmetry_Info()->get_slide_info() );
	TR << "Reading options..." << std::endl;
	if( option[ OptionKeys::symmetry::initialize_rigid_body_dofs ]() ) {
		TR << "initialize_rigid_body_dofs: true" << std::endl;
		rigid::RigidBodyDofSeqRandomizeMover mover( dofs );
		mover.apply( pose );
	}

	if( option[ OptionKeys::symmetry::perturb_rigid_body_dofs ].user() ) {
		TR << "perturb_rigid_body_dofs: true" << std::endl;
		/// read in dock_pert options from commandline.  the first value is the
		/// rotation magnitude and the second value is the translational value
		utility::vector1< Real > pert_mags = option[ OptionKeys::symmetry::perturb_rigid_body_dofs ]();
		TR << "option[ symmetry::perturb_rigid_body_dofs ]() rot=" << pert_mags[rot] << "  trans=" << pert_mags[trans] << std::endl;
		rigid::RigidBodyDofSeqPerturbMover mover( dofs, pert_mags[rot], pert_mags[trans] );
		mover.apply( pose );
	}
	// DO NOT do this for e.g. ligand docking
	if ( slide_ ) {
		if ( slide_info.get_slide_type() == SEQUENTIAL ) {
			SequentialSymmetrySlider symm_slider = SequentialSymmetrySlider( pose,
																																			 slide_info.get_SlideCriteriaType(),
																																			 slide_info.get_SlideCriteriaVal() );
			symm_slider.apply( pose );
		}
		if ( slide_info.get_slide_type() == ORDERED_SEQUENTIAL ) {
			OrderedSequentialSymmetrySlider symm_slider = OrderedSequentialSymmetrySlider( pose,
																																										 slide_info.get_SlideCriteriaType(),
																																										 slide_info.get_SlideCriteriaVal(),
																																										 slide_info.get_slide_order() );
			symm_slider.apply( pose );
		}
		if ( slide_info.get_slide_type() == RANDOM ) {
			RandomSymmetrySlider symm_slider = RandomSymmetrySlider( pose,
																															 slide_info.get_SlideCriteriaType(),
																															 slide_info.get_SlideCriteriaVal() );
			symm_slider.apply( pose );
		}

//		SymDockingSlideIntoContact slide( dofs );
//		slide.apply( pose );
	}
}

std::string
SymDockingInitialPerturbation::get_name() const {
	return "SymDockingInitialPerturbation";
}


	// default constructor
	SymDockingSlideIntoContact::SymDockingSlideIntoContact() : protocols::moves::Mover() {	}

	// constructor with arguments
	SymDockingSlideIntoContact::SymDockingSlideIntoContact(
		std::map< Size, core::conformation::symmetry::SymDof > dofs
	) : protocols::moves::Mover(),
			dofs_(dofs)
	{

		protocols::moves::Mover::type( "SymDockingSlideIntoContact" );
		core::scoring::symmetry::SymmetricScoreFunction scorefxn_sym ( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::CENTROID_WTS, core::scoring::DOCK_LOW_PATCH ) );
		scorefxn_ = new core::scoring::symmetry::SymmetricScoreFunction( scorefxn_sym );
	}

	SymDockingSlideIntoContact::~SymDockingSlideIntoContact(){}

void SymDockingSlideIntoContact::apply( core::pose::Pose & pose )
{
	using namespace moves;

	assert( core::pose::symmetry::is_symmetric( pose ));
	SymmetricConformation & symm_conf (
        dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

	Size num_slide_moves(1);
	std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );
	rigid::RigidBodyDofRandomTransMover mover( dofs );
	( *scorefxn_ )( pose );
	TR.Debug << "score " << pose.energies().total_energies()[ scoring::interchain_vdw ]  << std::endl;
	TR.Debug << "sliding into contact" << std::endl;
	TR.Debug << "Moving away" << std::endl;
	// first try moving away from each other
	while ( pose.energies().total_energies()[ scoring::interchain_vdw ] > 0.1 ) {
		mover.apply( pose );
		( *scorefxn_ )( pose );
		if ( ++num_slide_moves > 1000 ) {
			std::cerr << "To many slide moves. Subunits never touching..." << std::endl;
			utility_exit();
		}
	TR.Debug << "score away " << pose.energies().total_energies()[ scoring::interchain_vdw ]  << std::endl;
	}
	// then try moving towards each other
	TR.Debug << "Moving together" << std::endl;
	mover.trans_axis().negate();
	while ( pose.energies().total_energies()[ scoring::interchain_vdw ] < 0.1 ) {
		mover.apply( pose );
		( *scorefxn_ )( pose );
		if ( ++num_slide_moves > 1000 ) {
			std::cerr << "To many slide moves. Subunits never touching..." << std::endl;
			utility_exit();
		}
		TR.Debug << "score together " << pose.energies().total_energies()[ scoring::interchain_vdw ]  << std::endl;
	}
	// move away again until just touching
	mover.trans_axis().negate();
	mover.apply( pose );

}

std::string
SymDockingSlideIntoContact::get_name() const {
	return "SymDockingSlideIntoContact";
}


	FaSymDockingSlideTogether::FaSymDockingSlideTogether(
		std::map< Size, core::conformation::symmetry::SymDof > dofs
	) : protocols::moves::Mover(),
			dofs_(dofs),
			tolerance_(0.2)
	{
		protocols::moves::Mover::type( "FaSymDockingSlideTogether" );
		scorefxn_ = new core::scoring::symmetry::SymmetricScoreFunction();
		scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
	}

 FaSymDockingSlideTogether::~FaSymDockingSlideTogether(){}

void FaSymDockingSlideTogether::apply( core::pose::Pose & pose )
{
	using namespace core::scoring;

	assert( core::pose::symmetry::is_symmetric( pose ));
  SymmetricConformation & symm_conf (
        dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

	Size num_slide_moves(1);
  std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );

	// A very hacky way of guessing whether the components are touching:
	// if pushed together by 1A, does fa_rep change at all?
	// (The docking rb_* score terms aren't implemented as of this writing.)
	(*scorefxn_)( pose );
	core::Real const initial_fa_rep = pose.energies().total_energies()[ fa_rep ];
	bool are_touching = false;
	rigid::RigidBodyDofSeqTransMover trans_mover( dofs );

	//int i=1;
	// Take 2A steps till clash, then back apart one step.  Now you're within 2A of touching.
	// Repeat with 1A steps, 0.5A steps, 0.25A steps, etc until you're as close are you want.
	for( core::Real stepsize = 2.0; stepsize > tolerance_; stepsize /= 2.0 ) {
		trans_mover.trans_axis( trans_mover.trans_axis().negate() ); // now move together
		trans_mover.step_size(stepsize);
		do
		{
			trans_mover.apply( pose );
			(*scorefxn_)( pose );
			core::Real const push_together_fa_rep = pose.energies().total_energies()[ fa_rep ];
			//std::cout << "fa_rep = " << push_together_fa_rep << std::endl;
			are_touching = (std::abs(initial_fa_rep - push_together_fa_rep) > 1e-4);
			//std::ostringstream s;
			//s << "snapshot" << i << ".pdb";
			//pose.dump_pdb(s.str());
			//i += 1;
		} while( !are_touching );
			if ( ++num_slide_moves > 1000 ) {
			std::cerr << "To many slide moves. Subunits never touching..." << std::endl;
			utility_exit();
		}

		trans_mover.trans_axis( trans_mover.trans_axis().negate() ); // now move apart
		trans_mover.apply( pose );
	}
}

std::string
FaSymDockingSlideTogether::get_name() const {
	return "FaSymDockingSlideTogether";
}


SymmetrySlider::~SymmetrySlider(){}

SymmetrySlider::SymmetrySlider( core::pose::Pose & pose )
{
	setup( pose );
}

SymmetrySlider::SymmetrySlider(
	core::pose::Pose & pose,
	core::conformation::symmetry::SlideCriteriaType score_criteria,
	std::string SlideCriteriaVal
)
{
	setup( pose );
	runtime_assert( score_criteria >= 1 && score_criteria < TOTAL_NUM_CRITERIA );
	SlideCriteriaType_ = score_criteria;
	SlideThreshold_ = "AUTOMATIC";

	if ( SlideCriteriaType_ == CEN_DOCK_SCORE ){
		core::scoring::symmetry::SymmetricScoreFunction scorefxn_sym ( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::CENTROID_WTS, core::scoring::DOCK_LOW_PATCH ) );
		scorefxn_ = new core::scoring::symmetry::SymmetricScoreFunction( scorefxn_sym );
	} else if ( SlideCriteriaType_ == FA_REP_SCORE  ) {
		scorefxn_ = new core::scoring::symmetry::SymmetricScoreFunction();
		scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
		(*scorefxn_)( pose );
		core::Real const initial_fa_rep = pose.energies().total_energies()[ core::scoring::fa_rep ];
		std::ostringstream stream;
		stream << initial_fa_rep;
		SlideThreshold_ = stream.str();
	}

	if ( SlideCriteriaVal != "AUTOMATIC" ){
		SlideThreshold_ = SlideCriteriaVal;
	}
}

void SymmetrySlider::setup( core::pose::Pose & pose )
{
	using namespace core::scoring;
	using namespace moves;

  assert( core::pose::symmetry::is_symmetric( pose ));
  SymmetricConformation & symm_conf (
        dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

	std::map< Size, core::conformation::symmetry::SymDof > dofs = symm_conf.Symmetry_Info()->get_dofs();
	std::map< Size, core::conformation::symmetry::SymDof >::iterator jump_iterator;

	// Save jumps that are allowed to move and have a translation dof
  std::map< Size, core::conformation::symmetry::SymDof >::iterator it;
  std::map< Size, core::conformation::symmetry::SymDof >::iterator it_begin = dofs.begin();
  std::map< Size, core::conformation::symmetry::SymDof >::iterator it_end = dofs.end();
  for ( it = it_begin; it != it_end; ++it ) {
    int jump_nbr ( (*it).first );
    core::conformation::symmetry::SymDof dof ( (*it).second );
    if ( dof.allow_dof(1) || dof.allow_dof(2) || dof.allow_dof(3) ) {
			AllowSlideJumpMap_[jump_nbr] = true;
    }
  }
	current_jump_ = 0;
	reset_slide_ = false;
	SlideCriteriaType_ = CEN_DOCK_SCORE;
	SlideThreshold_ = "AUTOMATIC";

	core::scoring::symmetry::SymmetricScoreFunction scorefxn_sym ( ScoreFunctionFactory::create_score_function( CENTROID_WTS, DOCK_LOW_PATCH ) );
	scorefxn_ = new core::scoring::symmetry::SymmetricScoreFunction( scorefxn_sym );

	total_num_slides_ = 0;

	// initialize the InitialJumps_ map
	std::map< core::Size, bool >::const_iterator i_it;
  std::map< core::Size, bool >::const_iterator i_it_begin = AllowSlideJumpMap_.begin();
  std::map< core::Size, bool >::const_iterator i_it_end = AllowSlideJumpMap_.end();
  for ( i_it = i_it_begin; i_it != i_it_end; ++i_it ) {
		if ( (*i_it).second ) {
			TR.Debug << "Initial Jump (nbr): " << (*i_it).first << std::endl << pose.jump ( (*i_it).first ) << std::endl;
			InitialJumps_[ (*i_it).first ]=pose.jump( (*i_it).first );

			//fpd also find the correct direction to slide for each jump
			std::map< Size, core::conformation::symmetry::SymDof > dofs = symm_conf.Symmetry_Info()->get_dofs();
			std::map< Size, core::conformation::symmetry::SymDof >::iterator dof_iterator;

			dof_iterator = dofs.find( ( *i_it).first );
			rigid::RigidBodyDofTransMover dofmover( (*dof_iterator).second, (*i_it).first, step_size() );
			InvertJump_[ (*i_it).first ]=!dofmover_compresses( pose, dofmover );
		}
  }
}

bool SymmetrySlider::finished()
{
	bool finished = true;
	std::map< core::Size, bool >::const_iterator it;
	std::map< core::Size, bool >::const_iterator it_begin = AllowSlideJumpMap_.begin();
	std::map< core::Size, bool >::const_iterator it_end = AllowSlideJumpMap_.end();
	for ( it = it_begin; it != it_end; ++it ) {
		if ( ( *it).second ) finished = false;
	}
	return finished;
}

void SymmetrySlider::set_slide_criteria(std::string SlideCriteria )
{
		SlideThreshold_ = SlideCriteria;
}

void SymmetrySlider::set_current_jump(core::Size jump_nbr )
{
	TR.Debug << "Selecting Jump " << jump_nbr << std::endl;
	current_jump_ = jump_nbr;
}

// Start by sliding away. This is tricky if we are already in contact.
void SymmetrySlider::slide_away( core::pose::Pose & pose )
{
	using namespace moves;

	std::map< core::Size, bool >::const_iterator it;
	std::map< core::Size, bool >::const_iterator it_begin = AllowSlideJumpMap_.begin();
	std::map< core::Size, bool >::const_iterator it_end = AllowSlideJumpMap_.end();
	// slide away sequentially. Why? Becuase I say so.
	for ( it = it_begin; it != it_end; ++it ) {
		if ( ( *it).second ) {
			TR.Debug << "Sliding away..." << std::endl;
			// No point is sliding if we are already very_far_away
			if ( very_far_away( pose ) ) break;
			assert( core::pose::symmetry::is_symmetric( pose ));
			SymmetricConformation & symm_conf (
      dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

			std::map< Size, core::conformation::symmetry::SymDof > dofs = symm_conf.Symmetry_Info()->get_dofs();
			std::map< Size, core::conformation::symmetry::SymDof >::iterator dof_iterator;
			dof_iterator = dofs.find( ( *it).first );
			rigid::RigidBodyDofTransMover dofmover( (*dof_iterator).second, (*it).first, step_size() );
			if (!InvertJump_[ (*it).first ]) {
				dofmover.trans_axis().negate();
			}
			core::Real prev_score = slide_score(pose);
			core::Real new_score = prev_score;
			TR.Debug << "Sliding along " << (*it).first << " starting score: " << new_score << std::endl;
			// Slide away up to 200 steps
			for (int step=1; step< 200; ++step ){
				// stop if score did not change after 5 steps
				if ( step%5 == 0 ) {
					TR.Debug << "Slide away score: " << new_score << std::endl;
					if ( std::fabs( new_score - prev_score ) < 1e-1 ) {
						TR.Debug << "Done sliding, score is unchanged " << new_score << " " << prev_score << std::endl;
						break;
					}
					else prev_score = new_score;
				}
				dofmover.apply( pose );
				new_score = slide_score(pose);
				TR.Debug << "slide away score: " << new_score << " threshold: " << get_slide_threshold()  << std::endl;
				if ( new_score <= get_slide_threshold() ) {
					TR.Debug << "Done sliding " << new_score << std::endl;
					break;
				}
			}
		}
	}
}
// The threshold metric. We are using docking interchain_vdw for now
core::Real SymmetrySlider::slide_score( core::pose::Pose & pose )
{
	if ( SlideCriteriaType_ == CEN_DOCK_SCORE ) {
		(*scorefxn_)(pose);
		 return pose.energies().total_energies()[ scoring::interchain_vdw ];
	}	else if ( SlideCriteriaType_ == FA_REP_SCORE ) {
		(*scorefxn_)(pose);
		 return pose.energies().total_energies()[ core::scoring::fa_rep ];
	}

	return 0;
}

core::Real SymmetrySlider::step_size()
{
	return 1;
}

core::Real SymmetrySlider::get_slide_threshold()
{
	return std::atof( SlideThreshold_.c_str() );
}

void SymmetrySlider::set_slide_threshold( std::string threshold )
{
	SlideThreshold_ = threshold;
}

// If the slide threshold criteria has been met we stop sliding.
// Currently only score-based criteria implemented
bool SymmetrySlider::continue_slide( core::pose::Pose & pose )
{
	total_num_slides_++;

	if ( total_num_slides_ > 200 ) {
		total_num_slides_ = 0;
		reset_slide_ = true;
		TR.Debug << "Stop sliding because we have taken more than 200 steps..." << std::endl;
		return false;
	}
	if ( SlideCriteriaType_ >= CEN_DOCK_SCORE && SlideCriteriaType_ <= FA_REP_SCORE ) {
		if ( slide_score(pose) <= get_slide_threshold() ) {
			TR.Debug << "Continue sliding, step " << total_num_slides_ << " score: " << slide_score(pose) << std::endl;
			return true;
		}
		total_num_slides_ = 0;
	}
	return false;
}

core::Size SymmetrySlider::get_current_jump()
{
	return current_jump_;
}

bool SymmetrySlider::allowed_current_slide()
{
	if ( AllowSlideJumpMap_.find( current_jump_ ) == AllowSlideJumpMap_.end() ) {
		utility_exit_with_message( "[ERROR] slide jump not found..." );
	}
	return AllowSlideJumpMap_.find( current_jump_ )->second;
}

void SymmetrySlider::disallow_current_slide()
{
	if ( AllowSlideJumpMap_.find( current_jump_ ) == AllowSlideJumpMap_.end() ) {
		utility_exit_with_message( "[ERROR] slide jump not found..." );
	}
	AllowSlideJumpMap_[current_jump_] = false;
}

std::map< core::Size, bool >
SymmetrySlider::get_allow_slide_jump_map() const
{
	return AllowSlideJumpMap_;
}

// There is no point of sliding away if the chains are already really_far_away...
bool SymmetrySlider::very_far_away( core::pose::Pose & pose )
{
	SymmetricConformation & symm_conf (
        dynamic_cast<SymmetricConformation & > ( pose.conformation()) );
  SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	conformation::Residue anchor_res ( pose.residue( core::pose::symmetry::find_symmetric_basejump_anchor( pose ) ) );
	// find the residues that are connected to the virtual residues defined by the jump. Calculate the
	// distances of the CA's for these residues and use as a metric of distance
	core::Real min_distance (1000);
	for ( std::vector< Size>::const_iterator
          clone     = symm_info->bb_clones( anchor_res.seqpos() ).begin(),
          clone_end = symm_info->bb_clones( anchor_res.seqpos() ).end();
          clone != clone_end; ++clone ) {
		conformation::Residue anchor_clone = pose.residue( *clone );
		core::Real dist = pose.residue( anchor_res.seqpos() ).xyz("CA").distance( pose.residue( anchor_clone.seqpos() ).xyz("CA") );
		TR.Debug << "Distance residue  " << anchor_res.seqpos() << " to " <<  anchor_clone.seqpos() << " is " << dist << std::endl;
		if ( dist < min_distance ) min_distance = dist;
	}
	// 50 Angstroms is pretty far innit?
	TR.Debug << "Chains are " << min_distance << " away from each other" << std::endl;
	return ( min_distance > 50 ? true : false );

}

// Slide chains into contact. Calls pure virtual select_jump() so this class can not be used on its own
void SymmetrySlider::slide(core::pose::Pose & pose)
{
	using namespace moves;

	assert( core::pose::symmetry::is_symmetric( pose ));
	SymmetricConformation & symm_conf (
		dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

	 std::map< Size, core::conformation::symmetry::SymDof > dofs = symm_conf.Symmetry_Info()->get_dofs();

	// make sure we have a current_jump_
	select_jump();
	// start by sliding away
	slide_away( pose );
	// Go through all allowed translation jumps in the order defined by
	// the selection of current_jump_
	while ( SymmetrySlider::allowed_current_slide() ) {
		std::map< Size, core::conformation::symmetry::SymDof >::iterator dof_iterator;
			while ( SymmetrySlider::continue_slide ( pose ) ) {
				// Select a new slide
				// We slide until the slide criteria is satisfied
				select_jump();
				dof_iterator = dofs.find(SymmetrySlider::get_current_jump() );
				core::conformation::symmetry::SymDof dof ( (*dof_iterator).second );
				assert( dof_iterator != dofs.end() );

				// slide along the current jump
				rigid::RigidBodyDofTransMover dofmover( dof, SymmetrySlider::get_current_jump(), SymmetrySlider::step_size() );
				//fpd pick slide direction
				if ( InvertJump_[ (*dof_iterator).first ]) {
					dofmover.trans_axis().negate();
				}
				dofmover.apply( pose );
		}

		dof_iterator = dofs.find(SymmetrySlider::get_current_jump() );
		core::conformation::symmetry::SymDof dof ( (*dof_iterator).second );
		rigid::RigidBodyDofTransMover dofmover( dof, SymmetrySlider::get_current_jump(), SymmetrySlider::step_size() );
		if (!reset_slide_) {
			for ( core::Real step=step_size()/8; step<=step_size(); step *=2 ) {
				dofmover.step_size(step);
				dofmover.apply( pose );
				TR.Debug << "Refine slide position, step " << step << std::endl;
			if (SymmetrySlider::continue_slide ( pose ) ) break;
			dofmover.trans_axis().negate();
			dofmover.apply( pose );
			dofmover.trans_axis().negate();
			}
		}

		// if we did not find a contact during the slide then reset the rb position.
		if ( reset_slide_ ) {
			// we should really keep track of how many slides we do in each direction. 200 step reset is only correct if slide
			// in only one directions fails...
			TR.Debug << "Reset slide. No contact was found after 200 steps..." << std::endl;
			std::map< core::Size, core::kinematics::Jump > const & initial_jumps( InitialJumps_ );
			std::map< core::Size, core::kinematics::Jump >::const_iterator it;
			std::map< core::Size, core::kinematics::Jump >::const_iterator it_begin = initial_jumps.begin();
			std::map< core::Size, core::kinematics::Jump >::const_iterator it_end = initial_jumps.end();
			for ( it = it_begin; it != it_end; ++it ) {
				TR.Debug << "Reset Jump (nbr): " << (*it).first << std::endl << (*it).second  << std::endl;
				pose.set_jump( (*it).first, (*it).second );
			}
			reset_slide_ = false;
		}

		// Update the threshold due gained contact. Due to round-off errors we add a fudge factor of 0.1
		core::Real score( slide_score(pose) + 0.01 );
		// The threshold value is a string. We need to convert the real to a string. Akward...
		std::ostringstream stream;
		stream << score;
		std::string new_threshold ( stream.str() );
		set_slide_threshold( new_threshold );
		// Now we are done with sliding in this direction.
		TR.Debug << "Stop sliding along Jump " << get_current_jump() << std::endl;
		SymmetrySlider::disallow_current_slide();
		// Select a new slide
		select_jump();

	}
}

void SymmetrySlider::apply(core::pose::Pose & pose)
{
	slide(pose);
}

//fpd return the rg
core::Real SymmetrySlider::rg( core::pose::Pose const & pose )
{
	utility::vector1< numeric::xyzVector <core::Real> > cas;
	numeric::xyzVector <core::Real> sum_ca(0,0,0);
	core::Real rg = 0.0;

	for (int i=1; i<=(int)pose.total_residue(); ++i) {
		core::conformation::Residue const& rsd_i = pose.residue(i);
		if (!rsd_i.is_protein()) continue;
		cas.push_back( rsd_i.atom( pose.residue_type( i ).atom_index(" CA ") ).xyz() );
		sum_ca += cas[cas.size()];
	}
	sum_ca /= cas.size();
	for (int i=1; i<=(int)cas.size(); ++i) {
		rg += (sum_ca-cas[i]).length_squared();
	}
	return std::sqrt( rg/cas.size() );
}

bool SymmetrySlider::dofmover_compresses( core::pose::Pose & pose, protocols::rigid::RigidBodyDofTransMover & dofmover )
{
	core::Real radius_before = rg(pose);
	dofmover.apply(pose);
	core::Real radius_after =  rg(pose);

	// undo changes to the pose (cheaper than working in a pose copy)
	dofmover.trans_axis().negate();
	dofmover.apply( pose );
	dofmover.trans_axis().negate();

	if (radius_before > radius_after)
		return true;
	return false;
}


SequentialSymmetrySlider::SequentialSymmetrySlider( core::pose::Pose & pose ) :
	SymmetrySlider( pose )
{
	init();
}

SequentialSymmetrySlider::SequentialSymmetrySlider( SymmetrySlider const & Slide ) :
  SymmetrySlider( Slide )
{
	init();
}

SequentialSymmetrySlider::SequentialSymmetrySlider(
	core::pose::Pose & pose,
	core::conformation::symmetry::SlideCriteriaType score_criteria,
	std::string SlideCriteriaVal
) : SymmetrySlider( pose, score_criteria, SlideCriteriaVal )
{
	init();
}

void SequentialSymmetrySlider::init()
{
	std::map< core::Size, bool > const & allow_slide_jump_map( get_allow_slide_jump_map() );
	std::map< core::Size, bool >::const_iterator it;
  std::map< core::Size, bool >::const_iterator it_begin = allow_slide_jump_map.begin();
  std::map< core::Size, bool >::const_iterator it_end = allow_slide_jump_map.end();
  for ( it = it_begin; it != it_end; ++it ) {
    if ( ( *it).second ) {
			slide_order_.push_back( (*it).first );
		}
	}
	numeric::random::random_permutation( slide_order_, RG );
}

void SequentialSymmetrySlider::select_jump()
{
	std::map< core::Size, bool > const & allow_slide_jump_map( get_allow_slide_jump_map() );
	std::map< core::Size, bool >::const_iterator allow_it;

	std::vector< core::Size >::const_iterator it;
  std::vector< core::Size >::const_iterator it_begin = slide_order_.begin();
  std::vector< core::Size >::const_iterator it_end = slide_order_.end();
  for ( it = it_begin; it != it_end; ++it ) {
		TR.Debug << "JUMP select " << (*it) << std::endl;
		allow_it = allow_slide_jump_map.find( *it );
    if ( allow_it != allow_slide_jump_map.end() ) {
			if ( (*allow_it).second ) {
				set_current_jump( *it );
				break;
			}
		}
	}

/*	std::map< core::Size, bool > const & allow_slide_jump_map( get_allow_slide_jump_map() );
	std::map< core::Size, bool >::const_iterator it;
  std::map< core::Size, bool >::const_iterator it_begin = allow_slide_jump_map.begin();
  std::map< core::Size, bool >::const_iterator it_end = allow_slide_jump_map.end();
  for ( it = it_begin; it != it_end; ++it ) {
    if ( ( *it).second ) {
			set_current_jump( ( *it).first );
			break;
		}
	}*/
}

OrderedSequentialSymmetrySlider::OrderedSequentialSymmetrySlider( core::pose::Pose & pose, std::vector<core::Size> slide_order ) :
	SymmetrySlider( pose ),
	slide_order_( slide_order)
{}

OrderedSequentialSymmetrySlider::OrderedSequentialSymmetrySlider( SymmetrySlider const & Slide, std::vector<core::Size> slide_order ) :
  SymmetrySlider( Slide ),
	slide_order_( slide_order)
{}

OrderedSequentialSymmetrySlider::OrderedSequentialSymmetrySlider(
	core::pose::Pose & pose,
	core::conformation::symmetry::SlideCriteriaType score_criteria,
	std::string SlideCriteriaVal,
	std::vector<core::Size> slide_order
) : SymmetrySlider( pose, score_criteria, SlideCriteriaVal ),
		slide_order_( slide_order)
{}

void OrderedSequentialSymmetrySlider::select_jump()
{
	std::map< core::Size, bool > const & allow_slide_jump_map( get_allow_slide_jump_map() );
	std::map< core::Size, bool >::const_iterator allow_it;

	std::vector< core::Size >::const_iterator it;
  std::vector< core::Size >::const_iterator it_begin = slide_order_.begin();
  std::vector< core::Size >::const_iterator it_end = slide_order_.end();
  for ( it = it_begin; it != it_end; ++it ) {
		allow_it = allow_slide_jump_map.find( *it );
    if ( allow_it != allow_slide_jump_map.end() ) {
			if ( (*allow_it).second ) {
				set_current_jump( *it );
				break;
			}
		}
	}
}


RandomSymmetrySlider::RandomSymmetrySlider( core::pose::Pose & pose ) :
	SymmetrySlider( pose )
{}

RandomSymmetrySlider::RandomSymmetrySlider( SymmetrySlider const & Slide ) :
  SymmetrySlider( Slide )
{}

RandomSymmetrySlider::RandomSymmetrySlider(
	core::pose::Pose & pose,
	core::conformation::symmetry::SlideCriteriaType score_criteria,
	std::string SlideCriteriaVal
) : SymmetrySlider( pose, score_criteria, SlideCriteriaVal )
{}

void RandomSymmetrySlider::select_jump()
{
	std::map< core::Size, bool > const & allow_slide_jump_map( get_allow_slide_jump_map() );
	std::map< core::Size, bool >::const_iterator it;
  std::map< core::Size, bool >::const_iterator it_begin = allow_slide_jump_map.begin();
  std::map< core::Size, bool >::const_iterator it_end = allow_slide_jump_map.end();
	utility::vector1<core::Size> allowed;
	for ( it = it_begin; it != it_end; ++it ) {
    if ( ( *it).second ) {
			allowed.push_back( ( *it).first );
		}
	}
	if ( !allowed.empty() )
		set_current_jump( numeric::random::random_element( allowed ) );
}

}
} // namespace docking
} // namespace protocols
