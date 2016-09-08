// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/simple_moves/BackboneMover.cc
/// @brief  Method definitions for SmallMover and ShearMover

// Unit headers
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/BackboneMoverCreator.hh>

// Project Headers
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/id/TorsionID_Range.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

#include <protocols/rosetta_scripts/util.hh>

// Basic headers
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/NumericTraits.hh>
#include <numeric/random/random.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>

// C++ headers
#include <string>

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.BackboneMover" );

namespace protocols {
namespace simple_moves {


using namespace core;

//constructor
BackboneMover::BackboneMover() :
	protocols::canonical_sampling::ThermodynamicMover(),
	movemap_(),
	scorefxn_(),
	temperature_( 0.5 ),
	nmoves_( 1 ),
	resnum_( 0 ),
	big_angle_( 0 ),
	small_angle_( 0 ),
	pos_list_( ),
	old_phi_(0), new_phi_(0), old_psi_(0), new_psi_(0), old_rama_score_(0), new_rama_score_(0),
	preserve_detailed_balance_(false),
	selector_()
{
	moves::Mover::type( "BackboneMoverBase" );
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_bb( true ); // allow all backbone residues to move
	movemap_=movemap; // and make this a constant object
	angle_max( 'H', 0.0 ); // helix
	angle_max( 'L', 6.0 ); // other
	angle_max( 'E', 5.0 ); // strand
}

BackboneMover::BackboneMover(
	core::kinematics::MoveMapOP movemap_in,
	core::Real temperature_in,
	core::Size nmoves_in
) : protocols::canonical_sampling::ThermodynamicMover(),
	movemap_(std::move( movemap_in )),
	scorefxn_(),
	temperature_( temperature_in),
	nmoves_( nmoves_in ),
	resnum_( 0 ),
	big_angle_( 0 ),
	small_angle_( 0 ),
	pos_list_( ),
	old_phi_(0), new_phi_(0), old_psi_(0), new_psi_(0), old_rama_score_(0), new_rama_score_(0),
	preserve_detailed_balance_(false),
	selector_()
{
	moves::Mover::type( "BackboneMoverBase" );
	// set default values for angle_max since it is not passed in
	angle_max( 'H', 0.0 ); // helix
	angle_max( 'L', 6.0 ); // other
	angle_max( 'E', 5.0 ); // strand
}

/// @brief Copy constructor
///
BackboneMover::BackboneMover( BackboneMover const &src ) :
	protocols::canonical_sampling::ThermodynamicMover( src ),
	movemap_(), // Copied below
	scorefxn_(), // Cloned below
	temperature_( src.temperature_ ),
	nmoves_( src.nmoves_ ),
	angle_max_( src.angle_max_),
	resnum_( src.resnum_ ),
	big_angle_( src.big_angle_ ),
	small_angle_( src.small_angle_ ),
	pos_list_( src.pos_list_ ),
	old_phi_(src.old_phi_), new_phi_(src.new_phi_), old_psi_(src.old_psi_), new_psi_(src.new_psi_), old_rama_score_(src.old_rama_score_), new_rama_score_(src.new_rama_score_),
	preserve_detailed_balance_( src.preserve_detailed_balance_ ),
	selector_() // Cloned below
{
	if ( src.scorefxn_ ) {
		scorefxn_ = src.scorefxn_->clone();
	}
	if ( src.movemap_ ) {
		movemap_ = src.movemap_;
	} else {
		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		movemap->set_bb( true ); // allow all backbone residues to move
		movemap_=movemap; // and make this a constant object
	}
	if ( src.selector_ ) {
		selector_ = src.selector_->clone();
	}
}

//destructor
BackboneMover::~BackboneMover() = default;


void BackboneMover::temperature( core::Real const temperature_in ) { temperature_ = temperature_in; }

core::Real
BackboneMover::temperature() const {
	return temperature_;
}

void BackboneMover::nmoves( core::Size const nmoves_in ) { nmoves_ = nmoves_in; }

core::Size
BackboneMover::nmoves() const {
	return nmoves_;
}

core::kinematics::MoveMapCOP BackboneMover::movemap() const { return movemap_; }

void BackboneMover::movemap(core::kinematics::MoveMapOP new_movemap) { movemap_=new_movemap; }

core::select::residue_selector::ResidueSelectorCOP
BackboneMover::selector() const { return selector_; }

void BackboneMover::selector( core::select::residue_selector::ResidueSelectorCOP selector )
{
	selector_ = selector;
}

core::scoring::ScoreFunctionOP
BackboneMover::scorefxn() const { return scorefxn_; }

void
BackboneMover::scorefxn( core::scoring::ScoreFunctionOP sfx ) { scorefxn_ = sfx; }

void BackboneMover::angle_max( core::Real const angle )
{
	angle_max( 'H', angle ); // helix
	angle_max( 'L', angle ); // other
	angle_max( 'E', angle ); // strand
}

void BackboneMover::angle_max( char const type, core::Real const angle ) { angle_max_[ type ] = angle; }

void BackboneMover::angle_max( std::map< char, core::Real > angle_max_in ) { angle_max_.swap( angle_max_in ); }

core::Real
BackboneMover::get_angle_max(char const type) const
{
	return angle_max_.find(type)->second;
}


bool
BackboneMover::preserve_detailed_balance() const
{
	return preserve_detailed_balance_;
}

void
BackboneMover::set_preserve_detailed_balance(
	bool preserve_detailed_balance
)
{
	preserve_detailed_balance_ = preserve_detailed_balance;
}

core::Real BackboneMover::new_phi() { return new_phi_; }

core::Real BackboneMover::new_psi() { return new_psi_; }

/// @brief Set the ResidueSelector that this mover will use.
/// @details Clones the input.
void
BackboneMover::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector
) {
	runtime_assert_string_msg( selector, "Error in protocols::simple_moves::BackboneMover::set_residue_selector(): An invalid ResidueSelector owning pointer was passed to this function." );
	selector_ = selector->clone();
	return;
}

void
BackboneMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & pose
)
{
	if ( tag->hasOption( "residue_selector" ) ) {
		core::select::residue_selector::ResidueSelectorCOP res_selector =
			protocols::rosetta_scripts::parse_residue_selector( tag, data );
		if ( !res_selector ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "ResidueSelector passed to Shear or Small backbone mover could not be found." );
		}
		selector( res_selector );
	}
	if ( tag->hasOption( "scorefxn" ) ) {
		core::scoring::ScoreFunctionOP sfx =
			protocols::rosetta_scripts::parse_score_function( tag, data );
		if ( !sfx ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "Scorefxn passed to Shear or Small backbone mover could not be found." );
		}
		scorefxn( sfx->clone() );
	}

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	protocols::rosetta_scripts::parse_movemap( tag, pose, mm, data );
	movemap( mm );
	temperature( tag->getOption<core::Real>( "temperature", temperature_ ) );
	nmoves( tag->getOption<core::Size>( "nmoves", nmoves_ ) );
	angle_max( tag->getOption<core::Real>( "angle_max", 6.0 ) );
	set_preserve_detailed_balance( tag->getOption<bool>( "preserve_detailed_balance", preserve_detailed_balance_ ) );
}

void BackboneMover::apply( core::pose::Pose & pose )
{
	// clear everything to start from clean
	clear();

	//what if symmetry input
	// will not run for now, can change this to true if symmetry desired
	// potential problem with incoming modified symmetric movemaps, this will
	// reset any changes to the input movemap
	if ( false ) {
		if ( core::pose::symmetry::is_symmetric( pose )  )  {
			//TR << "Pose is symmetric, making backbone movemap symmetric."<< std::endl;
			core::pose::symmetry::make_symmetric_movemap( pose, *movemap_);
		}
	}

	// currently secstruct is not implemented and returns 'L' for all residues
	setup_list( pose );

	if ( pos_list_.empty() ) {
		Warning() << "no movable positions in " << type() << "!" << std::endl;
		return;
	}

	// how many moves to make
	int const num_iterations = std::max<Size>( Size(1), std::min<Size>( nmoves_, pos_list_.size()/2 ) );
	std::set< int > already_moved;
	int tries = 0;

	// now loop
	for ( int k=1; k<=num_iterations; ++k ) {
		while ( true ) {
			++tries;

			// give up after trying 10000 times
			if ( tries > 10000 ) {
				break;
			}

			//choose a random position from the list of positions previously chose to be candidate positions
			std::pair< int, Real > const & random_pos( pos_list_[ static_cast< int >( numeric::random::rg().uniform() * pos_list_.size() + 1 ) ] );
			resnum_ = random_pos.first;

			set_angles( random_pos.second );

			//  next three lines skip ends of structure !!
			//  fix a logic error here: the closer to the end, the higher the probability
			//  to skip it; and when it is 1 or total_residue, the probability should be
			//  zero; and the probability distribution should be symmetrical for two ends
			//  residue:    N-   1   2   3   4   5   6
			//  residue:    C-   t  t-1 t-2 t-3 t-4 t-5
			//  prob to skip:    1  0.8 0.6 0.4 0.2 0.0
			//  -- chu

			// need to add this back, prob want a pose.distance_to_chain_end(i)
			// function
			//
			//end = total_residue - 5;
			//if ( resnum <= 5 && static_cast< int >(ran3()*5+1) >= resnum ) goto L401;
			//if ( resnum > end && static_cast< int >(ran3()*5) + end <= resnum ) goto L401;

			// maybe we've already moved this position ?
			if ( already_moved.find( resnum_ ) != already_moved.end() ) {
				continue;
			}

			if ( !make_move( pose ) ) continue;
			already_moved.insert( resnum_ );
			break;
		} // while ( true )

	} // k=1,num
}

std::string
BackboneMover::get_name() const {
	return "BackboneMover";
}

void
BackboneMover::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Max angle for helices (H): " << get_angle_max('H') <<
		"\nMax angle for strands (E): " << get_angle_max('E') <<
		"\nMax angle for loops (L):   " << get_angle_max('L') <<
		"\nTemperature factor (kT):   " << temperature() <<
		"\nNumber of moves:           " << nmoves() << std::endl;
	output << "MoveMap:" << std::endl;
	movemap()->show(output);
}

void BackboneMover::clear() {
	pos_list_.erase(pos_list_.begin(), pos_list_.end());
}

bool BackboneMover::check_rama() {
	if ( preserve_detailed_balance_ ) return true;
	if ( new_rama_score_ > old_rama_score_ ) {
		Real const boltz_factor = (old_rama_score_-new_rama_score_)/temperature_;
		Real const probability = std::exp(std::max(Real(-40.0),boltz_factor) );
		if ( numeric::random::rg().uniform() >= probability ) return( false );
	}
	return( true );
}

core::select::residue_selector::ResidueSubset
BackboneMover::compute_selected_residues( core::pose::Pose const & pose ) const
{
	if ( selector_ ) {
		return selector_->apply( pose );
	} else {
		return core::select::residue_selector::ResidueSubset( pose.size(), true );
	}
}

std::ostream &operator<< (std::ostream &os, BackboneMover const &mover)
{
	mover.show(os);
	return os;
}


std::string
SmallMoverCreator::keyname() const {
	return SmallMoverCreator::mover_name();
}

protocols::moves::MoverOP
SmallMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SmallMover );
}

std::string
SmallMoverCreator::mover_name() {
	return "Small";
}


//constructor
SmallMover::SmallMover() : BackboneMover() { moves::Mover::type( "SmallMover" ); }

SmallMover::SmallMover(
	core::kinematics::MoveMapOP movemap_in,
	core::Real temperature_in,
	core::Size nmoves_in
)
: BackboneMover( movemap_in, temperature_in, nmoves_in )
{
	moves::Mover::type( "Small" );
}

//destructor
SmallMover::~SmallMover() = default;

std::string
SmallMover::get_name() const {
	return "Small";
}

protocols::moves::MoverOP
SmallMover::clone() const {
	return protocols::moves::MoverOP( new SmallMover(*this) );
}

protocols::moves::MoverOP
SmallMover::fresh_instance() const {
	return protocols::moves::MoverOP( new SmallMover() );
}

void SmallMover::setup_list( core::pose::Pose & pose )
{
	//Mask by ResidueSelector:
	core::select::residue_selector::ResidueSubset const subset = compute_selected_residues( pose );

	using namespace id;
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( !subset[i] ) continue; //Skip residues masked by the ResidueSelector.
		conformation::Residue const & rsd( pose.residue( i ) );

		//Checks if the residue is protein and has a free psi and phi angle as determined by the move map
		if ( rsd.is_protein() && movemap()->get( TorsionID( i, BB, phi_torsion ) ) &&
				movemap()->get( TorsionID( i, BB, psi_torsion ) ) ) {
			//Gets the secondary structure nature of the residue
			char const ss( pose.secstruct( i ) );
			if ( angle_max_.count( ss ) ) {
				//Checks the maximum angle deviation for the specific kind of secondary structure and if it is greater than 0 pushes back the residue and
				//the maximum angle deviation to the list of movable residues.
				Real const mx( angle_max_.find( ss )->second );
				if ( mx > 0.0 ) {
					pos_list_.push_back( std::make_pair( i, mx ) );
				}
			}
			// Check if the residue is a monosaccharide and has a free phi (BB 5) and psi (BB 4).
		} else if ( rsd.is_carbohydrate() && movemap()->get(TorsionID(i, BB, 4)) &&
				movemap()->get(TorsionID(i, BB, 5)) ) {
			// Carbohydrates are always considered loops for now.
			Real const mx = angle_max_.find('L')->second;
			if ( mx > 0.0 ) {
				pos_list_.push_back(std::make_pair(i, mx));
			}
		}
	}
}

utility::vector1<core::id::TorsionID_Range>
SmallMover::torsion_id_ranges(
	core::pose::Pose & //pose
)
{
	return utility::vector1<core::id::TorsionID_Range>();
}

utility::vector1<core::id::DOF_ID_Range>
SmallMover::dof_id_ranges(
	core::pose::Pose & pose
)
{
	Real static const pi(numeric::NumericTraits<Real>::pi());
	clear();
	setup_list( pose );

	utility::vector1<core::id::DOF_ID_Range> range_vector;

	for ( core::Size i = 1; i <= pos_list_.size(); ++i ) {

		core::id::TorsionID phi_torsion(pos_list_[i].first, core::id::BB, core::id::phi_torsion);
		core::id::DOF_ID phi_dof(pose.conformation().dof_id_from_torsion_id(phi_torsion));
		if ( phi_dof.valid() ) range_vector.push_back(core::id::DOF_ID_Range(phi_dof, -pi, pi));

		core::id::TorsionID psi_torsion(pos_list_[i].first, core::id::BB, core::id::psi_torsion);
		core::id::DOF_ID psi_dof(pose.conformation().dof_id_from_torsion_id(psi_torsion));
		if ( psi_dof.valid() ) range_vector.push_back(core::id::DOF_ID_Range(psi_dof, -pi, pi));
	}

	return range_vector;
}

void SmallMover::set_angles( core::Real angle_in ) {
	///< if anyone can think of how to explain this better, please do so!!
	big_angle_ = angle_in; ///< this number comes from max angle, it is the maximum
	///< total deviation
	small_angle_ = big_angle_/2.0; ///< this is max_angle/2, which is the deviation from the angle input

}

bool SmallMover::move_with_scorefxn( core::pose::Pose & pose )
{
	debug_assert( scorefxn() );
	old_rama_score_ = (*scorefxn())( pose );
	pose.set_phi( resnum_, new_phi_ );
	pose.set_psi( resnum_, new_psi_ );
	pose.energies().clear();
	new_rama_score_ = (*scorefxn())( pose );
	TR.Debug << "Using score function to evaluate ramas. Old=" << old_rama_score_ << " New=" << new_rama_score_
		<< " Temp=" << temperature() << std::endl;

	if ( check_rama() ) return true;

	pose.set_phi( resnum_, old_phi_ );
	pose.set_psi( resnum_, old_psi_ );
	pose.energies().clear();
	TR.Debug << "Reject: reverted score = " << (*scorefxn())( pose ) << std::endl;
	return false;
}

bool SmallMover::move_with_rama( core::pose::Pose & pose )
{
	scoring::Ramachandran const & rama( scoring::ScoringManager::get_instance()->get_Ramachandran() );
	conformation::Residue const & current_rsd( pose.residue(resnum_) );

	// Always accept carbohydrate moves for now....
	if ( current_rsd.is_protein() && current_rsd.aa() != chemical::aa_unk ) {
		old_rama_score_ = rama.eval_rama_score_residue( current_rsd.aa(), old_phi_, old_psi_ );
		new_rama_score_ = rama.eval_rama_score_residue( current_rsd.aa(), new_phi_, new_psi_ );

		// decide whether to accept the move
		if ( !check_rama() ) return false;
	}

	// set the new values for residue resnum
	pose.set_phi( resnum_, new_phi_ );
	pose.set_psi( resnum_, new_psi_ );

	return true;
}

bool SmallMover::make_move( core::pose::Pose & pose )
{
	old_phi_ = pose.phi(resnum_);
	new_phi_ = basic::periodic_range( old_phi_ - small_angle_ + numeric::random::rg().uniform() * big_angle_, 360.0 );

	old_psi_ = pose.psi(resnum_);
	new_psi_ = basic::periodic_range( old_psi_ - small_angle_ + numeric::random::rg().uniform() * big_angle_, 360.0 );

	if ( scorefxn() ) {
		return move_with_scorefxn( pose );
	} else {
		return move_with_rama( pose );
	}
}

void SmallMover::test_move( core::pose::Pose & pose)
{
	kinematics::MoveMapOP mmap( new kinematics::MoveMap() );
	mmap->set_chi( true );
	mmap->set_bb( true );

	movemap(mmap);

	apply(pose);
}

std::string
ShearMoverCreator::keyname() const {
	return protocols::simple_moves::ShearMoverCreator::mover_name();
}

protocols::moves::MoverOP
ShearMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new protocols::simple_moves::ShearMover );
}

std::string
ShearMoverCreator::mover_name() {
	return "Shear";
}

//constructor

ShearMover::ShearMover() : BackboneMover() { protocols::moves::Mover::type( "ShearMover" ); }

ShearMover::ShearMover(
	core::kinematics::MoveMapOP movemap_in,
	core::Real temperature_in,
	core::Size nmoves_in
)
: BackboneMover( movemap_in, temperature_in, nmoves_in )
{
	protocols::moves::Mover::type( "Shear" );
}

//destructor
ShearMover::~ShearMover() = default;

std::string
ShearMover::get_name() const {
	return "Shear";
}

protocols::moves::MoverOP
ShearMover::clone() const {
	return protocols::moves::MoverOP( new protocols::simple_moves::ShearMover(*this) );
}

protocols::moves::MoverOP
ShearMover::fresh_instance() const {
	return protocols::moves::MoverOP( new protocols::simple_moves::ShearMover() );
}

void protocols::simple_moves::ShearMover::setup_list( core::pose::Pose & pose )
{
	//Mask by ResidueSelector:
	core::select::residue_selector::ResidueSubset const subset = compute_selected_residues( pose );

	using namespace id;
	// Compare code below to that for SmallMover above.
	for ( core::Size i=2; i<=pose.size(); ++i ) {
		// both residue i and residue i-1 must be selected
		if ( !subset[i] ) continue;
		if ( !subset[i-1] ) continue;

		conformation::Residue const & rsd( pose.residue( i ) );
		if ( rsd.is_protein() && movemap()->get( TorsionID( i, BB, phi_torsion ) /*phi of i*/) &&
				movemap()->get( TorsionID( i-1, BB, psi_torsion ) /*psi of i-1*/) ) {
			char const ss( pose.secstruct( i ) );
			if ( angle_max_.count( ss ) ) {
				Real const mx( angle_max_.find( ss )->second );
				if ( mx > 0.0 ) {
					pos_list_.push_back( std::make_pair( i, mx ) );
				}
			}
			// Check if the residue is a monosaccharide and has a free phi (BB 5) and psi (BB 4).
		} else if ( rsd.is_carbohydrate() && movemap()->get(TorsionID(i - 1, BB, 4)) &&
				movemap()->get(TorsionID(i, BB, 5)) ) {
			// Carbohydrates are always considered loops for now.
			Real const mx = angle_max_.find('L')->second;
			if ( mx > 0.0 ) {
				pos_list_.push_back(std::make_pair(i, mx));
			}
		}
	}
}

utility::vector1<core::id::TorsionID_Range>
ShearMover::torsion_id_ranges(
	core::pose::Pose & //pose
)
{
	return utility::vector1<core::id::TorsionID_Range>();
}

utility::vector1<core::id::DOF_ID_Range>
ShearMover::dof_id_ranges(
	core::pose::Pose & pose
)
{
	Real static const pi(numeric::NumericTraits<Real>::pi());

	clear();
	setup_list( pose );

	utility::vector1<core::id::DOF_ID_Range> range_vector;

	for ( core::Size i = 1; i <= pos_list_.size(); ++i ) {

		core::id::TorsionID phi_torsion(pos_list_[i].first, core::id::BB, core::id::phi_torsion);
		core::id::DOF_ID phi_dof(pose.conformation().dof_id_from_torsion_id(phi_torsion));
		if ( phi_dof.valid() ) range_vector.push_back(core::id::DOF_ID_Range(phi_dof, -pi, pi));

		core::id::TorsionID psi_torsion(pos_list_[i].first-1, core::id::BB, core::id::psi_torsion);
		core::id::DOF_ID psi_dof(pose.conformation().dof_id_from_torsion_id(psi_torsion));
		if ( phi_dof.valid() ) range_vector.push_back(core::id::DOF_ID_Range(psi_dof, -pi, pi));
	}

	return range_vector;
}

// db (from rosetta++) set maximum deviation to be twice that of small move since shears are less perturbing
void protocols::simple_moves::ShearMover::set_angles( core::Real angle_in ) {
	///< if anyone can think of how to explain this, please do so!!
	big_angle_ = angle_in*2.0;
	small_angle_ = big_angle_/2.0;
}


bool
ShearMover::move_with_scorefxn( core::pose::Pose & pose )
{
	old_rama_score_ = (*scorefxn())( pose );
	pose.set_phi( resnum_, new_phi_ );
	pose.set_psi( resnum_-1, new_psi_ );
	pose.energies().clear();
	new_rama_score_ = (*scorefxn())( pose );

	if ( check_rama() ) return true;

	pose.set_phi( resnum_, old_phi_ );
	pose.set_psi( resnum_-1, old_psi_ );
	pose.energies().clear();
	return false;
}

bool
ShearMover::move_with_rama( core::pose::Pose & pose )
{
	scoring::Ramachandran const & rama( scoring::ScoringManager::get_instance()->get_Ramachandran() );
	// grab the residues
	conformation::Residue const & current_rsd( pose.residue(resnum_) );
	conformation::Residue const & prev_rsd( pose.residue(resnum_-1) );

	// Always accept carbohydrate moves for now....
	if ( current_rsd.is_protein() &&  current_rsd.aa() != chemical::aa_unk ) {
		// rama for phi of resnum and psi of resnum-1
		old_rama_score_ =  rama.eval_rama_score_residue( current_rsd.aa(), old_phi_, pose.psi(resnum_));
		new_rama_score_ =  rama.eval_rama_score_residue( current_rsd.aa(), new_phi_, pose.psi(resnum_));

		// decide whether to accept the move
		if ( !check_rama() ) return false;
	}
	if ( current_rsd.is_protein() &&  prev_rsd.aa() != chemical::aa_unk ) {
		// rama for residue resnum-1
		old_rama_score_ =  rama.eval_rama_score_residue( prev_rsd.aa(), pose.phi(resnum_-1), old_psi_ );
		new_rama_score_ =  rama.eval_rama_score_residue( prev_rsd.aa(), pose.phi(resnum_-1), new_psi_ );

		// decide whether to accept the move
		if ( !check_rama() ) return false;
	}

	// set the new values phi of resnum and psi of resnum-1
	pose.set_phi( resnum_, new_phi_ );
	pose.set_psi( resnum_-1, new_psi_ );

	return true;
}

bool protocols::simple_moves::ShearMover::make_move( core::pose::Pose & pose )
{
	Real shear_delta = small_angle_ - numeric::random::rg().uniform() * big_angle_;
	old_phi_ = pose.phi(resnum_);
	new_phi_ = basic::periodic_range( old_phi_ - shear_delta, 360.0 );
	old_psi_ = pose.psi(resnum_-1);
	new_psi_ = basic::periodic_range( old_psi_ + shear_delta, 360.0 );

	if ( scorefxn() ) {
		return move_with_scorefxn( pose );
	} else {
		return move_with_rama( pose );
	}
}

void protocols::simple_moves::ShearMover::test_move( core::pose::Pose & pose)
{
	kinematics::MoveMapOP mmap( new kinematics::MoveMap() );
	mmap->set_chi( true );
	mmap->set_bb( true );

	movemap(mmap);

	apply(pose);
}

}  // namespace moves
}  // namespace protocols
