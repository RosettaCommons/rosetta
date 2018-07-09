// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    src/protocols/simple_moves/BackboneMover.cc
/// @brief   Method definitions for SmallMover and ShearMover
/// @author  Labonte <JWLabonte@jhu.edu> (major refactoring to allow for non-peptide moves)


// Unit headers
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/BackboneMoverCreator.hh>

// Project Headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/carbohydrates/util.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>

#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/id/TorsionID_Range.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>

// Basic headers
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/NumericTraits.hh>
#include <numeric/random/random.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <string>


static basic::Tracer TR( "protocols.simple_moves.BackboneMover" );


namespace protocols {
namespace simple_moves {


using namespace core;

// Constructor
BackboneMover::BackboneMover() :
	protocols::canonical_sampling::ThermodynamicMover(),
	movemap_(),
	scorefxn_(),
	temperature_( 0.5 ),
	nmoves_( 1 ),
	resnum_( 0 ),
	big_angle_( 0 ),
	small_angle_( 0 ),
	moving_torsions_(),
	move_pos_list_( ),
	old_rama_score_(0), new_rama_score_(0),
	preserve_detailed_balance_( false ),
	selector_()
{
	moves::Mover::type( "BackboneMoverBase" );
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
	moving_torsions_(),
	move_pos_list_( ),
	old_rama_score_(0), new_rama_score_(0),
	preserve_detailed_balance_(false),
	selector_()
{
	moves::Mover::type( "BackboneMoverBase" );
	// Set default values for angle_max since it is not passed in.
	angle_max( 'H', 0.0 ); // helix
	angle_max( 'L', 6.0 ); // other
	angle_max( 'E', 5.0 ); // strand
}

/// @brief Copy constructor
BackboneMover::BackboneMover( BackboneMover const &src ) :
	protocols::canonical_sampling::ThermodynamicMover( src ),
	movemap_factory_( src.movemap_factory_ ),
	movemap_( src.movemap_ ),
	scorefxn_(), // Cloned below
	temperature_( src.temperature_ ),
	nmoves_( src.nmoves_ ),
	angle_max_( src.angle_max_),
	resnum_( src.resnum_ ),
	big_angle_( src.big_angle_ ),
	small_angle_( src.small_angle_ ),
	moving_torsions_( src.moving_torsions_ ),
	move_pos_list_( src.move_pos_list_ ),
	old_rama_score_(src.old_rama_score_), new_rama_score_(src.new_rama_score_),
	preserve_detailed_balance_( src.preserve_detailed_balance_ ),
	selector_() // Cloned below
{
	if ( src.scorefxn_ ) {
		scorefxn_ = src.scorefxn_->clone();
	}
	if ( src.selector_ ) {
		selector_ = src.selector_->clone();
	}
}

// Destructor.
BackboneMover::~BackboneMover() = default;


void
BackboneMover::angle_max( core::Real const angle )
{
	angle_max( 'H', angle ); // helix
	angle_max( 'L', angle ); // other
	angle_max( 'E', angle ); // strand
}

void
BackboneMover::angle_max( char const type, core::Real const angle )
{
	angle_max_[ type ] = angle;
}

void
BackboneMover::angle_max( std::map< char, core::Real > angle_max_in )
{
	angle_max_.swap( angle_max_in );
}

core::Real
BackboneMover::get_angle_max(char const type) const
{
	return angle_max_.find(type)->second;
}


/// @details Clones the input.
void
BackboneMover::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector )
{
	runtime_assert_string_msg( selector, "Error in protocols::simple_moves::BackboneMover::set_residue_selector(): "
		"An invalid ResidueSelector owning pointer was passed to this function." );
	selector_ = selector->clone();
}


void
BackboneMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const &
)
{
	if ( tag->hasOption( "residue_selector" ) ) {
		core::select::residue_selector::ResidueSelectorCOP res_selector =
			protocols::rosetta_scripts::parse_residue_selector( tag, data );
		if ( !res_selector ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "ResidueSelector passed to Shear or Small backbone mover could not be found." );
		}
		selector( res_selector );
	}
	if ( tag->hasOption( "scorefxn" ) ) {
		core::scoring::ScoreFunctionOP sfx =
			protocols::rosetta_scripts::parse_score_function( tag, data );
		if ( !sfx ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Scorefxn passed to Shear or Small backbone mover could not be found." );
		}
		scorefxn( sfx->clone() );
	}

	movemap_factory( protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data ) );
	temperature( tag->getOption<core::Real>( "temperature", temperature_ ) );
	nmoves( tag->getOption<core::Size>( "nmoves", nmoves_ ) );
	angle_max( tag->getOption<core::Real>( "angle_max", 6.0 ) );
	set_preserve_detailed_balance( tag->getOption<bool>( "preserve_detailed_balance", preserve_detailed_balance_ ) );
}

void
BackboneMover::apply( core::pose::Pose & pose )
{
	using namespace std;

	// Clear everything to start from clean.
	clear();

	setup_list( pose );

	if ( move_pos_list_.empty() ) {
		TR.Warning << "no movable positions in " << type() << "!" << endl;
		return;
	}

	// how many moves to make
	Size const num_iterations = std::max< Size >( Size( 1 ), std::min< Size >( nmoves_, move_pos_list_.size() ) );
	set< Size > already_moved;
	Size tries( 0 );

	// Now loop.
	for ( Size k = 1; k <= num_iterations; ++k ) {
		while ( true ) {
			++tries;

			// Give up after trying 10,000 times.
			if ( tries > 10000 ) { break; }

			// Choose a random position from the list of positions previously chosen to be candidate positions.
			core::Real max_angle;
			tuple< core::uint, utility::vector1< core::id::TorsionID >, Real > const & random_pos(
				move_pos_list_[ numeric::random::rg().uniform() * move_pos_list_.size() + 1 ] );
			tie( resnum_, moving_torsions_, max_angle ) = random_pos;

			set_angles( max_angle );

			// Maybe we've already moved this position ?
			if ( already_moved.find( resnum_ ) != already_moved.end() ) { continue; }

			if ( ! make_move( pose ) ) { continue; }
			already_moved.insert( resnum_ );
			break;
		}  // while ( true )
	}  // k = 1, num
}

bool
BackboneMover::make_move( core::pose::Pose & pose )
{
	if ( scorefxn() ) {
		return move_with_scorefxn( pose );
	} else {
		return move_with_rama( pose );
	}
}

void
BackboneMover::test_move( core::pose::Pose & pose )
{
	kinematics::MoveMapOP mmap( new kinematics::MoveMap() );
	mmap->set_chi( true );
	mmap->set_bb( true );

	movemap( mmap );

	apply( pose );
}

std::string
BackboneMover::get_name() const {
	return "BackboneMover";
}

void
BackboneMover::show( std::ostream & output ) const
{
	Mover::show( output );
	output << "Max angle for helices (H): " << get_angle_max( 'H' ) <<
		"\nMax angle for strands (E): " << get_angle_max( 'E' ) <<
		"\nMax angle for loops (L):   " << get_angle_max( 'L' ) <<
		"\nTemperature factor (kT):   " << temperature() <<
		"\nNumber of moves:           " << nmoves() << std::endl;
}

void
BackboneMover::clear()
{
	move_pos_list_.erase( move_pos_list_.begin(), move_pos_list_.end() );
}

bool
BackboneMover::check_rama()
{
	if ( preserve_detailed_balance_ ) { return true; }
	if ( new_rama_score_ > old_rama_score_ ) {
		Real const boltz_factor = ( old_rama_score_ - new_rama_score_ ) / temperature_;
		Real const probability = std::exp( std::max( Real( -40.0 ), boltz_factor ) );
		if ( numeric::random::rg().uniform() >= probability ) { return false; }
	}
	return true;
}

core::kinematics::MoveMapCOP
BackboneMover::movemap( core::pose::Pose const & pose ) {
	if ( movemap_ ) {
		return movemap_;
	} else if ( movemap_factory_ ) {
		return movemap_factory_->create_movemap_from_pose( pose );
	} else {
		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		movemap->set_bb( true ); // allow all backbone residues to move
		return movemap;
	}
}

utility::tag::XMLSchemaComplexTypeGeneratorOP
BackboneMover::complex_type_generator_for_backbone_mover( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	XMLSchemaComplexTypeGeneratorOP ct_gen( new XMLSchemaComplexTypeGenerator );
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute(
		"temperature", xsct_real,
		"What MC acceptance temperature to use (tests only the rama score, so not a full MC)" )
		+ XMLSchemaAttribute(
		"nmoves", xsct_non_negative_integer,
		"How many consecutive moves to make" )
		+ XMLSchemaAttribute::attribute_w_default(
		"angle_max", xsct_real,
		"By how much to perturb the backbone",
		"6.0" )
		+ XMLSchemaAttribute(
		"preserve_detailed_balance", xsct_rosetta_bool,
		"If set to true, does not test the MC acceptance criterion, and instead always accepts");

	rosetta_scripts::attributes_for_parse_score_function( attlist );
	core::select::residue_selector::attributes_for_parse_residue_selector(
		attlist, "residue_selector",
		"An optional, previously-defined ResidueSelector, specifying the subset "
		"of residues to which the mover will be applied. If not provided, "
		"the mover is applied to the whole pose. "
		"(Alternatively, a MoveMap may be used -- see below)" );
	//get subelement for parse movemap
	XMLSchemaSimpleSubelementList subelements;
	rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy( xsd, subelements );


	ct_gen->complex_type_naming_func( & protocols::moves::complex_type_name_for_mover )
		.set_subelements_repeatable( subelements )
		.add_attributes( attlist )
		.add_optional_name_attribute(
		"BackboneMover class has elements of the MC temperature to do repetitions of bb moves "
		"(small, shear, wobble, etc.)");

	return ct_gen;
}

// AMW: this is probably unnecessarily + a curse
std::string
BackboneMover::backbone_mover_complex_type_namer( std::string tag_name ){
	return "backbone_mover_" + tag_name + "_complex_type";
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

/// @brief   Return a list of TorsionIDs for the standard main-chain torsions of this residue.
/// @details Every residue family has its own set of standard (usually named) main-chain torsions, which have a variety
///          of TorsionIDs "under the hood".  This function gets these IDs from the residue.
/// @remark  Ultimately, I plan to finish a project to store each of these IDs in ResidueType.  When and if that
///          happens, this helper function can be replaced by a direct call to ResidueType::get_mainchain_TorsionIDs()
///          instead.  ~Labonte
/// @note    The Conformation is only needed for the carbohydrate case, as its main-chain torsions are not even a part
///          of the Residue being queried from Rosetta's point of view.
/// @author  Labonte <JWLabonte@jhu.edu>
utility::vector1< core::id::TorsionID >
BackboneMover::get_mainchain_TorsionIDs( core::conformation::Conformation const & conf, core::uint const seqpos ) const
{
	using namespace core::id;
	using namespace core::conformation;

	Residue const & rsd( conf.residue( seqpos ) );
	utility::vector1< core::id::TorsionID > torsions;
	if ( rsd.is_protein() && ! rsd.type().is_oligourea() ) {
		torsions.push_back( TorsionID( seqpos, BB, phi_torsion ) );
		torsions.push_back( TorsionID( seqpos, BB, psi_torsion ) );
	} else if ( rsd.type().is_oligourea() ) {
		torsions.push_back( TorsionID( seqpos, BB, phi_torsion_oligourea ) );
		torsions.push_back( TorsionID( seqpos, BB, psi_torsion_oligourea ) );
		torsions.push_back( TorsionID( seqpos, BB, theta_torsion_oligourea ) );
	} else if ( rsd.is_carbohydrate() ) {
		// Because a saccharide i's glycosidic torsions do not even belong to residue i from Rosetta's point of
		// view, we need a helper function to get the appropriate TorsionIDs and check if they are each movable.
		torsions = carbohydrates::get_glycosidic_TorsionIDs( conf, seqpos );
	}
	return torsions;
}

/// @brief    Is this set of torsions allowed to move by the MoveMap?
/// @details  Returns false if any of the torsions in the passed list are disallowed OR if the list is empty.
/// @author   Labonte <JWLabonte@jhu.edu>
bool
BackboneMover::are_torsions_allowed(
	utility::vector1< core::id::TorsionID > torsions,
	core::kinematics::MoveMapCOP mm ) const
{
	if ( ! torsions.size() ) { return false; }
	for ( auto torsion : torsions ) {
		if ( ! mm->get( torsion ) ) { return false; }
	}
	return true;
}

std::ostream &
operator<< ( std::ostream & os, BackboneMover const & mover )
{
	mover.show( os );
	return os;
}


// SmallMover //////////////////////////////////////////////////////////////////

//constructor
SmallMover::SmallMover() : BackboneMover()
{
	moves::Mover::type( "SmallMover" );  // Why is the type different for the different constructors?  ~Labonte
}

SmallMover::SmallMover(
	core::kinematics::MoveMapOP movemap_in,
	core::Real temperature_in,
	core::Size nmoves_in ) :
	BackboneMover( movemap_in, temperature_in, nmoves_in )
{
	moves::Mover::type( "Small" );
}

// Destructor
SmallMover::~SmallMover() = default;


protocols::moves::MoverOP
SmallMover::clone() const
{
	return protocols::moves::MoverOP( new SmallMover( *this ) );
}

protocols::moves::MoverOP
SmallMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new SmallMover() );
}


void
SmallMover::setup_list( core::pose::Pose & pose )
{
	using namespace core;

	// Mask by ResidueSelector.
	select::residue_selector::ResidueSubset const subset( compute_selected_residues( pose ) );

	for ( core::uint i = 1; i <= pose.size(); ++i ) {
		if ( ! subset[ i ] ) { continue; }  // Skip residues masked by the ResidueSelector.

		// Check the maximum angle deviation for the specific kind of secondary structure and, if it is greater
		// than 0, push back the residue, its moving torsions, and the maximum angle deviation to the list of
		// movable residues.
		char ss( pose.secstruct( i ) );
		if ( ! angle_max_.count( ss ) ) { ss = 'L'; }  // Can this ever happen? ~Labonte
		Real const mx( angle_max_.find( ss )->second );

		if ( mx > 0.0 ) {
			// Check if the residue has free main-chain torsions as determined by the MoveMap.
			utility::vector1< id::TorsionID > const torsions( get_mainchain_TorsionIDs( pose.conformation(), i ) );
			if ( are_torsions_allowed( torsions, movemap( pose ) ) ) {
				move_pos_list_.push_back( std::make_tuple( i, torsions, mx ) );
			}
		}
	}
}

void
SmallMover::set_angles( core::Real angle_in ) {
	// If anyone can think of how to explain this better, please do so!!
	big_angle_ = angle_in;  // <- This number comes from max angle, it is the maximum total deviation.
	small_angle_ = big_angle_ / 2.0;  // <- This is max_angle/2, which is the deviation from the angle input.
}

bool
SmallMover::move_with_scorefxn( core::pose::Pose & pose )
{
	debug_assert( scorefxn() );
	old_rama_score_ = ( *scorefxn() )( pose );

	utility::vector1< std::pair< core::id::TorsionID, Real > > old_angles;
	for ( auto torsion : moving_torsions_ ) {
		Real const old_angle( pose.torsion( torsion ) );
		old_angles.push_back( std::make_pair( torsion, old_angle ) );
		Real const new_angle(
			basic::periodic_range( old_angle - small_angle_ + numeric::random::rg().uniform() * big_angle_, 360.0 ) );
		pose.set_torsion( torsion, new_angle );
	}

	pose.energies().clear();
	new_rama_score_ = ( *scorefxn() )( pose );
	TR.Debug << "Using score function to evaluate ramas. Old=" << old_rama_score_ << " New=" << new_rama_score_
		<< " Temp=" << temperature() << std::endl;

	if ( check_rama() ) { return true; }

	for ( auto torsion : old_angles ) {
		pose.set_torsion( torsion.first, torsion.second );
	}

	pose.energies().clear();
	TR.Debug << "Reject: reverted score = " << ( *scorefxn() )( pose ) << std::endl;
	return false;
}

bool
SmallMover::move_with_rama( core::pose::Pose & pose )
{
	utility::vector1< std::pair< core::id::TorsionID, Real > > old_angles;
	utility::vector1< std::pair< core::id::TorsionID, Real > > new_angles;
	for ( auto torsion : moving_torsions_ ) {
		Real const old_angle( pose.torsion( torsion ) );
		old_angles.push_back( std::make_pair( torsion, old_angle ) );
		Real const new_angle(
			basic::periodic_range( old_angle - small_angle_ + numeric::random::rg().uniform() * big_angle_, 360.0 ) );
		new_angles.push_back( std::make_pair( torsion, new_angle ) );
	}

	conformation::Residue const & current_rsd( pose.residue( resnum_ ) );
	// Always accept carbohydrate moves for now....
	if ( current_rsd.is_protein() && ! current_rsd.type().is_oligourea() && current_rsd.aa() != chemical::aa_unk ) {
		scoring::Ramachandran const & rama( scoring::ScoringManager::get_instance()->get_Ramachandran() );

		debug_assert( new_angles[ 1 ].first.torsion() == core::id::phi_torsion );
		debug_assert( new_angles[ 2 ].first.torsion() == core::id::psi_torsion );
		old_rama_score_ =
			rama.eval_rama_score_residue( current_rsd.aa(), old_angles[ 1 ].second, old_angles[ 2 ].second );
		new_rama_score_ =
			rama.eval_rama_score_residue( current_rsd.aa(), new_angles[ 1 ].second, new_angles[ 2 ].second );

		// Decide whether to accept the move.
		if ( !check_rama() ) { return false; }
	}

	// Set the new values for residue resnum.
	for ( auto torsion : new_angles ) {
		pose.set_torsion( torsion.first, torsion.second );
	}

	return true;
}

utility::vector1< core::id::TorsionID_Range >
SmallMover::torsion_id_ranges( core::pose::Pose & /* pose */ )
{
	return utility::vector1< core::id::TorsionID_Range >();
}

utility::vector1< core::id::DOF_ID_Range >
SmallMover::dof_id_ranges( core::pose::Pose & pose )
{
	using namespace core::id;

	Real static const pi( numeric::NumericTraits< Real >::pi() );
	clear();
	setup_list( pose );

	utility::vector1< DOF_ID_Range > range_vector;

	for ( core::Size i = 1; i <= move_pos_list_.size(); ++i ) {
		TorsionID phi_id( std::get< 0 >( move_pos_list_[ i ] ), BB, phi_torsion );
		DOF_ID phi_dof( pose.conformation().dof_id_from_torsion_id( phi_id ) );
		if ( phi_dof.valid() ) {
			range_vector.push_back( DOF_ID_Range( phi_dof, -pi, pi ) );
		}

		TorsionID psi_id( std::get< 0 >( move_pos_list_[ i ] ), BB, psi_torsion );
		DOF_ID psi_dof( pose.conformation().dof_id_from_torsion_id( psi_id ) );
		if ( psi_dof.valid() ) {
			range_vector.push_back( DOF_ID_Range( psi_dof, -pi, pi ) );
		}
	}

	return range_vector;
}


std::string
SmallMover::get_name() const
{
	return mover_name();
}

std::string
SmallMover::mover_name()
{
	return "Small";
}

void
SmallMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaComplexTypeGeneratorOP ct_gen = BackboneMover::complex_type_generator_for_backbone_mover( xsd );
	ct_gen->element_name( mover_name() )
		.complex_type_naming_func( & moves::complex_type_name_for_mover )
		.description( "Small-move style backbone-torsion moves that, "
		"unlike shear, do not minimize downstream propagation." )
		.write_complex_type_to_schema( xsd );
}


// SmallMoverCreator ///////////////////////////////////////////////////////////

protocols::moves::MoverOP
SmallMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new SmallMover );
}

std::string
SmallMoverCreator::keyname() const
{
	return SmallMover::mover_name();
}

void
SmallMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SmallMover::provide_xml_schema( xsd );
}


// ShearMover //////////////////////////////////////////////////////////////////

// Constructor
ShearMover::ShearMover() : BackboneMover()
{
	protocols::moves::Mover::type( "ShearMover" );  // Why is the type different for the different constructors?  ~Labonte
}


ShearMover::ShearMover(
	core::kinematics::MoveMapOP movemap_in,
	core::Real temperature_in,
	core::Size nmoves_in ) :
	BackboneMover( movemap_in, temperature_in, nmoves_in )
{
	protocols::moves::Mover::type( "Shear" );
}

// Destructor
ShearMover::~ShearMover() = default;


protocols::moves::MoverOP
ShearMover::clone() const
{
	return protocols::moves::MoverOP( new protocols::simple_moves::ShearMover( *this ) );
}

protocols::moves::MoverOP
ShearMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new protocols::simple_moves::ShearMover() );
}

void
ShearMover::setup_list( core::pose::Pose & pose )
{
	// Compare code below to that for SmallMover above.

	using namespace std;
	using namespace utility;
	using namespace id;
	using namespace core::conformation;

	// Mask by ResidueSelector.
	core::select::residue_selector::ResidueSubset const subset( compute_selected_residues( pose ) );

	Size const n_rsd( pose.size() );
	for ( core::uint i = 2; i <= n_rsd; ++i ) {
		if ( ! subset[ i ] ) { continue; }

		char ss( pose.secstruct( i ) );
		if ( ! angle_max_.count( ss ) ) { ss = 'L'; }
		Real const mx( angle_max_.find( ss )->second );

		if ( mx > 0.0 ) {
			Conformation const & conf( pose.conformation() );
			Residue const & rsd( pose.residue( i ) );
			vector1< TorsionID > move_torsions, torsions_to_search;

			if ( ! rsd.is_carbohydrate() ) {
				// Technically, this block of code should handle any linear biopolymer, provided that appropriate code
				// is added to BackboneMover::get_mainchain_TorsionIDs().

				// First, get the torsions for the "move" from this residue.
				move_torsions = get_mainchain_TorsionIDs( conf, i );

				// Next, we need to generate a list of potential torsions for the "counter move".
				torsions_to_search = move_torsions;
				if ( ( ! rsd.is_lower_terminus() ) && ( subset[ i - 1 ] ) ) {
					torsions_to_search.append( get_mainchain_TorsionIDs( conf, i - 1 ) );  // Add parent's torsions.
				}
				if ( ( ! rsd.is_upper_terminus() ) && ( i < n_rsd ) && ( subset[ i + 1 ] ) ) {
					torsions_to_search.append( get_mainchain_TorsionIDs( conf, i + 1 ) );  // Add child's torsions.
				}

				if ( torsions_to_search.size() < 3 ) {
					// There are only two movable torsions, and they are adjacent and cannot possibly be near parallel.
					continue;
				}
			} else /*rsd is carbohydrate*/ {
				// Because of branching, carbohydrates require additional setup.
				// If other branching systems are added to Rosetta, a more general method might be feasible.
				tie( move_torsions, torsions_to_search ) = setup_list_for_saccharide_residue( pose, i );
			}

			for ( auto move_torsion : move_torsions ) {
				TorsionID const counter_move_torsion(
					find_bond_torsion_with_nearest_orientation( conf, torsions_to_search, move_torsion ) );
				vector1< TorsionID > const torsions = { move_torsion, counter_move_torsion };
				if ( are_torsions_allowed( torsions, movemap( pose ) ) ) {
					move_pos_list_.push_back( make_tuple( i, torsions, mx ) );
				}
			}
		}
	}
}

// db (from rosetta++) set maximum deviation to be twice that of small move since shears are less perturbing
void
ShearMover::set_angles( core::Real angle_in )
{
	// If anyone can think of how to explain this, please do so!!
	big_angle_ = angle_in * 2.0;
	small_angle_ = big_angle_ / 2.0;
}

bool
ShearMover::move_with_scorefxn( core::pose::Pose & pose )
{
	debug_assert( moving_torsions_.size() == 2 );
	old_rama_score_ = ( *scorefxn() )( pose );

	Real const old_angle1( pose.torsion( moving_torsions_[ 1 ] ) );
	Real const old_angle2( pose.torsion( moving_torsions_[ 2 ] ) );

	Real const shear_delta = small_angle_ - numeric::random::rg().uniform() * big_angle_;
	Real const new_angle1( basic::periodic_range( old_angle1 - shear_delta, 360.0 ) );
	Real const new_angle2( basic::periodic_range( old_angle2 + shear_delta, 360.0 ) );

	pose.set_torsion( moving_torsions_[ 1 ], new_angle1 );
	pose.set_torsion( moving_torsions_[ 2 ], new_angle2 );

	pose.energies().clear();
	new_rama_score_ = ( *scorefxn() )( pose );

	if ( check_rama() ) { return true; }

	pose.set_torsion( moving_torsions_[ 1 ], old_angle1 );
	pose.set_torsion( moving_torsions_[ 2 ], old_angle2 );

	pose.energies().clear();
	return false;
}

bool
ShearMover::move_with_rama( core::pose::Pose & pose )
{
	debug_assert( moving_torsions_.size() == 2 );

	Real const old_angle1( pose.torsion( moving_torsions_[ 1 ] ) );
	Real const old_angle2( pose.torsion( moving_torsions_[ 2 ] ) );

	Real const shear_delta = small_angle_ - numeric::random::rg().uniform() * big_angle_;
	Real const new_angle1( basic::periodic_range( old_angle1 - shear_delta, 360.0 ) );
	Real const new_angle2( basic::periodic_range( old_angle2 + shear_delta, 360.0 ) );

	// Grab the main residue.
	conformation::Residue const & current_rsd( pose.residue( resnum_ ) );
	// Always accept carbohydrate moves for now....
	if ( current_rsd.is_protein() ) {
		conformation::Residue const & prev_rsd( pose.residue( resnum_ - 1 ) );
		scoring::Ramachandran const & rama( scoring::ScoringManager::get_instance()->get_Ramachandran() );
		if ( ! current_rsd.type().is_oligourea() && current_rsd.aa() != chemical::aa_unk ) {
			// rama for phi of resnum and psi of resnum-1
			old_rama_score_ = rama.eval_rama_score_residue( current_rsd.aa(), old_angle1, pose.psi( resnum_ ) );
			new_rama_score_ = rama.eval_rama_score_residue( current_rsd.aa(), new_angle1, pose.psi( resnum_ ) );

			// Decide whether to accept the move.
			if ( ! check_rama() ) { return false; }
		}
		if ( prev_rsd.aa() != chemical::aa_unk && ! prev_rsd.type().is_oligourea() ) {
			// rama for residue resnum-1
			old_rama_score_ = rama.eval_rama_score_residue( prev_rsd.aa(), pose.phi( resnum_ - 1 ), old_angle2 );
			new_rama_score_ = rama.eval_rama_score_residue( prev_rsd.aa(), pose.phi( resnum_ - 1 ), new_angle2 );

			// Decide whether to accept the move.
			if ( ! check_rama() ) { return false; }
		}
	}

	// Set the new values.
	pose.set_torsion( moving_torsions_[ 1 ], new_angle1 );
	pose.set_torsion( moving_torsions_[ 2 ], new_angle2 );

	return true;
}

utility::vector1< core::id::TorsionID_Range >
ShearMover::torsion_id_ranges( core::pose::Pose & /* pose */ )
{
	return utility::vector1< core::id::TorsionID_Range >();
}

utility::vector1< core::id::DOF_ID_Range >
ShearMover::dof_id_ranges( core::pose::Pose & pose )
{
	using namespace core::id;

	Real static const pi( numeric::NumericTraits< Real >::pi() );

	clear();
	setup_list( pose );

	utility::vector1< DOF_ID_Range > range_vector;

	for ( core::Size i = 1; i <= move_pos_list_.size(); ++i ) {
		TorsionID phi_id( std::get< 0 >( move_pos_list_[ i ] ), BB, phi_torsion );
		DOF_ID phi_dof( pose.conformation().dof_id_from_torsion_id( phi_id ) );
		if ( phi_dof.valid() ) {
			range_vector.push_back( core::id::DOF_ID_Range( phi_dof, -pi, pi ) );
		}

		TorsionID psi_id( std::get< 0 >( move_pos_list_[ i ] ) - 1, BB, psi_torsion );
		DOF_ID psi_dof( pose.conformation().dof_id_from_torsion_id( psi_id ) );
		if ( phi_dof.valid() ) {
			range_vector.push_back( core::id::DOF_ID_Range( psi_dof, -pi, pi ) );
		}
	}

	return range_vector;
}

std::string
ShearMover::get_name() const
{
	return mover_name();
}

std::string
ShearMover::mover_name()
{
	return "Shear";
}

void
ShearMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaComplexTypeGeneratorOP ct_gen = BackboneMover::complex_type_generator_for_backbone_mover( xsd );
	ct_gen->element_name( mover_name() )
		.complex_type_naming_func( & moves::complex_type_name_for_mover )
		.description( "Shear style backbone-torsion moves that minimize downstream propagation" )
		.write_complex_type_to_schema( xsd );
}

std::pair< utility::vector1< core::id::TorsionID >, utility::vector1< core::id::TorsionID > >
ShearMover::setup_list_for_saccharide_residue( core::pose::Pose & pose, core::uint seq_pos )
{
	using namespace std;
	using namespace utility;
	using namespace id;
	using namespace core::conformation;

	Residue const & rsd( pose.residue( seq_pos ) );

	// First, get the torsions for the "move" from this residue.
	Conformation const & conf( pose.conformation() );
	vector1< TorsionID > move_torsions( carbohydrates::get_glycosidic_TorsionIDs( conf, seq_pos ) );

	// Next, we need to generate a list of potential torsions for the "counter move".
	// To do this, we must walk through the residues of the glycan and expand a list of torsions to search
	// to include nearby torsions that will NOT cause effects across branches.
	vector1< TorsionID > torsions_to_search( move_torsions );

	// We will always need to check the parent, assuming it exists.
	core::uint parent_pos( pose.glycan_tree_set()->get_parent( seq_pos ) );
	if ( parent_pos ) {
		Residue const & parent( pose.residue( parent_pos ) );
		if ( ! parent.is_branch_point() ) {
			// (If the parent IS a branch point, we cannot use any of its bonds for the counter move, because
			// it will move the other branch(es).  It does not matter whether the parent is another
			// saccharide or an AA in such a case.)

			// The parent is NOT a branch point, so any of its main-chain torsions could be used for a
			// counter move.
			if ( parent.is_carbohydrate() ) {  // TODO: Handle other cases, such as glycolipids.
				vector1< TorsionID > const parent_torsions(
					carbohydrates::get_glycosidic_TorsionIDs( conf, parent_pos ) );
				torsions_to_search.append( parent_torsions );
			}
		}
	}

	// Now deal with any children.
	if ( ! rsd.is_upper_terminus() ) {
		// (If the residue is an upper terminus, it has no children.)

		if ( ! rsd.is_branch_point() ) {
			// (If the residue IS a branch point, there is no torsion on any child that can act as a
			// counter move.)

			// The residue is NOT a branch point, so it only has one child, whose main-chain torsions
			// could be used for a counter move.
			core::uint child_pos( pose.glycan_tree_set()->get_mainchain_child( seq_pos ) );
			//core::uint child_pos( pose.glycan_tree_set()->get_node( seq_pos )->get_mainchain_child() );
			Residue const & child( pose.residue( child_pos ) );
			if ( child.is_carbohydrate() ) {
				vector1< TorsionID > const child_torsions(
					carbohydrates::get_glycosidic_TorsionIDs( conf, child_pos ) );
				torsions_to_search.append( child_torsions );
			}
		}
	}

	if ( torsions_to_search.size() < 3 ) {
		debug_assert( move_torsions.size() == torsions_to_search.size() );
		// There are only two movable torsions, and they are adjacent and cannot possibly be near parallel.
		return make_pair( vector1< TorsionID >(), vector1< TorsionID >() );  // empty vectors
	} else if ( torsions_to_search.size() == 3 ) {
		debug_assert( move_torsions.size() == torsions_to_search.size() );
		// Since there are always at least 2 glycosidic bonds, if there are exactly 3 movable torsions,
		// they must be from a single residue with a phi, psi, and omega.  A move of psi cannot be
		// countered by either phi or omega, since they are adjacent, but phi and omega are likely near-
		// parallel, because of a tendency of psi to maintain an anti rotamer.

		// Remove the 2nd element, which must be psi.
		// (This will force phi to pair with omega and vice versa below.)
		move_torsions.erase( move_torsions.begin() + 1 );
		torsions_to_search.erase( torsions_to_search.begin() + 1 );
	}

	return make_pair( move_torsions, torsions_to_search );
}


// ShearMoverCreator ///////////////////////////////////////////////////////////

protocols::moves::MoverOP
ShearMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new ShearMover );
}

std::string
ShearMoverCreator::keyname() const
{
	return ShearMover::mover_name();
}

void
ShearMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ShearMover::provide_xml_schema( xsd );
}

}  // namespace simple_moves
}  // namespace protocols
