// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/loops/loop_closure/ccd/CCDLoopClosureMover.cc
/// @brief   Method definitions for CCDLoopClosureMover
/// @author  Phil Bradley
/// @author  Oliver Lange
/// @author  Brian Weitzner
/// @author  Labonte <JWLabonte@jhu.edu>
/// @note    This file is the result of a refactor of code written by Phil and later wrapped in a Mover by Oliver.

// Unit Headers
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMoverCreator.hh>

// Package Headers
#include <protocols/loops/loop_closure/ccd/RamaCheck.hh>

// Project Headers
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/assert.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/py/PyAssert.hh>
#include <utility/tag/Tag.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/vector1.hh>

// Numeric Headers
#include <numeric/cyclic_coordinate_descent.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

// Basic header
#include <basic/prof.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.loops.loop_closure.ccd.CCDLoopClosureMover" );

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

using namespace core;

// Define constants.
// Indexed by ChainDirection
utility::vector1< SSize > const CCDLoopClosureMover::STEP_SIZE = utility::tools::make_vector1( 1, -1 );
Distance const CCDLoopClosureMover::MALARKEY = -1.0;
Real const CCDLoopClosureMover::BAD_SCORE = 500.0;

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Empty constructor
/// @details  By default, all residues within the given Loop will be allowed to move.
CCDLoopClosureMover::CCDLoopClosureMover() : Mover()
{
	using kinematics::MoveMap;
	using kinematics::MoveMapOP;

	Loop empty_loop;

	// Set default MoveMap.
	MoveMapOP default_movemap( new MoveMap() );
	default_movemap->set_bb( true );

	init( empty_loop, default_movemap );
}

// Copy constructor
CCDLoopClosureMover::CCDLoopClosureMover( CCDLoopClosureMover const & object_to_copy ) : Mover( object_to_copy )
{
	copy_data( *this, object_to_copy );
}

// Constructor with Loop input option
/// @details  By default, all residues within the given Loop will be allowed to move.
CCDLoopClosureMover::CCDLoopClosureMover( protocols::loops::Loop const & loop ) : Mover()
{
	using kinematics::MoveMap;
	using kinematics::MoveMapOP;

	// Set default MoveMap.
	MoveMapOP default_movemap( new MoveMap() );
	for ( core::uint resnum = loop.start(); resnum <= loop.stop(); ++resnum ) {
		default_movemap->set_bb( resnum, true );
	}

	init( loop, default_movemap );
}

// Constructor with Loop and MoveMap input options
CCDLoopClosureMover::CCDLoopClosureMover( protocols::loops::Loop const & loop, kinematics::MoveMapCOP mm ) :
	Mover()
{
	init( loop, mm );
}

// Destructor
CCDLoopClosureMover::~CCDLoopClosureMover() {}

// Assignment operator
CCDLoopClosureMover &
CCDLoopClosureMover::operator=( CCDLoopClosureMover const & object_to_copy )
{
	// Abort self-assignment.
	if ( this == &object_to_copy ) {
		return *this;
	}

	Mover::operator=( object_to_copy );
	copy_data( *this, object_to_copy );
	return *this;
}


// Standard Rosetta methods ///////////////////////////////////////////////////
// General methods
void
CCDLoopClosureMover::register_options()
{
	using namespace basic::options;
	option.add_relevant( OptionKeys::loops::ccd::max_torsion_delta_per_move );
	option.add_relevant( OptionKeys::loops::ccd::max_torsion_delta );
	option.add_relevant( OptionKeys::loops::ccd::tolerance );
	option.add_relevant( OptionKeys::loops::ccd::max_cycles );
	option.add_relevant( OptionKeys::loops::ccd::check_rama_scores );

	// call register_options() on all other Movers used by this class
	Mover::register_options(); // Mover's register_options doesn't do anything, it's really just here as an example.
	RamaCheckBase::register_options();
}

/// @details Print the CcdLoopMover's basic information, the loop it is working on, and details about the instance's
/// configuration.
void
CCDLoopClosureMover::show( std::ostream & output ) const
{
	using std::endl;

	Mover::show( output );  // name, type, tag
	loop_.show( output );
	output << "\nNumber of CCD cycles:    " << max_cycles_ <<
		"\nTolerance:               " << tolerance_ <<
		"\nbRama check:             " << ( check_rama_scores_ ? "True" : "False" ) <<
		"\nMax total delta helix:   " << max_total_torsion_delta_per_residue( helix ) <<
		"\nMax total delta strand:  " << max_total_torsion_delta_per_residue( strand ) <<
		"\nMax total delta loop:    " << max_total_torsion_delta_per_residue( coil ) << endl;
	output << "Movemap: " << endl;
	movemap_->show( output );

	if ( check_rama_scores() ) {
		output << "RamaCheck: " << endl;
		rama()->show( output );
	}
}

// Mover methods
std::string
CCDLoopClosureMover::get_name() const
{
	return type();
}

protocols::moves::MoverOP
CCDLoopClosureMover::clone() const
{
	return protocols::moves::MoverOP( new CCDLoopClosureMover( *this ) );
}

protocols::moves::MoverOP
CCDLoopClosureMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new CCDLoopClosureMover() );
}

void
CCDLoopClosureMover::apply( pose::Pose & pose )
{
	if ( verbose_ ) TR << "Closing loop " << loop_.start() << " " << loop_.stop() << std::endl;
	PROF_START( basic::CCD_CLOSE );

	using std::endl;
	using utility::vector1;

	// Reset status variables.
	average_change_in_torsion_angle_ = 0.0;
	average_change_in_rama_score_ = 0.0;

	// Save the starting pose for later comparison.
	pose::Pose starting_pose( pose );

	// Save the starting rama scores if check_rama_scores_ is true.
	if ( check_rama_scores() ) {
		rama()->initialize_starting_rama_scores( pose );
	}

	actual_cycles_ = 0;  // for reporting how many cycles we actually took

	// Now start cycling through.
	for ( core::uint cycle = 1; cycle <= max_cycles_; ++cycle ) {
		++actual_cycles_;

		// First forward, then backward.
		close_loop_in_single_direction ( pose, starting_pose, forward );
		close_loop_in_single_direction ( pose, starting_pose, backward );

		// Check overlap deviations to see if loop is closed every 5 cycles or on the last one.
		if ( ( cycle % 5 ) == 0 || cycle == max_cycles_ ) {
			// Calculate the RMSD of the virtual atoms to the atoms with which they are supposed to overlap

			if ( TR.Debug.visible() ) {
				TR.Debug << "cycle / deviation: " << cycle << " " << deviation_ << endl;
			}

			if ( success() ) {
				if ( TR.Debug.visible() ) {
					TR.Debug << "closed early: cycle= " << cycle << " " << tolerance_ << endl;
				}
				break;
			}
		}
	}  // next cycle

	compute_closure_metrics( pose, starting_pose );
	PROF_STOP( basic::CCD_CLOSE );
}

void CCDLoopClosureMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	Pose const & /*pose*/ )
{
	// max_torsion_delta_per_move: This is split into three tags.
	if ( tag->hasOption( "max_torsion_delta_per_move_H" ) ) {
		max_per_move_torsion_delta_[ helix ] = tag->getOption< Real >( "max_torsion_delta_per_move_H" );
	}
	if ( tag->hasOption( "max_torsion_delta_per_move_E" ) ) {
		max_per_move_torsion_delta_[ strand ] = tag->getOption< Real >( "max_torsion_delta_per_move_E" );
	}
	if ( tag->hasOption( "max_torsion_delta_per_move_L" ) ) {
		max_per_move_torsion_delta_[ coil ] = tag->getOption< Real >( "max_torsion_delta_per_move_L" );
	}

	// max_torsion_delta: This is split into three tags.
	if ( tag->hasOption( "max_torsion_delta_H" ) ) {
		max_total_torsion_delta_[ helix ] = tag->getOption< Real >( "max_torsion_delta_H" );
	}
	if ( tag->hasOption( "max_torsion_delta_E" ) ) {
		max_total_torsion_delta_[ strand ] = tag->getOption< Real >( "max_torsion_delta_E" );
	}
	if ( tag->hasOption( "max_torsion_delta_L" ) ) {
		max_total_torsion_delta_[ coil ] = tag->getOption< Real >( "max_torsion_delta_L" );
	}

	if ( tag->hasOption( "tolerance" ) ) {
		tolerance( tag->getOption< Real >( "tolerance" ) );
	}

	if ( tag->hasOption( "max_cycles" ) ) {
		max_cycles( tag->getOption< core::uint >( "max_cycles" ) );
	}

	if ( tag->hasOption( "check_rama_scores" ) ) {
		check_rama_scores( tag->getOption< bool >( "check_rama_scores" ) );
	}

	if ( tag->hasOption( "rama_2b" ) ) {
		use_rama_2B( tag->getOption< bool >( "rama_2b" ) );
	}

	// Have RamaCheck parse the tag if it has been enabled
	if ( check_rama_scores() ) {
		rama()->parse_my_tag( tag );
	}
}

// Accessors/Mutators /////////////////////////////////////////////////////////
//Get the Loop to be closed.
protocols::loops::Loop
CCDLoopClosureMover::loop() const
{
	return loop_;
}

// Set the Loop to be closed.
void
CCDLoopClosureMover::loop( protocols::loops::Loop new_loop )
{
	loop_ = new_loop;
}

// Get the current MoveMap.
kinematics::MoveMap
CCDLoopClosureMover::movemap() const {
	return *movemap_;
}

// Set the MoveMap.
void
CCDLoopClosureMover::movemap( kinematics::MoveMapCOP new_movemap )
{
	movemap_ = new_movemap;
}

/// @param   <sectruct>: 'H', 'E', or 'L'
/// @return  angle in degrees
Angle
CCDLoopClosureMover::max_per_move_torsion_delta_per_residue( char secstruct ) const
{
	using utility::excn::EXCN_KeyError;

	if ( ! sec_struc_char_to_enum_map().count(  secstruct ) ) {
		std::stringstream msg;
		msg << "CCDLoopClosureMover::max_per_move_delta_per_residue( char secstruct ): secstruct must be 'H', 'E', or 'L'. '";
		msg << secstruct << "' is not valid.";

		// Can python catch exceptions thrown from C++? If not, we had better keep this PyAssert!
		PyAssert( ( secstruct == 'H' ) || ( secstruct == 'E' ) || ( secstruct == 'L' ), msg.str() );
		throw EXCN_KeyError( msg.str() );
	}

	return max_per_move_torsion_delta_per_residue( sec_struc_char_to_enum_map().find( secstruct )->second );
}

/// @param   <sectruct>: 'H', 'E', or 'L'
/// @return  angle in degrees
Angle
CCDLoopClosureMover::max_total_torsion_delta_per_residue( char secstruct ) const
{
	using utility::excn::EXCN_KeyError;

	if ( ! sec_struc_char_to_enum_map().count(  secstruct ) ) {
		std::stringstream msg;
		msg << "CCDLoopClosureMover::max_total_delta_per_residue( char secstruct ): secstruct must be 'H', 'E', or 'L'. '";
		msg << secstruct << "' is not valid.";

		// Can python catch exceptions thrown from C++? If not, we had better keep this PyAssert!
		PyAssert( ( secstruct == 'H' ) || ( secstruct == 'E' ) || ( secstruct == 'L' ), msg.str() );
		throw EXCN_KeyError( msg.str() );
	}
	return max_total_torsion_delta_per_residue( sec_struc_char_to_enum_map().find( secstruct )->second );
}

// Other Public Methods ///////////////////////////////////////////////////////
bool
CCDLoopClosureMover::success() const
{
	return ( ( deviation() < tolerance_ ) );
}

// Private methods ////////////////////////////////////////////////////////////
// Initialize data members from arguments.
void
CCDLoopClosureMover::init( protocols::loops::Loop const & loop, kinematics::MoveMapCOP movemap )
{
	type("CCDLoopClosureMover");
	loop_ = loop;
	movemap_ = movemap;
	verbose_ = false;

	// These defaults are copied from Rosetta++: map_squence.cc::scored_frag_close()
	max_per_move_torsion_delta_.resize( n_secondary_structure_types );
	max_per_move_torsion_delta_[ helix ] = 1.0;
	max_per_move_torsion_delta_[ strand ] = 5.0;
	max_per_move_torsion_delta_[ coil ] = 10.0;
	max_total_torsion_delta_.resize( n_secondary_structure_types );
	max_total_torsion_delta_[ helix ] = 10.0;
	max_total_torsion_delta_[ strand ] = 50.0;
	max_total_torsion_delta_[ coil ] = 75.0;
	tolerance_ = 0.08;

	max_cycles_ = 100;
	check_rama_scores_ = true;
	use_rama_2b_ = false; // not tested, but probably a good move.

	deviation_ = MALARKEY;

	rama_ = NULL;

	init_options();
}

void CCDLoopClosureMover::init_options()
{
	// I like all of my options reading to happen at initialization (and then never again!).
	// For clarity, I only set options if the user has it on the command line so the defaults
	// are what you see in code, specifically in the init() method.
	using namespace basic::options;

	if ( option[ OptionKeys::loops::ccd::max_torsion_delta_per_move ].user() ) {
		utility::vector1< Angle > const opt = option[ OptionKeys::loops::ccd::max_torsion_delta_per_move ]();
		assert( opt.size() == 3 ); // this should be taken care of by the option system
		max_per_move_torsion_delta_per_residue( opt[ 1 ], opt[ 2 ], opt[ 3 ] );
	}

	if ( option[ OptionKeys::loops::ccd::max_torsion_delta ].user() ) {
		utility::vector1< Angle > const opt = option[ OptionKeys::loops::ccd::max_torsion_delta ]();
		assert( opt.size() == 3 ); // this should be taken care of by the option system
		max_total_torsion_delta_per_residue( opt[ 1 ], opt[ 2 ], opt[ 3 ] );
	}

	if ( option[ OptionKeys::loops::ccd::tolerance ].user() ) {
		tolerance( option[ OptionKeys::loops::ccd::tolerance ]() );
	}

	if ( option[ OptionKeys::loops::ccd::max_cycles ].user() ) {
		max_cycles( option[ OptionKeys::loops::ccd::max_cycles ]() );
	}

	if ( option[ OptionKeys::loops::ccd::check_rama_scores ].user() ) {
		check_rama_scores( option[ OptionKeys::loops::ccd::check_rama_scores ]() );
	}

	if ( option[ OptionKeys::loops::ccd::rama_2b ].user() ) {
		use_rama_2B( option[ OptionKeys::loops::ccd::rama_2b ]() );
	}
}

// Copy all data members from <from> to <to>.
void
CCDLoopClosureMover::copy_data( CCDLoopClosureMover & to, CCDLoopClosureMover const & from ) const
{
	to.loop_ = from.loop_;
	to.movemap_ = core::kinematics::MoveMapCOP( core::kinematics::MoveMapOP( new kinematics::MoveMap( * from.movemap_ ) ) );  // deep copy pointer

	to.max_per_move_torsion_delta_ = from.max_per_move_torsion_delta_;
	to.max_total_torsion_delta_ = from.max_total_torsion_delta_;
	to.tolerance_ = from.tolerance_;

	to.max_cycles_ = from.max_cycles_;
	to.check_rama_scores_ = from.check_rama_scores_;
	to.use_rama_2b_ = from.use_rama_2b_;

	to.deviation_ = from.deviation_;
	to.average_change_in_torsion_angle_ = from.average_change_in_torsion_angle_;

	to.rama_ = NULL; // correct rama_ is initialized lazily based on check_rama_scores_ and use_rama_2b_.

	to.average_change_in_rama_score_ = from.average_change_in_rama_score_;

	to.actual_cycles_ = from.actual_cycles_;
}

// Return the coordinates of the two or three atoms to be overlapped for a given residue.
utility::vector1< PointPosition >
CCDLoopClosureMover::get_anchors( conformation::Residue const & residue ) const
{
	using utility::vector1;

	vector1< PointPosition > anchors;

	// TODO: Edit Residue and ResidueType to avoid string look-ups.
	if ( residue.has_variant_type( chemical::CUTPOINT_LOWER ) ) {
		core::uint const last_lower_mainchain_atom( residue.mainchain_atom( residue.mainchain_atoms().size() ) );

		anchors.push_back( residue.atom( last_lower_mainchain_atom ).xyz() );
		anchors.push_back( residue.atom( "OVL1" ).xyz() );
		// If the residue does not have an OVL2, it means that OVL2 is not necessary for proper loop closure, i.e., the
		// bond being formed is freely rotating, i.e. does not have partial-double or greater bond order.

		// When comparing M and F later, use the min size of the two vectors to get the "magic number" 3 for
		// peptides and 2 for things like sugars. ~Labonte
		// Jason, I think that approach will have some problems.  I think it would be better to adjust the sizes of
		// the vectors here and rely on them being the same size everywhere else. ~BDW
		// Well, I made the change now, and all unit tests passed.  If it becomes a problem, we can try to figure out
		// another solution. ~Labonte
		if ( residue.has( "OVL2" ) ) {
			anchors.push_back( residue.atom( "OVL2" ).xyz() );
		}
	} else if ( residue.has_variant_type( chemical::CUTPOINT_UPPER ) ) {
		core::uint const first_upper_mainchain_atom( residue.mainchain_atom( 1 ) );
		core::uint const second_upper_mainchain_atom( residue.mainchain_atom( 2 ) );

		anchors.push_back( residue.atom( "OVU1" ).xyz() );
		anchors.push_back( residue.atom( first_upper_mainchain_atom ).xyz() );
		anchors.push_back( residue.atom( second_upper_mainchain_atom ).xyz() );
	} else {
		std::string const msg( "CCDLoopClosureMover::get_anchors( core::conformation::Residue const & residue ): "
			"Residue is not a cutpoint variant! You must add cutpoint variants before applying this Mover." );
#ifdef PYROSETTA
		  PyAssert( false, msg );
#endif
		throw utility::excn::EXCN_BadInput( msg );
	}

	return anchors;
}

// TODO: Is it possible to ever call this in cases where atom > ( 1  + number_mc_atoms )?
void
CCDLoopClosureMover::index_pair_in_range(
	core::uint & pos, core::uint & atom, Size const n_mainchain_atoms ) const
{
	while ( atom > n_mainchain_atoms ) {
		atom -= n_mainchain_atoms;
		pos += 1;
	}
	while ( atom < 1 ) {
		atom += n_mainchain_atoms;
		if ( pos == 0 ) {
			throw utility::excn::EXCN_RangeError( "pos is unsigned. Cannot make it negative!" );
		}
		pos -= 1;
	}
}

void
CCDLoopClosureMover::get_torsion_axis(
	pose::Pose const & pose,
	core::uint const seqpos,
	core::uint const torsion_num,
	Vector & axis_atom_coords,
	Vector & axis_unit_vector ) const
{
	Size const n_mainchain_atoms( pose.residue( seqpos ).mainchain_atoms().size() );
	core::uint const upstream_atom_resnum( seqpos ), upstream_mainchain_atom_num( torsion_num );
	core::uint downstream_atom_resnum( seqpos ), downstream_mainchain_atom_num( torsion_num + 1 );

	// Adjust atom and residue number for the last torsion, if applicable.
	index_pair_in_range( downstream_atom_resnum, downstream_mainchain_atom_num, n_mainchain_atoms );

	// TODO: Pull from the list of main-chain atoms, in case they aren't the 1st atoms in the residue.
	axis_atom_coords = pose.residue( downstream_atom_resnum ).xyz( downstream_mainchain_atom_num );
	axis_unit_vector = ( axis_atom_coords -
		pose.residue( upstream_atom_resnum ).xyz( upstream_mainchain_atom_num ) ).normalized();
}

/// @param <F>: the coordinates of the fixed target atoms
/// @param <M>: the coordinates of the moving positions to be overlapped with the target atoms
/// @param <coords>: all of the coordinates of mainchain atoms from a Pose
/// @param <seqpos>: the index of the residue for which we are calculating a torsion angle change
/// @param <torsion_num> : the index of the mainchain torsion for which we are calculating a change
/// @note  The math (including the variable names used) in this function comes from
/// Canutescu & Dunbrack (2003) Prot. Sci. 12, 963.
Angle
CCDLoopClosureMover::calculate_ccd_angle(
	pose::Pose const & pose,
	core::uint const seqpos,
	core::uint const torsion_num,
	ChainDirection const direction )
{
	using std::endl;
	using utility::vector1;
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	// Get the fixed and moving positions for this direction.
	// For the forward direction (e.g., N->C), the 1st 3 main-chain atoms of the residue after cut-point are fixed
	// For the backward direction (e.g., C->N), the last 3 main-chain atoms of the cut-point residue are fixed
	// We need to change it so that it is 2 instead of three for non-peptide cases. ~Labonte
	vector1< PointPosition > fixed_atoms, moving_atoms;
	switch ( direction ) {
	case forward :
		moving_atoms = get_anchors( pose.residue( loop_.cut() ) );
		fixed_atoms = get_anchors( pose.residue( loop_.cut() + 1 ) );
		break;
	case backward :
		moving_atoms = get_anchors( pose.residue( loop_.cut() + 1 ) );
		fixed_atoms = get_anchors( pose.residue( loop_.cut() ) );
		break;
	}
	Size number_of_atoms( numeric::min( fixed_atoms.size(), moving_atoms.size() ) );

	// Calculate the axis vector of the torsion angle.
	PointPosition axis_atom;
	Vector axis;
	get_torsion_axis( pose, seqpos, torsion_num, axis_atom, axis );

	axis *= STEP_SIZE[ direction ]; // Make sure theta_hat is pointing toward the cutpoint

	if ( TR.Debug.visible() ) {
		TR.Debug << "Initial deviation: " << deviation_ << endl;
	}

	Angle angle;
	Real dev;
	numeric::ccd_angle( fixed_atoms, moving_atoms, axis_atom, axis, angle, dev);

	// convert dev (sum of the squared deviations) to RMSD
	deviation_ = sqrt( dev / number_of_atoms );

	if ( TR.Debug.visible() ) {
		TR.Debug << "Final deviation: " << deviation_ << endl;
	}

	return angle;
}

void CCDLoopClosureMover::get_maximum_torsion_deltas_for_residue(
	core::pose::Pose const & pose,
	core::uint const seqpos,
	core::Real & per_move_allowed_delta,
	core::Real & total_allowed_delta ) const
{
	per_move_allowed_delta = max_per_move_torsion_delta_per_residue( pose.secstruct( seqpos ) );
	total_allowed_delta = max_total_torsion_delta_per_residue( pose.secstruct( seqpos ) );
}

// This method performs much of the work of the CCD closure.
// It operates on a single residue
void CCDLoopClosureMover::adjust_residue_to_minimize_deviation(
	core::pose::Pose & pose,
	core::pose::Pose const & starting_pose,
	core::uint const seqpos,
	ChainDirection const direction,
	Real const max_per_move_torsion_delta,
	Real const max_total_torsion_delta )
{
	using std::endl;
	using basic::subtract_degree_angles;
	using id::BB;
	using id::TorsionID;

	// Cycle through each movable torsion at the current residue position.
	core::conformation::Residue const & res( pose.residue( seqpos ) );
	Size const n_mainchain_torsions( res.mainchain_atoms().size() );
	for ( core::uint torsion_num = 1; torsion_num <= n_mainchain_torsions; ++torsion_num ) {

		if ( TR.Debug.visible() ) {
			TR.Debug << "Residue number: " << seqpos << " Torsion number: " << torsion_num << " Direction: "
				<< ( (direction == forward ) ? "forward" : "backward" ) << endl;
		}

		TorsionID const torsion_id( seqpos, BB, torsion_num );
		if ( ! movemap_->get( torsion_id ) ) { continue; }
		if ( is_mainchain_torsion_also_ring_torsion( res.type(), torsion_num ) ) { continue; }
		if ( seqpos == loop_.cut() && torsion_num == n_mainchain_torsions ) { continue; }  // TODO: Do we need this?

		Angle alpha = calculate_ccd_angle( pose, seqpos, torsion_num, direction );

		// Impose per-move maximum deltas.
		if ( alpha >  max_per_move_torsion_delta ) { alpha =  max_per_move_torsion_delta; }
		if ( alpha < -max_per_move_torsion_delta ) { alpha = -max_per_move_torsion_delta; }

		// Check for total movement during closure run.

		Angle const total_torsion_delta( subtract_degree_angles( starting_pose.torsion( torsion_id ),
			pose.torsion( torsion_id ) + alpha ) );

		if ( total_torsion_delta > max_total_torsion_delta ) {
			// This logic is a little tricky:
			// If adding alpha to the previous total_torsion_delta pushes us past 180 degrees from the starting
			// torsion angle, then it won't work, so check for that case.
			// (Note that if max_total_torsion_delta_ > 180 degrees, we won't ever get here.)
			assert( alpha + max_total_torsion_delta < 180 );
			if ( alpha > 0 ) {
				alpha -= ( total_torsion_delta - max_total_torsion_delta + 0.01 );
			} else {
				alpha += ( total_torsion_delta - max_total_torsion_delta + 0.01 );
			}
		}

		if ( ! check_rama_scores() || rama()->accept_new_conformation( pose, torsion_id, alpha ) ) {
			// Update the coordinates to reflect the torsion angle change.
			pose.set_torsion( torsion_id, pose.torsion( torsion_id ) + alpha );
		}
	} // next torsion_num
}

// This method performs most of the work of the CCD closure.
// It operates on the loop from the anchor to the cutpoint
void
CCDLoopClosureMover::close_loop_in_single_direction(
	pose::Pose & pose,
	pose::Pose const & starting_pose,
	ChainDirection const direction )
{
	using std::endl;

	core::uint start_pos, stop_pos;
	switch ( direction ) {
	case forward :
		start_pos = loop_.start();
		stop_pos = loop_.cut();
		break;
	case backward :
		start_pos = loop_.stop();
		stop_pos = loop_.cut() + 1;
		break;
	default :
		utility_exit_with_message("Unknown chain direction.");
		break;
	}

	// Cycle through each residue position in the loop in the proper direction.
	for ( core::uint seqpos = start_pos; seqpos * STEP_SIZE[ direction ] <= stop_pos * STEP_SIZE[ direction ];
			seqpos += STEP_SIZE[ direction ] ) {
		// Move to the next residue if the MoveMap prevents us from changing this residue.
		if ( ! movemap_->get_bb( seqpos ) ) {
			if ( TR.Debug.visible() ) {
				TR.Debug << "The MoveMap prevents movement of residue: " << seqpos << endl;
			}
			continue;
		}

		// Check maximums.
		Angle max_per_move_torsion_delta, max_total_torsion_delta;
		get_maximum_torsion_deltas_for_residue( pose, seqpos, max_per_move_torsion_delta, max_total_torsion_delta );

		// Move to the next residue if the angle maximums prevent movement of this residue.
		if ( max_total_torsion_delta <= 0.01 ) {
			if ( TR.Debug.visible() ) {
				TR.Debug << "Settings for maximum allowed changes in torsion angles prevent movement of residue: " <<
					seqpos << endl;
			}
			continue;
		}

		adjust_residue_to_minimize_deviation( pose, starting_pose, seqpos, direction, max_per_move_torsion_delta,
			max_total_torsion_delta );
	}  // next seqpos
}

void CCDLoopClosureMover::compute_closure_metrics(
	pose::Pose const & pose,
	pose::Pose const & starting_pose )
{
	using id::BB;
	using id::TorsionID;
	using basic::subtract_degree_angles;

	// Calculate torsion dev
	average_change_in_torsion_angle_ = 0.0;

	Angle total_change_in_torsion_angle( 0.0 );

	for ( core::uint i = loop_.start(); i <= loop_.stop(); ++i ) {
		Size n_mainchain_atoms( pose.residue( i ).mainchain_atoms().size() );
		for ( core::uint j = 1; j<= n_mainchain_atoms; ++j ) {
			TorsionID torsion_id( i, BB, j );
			total_change_in_torsion_angle +=
				std::abs( subtract_degree_angles( starting_pose.torsion( torsion_id ), pose.torsion( torsion_id ) ) );
		}
	}
	average_change_in_torsion_angle_ = total_change_in_torsion_angle / loop_.size();

	// Calculate change in Rama score if applicable
	average_change_in_rama_score_ = 0.0;
	if ( check_rama_scores() ) {
		average_change_in_rama_score_ = rama()->average_change_in_rama_score_over_range(
			pose, loop_.start(), loop_.stop() );
	}
}

RamaCheckBaseOP CCDLoopClosureMover::rama() const
{
	if ( rama_.get() == NULL ) {
		if ( ! use_rama_2B() ) {
			// one-body Ramachandran score
			rama_ = RamaCheckBaseOP( new RamaCheck1B );
		} else {
			// two-body Ramachandran score
			rama_ = RamaCheckBaseOP( new RamaCheck2B );
		}
	}

	return rama_;
}

std::map< char, SecondaryStructureType > const &
CCDLoopClosureMover::sec_struc_char_to_enum_map() {
	// This line will only ever be executed once.
	static std::map< char, SecondaryStructureType > * sec_struc_char_to_enum = 0;

	if ( sec_struc_char_to_enum == 0 ) {
		sec_struc_char_to_enum = new std::map< char, SecondaryStructureType >;
		( *sec_struc_char_to_enum )[ 'H' ] = helix;
		( *sec_struc_char_to_enum )[ 'E' ] = strand;
		( *sec_struc_char_to_enum )[ 'L' ] = coil;
		( *sec_struc_char_to_enum )[ 'X' ] = coil; // other
	}

	return *sec_struc_char_to_enum;
}

// Helper methods /////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that CCDLoopClosureMover can be "printed" in PyRosetta).
std::ostream &
operator<< ( std::ostream & os, CCDLoopClosureMover const & mover )
{
	mover.show( os );
	return os;
}


// Creator methods ////////////////////////////////////////////////////////////
std::string
CCDLoopClosureMoverCreator::keyname() const
{
	return CCDLoopClosureMoverCreator::mover_name();
}

protocols::moves::MoverOP
CCDLoopClosureMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new CCDLoopClosureMover );
}

std::string
CCDLoopClosureMoverCreator::mover_name()
{
	return "CCDLoopClosureMover";
}

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols
