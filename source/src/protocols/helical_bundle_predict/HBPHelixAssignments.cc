// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HBPHelixAssignments.cc
/// @brief A class for storing the helix assignments for a pose.  This can represent those proposed from an input file, or those at the current state of a trajectory.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Project headers
#include <protocols/helical_bundle_predict/HBPHelixAssignments.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/id/TorsionID.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <core/conformation/parametric/RealVectorValuedParameter.hh>
#include <core/conformation/parametric/SizeValuedParameter.hh>
#include <core/conformation/parametric/SizeVectorValuedParameter.hh>

// Protocols headers
#include <protocols/helical_bundle/util.hh>
#include <protocols/helical_bundle/BundleParametrizationCalculator.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// Utility headers
#include <utility/pointer/memory.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>

// Numeric headers
#include <numeric/angle.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/util.hh>

static basic::Tracer TR( "protocols.helical_bundle_predict.HBPHelixAssignments" );


namespace protocols {
namespace helical_bundle_predict {

////////////////////////////// HBPHelixParameters class //////////////////////////////

/// @brief Default constructor.
HBPHelixParameters::HBPHelixParameters() :
	utility::VirtualBase(),
	r0_range_( std::make_pair( 0.0, 0.0 ) ),
	omega0_range_( std::make_pair( 0.0, 0.0 ) ),
	delta_omega1_range_( std::make_pair( 0.0, 0.0 ) ),
	bundle_calculator_( nullptr )
{}

/// @brief Copy constructor.
HBPHelixParameters::HBPHelixParameters(
	HBPHelixParameters const & src
) :
	VirtualBase( src ),
	r0_range_( src.r0_range_ ),
	omega0_range_( src.omega0_range_ ),
	delta_omega1_range_( src.delta_omega1_range_ ),
	bundle_calculator_( src.bundle_calculator_ == nullptr ? nullptr : utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::BundleParametrizationCalculator >( src.bundle_calculator_->clone() ) )
{
	if ( src.bundle_calculator_ != nullptr ) {
		runtime_assert( bundle_calculator_ != nullptr );
	}
}

/// @brief Destructor.
HBPHelixParameters::~HBPHelixParameters() {}

/// @brief Create a copy of this object and return an owning pointer to the copy.
HBPHelixParametersOP
HBPHelixParameters::clone() const {
	return utility::pointer::make_shared< HBPHelixParameters >(*this);
}

/// @brief Set minimum r0 value.
void HBPHelixParameters::set_r0_min( core::Real const &setting ) { r0_range_.first = setting; }

/// @brief Set maximum r0 value.
void HBPHelixParameters::set_r0_max( core::Real const &setting ) { r0_range_.second = setting; }

/// @brief Set minimum omega0 value.
void HBPHelixParameters::set_omega0_min( core::Real const &setting ) { omega0_range_.first = setting; }

/// @brief Set maximum omega0 value.
void HBPHelixParameters::set_omega0_max( core::Real const &setting ) { omega0_range_.second = setting; }

/// @brief Set minimum delta_omega1 value.
void HBPHelixParameters::set_delta_omega1_min( core::Real const &setting ) { delta_omega1_range_.first = setting; }

/// @brief Set maximum delta_omega1 value.
void HBPHelixParameters::set_delta_omega1_max( core::Real const &setting ) { delta_omega1_range_.second = setting; }

/// @brief Set the bundle parameterization calculator.
/// @details Input is cloned.
void
HBPHelixParameters::set_calculator(
	protocols::helical_bundle::BundleParametrizationCalculatorCOP const & calculator_in
) {
	runtime_assert( calculator_in != nullptr );
	bundle_calculator_ = utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::BundleParametrizationCalculator >( calculator_in->clone() );
	runtime_assert( bundle_calculator_ != nullptr );
}

/// @brief Create a new bundle parameterization calculator, and initialize it from a Crick params file.
/// @details WARNING!  Triggers read from disk!
void
HBPHelixParameters::create_calculator_from_file(
	std::string const & filename
) {
	bundle_calculator_ = utility::pointer::make_shared< protocols::helical_bundle::BundleParametrizationCalculator >( false );
	bundle_calculator_->init_from_file( filename );
}

////////////////////////////// HBPHelix class //////////////////////////////

/// @brief Default constructor.
HBPHelix::HBPHelix() :
	utility::VirtualBase(),
	start_position_(0),
	end_position_(0),
	nucleation_prob_( 0.01 ),
	extension_prob_( 0.05 ),
	retraction_prob_( 0.03 ),
	parameters_(utility::pointer::make_shared< HBPHelixParameters >()),
	crick_params_filename_("")
{}

/// @brief Copy constructor.
/// @details Deep-copies paramters.
HBPHelix::HBPHelix(
	HBPHelix const & src
) :
	VirtualBase( src ),
	start_position_(src.start_position_),
	end_position_(src.end_position_),
	nucleation_prob_(src.nucleation_prob_),
	extension_prob_(src.extension_prob_),
	retraction_prob_(src.retraction_prob_),
	parameters_( src.parameters_ == nullptr ? nullptr : src.parameters_->clone() ),
	crick_params_filename_( src.crick_params_filename_ )
{}

/// @brief Destructor.
HBPHelix::~HBPHelix() {}

/// @brief Create a copy of this object and return an owning pointer to the copy.
HBPHelixOP
HBPHelix::clone() const {
	return utility::pointer::make_shared< HBPHelix >(*this);
}

/// @brief Get the parameters (nonconst access).
HBPHelixParametersOP
HBPHelix::parameters() {
	return parameters_;
}

/// @brief Is a sequence position within a helix?
/// @details Returns false if either of startpos or endpos is zero.
bool
HBPHelix::is_in_helix(
	core::Size const seqpos
) const {
	if ( start_position() == 0 || end_position() == 0 ) return false;
	return ( start_position() <= seqpos && seqpos <= end_position() );
}

/// @brief Get the number of residues per turn of this helix type.
core::Size
HBPHelix::residues_per_turn() const {
	debug_assert( parameters_ != nullptr );
	debug_assert( parameters_->bundle_calculator() != nullptr );
	return parameters_->bundle_calculator()->residues_per_turn();
}

/// @brief Given the start and end position of a helix, determine from the residue composition what the repeating unit offset should be.
/// @details For example, if the stretch had pattern alpha-alpha-beta-alpha, and the Crick params file expected alpha-alpha-alpha-beta, the offset
/// would be 3.
core::Size
HBPHelix::determine_repeating_unit_offset(
	core::pose::Pose const & ,//pose,
	core::Size const ,//start_position,
	core::Size const //end_position
) const {
	runtime_assert_string_msg( parameters_->bundle_calculator()->residues_per_repeat() == 1, "Error in HBPHelix::determine_repeating_unit_offset(): This function has not yet been generalized for repeating units of more than one residue." );
	return 0;
}

/// @brief Generate a map of (TorsionID->torsion value) for all mainchain torsions in the current helix.
/// @returns True for FAILURE (Crick parameters don't generate a helix), false for success.
bool
HBPHelix::get_torsions_for_helix(
	std::map< core::id::TorsionID, core::Real > & torsion_map_out,
	core::pose::Pose const & pose_for_reference
) const {
	utility::vector1< utility::vector1< numeric::xyzVector< core::Real > > > atom_positions;
	bool failed( false );
	protocols::helical_bundle::generate_atom_positions( atom_positions, pose_for_reference, start_position_, end_position_,
		parameters_->bundle_calculator()->real_parameter_cop( protocols::helical_bundle::BPC_r0 )->value(),
		parameters_->bundle_calculator()->real_parameter_cop( protocols::helical_bundle::BPC_omega0 )->value(),
		0.0, 0.0, 0.0, 0.0, 1.0, false,
		parameters_->bundle_calculator()->realvector_parameter_cop( protocols::helical_bundle::BPC_r1_peratom )->value(),
		parameters_->bundle_calculator()->real_parameter_cop( protocols::helical_bundle::BPC_omega1 )->value(),
		parameters_->bundle_calculator()->real_parameter_cop( protocols::helical_bundle::BPC_z1 )->value(),
		parameters_->bundle_calculator()->realvector_parameter_cop( protocols::helical_bundle::BPC_delta_omega1_peratom )->value(),
		parameters_->bundle_calculator()->real_parameter_cop( protocols::helical_bundle::BPC_delta_omega1 )->value(),
		parameters_->bundle_calculator()->realvector_parameter_cop( protocols::helical_bundle::BPC_delta_z1_peratom )->value(),
		parameters_->bundle_calculator()->residues_per_repeat(),
		parameters_->bundle_calculator()->sizevector_parameter_cop( protocols::helical_bundle::BPC_atoms_per_residue )->value(),
		determine_repeating_unit_offset( pose_for_reference, start_position_, end_position_ ),
		failed
	);

	if ( failed ) {
		if ( TR.visible() ) {
			TR << "The current parameter values failed to generate a sensible helix for span " << start_position_ << "-" << end_position_ << "." << std::endl;
		}
		return true; //Failed to generate a sensible helix
	}

	for ( core::Size ir(start_position_); ir<=end_position_; ++ir ) { //Loop through residues in helix.
		for ( core::Size itors(1), itorsmax( pose_for_reference.residue(ir).mainchain_torsions().size() ); itors<=itorsmax; ++itors ) { //Loop through mainchain torsions in residue
			core::id::TorsionID const torsid( ir, core::id::BB, itors );
			core::id::AtomID at1, at2, at3, at4;
			pose_for_reference.conformation().backbone_torsion_angle_atoms( torsid, at1, at2, at3, at4 );

			// Correct the residue indices to match the vector:
			core::Size const offsetval( start_position_ - 2 );
			at1.rsd() -= offsetval;
			at2.rsd() -= offsetval;
			at3.rsd() -= offsetval;
			at4.rsd() -= offsetval;

			debug_assert( torsion_map_out.count(torsid) == 0 );
			core::Real angle;
			numeric::dihedral_degrees( atom_positions[ at1.rsd() ][ at1.atomno()], atom_positions[ at2.rsd() ][at2.atomno()], atom_positions[ at3.rsd() ][at3.atomno()], atom_positions[ at4.rsd() ][at4.atomno()], angle );
			torsion_map_out[torsid] = angle;
		}
	}

	return false;
}

/// @brief Set the start position.
void
HBPHelix::set_start_position(
	core::Size const position_in
) {
	runtime_assert_string_msg( position_in > 0, "Error in HBPHelix::set_start_position(): The position must be greater than zero." );
	start_position_ = position_in;
	if ( end_position_ != 0 ) {
		runtime_assert_string_msg( end_position_ > start_position_, "Error in HBPHelix::set_start_position(): The end position must be greater than the start position." );
	}
}

/// @brief Set the end position.
void
HBPHelix::set_end_position(
	core::Size const position_in
) {
	runtime_assert_string_msg( position_in > 0, "Error in HBPHelix::set_end_position(): The position must be greater than zero." );
	end_position_ = position_in;
	if ( start_position_ != 0 ) {
		runtime_assert_string_msg( end_position_ > start_position_, "Error in HBPHelix::set_end_position(): The end position must be greater than the start position." );
	}
}

/// @brief Set the probability of nucleating a helix.
void
HBPHelix::set_nucleation_prob(
	core::Real const & setting
) {
	runtime_assert_string_msg( setting > 0, "Error in HBPHelix::set_nucleation_prob(): The probability of nucleating a helix must be greater than zero." );
	nucleation_prob_ = setting;
}

/// @brief Set the probability of extending an existing helix.
void
HBPHelix::set_extension_prob(
	core::Real const & setting
) {
	runtime_assert_string_msg( setting > 0, "Error in HBPHelix::set_extension_prob(): The probability of extending a helix must be greater than zero." );
	extension_prob_ = setting;
}

/// @brief Set the probability of shrinking an existing helix.
void
HBPHelix::set_retraction_prob(
	core::Real const & setting
) {
	runtime_assert_string_msg( setting > 0, "Error in HBPHelix::set_retraction_prob(): The probability of shrinking a helix must be greater than zero." );
	retraction_prob_ = setting;
}

/// @brief Set the filename for the Crick params for this helix.
/// @details Used to determine whether two helices are of a matching type.
void
HBPHelix::set_crick_params_filename(
	std::string const & name_in
) {
	crick_params_filename_ = name_in;
}

////////////////////////////// HBPHelixAssignments class //////////////////////////////

/// @brief Default constructor.
HBPHelixAssignments::HBPHelixAssignments() :
	utility::VirtualBase()
{
	clear(); //Initializes everything.
}

/// @brief Copy constructor.
/// @details Note explicit deep-copy.
HBPHelixAssignments::HBPHelixAssignments( HBPHelixAssignments const & src ) :
	VirtualBase( src ),
	helices_() /*Deep-copied below.*/ ,
	global_bundle_calculator_( utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::BundleParametrizationCalculator>( src.global_bundle_calculator_->clone() ) ),
	common_r0_(src.common_r0_),
	common_omega0_(src.common_omega0_),
	common_delta_omega1_(src.common_delta_omega1_),
	confine_to_user_defined_helices_(src.confine_to_user_defined_helices_),
	global_r0_min_(src.global_r0_min_),
	global_r0_max_(src.global_r0_max_),
	global_omega0_min_(src.global_omega0_min_),
	global_omega0_max_(src.global_omega0_max_),
	global_delta_omega1_min_(src.global_delta_omega1_min_),
	global_delta_omega1_max_(src.global_delta_omega1_max_),
	global_crick_params_file_(src.global_crick_params_file_),
	nucleation_prob_(src.nucleation_prob_),
	extension_prob_(src.extension_prob_),
	retraction_prob_(src.retraction_prob_),
	global_perturbation_prob_(src.global_perturbation_prob_),
	global_fractional_perturbation_magnitude_(src.global_fractional_perturbation_magnitude_)
{
	runtime_assert( global_bundle_calculator_ != nullptr );
	helices_.clear();
	helices_.reserve( src.helices_.size() );
	for ( core::Size i(1), imax(src.helices_.size()); i<=imax; ++i ) {
		helices_.push_back( src.helices_[i]->clone() );
	}
}

/// @brief Destructor.
HBPHelixAssignments::~HBPHelixAssignments(){}

/// @brief Create a copy of this object and return an owning pointer to the copy.
HBPHelixAssignmentsOP
HBPHelixAssignments::clone() const {
	return utility::pointer::make_shared< HBPHelixAssignments >( *this );
}

/// @brief Assignment operator.
/// @details Deep-clones helices.
HBPHelixAssignments &
HBPHelixAssignments::operator= (
	HBPHelixAssignments const &src
) {
	if ( &src == this ) return *this;

	helices_.clear();
	helices_.reserve( src.helices_.size() );
	for ( core::Size i(1), imax(src.helices_.size()); i<=imax; ++i ) {
		helices_.push_back( src.helices_[i]->clone() );
	}

	global_bundle_calculator_ = utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::BundleParametrizationCalculator>( src.global_bundle_calculator_->clone() );
	runtime_assert( global_bundle_calculator_ != nullptr );

	common_r0_ = src.common_r0_;
	common_omega0_ = src.common_omega0_;
	common_delta_omega1_ = src.common_delta_omega1_;
	confine_to_user_defined_helices_ = src.confine_to_user_defined_helices_;
	global_r0_min_ = src.global_r0_min_;
	global_r0_max_ = src.global_r0_max_;
	global_omega0_min_ = src.global_omega0_min_;
	global_omega0_max_ = src.global_omega0_max_;
	global_delta_omega1_min_ = src.global_delta_omega1_min_;
	global_delta_omega1_max_ = src.global_delta_omega1_max_;
	global_crick_params_file_ = src.global_crick_params_file_;
	nucleation_prob_ = src.nucleation_prob_;
	extension_prob_ = src.extension_prob_;
	retraction_prob_ = src.retraction_prob_;
	global_perturbation_prob_ = src.global_perturbation_prob_;
	global_fractional_perturbation_magnitude_ = src.global_fractional_perturbation_magnitude_;

	return *this;
}

////// PUBLIC FUNCTIONS //////

/// @brief Reset this object.
void
HBPHelixAssignments::clear() {
	helices_.clear();
	global_bundle_calculator_ = utility::pointer::make_shared< protocols::helical_bundle::BundleParametrizationCalculator >(false);
	common_r0_ = false;
	common_omega0_ = false;
	common_delta_omega1_ = false;
	confine_to_user_defined_helices_ = true;
	global_r0_min_ = 3.0;
	global_r0_max_ = 8.0;
	global_omega0_min_ = numeric::conversions::radians( -3.0 );
	global_omega0_max_ = numeric::conversions::radians( 3.0 );
	global_delta_omega1_min_ = numeric::conversions::radians( -180.0 );
	global_delta_omega1_max_ = numeric::conversions::radians( 180.0 );
	global_crick_params_file_ = "alpha_helix.crick_params";
	nucleation_prob_ = 0.01;
	extension_prob_ = 0.05;
	retraction_prob_ = 0.03;
	global_perturbation_prob_ = 0.1;
	global_fractional_perturbation_magnitude_ = 0.1;
}

/// @brief Set this object up from the contents of a file.
/// @details Warning!  Reads Crick params files from disk!
void
HBPHelixAssignments::initialize_from_file_contents(
	std::string const & file_contents
) {
	clear();
	initialize_global_options_from_file_contents(file_contents);
	initialize_helices_from_file_contents(file_contents);
}

/// @brief Given a sequence index, determine whether it is in a helix.
bool
HBPHelixAssignments::is_in_helix(
	core::Size const sequence_index
) const {
	for ( core::Size i(1), imax(helices_.size()); i<=imax; ++i ) {
		if ( helices_[i]->is_in_helix(sequence_index) ) return true;
	}
	return false;
}

/// @brief Roll a die and decide whether to nucleate a helix at the present position.
/// @details If we do it, then we add one turn of helix as a new helix (updating the helices_ list),
/// and we copy helix paramters from reference_helix_assignments.  The probabilities for nucleation are taken
/// from reference_helix_assignments.
/// @returns True if we nucleated here, false otherwise.
bool
HBPHelixAssignments::try_nucleating_helix(
	core::Size const seqpos,
	core::pose::Pose const &pose,
	HBPHelixAssignments const & reference_helix_assignments
) {
	debug_assert( seqpos > 0 );

	//Do not nucleate if this isn't an allowed region for a helix in the user assignments:
	if ( !(reference_helix_assignments.is_in_helix(seqpos)) ) return false; //We did not nucleate.

	//Do not nucleate if we're already in a helix:
	if ( is_in_helix(seqpos) ) return false; //We did not nucleate.

	//Get the user-defined helix:
	HBPHelixCOP reference_helix( reference_helix_assignments.helix_from_seqpos( seqpos ) );
	debug_assert( reference_helix != nullptr );

	//Get the residues per turn:
	core::Size const residues_per_turn( reference_helix->residues_per_turn() );

	//Do not nucleate if we're at the end of the pose, such that the nucleated helix would extend beyond
	//the pose:
	if ( seqpos > (pose.total_residue() - residues_per_turn + 1) ) return false; //We did not nucleate.

	//Do not nucleate if we would extend beyond the current reference helix and we're not supposed to:
	if ( reference_helix_assignments.confine_to_user_defined_helices() ) {
		if ( seqpos + residues_per_turn - 1 > reference_helix->end_position() ) return false; //We did not nucleate.
	}

	//Do not nucleate unless we roll a suitable random value:
	if ( numeric::random::uniform() > reference_helix->nucleation_prob() ) return false; //We did not nucleate.

	//If we reach here, we're committed to nucleating:
	HBPHelixOP new_helix( add_helix( seqpos, seqpos + residues_per_turn - 1, *reference_helix ) );
	HBPHelixParametersOP params( new_helix->parameters() );

	//Randomize parameters if we're not copying parameters from globals:
	if ( !common_r0_ ) {
		core::conformation::parametric::RealValuedParameterOP r0_param( params->bundle_calculator()->real_parameter( protocols::helical_bundle::BPC_r0 ) );
		r0_param->set_value( numeric::random::uniform() * ( params->r0_range().second - params->r0_range().first ) + params->r0_range().first );
	}
	if ( !common_omega0_ ) {
		core::conformation::parametric::RealValuedParameterOP omega0_param( params->bundle_calculator()->real_parameter( protocols::helical_bundle::BPC_omega0 ) );
		omega0_param->set_value( numeric::random::uniform() * ( params->omega0_range().second - params->omega0_range().first ) + params->omega0_range().first );
	}
	if ( !common_delta_omega1_ ) {
		core::conformation::parametric::RealValuedParameterOP delta_omega1_param( params->bundle_calculator()->real_parameter( protocols::helical_bundle::BPC_delta_omega1 ) );
		delta_omega1_param->set_value( numeric::random::uniform() * ( params->delta_omega1_range().second - params->delta_omega1_range().first ) + params->delta_omega1_range().first );
	}

	check_for_helix_merges( new_helix ); //Merge helices that overlap.

	return true; //We nucleated.
}

/// @brief Roll a die and decide whether to elongate a helix at the present position.
/// @details If we do it, then we add one residue to the helix (updating the helices_ list).
/// The probabilities for elongation are taken from the adjacent helix.
/// @note Does nothing if the present position is not adjacent to a helix.  Automatically merges
/// helices if the present position bridges two helices, and we decide to elongate.  In this case,
/// the parameters of one helix are randomly chosen as the parameters to keep.
/// @returns True if we elongated here, false otherwise.
bool
HBPHelixAssignments::try_elongating_helix(
	core::Size const seqpos,
	core::pose::Pose const &/*pose*/,
	HBPHelixAssignments const & reference_helix_assignments
) {
	if ( is_in_helix(seqpos) ) return false; //We can't elongate if we're already in a helix.
	if ( reference_helix_assignments.confine_to_user_defined_helices() && !reference_helix_assignments.is_in_helix(seqpos) ) return false; //We can't elongate beyond user-defined helix boundaries if this option is set.

	HBPHelixOP adjacent_helix( find_adjacent_helix( seqpos ) ); //Find the helix that we're next to.
	if ( adjacent_helix == nullptr ) return false; //This position isn't adjacent to a helix, so we don't elongate.

	//Roll the dice and decide whether to :
	if ( numeric::random::uniform() > adjacent_helix->extension_prob() ) return false; //We did not elongate.

	//If we reach here, we are committed to elongating.
	if ( adjacent_helix->start_position() == seqpos + 1 ) {
		adjacent_helix->set_start_position(seqpos);
	} else if ( adjacent_helix->end_position() == seqpos - 1 ) {
		adjacent_helix->set_end_position( seqpos );
	} else {
		utility_exit_with_message( "Program error in HBPHelixAssignments::try_elongating_helix().  Please contact Vikram K. Mulligan (vmulligan@flatironinstitute.org)." );
	}

	check_for_helix_merges(); //Merge helices that overlap.
	return true; //We elongated.
}

/// @brief Roll a die and decide whether to contract a helix at the present position.
/// @details If we do it, then we subtract one residue from the helix (updating the helices_ list).
/// The probabilities for retraction are taken from the adjacent helix.  If the helix shrinks to zero size,
/// then we remove it from the helices_ list.
/// @note Does nothing if the present position is not at one of the ends of a helix.
/// @returns True if we retracted here, false otherwise.
bool
HBPHelixAssignments::try_retracting_helix(
	core::Size const seqpos,
	core::pose::Pose const &/*pose*/,
	HBPHelixAssignments const &/*reference_helix_assignments*/
) {
	if ( !is_in_helix(seqpos) ) return false; //If we're not in a helix, do nothing.
	HBPHelixOP containing_helix( helix_from_seqpos( seqpos ) );
	debug_assert( containing_helix != nullptr );
	if ( seqpos != containing_helix->start_position() && seqpos != containing_helix->end_position() ) return false; //Can't contract if we're not at the end of a helix.

	//Roll a die and decide whether to retract:
	if ( numeric::random::uniform() > containing_helix->retraction_prob() ) return false; //We didn't contract.

	//If we reach here, we are committed to retracting.
	if ( containing_helix->start_position() >= containing_helix->end_position() - containing_helix->residues_per_turn() + 1 ) {
		//The helix is one turn, and should be removed.
		remove_helix( containing_helix );
		return true; //We retracted.
	}

	if ( containing_helix->start_position() == seqpos ) {
		containing_helix->set_start_position( seqpos + 1 );
	} else if ( containing_helix->end_position() == seqpos ) {
		containing_helix->set_end_position( seqpos - 1 );
	} else {
		utility_exit_with_message( "Program error in HBPHelixAssignments::try_retracting_helix().  Please contact Vikram K. Mulligan (vmulligan@flatironinstitute.org)." );
	}

	return true; //We retracted.
}

/// @brief Try perturbing the parameters for all helices.
/// @details Updates individual helices if parameters are perturbed.
bool
HBPHelixAssignments::try_perturbing_global_helix_parameters() {
	//Roll a die and decide whether to include a global perturbation move:
	if ( numeric::random::uniform() > global_perturbation_prob_ ) return false; //Do nothing unless random val is less than perturbation probability.

	core::conformation::parametric::RealValuedParameterOP r0_parameter( global_bundle_calculator_->real_parameter( protocols::helical_bundle::BPC_r0 ) );
	r0_parameter->set_value( numeric::clamp( r0_parameter->value() + numeric::random::gaussian() * global_fractional_perturbation_magnitude_ * ( global_r0_max_ - global_r0_min_ ), global_r0_min_, global_r0_max_ ) );
	core::conformation::parametric::RealValuedParameterOP omega0_parameter( global_bundle_calculator_->real_parameter( protocols::helical_bundle::BPC_omega0 ) );
	omega0_parameter->set_value( numeric::clamp( omega0_parameter->value() + numeric::random::gaussian() * global_fractional_perturbation_magnitude_ * ( global_omega0_max_ - global_omega0_min_ ), global_omega0_min_, global_omega0_max_ ) );
	core::conformation::parametric::RealValuedParameterOP delta_omega1_parameter( global_bundle_calculator_->real_parameter( protocols::helical_bundle::BPC_delta_omega1 ) );
	delta_omega1_parameter->set_value( numeric::clamp( delta_omega1_parameter->value() + numeric::random::gaussian() * global_fractional_perturbation_magnitude_ * ( global_delta_omega1_max_ - global_delta_omega1_min_ ), global_delta_omega1_min_, global_delta_omega1_max_ ) );

	if ( !(common_r0_ || common_omega0_ || common_delta_omega1_ ) ) return true; //We're done if we don't have to update individual helices.

	for ( core::Size i(1), imax(helices_.size()); i<=imax; ++i ) {
		protocols::helical_bundle::BundleParametrizationCalculatorOP curcalc( helices_[i]->parameters()->bundle_calculator() );

		if ( common_r0_ ) {
			curcalc->real_parameter( protocols::helical_bundle::BPC_r0 )->set_value( r0_parameter->value() );
		}

		if ( common_omega0_ ) {
			curcalc->real_parameter( protocols::helical_bundle::BPC_omega0 )->set_value( omega0_parameter->value() );
		}

		if ( common_delta_omega1_ ) {
			curcalc->real_parameter( protocols::helical_bundle::BPC_delta_omega1 )->set_value( delta_omega1_parameter->value() );
		}
	}
	return true; //We perturbed.
}

/// @brief Try perturbing the parameters of individual helices, that are not set to copy globals.
bool
HBPHelixAssignments::try_perturbing_local_helix_parameters() {
	if ( common_r0_ && common_omega0_ && common_delta_omega1_ ) return false; //We're done if we don't have to update these.
	for ( core::Size i(1), imax(helices_.size()); i<=imax; ++i ) {
		if ( numeric::random::uniform() > global_perturbation_prob_ ) continue; //Do nothing unless random val is less than perturbation probability.
		HBPHelixParametersOP params( helices_[i]->parameters() );
		protocols::helical_bundle::BundleParametrizationCalculatorOP curcalc( params->bundle_calculator() );

		if ( !common_r0_ ) {
			core::conformation::parametric::RealValuedParameterOP r0_parameter( curcalc->real_parameter( protocols::helical_bundle::BPC_r0 ) );
			r0_parameter->set_value( numeric::clamp( r0_parameter->value() + numeric::random::gaussian() * global_fractional_perturbation_magnitude_ * ( params->r0_range().second - params->r0_range().first ), params->r0_range().first, params->r0_range().second ) );
		}

		if ( !common_omega0_ ) {
			core::conformation::parametric::RealValuedParameterOP omega0_parameter( curcalc->real_parameter( protocols::helical_bundle::BPC_omega0 ) );
			omega0_parameter->set_value( numeric::clamp( omega0_parameter->value() + numeric::random::gaussian() * global_fractional_perturbation_magnitude_ * ( params->omega0_range().second - params->omega0_range().first ), params->omega0_range().first, params->omega0_range().second ) );
		}

		if ( !common_delta_omega1_ ) {
			core::conformation::parametric::RealValuedParameterOP delta_omega1_parameter( curcalc->real_parameter( protocols::helical_bundle::BPC_delta_omega1 ) );
			delta_omega1_parameter->set_value( numeric::clamp( delta_omega1_parameter->value() + numeric::random::gaussian() * global_fractional_perturbation_magnitude_ * ( params->delta_omega1_range().second - params->delta_omega1_range().first ), params->delta_omega1_range().first, params->delta_omega1_range().second ) );
		}
	}
	return true; //We perturbed.
}

/// @brief Generate a map of torsion ID -> torsion value for each helix, and return a vector of these maps.
/// @returns True for FAILURE, false for success.
bool
HBPHelixAssignments::generate_torsion_values_for_helices(
	utility::vector1< std::map < core::id::TorsionID, core::Real > > & torsion_values_out,
	core::pose::Pose const & pose_for_reference
) const {
	if ( helices_.size() == 0 ) {
		torsion_values_out.clear();
		return false;
	}
	torsion_values_out.resize( helices_.size() );
	for ( core::Size i(1), imax(helices_.size()); i<=imax; ++i ) {
		if ( helices_[i]->get_torsions_for_helix( torsion_values_out[i], pose_for_reference ) ) return true;
	}

	return false;
}

/// @brief Given a sequence position, get the helix index that contains it.
/// @details Returns zero if no helix contains this sequence position.
core::Size
HBPHelixAssignments::get_containing_helix_index(
	core::Size const seqpos
) const {
	if ( helices_.size() == 0 ) return 0;
	for ( core::Size i(1), imax(helices_.size()); i<=imax; ++i ) {
		if ( helices_[i]->start_position() <= seqpos && seqpos <= helices_[i]->end_position() ) {
			return i;
		}
	}
	return 0;
}

////// PRIVATE FUNCTIONS //////

/// @brief Get a nonconst owning pointer to the helix containing a given sequence position.
/// @details Returns nullptr if no helix contains that sequence position.
HBPHelixOP
HBPHelixAssignments::helix_from_seqpos(
	core::Size const position
) {
	for ( core::Size i(1), imax(helices_.size()); i<=imax; ++i ) {
		HBPHelixOP curhelix( helices_[i] );
		debug_assert( curhelix != nullptr );
		if ( curhelix->start_position() <= position && curhelix->end_position() >= position ) {
			return curhelix;
		}
	}
	return nullptr;
}

/// @brief Get a const owning pointer to the helix containing a given sequence position.
/// @details Returns nullptr if no helix contains that sequence position.
HBPHelixCOP
HBPHelixAssignments::helix_from_seqpos(
	core::Size const position
) const {
	for ( core::Size i(1), imax(helices_.size()); i<=imax; ++i ) {
		HBPHelixOP curhelix( helices_[i] );
		debug_assert( curhelix != nullptr );
		if ( curhelix->start_position() <= position && curhelix->end_position() >= position ) {
			return curhelix;
		}
	}
	return nullptr;
}

/// @brief Loop through all helices and find a helix that this sequence position is adjacent to.
/// @details.  If the position is adjacent to no helix, returns nullptr.  If the position is adjacent
/// to two helices (i.e. it bridges two helices separated by one residue), we pick one with 50/50 odds.
HBPHelixOP
HBPHelixAssignments::find_adjacent_helix(
	core::Size seqpos
) {
	if ( is_in_helix(seqpos) ) return nullptr;

	HBPHelixOP upper_helix(nullptr), lower_helix(nullptr);
	for ( core::Size i(1), imax(helices_.size()); i<=imax; ++i ) {
		if ( helices_[i]->start_position() == seqpos + 1 ) upper_helix = helices_[i];
		else if ( helices_[i]->end_position() == seqpos - 1 ) lower_helix = helices_[i];

		if ( upper_helix != nullptr && lower_helix != nullptr ) break;
	}

	if ( upper_helix != nullptr ) {
		if ( lower_helix != nullptr ) { // neither is null -- randomly select one to return.
			if ( numeric::random::uniform() > 0.5 ) {
				return upper_helix;
			} else {
				return lower_helix;
			}
		} else { // lower is null, upper is not
			return upper_helix;
		}
	} else { // upper is null
		if ( lower_helix != nullptr ) { // upper is null, lower is not
			return lower_helix;
		} // Both null handled below
	}

	return nullptr;
}

/// @brief Are two helices of the same type?
bool
HBPHelixAssignments::helix_types_match(
	HBPHelix const &first_helix,
	HBPHelix const &second_helix
) const {
	return first_helix.crick_params_filename() == second_helix.crick_params_filename();
}


/// @brief Check whether the addition or elongation of current_helix causes a new contiguous stretch of helix, and
/// if so, replace helix records in the helices_ list with a single contiguous helix.
/// @details When helices are merged, the helix contributing parameters is chosen randomly. If current_helix is not nullptr,
/// then we do not allow that helix to contribute the parameters -- i.e. they are chosen randomly from one of the other helices,
/// or we just use the other helix if there's only one.
void
HBPHelixAssignments::check_for_helix_merges(
	HBPHelixCOP const current_helix /*= nullptr*/
) {
	if ( helices_.size() < 2 ) return;

	while ( true ) { //Keep looping until we find no more merges.
		bool found_overlap(false);
		for ( core::Size i(2), imax(helices_.size()); i<=imax; ++i ) {
			for ( core::Size j(1); j<i; ++j ) {
				HBPHelixCOP first_helix( helices_[i]->start_position() < helices_[j]->start_position() ? helices_[i] : helices_[j] );
				HBPHelixCOP second_helix( helices_[i]->start_position() < helices_[j]->start_position() ? helices_[j] : helices_[i] );
				if ( first_helix->end_position() >= second_helix->start_position() - 1 ) {
					if ( helix_types_match( *first_helix, *second_helix ) ) {
						found_overlap = true;
						merge_helices(i, j, current_helix);
						break;
					} else { // Helix types DON'T match.
						found_overlap = true;
						resolve_mismatching_helix_overlap( i, j );
						break;
					}
				}
			}
			if ( found_overlap ) break;
		}
		if ( !found_overlap ) {
			break;
		}
	}
}

/// @brief Merge two helices together.
/// @details If current_helix is not nullptr, then we randomly pick one of the two helices to provide
/// the parameters.  If current_helix is nullptr and it matches one of the two helices, then the other
/// provides the parameters.
void
HBPHelixAssignments::merge_helices(
	core::Size const first_helix_index,
	core::Size const second_helix_index,
	HBPHelixCOP const current_helix /*= nullptr*/
) {
	debug_assert( first_helix_index > 0 && first_helix_index <= helices_.size() );
	debug_assert( second_helix_index > 0 && second_helix_index <= helices_.size() );
	debug_assert( first_helix_index != second_helix_index );

	HBPHelixOP helix1( helices_[first_helix_index] );
	HBPHelixOP helix2( helices_[second_helix_index] );

	debug_assert( helix1 != helix2 );
	debug_assert( helix1 != nullptr );
	debug_assert( helix2 != nullptr );

	HBPHelixOP helix_to_copy;
	if ( helix1 == current_helix ) helix_to_copy = helix2;
	else if ( helix2 == current_helix ) helix_to_copy = helix1;
	else { //Should cover the case of curent_helix == nullptr.
		if ( numeric::random::uniform() > 0.5 ) {
			helix_to_copy = helix1;
		} else {
			helix_to_copy = helix2;
		}
	}

	core::Size const start_res( std::min( helix1->start_position(), helix2->start_position() ) );
	core::Size const end_res( std::max( helix1->end_position(), helix2->end_position() ) );
	add_helix( start_res, end_res, *helix_to_copy );
	remove_helix(helix1);
	remove_helix(helix2);
}

/// @brief Given two overlapping helices of different types, prevent them from overlapping.
/// @details The cutoff point is randomly chosen within the overlap region.
void
HBPHelixAssignments::resolve_mismatching_helix_overlap(
	core::Size const helix1_index,
	core::Size const helix2_index
) {
	if ( helix1_index == helix2_index ) return; //Do nothing if the same helix.
	HBPHelixOP helix1( helices_[helix1_index] );
	HBPHelixOP helix2( helices_[helix2_index] );

	bool helix1_is_first( true );
	if ( helix1->start_position() >= helix2->start_position() ) {
		helix1_is_first = false;
	}
	HBPHelixOP earlier_helix( helix1_is_first ? helix1 : helix2 );
	HBPHelixOP later_helix( helix1_is_first ? helix2 : helix1 );

	if ( earlier_helix->end_position() < later_helix->start_position() ) return; //Do nothing if the helices don't overlap.
	core::Size const overlap_start( later_helix->start_position() ); //The first residue of the overlap region.
	core::Size const overlap_end( later_helix->end_position() < earlier_helix->end_position() ? later_helix->end_position() + 1 : earlier_helix->end_position() + 1 ); //The first residue PAST the overlap region.
	core::Size cutoff_residue( numeric::random::random_range( overlap_start, overlap_end ) );

	if ( cutoff_residue < earlier_helix->start_position() + earlier_helix->residues_per_turn() ) { //There's less than a turn of the earlier helix left:
		remove_helix( earlier_helix );
	} else {
		earlier_helix->set_end_position( cutoff_residue - 1 );
	}

	if ( cutoff_residue > later_helix->end_position() - later_helix->residues_per_turn() + 1 ) { //There's less than a turn of the later helix left:
		remove_helix(later_helix);
	} else {
		later_helix->set_start_position( cutoff_residue );
	}
}

/// @brief Given file contents, pull out the first BEGIN_GLOBALS...END_GLOBALS block, and throw an
/// error message if any other occurrences of this exist in the file.
void
HBPHelixAssignments::pull_out_global_block(
	std::string const &file_contents,
	std::string & global_block
) const {
	static std::string const errmsg( "Error when parsing helix definitions file in HBPHelixAssignments::pull_out_global_block(): " );
	std::stringstream file_ss( file_contents );
	std::string line, chunk;
	bool found( false ), in_block( false );
	std::stringstream outstr;

	while ( !file_ss.eof() ) {
		getline( file_ss, line );
		utility::strip( line, " \t\n" );
		if ( line.empty() ) continue;
		if ( line[0] == '#' ) continue;

		std::stringstream line_ss( line );
		utility::parse_out( line_ss, chunk );
		if ( !found && !in_block ) {
			if ( chunk == "BEGIN_GLOBALS" ) {
				found = true;
				in_block = true;
			} else if ( chunk == "END_GLOBALS" ) {
				utility_exit_with_message( errmsg + "\"END_GLOBALS\" found without \"BEGIN_GLOBALS\"." );
			}
		} else if ( found && in_block ) {
			if ( chunk == "END_GLOBALS" ) {
				in_block = false;
			} else if ( chunk == "BEGIN_GLOBALS" ) {
				utility_exit_with_message( errmsg + "Second \"BEGIN_GLOBALS\" found within a \"BEGIN_GLOBALS\"...\"END_GLOBALS\" block." );
			} else {
				outstr << line << "\n";
			}
		} else if ( found && !in_block ) {
			runtime_assert_string_msg( chunk != "BEGIN_GLOBALS", errmsg + "Second \"BEGIN_GLOBALS\" line found after first \"BEGIN_GLOBALS\"...\"END_GLOBALS\" block." );
			runtime_assert_string_msg( chunk != "END_GLOBALS", errmsg + "Stray \"END_GLOBALS\" line found after first \"BEGIN_GLOBALS\"...\"END_GLOBALS\" block." );
		}
	}
	runtime_assert_string_msg( !( found && in_block ), errmsg + "\"BEGIN_GLOBALS\" found without corresponding \"END_GLOBALS\"." );
	global_block = outstr.str();
}

/// @brief Read the BEGIN_GLOBALS...END_GLOBALS block in the file contents and set up global options.
/// @details Warning!  Reads Crick params files from disk!
void
HBPHelixAssignments::initialize_global_options_from_file_contents(
	std::string const &file_contents
) {
	static std::string const errmsg( "Error when parsing helix definitions file in HBPHelixAssignments::initialize_global_options_from_file_contents(): " );

	std::string global_block;
	pull_out_global_block(file_contents, global_block);
	if ( global_block.empty() ) return; //No global block, so we're done.
	std::stringstream global_ss( global_block );
	std::string line;
	while ( !global_ss.eof() ) {
		std::getline( global_ss, line );
		utility::strip( line, " \t\n");
		if ( line[0] == '#' ) continue; //Ignore comments.
		if ( line.empty() ) continue; //Ignore blank lines.
		std::stringstream line_ss( line );
		std::string chunk;
		utility::parse_out( line_ss, chunk, errmsg );
		if ( chunk == "COMMON_R0" ) {
			common_r0_ = utility::parse_boolean( line_ss, errmsg );
		} else if ( chunk == "COMMON_OMEGA0" ) {
			common_omega0_ = utility::parse_boolean( line_ss, errmsg );
		} else if ( chunk == "COMMON_DELTA_OMEGA1" ) {
			common_delta_omega1_ = utility::parse_boolean( line_ss, errmsg );
		} else if ( chunk == "CONFINE_TO_USER_DEFINED_HELICES" ) {
			confine_to_user_defined_helices_ = utility::parse_boolean( line_ss, errmsg );
		} else if ( chunk == "R0_MIN" ) {
			utility::parse_out( line_ss, global_r0_min_, errmsg );
		} else if ( chunk == "R0_MAX" ) {
			utility::parse_out( line_ss, global_r0_max_, errmsg );
		} else if ( chunk == "OMEGA0_MIN" ) {
			utility::parse_out( line_ss, global_omega0_min_, errmsg );
			global_omega0_min_ = numeric::conversions::radians( global_omega0_min_ );
		} else if ( chunk == "OMEGA0_MAX" ) {
			utility::parse_out( line_ss, global_omega0_max_, errmsg );
			global_omega0_max_ = numeric::conversions::radians( global_omega0_max_ );
		} else if ( chunk == "DELTA_OMEGA1_MIN" ) {
			utility::parse_out( line_ss, global_delta_omega1_min_, errmsg );
			global_delta_omega1_min_ = numeric::conversions::radians( global_delta_omega1_min_ );
		} else if ( chunk == "DELTA_OMEGA1_MAX" ) {
			utility::parse_out( line_ss, global_delta_omega1_max_, errmsg );
			global_delta_omega1_max_ = numeric::conversions::radians( global_delta_omega1_max_ );
		} else if ( chunk == "NUCLEATION_PROB" ) {
			utility::parse_out( line_ss, nucleation_prob_, errmsg );
		} else if ( chunk == "EXTENSION_PROB" ) {
			utility::parse_out( line_ss, extension_prob_, errmsg );
		} else if ( chunk == "RETRACTION_PROB" ) {
			utility::parse_out( line_ss, retraction_prob_, errmsg );
		} else if ( chunk == "GLOBAL_PERTURBATION_PROB" ) {
			utility::parse_out( line_ss, global_perturbation_prob_, errmsg );
		} else if ( chunk == "FRACTIONAL_PERTURBATION_MAGNITUDE" ) {
			utility::parse_out( line_ss, global_fractional_perturbation_magnitude_, errmsg );
		} else if ( chunk == "CRICK_PARAMS_FILE" ) {
			utility::parse_out( line_ss, global_crick_params_file_, errmsg );
		} else {
			utility_exit_with_message( errmsg + " Could not parse \"" + line_ss.str() + "\"." );
		}
	}
	global_bundle_calculator_->init_from_file( global_crick_params_file_ ); //Reads from file!

	// Set values to random initial values within appropriate parameter ranges:
	global_bundle_calculator_->real_parameter( protocols::helical_bundle::BPC_r0 )->set_value( numeric::random::uniform() * ( global_r0_max_ - global_r0_min_ ) + global_r0_min_ );
	global_bundle_calculator_->real_parameter( protocols::helical_bundle::BPC_omega0 )->set_value( numeric::random::uniform() * ( global_omega0_max_ - global_omega0_min_ ) + global_omega0_min_ );
	global_bundle_calculator_->real_parameter( protocols::helical_bundle::BPC_delta_omega1 )->set_value( numeric::random::uniform() * ( global_delta_omega1_max_ - global_delta_omega1_min_ ) + global_delta_omega1_min_ );
}

/// @brief Given file contents, pull out all of the stuff in between the BEGIN_HELIX and END_HELIX lines.
/// @details Omits the BEGIN_HELIX and END_HELIX lines.
void
HBPHelixAssignments::pull_out_helix_blocks(
	std::string const &file_contents,
	utility::vector1<std::string> &helix_blocks
) const {
	static std::string const errmsg( "Error when parsing helix definitions file in HBPHelixAssignments::pull_out_helix_blocks(): " );
	helix_blocks.clear();

	std::stringstream file_ss( file_contents );
	std::string line, chunk;

	bool in_block( false );
	std::string current_block("");

	while ( !file_ss.eof() ) {
		getline( file_ss, line );
		utility::strip( line, " \t\n");
		if ( line.empty() || line[0] == '#' ) continue;

		std::stringstream line_ss( line );
		utility::parse_out(line_ss, chunk, errmsg);
		if ( !in_block ) {
			runtime_assert_string_msg( chunk != "END_HELIX", errmsg + "Stray \"END_HELIX\" found without a preceding \"BEGIN_HELIX\"." );

			if ( chunk == "BEGIN_HELIX" ) {
				in_block = true;
			}
		} else { // We ARE in a block
			runtime_assert_string_msg( chunk != "BEGIN_GLOBALS" && chunk != "END_GLOBALS", errmsg + "Stray \"" + chunk + "\" found within a \"BEGIN_HELIX\"...\"END_HELIX\" block." );
			runtime_assert_string_msg( chunk != "BEGIN_HELIX", errmsg + "Stray \"BEGIN_HELIX\" found within a \"BEGIN_HELIX\"...\"END_HELIX\" block." );
			if ( chunk == "END_HELIX" ) {
				helix_blocks.push_back( current_block );
				current_block = "";
				in_block = false;
			} else {
				current_block += line;
				current_block += "\n";
			}
		}
	}
	runtime_assert_string_msg( !in_block, errmsg + "\"BEGIN_HELIX\" found with no corresponding \"END_HELIX\"." );
}

/// @brief Read BEGIN_HELIX...END_HELIX blocks in the file contents and set up local helix options.
/// @details Warning!  Triggers read from disk if Crick params file is set.
void
HBPHelixAssignments::initialize_helices_from_file_contents(
	std::string const &file_contents
) {
	static std::string const errmsg( "Error when parsing helix definitions file in HBPHelixAssignments::initialize_helices_from_file_contents(): " );

	utility::vector1< std::string > helix_blocks;
	pull_out_helix_blocks( file_contents, helix_blocks );

	for ( core::Size i(1), imax(helix_blocks.size()); i<=imax; ++i ) {
		HBPHelixOP curhelix( add_helix() );
		HBPHelixParametersOP params( curhelix->parameters() );

		std::stringstream curblock_ss( helix_blocks[i] );
		std::string line, chunk;

		std::string crick_params_filename;

		while ( !curblock_ss.eof() ) {
			getline( curblock_ss, line );
			utility::strip( line, " \t\n" );
			if ( line.empty() || line[0] == '#' ) continue;
			std::stringstream line_ss( line );
			utility::parse_out( line_ss, chunk, errmsg );
			if ( chunk == "START_RES" ) {
				core::Size start_res;
				utility::parse_out( line_ss, start_res, errmsg );
				curhelix->set_start_position(start_res);
			} else if ( chunk == "END_RES" ) {
				core::Size end_res;
				utility::parse_out( line_ss, end_res, errmsg );
				curhelix->set_end_position(end_res);
			} else if ( chunk == "R0_MIN" ) {
				core::Real val;
				utility::parse_out( line_ss, val, errmsg );
				runtime_assert_string_msg( !common_r0_, errmsg + "An \"R0_MIN\" line was found in a helix definition, but \"COMMON_R0 TRUE\" was set globally." );
				params->set_r0_min( val );
			} else if ( chunk == "R0_MAX" ) {
				core::Real val;
				utility::parse_out( line_ss, val, errmsg );
				runtime_assert_string_msg( !common_r0_, errmsg + "An \"R0_MAX\" line was found in a helix definition, but \"COMMON_R0 TRUE\" was set globally." );
				params->set_r0_max( val );
			} else if ( chunk == "OMEGA0_MIN" ) {
				core::Real val;
				utility::parse_out( line_ss, val, errmsg );
				runtime_assert_string_msg( !common_omega0_, errmsg + "An \"OMEGA0_MIN\" line was found in a helix definition, but \"COMMON_OMEGA0 TRUE\" was set globally." );
				params->set_omega0_min( numeric::conversions::radians( val ) );
			} else if ( chunk == "OMEGA0_MAX" ) {
				core::Real val;
				utility::parse_out( line_ss, val, errmsg );
				runtime_assert_string_msg( !common_omega0_, errmsg + "An \"OMEGA0_MAX\" line was found in a helix definition, but \"COMMON_OMEGA0 TRUE\" was set globally." );
				params->set_omega0_max( numeric::conversions::radians( val ) );
			} else if ( chunk == "DELTA_OMEGA1_MIN" ) {
				core::Real val;
				utility::parse_out( line_ss, val, errmsg );
				runtime_assert_string_msg( !common_delta_omega1_, errmsg + "An \"DELTA_OMEGA1_MIN\" line was found in a helix definition, but \"COMMON_DELTA_OMEGA1 TRUE\" was set globally." );
				params->set_delta_omega1_min( numeric::conversions::radians( val ) );
			} else if ( chunk == "DELTA_OMEGA1_MAX" ) {
				core::Real val;
				utility::parse_out( line_ss, val, errmsg );
				runtime_assert_string_msg( !common_delta_omega1_, errmsg + "An \"DELTA_OMEGA1_MAX\" line was found in a helix definition, but \"COMMON_DELTA_OMEGA1 TRUE\" was set globally." );
				params->set_delta_omega1_max( numeric::conversions::radians( val ) );
			} else if ( chunk == "NUCLEATION_PROB" ) {
				core::Real val;
				utility::parse_out( line_ss, val, errmsg );
				curhelix->set_nucleation_prob(val);
			} else if ( chunk == "EXTENSION_PROB" ) {
				core::Real val;
				utility::parse_out( line_ss, val, errmsg );
				curhelix->set_extension_prob(val);
			} else if ( chunk == "RETRACTION_PROB" ) {
				core::Real val;
				utility::parse_out( line_ss, val, errmsg );
				curhelix->set_retraction_prob(val);
			} else if ( chunk == "CRICK_PARAMS_FILE" ) {
				utility::parse_out( line_ss, crick_params_filename, errmsg );
			} else {
				utility_exit_with_message( errmsg + "Could not parse line \"" + line + "\"." );
			}
		}
		if ( !crick_params_filename.empty() ) {
			params->create_calculator_from_file( crick_params_filename ); //Reads from disk!
			curhelix->set_crick_params_filename( crick_params_filename );
		} // Else, this remains a default calculator that was cloned from the global calculator when the helix was created, above.
	}
}

/// @brief Add a helix and return an owning pointer to it.
/// @details Parameters are initialized from globals.
HBPHelixOP
HBPHelixAssignments::add_helix() {
	HBPHelixOP new_helix( utility::pointer::make_shared< HBPHelix >() );
	new_helix->set_nucleation_prob( nucleation_prob_ );
	new_helix->set_extension_prob( extension_prob_ );
	new_helix->set_retraction_prob( retraction_prob_ );
	new_helix->set_crick_params_filename( global_crick_params_file_ );

	HBPHelixParametersOP params( new_helix->parameters() );
	params->set_r0_min( global_r0_min_ );
	params->set_r0_max( global_r0_max_ );
	params->set_omega0_min( global_omega0_min_ );
	params->set_omega0_max( global_omega0_max_ );
	params->set_delta_omega1_min( global_delta_omega1_min_ );
	params->set_delta_omega1_max( global_delta_omega1_max_ );
	params->set_calculator( global_bundle_calculator_ ); //Clones this calculator

	helices_.push_back( new_helix );
	return new_helix;
}

/// @brief Add a helix and return an owning pointer to it, initializing its start, end, and parameters.
/// @details Parameters are initialized from reference_helix.
HBPHelixOP
HBPHelixAssignments::add_helix(
	core::Size const startpos,
	core::Size const endpos,
	HBPHelix const &reference_helix
) {
	HBPHelixOP new_helix( utility::pointer::make_shared< HBPHelix >(reference_helix) );
	new_helix->set_start_position(startpos);
	new_helix->set_end_position(endpos);
	helices_.push_back(new_helix);
	return new_helix;
}

/// @brief Remove a helix.
void
HBPHelixAssignments::remove_helix(
	HBPHelixOP helix_to_remove
) {
	for ( utility::vector1< HBPHelixOP >::iterator it( helices_.begin() ); it!=helices_.end(); ++it ) {
		if ( *it == helix_to_remove ) {
			helices_.erase(it);
			break;
		}
	}
}

} //protocols
} //helical_bundle_predict






