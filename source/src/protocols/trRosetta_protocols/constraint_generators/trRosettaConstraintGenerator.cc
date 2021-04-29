// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta_protocols/constraint_generators/trRosettaConstraintGenerator.cc
/// @brief A module that runs a trRosetta neural network on an input multiple sequence alignment and
/// uses the output to apply distance and/or angle constraints to a pose for subsequent structure
/// prediction or refinement.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Unit headers
#include <protocols/trRosetta_protocols/constraint_generators/trRosettaConstraintGenerator.hh>
#include <protocols/trRosetta_protocols/constraint_generators/trRosettaConstraintGeneratorCreator.hh>
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>

// Protocols headers
#include <protocols/trRosetta/trRosettaOutputs_v1.hh>
#include <protocols/trRosetta/trRosettaProtocol_v1.hh>

// Core headers
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/func/SplineFunc.hh>

// Numeric headers
#include <numeric/constants.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/tensorflow_manager/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/trRosetta.OptionKeys.gen.hh>
#include <basic/citation_manager/CitationManager.hh>
#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/pointer/memory.hh>
#include <utility/file/file_sys_util.hh>

static basic::Tracer TR( "protocols.trRosetta_protocols.constraint_generators.trRosettaConstraintGenerator" );

namespace protocols {
namespace trRosetta_protocols {
namespace constraint_generators {

protocols::constraint_generator::ConstraintGeneratorOP
trRosettaConstraintGeneratorCreator::create_constraint_generator() const
{
	return utility::pointer::make_shared< trRosettaConstraintGenerator >();
}

std::string
trRosettaConstraintGeneratorCreator::keyname() const
{
	return trRosettaConstraintGenerator::class_name();
}

/// @brief Default constructor.
trRosettaConstraintGenerator::trRosettaConstraintGenerator() :
	protocols::constraint_generator::ConstraintGenerator( trRosettaConstraintGenerator::class_name() )
{
#ifdef USE_TENSORFLOW
	init_from_commandline( basic::options::option );
#endif //USE_TENSORFLOW
}

/// @brief OptionsCollection constructor.
trRosettaConstraintGenerator::trRosettaConstraintGenerator(
#ifdef USE_TENSORFLOW
	utility::options::OptionCollection const & opts
#else
	utility::options::OptionCollection const &
#endif
) :
	protocols::constraint_generator::ConstraintGenerator( trRosettaConstraintGenerator::class_name() )
{
#ifdef USE_TENSORFLOW
	init_from_commandline(opts);
#endif //USE_TENSORFLOW
}

#ifdef USE_TENSORFLOW
/// @brief Options constructor.
trRosettaConstraintGenerator::trRosettaConstraintGenerator(
	std::string const & msa_file
) :
	protocols::constraint_generator::ConstraintGenerator( trRosettaConstraintGenerator::class_name() )
{
	init_from_commandline( basic::options::option );
	set_msa_file(msa_file);
}
#endif //USE_TENSORFLOW

/// @brief Destructor.
trRosettaConstraintGenerator::~trRosettaConstraintGenerator() = default;

/// @brief Copy this object and return owning pointer to copy.
protocols::constraint_generator::ConstraintGeneratorOP
trRosettaConstraintGenerator::clone() const
{
	return utility::pointer::make_shared< trRosettaConstraintGenerator >( *this );
}

std::string
trRosettaConstraintGenerator::class_name()
{
	return "trRosettaConstraintGenerator";
}

void
trRosettaConstraintGenerator::parse_tag(
#ifdef USE_TENSORFLOW
	utility::tag::TagCOP tag,
#else
	utility::tag::TagCOP,
#endif
	basic::datacache::DataMap &
) {
#ifdef USE_TENSORFLOW
	if ( tag->hasOption("msa_file") ) {
		set_msa_file( tag->getOption<std::string>("msa_file") );
	}
	if ( tag->hasOption("generate_distance_constraints") ) {
		set_generate_dist_constraints( tag->getOption< bool >("generate_distance_constraints") );
	}
	if ( tag->hasOption("generate_theta_constraints") ) {
		set_generate_theta_constraints( tag->getOption< bool >("generate_theta_constraints") );
	}
	if ( tag->hasOption("generate_phi_constraints") ) {
		set_generate_phi_constraints( tag->getOption< bool >("generate_phi_constraints") );
	}
	if ( tag->hasOption("generate_omega_constraints") ) {
		set_generate_omega_constraints( tag->getOption< bool >("generate_omega_constraints") );
	}

	set_prob_cutoffs(
		tag->getOption< core::Real >( "distance_constraint_prob_cutoff", dist_prob_cutoff() ),
		tag->getOption< core::Real >( "omega_constraint_prob_cutoff", omega_prob_cutoff() ),
		tag->getOption< core::Real >( "theta_constraint_prob_cutoff", theta_prob_cutoff() ),
		tag->getOption< core::Real >( "phi_constraint_prob_cutoff", phi_prob_cutoff() )
	);

	set_constraint_weights(
		tag->getOption< core::Real >( "distance_constraint_Weight", distance_constraint_weight() ),
		tag->getOption< core::Real >( "omega_constraint_Weight", omega_constraint_weight() ),
		tag->getOption< core::Real >( "theta_constraint_Weight", theta_constraint_weight() ),
		tag->getOption< core::Real >( "phi_constraint_Weight", phi_constraint_weight() )
	);
#else // !USE_TENSORFLOW
	utility_exit_with_message(
		"Error in trRosettaConstraintGenerator::parse_tag(): The trRosettaConstraintGenerator requires compilation with Tensorflow support.\n\n"
		+ basic::tensorflow_manager::get_tensorflow_compilation_instructions( "trRosettaConstraintGenerator", false )
	);
#endif // USE_TENSORFLOW
}

core::scoring::constraints::ConstraintCOPs
trRosettaConstraintGenerator::apply(
#ifdef USE_TENSORFLOW
	core::pose::Pose const & pose
#else
	core::pose::Pose const &
#endif
) const {
#ifdef USE_TENSORFLOW
	// Check whether we have already run the neural network, and run it
	// if we have not:
	if ( trRosetta_outputs_ == nullptr ) {
		load_msa_and_run_neural_network();
	}

	// Generate the constraints from the output of the network:
	return generate_constraints(pose);
#else // !USE_TENSORFLOW
	utility_exit_with_message(
		"Error in trRosettaConstraintGenerator::apply(): The trRosettaConstraintGenerator requires compilation with Tensorflow support.\n\n"
		+ basic::tensorflow_manager::get_tensorflow_compilation_instructions( "trRosettaConstraintGenerator", false )
	);
	return utility::vector1< core::scoring::constraints::ConstraintCOP >(); //To keep the compiler happy.
#endif // USE_TENSORFLOW
}

void
trRosettaConstraintGeneratorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	trRosettaConstraintGenerator::provide_xml_schema( xsd );
}

void
trRosettaConstraintGenerator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default( "msa_file", xs_string,
		"Filename for a multiple sequence alignment file, in a3m format.  "
		"Dashes indicate gap sequences, and lowercase characters will be removed "
		"(and flanking regions ligated).  If not provided, the commandline option "
		"-trRosetta:msa_file will be used.  One or the other is required.", "" )

		+ XMLSchemaAttribute::attribute_w_default( "generate_distance_constraints", xsct_rosetta_bool,
		"Set whether this will generate distance constraints.  Distance constraints constrain "
		"the distance between pairs of amino acids.  These are symmetric, and are only "
		"generated once per amino acid pair (since dist(a,b) == dist(b,a)).  Defaults to "
		"commandline setting -trRosetta:use_distance_constraints.", "true"
		)

		+ XMLSchemaAttribute::attribute_w_default( "generate_omega_constraints", xsct_rosetta_bool,
		"Set whether this will generate omega dihedral constraints.  Omega constraints "
		"constrain the dihedral between CA1-CB1-CB2-CA2 in pairs of amino acids.  These "
		"are symmetric, and are only generated once per amino acid pair (since omega(a,b) "
		"== omega(b,a)).  Note that this is NOT omega the backbone dihedral torsion!  Defaults to "
		"commandline setting -trRosetta:use_omega_constraints.", "true"
		)

		+ XMLSchemaAttribute::attribute_w_default( "generate_theta_constraints", xsct_rosetta_bool,
		"Set whether this will generate theta dihedral constraints.  Theta constraints "
		"constrain the dihedral between N1-CA1-CB1-CB2 in pairs of amino acids.  These are "
		"asymmetric (i.e. theta(a,b)!=theta(b,a)), so there are two per amino acid pair "
		"(unless a == b, which is skipped).  Defaults to commandline setting "
		"-trRosetta:use_theta_constraints.", "true"
		)

		+ XMLSchemaAttribute::attribute_w_default( "generate_phi_constraints", xsct_rosetta_bool,
		"Set whether this will generate phi angle constraints.  Phi constraints constrain "
		"the angle between CA1-CB1-CB2 in pairs of amino acids.  These are asymmetric "
		"(i.e. phi(a,b)!=phi(b,a)), so there are two per amino acid pair (unless a == b, "
		"which is skipped).  Note that this is NOT phi the backbone dihedral torsion!  Defaults to "
		"commandline setting -trRosetta:use_phi_constraints.", "true"
		)

		+ XMLSchemaAttribute::attribute_w_default( "distance_constraint_prob_cutoff", xsct_real,
		"Set the probability cutoff below which we omit a distance constraint.  Default 0.05, or "
		"whatever is set on the commandline with the -trRosetta:distance_constraint_prob_cutoff "
		"commandline option.",
		"0.05" )

		+ XMLSchemaAttribute::attribute_w_default( "omega_constraint_prob_cutoff", xsct_real,
		"Set the probability cutoff below which we omit a omega dihedral constraint.  Default 0.55, or "
		"whatever is set on the commandline with the -trRosetta:omega_constraint_prob_cutoff "
		"commandline option.",
		"0.55" )

		+ XMLSchemaAttribute::attribute_w_default( "theta_constraint_prob_cutoff", xsct_real,
		"Set the probability cutoff below which we omit a theta dihedral constraint.  Default 0.55, or "
		"whatever is set on the commandline with the -trRosetta:theta_constraint_prob_cutoff "
		"commandline option.",
		"0.55" )

		+ XMLSchemaAttribute::attribute_w_default( "phi_constraint_prob_cutoff", xsct_real,
		"Set the probability cutoff below which we omit a phi angle constraint.  Default 0.65, or "
		"whatever is set on the commandline with the -trRosetta:phi_constraint_prob_cutoff "
		"commandline option.",
		"0.65" )

		+ XMLSchemaAttribute::attribute_w_default( "distance_constraint_weight", xsct_real,
		"Set the weight for trRosetta-generated distance constraints.  Defaults to 1.0, or whatever "
		"was set on the commandline with the -trRosetta:distance_constraint_weight commandline option.",
		"1.0" )

		+ XMLSchemaAttribute::attribute_w_default( "omega_constraint_weight", xsct_real,
		"Set the weight for trRosetta-generated omega dihedral constraints.  Defaults to 1.0, or whatever "
		"was set on the commandline with the -trRosetta:omega_constraint_weight commandline option.",
		"1.0" )

		+ XMLSchemaAttribute::attribute_w_default( "theta_constraint_weight", xsct_real,
		"Set the weight for trRosetta-generated theta dihedral constraints.  Defaults to 1.0, or whatever "
		"was set on the commandline with the -trRosetta:theta_constraint_weight commandline option.",
		"1.0" )

		+ XMLSchemaAttribute::attribute_w_default( "phi_constraint_weight", xsct_real,
		"Set the weight for trRosetta-generated phi angle constraints.  Defaults to 1.0, or whatever "
		"was set on the commandline with the -trRosetta:phi_constraint_weight commandline option.",
		"1.0" );

	protocols::constraint_generator::ConstraintGeneratorFactory::xsd_constraint_generator_type_definition_w_attributes(
		xsd,
		class_name(),
		"The trRosettaConstraintGenerator takes as input a file containing a multiple sequence alignment, feeds this "
		"to the trRosetta neural network, and uses the output to generate distance and angle constraints between "
		"pairs of residues as described in Yang et al. (2020) Improved protein "
		"structure prediction using predicted interresidue orientations. Proc. Natl. Acad. Sci. USA "
		"117(3):1496-503. https://doi.org/10.1073/pnas.1914677117.\n\n"
		+ basic::tensorflow_manager::get_tensorflow_compilation_instructions( "trRosettaConstraintGenerator", true ),
		attlist );
}

/// @brief Provide citations to the passed CitationCollectionList.
/// This allows the constraint generator to provide citations for itself
/// and for any modules that it invokes.
/// @details The default implementation of this does nothing.  Should be overridden
/// by derived classes.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
trRosettaConstraintGenerator::provide_citation_info(
	basic::citation_manager::CitationCollectionList & citations
) const {
	provide_citation_info_static( citations );
}

/// @brief Provide citations for hte trRosettaConstraintGenerator, without needing an instance.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
/*static*/
void
trRosettaConstraintGenerator::provide_citation_info_static(
	basic::citation_manager::CitationCollectionList & citations
) {
	using namespace basic::citation_manager;
	citations.add(
		utility::pointer::make_shared< UnpublishedModuleInfo >(
		class_name(), CitedModuleType::ConstraintGenerator,
		"Vikram K. Mulligan", "Systems Biology, Center for Computational Biology, Flatiron Institute",
		"vmulligan@flatironinstitute.org",
		"Integrated trRosetta into Rosetta, and wrote the trRosettaConstraintGenerator."
		)
	);

	citations.add( protocols::trRosetta::trRosettaProtocol_v1::get_trRosetta_neural_net_citation() );
}


#ifdef USE_TENSORFLOW

/// @brief Set the name of the multiple sequence alignment file.
/// @details Calls reset().
void
trRosettaConstraintGenerator::set_msa_file(
	std::string const & filename_in
) {
	msa_file_ = filename_in;
	reset();
}

/// @brief Reset (delete) the trRosetta neural network output.
/// @details Forces a re-run of the trRosetta neural network when
/// the apply() function is called.
void
trRosettaConstraintGenerator::reset() {
	trRosetta_outputs_ = nullptr;
}

/// @brief Set whether this will generate distance constraints.
/// @details Distance constraints constrain the distance between pairs of amino acids.
/// These are symmetric, and are only generated once per amino acid pair (since
/// dist(a,b) == dist(b,a)).
void
trRosettaConstraintGenerator::set_generate_dist_constraints(
	bool const setting
) {
	TR << "Configuring trRosettaConstraintGenerator " << (setting ? " to set " : " not to set ") << " distance constraints." << std::endl;
	generate_dist_constraints_ = setting;
}

/// @brief Set whether this will generate omega dihedral constraints.
/// @details Omega constraints constrain the dihedral between CA1-CB1-CB2-CA2 in pairs
/// of amino acids.  These are symmetric, and are only generated once per amino acid
/// pair (since omega(a,b) == omega(b,a)).
/// @note This is NOT omega the backbone dihedral torsion!
void
trRosettaConstraintGenerator::set_generate_omega_constraints(
	bool const setting
) {
	TR << "Configuring trRosettaConstraintGenerator " << (setting ? " to set " : " not to set ") << " omega dihedral constraints." << std::endl;
	generate_omega_constraints_ = setting;
}

/// @brief Set whether this will generate theta dihedral constraints.
/// @details Theta constraints constrain the dihedral between N1-CA2-CB1-CB2 in pairs
/// of amino acids.  These are asymmetric (i.e. theta(a,b)!=theta(b,a)), so there are
/// two per amino acid pair (unless a == b, which is skipped).
void
trRosettaConstraintGenerator::set_generate_theta_constraints(
	bool const setting
) {
	TR << "Configuring trRosettaConstraintGenerator " << (setting ? " to set " : " not to set ") << " theta dihedral constraints." << std::endl;
	generate_theta_constraints_ = setting;
}

/// @brief Set whether this will generate phi angle constraints.
/// @details Phi constraints constrain the angle between CA1-CB1-CB2 in pairs
/// of amino acids.  These are asymmetric (i.e. phi(a,b)!=phi(b,a)), so there are
/// two per amino acid pair (unless a == b, which is skipped).
/// @note This is NOT phi the backbone dihedral torsion!
void
trRosettaConstraintGenerator::set_generate_phi_constraints(
	bool const setting
) {
	TR << "Configuring trRosettaConstraintGenerator " << (setting ? " to set " : " not to set ") << " phi angle constraints." << std::endl;
	generate_phi_constraints_ = setting;
}

/// @brief Set the probability cutoffs for distance, omega, theta, and phi.
void
trRosettaConstraintGenerator::set_prob_cutoffs(
	core::Real const dist_cutoff,
	core::Real const omega_cutoff,
	core::Real const theta_cutoff,
	core::Real const phi_cutoff
) {
	std::string const errmsg("Error in trRosettaConstraintGenerator::set_prob_cutoffs(): ");
	runtime_assert_string_msg( dist_cutoff <=1.0 && dist_cutoff >= 0.0, errmsg + "The distance probability cutoff must be between 0.0 and 1.0!" );
	runtime_assert_string_msg( omega_cutoff <=1.0 && omega_cutoff >= 0.0, errmsg + "The omega dihedral probability cutoff must be between 0.0 and 1.0!" );
	runtime_assert_string_msg( theta_cutoff <=1.0 && theta_cutoff >= 0.0, errmsg + "The theta dihedral probability cutoff must be between 0.0 and 1.0!" );
	runtime_assert_string_msg( phi_cutoff <=1.0 && phi_cutoff >= 0.0, errmsg + "The phi angle probability cutoff must be between 0.0 and 1.0!" );
	dist_prob_cutoff_ = dist_cutoff;
	omega_prob_cutoff_ = omega_cutoff;
	theta_prob_cutoff_ = theta_cutoff;
	phi_prob_cutoff_ = phi_cutoff;
}

/// @brief Set the constraint weights for distance, omega, theta, and phi.
void
trRosettaConstraintGenerator::set_constraint_weights(
	core::Real const dist_weight,
	core::Real const omega_weight,
	core::Real const theta_weight,
	core::Real const phi_weight
) {
	std::string const errmsg( "Error in trRosettaConstraintGenerator::set_constraint_weights(): " );
	runtime_assert_string_msg( dist_weight >= 0.0, "Distance constraint weight must be positive!" );
	runtime_assert_string_msg( omega_weight >= 0.0, "Omega dihedral constraint weight must be positive!" );
	runtime_assert_string_msg( theta_weight >= 0.0, "Theta dihedral constraint weight must be positive!" );
	runtime_assert_string_msg( phi_weight >= 0.0, "Phi angle constraint weight must be positive!" );
	distance_constraint_weight_= dist_weight;
	omega_constraint_weight_ = omega_weight;
	theta_constraint_weight_ = theta_weight;
	phi_constraint_weight_ = phi_weight;
}

#endif //USE_TENSORFLOW

////////////////////////////////////////
// PRIVATE MEMBER FUNCTIONS
////////////////////////////////////////

#ifdef USE_TENSORFLOW

/// @brief Initialize this constraint generator from an options collection.
/// @details Called once by constructor.
void
trRosettaConstraintGenerator::init_from_commandline(
	utility::options::OptionCollection const & opts
) {
	set_msa_file( opts[basic::options::OptionKeys::trRosetta::msa_file]() );
	set_prob_cutoffs(
		opts[basic::options::OptionKeys::trRosetta::distance_constraint_prob_cutoff](),
		opts[basic::options::OptionKeys::trRosetta::omega_constraint_prob_cutoff](),
		opts[basic::options::OptionKeys::trRosetta::theta_constraint_prob_cutoff](),
		opts[basic::options::OptionKeys::trRosetta::phi_constraint_prob_cutoff]()
	);
	set_constraint_weights(
		opts[ basic::options::OptionKeys::trRosetta::distance_constraint_weight ](),
		opts[ basic::options::OptionKeys::trRosetta::omega_constraint_weight ](),
		opts[ basic::options::OptionKeys::trRosetta::theta_constraint_weight ](),
		opts[ basic::options::OptionKeys::trRosetta::phi_constraint_weight ]()
	);
	set_generate_dist_constraints( opts[basic::options::OptionKeys::trRosetta::use_distance_constraints]() );
	set_generate_omega_constraints( opts[basic::options::OptionKeys::trRosetta::use_omega_constraints]() );
	set_generate_theta_constraints( opts[basic::options::OptionKeys::trRosetta::use_theta_constraints]() );
	set_generate_phi_constraints( opts[basic::options::OptionKeys::trRosetta::use_phi_constraints]() );
}

/// @brief Load the multiple sequence alignment file and run the neural network.
/// @details Assumes that trRosetta_outputs_ is nullptr; throws if msa_file_ is empty.
/// @note TRIGGERS READ FROM DISK.
void
trRosettaConstraintGenerator::load_msa_and_run_neural_network() const {
	using namespace protocols::trRosetta;

	//Initial checks:
	std::string const errmsg( "Error in trRosettaConstraintGenerator::load_msa_and_run_neural_network(): " );
	debug_assert( trRosetta_outputs_ == nullptr ); //Should be true.
	runtime_assert_string_msg( !msa_file_.empty(), errmsg + "No multiple sequence alignment file was provided.  A "
		"multiple sequence alignment file must be specified at the commandline with the \"-trRosetta:msa_file\" option, "
		"in RosettaScripts XML with the \"msa_file\" option, or from Python or C++ by calling the "
		"\"trRosettaConstraintGenerator::set_msa_file()\" function."
	);
	runtime_assert_string_msg( utility::file::file_exists( msa_file_ ), errmsg + "The file \"" + msa_file_
		+ "\" could not be found!"
	);

	//Initialize the model -- TRIGGERS READ FROM DISK:
	TR << "Loading trRosetta model version " << trRosetta_model_version_ << "." << std::endl;
	trRosettaProtocolBaseOP model;
	switch( trRosetta_model_version_ ) {
		case 1:
			//Model version 1.
			model = utility::pointer::make_shared<trRosettaProtocol_v1 >();
			break;
		default:
			utility_exit_with_message( errmsg + "Model version was set to " + std::to_string( trRosetta_model_version_ )
				+ ", but this value is not a valid trRosetta model version number!"
			);
	};
	model->set_input_msa_file(msa_file_);

	//Run the model and store the output:
	trRosetta_outputs_ = model->run();
}

/// @brief Given a pose, generate the constraints.
/// @details Assumes that trRosetta_outputs_ is not a nullptr.  Throws if it is, in debug mode.
core::scoring::constraints::ConstraintCOPs
trRosettaConstraintGenerator::generate_constraints(
	core::pose::Pose const & pose
) const {
	debug_assert( trRosetta_outputs_ != nullptr );

	switch( trRosetta_model_version_ ) {
		case 1:
			return generate_constraints_v1( pose );
		default:
			utility_exit_with_message( "Error in trRosettaConstraintGenerator::generate_constraints(): Model version was set to "
				+ std::to_string( trRosetta_model_version_ )
				+ ", but this value is not a valid trRosetta model version number!"
			);
	};

	return utility::vector1< core::scoring::constraints::ConstraintCOP >(); //Keep the compiler happy.
}

/// @brief Given a pose, generate the constraints using the version 1 model.
/// @details Assumes that trRosetta_outputs_ is not a nullptr.  Throws if it is, in debug mode.
core::scoring::constraints::ConstraintCOPs
trRosettaConstraintGenerator::generate_constraints_v1(
	core::pose::Pose const & pose
) const {
	using namespace protocols::trRosetta;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	utility::vector1< core::scoring::constraints::ConstraintCOP > outputvec;

	core::Size const seqlength( pose.total_residue() );

	//Parameters for v1 model:
	core::Size const n_dist_bins_model(37); //In neural net.
	core::Size const n_omega_bins_model(25); //In neural net.
	core::Size const n_theta_bins_model(25); //In neural net.
	core::Size const n_phi_bins_model(13); //In neural net.
	core::Size const n_dist_bins(35); //In constraints spline.
	core::Size const n_omega_bins(28); //In constraints spline.
	core::Size const n_theta_bins(28); //In constraints spline.
	core::Size const n_phi_bins(16); //In constraints spline.
	core::Real const dist_bins_denominator( 19.5 );
	core::Real const dist_bins_exponent( 1.57 );
	utility::vector1< core::Real > dist_bins_vect(n_dist_bins);
	utility::vector1< core::Real > dist_bins_background( n_dist_bins - 3 );
	utility::vector1< core::Real > omega_bins_vect(n_omega_bins);
	utility::vector1< core::Real > phi_bins_vect( n_phi_bins );
	dist_bins_vect[1] = 0.0;
	dist_bins_vect[2] = 2.0;
	dist_bins_vect[3] = 3.5;
	for( core::Size ii(4); ii<=n_dist_bins; ++ii ) {
		dist_bins_vect[ii] = 4.25 + static_cast<core::Real>(ii-4)*0.5;
		dist_bins_background[ii-3] = std::pow( dist_bins_vect[ii]/dist_bins_denominator, dist_bins_exponent );
	}
	for( core::Size ii(1); ii<=n_omega_bins; ++ii ) {
		omega_bins_vect[ii] = (1.125*numeric::constants::d::pi)*( -1 + 2*(ii-1)/( static_cast<core::Real>(n_omega_bins-1) ) );
	}
	for( core::Size ii(1); ii<=n_phi_bins; ++ii ) {
		phi_bins_vect[ii] = ( (static_cast<core::Real>(ii-1)/static_cast<core::Real>(n_phi_bins-1)) * (1.0 + 3.0*(15.0/180.0)) - 1.5 * 15.0/180.0 ) * numeric::constants::d::pi;
	}
	utility::vector1< core::Real > const & theta_bins_vect(omega_bins_vect);

#ifdef NDEBUG
	trRosettaOutputs_v1COP outputs( utility::pointer::static_pointer_cast< trRosettaOutputs_v1 const >( trRosetta_outputs_ ) );
#else
	trRosettaOutputs_v1COP outputs( utility::pointer::dynamic_pointer_cast< trRosettaOutputs_v1 const >( trRosetta_outputs_ ) );
	debug_assert( outputs != nullptr );
#endif

	bool const pose_is_centroid( pose.residue_type_set_for_pose()->mode() == core::chemical::CENTROID_t );

	//Looping over every pair of residues
	for( core::Size ires(1); ires<seqlength; ++ires ) {
		core::chemical::ResidueType const & restype_i( pose.residue_type(ires) );
		if( !restype_i.is_canonical_aa() ) continue;
		bool res_i_is_gly( restype_i.aa() == core::chemical::aa_gly );
		core::id::AtomID const cb_atom_i( restype_i.atom_index( res_i_is_gly ? "CA" : "CB" ), ires );
		core::id::AtomID const ca_atom_i( restype_i.atom_index( "CA" ), ires );
		core::id::AtomID const n_atom_i( restype_i.atom_index( "N" ), ires );
		core::id::AtomID const dihedral_cb_atom_i( restype_i.atom_index( res_i_is_gly ? ( pose_is_centroid ? "CEN" : "1HA" ) : "CB" ), ires );

		for( core::Size jres(ires+1); jres<=seqlength; ++jres ) {
			core::chemical::ResidueType const & restype_j( pose.residue_type(jres) );
			if( !restype_j.is_canonical_aa() ) continue;
			bool res_j_is_gly( restype_j.aa() == core::chemical::aa_gly );
			core::id::AtomID const cb_atom_j( restype_j.atom_index( res_j_is_gly ? "CA" : "CB" ), jres );
			core::id::AtomID const ca_atom_j( restype_j.atom_index( "CA" ), jres );
			core::id::AtomID const n_atom_j( restype_j.atom_index( "N" ), jres );
			core::id::AtomID const dihedral_cb_atom_j( restype_j.atom_index( res_j_is_gly ? ( pose_is_centroid ? "CEN" : "1HA" ) : "CB" ), jres );

			////////////////////
			// DISTANCE:
			////////////////////
			if( generate_dist_constraints_ ) {
				generate_dist_constraints_v1(
					outputs, outputvec, dist_bins_background, dist_bins_vect,
					cb_atom_i, cb_atom_j,
					n_dist_bins, n_dist_bins_model, ires, jres, dist_prob_cutoff()
				);
			}
			////////////////////
			// END DISTANCE
			////////////////////

			////////////////////
			// OMEGA (DIHEDRAL):
			////////////////////
			if( generate_omega_constraints_ ) {
				generate_omega_or_theta_constraints_v1(
					false,
					outputs, outputvec, omega_bins_vect,
					ca_atom_i, dihedral_cb_atom_i, dihedral_cb_atom_j, ca_atom_j,
					n_omega_bins, n_omega_bins_model, ires, jres, omega_prob_cutoff()
				);
			}
			////////////////////
			// END OMEGA
			// (DIHEDRAL)
			////////////////////

			////////////////////
			// THETA (DIHEDRAL)
			////////////////////
			if( generate_theta_constraints_ ) {
				generate_omega_or_theta_constraints_v1(
					true,
					outputs, outputvec, theta_bins_vect,
					n_atom_i, ca_atom_i, dihedral_cb_atom_i, dihedral_cb_atom_j,
					n_theta_bins, n_theta_bins_model, ires, jres, theta_prob_cutoff()
				);
				generate_omega_or_theta_constraints_v1(
					true,
					outputs, outputvec, theta_bins_vect,
					n_atom_j, ca_atom_j, dihedral_cb_atom_j, dihedral_cb_atom_i,
					n_theta_bins, n_theta_bins_model, jres, ires, theta_prob_cutoff()
				);
			}
			////////////////////
			// END THETA
			// (DIHEDRAL)
			////////////////////

			////////////////////
			// PHI (ANGLE)
			////////////////////
			if( generate_phi_constraints_ ) {
				generate_phi_constraints_v1(
					outputs, outputvec, phi_bins_vect,
					ca_atom_i, dihedral_cb_atom_i, dihedral_cb_atom_j,
					n_phi_bins, n_phi_bins_model, ires, jres, phi_prob_cutoff()
				);
				generate_phi_constraints_v1(
					outputs, outputvec, phi_bins_vect,
					ca_atom_j, dihedral_cb_atom_j, dihedral_cb_atom_i,
					n_phi_bins, n_phi_bins_model, jres, ires, phi_prob_cutoff()
				);
			}
			////////////////////
			// END PHI
			// (ANGLE)
			////////////////////

		}
	} //End loop over every pair of residues.

	return outputvec;
}

/// @brief Generate the distance constraints using the version 1 model.
/// @details New constraints are appended to outputvec.  Outputvec is not cleared
/// by this operation.
void
trRosettaConstraintGenerator::generate_dist_constraints_v1(
	protocols::trRosetta::trRosettaOutputs_v1COP trRosetta_outputs,
	utility::vector1 <core::scoring::constraints::ConstraintCOP> & outputvec,
	utility::vector1< core::Real > const & dist_bins_background,
	utility::vector1< core::Real > const & dist_bins_vect,
	core::id::AtomID const & cb_atom_i,
	core::id::AtomID const & cb_atom_j,
	core::Size const n_dist_bins,
	core::Size const n_dist_bins_model,
	core::Size const ires,
	core::Size const jres,
	core::Real const prob_cutoff
) const {
	using namespace protocols::trRosetta;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	core::Real dist_probsum(0.0);
	for( core::Size ibin(6); ibin<=n_dist_bins_model; ++ibin ) {
		dist_probsum += trRosetta_outputs->dist(ires, jres, ibin);
	}

	if( dist_probsum <= prob_cutoff ) return; //Can stop here if below probability cutoff.

	if( TR.Debug.visible() ) {
		TR.Debug << "Generating distance constraints between residues " << ires << " and " << jres << "." << std::endl;
	}

	utility::vector1< core::Real > dist_attractive_repulsive(n_dist_bins);
	for( core::Size ibin(6); ibin<=n_dist_bins_model; ++ibin ) {
		dist_attractive_repulsive[ibin-2] = -std::log( (trRosetta_outputs->dist(ires,jres,ibin) + 0.0001) / (trRosetta_outputs->dist(ires,jres,n_dist_bins_model) * dist_bins_background[ibin-5] ) ) - 0.5;
	}
	core::Real const addval( std::max(0.0, dist_attractive_repulsive[4]) );
	dist_attractive_repulsive[1] = addval + 10.0;
	dist_attractive_repulsive[2] = addval + 3.0;
	dist_attractive_repulsive[3] = addval + 0.5;

	if( TR.Debug.visible() ) {
		TR.Debug << "x_axis\t";
		for( core::Size ii(1); ii<=dist_bins_vect.size(); ++ii) {
			TR.Debug << dist_bins_vect[ii] << "\t";
		}
		TR.Debug << std::endl;
		TR.Debug << "y_axis\t";
		for( core::Size ii(1); ii<=dist_attractive_repulsive.size(); ++ii) {
			TR.Debug << dist_attractive_repulsive[ii] << "\t";
		}
		TR.Debug << std::endl;
	}

	SplineFuncOP splinefunc(
		utility::pointer::make_shared< SplineFunc >(
			"TAG", distance_constraint_weight(), 1.0, 0.5,
			dist_bins_vect, dist_attractive_repulsive
		)
	);
	outputvec.push_back( utility::pointer::make_shared< AtomPairConstraint >( cb_atom_i, cb_atom_j, splinefunc ) );
}

/// @brief Generate the omega or theta dihedral constraints using the version 1 model.
/// @details New constraints are appended to outputvec.  Outputvec is not cleared
/// by this operation.
void
trRosettaConstraintGenerator::generate_omega_or_theta_constraints_v1(
	bool const generate_theta_constraints,
	protocols::trRosetta::trRosettaOutputs_v1COP trRosetta_outputs,
	utility::vector1 <core::scoring::constraints::ConstraintCOP> & outputvec,
	utility::vector1< core::Real > const & bins_vect,
	core::id::AtomID const & at1,
	core::id::AtomID const & at2,
	core::id::AtomID const & at3,
	core::id::AtomID const & at4,
	core::Size const n_bins,
	core::Size const n_bins_model,
	core::Size const ires,
	core::Size const jres,
	core::Real const prob_cutoff
) const {
	using namespace protocols::trRosetta;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	core::Real probsum(0.0);
	for( core::Size ibin(2); ibin<=n_bins_model; ++ibin ) {
		probsum += (generate_theta_constraints ? trRosetta_outputs->theta(ires, jres, ibin) : trRosetta_outputs->omega(ires, jres, ibin));
	}
	if( probsum <= prob_cutoff ) return; //Skip omega or theta constraints if below probability cutoff.

	if( TR.Debug.visible() ) {
		TR.Debug << "Generating " << ( generate_theta_constraints ? "theta" : "omega" ) << " (inter-residue dihedral) constraints between residues " << ires << " and " << jres << "." << std::endl;
	}
	utility::vector1< core::Real > logvect( n_bins );
	core::Real const lastbin( ( generate_theta_constraints ? trRosetta_outputs->theta(ires, jres, n_bins_model) : trRosetta_outputs->omega(ires, jres, n_bins_model) ) + 0.0001 );
	for( core::Size ibin(1); ibin<=n_bins_model; ++ibin ) {
		logvect[ibin+1] = -std::log( ( ( generate_theta_constraints ? trRosetta_outputs->theta(ires, jres, ibin ) : trRosetta_outputs->omega(ires, jres, ibin ) ) + 0.0001 ) / ( lastbin ) );
	}
	logvect[n_bins_model+2] = logvect[3];
	logvect[n_bins_model+3] = logvect[4];
	logvect[1] = logvect[n_bins_model];
	logvect[2] = logvect[n_bins_model+1];
	if( TR.Debug.visible() ) {
		TR.Debug << "x_axis\t";
		for( core::Size ii(1); ii<=bins_vect.size(); ++ii) {
			TR.Debug << bins_vect[ii] << "\t";
		}
		TR.Debug << std::endl;
		TR.Debug << "y_axis\t";
		for( core::Size ii(1); ii<=logvect.size(); ++ii) {
			TR.Debug << logvect[ii] << "\t";
		}
		TR.Debug << std::endl;
	}

	SplineFuncOP splinefunc(
		utility::pointer::make_shared< SplineFunc >(
			"TAG",
			generate_theta_constraints ? theta_constraint_weight() : omega_constraint_weight(),
			1.0, bins_vect[2]-bins_vect[1],
			bins_vect, logvect
		)
	);

	outputvec.push_back(
		utility::pointer::make_shared< DihedralConstraint >(
			at1, at2, at3, at4,
			splinefunc
		)
	);
}

/// @brief Generate the phi angle constraints using the version 1 model.
/// @details New constraints are appended to outputvec.  Outputvec is not cleared
/// by this operation.
void
trRosettaConstraintGenerator::generate_phi_constraints_v1(
	protocols::trRosetta::trRosettaOutputs_v1COP trRosetta_outputs,
	utility::vector1 <core::scoring::constraints::ConstraintCOP> & outputvec,
	utility::vector1< core::Real > const & phi_bins_vect,
	core::id::AtomID const & at1,
	core::id::AtomID const & at2,
	core::id::AtomID const & at3,
	core::Size const n_phi_bins,
	core::Size const n_phi_bins_model,
	core::Size const ires,
	core::Size const jres,
	core::Real const prob_cutoff
) const {
	using namespace protocols::trRosetta;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	core::Real probsum(0.0);
	for( core::Size ibin(2); ibin <= n_phi_bins_model; ++ibin ) {
		probsum += trRosetta_outputs->phi(ires,jres,ibin);
	}
	if( probsum <= prob_cutoff ) return;

	utility::vector1< core::Real > logvect( n_phi_bins );
	core::Real const denom( trRosetta_outputs->phi(ires,jres,n_phi_bins_model) + 0.0001 );
	for( core::Size ibin(2); ibin<=n_phi_bins_model; ++ibin ) {
		logvect[ibin + 1] = -std::log( (trRosetta_outputs->phi(ires,jres,ibin) + 0.0001) / denom );
	}
	logvect[n_phi_bins_model+2] = logvect[n_phi_bins_model+1];
	logvect[n_phi_bins_model+3] = logvect[n_phi_bins_model];
	logvect[1] = logvect[4];
	logvect[2] = logvect[3];

	if( TR.Debug.visible() ) {
		TR.Debug << "x_axis\t";
		for( core::Size ii(1); ii<=phi_bins_vect.size(); ++ii) {
			TR.Debug << phi_bins_vect[ii] << "\t";
		}
		TR.Debug << std::endl;
		TR.Debug << "y_axis\t";
		for( core::Size ii(1); ii<=logvect.size(); ++ii) {
			TR.Debug << logvect[ii] << "\t";
		}
		TR.Debug << std::endl;
	}

	SplineFuncOP splinefunc(
		utility::pointer::make_shared< SplineFunc >(
			"TAG", phi_constraint_weight(), 1.0, phi_bins_vect[2]-phi_bins_vect[1],
			phi_bins_vect, logvect
		)
	);

	outputvec.push_back(
		utility::pointer::make_shared< AngleConstraint >(
			at1, at2, at3,
			splinefunc
		)
	);
}

#endif //USE_TENSORFLOW

} //constraint_generators
} //trRosetta_protocols
} //protocols
