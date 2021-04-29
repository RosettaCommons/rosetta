// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file trRosetta_protocols/movers/trRosettaProtocolMover.cc
/// @brief The full trRosetta structure prediction protocol from Yang et al, converted to
/// C++ and implemented as a mover.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Unit headers
#include <protocols/trRosetta_protocols/movers/trRosettaProtocolMover.hh>
#include <protocols/trRosetta_protocols/movers/trRosettaProtocolMoverCreator.hh>

// Core headers
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AA.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/pose/variant_util.hh>
#include <core/scoring/RamaPrePro.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/simple_metrics/metrics/RMSDMetric.hh>
#include <core/simple_metrics/metrics/TotalEnergyMetric.hh>
#include <core/simple_metrics/metrics/TimingProfileMetric.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// Protocols headers
#include <protocols/trRosetta/trRosettaProtocol_v1.hh>
#include <protocols/trRosetta_protocols/constraint_generators/trRosettaConstraintGenerator.hh>
#include <protocols/jd2/util.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/residue_selectors/StoredResidueSubsetSelector.hh>
#include <protocols/residue_selectors/StoreResidueSubsetMover.hh>
#include <protocols/simple_moves/bin_transitions/InitializeByBins.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/relax/FastRelax.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/trRosetta.OptionKeys.gen.hh>
#include <basic/tensorflow_manager/util.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// Numeric headers
#include <numeric/random/random.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <utility/vector1.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/CitationManager.hh>

static basic::Tracer TR( "protocols.trRosetta_protocols.movers.trRosettaProtocolMover" );

namespace protocols {
namespace trRosetta_protocols {
namespace movers {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
/// @details Initializes from the global options system.
trRosettaProtocolMover::trRosettaProtocolMover():
	protocols::moves::Mover( trRosettaProtocolMover::mover_name() )
{
#ifdef USE_TENSORFLOW
	init_from_options( basic::options::option ); //Init from global options.
#endif //USE_TENSORFLOW
}

/// @brief Options constructor.
/// @details Initializes from a local options collection.
trRosettaProtocolMover::trRosettaProtocolMover(
#ifdef USE_TENSORFLOW
	utility::options::OptionCollection const & options
#else // !USE_TENSORFLOW
	utility::options::OptionCollection const &
#endif // USE_TENSORFLOW
) :
	protocols::moves::Mover( trRosettaProtocolMover::mover_name() )
{
#ifdef USE_TENSORFLOW
	init_from_options(options);
#endif //USE_TENSORFLOW
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
trRosettaProtocolMover::~trRosettaProtocolMover() = default;

/// @brief Indicate commandline flags that are relevant to this mover.
void
trRosettaProtocolMover::register_options( bool const only_constraint_generator_options /*= false*/ ) {
	using namespace basic::options;
	option.add_relevant( basic::options::OptionKeys::in::file::fasta );
	option.add_relevant( basic::options::OptionKeys::in::file::native );
	option.add_relevant( basic::options::OptionKeys::trRosetta::msa_file );

	option.add_relevant( basic::options::OptionKeys::trRosetta::use_distance_constraints );
	option.add_relevant( basic::options::OptionKeys::trRosetta::use_omega_constraints );
	option.add_relevant( basic::options::OptionKeys::trRosetta::use_theta_constraints );
	option.add_relevant( basic::options::OptionKeys::trRosetta::use_phi_constraints );

	option.add_relevant( basic::options::OptionKeys::trRosetta::distance_constraint_prob_cutoff );
	option.add_relevant( basic::options::OptionKeys::trRosetta::omega_constraint_prob_cutoff );
	option.add_relevant( basic::options::OptionKeys::trRosetta::theta_constraint_prob_cutoff );
	option.add_relevant( basic::options::OptionKeys::trRosetta::phi_constraint_prob_cutoff );

	option.add_relevant( basic::options::OptionKeys::trRosetta::distance_constraint_weight );
	option.add_relevant( basic::options::OptionKeys::trRosetta::omega_constraint_weight );
	option.add_relevant( basic::options::OptionKeys::trRosetta::theta_constraint_weight );
	option.add_relevant( basic::options::OptionKeys::trRosetta::phi_constraint_weight );

	if ( !only_constraint_generator_options ) {
		option.add_relevant( basic::options::OptionKeys::trRosetta::backbone_randomization_mode );
		option.add_relevant( basic::options::OptionKeys::trRosetta::backbone_minimization_mode );

		option.add_relevant( basic::options::OptionKeys::trRosetta::cis_peptide_prob_non_prepro );
		option.add_relevant( basic::options::OptionKeys::trRosetta::cis_peptide_prob_prepro );
		option.add_relevant( basic::options::OptionKeys::trRosetta::scorefxn0 );
		option.add_relevant( basic::options::OptionKeys::trRosetta::scorefxn1 );
		option.add_relevant( basic::options::OptionKeys::trRosetta::scorefxn2 );
		option.add_relevant( basic::options::OptionKeys::trRosetta::scorefxn3 );
		option.add_relevant( basic::options::OptionKeys::trRosetta::scorefxn_fullatom );

		option.add_relevant( basic::options::OptionKeys::trRosetta::mutate_gly_to_ala );
		option.add_relevant( basic::options::OptionKeys::trRosetta::fullatom_refinement );
	}
}

/// @brief Given a backbone randomization mode enum, get a string corresponding to the enum.
/// @details This function must be updated if entries are added to the trRosettaProtocolBackboneRandomizationMode enum class!
/*static*/
std::string
trRosettaProtocolMover::get_randomization_mode_string_from_enum(
	trRosettaProtocolBackboneRandomizationMode const mode_enum
) {
	switch( mode_enum ) {
	case trRosettaProtocolBackboneRandomizationMode::classic :
		return "classic";
	case trRosettaProtocolBackboneRandomizationMode::ramachandran :
		return "ramachandran";
	case trRosettaProtocolBackboneRandomizationMode::bins :
		return "bins";
	default :
		utility_exit_with_message( "Error in trRosettaProtocolMover::get_randomization_mode_string_from_enum(): Invalid enum value \"" + std::to_string( static_cast< int >(mode_enum) ) + "\" provided to this function!" );

	};
	return ""; //We'll never reach here, but this keeps the compiler happy.
}

/// @brief Given a backbone randomization mode name, get the enum corresponding to the name.
/// @details Throws if an invalid name is provided.
/*static*/
trRosettaProtocolBackboneRandomizationMode
trRosettaProtocolMover::get_randomization_mode_enum_from_string(
	std::string const & mode_name
) {
	for ( core::Size i(1); i<=static_cast<core::Size>(trRosettaProtocolBackboneRandomizationMode::NUM_ENTRIES); ++i ) {
		if ( mode_name == get_randomization_mode_string_from_enum( static_cast< trRosettaProtocolBackboneRandomizationMode >(i) ) ) {
			return static_cast< trRosettaProtocolBackboneRandomizationMode >(i);
		}
	}
	utility_exit_with_message( "Error in trRosettaProtocolMover::get_randomization_mode_enum_from_string(): Could not parse \"" + mode_name + "\" as a valid backbone randomization mode for the trRosettaProtocol mover." );
	return trRosettaProtocolBackboneRandomizationMode::classic; //Keep the compiler happy.
}

/// @brief Given a backbone minimization mode enum, get a string corresponding to the enum.
/// @details This function must be updated if entries are added to the trRosettaProtocolBackboneMinimizationMode enum class!
/*static*/
std::string
trRosettaProtocolMover::get_minimization_mode_string_from_enum(
	trRosettaProtocolBackboneMinimizationMode const mode_enum
) {
	switch( mode_enum ) {
	case trRosettaProtocolBackboneMinimizationMode::classic0 :
		return "classic0";
	case trRosettaProtocolBackboneMinimizationMode::classic1 :
		return "classic1";
	case trRosettaProtocolBackboneMinimizationMode::classic2 :
		return "classic2";
	default :
		utility_exit_with_message( "Error in trRosettaProtocolMover::get_minimization_mode_string_from_enum(): Invalid enum value \"" + std::to_string( static_cast< int >(mode_enum) ) + "\" provided to this function!" );
	};
	return ""; //We'll never reach here, but this keeps the compiler happy.
}

/// @brief Given a backbone minimization mode name, get the enum corresponding to the name.
/// @details Throws if an invalid name is provided.
/*static*/
trRosettaProtocolBackboneMinimizationMode
trRosettaProtocolMover::get_minimization_mode_enum_from_string(
	std::string const & mode_name
) {
	for ( core::Size i(1); i<=static_cast<core::Size>(trRosettaProtocolBackboneMinimizationMode::NUM_ENTRIES); ++i ) {
		if ( mode_name == get_minimization_mode_string_from_enum( static_cast< trRosettaProtocolBackboneMinimizationMode >(i) ) ) {
			return static_cast< trRosettaProtocolBackboneMinimizationMode >(i);
		}
	}
	utility_exit_with_message( "Error in trRosettaProtocolMover::get_minimization_mode_enum_from_string(): Could not parse \"" + mode_name + "\" as a valid backbone minimization mode for the trRosettaProtocol mover." );
	return trRosettaProtocolBackboneMinimizationMode::classic2; //Keep the compiler happy.
}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
trRosettaProtocolMover::apply(
#ifdef USE_TENSORFLOW
	core::pose::Pose & pose
#else // !USE_TENSORFLOW
	core::pose::Pose &
#endif // USE_TENSORFLOW
) {
#ifdef USE_TENSORFLOW

	using namespace protocols::minimization_packing;

	TR << "Starting trRosetta structure prediction protocol." << std::endl;

	core::simple_metrics::metrics::TimingProfileMetric timing_metric_centroid;

	std::string const errmsg( "Error in trRosettaProtocolMover::apply(): " );
	runtime_assert_string_msg( !msa_file_.empty(), errmsg + "A multiple sequence alignment file must be set either "
		"at the commandline with the -trRosetta:msa_file option, from XML with the msa_file option, "
		"or from C++ or Python code with the trRosettaProtocolMover::set_msa_file() method."
	);

	//Get the sequence:
	if ( sequence().empty() ) {
		runtime_assert_string_msg( !(fasta_file().empty()), errmsg + "No sequence or FASTA file was provided to the trRosettaProtocol mover!" );
		utility::vector1< core::sequence::SequenceOP > const seqs( core::sequence::read_fasta_file( fasta_file() ) );
		runtime_assert_string_msg( seqs.size() == 1, errmsg + "The FASTA file \"" + fasta_file() + "\" contained " + std::to_string(seqs.size()) + " sequences.  Exactly one sequence must be supplied to the trRosettaProtocol mover." );
		sequence_ /*mutable*/ = seqs[1]->ungapped_sequence();
	}
	validate_sequence();

	if ( TR.Warning.visible() && !( use_distance_constraints() || use_omega_constraints() || use_theta_constraints() || use_phi_constraints() ) ) {
		TR.Warning << "The trRosettaProtocol mover is set to use none of the constraints from the trRosetta constraint generator.  ";
		TR.Warning << "Unconstrained minimization from a random starting conformation is not expected to yield a meaningful structure!" << std::endl;
	}

	//Create MinMovers that we'll use:
	core::kinematics::MoveMapOP movemap( utility::pointer::make_shared< core::kinematics::MoveMap >() );
	movemap->set_bb(true);
	movemap->set_chi(false);
	movemap->set_jump(true); //Should really be false, but it doesn't matter.  There are no jumps, since below we clear the pose and build a new single-chain pose.
	MinMoverOP minmover0( configure_minmover( movemap, sfxn0_, 1000, false ) );
	MinMoverOP minmover1( configure_minmover( movemap, sfxn1_, 1000, false ) );
	MinMoverOP minmover2( configure_minmover( movemap, sfxn2_, 500, false ) );
	MinMoverOP minmover3( configure_minmover( movemap, sfxn3_, 1000, true ) );
	protocols::moves::RepeatMoverOP repeat_minmover0( utility::pointer::make_shared< protocols::moves::RepeatMover >(minmover0, 3) );

	//Clear the pose and create a new one:
	TR << "Building pose from sequence." << std::endl;
	pose.clear();
	pose = make_centroid_pose_from_sequence();

	//Mutate glycines to alanine (and get back a StoredResidueSubsetSelector that stores the gly positions, to allow them to be mutated back to gly later):
	core::select::residue_selector::ResidueSelectorCOP glycine_selector(
		mutate_gly_to_ala() ?
		do_mutate_gly_to_ala( pose ) :
		nullptr
	);

	//Randomize the backbone dihedrals:
	randomize_backbone_dihedrals(pose);

	//Remove clashes:
	remove_clashes(pose, *minmover2, *sfxn2_);

	// Generate the trRosetta constraints:
	utility::vector1< core::scoring::constraints::ConstraintCOP > trRosetta_constraints(
		generate_trRosetta_constraints( pose )
	);

	// Do the actual minimization protocol:
	TR << "Performing centroid-mode backbone optimization protocol \"" << get_backbone_minimization_mode_name() << "\"." << std::endl;
	switch( backbone_minimization_mode_ ) {
	case trRosettaProtocolBackboneMinimizationMode::classic0 :
		perform_classic0_minimization_protocol( pose, trRosetta_constraints, repeat_minmover0, minmover1, minmover3 );
		break;
	case trRosettaProtocolBackboneMinimizationMode::classic1 :
		perform_classic1_minimization_protocol( pose, trRosetta_constraints, repeat_minmover0, minmover1, minmover3 );
		break;
	case trRosettaProtocolBackboneMinimizationMode::classic2 :
		perform_classic2_minimization_protocol( pose, trRosetta_constraints, repeat_minmover0, minmover1, minmover3 );
		break;
	default :
		utility_exit_with_message(errmsg + "Unknown minimization mode specified!");
	};

	store_constraints_score(pose, "after_centroid_phase");
	store_rmsd( pose, get_native_pose(), "after_centroid_phase" );

	if ( mutate_gly_to_ala() ) {
		remove_trRosetta_constraints(trRosetta_constraints, pose);
		do_mutate_ala_to_gly( glycine_selector, pose );
	}
	(*sfxn1_)(pose);

	timing_metric_centroid.apply( "time_for_centroid_phase_in_minutes", pose );

	if ( fullatom_refinement() ) {
		core::simple_metrics::metrics::TimingProfileMetric timing_metric_fullatom;
		TR << "Carrying out fullatom refinement." << std::endl;
		convert_to_fullatom(pose);
		trRosetta_constraints = generate_trRosetta_constraints( pose, 0.1 );
		add_constraints_to_pose( pose, trRosetta_constraints, 1, pose.total_residue() );
		do_fullatom_refinement(pose);
		store_constraints_score(pose, "after_fullatom_refinement");
		store_rmsd( pose, get_native_pose(), "after_fullatom_refinement" );
		remove_trRosetta_constraints( trRosetta_constraints, pose );
		(*sfxn_fullatom_)(pose);
		timing_metric_fullatom.apply( "time_for_fullatom_phase_in_minutes", pose );
	}

	TR << "Completed trRosetta structure prediction protocol." << std::endl;
	TR.flush();
#else // !USE_TENSORFLOW
	utility_exit_with_message(
		"Error in trRosettaProtocolMover::apply():  This mover can only be used in the Tensorflow-enabled Rosetta compilation.\n\n"
		+ basic::tensorflow_manager::get_tensorflow_compilation_instructions( "trRosettaProtocolMover" )
	);
#endif //USE_TENSORFLOW
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
trRosettaProtocolMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// RosettaScripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
trRosettaProtocolMover::parse_my_tag(
#ifdef USE_TENSORFLOW
	utility::tag::TagCOP tag,
#else // !USE_TENSORFLOW
	utility::tag::TagCOP,
#endif // USE_TENSORFLOW
	basic::datacache::DataMap &
) {
	std::string const errmsg( "Error in trRosettaProtocolMover::parse_my_tag(): " );

#ifdef USE_TENSORFLOW

	runtime_assert_string_msg(
		!(tag->hasOption("sequence") && tag->hasOption("fasta_file")),
		errmsg + "One or the other of the \"sequence\" and \"fasta_file\" options may be provided.  Both were detected!"
	);

	set_fasta_file( tag->getOption<std::string>("fasta_file", fasta_file()) );
	set_sequence( tag->getOption<std::string>("sequence", sequence()) );
	set_msa_file( tag->getOption<std::string>("msa_file", msa_file()) );
	set_use_distance_constraints( tag->getOption<bool>("use_distance_constraints", use_distance_constraints()));
	set_use_omega_constraints( tag->getOption<bool>("use_omega_constraints", use_omega_constraints()));
	set_use_theta_constraints( tag->getOption<bool>("use_theta_constraints", use_theta_constraints()));
	set_use_phi_constraints( tag->getOption<bool>("use_phi_constraints", use_phi_constraints()));

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

	set_backbone_randomization_mode( tag->getOption<std::string>("backbone_randomization_mode", get_backbone_randomization_mode_name()) );
	set_ramachandran_mode_cis_probabilities(
		tag->getOption< core::Real >( "cis_peptide_prob_non_prepro", ramachandran_mode_cis_probability_non_prepro() ),
		tag->getOption< core::Real >( "cis_peptide_prob_prepro", ramachandran_mode_cis_probability_prepro() )
	);

	if ( tag->hasOption("scorefxn0") ) {
		set_scorefunction_for_minimization_stage( tag->getOption< std::string >("scorefxn0"), 0 );
	}
	if ( tag->hasOption("scorefxn1") ) {
		set_scorefunction_for_minimization_stage( tag->getOption< std::string >("scorefxn1"), 1 );
	}
	if ( tag->hasOption("scorefxn2") ) {
		set_scorefunction_for_minimization_stage( tag->getOption< std::string >("scorefxn2"), 2 );
	}
	if ( tag->hasOption("scorefxn3") ) {
		set_scorefunction_for_minimization_stage( tag->getOption< std::string >("scorefxn3"), 3 );
	}
	if ( tag->hasOption("scorefxn_fullatom") ) {
		set_scorefunction_for_fullatom_refinement( tag->getOption< std::string >("scorefxn_fullatom") );
	}

	set_fullatom_refinement( tag->getOption<bool>("fullatom_refinement", fullatom_refinement() ) );

	set_backbone_minimization_mode( tag->getOption<std::string>("backbone_minimization_mode", get_backbone_minimization_mode_name()) );
	set_mutate_gly_to_ala( tag->getOption<bool>("mutate_gly_to_ala", mutate_gly_to_ala() ) );

#else // !USE_TENSORFLOW
	utility_exit_with_message(
		errmsg + "The trRosettaProtocolMover cannot be used in this compilation of Rosetta.\n\n"
		+ basic::tensorflow_manager::get_tensorflow_compilation_instructions("trRosettaProtocolMover")
	);
#endif //USE_TENSORFLOW
}


void trRosettaProtocolMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default( "msa_file", xs_string,
		"Filename for a multiple sequence alignment file, in a3m format.  "
		"Dashes indicate gap sequences, and lowercase characters will be removed "
		"(and flanking regions ligated).  If not provided, the commandline option "
		"-trRosetta:msa_file will be used.  One or the other is required.", "" )

		+ XMLSchemaAttribute::attribute_w_default( "use_distance_constraints", xsct_rosetta_bool,
		"Set whether inter-residue distance constraints generated by the trRosetta neural "
		"network should be used for structure prediction.  True by default, unless a default is set "
		"at the commandline with the -trRosetta:use_distance_constraints flag.",
		"true" )

		+ XMLSchemaAttribute::attribute_w_default( "use_omega_constraints", xsct_rosetta_bool,
		"Set whether inter-residue omega dihedral constraints generated by the trRosetta neural "
		"network should be used for structure prediction.  Note that this is NOT the omega backbone "
		"dihedral angle, but an inter-residue dihedral defined by CA1-CB1-CB2-CA2.  True by default, "
		"unless a default is set at the commandline with the -trRosetta:use_omega_constraints flag.",
		"true" )

		+ XMLSchemaAttribute::attribute_w_default( "use_theta_constraints", xsct_rosetta_bool,
		"Set whether inter-residue theta dihedral constraints generated by the trRosetta neural "
		"network should be used for structure prediction.  Note that this is NOT a backbone "
		"dihedral angle, but an inter-residue dihedral defined by N1-CA1-CB1-CB2.  True by default, "
		"unless a default is set at the commandline with the -trRosetta:use_theta_constraints flag.",
		"true" )

		+ XMLSchemaAttribute::attribute_w_default( "use_phi_constraints", xsct_rosetta_bool,
		"Set whether inter-residue phi angle constraints generated by the trRosetta neural "
		"network should be used for structure prediction.  Note that this is NOT the phi backbone "
		"dihedral angle, but an inter-residue angle defined by CA1-CB1-CB2.  True by default, "
		"unless a default is set at the commandline with the -trRosetta:use_phi_constraints flag.",
		"true" )

		+ XMLSchemaAttribute::attribute_w_default( "distance_constraint_prob_cutoff", xsct_real,
		"Set the probability cutoff below which we omit a distance constraint.  Default 0.05, or "
		"whatever is set on the commandline with the -trRosetta::distance_constraint_prob_cutoff "
		"commandline option.",
		"0.05" )

		+ XMLSchemaAttribute::attribute_w_default( "omega_constraint_prob_cutoff", xsct_real,
		"Set the probability cutoff below which we omit a omega dihedral constraint.  Default 0.55, or "
		"whatever is set on the commandline with the -trRosetta::omega_constraint_prob_cutoff "
		"commandline option.",
		"0.55" )

		+ XMLSchemaAttribute::attribute_w_default( "theta_constraint_prob_cutoff", xsct_real,
		"Set the probability cutoff below which we omit a theta dihedral constraint.  Default 0.55, or "
		"whatever is set on the commandline with the -trRosetta::theta_constraint_prob_cutoff "
		"commandline option.",
		"0.55" )

		+ XMLSchemaAttribute::attribute_w_default( "phi_constraint_prob_cutoff", xsct_real,
		"Set the probability cutoff below which we omit a phi angle constraint.  Default 0.65, or "
		"whatever is set on the commandline with the -trRosetta::phi_constraint_prob_cutoff "
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
		"1.0" )

		+ XMLSchemaAttribute::attribute_w_default( "sequence", xs_string,
		"The amino acid sequence to predict.  EITHER this OR a FASTA file must be provided. "
		"Sequences must be single-letter amino acid codes, and must contain "
		"only the 20 canonical amino acids.",
		"" )

		+ XMLSchemaAttribute::attribute_w_default( "fasta_file", xs_string ,
		"A FASTA file containing a single sequence, the amino acid sequence to predict.  EITHER this OR "
		"a sequence must be provided.  Sequences must be single-letter amino acid codes, and must contain "
		"only the 20 canonical amino acids.  A FASTA file can also be set with the -in:file:fasta commandline "
		"flag, which sets the default for this mover (overrideable either with the fasta_file option or the "
		"sequence option).",
		"" )

		+ XMLSchemaAttribute::attribute_w_default( "backbone_randomization_mode", xs_string ,
		"The manner in wihch the polypeptide backbone will be initially randomized.  Options are 'classic' "
		"(the manner used in the original Yang et al. PyRosetta protocol, which randomly selects from one of "
		"six phi/psi pairs for each residue), 'ramachandran' (randomizing biased by the Ramachandran preferences "
		"of each amino acid type), or 'bins' (randomizing biased by the probabilities of residue type i being in "
		"backbone bin X and residue type i+1 being in backbone bin Y).  Defaults to 'classic', or whatever is set "
		"at the commandline with the -trRosetta::backbone_randomization_mode commandline option.", "classic" )

		+ XMLSchemaAttribute::attribute_w_default( "backbone_minimization_mode", xs_string ,
		"The manner in wihch the polypeptide backbone will be minimized using the constraints from the trRosetta "
		"neural network.  Options are: 'classic0' (minimize using short-range constraints, then minimize using "
		"medium-range constraints, then minimize using long-range constraints), 'classic1' (minimize using "
		"short- and medium-range constraints, then minimize using long-range constraints), or 'classic2' (minimize "
		"using all constraints).  Defaults to 'classic2', or whatever is set at the commandline with the "
		"-trRosetta::backbone_minimization_mode commandline option.", "classic2" )

		+ XMLSchemaAttribute::attribute_w_default( "cis_peptide_prob_non_prepro", xsct_real,
		"The probability of sampling a cis peptide bond at a position that is NOT followed by a proline when "
		"'ramachandran' backbone randomization mode is used.  Defaults to 0.0005 (or a setting provided at the "
		"commandline with the -trRosetta:cis_peptide_prob_non_prepro flag).  Ignored for 'classic' or "
		"'bins' modes.", "0.0005" )

		+ XMLSchemaAttribute::attribute_w_default( "cis_peptide_prob_prepro", xsct_real,
		"The probability of sampling a cis peptide bond at a position that IS followed by a proline when "
		"'ramachandran' backbone randomization mode is used.  Defaults to 0.05 (or a setting provided at the "
		"commandline with the -trRosetta:cis_peptide_prob_prepro flag).  Ignored for 'classic' or "
		"'bins' modes.", "0.05" )

		+ XMLSchemaAttribute::attribute_w_default( "scorefxn0", xs_string,
		"The scoring function used for stage 0 energy minimization.  Defaults to trRosetta_cen0 (or to "
		"whatever is set on the commandline with the -trRosetta:scorefxn0 commandline option).",
		"trRosetta_cen0" )

		+ XMLSchemaAttribute::attribute_w_default( "scorefxn1", xs_string,
		"The scoring function used for stage 1 energy minimization.  Defaults to trRosetta_cen1 (or to "
		"whatever is set on the commandline with the -trRosetta:scorefxn1 commandline option).",
		"trRosetta_cen1" )

		+ XMLSchemaAttribute::attribute_w_default( "scorefxn2", xs_string,
		"The scoring function used for stage 2 (Van der Waals) energy minimization.  Defaults to trRosetta_cen2 (or to "
		"whatever is set on the commandline with the -trRosetta:scorefxn2 commandline option).",
		"trRosetta_cen2" )

		+ XMLSchemaAttribute::attribute_w_default( "scorefxn3", xs_string,
		"The scoring function used for stage 3 (Cartesian) energy minimization.  Defaults to trRosetta_cart (or to "
		"whatever is set on the commandline with the -trRosetta:scorefxn3 commandline option).",
		"trRosetta_cart" )

		+ XMLSchemaAttribute::attribute_w_default("mutate_gly_to_ala", xsct_rosetta_bool,
		"If true, glycine residues are mutated to alanine during the initial centroid phases of "
		"minimization to match the original PyRosetta trRosetta protocol (then mutated back to "
		"glycine for fullatom refinement).  True by default.", "true")

		+ XMLSchemaAttribute::attribute_w_default("fullatom_refinement", xsct_rosetta_bool,
		"If true, we do fullatom refinement at the end with the FastRelax protocol, using the scoring function "
		"specified with the scorefxn_fullatom option.  If the atom_pair, "
		"dihedral, and angle constraint scoreterms are not on, they are turned on.  True by default.", "true")

		+ XMLSchemaAttribute::attribute_w_default("scorefxn_fullatom", xs_string,
		"Weights file for scorefunction used for fullatom refinement with FastRelax.  If atom-pair_constraint, "
		"dihedral_constriant, or angle_constraint terms are zero, they will be set to 5.0, 1.0, and 1.0 "
		"respectively.  If empty (the default), then the scoring function specified with -score:weights is used "
		"instead.", "");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Implements the full trRosetta protocol, as described in Yang et al. (2020) Improved protein "
		"structure prediction using predicted interresidue orientations. Proc. Natl. Acad. Sci. USA "
		"117(3):1496-503. https://doi.org/10.1073/pnas.1914677117.  This mover takes as input a multiple "
		"sequence alignment, runs the trRosetta neural network, generates distance and angle constraints "
		"between pairs of residues, and carries out energy-minimization to produce a structure.  Note that "
		"this mover deletes and replaces the input structure.  If a native structure is provided, the "
		"mover tags the output structure with the RMSD to native.\n\n"
		+ basic::tensorflow_manager::get_tensorflow_compilation_instructions( mover_name() + " mover", true ),
		attlist
	);
}

////////////////////////////////////////////////////////////////////////////////
// PUBLIC FUNCTIONS -- SETTERS
////////////////////////////////////////////////////////////////////////////////

#ifdef USE_TENSORFLOW

/// @brief Set the amino acid sequence of the protein being predicted.
/// @details Must be one-letter codes; no NCAAs allowed.  Overrides any setting for
/// fasta_file_ (i.e. sets fasta_file_ to "").
void
trRosettaProtocolMover::set_sequence(
	std::string const & seq_in
) {
	sequence_ = seq_in;
	fasta_file_ = "";
	if ( !sequence_.empty() ) {
		validate_sequence();
	}
}

/// @brief Checks that the current sequence is valid, and throws if it is not.
void
trRosettaProtocolMover::validate_sequence() const {
	std::string const errmsg( "Error in trRosettaProtocolMover::validate_sequence(): " );
	runtime_assert_string_msg( !(sequence_.empty()), errmsg + "The sequence provided to the trRosettaProtocol was empty!" );
	for ( core::Size i(0), imax(sequence_.length()); i<imax; ++i ) {
		runtime_assert_string_msg(
			core::chemical::is_canonical_L_aa_or_gly( core::chemical::aa_from_oneletter_code( sequence_[i] ) ),
			errmsg + "Amino acid \"" + sequence_[i] + "\" at position " + std::to_string(i+1) + " is not recognized as a canonical amino acid!"
		);
	}
}

/// @brief Set the FASTA file containing the sequence of the protein to predict.
/// @details Overrides any setting for sequence_ (i.e. sets sequence_ to "").
void
trRosettaProtocolMover::set_fasta_file(
	std::string const & filename
) {
	fasta_file_ = filename;
	sequence_ = "";
}


/// @brief Set the multiple sequence alignment filename.
/// @details Resets the ConstraintGenerator, if already loaded.
void
trRosettaProtocolMover::set_msa_file(
	std::string const & filename
) {
	msa_file_ = filename;
	constraint_generator_ = nullptr;
}

/// @brief Set whether we are using the constraint generator to set
/// distance constraints.
void
trRosettaProtocolMover::set_use_distance_constraints(
	bool const setting
) {
	use_distance_constraints_ = setting;
}

/// @brief Set whether we are using the constraint generator to set
/// omega dihedral constraints.
/// @details This is NOT the backbone omega dihedral.  It is the dihedral angle
/// between CA1-CB1-CB2-CA2.
void
trRosettaProtocolMover::set_use_omega_constraints(
	bool const setting
) {
	use_omega_constraints_ = setting;
}

/// @brief Set whether we are using the constraint generator to set
/// theta dihedral constraints.
/// @details This is the dihedral between N1-CA1-CB1-CB2.
void
trRosettaProtocolMover::set_use_theta_constraints(
	bool const setting
) {
	use_theta_constraints_ = setting;
}

/// @brief Set whether we are using the constraint generator to set
/// phi angle constraints.
/// @details This is NOT the backbone phi dihedral.  It is the angle between
/// CA1-CB1-CB2.
void
trRosettaProtocolMover::set_use_phi_constraints(
	bool const setting
) {
	use_phi_constraints_ = setting;
}

/// @brief Set the backbone randomization mode.
void
trRosettaProtocolMover::set_backbone_randomization_mode(
	trRosettaProtocolBackboneRandomizationMode const mode_in
) {
	debug_assert( static_cast<core::Size>(mode_in) > 0 && static_cast<core::Size>(mode_in) <= static_cast<core::Size>(trRosettaProtocolBackboneRandomizationMode::NUM_ENTRIES) );
	backbone_randomization_mode_ = mode_in;
}

/// @brief Set the backbone randomization mode, using the mode name string.
void
trRosettaProtocolMover::set_backbone_randomization_mode(
	std::string const & mode_in
) {
	set_backbone_randomization_mode( get_randomization_mode_enum_from_string(mode_in) );
}

/// @brief Set the backbone minimization mode.
void
trRosettaProtocolMover::set_backbone_minimization_mode(
	trRosettaProtocolBackboneMinimizationMode const mode_in
) {
	debug_assert( static_cast<core::Size>(mode_in) > 0 && static_cast<core::Size>(mode_in) <= static_cast<core::Size>(trRosettaProtocolBackboneMinimizationMode::NUM_ENTRIES) );
	backbone_minimization_mode_ = mode_in;
}

/// @brief Set the backbone minimization mode, using the mode name string.
void
trRosettaProtocolMover::set_backbone_minimization_mode(
	std::string const & mode_in
) {
	set_backbone_minimization_mode( get_minimization_mode_enum_from_string( mode_in ) );
}


/// @brief Set the probability of sampling a cis peptide bond in 'ramachandran' mode,
/// at non-pre-proline and pre-proline positions.
void
trRosettaProtocolMover::set_ramachandran_mode_cis_probabilities(
	core::Real const non_prepro_prob,
	core::Real const prepro_prob
) {
	std::string const errmsg( "Error in trRosettaProtocolMover::set_ramachandran_mode_cis_probabilities(): " );
	runtime_assert_string_msg( non_prepro_prob >= 0.0 && non_prepro_prob <= 1.0, errmsg + "The probability of sampling a cis peptide bond at a position that is not followed by proline must be between 0.0 and 1.0.");
	runtime_assert_string_msg( prepro_prob >= 0.0 && prepro_prob <= 1.0, errmsg + "The probability of sampling a cis peptide bond at a position that is followed by proline must be between 0.0 and 1.0.");
	ramachandran_mode_cis_probability_non_prepro_ = non_prepro_prob;
	ramachandran_mode_cis_probability_prepro_ = prepro_prob;
}

/// @brief Set the weights file to use for minimization stage 0, 1, 2, or 3.
/// @details Triggers read of weights from disk!
void
trRosettaProtocolMover::set_scorefunction_for_minimization_stage(
	std::string const & weights_filename,
	core::Size const stage
) {
	switch( stage ) {
	case 0 :
		sfxn0_ = core::scoring::ScoreFunctionFactory::create_score_function( weights_filename );
		set_emethod_options(*sfxn0_);
		break;
	case 1 :
		sfxn1_ = core::scoring::ScoreFunctionFactory::create_score_function( weights_filename );
		set_emethod_options(*sfxn1_);
		break;
	case 2 :
		sfxn2_ = core::scoring::ScoreFunctionFactory::create_score_function( weights_filename );
		set_emethod_options(*sfxn2_);
		break;
	case 3 :
		sfxn3_ = core::scoring::ScoreFunctionFactory::create_score_function( weights_filename );
		set_emethod_options(*sfxn3_);
		break;
	default :
		utility_exit_with_message( "Error in trRosettaProtocolMover::set_scorefunction_for_minimization_stage(): Stage must be 0, 1, 2, or 3!" );
	};
	TR << "Set scorefunction for centroid stage " << stage  << " to " << weights_filename << "." << std::endl;
}

/// @brief Set the scoring function used for fullatom refinement.
/// @details Triggers read of weights from disk!
/// @note Turns on atom_pair_constraint, dihedral_constraint, and angle_constraint terms if they
/// are not on in the scorefunction that is loaded.  Weights filename is deliberately passed by
/// copy.
void
trRosettaProtocolMover::set_scorefunction_for_fullatom_refinement(
	std::string weights_filename
) {
	if ( weights_filename.empty() ) {
		weights_filename = basic::options::option[basic::options::OptionKeys::score::weights](); //Get default if this is empty.
	}
	sfxn_fullatom_ = core::scoring::ScoreFunctionFactory::create_score_function( weights_filename );
	TR << "Set scorefunction for full-atom refinement to " << weights_filename << "." << std::endl;
	if ( sfxn_fullatom_->has_zero_weight( core::scoring::atom_pair_constraint ) ) {
		sfxn_fullatom_->set_weight( core::scoring::atom_pair_constraint, 5.0 );
		TR << "Setting full-atom refinement scorefunction atom_pair_constraint weight to 5.0 (was 0.0)." << std::endl;
	}
	if ( sfxn_fullatom_->has_zero_weight( core::scoring::dihedral_constraint ) ) {
		sfxn_fullatom_->set_weight( core::scoring::dihedral_constraint, 1.0 );
		TR << "Setting full-atom refinement scorefunction dihedral_constraint weight to 1.0 (was 0.0)." << std::endl;
	}
	if ( sfxn_fullatom_->has_zero_weight( core::scoring::angle_constraint ) ) {
		sfxn_fullatom_->set_weight( core::scoring::angle_constraint, 1.0 );
		TR << "Setting full-atom refinement scorefunction angle_constraint weight to 1.0 (was 0.0)." << std::endl;
	}
}

/// @brief Set the probability cutoffs for distance, omega, theta, and phi.
void
trRosettaProtocolMover::set_prob_cutoffs(
	core::Real const dist_cutoff,
	core::Real const omega_cutoff,
	core::Real const theta_cutoff,
	core::Real const phi_cutoff
) {
	std::string const errmsg("Error in trRosettaProtocolMover::set_prob_cutoffs(): ");
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
trRosettaProtocolMover::set_constraint_weights(
	core::Real const dist_weight,
	core::Real const omega_weight,
	core::Real const theta_weight,
	core::Real const phi_weight
) {
	std::string const errmsg( "Error in trRosettaProtocolMover::set_constraint_weights(): " );
	runtime_assert_string_msg( dist_weight >= 0.0, "Distance constraint weight must be positive!" );
	runtime_assert_string_msg( omega_weight >= 0.0, "Omega dihedral constraint weight must be positive!" );
	runtime_assert_string_msg( theta_weight >= 0.0, "Theta dihedral constraint weight must be positive!" );
	runtime_assert_string_msg( phi_weight >= 0.0, "Phi angle constraint weight must be positive!" );
	distance_constraint_weight_= dist_weight;
	omega_constraint_weight_ = omega_weight;
	theta_constraint_weight_ = theta_weight;
	phi_constraint_weight_ = phi_weight;
}

/// @brief Sets whether glycine residues should be mutated to alanine during the centroid phase.
/// @details True by default to match the original PyRosetta protocol.
void
trRosettaProtocolMover::set_mutate_gly_to_ala(
	bool const setting
) {
	mutate_gly_to_ala_ = setting;
}

/// @brief Set whether we do fullatom refinement (with FastRelax) at the end.  Default true.
void
trRosettaProtocolMover::set_fullatom_refinement(
	bool const setting
) {
	fullatom_refinement_ = setting;
}

/// @brief Get the current backbone randomization mode's name, as a string.
std::string
trRosettaProtocolMover::get_backbone_randomization_mode_name() const {
	return get_randomization_mode_string_from_enum( backbone_randomization_mode_ );
}

/// @brief Get the current backbone minimization mode's name, as a string.
std::string
trRosettaProtocolMover::get_backbone_minimization_mode_name() const {
	return get_minimization_mode_string_from_enum( backbone_minimization_mode_ );
}

#endif //USE_TENSORFLOW

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
trRosettaProtocolMover::fresh_instance() const
{
	return utility::pointer::make_shared< trRosettaProtocolMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
trRosettaProtocolMover::clone() const
{
	return utility::pointer::make_shared< trRosettaProtocolMover >( *this );
}

std::string trRosettaProtocolMover::get_name() const {
	return mover_name();
}

std::string trRosettaProtocolMover::mover_name() {
	return "trRosettaProtocol";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
trRosettaProtocolMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< trRosettaProtocolMover >();
}

std::string
trRosettaProtocolMoverCreator::keyname() const
{
	return trRosettaProtocolMover::mover_name();
}

void trRosettaProtocolMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	trRosettaProtocolMover::provide_xml_schema( xsd );
}

/// @brief This mover is unpublished.  It returns Vikram K. Mulligan as its author.
/// It also adds Yang et al. (2020) as a citiation for trRosetta itself.
void
trRosettaProtocolMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"trRosettaProtocol", basic::citation_manager::CitedModuleType::Mover,
		"Vikram K. Mulligan",
		"Systems Biology, Center for Computational Biology, Flatiron Institute",
		"vmulligan@flatironinstitute.org",
		"Converted the Python trRosetta protocol from Yang et al. (2020) to C++, and implemented it as the trRosettaProtocolMover."
		)
	);

	citations.add( protocols::trRosetta::trRosettaProtocol_v1::get_trRosetta_neural_net_citation() );

#ifdef USE_TENSORFLOW
	if ( fullatom_refinement() )
#endif
	{
		protocols::relax::FastRelax frlx;
		frlx.provide_citation_info(citations);

		core::simple_metrics::metrics::RMSDMetric rmsdmetric;
		rmsdmetric.provide_citation_info(citations);
	}

	core::simple_metrics::metrics::TotalEnergyMetric totalE;
	totalE.provide_citation_info(citations);

	core::simple_metrics::metrics::TimingProfileMetric timing;
	timing.provide_citation_info(citations);
}


////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

#ifdef USE_TENSORFLOW

/// @brief Initialize this mover from an options collection.
/// @details Intended to be called once and once only, from the
/// constructor for this mover.
/// @note Read from disk of MSA file is deferred to apply time.
void
trRosettaProtocolMover::init_from_options(
	utility::options::OptionCollection const & opts
) {
	using namespace basic::options::OptionKeys;

	std::string const errmsg( "Error in trRosettaProtocolMover::init_from_options(): " );

	// Import the native pose:
	if ( opts[in::file::native].user() ) {
		TR << "Setting native structure in trRosettaProtocol mover." << std::endl;
		if ( protocols::jd2::jd2_used() ) {
			protocols::jd2::set_native_in_mover( *this, opts );
		} else {
			core::pose::PoseOP native_pose( utility::pointer::make_shared< core::pose::Pose >() );
			core::chemical::ResidueTypeSetCOP rsd_set;
			if ( opts[ in::file::fullatom ]() ) {
				rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
			} else {
				rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
			}
			core::import_pose::pose_from_file( *native_pose, *rsd_set, opts[ in::file::native ]() , core::import_pose::PDB_file);
			set_native_pose(native_pose);
		}
		runtime_assert_string_msg( get_native_pose() != nullptr, errmsg + "Could not load native pose!" );
	} else {
		TR << "No native structure specified.  Not setting native structure in trRosettaProtocol mover." << std::endl;
	}

	// Set the MSA file.  May be overridden later.
	set_msa_file( opts[basic::options::OptionKeys::trRosetta::msa_file]() );

	// Set the FASTA file. May be overriden later.
	if ( opts[basic::options::OptionKeys::in::file::fasta].user() ) {
		runtime_assert_string_msg( opts[basic::options::OptionKeys::in::file::fasta].size() == 1, errmsg + "Please pass exactly one FASTA file with the -in:file:fasta commandline option." );
		set_fasta_file( opts[basic::options::OptionKeys::in::file::fasta]()[1] );
	}

	// Set options for the constraints to be used.  May be overridden later.
	set_use_distance_constraints( opts[basic::options::OptionKeys::trRosetta::use_distance_constraints]() );
	set_use_omega_constraints( opts[basic::options::OptionKeys::trRosetta::use_omega_constraints]() );
	set_use_theta_constraints( opts[basic::options::OptionKeys::trRosetta::use_theta_constraints]() );
	set_use_phi_constraints( opts[basic::options::OptionKeys::trRosetta::use_phi_constraints]() );

	// Set backbone randomization mode.  May be overridden later.
	set_backbone_randomization_mode( opts[basic::options::OptionKeys::trRosetta::backbone_randomization_mode]() );

	// Set probabilities for sampling cis peptide bonds.  May be overridden later.
	set_ramachandran_mode_cis_probabilities(
		opts[basic::options::OptionKeys::trRosetta::cis_peptide_prob_non_prepro](),
		opts[basic::options::OptionKeys::trRosetta::cis_peptide_prob_prepro]()
	);

	// Set scorefunctions.  May be overrridden later.
	set_scorefunction_for_minimization_stage( opts[basic::options::OptionKeys::trRosetta::scorefxn0](), 0 );
	set_scorefunction_for_minimization_stage( opts[basic::options::OptionKeys::trRosetta::scorefxn1](), 1 );
	set_scorefunction_for_minimization_stage( opts[basic::options::OptionKeys::trRosetta::scorefxn2](), 2 );
	set_scorefunction_for_minimization_stage( opts[basic::options::OptionKeys::trRosetta::scorefxn3](), 3 );
	set_scorefunction_for_fullatom_refinement( opts[basic::options::OptionKeys::trRosetta::scorefxn_fullatom]() );

	// Set backbone minimization mode.  May be overridden later.
	set_backbone_minimization_mode( opts[basic::options::OptionKeys::trRosetta::backbone_minimization_mode]() );

	// Set probability cutoffs.  May be overridden later.
	set_prob_cutoffs(
		opts[basic::options::OptionKeys::trRosetta::distance_constraint_prob_cutoff](),
		opts[basic::options::OptionKeys::trRosetta::omega_constraint_prob_cutoff](),
		opts[basic::options::OptionKeys::trRosetta::theta_constraint_prob_cutoff](),
		opts[basic::options::OptionKeys::trRosetta::phi_constraint_prob_cutoff]()
	);

	// Set the constraint weights.
	set_constraint_weights(
		opts[ basic::options::OptionKeys::trRosetta::distance_constraint_weight ](),
		opts[ basic::options::OptionKeys::trRosetta::omega_constraint_weight ](),
		opts[ basic::options::OptionKeys::trRosetta::theta_constraint_weight ](),
		opts[ basic::options::OptionKeys::trRosetta::phi_constraint_weight ]()
	);

	// Set whether we mutate glycines to alanine during centroid phase.
	set_mutate_gly_to_ala( opts[basic::options::OptionKeys::trRosetta::mutate_gly_to_ala]() );
}

/// @brief Create the trRosettaConstraintGenerator, apply it to the pose, and return a vector
/// of constraint objects.
/// @details If probability_offset is provided, the probability thresholds are all shifted by
/// this amount.  Used during fullatom refinement.
utility::vector1< core::scoring::constraints::ConstraintCOP >
trRosettaProtocolMover::generate_trRosetta_constraints(
	core::pose::Pose const & pose,
	core::Real const probability_offset /*=0.0*/
) const {
	using namespace protocols::trRosetta_protocols::constraint_generators;

	//Create the constraint generator, if needed:
	if ( constraint_generator_ == nullptr ) {
		// Initialize the constraint generator:
		TR << "Initializing constraint generator from multiple sequence alignment file \"" << msa_file_ << "\"." << std::endl;
		constraint_generator_ = utility::pointer::make_shared< trRosettaConstraintGenerator >(msa_file_);
	} else {
		if ( constraint_generator_->msa_file() != msa_file_ ) {
			TR.Warning << "Constraint generator was configured, but its MSA file was set to \"" << constraint_generator_->msa_file() << "\", while the trRosettaProtocol mover was set to use msa file \"" << msa_file_ << "\".  Instantiating new trRosetta constraint generator to match mover configuration." << std::endl;
			constraint_generator_ = utility::pointer::make_shared< trRosettaConstraintGenerator >(msa_file_);
		}
	}

	//Configure the constraint generator:
	constraint_generator_->set_generate_dist_constraints( use_distance_constraints() );
	constraint_generator_->set_generate_omega_constraints( use_omega_constraints() );
	constraint_generator_->set_generate_theta_constraints( use_theta_constraints() );
	constraint_generator_->set_generate_phi_constraints( use_phi_constraints() );
	constraint_generator_->set_prob_cutoffs(
		dist_prob_cutoff() + probability_offset,
		omega_prob_cutoff() + probability_offset,
		theta_prob_cutoff() + probability_offset,
		phi_prob_cutoff() + probability_offset
	);
	constraint_generator_->set_constraint_weights(
		distance_constraint_weight(),
		omega_constraint_weight(),
		theta_constraint_weight(),
		phi_constraint_weight()
	);

	TR << "Generating trRosetta constraints." << std::endl;

	return constraint_generator_->apply(pose);
}

/// @brief Remove previously-added trRosetta constraints from the pose.
void
trRosettaProtocolMover::remove_trRosetta_constraints(
	utility::vector1< core::scoring::constraints::ConstraintCOP > const & constraints,
	core::pose::Pose & pose
) const {
	runtime_assert_string_msg( pose.remove_constraints( constraints, true ), "Program error in trRosettaProtocolMover::remove_trRosetta_constraints(): Something went wrong when trying to remove trRosetta constraints." );
}

/// @brief Add constraints from a list to the pose, if the constraints are between residues that are
/// separated by at least min_seqsep but less than max_seqsep residues.
/// @details This does not clear constraints from the pose.
/// @note It is "less than" max_seqsep, not "less than or equal to".
void
trRosettaProtocolMover::add_constraints_to_pose(
	core::pose::Pose & pose,
	utility::vector1< core::scoring::constraints::ConstraintCOP > const & trRosetta_constraints,
	core::Size const min_seqsep,
	core::Size const max_seqsep,
	bool const skip_glycine_positions /*=false*/
) const {
	for ( auto const & cst : trRosetta_constraints ) {
		core::Size res1, res2;
		get_residues_from_constraint( res1, res2, cst );
		if ( skip_glycine_positions ) {
			if (
					pose.residue_type(res1).aa() == core::chemical::aa_gly ||
					pose.residue_type(res2).aa() == core::chemical::aa_gly
					) {
				continue;
			}
		}
		debug_assert(res2 > res1);
		core::Size const seqdist( res2 - res1 ); //Assumes res2>res1
		if ( seqdist >= min_seqsep && seqdist < max_seqsep ) {
			pose.add_constraint( cst );
		}
	}
}

/// @brief Given a constraint, determine if it is an AtomPairConstraint, an
/// AngleConstraint, or a DihedralConstraint, pull out the pair of residues
/// that are constrained, and return the pair.
/// @details Throws if type is unrecognized or if more than two residues are
/// constrained.  Values of res1 and res2 are overwritten by this operation.
void
trRosettaProtocolMover::get_residues_from_constraint(
	core::Size & res1,
	core::Size & res2,
	core::scoring::constraints::ConstraintCOP const & cst
) const {
	if ( get_residues_from_atom_pair_constraint( res1, res2, cst ) ) {
		return;
	}
	if ( get_residues_from_angle_constraint( res1, res2, cst ) ) {
		return;
	}
	if ( get_residues_from_dihedral_constraint( res1, res2, cst ) ) {
		return;
	}
	utility_exit_with_message( "trRosettaProtocolMover::get_residues_from_constraint(): Selected constraint was not an AtomPairConstraint, an AngleConstraint, or a DihedralConstraint!" );
}

/// @brief Given a constraint, determine if it is an AtomPairConstraint pull
/// out the pair of residues that are constrained, and return the pair.
/// @details If successful, values of res1 and res2 are overwritten by this
/// operation, and the function returns "true".  Otherwise, values are not
/// altered, and the function returns "false".
bool
trRosettaProtocolMover::get_residues_from_atom_pair_constraint(
	core::Size & res1,
	core::Size & res2,
	core::scoring::constraints::ConstraintCOP const & cst
) const {
	using namespace core::scoring::constraints;
	AtomPairConstraintCOP atpair_cst( utility::pointer::dynamic_pointer_cast< AtomPairConstraint const >(cst) );
	if ( atpair_cst == nullptr ) return false;
	res1 = atpair_cst->atom1().rsd();
	res2 = atpair_cst->atom2().rsd();
	if ( res1 > res2 ) {
		std::swap(res1, res2);
	}
	return true;
}

/// @brief Given a constraint, determine if it is an AngleConstraint, pull
/// out the pair of residues that are constrained, and return the pair.
/// @details If successful, values of res1 and res2 are overwritten by this
/// operation, and the function returns "true".  Otherwise, values are not
/// altered, and the function returns "false".
/// @note Throws if the middle residue doesn't match either the first or
/// last.
bool
trRosettaProtocolMover::get_residues_from_angle_constraint(
	core::Size & res1,
	core::Size & res2,
	core::scoring::constraints::ConstraintCOP const & cst
) const {
	using namespace core::scoring::constraints;
	AngleConstraintCOP angle_cst( utility::pointer::dynamic_pointer_cast< AngleConstraint const >(cst) );
	if ( angle_cst == nullptr ) return false;
	res1 = angle_cst->atom(1).rsd();
	res2 = angle_cst->atom(3).rsd();
	core::Size const midres( angle_cst->atom(2).rsd() );
	runtime_assert_string_msg(
		midres == res1 || midres == res2,
		"Error in trRosettaProtocolMover::get_residues_from_angle_constraint(): This AngleConstraint constrains three different residues!"
	);
	if ( res1 > res2 ) {
		std::swap(res1, res2);
	}
	return true;
}

/// @brief Given a constraint, determine if it is a DihedralConstraint, pull
/// out the pair of residues that are constrained, and return the pair.
/// @details If successful, values of res1 and res2 are overwritten by this
/// operation, and the function returns "true".  Otherwise, values are not
/// altered, and the function returns "false".
/// @note Throws if the middle residues doesn't match either the first or
/// last.
bool
trRosettaProtocolMover::get_residues_from_dihedral_constraint(
	core::Size & res1,
	core::Size & res2,
	core::scoring::constraints::ConstraintCOP const & cst
) const {
	using namespace core::scoring::constraints;
	DihedralConstraintCOP dihedral_cst( utility::pointer::dynamic_pointer_cast< DihedralConstraint const >(cst) );
	if ( dihedral_cst == nullptr ) return false;
	res1 = dihedral_cst->atom(1).rsd();
	res2 = dihedral_cst->atom(4).rsd();
	core::Size const midres1( dihedral_cst->atom(2).rsd() );
	core::Size const midres2( dihedral_cst->atom(3).rsd() );
	runtime_assert_string_msg(
		(midres1 == res1 || midres1 == res2) && (midres2 == res1 || midres2 == res2),
		"Error in trRosettaProtocolMover::get_residues_from_dihedral_constraint(): This DihedralConstraint constrains more than two different residues!"
	);
	if ( res1 > res2 ) {
		std::swap(res1, res2);
	}
	return true;
}

/// @brief Perform minimization using short-range constraints, then medium-range, then long-range.
/// @details Note that constraints are cumulative (short-range only in round 1, short- and medium-range
/// in round 2, and short-, medium-, and long-range in round 3).
void
trRosettaProtocolMover::perform_classic0_minimization_protocol(
	core::pose::Pose & pose,
	utility::vector1< core::scoring::constraints::ConstraintCOP > const & trRosetta_constraints,
	protocols::moves::MoverOP repeat_minmover0,
	protocols::moves::MoverOP minmover1,
	protocols::moves::MoverOP minmover3
) const {
	TR << "Performing classic0 minimization protocol" << std::endl;
	TR << "Adding short-range sequence constraints (residues less than " << medium_range_seqsep_cutoff_ << " residues apart in linear sequence)." << std::endl;
	add_constraints_to_pose( pose, trRosetta_constraints, 1, medium_range_seqsep_cutoff_ );
	perform_classic_minsteps( pose, repeat_minmover0, minmover1, minmover3 );

	TR << "Adding medium-range sequence constraints (residues less than " << long_range_seqsep_cutoff_ << " residues apart in linear sequence)." << std::endl;
	add_constraints_to_pose( pose, trRosetta_constraints, medium_range_seqsep_cutoff_, long_range_seqsep_cutoff_ );
	perform_classic_minsteps( pose, repeat_minmover0, minmover1, minmover3 );

	TR << "Adding long-range sequence constraints (all residues)." << std::endl;
	add_constraints_to_pose( pose, trRosetta_constraints, long_range_seqsep_cutoff_, pose.total_residue() );
	perform_classic_minsteps( pose, repeat_minmover0, minmover1, minmover3 );
}

/// @brief Perform minimization using short-range and medium-range constraints, then long-range.
/// @details Note that constraints are cumulative (short-range and medium-range
/// in round 1, and short-, medium-, and long-range in round 2).
/// @note For some reason, this protocol does not constrain residues closer than 3 amino acids apart.
void
trRosettaProtocolMover::perform_classic1_minimization_protocol(
	core::pose::Pose & pose,
	utility::vector1< core::scoring::constraints::ConstraintCOP > const & trRosetta_constraints,
	protocols::moves::MoverOP repeat_minmover0,
	protocols::moves::MoverOP minmover1,
	protocols::moves::MoverOP minmover3
) const {
	TR << "Performing classic1 minimization protocol" << std::endl;
	TR << "Adding short- and medium-range sequence constraints (residues less than " << long_range_seqsep_cutoff_ << " residues apart in linear sequence)." << std::endl;
	add_constraints_to_pose( pose, trRosetta_constraints, 3, long_range_seqsep_cutoff_ );
	perform_classic_minsteps( pose, repeat_minmover0, minmover1, minmover3 );

	TR << "Adding long-range sequence constraints (all residues)." << std::endl;
	add_constraints_to_pose( pose, trRosetta_constraints, long_range_seqsep_cutoff_, pose.total_residue() );
	perform_classic_minsteps( pose, repeat_minmover0, minmover1, minmover3 );
}


/// @brief Perform minimization using short-, medium-, and long-range constraints all in one go.
void
trRosettaProtocolMover::perform_classic2_minimization_protocol(
	core::pose::Pose & pose,
	utility::vector1< core::scoring::constraints::ConstraintCOP > const & trRosetta_constraints,
	protocols::moves::MoverOP repeat_minmover0,
	protocols::moves::MoverOP minmover1,
	protocols::moves::MoverOP minmover3
) const {
	TR << "Performing classic2 minimization protocol" << std::endl;
	TR << "Adding long-range sequence constraints (all residues)." << std::endl;
	add_constraints_to_pose( pose, trRosetta_constraints, 1, pose.total_residue() );
	perform_classic_minsteps( pose, repeat_minmover0, minmover1, minmover3 );
}

/// @brief Perform the minimization steps common to the above three protocols.
void
trRosettaProtocolMover::perform_classic_minsteps(
	core::pose::Pose & pose,
	protocols::moves::MoverOP repeat_minmover0,
	protocols::moves::MoverOP minmover1,
	protocols::moves::MoverOP minmover3
) const {
	TR << "Performing torsion-space minimization with scorefunction 0." << std::endl;
	repeat_minmover0->apply(pose);
	TR << "Performing Cartesian minimization with scorefunction 3." << std::endl;
	minmover3->apply(pose);
	TR << "Performing torsion-space clash removal minimization with scorefunction 1 (validating with scorefunction 2)." << std::endl;
	remove_clashes( pose, *minmover1, *sfxn2_ );
}

/// @brief Set up a MinMover for use by this protocol.
protocols::minimization_packing::MinMoverOP
trRosettaProtocolMover::configure_minmover(
	core::kinematics::MoveMapOP movemap,
	core::scoring::ScoreFunctionCOP sfxn,
	core::Size const maxiters,
	bool const do_cartesian
) const {
	using namespace protocols::minimization_packing;
	MinMoverOP minmover( utility::pointer::make_shared< MinMover >() );
	minmover->score_function( sfxn );
	minmover->set_movemap(movemap);
	minmover->tolerance(0.0001);
	minmover->nb_list(true);
	minmover->max_iter(maxiters);
	minmover->cartesian(do_cartesian);
	return minmover;
}

/// @brief Convert a pose to fullatom from centroid.
void
trRosettaProtocolMover::convert_to_fullatom(
	core::pose::Pose & pose
) const {
	core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
}

/// @brief Do all-atom FastRelax refinement.
void
trRosettaProtocolMover::do_fullatom_refinement(
	core::pose::Pose & pose
) const {
	debug_assert( sfxn_fullatom_ != nullptr );
	TR << "Performing all-atom refinement using the FastRelax protocol with DualSpace." << std::endl;

	core::kinematics::MoveMapOP movemap( utility::pointer::make_shared< core::kinematics::MoveMap >() );
	movemap->set_bb(true);
	movemap->set_chi(true);
	movemap->set_jump(true);

	protocols::relax::FastRelaxOP frlx( utility::pointer::make_shared< protocols::relax::FastRelax >() );
	frlx->set_scorefxn(sfxn_fullatom_);
	frlx->dualspace(true);
	frlx->set_movemap(movemap);

	frlx->apply(pose);
}

/// @brief Store the constraints score in the pose.
void
trRosettaProtocolMover::store_constraints_score(
	core::pose::Pose & pose,
	std::string const & metric_suffix
) const {
	core::scoring::ScoreFunctionOP sfxn( utility::pointer::make_shared< core::scoring::ScoreFunction >() );
	sfxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	sfxn->set_weight( core::scoring::angle_constraint, 1.0 );
	sfxn->set_weight( core::scoring::dihedral_constraint, 1.0 );

	core::simple_metrics::metrics::TotalEnergyMetric totalE;
	totalE.set_scorefunction(sfxn);
	totalE.apply( "Total_constraints_score_" + metric_suffix, pose );

	core::simple_metrics::metrics::TotalEnergyMetric atpairE;
	atpairE.set_scorefunction(sfxn);
	atpairE.set_scoretype( core::scoring::atom_pair_constraint );
	atpairE.apply( "Atom_pair_constraints_score_" + metric_suffix, pose );

	core::simple_metrics::metrics::TotalEnergyMetric angE;
	angE.set_scorefunction(sfxn);
	angE.set_scoretype( core::scoring::angle_constraint );
	angE.apply( "Angle_constraints_score_" + metric_suffix, pose );

	core::simple_metrics::metrics::TotalEnergyMetric dihedE;
	dihedE.set_scorefunction(sfxn);
	dihedE.set_scoretype( core::scoring::dihedral_constraint );
	dihedE.apply( "Dihedral_constraints_score_" + metric_suffix, pose );
}

/// @brief Store the RMSD to native in the pose.
/// @details Does nothing if native_pose == nullptr.
void
trRosettaProtocolMover::store_rmsd(
	core::pose::Pose & pose,
	core::pose::PoseCOP const & native_pose,
	std::string const & metric_suffix
) const {
	if ( native_pose != nullptr  ) {
		TR << "Computing RMSD to native..." << std::endl;
		core::simple_metrics::metrics::RMSDMetric rmsd_metric( native_pose );
		rmsd_metric.set_run_superimpose(true);
		rmsd_metric.set_rmsd_type( core::scoring::rmsd_protein_bb_ca );
		rmsd_metric.apply( "Backbone_CA_RMSD_to_native_" + metric_suffix, pose );

		core::simple_metrics::metrics::RMSDMetric rmsd_metric2( native_pose );
		rmsd_metric2.set_run_superimpose(true);
		rmsd_metric2.set_rmsd_type( core::scoring::rmsd_protein_bb_heavy_including_O );
		rmsd_metric2.apply( "Backbone_heavyatom_RMSD_to_native_" + metric_suffix, pose );
	} else {
		TR << "No native structure supplied.  Skipping RMSD to native." << std::endl;
	}
}


/// @brief Make a centroid pose from the sequence stored in sequence_.
/// @details Does not check that the sequence is valid first!
core::pose::Pose
trRosettaProtocolMover::make_centroid_pose_from_sequence() const {
	core::pose::Pose pose;
	core::chemical::ResidueTypeSetCOP restypeset(
		core::chemical::ChemicalManager::get_instance()->residue_type_set(
		core::chemical::CENTROID_t
		)
	);
	core::pose::make_pose_from_sequence( pose, sequence_, restypeset, true );
	return pose; //Return value optimization should prevent an expensive pose copying step here.
}

/// @brief Set or randomize the backbone phi, psi, and omega angles based on the backbone_randomization_mode_ setting.
/// @details Calls one of randomize_backbone_dihedrals_classic(), randomize_backbone_dihedrals_rama_prepro(), or
/// randomize_backbone_dihedrals_by_bins().
void
trRosettaProtocolMover::randomize_backbone_dihedrals(
	core::pose::Pose & pose
) const {
	TR << "Randomizing backbone dihedrals using protocol \"" << get_backbone_randomization_mode_name() << "\"." << std::endl;
	switch( backbone_randomization_mode_ ) {
	case trRosettaProtocolBackboneRandomizationMode::classic :
		randomize_backbone_dihedrals_classic(pose);
		break;
	case trRosettaProtocolBackboneRandomizationMode::ramachandran :
		randomize_backbone_dihedrals_rama_prepro(pose);
		break;
	case trRosettaProtocolBackboneRandomizationMode::bins :
		randomize_backbone_dihedrals_by_bins(pose);
		break;
	default :
		utility_exit_with_message("trRosettaProtocolMover::randomize_backbone_dihedrals(): Unrecognized backbone randomization mode!");
	};
}

/// @brief Set or randomize the backbone phi, psi, and omega angles using the method described in Yang et al.
/// @details This sets each residue's dihedrals to one of the following phi/psi combinations, chosen randomly:
/// -140  153 180 0.135 B
///  -72  145 180 0.155 B
/// -122  117 180 0.073 B
///  -82  -14 180 0.122 A
///  -61  -41 180 0.497 A
///   57   39 180 0.018 L
void
trRosettaProtocolMover::randomize_backbone_dihedrals_classic(
	core::pose::Pose & pose
) const {
	for ( core::Size ir(1), irmax(pose.total_residue()); ir<=irmax; ++ir ) {
		pose.set_omega( ir, 180.0 ); //Always set omega to trans.
		core::Real const randval( numeric::random::uniform() );
		if ( randval <= 0.135 ) {
			pose.set_phi( ir, -140.0 );
			pose.set_psi( ir,  153.0 );
		} else if ( randval <= 0.290 ) {
			pose.set_phi( ir,  -72.0 );
			pose.set_psi( ir,  145.0 );
		} else if ( randval <= 0.363 ) {
			pose.set_phi( ir, -122.0 );
			pose.set_psi( ir,  117.0 );
		} else if ( randval <= 0.485 ) {
			pose.set_phi( ir,  -82.0 );
			pose.set_psi( ir,  -14.0 );
		} else if ( randval <= 0.982 ) {
			pose.set_phi( ir,  -61.0 );
			pose.set_psi( ir,  -41.0 );
		} else {
			pose.set_phi( ir,   57.0 );
			pose.set_psi( ir,   39.0 );
		}
	}
	pose.update_residue_neighbors();
}

/// @brief Set or randomize the backbone phi, psi, and omega angles based on the Ramachandran preferences of
/// each amino acid type (using the rama_prepro scoreterm).
/// @details Omega is set to cis with 99.95% probability except at pre-proline positions, where it is set to cis
/// with 95% probability.  This crudely matches frequencies in PDB structures.
void
trRosettaProtocolMover::randomize_backbone_dihedrals_rama_prepro(
	core::pose::Pose & pose
) const {
	core::chemical::ResidueTypeCOP ala_type(
		core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID_t )->get_representative_type_aa( core::chemical::aa_ala )
	);
	core::scoring::RamaPrePro const rpp; //Genrator for random phi/psi angles.
	core::Real const prepro_trans_prob(1.0 - ramachandran_mode_cis_probability_prepro_);
	core::Real const non_prepro_trans_prob(1.0 - ramachandran_mode_cis_probability_non_prepro_);
	for ( core::Size ir(1), irmax(pose.total_residue()); ir<=irmax; ++ir ) {
		core::Real const randval( numeric::random::uniform() );
		bool do_cis = false;
		if ( ir != irmax && pose.residue_type(ir+1).aa() == core::chemical::aa_pro ) {
			if ( randval > prepro_trans_prob ) {
				do_cis = true;
			}
		} else {
			if ( randval > non_prepro_trans_prob ) {
				do_cis = true;
			}
		}
		pose.set_omega( ir, do_cis ? 0.0 : 180.0 );

		utility::vector1< core::Real > torsions(2);
		rpp.random_mainchain_torsions( pose.residue_type(ir).aa(), (ir < irmax ? pose.residue_type_ptr(ir+1) : ala_type), torsions );
		pose.set_phi( ir, torsions[1] );
		pose.set_psi( ir, torsions[2] );
	}
	pose.update_residue_neighbors();
}

/// @brief Set or randomize the backbone phi, psi, and omega angles based on the probability of observing
/// residue type i in backbone bin X and residue type i+1 in backbone bin Y.
void
trRosettaProtocolMover::randomize_backbone_dihedrals_by_bins(
	core::pose::Pose & pose
) const {
	protocols::simple_moves::bin_transitions::InitializeByBins init_by_bins;
	init_by_bins.set_binfile_and_load( "ABBA" );
	init_by_bins.apply(pose);
}

/// @brief Mutate all the glycines to alanine.
/// @returns A StoredResidueSubsetSelector that preserves the glycine positions, for later mutating
/// them all back to alanines.
core::select::residue_selector::ResidueSelectorCOP
trRosettaProtocolMover::do_mutate_gly_to_ala(
	core::pose::Pose & pose
) const {
	using namespace core::select::residue_selector;
	using namespace protocols::residue_selectors;
	using namespace protocols::simple_moves;

	TR << "Mutating glycine positions to alanine for centroid phase." << std::endl;

	// Select all the gly positions:
	ResidueIndexSelectorOP gly_selector(
		utility::pointer::make_shared< ResidueIndexSelector >()
	);

	TR << "Glycine positions found: ";
	bool first (true);
	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
		if ( pose.residue_type(i).aa() == core::chemical::aa_gly ) {
			gly_selector->append_index(i);
			if ( first ) {
				first = false;
			} else {
				TR << ", ";
			}
			TR << i;
		}
	}
	if ( first ) {
		TR << "[none found]." << std::endl;
	} else {
		TR << "." << std::endl;
	}

	// Store them in the pose:
	StoreResidueSubsetMover store_gly_positions( gly_selector, std::string("glycine_positions"), true );
	store_gly_positions.apply(pose);

	// Mutate them all to alanine:
	MutateResidue mutate_to_ala;
	mutate_to_ala.set_selector( gly_selector );
	mutate_to_ala.set_res_name( "ALA" );
	mutate_to_ala.apply(pose);

	// Make sure that we still have terminal types:
	add_termini( pose );

	// Return a selector that selects the stored glycine positions:
	return utility::pointer::make_shared< StoredResidueSubsetSelector >( "glycine_positions" );
}

/// @brief Mutate the alanines that were previously glycines back to glycine.
void
trRosettaProtocolMover::do_mutate_ala_to_gly(
	core::select::residue_selector::ResidueSelectorCOP gly_selector,
	core::pose::Pose & pose
) const {
	TR << "Mutating positions that were formerly glycine back to glycine." << std::endl;
	protocols::simple_moves::MutateResidue mutate_to_gly;
	mutate_to_gly.set_selector( gly_selector );
	mutate_to_gly.set_res_name( "GLY" );
	mutate_to_gly.apply(pose);
	add_termini(pose);
}

/// @brief Ensure that the termini have the N-terminal and C-terminal protein types.
void
trRosettaProtocolMover::add_termini(
	core::pose::Pose & pose
) const {
	using namespace core::pose;
	add_lower_terminus_type_to_pose_residue( pose, 1 );
	add_upper_terminus_type_to_pose_residue( pose, pose.total_residue() );
	pose.update_residue_neighbors();
}

/// @brief Set energy method options for a scorefunction (to turn on hb_cen_soft).
void
trRosettaProtocolMover::set_emethod_options(
	core::scoring::ScoreFunction & sfxn
) const {
	core::scoring::methods::EnergyMethodOptions opts( sfxn.energy_method_options() );
	opts.hb_cen_soft(true);
	sfxn.set_energy_method_options(opts);
}

/// @brief Perform up to 5 rounds of minimization, or until the score is less than 10.0.
void
trRosettaProtocolMover::remove_clashes(
	core::pose::Pose & pose,
	protocols::moves::Mover & minmover,
	core::scoring::ScoreFunction const & sfxn
) const {
	TR << "Removing clashes." << std::endl;
	for ( core::Size i(1); i<=5; ++i ) {
		if ( (sfxn)(pose) < 10.0 ) break;
		minmover.apply(pose);
	}
}

#endif //USE_TENSORFLOW

std::ostream &
operator<<( std::ostream & os, trRosettaProtocolMover const & mover )
{
	mover.show(os);
	return os;
}


} //movers
} //trRosetta_protocols
} //protocols
