// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/MakeBarrel.cc
/// @brief Builds a beta-barrel using the Crick parameters.
/// @details The barrel is centred on the origin, with the barrel axis pointing along the
/// global z-axis.  This mover calls the MakeBarrelStrand mover.
/// @author Andy Watkins

// Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/beta_barrel/MakeBarrel.hh>
#include <protocols/beta_barrel/MakeBarrelCreator.hh>
#include <protocols/beta_barrel/BarrelParametrizationCalculator.hh>
#include <protocols/beta_barrel/util.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <utility/tag/Tag.hh>

#include <numeric/constants.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <core/types.hh>
#include <core/conformation/parametric/Parameter.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <core/conformation/parametric/BooleanValuedParameter.hh>
#include <utility/pointer/memory.hh>

//Auto Headers
#include <core/pose/Pose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <protocols/beta_barrel/MakeBarrelStrand.hh>
#include <protocols/helical_bundle/util.hh>
#include <core/conformation/Conformation.hh>

namespace protocols {
namespace beta_barrel {

static basic::Tracer TR("protocols.beta_barrel.MakeBarrel");

/// @brief Constructor for MakeBarrel mover.
MakeBarrel::MakeBarrel():
	Mover("MakeBarrel"),
	reset_pose_(true),
	default_calculator_( utility::pointer::make_shared< BarrelParametrizationCalculator >()  ),
	make_barrel_strand_movers_(),
	n_strands_(8),
	shear_number_(8),
	antiparallel_(true),
	default_crick_params_file_("beta_strand"),
	default_residue_name_(),
	default_strand_length_(7),
	use_degrees_(false),
	last_apply_failed_(false),
	defaults_set_(false)
{
	default_residue_name_ = utility::vector1<std::string>({ "ALA" });
	initialize_default_calculator_from_default_crick_params_file();
}


/// @brief Copy constructor for MakeBarrel mover.
MakeBarrel::MakeBarrel( MakeBarrel const & src ):
	protocols::moves::Mover(src),
	reset_pose_(src.reset_pose_),
	default_calculator_( utility::pointer::static_pointer_cast< BarrelParametrizationCalculator >(src.default_calculator_->clone()) ), //Cloned (deep-copied)
	make_barrel_strand_movers_(), //copy this below
	n_strands_(src.n_strands_),
	shear_number_(src.shear_number_),
	antiparallel_(src.antiparallel_),
	default_crick_params_file_(src.default_crick_params_file_),
	default_residue_name_(src.default_residue_name_),
	default_strand_length_(src.default_strand_length_),
	use_degrees_(src.use_degrees_),
	last_apply_failed_(src.last_apply_failed_),
	defaults_set_(src.defaults_set_)
{
	make_barrel_strand_movers_.clear();
	for ( core::Size i=1, imax=src.make_barrel_strand_movers_.size(); i<=imax; ++i ) {
		make_barrel_strand_movers_.push_back( utility::pointer::static_pointer_cast< MakeBarrelStrand >( src.make_barrel_strand_movers_[i]->clone() ) );
	}
}


/// @brief Destructor for MakeBarrel mover.
MakeBarrel::~MakeBarrel() = default;


/// @brief Clone operator to create a pointer to a fresh MakeBarrel object that copies this one.
protocols::moves::MoverOP MakeBarrel::clone() const {
	return utility::pointer::make_shared< MakeBarrel >( *this );
}


/// @brief Fresh_instance operator to create a pointer to a fresh MakeBarrel object that does NOT copy this one.
protocols::moves::MoverOP MakeBarrel::fresh_instance() const {
	return utility::pointer::make_shared< MakeBarrel >();
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void MakeBarrel::apply (core::pose::Pose & pose)
{
	if ( TR.visible() ) TR << "Building a beta-barrel using the Crick equations." << std::endl;

	bool failed=false;

	core::pose::Pose newpose;

	//Compute per-strand positions from global barrel parameters.
	//We need z1 from the default calculator to compute the axial stagger.
	core::Real const z1_val( default_calculator_->real_parameter_cop( BBPC_z1 )->value() );

	utility::vector1< core::Real > delta_omega0_values;
	utility::vector1< core::Real > delta_z0_values;
	utility::vector1< bool > invert_values;
	compute_strand_positions( n_strands_, shear_number_, z1_val, antiparallel_,
		delta_omega0_values, delta_z0_values, invert_values );

	//The BarrelParametersSet object that will hold the Crick parameters and will be passed into the Conformation object of the pose on which we will be operating.
	BarrelParametersSetOP output_parameters_set( utility::pointer::make_shared< BarrelParametersSet >() );
	output_parameters_set->set_n_strands( n_strands_ );
	output_parameters_set->set_shear_number( shear_number_ );
	output_parameters_set->set_antiparallel( antiparallel_ );
	output_parameters_set->set_barrel_radius( default_calculator_->real_parameter_cop( BBPC_r0 )->value() );

	for ( core::Size istrand=1; istrand<=n_strands_; ++istrand ) { //Loop through all strands.
		// Create a MakeBarrelStrand mover for this strand, either from the per-strand override or defaults.
		protocols::beta_barrel::MakeBarrelStrandOP make_strand;
		if ( istrand <= make_barrel_strand_movers_.size() ) {
			make_strand = utility::pointer::static_pointer_cast<protocols::beta_barrel::MakeBarrelStrand>( make_barrel_strand_movers_[istrand]->clone() );
		} else {
			// No per-strand override; create from defaults.
			make_strand = utility::pointer::make_shared< MakeBarrelStrand >( default_calculator_ );
			make_strand->set_residue_name( default_residue_name_ );
			make_strand->set_strand_length( default_strand_length_ );
		}

		make_strand->set_reset_pose(true);
		make_strand->copy_unset_params_from_globals( default_calculator_ );

		// Set per-strand computed values (delta_omega0, delta_z0, invert).
		make_strand->calculator_op()->real_parameter( BBPC_delta_omega0 )->set_value( delta_omega0_values[istrand] );
		make_strand->calculator_op()->real_parameter( BBPC_delta_z0 )->set_value( delta_z0_values[istrand] );
		make_strand->calculator_op()->boolean_parameter( BBPC_invert_strand )->set_value( invert_values[istrand] );

		core::pose::Pose strandpose; //A temporary pose.
		make_strand->apply(strandpose);

		failed=make_strand->last_apply_failed(); //Did the apply fail?
		if ( failed ) break;

		//Add the newly-generated strand to the pose for the barrel.
		if ( newpose.empty() ) newpose=strandpose;
		else newpose.append_pose_by_jump(strandpose,1); //This will also append the ParametersSet object.
	}

	if ( failed ) {
		set_last_apply_failed(true);
		if ( TR.visible() ) TR << "Build of beta-barrel failed, likely due to bad Crick parameters resulting in nonsensical geometry.  Returning input pose." << std::endl;
	} else {
		set_last_apply_failed(false);

		//At this point, newpose has a bunch of ParametersSet objects, each with one Parameters object.  We want a single ParametersSet object with all of the Parameters objects.
		//So we do some shuffling around:
		for ( core::Size i=1, imax=newpose.conformation().n_parameters_sets(); i<=imax; ++i ) { //Loop through all ParametersSets
			output_parameters_set->add_parameters( newpose.conformation().parameters_set(i)->parameters(1) ); //Put all of the Parameters objects into this ParametersSet object.
		}
		newpose.conformation().clear_parameters_set_list(); //Delete all of the ParametersSet objects in newpose.
		newpose.conformation().add_parameters_set( output_parameters_set );

		if ( reset_pose() ) {
			if ( TR.Debug.visible() ) TR.Debug << "Replacing input pose with beta-barrel." << std::endl;
			pose.clear();
			pose=newpose;
		} else {
			if ( TR.Debug.visible() ) TR.Debug << "Appending beta-barrel to pose." << std::endl;
			pose.append_pose_by_jump(newpose, 1);
		}
	}

	if ( TR.Debug.visible() ) TR.Debug << "Finished apply function." << std::endl;

	cyclic_peptide::PeptideStubMover::assign_chain_ids( pose );
}

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////

/// @brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
///
void
MakeBarrel::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data_map*/
) {
	runtime_assert_string_msg( tag->getName() == "MakeBarrel", "This should be impossible -- the tag name does not match the mover name.");

	if ( TR.visible() ) TR << "Parsing options for MakeBarrel (\"" << tag->getOption<std::string>("name" ,"") << "\") mover." << std::endl;

	//Get whether we're using radians or degrees:
	set_use_degrees( tag->getOption<bool>("use_degrees", false) );
	if ( TR.visible() ) TR << "Reading " << (use_degrees() ? "degrees" : "radians") << " for all input angle measurements.  (Radians are used internally, and output is in radians.)" << std::endl;

	//Set the default Crick params file:
	if ( tag->hasOption("crick_params_file") ) {
		set_default_crick_params_file( tag->getOption<std::string>("crick_params_file") );
		if ( TR.visible() ) TR << "Default Crick params file set to " << default_crick_params_file() << "." << std::endl;
	}

	//Set the reset option:
	if ( tag->hasOption("reset") ) {
		set_reset_pose( tag->getOption<bool>("reset") );
		if ( TR.visible() ) TR << "Reset mode set to " << (reset_pose() ? "true" : "false") << "." << std::endl;
	}

	//Set the number of strands:
	if ( tag->hasOption("n_strands") ) {
		set_n_strands( tag->getOption<core::Size>("n_strands") );
		if ( TR.visible() ) TR << "Number of strands set to " << n_strands() << "." << std::endl;
	}

	//Set the shear number:
	if ( tag->hasOption("shear_number") ) {
		set_shear_number( tag->getOption<core::Size>("shear_number") );
		if ( TR.visible() ) TR << "Shear number set to " << shear_number() << "." << std::endl;
	}

	//Set antiparallel:
	if ( tag->hasOption("antiparallel") ) {
		set_antiparallel( tag->getOption<bool>("antiparallel") );
		if ( TR.visible() ) TR << "Antiparallel set to " << (antiparallel() ? "true" : "false") << "." << std::endl;
	}

	//Set defaults for those parameters that can be set from the tag:
	runtime_assert_string_msg( !defaults_set_, "Error in MakeBarrel::parse_my_tag(): The defaults have already been set, and cannot be set again from the tag!" );
	for ( core::Size i(1); static_cast< BBPC_Parameters>(i) < BBPC_end_of_list; ++i ) { //Process all defaults that can be set:
		core::conformation::parametric::ParameterOP curparam( default_calculator_->parameter(i) );
		if ( curparam->can_be_set() ) {
			curparam->parse_setting( tag, curparam->can_be_set(), false, false, false );
		}
	}

	//Set residue name:
	if ( tag->hasOption("residue_name") ) {
		std::string resname_string( tag->getOption<std::string>("residue_name") );
		utility::vector1 < std::string > resname_vect;
		helical_bundle::parse_resnames( resname_string, resname_vect );
		set_default_residue_name( resname_vect );
		if ( TR.visible() ) {
			if ( resname_vect.size()==1 ) TR << "Default residue name set to " << default_residue_name(1) << "." << std::endl;
			else {
				TR << "Default residue names set to ";
				for ( core::Size i=1, imax=resname_vect.size(); i<=imax; ++i ) {
					TR << default_residue_name(i);
					if ( i<imax ) TR << ", ";
				}
				TR << "." << std::endl;
			}
		}
	}

	//Set strand length:
	if ( tag->hasOption("strand_length") ) {
		set_default_strand_length( tag->getOption<core::Size>("strand_length") );
		if ( TR.visible() ) TR << "Default strand length set to " << default_strand_length() << "." << std::endl;
	}

	//AT THIS POINT, DEFAULTS HAVE BEEN SET.
	defaults_set_ = true;

	//Parse sub-tags, setting per-strand params.
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	for ( utility::vector1< utility::tag::TagCOP >::const_iterator branch_tag( branch_tags.begin() ); branch_tag != branch_tags.end(); ++branch_tag ) {
		if ( (*branch_tag)->getName() == "Strand" ) { //A strand has been added.  Add it, and parse its options.
			add_strand(); //This will set all of the parameters to defaults

			core::Size const strand_index( make_barrel_strand_movers_.size() ); //The index of the current strand
			if ( TR.visible() ) TR << "Added a strand." << std::endl;

			//First, set Crick params file for this strand (if specified):
			if ( (*branch_tag)->hasOption( "crick_params_file" ) ) {
				std::string const paramsfile( (*branch_tag)->getOption<std::string>( "crick_params_file", "" ) );
				set_crick_params_file_for_strand( paramsfile, strand_index );
				if ( TR.visible() ) TR << "\tRead minor helix parameters from " << paramsfile << "." << std::endl;
			}

			//Set residue_name for this strand:
			if ( (*branch_tag)->hasOption( "residue_name" ) ) {
				std::string const resname( (*branch_tag)->getOption<std::string>("residue_name", "ALA") );
				utility::vector1< std::string > resnames;
				helical_bundle::parse_resnames(resname, resnames);
				set_residue_name_for_strand( resnames, strand_index );
				if ( TR.visible() ) {
					if ( resnames.size()==1 ) {
						TR << "\tSet the residue type for the strand to " << resnames[1] << "." << std::endl;
					} else {
						TR << "\tSet the following residue types for the repeating unit in this strand: ";
						for ( core::Size i=1, imax=resnames.size(); i<=imax; ++i ) {
							TR << resnames[i];
							if ( i<imax ) TR << ", ";
						}
						TR << "." << std::endl;
					}
				}
			}

			//Set strand length for this strand:
			if ( (*branch_tag)->hasOption( "strand_length" ) ) {
				core::Size const strandlength( (*branch_tag)->getOption<core::Size>("strand_length", core::Size(0)) );
				set_strand_length_for_strand( strandlength, strand_index );
				if ( TR.visible() ) TR << "\tSet strand length to " << strandlength << "." << std::endl;
			}

			//Set parameters for this strand:
			for ( core::Size i(1); static_cast< BBPC_Parameters>(i) < BBPC_end_of_list; ++i ) { //Process all that can be set:
				core::conformation::parametric::ParameterOP curparam( strand(strand_index)->calculator_op()->parameter(i) );
				if ( ( !curparam->global_for_parameters_set() ) && ( curparam->can_be_set() || curparam->can_be_copied() ) ) {
					curparam->parse_setting( (*branch_tag), curparam->can_be_set(), false, false, curparam->can_be_copied() );
				}
			}
		}
	}

	// Note: unlike MakeBundle, we do not require at least one Strand sub-tag, because
	// the barrel auto-generates n_strands strands from the global parameters if no
	// per-strand overrides are provided.

	return;
} //parse_my_tag


/// @brief Function to add a strand.
/// @details This creates a MakeBarrelStrand mover that will be called at apply time.
/// Note that this function assumes that defaults have been set already.
void MakeBarrel::add_strand() {
	defaults_set_ = true; //At this point, we assume that defaults are set.

	make_barrel_strand_movers_.push_back( utility::pointer::make_shared< protocols::beta_barrel::MakeBarrelStrand >( default_calculator_ ) );

	core::Size const newindex( make_barrel_strand_movers_.size() );
	strand(newindex)->set_residue_name( default_residue_name() );
	strand(newindex)->set_strand_length( default_strand_length() );
	return;
}

/// @brief Set the default Crick params file name.
/// @details Triggers a read from disk!
void
MakeBarrel::set_default_crick_params_file(
	std::string const &input_string
) {
	runtime_assert_string_msg( !defaults_set_, "Error in MakeBarrel::set_default_crick_params_file(): An attempt was made to set the default Crick parameters file name, but the defaults have already been set!" );
	default_crick_params_file_=input_string;
	initialize_default_calculator_from_default_crick_params_file();
}

/// @brief Initialize the default calculator from the default Crick params file.
/// @details Triggers a read from disk!
void
MakeBarrel::initialize_default_calculator_from_default_crick_params_file() {
	default_calculator_->init_from_file( default_crick_params_file_ );
}

/// @brief Set the Crick params file for a particular strand.
/// @details Triggers a read from disk!
void
MakeBarrel::set_crick_params_file_for_strand(
	std::string const &filename,
	core::Size const strand_index
) {
	runtime_assert_string_msg( strand_index > 0 && strand_index <= make_barrel_strand_movers_.size(), "Error in MakeBarrel::set_crick_params_file_for_strand(): The strand index provided is outside of the range of strand indices configured for this mover." );
	strand(strand_index)->set_minor_helix_params_from_file( filename );
}

/// @brief Returns the default Crick params file name.
///
std::string const & MakeBarrel::default_crick_params_file() const {
	return default_crick_params_file_;
}

/// @brief Set the default residue name
///
void
MakeBarrel::set_default_residue_name(
	utility::vector1< std::string > const &names
) {
	runtime_assert_string_msg( !defaults_set_, "Error in MakeBarrel::set_default_residue_name(): An attempt was made to set the default residue name(s), but the defaults have already been set!" );
	runtime_assert_string_msg( names.size() > 0, "Error in MakeBarrel::set_default_residue_name():  The residue names vector passed to this function was empty!" );
	default_residue_name_=names;
}


/// @brief Returns the default residue name for a particular index in the repeating unit.
///
std::string const &
MakeBarrel::default_residue_name( core::Size const index_in_repeating_unit ) const {
	runtime_assert_string_msg( index_in_repeating_unit <= default_residue_name_.size() && index_in_repeating_unit > 0, "In protocols::beta_barrel::MakeBarrel::default_residue_name(): The index provided to this function is out of range." );
	return default_residue_name_[index_in_repeating_unit];
}

/// @brief Returns the default residue name vector.
///
utility::vector1 < std::string > const & MakeBarrel::default_residue_name() const {
	return default_residue_name_;
}

/// @brief Set the residue names for a given strand.
void
MakeBarrel::set_residue_name_for_strand(
	utility::vector1< std::string > const & names,
	core::Size const strand_index
) {
	runtime_assert_string_msg( strand_index > 0 && strand_index <= make_barrel_strand_movers_.size(), "Error in MakeBarrel::set_residue_name_for_strand(): The strand index provided is outside of the range of strand indices configured for this mover." );
	runtime_assert_string_msg( names.size() > 0, "Error in MakeBarrel::set_residue_name_for_strand():  The residue names vector passed to this function was empty!" );
	strand(strand_index)->set_residue_name(names);
}

/// @brief Set the default number of residues per strand
///
void
MakeBarrel::set_default_strand_length(
	core::Size const &val
) {
	runtime_assert_string_msg( !defaults_set_, "Error in MakeBarrel::set_default_strand_length(): An attempt was made to set the default strand length, but the defaults have already been set!" );
	runtime_assert_string_msg( val > 1, "Error in MakeBarrel::set_default_strand_length(): The strand length must be greater than 1." );
	default_strand_length_=val;
}

/// @brief Returns the default number of residues in each strand
///
core::Size MakeBarrel::default_strand_length() const {
	return default_strand_length_;
}

/// @brief Set the strand length, in residues, for the Nth strand.
void
MakeBarrel::set_strand_length_for_strand(
	core::Size const strand_length,
	core::Size const strand_index
) {
	runtime_assert_string_msg( strand_length > 1, "Error in MakeBarrel::set_strand_length_for_strand(): The strand length must be greater than 1." );
	runtime_assert_string_msg( strand_index > 0 && strand_index <= make_barrel_strand_movers_.size(), "Error in MakeBarrel::set_strand_length_for_strand(): The strand index provided is outside of the range of strand indices configured for this mover." );
	strand(strand_index)->set_strand_length(strand_length);
}

/// @brief Set whether we're using degrees (true) or radians (false).
///
void
MakeBarrel::set_use_degrees( bool const val/*=true*/ ) {
	use_degrees_=val;
	default_calculator_->set_use_degrees( use_degrees() );
}

std::string MakeBarrel::get_name() const {
	return mover_name();
}

std::string MakeBarrel::mover_name() {
	return "MakeBarrel";
}

void MakeBarrel::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	BarrelParametrizationCalculatorOP default_calculator( utility::pointer::make_shared< BarrelParametrizationCalculator >() );

	using namespace utility::tag;
	AttributeList attlist;

	for ( core::Size i(1); static_cast< BBPC_Parameters >(i) < BBPC_end_of_list; ++i ) { //Process all defaults that can be set:
		core::conformation::parametric::ParameterOP curparam( default_calculator->parameter(i) );
		if ( curparam->can_be_set() ) {
			curparam->provide_xsd_information( attlist, curparam->can_be_set(), false, false, false );
		}
	}

	attlist + XMLSchemaAttribute::attribute_w_default( "use_degrees", xsct_rosetta_bool, "Input values in degrees, instead of radians", "false" );
	attlist + XMLSchemaAttribute::attribute_w_default( "n_strands", xsct_non_negative_integer, "Number of strands in the barrel", "8" );
	attlist + XMLSchemaAttribute::attribute_w_default( "shear_number", xsct_non_negative_integer, "Shear number of the barrel", "8" );
	attlist + XMLSchemaAttribute::attribute_w_default( "antiparallel", xsct_rosetta_bool, "If true, adjacent strands run in opposite directions", "true" );
	attlist + XMLSchemaAttribute::attribute_w_default( "reset", xsct_rosetta_bool, "Reset the input pose, instead of appending the barrel to it", "true" );
	attlist + XMLSchemaAttribute::attribute_w_default( "crick_params_file", xs_string, "Crick params file for strand geometry", "beta_strand" );
	attlist + XMLSchemaAttribute::attribute_w_default( "residue_name", xs_string, "Residue type name (full name, not 3-letter code)", "ALA" );
	attlist + XMLSchemaAttribute::attribute_w_default( "strand_length", xsct_non_negative_integer, "Default number of residues per strand", "7" );

	AttributeList subtag_attributes;
	subtag_attributes + XMLSchemaAttribute( "crick_params_file", xs_string, "Crick params file for this strand" );
	subtag_attributes + XMLSchemaAttribute( "residue_name", xs_string, "Residue type for this strand" );
	subtag_attributes + XMLSchemaAttribute( "strand_length", xsct_non_negative_integer, "Number of residues in this strand" );

	for ( core::Size i(1); static_cast< BBPC_Parameters >(i) < BBPC_end_of_list; ++i ) { //Process all defaults that can be set:
		core::conformation::parametric::ParameterOP curparam( default_calculator->parameter(i) );
		if ( !curparam->global_for_parameters_set() && ( curparam->can_be_set() || curparam->can_be_copied() ) ) {
			curparam->provide_xsd_information( subtag_attributes, curparam->can_be_set(), curparam->can_be_copied(), false, false );
		}
	}


	utility::tag::XMLSchemaSimpleSubelementList ssl;
	ssl.add_simple_subelement( "Strand", subtag_attributes, "Tags describing individual strands in the barrel" );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "The MakeBarrel mover builds a beta-barrel parametrically, using the Crick parameterization, given a set of Crick parameter values.  It arranges strands around a barrel axis with specified radius, shear number, and antiparallel topology.", attlist, ssl );
}

/// @brief Provide the citation.
void
MakeBarrel::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	using namespace basic::citation_manager;
	citations.add(
		utility::pointer::make_shared< UnpublishedModuleInfo >(
		"MakeBarrel", CitedModuleType::Mover,
		"Andy Watkins",
		"Department of Chemical Engineering, Stanford University",
		"andy.watkins2@gmail.com"
		)
	);
}

std::string MakeBarrelCreator::keyname() const {
	return MakeBarrel::mover_name();
}

protocols::moves::MoverOP
MakeBarrelCreator::create_mover() const {
	return utility::pointer::make_shared< MakeBarrel >();
}

void MakeBarrelCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MakeBarrel::provide_xml_schema( xsd );
}


////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////


} //namespace beta_barrel
} //namespace protocols
