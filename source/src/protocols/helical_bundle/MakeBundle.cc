// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/MakeBundle.cc
/// @brief Builds a helical bundle using the Crick parameters.
/// @details The bundle is centred on the origin, with the outer helix axis pointing along the
/// global z-axis.  This mover calls the MakeBundleHelix mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/helical_bundle/MakeBundle.hh>
#include <protocols/helical_bundle/MakeBundleCreator.hh>
#include <protocols/helical_bundle/BundleParametrizationCalculator.hh>
#include <numeric/crick_equations/BundleParams.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <utility/tag/Tag.hh>

#include <numeric/constants.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <numeric/random/random.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <core/scoring/rms_util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/parametric/Parameter.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <utility/pointer/memory.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace helical_bundle {

static basic::Tracer TR("protocols.helical_bundle.MakeBundle");

/// @brief Creator for MakeBundle mover.
MakeBundle::MakeBundle():
	Mover("MakeBundle"),
	reset_pose_(true),
	default_calculator_( utility::pointer::make_shared< BundleParametrizationCalculator >()  ),
	make_bundle_helix_movers_(),
	bundle_symmetry_(0),
	bundle_symmetry_copies_(0),
	default_crick_params_file_("alpha_helix.crick_params"),
	default_residue_name_(),
	default_helix_length_(10),
	default_helix_length_set_(false),
	use_degrees_(false),
	last_apply_failed_(false),
	defaults_set_(false)
{
	default_residue_name_ = utility::vector1<std::string>({ "ALA" });
	initialize_default_calculator_from_default_crick_params_file();
}


/// @brief Copy constructor for MakeBundle mover.
///
/// @brief Creator for MakeBundle mover.
MakeBundle::MakeBundle( MakeBundle const & src ):
	protocols::moves::Mover(src),
	reset_pose_(src.reset_pose_),
	default_calculator_( utility::pointer::static_pointer_cast< BundleParametrizationCalculator >(src.default_calculator_->clone()) ), //Cloned (deep-copied)
	make_bundle_helix_movers_(), //copy this below
	bundle_symmetry_(src.bundle_symmetry_),
	bundle_symmetry_copies_(src.bundle_symmetry_copies_),
	default_crick_params_file_(src.default_crick_params_file_),
	default_residue_name_(src.default_residue_name_),
	default_helix_length_(src.default_helix_length_),
	default_helix_length_set_(src.default_helix_length_set_),
	use_degrees_(src.use_degrees_),
	last_apply_failed_(src.last_apply_failed_),
	defaults_set_(src.defaults_set_)
{
	make_bundle_helix_movers_.clear();
	for ( core::Size i=1, imax=src.make_bundle_helix_movers_.size(); i<=imax; ++i ) {
		make_bundle_helix_movers_.push_back( utility::pointer::static_pointer_cast< MakeBundleHelix >( src.make_bundle_helix_movers_[i]->clone() ) );
	}
}


/// @brief Destructor for MakeBundle mover.
MakeBundle::~MakeBundle() = default;


/// @brief Clone operator to create a pointer to a fresh MakeBundle object that copies this one.
protocols::moves::MoverOP MakeBundle::clone() const {
	return protocols::moves::MoverOP( new MakeBundle( *this ) );
}


/// @brief Fresh_instance operator to create a pointer to a fresh MakeBundle object that does NOT copy this one.
protocols::moves::MoverOP MakeBundle::fresh_instance() const {
	return protocols::moves::MoverOP( new MakeBundle );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void MakeBundle::apply (core::pose::Pose & pose)
{
	if ( TR.visible() ) TR << "Building a helical bundle using the Crick equations." << std::endl;

	bool failed=false;

	core::pose::Pose newpose;

	core::Size const total_repeats(symmetry()<2 ? 1 : symmetry());
	core::Size const repeats(symmetry_copies()< 1 ? total_repeats : symmetry_copies());
	core::Real const delta_omega0_offset_increment((use_degrees() ? 360.0 : numeric::constants::d::pi_2) / static_cast<core::Real>(total_repeats));
	core::Real delta_omega0_offset(0.0);

	//The BundleParametersSet object that will hold the Crick parameters and will be passed into the Conformation object of the pose on which we will be operating.
	BundleParametersSetOP output_parameters_set(BundleParametersSetOP( new BundleParametersSet ));
	output_parameters_set->set_bundle_symmetry( symmetry() ); //Store the symmetry of this bundle.
	output_parameters_set->set_bundle_symmetry_copies( symmetry_copies() ); //Store the number of symmetry copies to actually generate.
	output_parameters_set->set_n_helices( n_helices() ); //Store the number of helices that are defined for each symmetry repeat.

	for ( core::Size irepeat=1; irepeat<=repeats; ++irepeat ) { //Loop through all of the symmetry copies.
		for ( core::Size ihelix=1, ihelixmax=n_helices(); ihelix<=ihelixmax; ++ihelix ) { //Loop through all of the helices defined for each symmetry repeat.
			protocols::helical_bundle::MakeBundleHelixOP make_helix( utility::pointer::static_pointer_cast<protocols::helical_bundle::MakeBundleHelix>(helix(ihelix)->clone()) );
			make_helix->set_reset_pose(true);
			make_helix->copy_unset_params_from_globals( default_calculator_ );
			if ( make_helix->copy_params_from_previous_helices( newpose ) ) {
				failed = true;
				break;
			}
			if ( irepeat > 1 ) make_helix->calculator_op()->real_parameter( BPC_delta_omega0 )->set_value( make_helix->calculator_cop()->real_parameter_cop(BPC_delta_omega0)->value() * ( use_degrees() ? 180.0 / numeric::constants::d::pi : 1.0 ) + delta_omega0_offset ); //Offset omega0 for symmetry copies.
			core::pose::Pose helixpose; //A temporary pose.
			make_helix->apply(helixpose);

			failed=make_helix->last_apply_failed(); //Did the apply fail?
			if ( failed ) break;

			//Add the newly-generated helix to the pose for the bundle.
			if ( newpose.empty() ) newpose=helixpose;
			else newpose.append_pose_by_jump(helixpose,1); //This will also append the ParametersSet object -- but this creates a list of ParametersSets, rather than one ParametersSet with a list of Parameters.
		}
		if ( failed ) break;
		if ( repeats>1 ) delta_omega0_offset += delta_omega0_offset_increment;
	}

	if ( failed ) {
		set_last_apply_failed(true);
		if ( TR.visible() ) TR << "Build of helical bundle failed, likely due to bad Crick parameters resulting in nonsensical geometry.  Returning input pose." << std::endl;
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
			if ( TR.Debug.visible() ) TR.Debug << "Replacing input pose with helical bundle." << std::endl;
			pose.clear();
			pose=newpose;
		} else {
			if ( TR.Debug.visible() ) TR.Debug << "Appending helical bundle to pose." << std::endl;
			pose.append_pose_by_jump(newpose, 1);
		}
	}

	if ( TR.Debug.visible() ) TR.Debug << "Finished apply function." << std::endl;

	return;
}

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////

/// @brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
///
void
MakeBundle::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data_map*/,
	protocols::filters::Filters_map const &/*filters*/,
	protocols::moves::Movers_map const &/*movers*/,
	core::pose::Pose const & /*pose*/
) {
	runtime_assert_string_msg( tag->getName() == "MakeBundle", "This should be impossible -- the tag name does not match the mover name.");

	if ( TR.visible() ) TR << "Parsing options for MakeBundle (\"" << tag->getOption<std::string>("name" ,"") << "\") mover." << std::endl;

	//Get whether we're using radians or degrees:
	set_use_degrees( tag->getOption<bool>("use_degrees", false) );
	if ( TR.visible() ) TR << "Reading " << (use_degrees() ? "degrees" : "radians") << " for all input angle measurements.  (Radians are used internally, and output is in radians.)" << std::endl;

	//Set the default Crick params file:
	if ( tag->hasOption("crick_params_file") ) {
		set_default_crick_params_file( tag->getOption<std::string>("crick_params_file") );
		if ( TR.visible() ) TR << "Default Crick params file set to " << default_crick_params_file() << "." << std::endl;
	}

	//Set symmetry options for the MakeBundle mover:
	set_symmetry_options_from_tag(tag);

	//Set the reset option:
	if ( tag->hasOption("reset") ) {
		set_reset_pose( tag->getOption<bool>("reset") );
		if ( TR.visible() ) TR << "Reset mode set to " << (reset_pose() ? "true" : "false") << "." << std::endl;
	}

	//Set defaults for those parameters that can be set from the tag:
	runtime_assert_string_msg( !defaults_set_, "Error in MakeBundle::parse_my_tag(): The defaults have already been set, and cannot be set again from the tag!" );
	for ( core::Size i(1); static_cast< BPC_Parameters>(i) < BPC_end_of_list; ++i ) { //Process all defaults that can be set:
		core::conformation::parametric::ParameterOP curparam( default_calculator_->parameter(i) );
		if ( curparam->can_be_set() ) {
			curparam->parse_setting( tag, curparam->can_be_set(), false, false, false );
		}
	}

	//Set residue name:
	if ( tag->hasOption("residue_name") ) {
		std::string resname_string( tag->getOption<std::string>("residue_name") );
		utility::vector1 < std::string > resname_vect;
		parse_resnames( resname_string, resname_vect );
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

	//Set helix length:
	if ( tag->hasOption("helix_length") ) {
		set_default_helix_length( tag->getOption<core::Size>("helix_length") );
		if ( TR.visible() ) TR << "Default helix length set to " << default_helix_length() << "." << std::endl;
	}

	//AT THIS POINT, DEFAULTS HAVE BEEN SET.
	defaults_set_ = true;

	//Check that at least one helix is defined:
	bool at_least_one_helix = false;

	//Parse sub-tags, setting major and minor helix params.  (Note that the add_helix function will automatically set the defaults).
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	for ( utility::vector1< utility::tag::TagCOP >::const_iterator branch_tag( branch_tags.begin() ); branch_tag != branch_tags.end(); ++branch_tag ) {
		if ( (*branch_tag)->getName() == "Helix" ) { //A helix or strand (synonymous concepts, here) has been added.  Add it, and parse its options.
			at_least_one_helix = true;
			add_helix(); //This will set all of the parameters to defaults

			core::Size const helix_index( n_helices() ); //The index of the current helix
			if ( TR.visible() ) TR << "Added a helix." << std::endl;

			//First, set Crick params file for this helix (if specified):
			if ( (*branch_tag)->hasOption( "crick_params_file" ) ) {
				std::string const paramsfile( (*branch_tag)->getOption<std::string>( "crick_params_file", "" ) );
				set_crick_params_file_for_helix( paramsfile, helix_index );
				if ( TR.visible() ) TR << "\tRead minor helix parameters from " << paramsfile << "." << std::endl;
			}

			//Set residue_name for this helix:
			if ( (*branch_tag)->hasOption( "residue_name" ) ) {
				std::string const resname( (*branch_tag)->getOption<std::string>("residue_name", "ALA") );
				utility::vector1< std::string > resnames;
				parse_resnames(resname, resnames);
				set_residue_name_for_helix( resnames, helix_index );
				if ( TR.visible() ) {
					if ( resnames.size()==1 ) {
						TR << "\tSet the residue type for the helix to " << resnames[1] << "." << std::endl;
					} else {
						TR << "\tSet the following residue types for the repeating unit in this helix: ";
						for ( core::Size i=1, imax=resnames.size(); i<=imax; ++i ) {
							TR << resnames[i];
							if ( i<imax ) TR << ", ";
						}
						TR << "." << std::endl;
					}
				}
			}

			//Set helix length for this helix:
			if ( (*branch_tag)->hasOption( "helix_length" ) ) {
				core::Size const helixlength( (*branch_tag)->getOption<core::Size>("helix_length", 0) );
				set_helix_length_for_helix( helixlength, helix_index );
				if ( TR.visible() ) TR << "\tSet helix length to " << helixlength << "." << std::endl;
			}

			//Set parameters for this helix:
			for ( core::Size i(1); static_cast< BPC_Parameters>(i) < BPC_end_of_list; ++i ) { //Process all defaults that can be set:
				core::conformation::parametric::ParameterOP curparam( helix(helix_index)->calculator_op()->parameter(i) );
				if ( ( !curparam->global_for_parameters_set() ) && ( curparam->can_be_set() || curparam->can_be_copied() ) ) {
					curparam->parse_setting( (*branch_tag), curparam->can_be_set(), false, false, curparam->can_be_copied() );
				}
			}
		}
	}

	runtime_assert_string_msg(at_least_one_helix, "In protocols::helical_bundle::MakeBundle::parse_my_tag(): At least one helix must be defined using a <Helix ...> sub-tag!");

	return;
} //parse_my_tag

/// @brief Set symmetry and symmetry_copies options based on an input tag.
///
void MakeBundle::set_symmetry_options_from_tag( utility::tag::TagCOP tag )
{
	//Set options for the MakeBundle mover:
	if ( tag->hasOption("symmetry") ) {
		set_symmetry( tag->getOption<core::Size>("symmetry", 0) );
		if ( TR.visible() ) {
			if ( symmetry()<2 ) TR << "Symmetry mode set to false." << std::endl;
			else TR << symmetry() << "-fold symmetry set." << std::endl;
		}
	}
	if ( tag->hasOption("symmetry_copies") ) {
		set_symmetry_copies( tag->getOption<core::Size>("symmetry_copies", 0) );
		if ( TR.visible() ) {
			if ( symmetry_copies()<1 ) TR << "All symmetry copies will be generated." << std::endl;
			if ( symmetry_copies()==1 ) TR << "Only the first symmetry copy will be generated." << std::endl;
			else {
				TR << symmetry_copies() << " symmetry copies will be generated." << std::endl;
				if ( symmetry_copies() > symmetry() ) {
					TR << "Note!  The number of symmetry copies is greater than the bundle symmetry.  This only really makes sense for helical repeat symmetries (i.e. symmetries with z-offset)." << std::endl;
				}
			}
		}
	}

	return;
} //set_global_options_from_tag


/// @brief Function to add a helix.
/// @details This creates a MakeBundleHelix mover that will be called at apply time.  Note that this function assumes that defaults have been set already.  They
/// cannot be set after calling this function.
void MakeBundle::add_helix() {
	defaults_set_ = true; //At this point, we assume that defaults are set.  Setting this to true prevents defaults from being set later.

	make_bundle_helix_movers_.push_back( protocols::helical_bundle::MakeBundleHelixOP(new protocols::helical_bundle::MakeBundleHelix( default_calculator_ )) );

	core::Size const newindex( make_bundle_helix_movers_.size() );
	helix(newindex)->set_residue_name( default_residue_name() );
	helix(newindex)->set_helix_length( default_helix_length() );
	return;
}

/// @brief Set the default Crick params file name.
/// @details Triggers a read from disk!
void
MakeBundle::set_default_crick_params_file(
	std::string const &input_string
) {
	runtime_assert_string_msg( !defaults_set_, "Error in MakeBundle::set_default_crick_params_file(): An attempt was made to set the default Crick parameters file name, but the defaults have already been set!" );
	default_crick_params_file_=input_string;
	initialize_default_calculator_from_default_crick_params_file();
}

/// @brief Initialize the default calculator from the default Crick params file.
/// @details Triggers a read from disk!
void
MakeBundle::initialize_default_calculator_from_default_crick_params_file() {
	default_calculator_->init_from_file( default_crick_params_file_ );
}

/// @brief Set the Crick params file for a particular helix.
/// @details Triggers a read from disk!
void
MakeBundle::set_crick_params_file_for_helix(
	std::string const &filename,
	core::Size const helix_index
) {
	runtime_assert_string_msg( helix_index > 0 && helix_index <= n_helices(), "Error in MakeBundle::set_crick_params_file_for_helix(): The helix index provided is outside of the range of helix indices configured for this mover." );
	helix(helix_index)->set_minor_helix_params_from_file( filename );
}

/// @brief Returns the default Crick params file name.
///
std::string const & MakeBundle::default_crick_params_file() const {
	return default_crick_params_file_;
}

/// @brief Set the default residue name
///
void
MakeBundle::set_default_residue_name(
	utility::vector1< std::string > const &names
) {
	runtime_assert_string_msg( !defaults_set_, "Error in MakeBundle::set_default_residue_name(): An attempt was made to set the default residue name(s), but the defaults have already been set!" );
	runtime_assert_string_msg( names.size() > 0, "Error in MakeBundle::set_default_residue_name():  The residue names vector passed to this function was empty!" );
	default_residue_name_=names;
}


/// @brief Returns the default residue name for a particular index in the repeating unit.
///
std::string const &
MakeBundle::default_residue_name( core::Size const index_in_repeating_unit ) const {
	runtime_assert_string_msg( index_in_repeating_unit <= default_residue_name_.size() && index_in_repeating_unit > 0, "In protocols::helical_bundle::MakeBundle::default_residue_name(): The index provided to this function is out of range." );
	return default_residue_name_[index_in_repeating_unit];
}

/// @brief Returns the default residue name vector.
///
utility::vector1 < std::string > const & MakeBundle::default_residue_name() const {
	return default_residue_name_;
}

/// @brief Set the residue names for a given helix.
void
MakeBundle::set_residue_name_for_helix(
	utility::vector1< std::string > const & names,
	core::Size const helix_index
) {
	runtime_assert_string_msg( helix_index > 0 && helix_index <= n_helices(), "Error in MakeBundle::set_residue_name_for_helix(): The helix index provided is outside of the range of helix indices configured for this mover." );
	runtime_assert_string_msg( names.size() > 0, "Error in MakeBundle::set_residue_name_for_helix():  The residue names vector passed to this function was empty!" );
	helix(helix_index)->set_residue_name(names);
}

/// @brief Set the default number of residues per helix
///
void
MakeBundle::set_default_helix_length(
	core::Size const &val
) {
	runtime_assert_string_msg( !defaults_set_, "Error in MakeBundle::set_default_helix_length(): An attempt was made to set the default helix length, but the defaults have already been set!" );
	runtime_assert_string_msg( val > 1, "Error in MakeBundle::set_default_helix_length(): The helix length must be greater than 1." );
	default_helix_length_=val;
}

/// @brief Returns the default number of residues in each helix
///
core::Size MakeBundle::default_helix_length() const {
	return default_helix_length_;
}

/// @brief Set the helix length, in residues, for the Nth helix.
void
MakeBundle::set_helix_length_for_helix(
	core::Size const helix_length,
	core::Size const helix_index
) {
	runtime_assert_string_msg( helix_length > 1, "Error in MakeBundle::set_helix_length_for_helix(): The helix length must be greater than 1." );
	runtime_assert_string_msg( helix_index > 0 && helix_index <= n_helices(), "Error in MakeBundle::set_helix_length_for_helix(): The helix index provided is outside of the range of helix indices configured for this mover." );
	helix(helix_index)->set_helix_length(helix_length);
}

/// @brief Set whether we're using degrees (true) or radians (false).
///
void
MakeBundle::set_use_degrees( bool const val/*=true*/ ) {
	use_degrees_=val;
	default_calculator_->set_use_degrees( use_degrees() );
}

std::string MakeBundle::get_name() const {
	return mover_name();
}

std::string MakeBundle::mover_name() {
	return "MakeBundle";
}

void MakeBundle::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	BundleParametrizationCalculatorOP default_calculator( utility::pointer::make_shared< BundleParametrizationCalculator >() );

	using namespace utility::tag;
	AttributeList attlist;

	for ( core::Size i(1); static_cast< BPC_Parameters >(i) < BPC_end_of_list; ++i ) { //Process all defaults that can be set:
		core::conformation::parametric::ParameterOP curparam( default_calculator->parameter(i) );
		if ( curparam->can_be_set() ) {
			curparam->provide_xsd_information( attlist, curparam->can_be_set(), false, false, false );
		}
	}

	attlist + XMLSchemaAttribute::attribute_w_default( "use_degrees", xsct_rosetta_bool, "Input values in degrees, instead of radians", "false" );
	add_attributes_for_make_bundle_symmetry( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "reset", xsct_rosetta_bool, "Reset the input pose, instead of appending the bundle to it", "false" );
	add_attributes_for_make_bundle_minorhelix_defaults( attlist );
	add_attributes_for_other_helix_params( attlist );

	AttributeList subtag_attributes;
	add_attributes_for_helix_params( subtag_attributes );
	add_attributes_for_minor_helix_params( subtag_attributes );
	add_attributes_for_other_helix_params( subtag_attributes );

	for ( core::Size i(1); static_cast< BPC_Parameters >(i) < BPC_end_of_list; ++i ) { //Process all defaults that can be set:
		core::conformation::parametric::ParameterOP curparam( default_calculator->parameter(i) );
		if ( !curparam->global_for_parameters_set() && ( curparam->can_be_set() || curparam->can_be_copied() ) ) {
			curparam->provide_xsd_information( subtag_attributes, curparam->can_be_set(), curparam->can_be_copied(), false, false );
		}
	}


	utility::tag::XMLSchemaSimpleSubelementList ssl;
	ssl.add_simple_subelement( "Helix", subtag_attributes, "Tags describing individual helices in the bundle"/*, 0 minoccurs*/ );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "The MakeBundle mover builds a helical bundle parametrically, using the Crick parameterization, given a set of Crick parameter values.  Note that the Crick parameterization is compatible with arbitrary helices (including strands, which are special cases of helices in which the turn per residue is about 180 degrees).", attlist, ssl );
}

std::string MakeBundleCreator::keyname() const {
	return MakeBundle::mover_name();
}

protocols::moves::MoverOP
MakeBundleCreator::create_mover() const {
	return protocols::moves::MoverOP( new MakeBundle );
}

void MakeBundleCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MakeBundle::provide_xml_schema( xsd );
}


////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////


} //namespace helical_bundle
} //namespace protocols
