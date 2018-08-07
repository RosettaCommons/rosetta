// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/PerturbBundle.cc
/// @brief  Perturbs a helical bundle by altering the Crick parameters.
/// @details The bundle is centred on the origin, with the outer helix axis pointing along the
/// global z-axis.  This mover calls the PerturbBundleHelix mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/helical_bundle/PerturbBundle.hh>
#include <protocols/helical_bundle/PerturbBundleCreator.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
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
#include <protocols/helical_bundle/BundleParametrizationCalculator.hh>
#include <protocols/helical_bundle/PerturbBundleHelix.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <utility/pointer/memory.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// STL headers:
#include <iostream>

using basic::Error;
using basic::Warning;

namespace protocols {
namespace helical_bundle {

static basic::Tracer TR("protocols.helical_bundle.PerturbBundle");


/// @brief Creator for PerturbBundle mover.
PerturbBundle::PerturbBundle():
	Mover("PerturbBundle"),
	default_calculator_( new BundleParametrizationCalculator(false) ),
	individual_helix_calculators_(),
	bundleparametersset_index_(1)
{}


/// @brief Copy constructor for PerturbBundle mover.
PerturbBundle::PerturbBundle( PerturbBundle const & src ):
	protocols::moves::Mover( src ),
	default_calculator_( utility::pointer::static_pointer_cast< BundleParametrizationCalculator >( src.default_calculator_->clone() ) ),
	individual_helix_calculators_( /*Cloned below*/),
	bundleparametersset_index_(src.bundleparametersset_index_)
{
	// Deep-cloning the individual helix calculators:
	for ( core::Size i(1), imax(src.individual_helix_calculators_.size()); i<=imax; ++i ) {
		individual_helix_calculators_.push_back( std::pair<core::Size, BundleParametrizationCalculatorOP>( src.individual_helix_calculators_[i].first, utility::pointer::static_pointer_cast< BundleParametrizationCalculator >( src.individual_helix_calculators_[i].second->clone() ) ) );
	}
}


/// @brief Destructor for PerturbBundle mover.
PerturbBundle::~PerturbBundle() = default;

/// @brief Clone operator to create a pointer to a fresh PerturbBundle object that copies this one.
protocols::moves::MoverOP PerturbBundle::clone() const {
	return protocols::moves::MoverOP( new PerturbBundle( *this ) );
}


/// @brief Fresh_instance operator to create a pointer to a fresh PerturbBundle object that does NOT copy this one.
protocols::moves::MoverOP PerturbBundle::fresh_instance() const {
	return protocols::moves::MoverOP( new PerturbBundle );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void PerturbBundle::apply( core::pose::Pose & pose )
{
	bool failed=false;

	if ( TR.visible() ) TR << "Copying pose." << std::endl;
	core::pose::Pose pose_copy(pose);

	if ( TR.visible() ) TR << "Finding BundleParametersSet object in pose." << std::endl;
	BundleParametersSetOP params_set(nullptr);
	core::Size params_set_index(0);
	core::Size n_encountered(0);
	bool breaknow(false);
	for ( core::Size i=1, imax=pose_copy.conformation().n_parameters_sets(); i<=imax; ++i ) {
		BundleParametersSetOP cur_set( utility::pointer::dynamic_pointer_cast< BundleParametersSet >( pose_copy.conformation().parameters_set(i) ) );
		if ( cur_set != nullptr ) {
			++n_encountered; //Increment the number of parameterssets encountered
			if ( n_encountered==bundleparametersset_index() ) {
				params_set_index=i;
				params_set=cur_set; //Assign the current owning pointer to be the params_set owning pointer.
				breaknow=true;
			}
		}
		if ( breaknow ) break;
	}
	runtime_assert_string_msg(params_set_index!=0 && params_set != nullptr, "In protocols::helical_bundle::PerturbBundle::apply() function: BundleparametersSet object with given index not found in pose!");

	write_report( params_set, true); //Write a pre-perturbation report summarizing the initial Crick parameter values.

	if ( TR.visible() ) TR << "Perturbing Crick equation values." << std::endl;
	if ( !perturb_values( params_set ) ) {
		if ( TR.visible() ) TR << "Perturbation failed -- senseless values resulted.  Returning input pose." << std::endl;
		failed=true;
	}

	if ( !failed ) {
		if ( TR.visible() ) TR << "Rebuilding the helical bundle conformation using the Crick parameters stored in the pose." << std::endl;
		rebuild_conformation(pose_copy, params_set, params_set_index, failed);
	}

	if ( !failed ) {
		if ( TR.visible() ) TR << "Perturbation successful.  Copying result to the input pose." << std::endl;
		pose = pose_copy;
		write_report( utility::pointer::dynamic_pointer_cast< BundleParametersSet >( pose.conformation().parameters_set(params_set_index) ), false); //Write a post-perturbation report summarizing the final Crick parameter values.
	} else {
		if ( TR.visible() ) TR << "The current attempt generated Crick parameters that did not permit sensible geometry.  Returning input pose." << std::endl;
	}

	if ( TR.Debug.visible() ) TR.Debug << "Finished apply function." << std::endl;

	TR.flush();
	TR.Debug.flush();

	return;
}

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////

/// @brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
///
void
PerturbBundle::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data_map*/,
	protocols::filters::Filters_map const &/*filters*/,
	protocols::moves::Movers_map const &/*movers*/,
	core::pose::Pose const & /*pose*/
) {

	if ( tag->getName() != "PerturbBundle" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible -- the tag name does not match the mover name.");
	}

	if ( TR.visible() ) TR << "Parsing options for PerturbBundle (\"" << tag->getOption<std::string>("name" ,"") << "\") mover." << std::endl;

	debug_assert( default_calculator_ != nullptr ); //Should have been created.

	//Determine whether input is in degrees or radians
	default_calculator_->set_use_degrees( tag->getOption<bool>( "use_degrees", false ) );
	if ( TR.visible() ) TR << "Interpreting user-input angles as being in " << ( default_calculator_->use_degrees() ? "degrees" : "radians") << ".  (Internally, radians are always used, and output will be in radians.)" << std::endl;

	//Set a default perturbation type:
	if ( tag->hasOption("default_perturbation_type") ) {
		std::string const perttype( tag->getOption<std::string>("default_perturbation_type", "") );
		default_calculator_->set_perturbation_type_globally( perttype );
	}

	//Set defaults for the various perturbable degrees of freedom:
	for ( core::Size i(1); i<=static_cast<core::Size>( BPC_last_parameter_to_be_sampled ); ++i ) {
		core::conformation::parametric::RealValuedParameterOP curparam( default_calculator_->real_parameter( i ) );
		if ( curparam == nullptr ) continue;
		if ( curparam->can_be_set() || curparam->can_be_perturbed() ) {
			curparam->parse_setting( tag, curparam->can_be_set(), false, curparam->can_be_perturbed(), false );
		}
	}

	//Parse options for specific helices:
	reset_helices(); //Make sure we have no defined helices.
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	for ( auto const & branch_tag : branch_tags ) {
		runtime_assert_string_msg( branch_tag->getName() == "Helix",
			"Error in PerturbBundle::parse_my_tag(): Sub-tags for the PerturbBundle mover must be \"Helix\".  Could not parse \"" + branch_tag->getName() + "\"." );
		core::Size helix_index( branch_tag->getOption<core::Size>("helix_index", 0) );
		runtime_assert_string_msg(helix_index>0, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was added, but its index was set to 0.  This is not allowed." );
		BundleParametrizationCalculatorOP this_calculator( add_helix(helix_index) ); //Add and initialize this helix.  By default, no degrees of freedom may be perturbed.  Store the current helix index in this_helix.
		for ( core::Size i(1); i<=static_cast<core::Size>( BPC_last_parameter_to_be_sampled ); ++i ) {
			core::conformation::parametric::RealValuedParameterOP curparam( this_calculator->real_parameter( i ) );
			if ( curparam == nullptr ) continue;
			if ( curparam->can_be_set() || curparam->can_be_perturbed() || curparam->can_be_copied() ) {
				curparam->parse_setting( branch_tag, curparam->can_be_set(), false, curparam->can_be_perturbed(), curparam->can_be_copied() );
			}
		}
	}

} //parse_my_tag

////////////////////////////////////////////////////////////////////////////////
//          PUBLIC FUNCTIONS                                                  //
////////////////////////////////////////////////////////////////////////////////

/// @brief Clear the list of helices.
/// @details Clears the individual_helix_calculators_ list.
void
PerturbBundle::reset_helices() {
	individual_helix_calculators_.clear();
}

/// @brief Add options for a new helix
///
BundleParametrizationCalculatorOP
PerturbBundle::add_helix( core::Size const helix_index ) {
	debug_assert( default_calculator_ != nullptr ); //Shouldn't be possible.
	runtime_assert_string_msg(helix_not_defined( helix_index ), "In protocols::helical_bundle::PerturbBundle::add_helix() function: could not add helix.  Helix defined multiple times!" );
	individual_helix_calculators_.push_back( std::make_pair( helix_index, utility::pointer::static_pointer_cast< BundleParametrizationCalculator >( default_calculator_->clone() ) ) );

	return individual_helix_calculators_[individual_helix_calculators_.size()].second;
}

/// @brief Access the calculator for a given helix (const access).
BundleParametrizationCalculatorCOP
PerturbBundle::individual_helix_calculator_cop(
	core::Size const helix_calculator_index
) const {
	debug_assert( helix_calculator_index > 0 );
	debug_assert( helix_calculator_index <= individual_helix_calculators_.size() );
	for ( core::Size i(1), imax(individual_helix_calculators_.size()); i<=imax; ++i ) {
		if ( individual_helix_calculators_[i].first == helix_calculator_index ) {
			return individual_helix_calculators_[i].second;
		}
	}
	utility_exit_with_message("Error in PerturbBundle::individual_helix_calculator_cop(): Helix with given index not found!");
	return BundleParametrizationCalculatorCOP( nullptr ); //To keep compiler happy.
}

/// @brief Access the calculator for a given helix (nonconst access).
BundleParametrizationCalculatorOP
PerturbBundle::individual_helix_calculator(
	core::Size const helix_calculator_index
) {
	debug_assert( helix_calculator_index > 0 );
	debug_assert( helix_calculator_index <= individual_helix_calculators_.size() );
	for ( core::Size i(1), imax(individual_helix_calculators_.size()); i<=imax; ++i ) {
		if ( individual_helix_calculators_[i].first == helix_calculator_index ) {
			return individual_helix_calculators_[i].second;
		}
	}
	utility_exit_with_message("Error in PerturbBundle::individual_helix_calculator(): Helix with given index not found!");
	return BundleParametrizationCalculatorOP( nullptr ); //To keep compiler happy.
}

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////


/// @brief Confirms that a helix has not yet been defined.  Returns "true" if the helix
/// has NOT been defined, false otherwise.
bool
PerturbBundle::helix_not_defined(
	core::Size const helix_index
) const {
	for ( core::Size i(1), imax(individual_helix_calculators_.size()); i<=imax; ++i ) {
		if ( individual_helix_calculators_[i].first == helix_index ) return false;
	}
	return true;
}

/// @brief Get a calculator for a particular helix.
/// @details Returns nullptr if no calculator for this helix has been defined.
BundleParametrizationCalculatorCOP
PerturbBundle::get_calculator_for_helix(
	core::Size const helix_index
) const {
	for ( core::Size i(1), imax(individual_helix_calculators_.size()); i<=imax; ++i ) {
		if ( individual_helix_calculators_[i].first == helix_index ) return individual_helix_calculators_[i].second;
	}
	return nullptr;
}

/// @brief Perturb the helical bundle parameter values in the pose, subject to the options already set.
/// @details Called by the apply() function.  Returns true for success, false for failure.
bool PerturbBundle::perturb_values( BundleParametersSetOP params_set) const {

	//Get the symmetry of the bundle:
	core::Size const symmetry( params_set->bundle_symmetry() < 2 ? 1 : params_set->bundle_symmetry() );
	//Get the number of symmetry copies in the bundle:
	core::Size const symmetry_copies( params_set->bundle_symmetry_copies()==0 ? symmetry : params_set->bundle_symmetry_copies() );
	//Get the number of helices defined in each symmetry copy in the bundle:
	core::Size const n_helices( params_set->n_helices() );
	runtime_assert_string_msg( n_helices > 0,
		"In protocols::helical_bundle::PerturbBundle::perturb_values() function: the number of helices in the pose is 0.  Unable to proceed -- nothing to perturb." );

	if ( TR.Debug.visible() ) TR.Debug << "params_set->n_parameters(): " << params_set->n_parameters() << std::endl; //DELETE ME
	if ( TR.Debug.visible() ) TR.Debug << "n_helices: " << n_helices << std::endl; //DELETE ME
	if ( TR.Debug.visible() ) TR.Debug << "symmetry_copies: " << symmetry_copies << std::endl; //DELETE ME
	runtime_assert_string_msg( params_set->n_parameters() == n_helices*symmetry_copies,
		"In protocols::helical_bundle::PerturbBundle::perturb_values() function: the pose has corrupted BundleParametersSet data.  Unable to proceed." );

	core::Size helix_index(0); //Counter for the index of the current helix.
	core::Real delta_omega0_offset(0.0); //In radians.
	core::Real delta_omega0_offset_increment( symmetry_copies > 1 ? numeric::constants::d::pi_2 / static_cast<core::Real>(symmetry) : 0.0 );

	bool loopthrough_failed(false); //Used to break from the nested loops.
	for ( core::Size isym=1; isym<=symmetry_copies; ++isym ) { //Loop through all of the symmetry copies.
		for ( core::Size ihelix=1; ihelix<=n_helices; ++ihelix ) { //Loop through all of the helices defined for this symmetry copy
			++helix_index; //Increment the index of the current helix.
			BundleParametersOP params( utility::pointer::dynamic_pointer_cast< BundleParameters >( params_set->parameters(helix_index) ) );
			runtime_assert_string_msg( params != nullptr,
				"In protocols::helical_bundle::PerturbBundle::perturb_values() function: unable to get an owning pointer for the BundleParameters object." );

			if ( helix_index <= n_helices ) { //If this is part of the first symmetry repeat, perturb the values
				protocols::helical_bundle::BundleParametrizationCalculatorCOP curcalculator( get_calculator_for_helix(helix_index) ); //Will be nullptr if helix not defined.
				if ( curcalculator == nullptr ) curcalculator = default_calculator_;
				for ( core::Size iparam(1); iparam<=static_cast<core::Size>(BPC_last_parameter_to_be_sampled); ++iparam ) {
					core::conformation::parametric::RealValuedParameterCOP curparam_calculator( curcalculator->real_parameter_cop( iparam ) );
					if ( curparam_calculator == nullptr ) continue;
					if ( !curparam_calculator->can_be_perturbed() && !curparam_calculator->can_be_copied() && !curparam_calculator->can_be_set() ) continue;
					if ( !curparam_calculator->perturbation_set() && !curparam_calculator->copying_information_was_set() && !curparam_calculator->value_was_set() ) continue;
#ifdef NDEBUG
					//Release mode: no dynamic_cast
					core::conformation::parametric::RealValuedParameterOP curparam( utility::pointer::static_pointer_cast<core::conformation::parametric::RealValuedParameter>( params->parameter_op( iparam ) ) );
#else
					//Debug mode: use dynamic_cast and check.
					core::conformation::parametric::RealValuedParameterOP curparam( utility::pointer::dynamic_pointer_cast<core::conformation::parametric::RealValuedParameter>( params->parameter_op( iparam ) ) );
					debug_assert(curparam != nullptr);
#endif
					if ( curparam_calculator->can_be_perturbed() && curparam_calculator->perturbation_set() ) {
						curparam->set_value( curparam_calculator->generate_perturbed_value( curparam->value() ), true );
					} else if ( curparam_calculator->can_be_copied() && curparam_calculator->copying_information_was_set() ) {
						runtime_assert_string_msg( curparam_calculator->copy_from_parameters_index() < helix_index, "Error in PerturbBundle::perturb_values(): The index of the helix from which we should be copying is not less than the index of the helix to which we're copying parameter values." );
						core::conformation::parametric::RealValuedParameterOP curparam_calculator_copy( utility::pointer::static_pointer_cast< core::conformation::parametric::RealValuedParameter >( curparam_calculator->clone() ) );
						loopthrough_failed = protocols::helical_bundle::BundleParametrizationCalculator::copy_params_from_previous_helices_perturbbundle_style( params_set, curparam_calculator->copy_from_parameters_index(), static_cast<BPC_Parameters>(iparam), curparam_calculator_copy, params );
						if ( loopthrough_failed ) break;
						params->replace_parameter_via_clone( iparam, curparam_calculator_copy );
					} else if ( curparam_calculator->can_be_set() && curparam_calculator->value_was_set() ) {
						curparam->set_value( curparam_calculator->value(), true );
					}
				}
			} else { //If this is part of a later symmetry repeat, just copy the values from the first symmetry repeat.
#ifdef NDEBUG
				BundleParametersOP ref_params( utility::pointer::static_pointer_cast< BundleParameters >( params_set->parameters(ihelix) ) );
				BundleParametersOP cur_params( utility::pointer::static_pointer_cast< BundleParameters >( params_set->parameters(helix_index) ) );
#else
				BundleParametersOP ref_params( utility::pointer::dynamic_pointer_cast< BundleParameters >( params_set->parameters(ihelix) ) );
				BundleParametersOP cur_params( utility::pointer::dynamic_pointer_cast< BundleParameters >( params_set->parameters(helix_index) ) );
				runtime_assert_string_msg( ref_params,
					"In protocols::helical_bundle::PerturbBundle::perturb_values() function: Unable to get an owning pointer for the reference BundleParameters object.  This is odd -- it should not happen." );
				runtime_assert_string_msg( cur_params,
					"In protocols::helical_bundle::PerturbBundle::perturb_values() function: Unable to get an owning pointer for the symmetry copy of a BundleParameters object.  This is odd -- it should not happen." );
#endif
				// Copy parameters:
				for ( core::Size iparam(1); iparam<=static_cast<core::Size>(BPC_last_parameter_to_be_sampled); ++iparam ) {
					core::conformation::parametric::RealValuedParameterOP refparam( utility::pointer::dynamic_pointer_cast<core::conformation::parametric::RealValuedParameter>( ref_params->parameter_op( iparam ) ) );
					debug_assert(refparam != nullptr);
#ifdef NDEBUG
					core::conformation::parametric::RealValuedParameterOP curparam( utility::pointer::static_pointer_cast<core::conformation::parametric::RealValuedParameter>( cur_params->parameter_op( iparam ) ) );
#else
					core::conformation::parametric::RealValuedParameterOP curparam( utility::pointer::dynamic_pointer_cast<core::conformation::parametric::RealValuedParameter>( cur_params->parameter_op( iparam ) ) );
					debug_assert(curparam != nullptr);
#endif
					if ( iparam == static_cast<core::Size>(BPC_delta_omega0) ) {
						curparam->set_value( refparam->value() + ( delta_omega0_offset ), true );
					} else {
						curparam->set_value( refparam->value(), true );
					}
				}
			}
			if ( loopthrough_failed ) break;
		}
		if ( loopthrough_failed ) break;
		delta_omega0_offset += delta_omega0_offset_increment;
	}

	return (!loopthrough_failed);
}

/// @brief Rebuild the helical bundle conformation using the bundle parameter values in the pose.
///
void PerturbBundle::rebuild_conformation(
	core::pose::Pose &pose,
	BundleParametersSetOP params_set,
	core::Size const params_set_index,
	bool &failed
) const {

	core::Size const n_params( params_set->n_parameters() ); //Get the number of parameters objects (number of helices that we're rebuilding).

	for ( core::Size ihelix=1; ihelix<=n_params; ++ihelix ) {

		PerturbBundleHelix perthelix; //Construct the mover to perturb this helix.

		perthelix.set_parameters_set_index(params_set_index);
		perthelix.set_parameters_index(ihelix);

		perthelix.apply( pose );

		failed = perthelix.last_apply_failed();
		if ( failed ) break;

	}

	return;
}

/// @brief Write out the perturbed Crick parameters.
/// @details The "before" parameter determines whether this is a pre-perturbation
/// report or a post-perturbation report.
void PerturbBundle::write_report(
	BundleParametersSetOP params_set,
	bool const before
) const {
	if ( !TR.visible() ) return; //Do nothing if the tracer isn't visible.

	if ( before ) {
		TR << "*** Crick parameter values prior to perturbation ***" << std::endl;
	} else {
		TR << "*** Crick parameter values after perturbation ***" << std::endl;
	}

	core::Size const nparams(params_set->n_parameters());
	if ( nparams==0 ) {
		TR << "No Parameters objects in the ParametersSet object!" << std::endl;
	}

	for ( core::Size iparams=1; iparams<=nparams; ++iparams ) { //Loop through the params objects
		BundleParametersOP params( utility::pointer::dynamic_pointer_cast<BundleParameters>(params_set->parameters(iparams)) );
		if ( params == nullptr ) continue;
		TR << "Helix " << iparams << ":" << std::endl;
		std::stringstream remark;
		params->get_pdb_remark( remark );
		TR << remark.str() << std::endl;
	}

	TR << "*** END SUMMARY ***" << std::endl;
	TR.flush();

	return;
}

std::string PerturbBundle::get_name() const {
	return mover_name();
}

std::string PerturbBundle::mover_name() {
	return "PerturbBundle";
}

void PerturbBundle::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	BundleParametrizationCalculatorOP default_calculator( utility::pointer::make_shared< BundleParametrizationCalculator >() );

	XMLSchemaRestriction pert_type;
	pert_type.name( "pert_type" );
	pert_type.base_type( xs_string );
	pert_type.add_restriction( xsr_enumeration, "gaussian" );
	pert_type.add_restriction( xsr_enumeration, "uniform" );
	pert_type.add_restriction( xsr_enumeration, "" ); // it can be the empty string
	xsd.add_top_level_element( pert_type );

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "use_degrees", xsct_rosetta_bool, "Interpret user-supplied angles as degrees rather than radians", "false" )
		+ XMLSchemaAttribute( "default_perturbation_type", "pert_type", "Default type for perturbations to the bundle parameters, either uniform or gaussian" );

	for ( core::Size i(1); static_cast< BPC_Parameters >(i) < BPC_end_of_list; ++i ) { //Process all defaults that can be set:
		core::conformation::parametric::ParameterCOP curparam( default_calculator->parameter_cop(i) );
		if ( curparam->can_be_set() || curparam->can_be_perturbed() ) {
			curparam->provide_xsd_information( attlist, curparam->can_be_set(), false, false, curparam->can_be_perturbed() );
		}
	}

	AttributeList subtag_attributes;

	subtag_attributes + XMLSchemaAttribute::required_attribute( "helix_index", xsct_positive_integer, "Numerical index for this particular helix" );


	for ( core::Size i(1); static_cast< BPC_Parameters >(i) < BPC_end_of_list; ++i ) { //Process all defaults that can be set:
		core::conformation::parametric::ParameterCOP curparam( default_calculator->parameter_cop(i) );
		if ( !curparam->global_for_parameters_set() && ( curparam->can_be_set() || curparam->can_be_copied() || curparam->can_be_perturbed() ) ) {
			curparam->provide_xsd_information( subtag_attributes, curparam->can_be_set(), curparam->can_be_copied(), false, curparam->can_be_perturbed() );
		}
	}

	utility::tag::XMLSchemaSimpleSubelementList ssl;
	ssl.add_simple_subelement( "Helix", subtag_attributes, "Tags describing the perturbation of individual helices in the bundle.");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "Perturb helical bundles by direct manipulation of their bundle parameters", attlist, ssl );
}

std::string PerturbBundleCreator::keyname() const {
	return PerturbBundle::mover_name();
}

protocols::moves::MoverOP
PerturbBundleCreator::create_mover() const {
	return protocols::moves::MoverOP( new PerturbBundle );
}

void PerturbBundleCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PerturbBundle::provide_xml_schema( xsd );
}



} //namespace helical_bundle
} //namespace protocols
