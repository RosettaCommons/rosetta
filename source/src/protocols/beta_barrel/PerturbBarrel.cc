// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/PerturbBarrel.cc
/// @brief  Perturbs a beta-barrel by altering the barrel parameters.
/// @details This mover calls the PerturbBarrelStrand mover.
/// @author Andy Watkins

// Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/beta_barrel/PerturbBarrel.hh>
#include <protocols/beta_barrel/PerturbBarrelCreator.hh>
#include <utility/tag/Tag.hh>

#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <core/types.hh>
#include <protocols/beta_barrel/BarrelParametrizationCalculator.hh>
#include <protocols/beta_barrel/PerturbBarrelStrand.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <utility/pointer/memory.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// STL headers:
#include <sstream>

#include <core/conformation/Conformation.hh>


namespace protocols {
namespace beta_barrel {

static basic::Tracer TR("protocols.beta_barrel.PerturbBarrel");


/// @brief Constructor for PerturbBarrel mover.
PerturbBarrel::PerturbBarrel():
	Mover("PerturbBarrel"),
	default_calculator_( new BarrelParametrizationCalculator(false) ),
	individual_strand_calculators_(),
	barrelparametersset_index_(1)
{}


/// @brief Copy constructor for PerturbBarrel mover.
PerturbBarrel::PerturbBarrel( PerturbBarrel const & src ):
	protocols::moves::Mover( src ),
	default_calculator_( utility::pointer::static_pointer_cast< BarrelParametrizationCalculator >( src.default_calculator_->clone() ) ),
	individual_strand_calculators_( /*Cloned below*/),
	barrelparametersset_index_(src.barrelparametersset_index_)
{
	// Deep-cloning the individual strand calculators:
	for ( core::Size i(1), imax(src.individual_strand_calculators_.size()); i<=imax; ++i ) {
		individual_strand_calculators_.push_back( std::pair<core::Size, BarrelParametrizationCalculatorOP>( src.individual_strand_calculators_[i].first, utility::pointer::static_pointer_cast< BarrelParametrizationCalculator >( src.individual_strand_calculators_[i].second->clone() ) ) );
	}
}


/// @brief Destructor for PerturbBarrel mover.
PerturbBarrel::~PerturbBarrel() = default;

/// @brief Clone operator to create a pointer to a fresh PerturbBarrel object that copies this one.
protocols::moves::MoverOP PerturbBarrel::clone() const {
	return utility::pointer::make_shared< PerturbBarrel >( *this );
}


/// @brief Fresh_instance operator to create a pointer to a fresh PerturbBarrel object that does NOT copy this one.
protocols::moves::MoverOP PerturbBarrel::fresh_instance() const {
	return utility::pointer::make_shared< PerturbBarrel >();
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void PerturbBarrel::apply( core::pose::Pose & pose )
{
	bool failed=false;

	if ( TR.visible() ) TR << "Copying pose." << std::endl;
	core::pose::Pose pose_copy(pose);

	if ( TR.visible() ) TR << "Finding BarrelParametersSet object in pose." << std::endl;
	BarrelParametersSetOP params_set(nullptr);
	core::Size params_set_index(0);
	core::Size n_encountered(0);
	bool breaknow(false);
	for ( core::Size i=1, imax=pose_copy.conformation().n_parameters_sets(); i<=imax; ++i ) {
		BarrelParametersSetOP cur_set( utility::pointer::dynamic_pointer_cast< BarrelParametersSet >( pose_copy.conformation().parameters_set(i) ) );
		if ( cur_set != nullptr ) {
			++n_encountered; //Increment the number of parameterssets encountered
			if ( n_encountered==barrelparametersset_index() ) {
				params_set_index=i;
				params_set=cur_set; //Assign the current owning pointer to be the params_set owning pointer.
				breaknow=true;
			}
		}
		if ( breaknow ) break;
	}
	runtime_assert_string_msg(params_set_index!=0 && params_set != nullptr, "In protocols::beta_barrel::PerturbBarrel::apply() function: BarrelParametersSet object with given index not found in pose!");

	write_report( params_set, true); //Write a pre-perturbation report summarizing the initial barrel parameter values.

	if ( TR.visible() ) TR << "Perturbing barrel parameter values." << std::endl;
	if ( !perturb_values( params_set ) ) {
		if ( TR.visible() ) TR << "Perturbation failed -- nonsensical values resulted.  Returning input pose." << std::endl;
		failed=true;
	}

	if ( !failed ) {
		if ( TR.visible() ) TR << "Rebuilding the beta-barrel conformation using the barrel parameters stored in the pose." << std::endl;
		rebuild_conformation(pose_copy, params_set, params_set_index, failed);
	}

	if ( !failed ) {
		if ( TR.visible() ) TR << "Perturbation successful.  Copying result to the input pose." << std::endl;
		pose = pose_copy;
		write_report( utility::pointer::dynamic_pointer_cast< BarrelParametersSet >( pose.conformation().parameters_set(params_set_index) ), false); //Write a post-perturbation report summarizing the final barrel parameter values.
	} else {
		if ( TR.visible() ) TR << "The current attempt generated barrel parameters that did not permit sensible geometry.  Returning input pose." << std::endl;
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
PerturbBarrel::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data_map*/
) {

	if ( tag->getName() != "PerturbBarrel" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible -- the tag name does not match the mover name.");
	}

	if ( TR.visible() ) TR << "Parsing options for PerturbBarrel (\"" << tag->getOption<std::string>("name" ,"") << "\") mover." << std::endl;

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
	for ( core::Size i(1); i<=static_cast<core::Size>( BBPC_last_parameter_to_be_sampled ); ++i ) {
		core::conformation::parametric::RealValuedParameterOP curparam( default_calculator_->real_parameter( i ) );
		if ( curparam == nullptr ) continue;
		if ( curparam->can_be_set() || curparam->can_be_perturbed() ) {
			curparam->parse_setting( tag, curparam->can_be_set(), false, curparam->can_be_perturbed(), false );
		}
	}

	//Parse options for specific strands:
	reset_strands(); //Make sure we have no defined strands.
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	for ( auto const & branch_tag : branch_tags ) {
		runtime_assert_string_msg( branch_tag->getName() == "Strand",
			"Error in PerturbBarrel::parse_my_tag(): Sub-tags for the PerturbBarrel mover must be \"Strand\".  Could not parse \"" + branch_tag->getName() + "\"." );
		core::Size strand_index( branch_tag->getOption<core::Size>("strand_index", 0) );
		runtime_assert_string_msg(strand_index>0, "In protocols::beta_barrel::PerturbBarrel::parse_my_tag() function: a strand was added, but its index was set to 0.  This is not allowed." );
		BarrelParametrizationCalculatorOP this_calculator( add_strand(strand_index) ); //Add and initialize this strand.  By default, no degrees of freedom may be perturbed.  Store the current strand index in this_strand.
		for ( core::Size i(1); i<=static_cast<core::Size>( BBPC_last_parameter_to_be_sampled ); ++i ) {
			core::conformation::parametric::RealValuedParameterOP curparam( this_calculator->real_parameter( i ) );
			if ( curparam == nullptr ) continue;
			if ( curparam->can_be_set() || curparam->can_be_perturbed() ) {
				curparam->parse_setting( branch_tag, curparam->can_be_set(), false, curparam->can_be_perturbed(), false );
			}
		}
	}

} //parse_my_tag

////////////////////////////////////////////////////////////////////////////////
//          PUBLIC FUNCTIONS                                                  //
////////////////////////////////////////////////////////////////////////////////

/// @brief Clear the list of strands.
/// @details Clears the individual_strand_calculators_ list.
void
PerturbBarrel::reset_strands() {
	individual_strand_calculators_.clear();
}

/// @brief Add options for a new strand
///
BarrelParametrizationCalculatorOP
PerturbBarrel::add_strand( core::Size const strand_index ) {
	debug_assert( default_calculator_ != nullptr ); //Shouldn't be possible.
	runtime_assert_string_msg(strand_not_defined( strand_index ), "In protocols::beta_barrel::PerturbBarrel::add_strand() function: could not add strand.  Strand defined multiple times!" );
	individual_strand_calculators_.push_back( std::make_pair( strand_index, utility::pointer::static_pointer_cast< BarrelParametrizationCalculator >( default_calculator_->clone() ) ) );

	return individual_strand_calculators_[individual_strand_calculators_.size()].second;
}

/// @brief Access the calculator for a given strand (const access).
BarrelParametrizationCalculatorCOP
PerturbBarrel::individual_strand_calculator_cop(
	core::Size const strand_calculator_index
) const {
	debug_assert( strand_calculator_index > 0 );
	debug_assert( strand_calculator_index <= individual_strand_calculators_.size() );
	for ( core::Size i(1), imax(individual_strand_calculators_.size()); i<=imax; ++i ) {
		if ( individual_strand_calculators_[i].first == strand_calculator_index ) {
			return individual_strand_calculators_[i].second;
		}
	}
	utility_exit_with_message("Error in PerturbBarrel::individual_strand_calculator_cop(): Strand with given index not found!");
	return BarrelParametrizationCalculatorCOP( nullptr ); //To keep compiler happy.
}

/// @brief Access the calculator for a given strand (nonconst access).
BarrelParametrizationCalculatorOP
PerturbBarrel::individual_strand_calculator(
	core::Size const strand_calculator_index
) {
	debug_assert( strand_calculator_index > 0 );
	debug_assert( strand_calculator_index <= individual_strand_calculators_.size() );
	for ( core::Size i(1), imax(individual_strand_calculators_.size()); i<=imax; ++i ) {
		if ( individual_strand_calculators_[i].first == strand_calculator_index ) {
			return individual_strand_calculators_[i].second;
		}
	}
	utility_exit_with_message("Error in PerturbBarrel::individual_strand_calculator(): Strand with given index not found!");
	return BarrelParametrizationCalculatorOP( nullptr ); //To keep compiler happy.
}

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////


/// @brief Confirms that a strand has not yet been defined.  Returns "true" if the strand
/// has NOT been defined, false otherwise.
bool
PerturbBarrel::strand_not_defined(
	core::Size const strand_index
) const {
	for ( core::Size i(1), imax(individual_strand_calculators_.size()); i<=imax; ++i ) {
		if ( individual_strand_calculators_[i].first == strand_index ) return false;
	}
	return true;
}

/// @brief Get a calculator for a particular strand.
/// @details Returns nullptr if no calculator for this strand has been defined.
BarrelParametrizationCalculatorCOP
PerturbBarrel::get_calculator_for_strand(
	core::Size const strand_index
) const {
	for ( core::Size i(1), imax(individual_strand_calculators_.size()); i<=imax; ++i ) {
		if ( individual_strand_calculators_[i].first == strand_index ) return individual_strand_calculators_[i].second;
	}
	return nullptr;
}

/// @brief Perturb the beta-barrel parameter values in the pose, subject to the options already set.
/// @details Called by the apply() function.  Returns true for success, false for failure.
bool PerturbBarrel::perturb_values( BarrelParametersSetOP params_set) const {

	core::Size const n_strands_in_set( params_set->n_parameters() );
	runtime_assert_string_msg( n_strands_in_set > 0,
		"In protocols::beta_barrel::PerturbBarrel::perturb_values() function: the number of strands in the pose is 0.  Unable to proceed -- nothing to perturb." );

	if ( TR.Debug.visible() ) TR.Debug << "params_set->n_parameters(): " << n_strands_in_set << std::endl;

	for ( core::Size istrand=1; istrand<=n_strands_in_set; ++istrand ) {
		BarrelParametersOP params( utility::pointer::dynamic_pointer_cast< BarrelParameters >( params_set->parameters(istrand) ) );
		runtime_assert_string_msg( params != nullptr,
			"In protocols::beta_barrel::PerturbBarrel::perturb_values() function: unable to get an owning pointer for the BarrelParameters object." );

		protocols::beta_barrel::BarrelParametrizationCalculatorCOP curcalculator( get_calculator_for_strand(istrand) ); //Will be nullptr if strand not defined.
		if ( curcalculator == nullptr ) curcalculator = default_calculator_;
		for ( core::Size iparam(1); iparam<=static_cast<core::Size>(BBPC_last_parameter_to_be_sampled); ++iparam ) {
			core::conformation::parametric::RealValuedParameterCOP curparam_calculator( curcalculator->real_parameter_cop( iparam ) );
			if ( curparam_calculator == nullptr ) continue;
			if ( !curparam_calculator->can_be_perturbed() && !curparam_calculator->can_be_set() ) continue;
			if ( !curparam_calculator->perturbation_set() && !curparam_calculator->value_was_set() ) continue;
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
			} else if ( curparam_calculator->can_be_set() && curparam_calculator->value_was_set() ) {
				curparam->set_value( curparam_calculator->value(), true );
			}
		}
	}

	return true;
}

/// @brief Rebuild the beta-barrel conformation using the barrel parameter values in the pose.
///
void PerturbBarrel::rebuild_conformation(
	core::pose::Pose &pose,
	BarrelParametersSetOP /*params_set*/,
	core::Size const params_set_index,
	bool &failed
) const {

	core::Size const n_params( pose.conformation().parameters_set( params_set_index )->n_parameters() ); //Get the number of parameters objects (number of strands that we're rebuilding).

	for ( core::Size istrand=1; istrand<=n_params; ++istrand ) {

		PerturbBarrelStrand pertstrand; //Construct the mover to perturb this strand.

		pertstrand.set_parameters_set_index(params_set_index);
		pertstrand.set_parameters_index(istrand);

		pertstrand.apply( pose );

		failed = pertstrand.last_apply_failed();
		if ( failed ) break;

	}

	return;
}

/// @brief Write out the perturbed barrel parameters.
/// @details The "before" parameter determines whether this is a pre-perturbation
/// report or a post-perturbation report.
void PerturbBarrel::write_report(
	BarrelParametersSetOP params_set,
	bool const before
) const {
	if ( !TR.visible() ) return; //Do nothing if the tracer isn't visible.

	if ( before ) {
		TR << "*** Barrel parameter values prior to perturbation ***" << std::endl;
	} else {
		TR << "*** Barrel parameter values after perturbation ***" << std::endl;
	}

	core::Size const nparams(params_set->n_parameters());
	if ( nparams==0 ) {
		TR << "No Parameters objects in the ParametersSet object!" << std::endl;
	}

	for ( core::Size iparams=1; iparams<=nparams; ++iparams ) { //Loop through the params objects
		BarrelParametersOP params( utility::pointer::dynamic_pointer_cast<BarrelParameters>(params_set->parameters(iparams)) );
		if ( params == nullptr ) continue;
		TR << "Strand " << iparams << ":" << std::endl;
		std::stringstream remark;
		params->get_pdb_remark( remark );
		TR << remark.str() << std::endl;
	}

	TR << "*** END SUMMARY ***" << std::endl;
	TR.flush();

	return;
}

std::string PerturbBarrel::get_name() const {
	return mover_name();
}

std::string PerturbBarrel::mover_name() {
	return "PerturbBarrel";
}

void PerturbBarrel::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	BarrelParametrizationCalculatorOP default_calculator( utility::pointer::make_shared< BarrelParametrizationCalculator >() );

	XMLSchemaRestriction pert_type;
	pert_type.name( "barrel_pert_type" );
	pert_type.base_type( xs_string );
	pert_type.add_restriction( xsr_enumeration, "gaussian" );
	pert_type.add_restriction( xsr_enumeration, "uniform" );
	pert_type.add_restriction( xsr_enumeration, "" ); // it can be the empty string
	xsd.add_top_level_element( pert_type );

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "use_degrees", xsct_rosetta_bool, "Interpret user-supplied angles as degrees rather than radians", "false" )
		+ XMLSchemaAttribute( "default_perturbation_type", "barrel_pert_type", "Default type for perturbations to the barrel parameters, either uniform or gaussian" );

	for ( core::Size i(1); static_cast< BBPC_Parameters >(i) < BBPC_end_of_list; ++i ) { //Process all defaults that can be set:
		core::conformation::parametric::ParameterCOP curparam( default_calculator->parameter_cop(i) );
		if ( curparam->can_be_set() || curparam->can_be_perturbed() ) {
			curparam->provide_xsd_information( attlist, curparam->can_be_set(), false, false, curparam->can_be_perturbed() );
		}
	}

	AttributeList subtag_attributes;

	subtag_attributes + XMLSchemaAttribute::required_attribute( "strand_index", xsct_positive_integer, "Numerical index for this particular strand" );

	for ( core::Size i(1); static_cast< BBPC_Parameters >(i) < BBPC_end_of_list; ++i ) { //Process all defaults that can be set:
		core::conformation::parametric::ParameterCOP curparam( default_calculator->parameter_cop(i) );
		if ( curparam->can_be_set() || curparam->can_be_perturbed() ) {
			curparam->provide_xsd_information( subtag_attributes, curparam->can_be_set(), false, false, curparam->can_be_perturbed() );
		}
	}

	utility::tag::XMLSchemaSimpleSubelementList ssl;
	ssl.add_simple_subelement( "Strand", subtag_attributes, "Tags describing the perturbation of individual strands in the barrel.");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "Perturb beta-barrels by direct manipulation of their barrel parameters", attlist, ssl );
}

/// @brief Provide the citation.
void
PerturbBarrel::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"PerturbBarrel", basic::citation_manager::CitedModuleType::Mover,
		"Andy Watkins",
		"Department of Biochemistry, Stanford University",
		"watkina6@stanford.edu"
		)
	);
}

std::string PerturbBarrelCreator::keyname() const {
	return PerturbBarrel::mover_name();
}

protocols::moves::MoverOP
PerturbBarrelCreator::create_mover() const {
	return utility::pointer::make_shared< PerturbBarrel >();
}

void PerturbBarrelCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PerturbBarrel::provide_xml_schema( xsd );
}



} //namespace beta_barrel
} //namespace protocols
