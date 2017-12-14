// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ResLvlTaskOperations.cc
/// @brief
/// @author ashworth

// Unit Headers
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/ResLvlTaskOperationCreators.hh>

#include <core/pack/task/operation/task_op_schemas.hh>

#include <core/pack/task/PackerTask.hh>

#include <core/chemical/AA.hh>

#include <basic/Tracer.hh>

#include <utility>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <boost/foreach.hpp>

#include <sstream>


namespace core {
namespace pack {
namespace task {
namespace operation {

static basic::Tracer TR( "core.pack.task.operation.ResLvlTaskOperations", basic::t_info );

RestrictToRepackingRLT::~RestrictToRepackingRLT() = default;

ResLvlTaskOperationOP
RestrictToRepackingRLT::clone() const { return ResLvlTaskOperationOP( new RestrictToRepackingRLT( *this ) ); }

void RestrictToRepackingRLT::apply( ResidueLevelTask & rlt ) const
{
	rlt.restrict_to_repacking();
}

std::string RestrictToRepackingRLT::keyname() { return "RestrictToRepackingRLT"; }

void RestrictToRepackingRLT::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_lvl_task_op_schema_empty(
		xsd, keyname(),
		"Turn off design on the positions selected by the accompanying ResFilter.");
}

ResLvlTaskOperationOP
RestrictToRepackingRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new RestrictToRepackingRLT );
}

std::string RestrictToRepackingRLTCreator::keyname() const { return RestrictToRepackingRLT::keyname(); }

void RestrictToRepackingRLTCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	RestrictToRepackingRLT::provide_xml_schema( xsd );
}


RestrictAbsentCanonicalAASRLT::RestrictAbsentCanonicalAASRLT()
: canonical_aas_to_keep_( chemical::num_canonical_aas, false )
{}

RestrictAbsentCanonicalAASRLT::~RestrictAbsentCanonicalAASRLT() = default;


ResLvlTaskOperationOP
RestrictAbsentCanonicalAASRLT::clone() const { return ResLvlTaskOperationOP( new RestrictAbsentCanonicalAASRLT( *this ) ); }

void RestrictAbsentCanonicalAASRLT::apply( ResidueLevelTask & rlt ) const
{
	rlt.restrict_absent_canonical_aas( canonical_aas_to_keep_ );
}

// if an amino acid is not present (false) in the boolean vector, then do not allow it at this position.  The boolean vector is a 20-length vector in alphabetical order by one-letter code.
void RestrictAbsentCanonicalAASRLT::aas_to_keep( utility::vector1< bool > const & aa_flags )
{
	runtime_assert( aa_flags.size() == chemical::num_canonical_aas );
	canonical_aas_to_keep_ = aa_flags;
}

void RestrictAbsentCanonicalAASRLT::aas_to_keep( std::string const & aastring )
{
	using namespace chemical;
	runtime_assert( canonical_aas_to_keep_.size() == num_canonical_aas );
	for ( char const code : aastring ) {
		if ( oneletter_code_specifies_aa( code ) ) {
			canonical_aas_to_keep_[ aa_from_oneletter_code( code ) ] = true;
		} else {
			std::ostringstream os;
			os << "aa letter " << code << " does not not correspond to a canonical AA";
			utility_exit_with_message( os.str() );
		}
	}
}

void RestrictAbsentCanonicalAASRLT::parse_tag( TagCOP tag )
{
	runtime_assert( tag != nullptr );
	if ( tag->hasOption("aas") ) aas_to_keep( tag->getOption<std::string>("aas") );
	else utility_exit_with_message("no aas tag option by which restrict absent canonical aas.");
}

std::string RestrictAbsentCanonicalAASRLT::keyname() { return "RestrictAbsentCanonicalAASRLT"; }

void RestrictAbsentCanonicalAASRLT::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes;
	attributes + utility::tag::XMLSchemaAttribute::required_attribute(
		"aas", utility::tag::xs_string,
		"list of one letter codes of permitted amino acids, with no separator. "
		"(e.g. aas=HYFW for only aromatic amino acids.)" );

	res_lvl_task_op_schema_w_attributes(
		xsd, keyname(), attributes,
		"Do not allow design to amino acid identities that are not listed (i.e. permit only those listed) "
		"at the positions selected by the accompanying ResFilter.");
}

void RestrictAbsentCanonicalAASRLTCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	RestrictAbsentCanonicalAASRLT::provide_xml_schema( xsd );
}

std::string RestrictAbsentCanonicalAASRLTCreator::keyname() const { return RestrictAbsentCanonicalAASRLT::keyname(); }

ResLvlTaskOperationOP
RestrictAbsentCanonicalAASRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new RestrictAbsentCanonicalAASRLT );
}

//BEGIN DisallowIfNonnativeRLTRLT
DisallowIfNonnativeRLT::DisallowIfNonnativeRLT():
	disallowed_aas_ ( chemical::num_canonical_aas, false ),
	allowed_aas_( invert_vector( disallowed_aas_ ) )
{}
DisallowIfNonnativeRLT::DisallowIfNonnativeRLT( utility::vector1< bool > disallowed_aas):
	disallowed_aas_( disallowed_aas ),
	allowed_aas_( invert_vector(disallowed_aas) )
{}

DisallowIfNonnativeRLT::~DisallowIfNonnativeRLT()= default;

ResLvlTaskOperationOP DisallowIfNonnativeRLT::clone() const
{
	return ResLvlTaskOperationOP( new DisallowIfNonnativeRLT( *this ) );
}

void DisallowIfNonnativeRLT::clear(){
	allowed_aas_.clear();
	disallowed_aas_.clear();
}
//private function to invert disallowed aas into allowed aas
utility::vector1< bool >
DisallowIfNonnativeRLT::invert_vector( utility::vector1< bool > disallowed_aas){
	utility::vector1< bool > inverted_vec;
	for ( bool const pos : disallowed_aas ) {
		inverted_vec.push_back( ! pos );
	}
	return inverted_vec;
}

void DisallowIfNonnativeRLT::apply( ResidueLevelTask & rlt ) const
{
	runtime_assert( allowed_aas_.size() == chemical::num_canonical_aas );
	rlt.restrict_nonnative_canonical_aas( allowed_aas_ );
}

void DisallowIfNonnativeRLT::disallow_aas( utility::vector1< bool > const & cannonical_disallowed ){
	runtime_assert( cannonical_disallowed.size() == chemical::num_canonical_aas );
	disallowed_aas_ = cannonical_disallowed;
	allowed_aas_ = invert_vector( disallowed_aas_ );
}
void DisallowIfNonnativeRLT::disallow_aas( std::string const & aa_string ){
	using namespace chemical;
	utility::vector1< bool > aa_vector ( chemical::num_canonical_aas, false );
	for ( char const code : aa_string ) {
		if ( oneletter_code_specifies_aa( code ) ) {
			aa_vector[ aa_from_oneletter_code( code ) ] = true;
		} else {
			std::ostringstream os;
			os << "aa letter " << code << " does not not correspond to a canonical AA";
			utility_exit_with_message( os.str() );
		}
	}
	disallowed_aas_ = aa_vector;
	allowed_aas_ = invert_vector( disallowed_aas_ );
}

void DisallowIfNonnativeRLT::parse_tag( TagCOP tag )
{
	runtime_assert( tag != nullptr );
	if ( tag->hasOption("disallow_aas") ) {
		disallow_aas( tag->getOption< std::string >( "disallow_aas" ) );
	} else utility_exit_with_message("no aas tag option by which restrict absent canonical aas.");
}

std::string DisallowIfNonnativeRLT::keyname() { return "DisallowIfNonnativeRLT"; }

void DisallowIfNonnativeRLT::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes;
	attributes + utility::tag::XMLSchemaAttribute::required_attribute(
		"disallow_aas", utility::tag::xs_string ,
		"takes a string of one letter amino acid codes, no separation needed. "
		"For example disallow_aas=GCP would prevent Gly, Cys, and Pro from being designed "
		"unless they were the native amino acid at a position."
		"This task is useful when you are designing in a region that has Gly and Pro and "
		"you do not want to include them at other positions that aren't already Gly or Pro.");

	res_lvl_task_op_schema_w_attributes(
		xsd, keyname(), attributes,
		"Restrict design to not include a residue as an possibility in the task at a position unless it is the starting residue. "
		"If resnum is left as 0, the restriction will apply throughout the pose.");
}

ResLvlTaskOperationOP
DisallowIfNonnativeRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new DisallowIfNonnativeRLT );
}

std::string DisallowIfNonnativeRLTCreator::keyname() const { return DisallowIfNonnativeRLT::keyname(); }

void DisallowIfNonnativeRLTCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	DisallowIfNonnativeRLT::provide_xml_schema( xsd );
}

//Begin PreventRepackingRLT
PreventRepackingRLT::~PreventRepackingRLT() = default;

ResLvlTaskOperationOP
PreventRepackingRLT::clone() const { return ResLvlTaskOperationOP( new PreventRepackingRLT( *this ) ); }

void PreventRepackingRLT::apply( ResidueLevelTask & rlt ) const
{
	rlt.prevent_repacking();
}

std::string PreventRepackingRLT::keyname() { return "PreventRepackingRLT"; }

void PreventRepackingRLT::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_lvl_task_op_schema_empty(
		xsd, keyname(),
		"Turn off design and repacking on the positions selected by the accompanying ResFilter.");
}

ResLvlTaskOperationOP
PreventRepackingRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new PreventRepackingRLT );
}

std::string PreventRepackingRLTCreator::keyname() const { return PreventRepackingRLT::keyname(); }

void PreventRepackingRLTCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PreventRepackingRLT::provide_xml_schema( xsd );
}

AddBehaviorRLT::AddBehaviorRLT() : parent() {}

AddBehaviorRLT::AddBehaviorRLT( std::string const & behavior )
: parent(),
	behavior_( behavior )
{}

AddBehaviorRLT::~AddBehaviorRLT() = default;

ResLvlTaskOperationOP
AddBehaviorRLT::clone() const { return ResLvlTaskOperationOP( new AddBehaviorRLT( *this ) ); }

void AddBehaviorRLT::apply( ResidueLevelTask & rlt ) const
{
	runtime_assert( ! behavior_.empty() );
	rlt.add_behavior( behavior_ );
}

void AddBehaviorRLT::parse_tag( TagCOP tag )
{
	runtime_assert( tag != nullptr );
	if ( tag->hasOption("behavior") ) behavior_ = tag->getOption<std::string>("behavior");
	else utility_exit_with_message("AddBehaviorRLT tag needs to define option \"behavior\".");
}

std::string AddBehaviorRLT::keyname() { return "AddBehaviorRLT"; }

void AddBehaviorRLT::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes;
	attributes + utility::tag::XMLSchemaAttribute::required_attribute(
		"behavior", utility::tag::xs_string,
		"Behavior string. These are protocol-specific. "
		"Consult the protocol documentation for if it responds to behavior strings.");
	res_lvl_task_op_schema_w_attributes(
		xsd, keyname(), attributes,
		"Add the given \"behavior\" to the positions selected by the accompanying ResFilter.");
}

ResLvlTaskOperationOP
AddBehaviorRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new AddBehaviorRLT );
}

std::string AddBehaviorRLTCreator::keyname() const { return AddBehaviorRLT::keyname(); }

void AddBehaviorRLTCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	AddBehaviorRLT::provide_xml_schema( xsd );
}

IncludeCurrentRLT::IncludeCurrentRLT() = default;
IncludeCurrentRLT::~IncludeCurrentRLT() = default;
ResLvlTaskOperationOP IncludeCurrentRLT::clone() const
{
	return ResLvlTaskOperationOP( new IncludeCurrentRLT );
}

void IncludeCurrentRLT::apply( ResidueLevelTask & rlt ) const
{
	rlt.or_include_current( true );
}

std::string IncludeCurrentRLT::keyname() { return "IncludeCurrentRLT"; }

void IncludeCurrentRLT::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_lvl_task_op_schema_empty(
		xsd, keyname(),
		"Includes current rotamers (eg - from input pdb) in the rotamer set. ");
}

ResLvlTaskOperationOP
IncludeCurrentRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new IncludeCurrentRLT );
}

std::string IncludeCurrentRLTCreator::keyname() const { return IncludeCurrentRLT::keyname(); }

void IncludeCurrentRLTCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	IncludeCurrentRLT::provide_xml_schema( xsd );
}

PreserveCBetaRLT::PreserveCBetaRLT() = default;
PreserveCBetaRLT::~PreserveCBetaRLT() = default;
ResLvlTaskOperationOP PreserveCBetaRLT::clone() const
{
	return ResLvlTaskOperationOP( new PreserveCBetaRLT );
}
void PreserveCBetaRLT::apply( ResidueLevelTask & rlt ) const
{
	rlt.or_preserve_c_beta( true );
}

std::string PreserveCBetaRLT::keyname() { return "PreserveCBetaRLT"; }

void PreserveCBetaRLT::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_lvl_task_op_schema_empty(
		xsd, keyname(),
		"preserves c-beta during rotamer building for all residues. "
		"Under development and untested. Use at your own risk.");
}

ResLvlTaskOperationOP
PreserveCBetaRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new PreserveCBetaRLT );
}

std::string PreserveCBetaRLTCreator::keyname() const { return PreserveCBetaRLT::keyname(); }

void PreserveCBetaRLTCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PreserveCBetaRLT::provide_xml_schema( xsd );
}

ExtraChiCutoffRLT::ExtraChiCutoffRLT() :
	extrachi_cutoff_( EXTRACHI_CUTOFF_LIMIT )
{}

ExtraChiCutoffRLT::~ExtraChiCutoffRLT() = default;
ResLvlTaskOperationOP ExtraChiCutoffRLT::clone() const
{
	return ResLvlTaskOperationOP( new ExtraChiCutoffRLT( *this ));
}
void ExtraChiCutoffRLT::apply( ResidueLevelTask & rlt ) const
{
	rlt.and_extrachi_cutoff( extrachi_cutoff_ );
}

void ExtraChiCutoffRLT::parse_tag( TagCOP tag )
{
	extrachi_cutoff_ = tag->getOption< core::Size >( "extrachi_cutoff", EXTRACHI_CUTOFF_LIMIT );
}

std::string ExtraChiCutoffRLT::keyname() { return "ExtraChiCutoffRLT"; }

void ExtraChiCutoffRLT::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes + XMLSchemaAttribute::attribute_w_default(
		"extrachi_cutoff", xsct_non_negative_integer,
		"lower extrachi_cutoff to given value; do nothing if not a decrease",
		utility::to_string( EXTRACHI_CUTOFF_LIMIT  ) );

	res_lvl_task_op_schema_w_attributes(
		xsd, keyname(), attributes,
		"Move only toward a lower cutoff for #neighbors w/i 10A that qualify a residue to be considered buried.");
}

ResLvlTaskOperationOP
ExtraChiCutoffRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new ExtraChiCutoffRLT );
}

std::string ExtraChiCutoffRLTCreator::keyname() const { return ExtraChiCutoffRLT::keyname(); }

void ExtraChiCutoffRLTCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ExtraChiCutoffRLT::provide_xml_schema( xsd );
}

ExtraRotamersGenericRLT::ExtraRotamersGenericRLT() = default;
ExtraRotamersGenericRLT::~ExtraRotamersGenericRLT() = default;
ResLvlTaskOperationOP ExtraRotamersGenericRLT::clone() const
{
	return ResLvlTaskOperationOP( new ExtraRotamersGenericRLT( *this ));
}
void ExtraRotamersGenericRLT::apply( ResidueLevelTask & rlt ) const
{
	set_rotamer_sampling_data_for_RLT( sampling_data_, rlt );
}
void ExtraRotamersGenericRLT::parse_tag( TagCOP tag )
{
	parse_rotamer_sampling_data( tag, sampling_data_ );
}

void ExtraRotamersGenericRLT::ex1( bool value ) {
	sampling_data_.ex1_ = value;
}
void ExtraRotamersGenericRLT::ex2( bool value ) {
	sampling_data_.ex2_ = value;
}
void ExtraRotamersGenericRLT::ex3( bool value ) {
	sampling_data_.ex3_ = value;
}
void ExtraRotamersGenericRLT::ex4( bool value ) {
	sampling_data_.ex4_ = value;
}
void ExtraRotamersGenericRLT::ex1aro( bool value ) {
	sampling_data_.ex1aro_ = value;
}
void ExtraRotamersGenericRLT::ex2aro( bool value ) {
	sampling_data_.ex2aro_ = value;
}
void ExtraRotamersGenericRLT::ex1aro_exposed( bool value ) {
	sampling_data_.ex1aro_exposed_ = value;
}
void ExtraRotamersGenericRLT::ex2aro_exposed( bool value ) {
	sampling_data_.ex2aro_exposed_ = value;
}
void ExtraRotamersGenericRLT::ex1_sample_level( ExtraRotSample value ) {
	sampling_data_.ex1_sample_level_ = value;
}
void ExtraRotamersGenericRLT::ex2_sample_level( ExtraRotSample value ) {
	sampling_data_.ex2_sample_level_ = value;
}
void ExtraRotamersGenericRLT::ex3_sample_level( ExtraRotSample value ) {
	sampling_data_.ex3_sample_level_ = value;
}
void ExtraRotamersGenericRLT::ex4_sample_level( ExtraRotSample value ) {
	sampling_data_.ex4_sample_level_ = value;
}
void ExtraRotamersGenericRLT::ex1aro_sample_level( ExtraRotSample value ) {
	sampling_data_.ex1aro_sample_level_ = value;
}
void ExtraRotamersGenericRLT::ex2aro_sample_level( ExtraRotSample value ) {
	sampling_data_.ex2aro_sample_level_ = value;
}
void ExtraRotamersGenericRLT::ex1aro_exposed_sample_level( ExtraRotSample value ) {
	sampling_data_.ex1aro_exposed_sample_level_ = value;
}
void ExtraRotamersGenericRLT::ex2aro_exposed_sample_level( ExtraRotSample value ) {
	sampling_data_.ex2aro_exposed_sample_level_ = value;
}
void ExtraRotamersGenericRLT::exdna_sample_level( ExtraRotSample value ) {
	sampling_data_.exdna_sample_level_ = value;
}
void ExtraRotamersGenericRLT::extrachi_cutoff( Size value ) {
	sampling_data_.extrachi_cutoff_ = value;
}

ExtraRotamerSamplingData const &
ExtraRotamersGenericRLT::sampling_data() const
{
	return sampling_data_;
}

std::string ExtraRotamersGenericRLT::keyname() { return "ExtraRotamersGenericRLT"; }

void ExtraRotamersGenericRLT::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes = rotamer_sampling_data_xml_schema_attributes( xsd );
	res_lvl_task_op_schema_w_attributes(
		xsd, keyname(), attributes,
		"During packing, extra rotamers can be used to increase sampling. "
		"Use this TaskOperation to specify for all residues at once what extra rotamers should be used. "
		"Note: The extrachi_cutoff is used to determine how many neighbors a residue "
		"must have before the extra rotamers are applied. For example of you want "
		"to apply extra rotamers to all residues, set extrachi_cutoff=0. ");
}

ResLvlTaskOperationOP
ExtraRotamersGenericRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new ExtraRotamersGenericRLT );
}

std::string ExtraRotamersGenericRLTCreator::keyname() const { return ExtraRotamersGenericRLT::keyname(); }

void ExtraRotamersGenericRLTCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ExtraRotamersGenericRLT::provide_xml_schema( xsd );
}



} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
