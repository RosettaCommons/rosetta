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

static THREAD_LOCAL basic::Tracer TR( "core.pack.task.operation.ResLvlTaskOperations", basic::t_info );

RestrictToRepackingRLT::~RestrictToRepackingRLT() {}

ResLvlTaskOperationOP
RestrictToRepackingRLT::clone() const { return ResLvlTaskOperationOP( new RestrictToRepackingRLT( *this ) ); }

void RestrictToRepackingRLT::apply( ResidueLevelTask & rlt ) const
{
	rlt.restrict_to_repacking();
}

std::string RestrictToRepackingRLT::keyname() { return "RestrictToRepackingRLT"; }

void RestrictToRepackingRLT::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_lvl_task_op_schema_empty( xsd, keyname() );
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

RestrictAbsentCanonicalAASRLT::~RestrictAbsentCanonicalAASRLT() {}


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
	for ( std::string::const_iterator it( aastring.begin() ), end( aastring.end() );
			it != end; ++it ) {
		if ( oneletter_code_specifies_aa( *it ) ) {
			canonical_aas_to_keep_[ aa_from_oneletter_code( *it ) ] = true;
		} else {
			std::ostringstream os;
			os << "aa letter " << *it << " does not not correspond to a canonical AA";
			utility_exit_with_message( os.str() );
		}
	}
}

void RestrictAbsentCanonicalAASRLT::parse_tag( TagCOP tag )
{
	runtime_assert( tag != 0 );
	if ( tag->hasOption("aas") ) aas_to_keep( tag->getOption<std::string>("aas") );
	else utility_exit_with_message("no aas tag option by which restrict absent canonical aas.");
}

std::string RestrictAbsentCanonicalAASRLT::keyname() { return "RestrictAbsentCanonicalAASRLT"; }

void RestrictAbsentCanonicalAASRLT::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes;
	attributes + utility::tag::XMLSchemaAttribute::required_attribute( "aas", utility::tag::xs_string );
	res_lvl_task_op_schema_w_attributes( xsd, keyname(), attributes );
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

DisallowIfNonnativeRLT::~DisallowIfNonnativeRLT(){}

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
	for ( core::Size ii=1; ii<=disallowed_aas_.size(); ii++ ) {
		inverted_vec.push_back( ! disallowed_aas[ii] );
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
	for ( std::string::const_iterator it( aa_string.begin() ), end( aa_string.end() );
			it != end; ++it ) {
		if ( oneletter_code_specifies_aa( *it ) ) {
			aa_vector[ aa_from_oneletter_code( *it ) ] = true;
		} else {
			std::ostringstream os;
			os << "aa letter " << *it << " does not not correspond to a canonical AA";
			utility_exit_with_message( os.str() );
		}
	}
	disallowed_aas_ = aa_vector;
	allowed_aas_ = invert_vector( disallowed_aas_ );
}

void DisallowIfNonnativeRLT::parse_tag( TagCOP tag )
{
	runtime_assert( tag != 0 );
	if ( tag->hasOption("disallow_aas") ) {
		disallow_aas( tag->getOption< std::string >( "disallow_aas" ) );
	} else utility_exit_with_message("no aas tag option by which restrict absent canonical aas.");
}

std::string DisallowIfNonnativeRLT::keyname() { return "DisallowIfNonnativeRLT"; }

void DisallowIfNonnativeRLT::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes;
	attributes + utility::tag::XMLSchemaAttribute::required_attribute( "disallow_aas", utility::tag::xs_string );
	res_lvl_task_op_schema_w_attributes( xsd, keyname(), attributes );
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
PreventRepackingRLT::~PreventRepackingRLT() {}

ResLvlTaskOperationOP
PreventRepackingRLT::clone() const { return ResLvlTaskOperationOP( new PreventRepackingRLT( *this ) ); }

void PreventRepackingRLT::apply( ResidueLevelTask & rlt ) const
{
	rlt.prevent_repacking();
}

std::string PreventRepackingRLT::keyname() { return "PreventRepackingRLT"; }

void PreventRepackingRLT::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_lvl_task_op_schema_empty( xsd, keyname());
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

AddBehaviorRLT::~AddBehaviorRLT() {}

ResLvlTaskOperationOP
AddBehaviorRLT::clone() const { return ResLvlTaskOperationOP( new AddBehaviorRLT( *this ) ); }

void AddBehaviorRLT::apply( ResidueLevelTask & rlt ) const
{
	runtime_assert( ! behavior_.empty() );
	rlt.add_behavior( behavior_ );
}

void AddBehaviorRLT::parse_tag( TagCOP tag )
{
	runtime_assert( tag != 0 );
	if ( tag->hasOption("behavior") ) behavior_ = tag->getOption<std::string>("behavior");
	else utility_exit_with_message("AddBehaviorRLT tag needs to define option \"behavior\".");
}

std::string AddBehaviorRLT::keyname() { return "AddBehaviorRLT"; }

void AddBehaviorRLT::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes;
	attributes + utility::tag::XMLSchemaAttribute::required_attribute( "behavior", utility::tag::xs_string );
	res_lvl_task_op_schema_w_attributes( xsd, keyname(), attributes );
}

ResLvlTaskOperationOP
AddBehaviorRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new AddBehaviorRLT );
}

std::string AddBehaviorRLTCreator::keyname() const { return AddBehaviorRLT::keyname(); }

void AddBehaviorRLTCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	AddBehaviorRLT::provide_xml_schema( xsd );
}

IncludeCurrentRLT::IncludeCurrentRLT() {}
IncludeCurrentRLT::~IncludeCurrentRLT() {}
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
	res_lvl_task_op_schema_empty( xsd, keyname() );
}

ResLvlTaskOperationOP
IncludeCurrentRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new IncludeCurrentRLT );
}

std::string IncludeCurrentRLTCreator::keyname() const { return IncludeCurrentRLT::keyname(); }

void IncludeCurrentRLTCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	IncludeCurrentRLT::provide_xml_schema( xsd );
}

PreserveCBetaRLT::PreserveCBetaRLT() {}
PreserveCBetaRLT::~PreserveCBetaRLT() {}
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
	res_lvl_task_op_schema_empty( xsd, keyname() );
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

ExtraChiCutoffRLT::~ExtraChiCutoffRLT() {}
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
	attributes + XMLSchemaAttribute::attribute_w_default( "extrachi_cutoff", xsct_non_negative_integer, utility::to_string( EXTRACHI_CUTOFF_LIMIT ) );
	res_lvl_task_op_schema_w_attributes( xsd, keyname(), attributes );
}

ResLvlTaskOperationOP
ExtraChiCutoffRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new ExtraChiCutoffRLT );
}

std::string ExtraChiCutoffRLTCreator::keyname() const { return ExtraChiCutoffRLT::keyname(); }

void ExtraChiCutoffRLTCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ExtraChiCutoffRLT::provide_xml_schema( xsd );
}

ExtraRotamersGenericRLT::ExtraRotamersGenericRLT() {}
ExtraRotamersGenericRLT::~ExtraRotamersGenericRLT() {}
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
	res_lvl_task_op_schema_w_attributes( xsd, keyname(), attributes );
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
