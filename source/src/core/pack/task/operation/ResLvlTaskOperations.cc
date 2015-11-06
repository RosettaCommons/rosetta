// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/ResLvlTaskOperations.cc
/// @brief
/// @author ashworth

// Unit Headers
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/ResLvlTaskOperationCreators.hh>

#include <core/pack/task/PackerTask.hh>

#include <core/chemical/AA.hh>

#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/exit.hh>
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
RestrictToRepackingRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new RestrictToRepackingRLT );
}

ResLvlTaskOperationOP
RestrictToRepackingRLT::clone() const { return ResLvlTaskOperationOP( new RestrictToRepackingRLT( *this ) ); }

void RestrictToRepackingRLT::apply( ResidueLevelTask & rlt ) const
{
	rlt.restrict_to_repacking();
}

RestrictAbsentCanonicalAASRLT::RestrictAbsentCanonicalAASRLT()
: canonical_aas_to_keep_( chemical::num_canonical_aas, false )
{}

RestrictAbsentCanonicalAASRLT::~RestrictAbsentCanonicalAASRLT() {}

ResLvlTaskOperationOP
RestrictAbsentCanonicalAASRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new RestrictAbsentCanonicalAASRLT );
}

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

ResLvlTaskOperationOP
DisallowIfNonnativeRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new DisallowIfNonnativeRLT );
}

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

//Begin PreventRepackingRLT
PreventRepackingRLT::~PreventRepackingRLT() {}

ResLvlTaskOperationOP
PreventRepackingRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new PreventRepackingRLT );
}

ResLvlTaskOperationOP
PreventRepackingRLT::clone() const { return ResLvlTaskOperationOP( new PreventRepackingRLT( *this ) ); }

void PreventRepackingRLT::apply( ResidueLevelTask & rlt ) const
{
	rlt.prevent_repacking();
}

AddBehaviorRLT::AddBehaviorRLT() : parent() {}

AddBehaviorRLT::AddBehaviorRLT( std::string const & behavior )
: parent(),
	behavior_( behavior )
{}

AddBehaviorRLT::~AddBehaviorRLT() {}

ResLvlTaskOperationOP
AddBehaviorRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new AddBehaviorRLT );
}

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

ResLvlTaskOperationOP
IncludeCurrentRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new IncludeCurrentRLT );
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

ResLvlTaskOperationOP
PreserveCBetaRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new PreserveCBetaRLT );
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

ResLvlTaskOperationOP
ExtraChiCutoffRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new ExtraChiCutoffRLT );
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


ResLvlTaskOperationOP
ExtraRotamersGenericRLTCreator::create_res_level_task_operation() const {
	return ResLvlTaskOperationOP( new ExtraRotamersGenericRLT );
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



} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
