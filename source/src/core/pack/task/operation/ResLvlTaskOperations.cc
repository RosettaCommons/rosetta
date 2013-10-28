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
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

#include <sstream>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {

RestrictToRepackingRLT::~RestrictToRepackingRLT() {}

ResLvlTaskOperationOP
RestrictToRepackingRLTCreator::create_res_level_task_operation() const {
	return new RestrictToRepackingRLT;
}

ResLvlTaskOperationOP
RestrictToRepackingRLT::clone() const { return new RestrictToRepackingRLT( *this ); }

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
	return new RestrictAbsentCanonicalAASRLT;
}

ResLvlTaskOperationOP
RestrictAbsentCanonicalAASRLT::clone() const { return new RestrictAbsentCanonicalAASRLT( *this ); }

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
	runtime_assert( tag );
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
	return new DisallowIfNonnativeRLT;
}

ResLvlTaskOperationOP DisallowIfNonnativeRLT::clone() const
{
	return new DisallowIfNonnativeRLT( *this );
}

void DisallowIfNonnativeRLT::clear(){
	allowed_aas_.clear();
	disallowed_aas_.clear();
}
//private function to invert disallowed aas into allowed aas
utility::vector1< bool >
DisallowIfNonnativeRLT::invert_vector( utility::vector1< bool > disallowed_aas){
	utility::vector1< bool > inverted_vec;
	for(core::Size ii=1; ii<=disallowed_aas_.size(); ii++ ){
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
	runtime_assert( tag );
	if ( tag->hasOption("disallow_aas") )
		disallow_aas( tag->getOption< std::string >( "disallow_aas" ) );
	else utility_exit_with_message("no aas tag option by which restrict absent canonical aas.");
}

	//Begin PreventRepackingRLT
PreventRepackingRLT::~PreventRepackingRLT() {}

ResLvlTaskOperationOP
PreventRepackingRLTCreator::create_res_level_task_operation() const {
	return new PreventRepackingRLT;
}

ResLvlTaskOperationOP
PreventRepackingRLT::clone() const { return new PreventRepackingRLT( *this ); }

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
	return new AddBehaviorRLT;
}

ResLvlTaskOperationOP
AddBehaviorRLT::clone() const { return new AddBehaviorRLT( *this ); }

void AddBehaviorRLT::apply( ResidueLevelTask & rlt ) const
{
	runtime_assert( ! behavior_.empty() );
	rlt.add_behavior( behavior_ );
}

void AddBehaviorRLT::parse_tag( TagCOP tag )
{
	runtime_assert( tag );
	if ( tag->hasOption("behavior") ) behavior_ = tag->getOption<std::string>("behavior");
	else utility_exit_with_message("AddBehaviorRLT tag needs to define option \"behavior\".");
}

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
