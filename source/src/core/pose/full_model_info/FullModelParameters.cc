// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/full_model_info/FullModelParameters.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "core.pose.full_model_info.FullModelParameters" );

namespace core {
namespace pose {
namespace full_model_info {

	//Constructor
	FullModelParameters::FullModelParameters( std::string const full_sequence ):
		full_sequence_( full_sequence )
	{
		initialize_parameters( *this );
		for ( Size n = 1; n <= full_sequence.size(); n++ ) conventional_numbering_.push_back( n );
	}

	//Constructor
	FullModelParameters::FullModelParameters( std::string const full_sequence,
																						utility::vector1< Size > const & cutpoint_open_in_full_model,
																						utility::vector1< Size > const & res_numbers_in_pose	):
		full_sequence_( full_sequence )
	{
		initialize_parameters( *this );
		set_parameter( CUTPOINT_OPEN, convert_to_parameter_values_at_res( cutpoint_open_in_full_model ) );
		utility::vector1< Size > fixed_domain_map;
		for ( Size n = 1; n <= full_sequence.size(); n++ ){
			if ( res_numbers_in_pose.has( n ) ){
				fixed_domain_map.push_back( 1 );
			} else {
				fixed_domain_map.push_back( 0 );
			}
			conventional_numbering_.push_back( static_cast<int>( n ) );
		}
		set_parameter( FIXED_DOMAIN, fixed_domain_map );
	}

	//Constructor
	FullModelParameters::FullModelParameters( pose::Pose const & pose,
																						utility::vector1< Size > const & res_list ) {

		initialize_parameters( *this );
		get_sequence_with_gaps_filled_with_n( pose, full_sequence_, conventional_numbering_, res_list );
		set_parameter( CUTPOINT_OPEN, convert_to_parameter_values_at_res( get_cutpoint_open_from_pdb_info( pose ) ) );

		// not sure what's best here -- for now setting that the pose's residues are 'fixed' within the domain map.
		utility::vector1< Size > fixed_domain_map_( full_sequence_.size(), 0 );
		for ( Size n = 1; n <= res_list.size(); n++ ) {
			Size const & res_num = res_list[ n ];
			runtime_assert( conventional_numbering_.has_value( static_cast<int>( res_num ) ) );
			fixed_domain_map_[ conventional_numbering_.index( static_cast<int>( res_num ) ) ] = 1;
		}
		set_parameter( FIXED_DOMAIN, fixed_domain_map_ );

	}

	// copy
	FullModelParameters::FullModelParameters( FullModelParameters const & src ) :
		ReferenceCount( src ),
		parameter_values_at_res_( src.parameter_values_at_res_ ),
		parameter_values_as_res_lists_( src.parameter_values_as_res_lists_ ),
		full_sequence_( src.full_sequence_ ),
		conventional_numbering_( src.conventional_numbering_ )
	{
	}

	//Destructor
	FullModelParameters::~FullModelParameters()
	{}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::map< Size, utility::vector1< Size > >
	FullModelParameters::convert_to_res_lists_by_value( utility::vector1< Size > const & parameter_values_at_res ){
		std::map< Size, utility::vector1< Size > > res_lists_by_value;
		runtime_assert( size() == parameter_values_at_res.size() );

		// basic initialization.
		utility::vector1< Size > blank_vector;
		res_lists_by_value[ 0 ] = blank_vector;
		res_lists_by_value[ 1 ] = blank_vector;

		for ( Size n = 1; n <= parameter_values_at_res.size(); n++ ){
			Size const & value = parameter_values_at_res[ n ];
			if ( res_lists_by_value.find( value ) == res_lists_by_value.end() ){
				res_lists_by_value[ value ] = blank_vector;
			}
			res_lists_by_value[ value ].push_back( n );
		}

		return res_lists_by_value;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	FullModelParameters::convert_to_parameter_values_at_res( utility::vector1< Size > const & res_list ){
		utility::vector1< Size > parameter_values_at_res( size(), 0 );
		for ( Size n = 1; n <= res_list.size(); n++ ){
			runtime_assert ( res_list[n] > 0 && res_list[n] < size() );
			parameter_values_at_res[ res_list[n] ] = 1;
		}
		return parameter_values_at_res;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size > const &
	FullModelParameters::get_res_list( FullModelParameterType const type, Size const value ) const {
		runtime_assert( parameter_values_as_res_lists_.find( type ) != parameter_values_as_res_lists_.end() );
		std::map< Size, utility::vector1< Size > > const & res_list_map = parameter_values_as_res_lists_.find( type )->second;

		runtime_assert( res_list_map.find( value ) != res_list_map.end() );
		return res_list_map.find( value )->second;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size > const &
	FullModelParameters::get_parameter( FullModelParameterType const type ) const {
		runtime_assert( parameter_values_at_res_.find( type ) != parameter_values_at_res_.end() );
		return parameter_values_at_res_.find( type )->second;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	FullModelParameters::get_sequence_with_gaps_filled_with_n( pose::Pose const & pose,
																														 std::string & sequence,
																														 utility::vector1< int > & full_numbering,
																														 utility::vector1< Size > const & res_list ) const {
		// should also be smart about not filling in n's between chains.
		// anyway. this is a quick hack for now.
		sequence = "";
		full_numbering.clear();

		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			Size const prev_res_num    = (n > 1 ) ? res_list[ n-1 ] : res_list[ n ];
			Size const current_res_num = res_list[ n ];
			for ( Size i = prev_res_num+1; i < current_res_num; i++ ) {
				sequence.push_back( 'n' );
				full_numbering.push_back( static_cast<int>( i ) );
			}
			sequence.push_back( pose.sequence()[ n-1 ] );
			full_numbering.push_back( static_cast<int>( current_res_num ) );
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	FullModelParameters::get_cutpoint_open_from_pdb_info( pose::Pose const & pose ) const {

		PDBInfoCOP pdb_info = pose.pdb_info();
		utility::vector1< Size > cutpoint_open;

		for ( Size n = 1; n < pose.total_residue(); n++ ){

			if ( pdb_info &&  (pdb_info->chain( n ) != pdb_info->chain( n+1 )) )	{
				cutpoint_open.push_back( n );
				continue;
			}

		}

		return cutpoint_open;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	FullModelParameters::chains_in_full_model() const {
		utility::vector1< Size > const & cutpoint_open_in_full_model = get_parameter( CUTPOINT_OPEN );
		utility::vector1< Size > chains;
		Size chain_number( 1 );
		for ( Size n = 1; n <= full_sequence_.size(); n++ ){
			chains.push_back( chain_number );
			if ( cutpoint_open_in_full_model.has_value( n ) ) chain_number++;
		}
		return chains;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	FullModelParameters::set_parameter( FullModelParameterType const type,
																			utility::vector1< Size > const & setting ){

		parameter_values_at_res_[ type ] = setting;
		parameter_values_as_res_lists_[ type ] = convert_to_res_lists_by_value( setting );

		keep_chain_and_cutpoint_open_matched( type );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	FullModelParameters::set_parameter_as_res_list( FullModelParameterType const type,
																									utility::vector1< Size > const & setting ){
		set_parameter( type, convert_to_parameter_values_at_res( setting ) );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	FullModelParameters::keep_chain_and_cutpoint_open_matched( FullModelParameterType const & type ){
		// kind of special -- CHAIN is a convenient thing to hold, but should match info in CUTPOINT_OPEN
		if ( type == CUTPOINT_OPEN ) set_parameter( CHAIN, chains_in_full_model() );
		if ( type == CHAIN ) runtime_assert( chains_in_full_model() == get_parameter( CHAIN ) );
	}

} //full_model_info
} //pose
} //core
