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
#include <core/pose/full_model_info/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>

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
																						utility::vector1< Size > & res_list ) {

		initialize_parameters( *this );
		get_sequence_with_gaps_filled_with_n( pose, full_sequence_, conventional_numbering_, res_list );
		set_parameter( CUTPOINT_OPEN, convert_to_parameter_values_at_res( get_cutpoint_open_from_pdb_info( pose ) ) );

		// not sure what's best here -- for now setting that the pose's residues are 'fixed' within the domain map.
		utility::vector1< Size > fixed_domain_map_( full_sequence_.size(), 0 );
		for ( Size n = 1; n <= res_list.size(); n++ ) {
			Size const & res_num = res_list[ n ];
			fixed_domain_map_[ res_num ] = 1;
		}
		set_parameter( FIXED_DOMAIN, fixed_domain_map_ );

	}

	// copy
	FullModelParameters::FullModelParameters( FullModelParameters const & src ) :
		ReferenceCount( src ),
		parameter_values_at_res_( src.parameter_values_at_res_ ),
		parameter_values_as_res_lists_( src.parameter_values_as_res_lists_ ),
		full_sequence_( src.full_sequence_ ),
		conventional_numbering_( src.conventional_numbering_ ),
		conventional_chains_( src.conventional_chains_ )
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
			runtime_assert ( res_list[n] >= 1 && res_list[n] <= size() );
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
																														 utility::vector1< int > & conventional_numbering,
																														 utility::vector1< Size > & res_list
																														 ) const {
		// should also be smart about not filling in n's between chains.
		// anyway. this is a quick hack for now.
		sequence = "";
		conventional_numbering.clear();

		utility::vector1< int > const pdb_res_list = get_res_num_from_pdb_info( pose );
		Size count( 0 );
		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			int const prev_res_num    = ( n > 1 && pose.chain( n-1 ) == pose.chain( n ) ) ? pdb_res_list[ n-1 ] : pdb_res_list[ n ];
			int const current_res_num = pdb_res_list[ n ];
			for ( int i = prev_res_num+1; i < current_res_num; i++ ) {
				sequence.push_back( 'n' );
				conventional_numbering.push_back( i );
				count++;
			}
			sequence.push_back( pose.sequence()[ n-1 ] );
			conventional_numbering.push_back( current_res_num );
			res_list.push_back( ++count );
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
			if ( pose.residue_type( n ).is_protein() != pose.residue_type( n+1 ).is_protein() ) {
				cutpoint_open.push_back( n );
				continue;
			}
			if ( pose.residue_type( n ).is_RNA() != pose.residue_type( n+1 ).is_RNA() ) {
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

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	FullModelParameters::conventional_to_full( utility::vector1< int > const & res_list ) const {
		utility::vector1< Size > res_list_in_full_numbering;
		for ( Size n = 1; n <= res_list.size(); n++ ){
			res_list_in_full_numbering.push_back( conventional_to_full( res_list[n] ) );
		}
		return res_list_in_full_numbering;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	FullModelParameters::conventional_to_full( int const res_num ) const {
		bool found_match( false );
		Size res_num_in_full_numbering( 0 );
		for ( Size n = 1; n <= conventional_numbering_.size() ; n++ ){
			if ( res_num == conventional_numbering_[ n ] ){
				if( found_match ) utility_exit_with_message( "ambiguous res_num (maybe supply chain?) "+ObjexxFCL::string_of(res_num) );
				res_num_in_full_numbering = n;
				found_match = true;
			}
		}
		runtime_assert( found_match );
		return res_num_in_full_numbering;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	FullModelParameters::has_conventional_residue( int const res_num ) const {
		for ( Size n = 1; n <= conventional_numbering_.size() ; n++ ){
			if ( res_num == conventional_numbering_[ n ] )	return true;
		}
		return false;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	FullModelParameters::conventional_to_full( std::pair< utility::vector1< int >, utility::vector1< char > > const & resnum_and_chain ) const {
		utility::vector1< Size > res_list_in_full_numbering;
		utility::vector1< int  > const & resnum = resnum_and_chain.first;
		utility::vector1< char > const & chain  = resnum_and_chain.second;
		runtime_assert( resnum.size() == chain.size() );
		for ( Size n = 1; n <= resnum.size(); n++ ){
			res_list_in_full_numbering.push_back( conventional_to_full( resnum[n], chain[n] ) );
		}
		return res_list_in_full_numbering;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	FullModelParameters::conventional_to_full( int const res_num, char const chain ) const {
		bool found_match( false );
		Size res_num_in_full_numbering( 0 );
		for ( Size n = 1; n <= conventional_numbering_.size() ; n++ ){
			if ( res_num != conventional_numbering_[ n ] ) continue;
			if ( chain != ' ' && conventional_chains_.size() > 0 && conventional_chains_[ n ] != ' ' && chain != conventional_chains_[ n ] ) continue;
			if( found_match ) utility_exit_with_message( "ambiguous res_num & chain "+ObjexxFCL::string_of(res_num)+" "+ObjexxFCL::string_of(chain) );
			res_num_in_full_numbering = n;
			found_match = true;
		}
		runtime_assert( found_match );
		return res_num_in_full_numbering;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	FullModelParameters::has_conventional_residue( int const res_num, char const chain ) const {
		for ( Size n = 1; n <= conventional_numbering_.size() ; n++ ){
			if ( res_num != conventional_numbering_[ n ] ) continue;
			if ( chain != ' ' && conventional_chains_.size() > 0 && conventional_chains_[ n ] != ' ' && chain != conventional_chains_[ n ] ) continue;
			return true;
		}
		return false;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< int >
	FullModelParameters::full_to_conventional( utility::vector1< Size > const & res_list ) const {
		utility::vector1< int > res_list_in_conventional_numbering;
		for ( Size n = 1; n <= res_list.size(); n++ ){
			res_list_in_conventional_numbering.push_back( full_to_conventional( res_list[n] ) );
		}
		return res_list_in_conventional_numbering;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int
	FullModelParameters::full_to_conventional( Size const res_num ) const {
		runtime_assert( res_num >= 1 && res_num <= conventional_numbering_.size() );
		return conventional_numbering_[ res_num ];
	}

} //full_model_info
} //pose
} //core
