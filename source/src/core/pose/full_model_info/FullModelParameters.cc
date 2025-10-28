// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/full_model_info/FullModelParameters.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/rna/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/stream_util.hh> // AUTO IWYU For operator<<

static basic::Tracer TR( "core.pose.full_model_info.FullModelParameters" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace full_model_info {

/////////////////////////////////////////////////////////////////////
// Get rid of this commented code when it is incorporated into a unit test.
/////////////////////////////////////////////////////////////////////
// TR << *const_full_model_info( pose ).full_model_parameters() << std::endl;

// // movnig in and out.
// std::ostringstream os;
// os << *const_full_model_info( pose ).full_model_parameters() << std::endl;
// std::istringstream is( os.str() );

// TR << "BACK OUT" << std::endl;
// FullModelParameters new_full_model_parameters;
// is >> new_full_model_parameters;
// TR << new_full_model_parameters << std::endl;
// return;
/////////////////////////////////////////////////////////////////////

//Constructor
FullModelParameters::FullModelParameters():
	full_sequence_( "" )
{
	initialize_parameters( *this );
}

//Constructor
FullModelParameters::FullModelParameters( std::string const & full_sequence ):
	full_sequence_( full_sequence )
{
	initialize_parameters( *this );
	for ( Size n = 1; n <= core::pose::rna::remove_bracketed( full_sequence ).size(); n++ ) {
		conventional_numbering_.push_back( n );
		conventional_chains_.push_back( " " );
		conventional_segids_.push_back( "    " );
	}
}

//Constructorx
FullModelParameters::FullModelParameters( std::string const & full_sequence,
	utility::vector1< Size > const & cutpoint_open_in_full_model,
	utility::vector1< Size > const & res_numbers_in_pose ):
	full_sequence_( full_sequence )
{
	initialize_parameters( *this );
	set_parameter( CUTPOINT_OPEN, convert_to_parameter_values_at_res( cutpoint_open_in_full_model ) );
	utility::vector1< Size > fixed_domain_map;
	for ( Size n = 1; n <= core::pose::rna::remove_bracketed( full_sequence ).size(); n++ ) {
		if ( res_numbers_in_pose.has( n ) ) {
			fixed_domain_map.push_back( 1 );
		} else {
			fixed_domain_map.push_back( 0 );
		}
		conventional_numbering_.push_back( static_cast<int>( n ) );
		conventional_chains_.push_back( " " );
		conventional_segids_.push_back( "    " );
	}
	set_parameter( FIXED_DOMAIN, fixed_domain_map );
	set_parameter( INPUT_DOMAIN, fixed_domain_map );
	set_parameter( WORKING,      fixed_domain_map );
}

//Constructor
FullModelParameters::FullModelParameters( pose::Pose const & pose,
	utility::vector1< Size > & res_list ) {

	get_sequence_with_gaps_filled_with_n( pose, full_sequence_,
		conventional_numbering_,
		conventional_chains_,
		conventional_segids_,
		res_list );

	initialize_parameters( *this );

	set_parameter( CUTPOINT_OPEN, convert_to_parameter_values_at_res( get_cutpoint_open_from_pdb_info( pose, res_list ) ) );

	// not sure what's best here -- for now setting that the pose's residues are 'fixed' within the domain map.
	utility::vector1< Size > fixed_domain_map( full_sequence_.size(), 0 );
	for ( Size n = 1; n <= res_list.size(); n++ ) {
		Size const & res_num = res_list[ n ];
		fixed_domain_map[ res_num ] = 1;
	}
	set_parameter( FIXED_DOMAIN, fixed_domain_map );
	set_parameter( INPUT_DOMAIN, fixed_domain_map );
	set_parameter( WORKING,      fixed_domain_map );

}

// copy
FullModelParameters::FullModelParameters( FullModelParameters const & /*src*/ ) = default;

//Destructor
FullModelParameters::~FullModelParameters() = default;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::map< Size, utility::vector1< Size > >
FullModelParameters::convert_to_res_lists_by_value( utility::vector1< Size > const & parameter_values_at_res ){
	std::map< Size, utility::vector1< Size > > res_lists_by_value;
	runtime_assert( size() == parameter_values_at_res.size() );

	// basic initialization.
	utility::vector1< Size > blank_vector;
	res_lists_by_value[ 0 ] = blank_vector;
	res_lists_by_value[ 1 ] = blank_vector;

	for ( Size n = 1; n <= parameter_values_at_res.size(); n++ ) {
		Size const & value = parameter_values_at_res[ n ];
		if ( res_lists_by_value.find( value ) == res_lists_by_value.end() ) {
			res_lists_by_value[ value ] = blank_vector;
		}
		res_lists_by_value[ value ].push_back( n );
	}

	return res_lists_by_value;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FullModelParameters::fill_parameter_values( utility::vector1< Size > & parameter_values_at_res,
	Size const idx, utility::vector1< Size > const & res_list ) const {
	for ( Size const resi : res_list ) {
		runtime_assert ( resi >= 1 && resi <= size() );
		parameter_values_at_res[ resi ] = idx;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelParameters::convert_to_parameter_values_at_res( utility::vector1< Size > const & res_list ){
	utility::vector1< Size > parameter_values_at_res( size(), 0 );
	fill_parameter_values( parameter_values_at_res, 1, res_list );
	return parameter_values_at_res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelParameters::convert_to_parameter_values_at_res( std::map< Size, utility::vector1< Size > > const & res_lists ){
	utility::vector1< Size > parameter_values_at_res( size(), 0 );
	for ( auto const & elem : res_lists ) {
		fill_parameter_values( parameter_values_at_res, elem.first, elem.second );
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
utility::vector1< std::pair< Size, Size > >
FullModelParameters::get_res_list_as_pairs( FullModelParameterType const type ) const {
	utility::vector1< std::pair< Size, Size > > res_list_as_pairs;
	std::map< Size, utility::vector1< Size > > const & res_lists = get_parameter_as_res_lists( type );
	for ( auto const & elem : res_lists ) {
		if ( elem.first == 0 ) continue;
		if ( elem.second.size() == 0 ) continue;
		runtime_assert( elem.second.size() == 2 );
		Size const res1_full =  elem.second[1];
		Size const res2_full =  elem.second[2];
		res_list_as_pairs.push_back( std::make_pair( res1_full, res2_full ) );
	}
	return res_list_as_pairs;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size > const &
FullModelParameters::get_parameter( FullModelParameterType const type ) const {
	runtime_assert( parameter_values_at_res_.find( type ) != parameter_values_at_res_.end() );
	return parameter_values_at_res_.find( type )->second;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::map< Size, utility::vector1< Size > > const &
FullModelParameters::get_parameter_as_res_lists( FullModelParameterType const type ) const {
	runtime_assert( parameter_values_as_res_lists_.find( type ) != parameter_values_as_res_lists_.end() );
	return parameter_values_as_res_lists_.find( type )->second;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FullModelParameters::get_sequence_with_gaps_filled_with_n( pose::Pose const & pose,
	std::string & sequence,
	utility::vector1< int >  & conventional_numbering,
	utility::vector1< std::string > & conventional_chains,
	utility::vector1< std::string > & conventional_segids,
	utility::vector1< Size > & res_list
) const {
	// should also be smart about not filling in n's between chains.
	// anyway. this is a quick hack for now.
	sequence = "";
	conventional_numbering.clear();
	conventional_chains.clear();
	conventional_segids.clear();

	utility::vector1< int > const pdb_res_list = get_res_num_from_pdb_info( pose );
	Size count( 0 );
	for ( Size n = 1; n <= pose.size(); n++ ) {
		int const prev_res_num    = ( n > 1 && pose.chain( n-1 ) == pose.chain( n ) ) ? pdb_res_list[ n-1 ] : pdb_res_list[ n ];
		int const current_res_num = pdb_res_list[ n ];
		for ( int i = prev_res_num+1; i < current_res_num; i++ ) {
			sequence.push_back( 'n' );
			conventional_numbering.push_back( i );
			if ( pose.pdb_info() ) {
				conventional_chains.push_back( pose.pdb_info()->chain( n ) );
				conventional_segids.push_back( pose.pdb_info()->segmentID( n ) );
			} else {
				conventional_chains.push_back( " " );
				conventional_segids.push_back( "    " );
			}
			count++;
		}
		sequence.push_back( pose.sequence()[ n-1 ] );
		conventional_numbering.push_back( current_res_num );
		if ( pose.pdb_info() ) {
			conventional_chains.push_back( pose.pdb_info()->chain( n ) );
			conventional_segids.push_back( pose.pdb_info()->segmentID( n ) );
		} else {
			conventional_chains.push_back( " " );
			conventional_segids.push_back( "    " );
		}
		res_list.push_back( ++count );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelParameters::get_cutpoint_open_from_pdb_info( pose::Pose const & pose, utility::vector1< Size > const & res_list ) const {

	PDBInfoCOP pdb_info = pose.pdb_info();
	utility::vector1< Size > cutpoint_open;

	for ( Size n = 1; n < pose.size(); n++ ) {

		if ( ( pdb_info &&  (pdb_info->chain( n ) != pdb_info->chain( n+1 ) ) ) ||
				( pose.residue_type( n ).is_protein() != pose.residue_type( n+1 ).is_protein() ) ||
				( pose.residue_type( n ).is_RNA() != pose.residue_type( n+1 ).is_RNA() ) ||
				( !pose.residue_type( n ).is_polymer() ) ) {
			cutpoint_open.push_back( res_list[ n ] );
		}
	}
	return cutpoint_open;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelParameters::chains_in_full_model() const {
	return get_chains_from_cutpoint_open( get_parameter( CUTPOINT_OPEN ), size() );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FullModelParameters::set_parameter( FullModelParameterType const type,
	utility::vector1< Size > const & setting ){

	parameter_values_at_res_[ type ] = setting;
	parameter_values_as_res_lists_[ type ] = convert_to_res_lists_by_value( setting );

	//  keep_chain_and_cutpoint_open_matched( type );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FullModelParameters::set_parameter_as_res_list( FullModelParameterType const type,
	utility::vector1< Size > const & setting ){
	set_parameter( type, convert_to_parameter_values_at_res( setting ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FullModelParameters::set_parameter_as_res_lists( FullModelParameterType const type,
	std::map< Size, utility::vector1< Size > > const & setting ){
	set_parameter( type, convert_to_parameter_values_at_res( setting ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelParameters::conventional_to_full( utility::vector1< int > const & res_list ) const {
	utility::vector1< Size > res_list_in_full_numbering;
	for ( Size n = 1; n <= res_list.size(); n++ ) {
		res_list_in_full_numbering.push_back( conventional_to_full( res_list[n] ) );
	}
	return res_list_in_full_numbering;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
FullModelParameters::conventional_to_full( int const res_num ) const {
	bool found_match( false );
	Size res_num_in_full_numbering( 0 );
	for ( Size n = 1; n <= conventional_numbering_.size() ; n++ ) {
		if ( res_num == conventional_numbering_[ n ] ) {
			if ( found_match ) utility_exit_with_message( "ambiguous res_num (maybe supply chain?) "+ObjexxFCL::string_of(res_num) );
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
	for ( Size n = 1; n <= conventional_numbering_.size() ; n++ ) {
		if ( res_num == conventional_numbering_[ n ] ) return true;
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelParameters::conventional_to_full( std::tuple< utility::vector1< int >, utility::vector1< std::string >, utility::vector1< std::string > > const & resnum_and_chain_and_segid ) const {
	utility::vector1< Size > res_list_in_full_numbering;
	utility::vector1< int  > const & resnum = std::get< 0 >( resnum_and_chain_and_segid );
	utility::vector1< std::string > const & chain  = std::get< 1 >( resnum_and_chain_and_segid );
	utility::vector1< std::string > const & segid  = std::get< 2 >( resnum_and_chain_and_segid );
	// AMW TODO: There is an issue in how segid is sometimes empty here. It should probably
	// be filled by all clients of this function, so this function can check for that...
	runtime_assert( resnum.size() == chain.size() );
	for ( Size n = 1; n <= resnum.size(); n++ ) {
		if ( !segid.empty() ) {
			res_list_in_full_numbering.push_back( conventional_to_full( resnum[n], chain[n], segid[n] ) );
		} else {
			res_list_in_full_numbering.push_back( conventional_to_full( resnum[n], chain[n], "    " ) );
		}
	}
	return res_list_in_full_numbering;
}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// utility::vector1< Size >
// FullModelParameters::conventional_to_full( std::pair< std::vector< int >, std::vector< char > > const & resnum_and_chain ) const {
//  utility::vector1< Size > res_list_in_full_numbering;
//  std::vector< int  > const & resnum = resnum_and_chain.first;
//  std::vector< char > const & chain  = resnum_and_chain.second;
//  runtime_assert( resnum.size() == chain.size() );
//  for ( Size n = 0; n < resnum.size(); n++ ) {
//   res_list_in_full_numbering.push_back( conventional_to_full( resnum[n], chain[n] ) );
//  }
//  return res_list_in_full_numbering;
// }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
FullModelParameters::conventional_to_full( int const res_num, std::string const & chain, std::string const & segid ) const {
	using namespace ObjexxFCL;
	bool found_match( false );
	Size res_num_in_full_numbering( 0 );
	//TR << "Chains: " << conventional_chains_ << std::endl;
	//TR << "Segid:  " << conventional_segids_ << std::endl;
	//TR << "Resnum: " << conventional_numbering_ << std::endl;

	for ( Size n = 1; n <= conventional_numbering_.size() ; n++ ) {
		//std::cout << " eval " << n << " " << conventional_numbering_[ n ] << " " << conventional_chains_[ n ] << " \"" << conventional_segids_[ n ] << "\"" << std::endl;
		if ( res_num != conventional_numbering_[ n ] ) continue;
		if ( chain != " " && conventional_chains_.size() > 0 && conventional_chains_[ n ] != " " && chain != conventional_chains_[ n ] ) continue;
		if ( segid != "    " && conventional_segids_.size() > 0 && conventional_segids_[ n ] != "    " && segid != conventional_segids_[ n ] ) continue;
		if ( found_match ) utility_exit_with_message( "ambiguous res_num & chain & segid "+ string_of(res_num)+" "+string_of(chain)+" "+segid );
		res_num_in_full_numbering = n;
		found_match = true;
	}
	if ( !found_match ) {
		TR << "Chains: " << conventional_chains_ << std::endl;
		TR << "Segid:  " << conventional_segids_ << std::endl;
		TR << "Resnum: " << conventional_numbering_ << std::endl;
		utility_exit_with_message( "Could not match residue number  " + string_of( res_num ) + " and chain " + '"' + string_of(chain) + '"' + " and segid " + '"' + segid + '"' );
	}
	return res_num_in_full_numbering;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
FullModelParameters::has_conventional_residue( int const res_num, std::string const & chain, std::string const & segid ) const {
	for ( Size n = 1; n <= conventional_numbering_.size() ; n++ ) {
		if ( res_num != conventional_numbering_[ n ] ) continue;
		if ( chain != " " && conventional_chains_.size() > 0 && conventional_chains_[ n ] != " " && chain != conventional_chains_[ n ] ) continue;
		if ( segid != "    " && conventional_segids_.size() > 0 && conventional_segids_[ n ] != "    " && segid != conventional_segids_[ n ] ) continue;
		return true;
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< int >
FullModelParameters::full_to_conventional( utility::vector1< Size > const & res_list ) const {
	utility::vector1< int > res_list_in_conventional_numbering;
	for ( Size n = 1; n <= res_list.size(); n++ ) {
		res_list_in_conventional_numbering.push_back( full_to_conventional( res_list[n] ) );
	}
	return res_list_in_conventional_numbering;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::tuple< utility::vector1< int >, utility::vector1< std::string >, utility::vector1< std::string > >
FullModelParameters::full_to_conventional_resnum_and_chain_and_segid( utility::vector1< Size > const & res_list ) const {
	utility::vector1< int > conventional_numbering_in_res_list;
	utility::vector1< std::string > conventional_chains_in_res_list;
	utility::vector1< std::string > conventional_segids_in_res_list;
	for ( Size n = 1; n <= res_list.size(); n++ ) {
		Size const & res_num = res_list[ n ];
		runtime_assert( res_num >= 1 && res_num <= conventional_numbering_.size() );
		runtime_assert( res_num <= conventional_chains_.size() || conventional_chains_.size() == 0 );
		runtime_assert( res_num <= conventional_segids_.size() || conventional_segids_.size() == 0 );
		conventional_numbering_in_res_list.push_back( conventional_numbering_[ res_num ] );
		if ( conventional_chains_.size() > 0 ) {
			conventional_chains_in_res_list.push_back( conventional_chains_[ res_num ] );
		} else {
			conventional_chains_in_res_list.push_back( " " );
		}
		if ( conventional_segids_.size() > 0 ) {
			conventional_segids_in_res_list.push_back( conventional_segids_[ res_num ] );
		} else {
			conventional_segids_in_res_list.push_back( "    " );
		}
	}
	return std::make_tuple( conventional_numbering_in_res_list, conventional_chains_in_res_list, conventional_segids_in_res_list );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int
FullModelParameters::full_to_conventional( Size const res_num ) const {
	runtime_assert( res_num >= 1 && res_num <= conventional_numbering_.size() );
	return conventional_numbering_[ res_num ];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::tuple< int, std::string, std::string >
FullModelParameters::full_to_conventional_resnum_and_chain_and_segid( Size const res_num ) const {
	runtime_assert( res_num >= 1 && res_num <= conventional_numbering_.size() );
	runtime_assert( res_num <= conventional_chains_.size() || conventional_chains_.size() == 0 );
	runtime_assert( res_num <= conventional_segids_.size() || conventional_segids_.size() == 0 );
	if ( conventional_chains_.size() > 0 ) {
		if ( conventional_segids_.size() > 0 ) {
			return std::make_tuple( conventional_numbering_[ res_num ], conventional_chains_[ res_num ], conventional_segids_[ res_num ] );
		} else {
			return std::make_tuple( conventional_numbering_[ res_num ], conventional_chains_[ res_num ], "    " );
		}
	}
	if ( conventional_segids_.size() > 0 ) {
		return std::make_tuple( conventional_numbering_[ res_num ], " ", conventional_segids_[ res_num ] );
	} else {
		return std::make_tuple( conventional_numbering_[ res_num ], " ", "    " );
	}
}

Size
FullModelParameters::size() const {
	return core::pose::rna::remove_bracketed( full_sequence_ ).size();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details  Nice one-line format summarizing everything in FullModelParameters.
std::ostream &
operator <<( std::ostream & os, FullModelParameters const & t )
{
	os << "FULL_MODEL_PARAMETERS";

	os << "  FULL_SEQUENCE " << t.full_sequence();
	if ( t.global_sequence() != "" ) {
		os << "  GLOBAL_SEQUENCE " << t.global_sequence();
		os << "  GLOBAL_MAPPING ";
		for ( Size i = 1; i <= t.global_mapping_.size() - 1; i++ ) {
			os << t.global_mapping_[i] << ',';
		}
		if ( t.global_mapping_.size() > 0 ) os << t.global_mapping_.back();
	}

	if ( t.conventional_chains().size() > 0 && t.conventional_segids().size() > 0 ) {
		os << "  CONVENTIONAL_RES_CHAIN "  << make_tag_with_dashes( t.conventional_numbering(), t.conventional_chains(), t.conventional_segids(), ',' );
	} else {
		os << "  CONVENTIONAL_RES_CHAIN "  << make_tag_with_dashes( t.conventional_numbering(), ',' );
	}

	for ( Size n = 1; n < LAST_TYPE; n++ ) {
		auto type = static_cast< FullModelParameterType >( n );
		std::map< Size, utility::vector1< Size > > const & res_lists = t.get_parameter_as_res_lists( type );
		runtime_assert( res_lists.size() > 0 );

		bool has_domain_higher_than_one( false );
		for ( auto const & res_list : res_lists ) {
			if ( res_list.first > 1 ) has_domain_higher_than_one = true;
		}

		std::ostringstream os_local;
		for ( auto const & res_list : res_lists ) {
			if ( res_list.first == 0 )         continue; // don't bother with 0.
			if ( res_list.second.size() == 0 ) continue;
			os_local << ' ';
			if ( has_domain_higher_than_one ) os_local << res_list.first << ':'; // give index.
			os_local << make_tag_with_dashes( res_list.second, ',' );
		}
		if ( os_local.str().size() == 0 ) continue;
		os << "  " << to_string( type ) << os_local.str();
	}

	// Heavens! We forgot about the non_standard_residue_map!
	// This is necessary because if only sequence is there, we can't re-parse it
	// to get the NSRM (and not every code path will ensure that consistency).
	// Obviously it would be IDEAL for these things to just be fundamental, managed
	// pose data but at the moment we should be careful to explicitly read them out/in.
	if ( !t.non_standard_residue_map().empty() ) {
		os << "  NON_STANDARD_RESIDUE_MAP ";
		for ( auto const & elem : t.non_standard_residue_map() ) {
			os << elem.first << ":" << elem.second << " ";
		}
	}

	// constraints are a special block -- they have atom_names and do not fit with above format.
	using namespace core::scoring::constraints;
	if ( t.cst_string_.size() > 0 ) os << "   CONSTRAINTS " << utility::replace_in( t.cst_string_, "\n", "; " );
	return os;
}


/////////////////////////////////////////////////////////////////////////////
/// @details Read in of one-line format for FullModelParameters -- better be exact reverse of <<
std::istream &
operator >>( std::istream & is, FullModelParameters & t )
{
	using namespace utility;
	std::string tag;

	is >> tag;
	runtime_assert ( !is.fail() && tag == "FULL_MODEL_PARAMETERS" );

	is >> tag;
	runtime_assert ( !is.fail() && tag == "FULL_SEQUENCE" );
	is >> t.full_sequence_;

	initialize_parameters( t ); // depends on size of full_sequence

	is >> tag;
	if ( !is.fail() && tag == "GLOBAL_SEQUENCE" ) {
		is >> t.global_sequence_;
		is >> tag;
		runtime_assert( !is.fail() && (tag == "GLOBAL_MAPPING") );
		is >> tag;
		utility::vector1< std::string > global_idxs = string_split( tag, ',');
		for ( Size i = 1; i <= global_idxs.size(); i++ ) {
			t.global_mapping_.push_back( string2int(global_idxs[i]) );
		}
		is >> tag;
	}
	std::tuple< utility::vector1<int>, utility::vector1<std::string>, utility::vector1< std::string > > resnum_chain;
	runtime_assert ( !is.fail() && ( tag == "CONVENTIONAL_RES_CHAIN_SEGID" || tag == "CONVENTIONAL_RES_CHAIN" ) );
	for ( bool ok = true; ok ; ) {
		is >> tag;
		resnum_chain = get_resnum_and_chain_and_segid( tag, ok );
		if ( ok ) {
			t.conventional_numbering_ = std::get< 0 >( resnum_chain );
			t.conventional_chains_    = std::get< 1 >( resnum_chain );
			t.conventional_segids_    = std::get< 2 >( resnum_chain );
		}
	}

	while ( !is.fail() && tag != "CONSTRAINTS" && tag != "NON_STANDARD_RESIDUE_MAP" /*special*/ ) {
		FullModelParameterType type = full_model_parameter_type_from_string( tag );
		utility::vector1< Size > parameter_values_at_res( t.size(), 0 );

		bool ok( true );
		while ( ok )  {
			// the 'chain' will be any string before a ':'. Could be blank.
			is >> tag;
			if ( tag == "NON_STANDARD_RESIDUE_MAP" ) break;
			TR.Debug << "Found tag " << tag << std::endl;

			utility::vector1< std::string > cols = string_split( tag, ':');
			runtime_assert( cols.size() == 1 || cols.size() == 2 );
			Size idx = ( cols.size() == 1 ) ? 1 : ObjexxFCL::int_of( cols[1] );
			std::vector< int > const resnum =  ObjexxFCL::ints_of( cols[ cols.size() ], ok );
			if ( ok ) {
				for ( int q : resnum ) parameter_values_at_res[ q ] = idx;
			}
			if ( is.fail() ) break;
		}
		t.set_parameter( type, parameter_values_at_res );
	}

	// don't do this -- we just got it. is >> tag;
	if ( tag == "NON_STANDARD_RESIDUE_MAP" ) {
		TR.Debug << "Found tag NON_STANDARD_RESIDUE_MAP " << std::endl;
		bool ok = true;
		while ( ok ) {
			// Read in size : string;
			is >> tag;

			TR.Debug << "Found tag " << tag << std::endl;
			utility::vector1< std::string > cols = string_split( tag, ':');
			// We used to assert there were == two columns. That isn't right: this could
			// have patched residues. (For example, in DNA-as-RNA stepwise, where we
			// add deoxxy residues.)
			runtime_assert( cols.size() >= 2 );
			Size pos = ObjexxFCL::int_of( cols[1] );
			if ( !pos ) break; // tag may be loaded up with CONSTRAINTS or done entirely
			t.non_standard_residue_map_nonconst()[ pos ] = cols[2];
			if ( is.fail() ) break;
		}
	}

	is >> tag;
	if ( tag == "CONSTRAINTS" ) {
		getline( is, t.cst_string_ ); // assume rest of line is taken up by constraints.
		t.cst_string_ = utility::replace_in( t.cst_string_, "; ", "\n"); // convert separator
	}
	return is;
}

/// @brief equal to operator
bool
operator==(
	FullModelParameters const & a,
	FullModelParameters const & b
){
	std::ostringstream ss_a, ss_b;
	ss_a << a;
	ss_b << b;
	return ( ss_a.str() == ss_b.str() );
}

/// @brief equal to operator
bool
operator!=(
	FullModelParameters const & a,
	FullModelParameters const & b
){
	return !( a == b );
}

///////////////////////////////////////////////////////////////////////////////////
// this is pretty clumsy -- keeping the constraints in a string format.
//  however... otherwise we need to instantiate a pose, and we don't always
//  know what the residue type set is.
void
FullModelParameters::read_cst_file( std::string const & cst_file ) {
	cst_string_ = "";
	full_model_pose_for_constraints_ = nullptr;
	cst_set_ = nullptr;
	runtime_assert( utility::file::file_exists( cst_file ) );
	utility::io::izstream data( cst_file.c_str() );
	while ( data.good() ) {
		std::string line;
		getline( data, line );
		if ( line[0] == '#' || line[0] == '\n' || line.size() == 0 ) continue;
		if ( cst_string_.size() > 0 ) cst_string_ += '\n';
		cst_string_ += line; //.substr( 0, line.size() - 1  ); // strip off endline.
	}
}


/////////////////////////////////////////////////////////////////////
void
FullModelParameters::update_pose_and_cst_set_from_cst_string( chemical::ResidueTypeSet const & rsd_type_set ) const {
	using namespace core::scoring::constraints;

	if ( full_model_pose_for_constraints_ != nullptr && cst_set_ != nullptr ) return;

	std::istringstream data( cst_string_ );

	PoseOP full_model_pose( new Pose );
	make_pose_from_sequence( *full_model_pose, full_sequence_, rsd_type_set, false /*auto termini*/ );

	PDBInfoOP pdb_info( new PDBInfo( *full_model_pose ) );
	pdb_info->set_numbering( conventional_numbering_ );
	pdb_info->set_chains( conventional_chains_ );
	full_model_pose->pdb_info( pdb_info );
	full_model_pose_for_constraints_ = full_model_pose;


	ConstraintSetOP full_model_cst_set( new ConstraintSet );
	ConstraintIO::read_constraints_new( data, full_model_cst_set, *full_model_pose );
	cst_set_ = full_model_cst_set;
}

////////////////////////////////////////////////////
scoring::constraints::ConstraintSetCOP
FullModelParameters::cst_set() const {
	runtime_assert( cst_set_ != nullptr ); // make sure that update_pose_and_cst_set_from_cst_string() has been called;
	return cst_set_;
}

////////////////////////////////////////////////////
Pose const &
FullModelParameters::full_model_pose_for_constraints() const {
	runtime_assert( full_model_pose_for_constraints_ != nullptr ); // make sure that update_pose_and_cst_set_from_cst_string() has been called;
	return *full_model_pose_for_constraints_;
}

////////////////////////////////////////////////////
void
FullModelParameters::read_global_seq_info( std::string const & global_seq_file ) {
	if ( global_seq_file.size() == 0 ) {
		return;
	}

	// Read as FASTA
	utility::vector1< core::sequence::SequenceOP > fasta_sequences = core::sequence::read_fasta_file( global_seq_file );
	if ( fasta_sequences.size() > 1 ) {
		TR.Warning << "Can't handle more than one global sequence for now." << std::endl;
	}

	// Fill in global sequence
	global_sequence_ = fasta_sequences[ 1 ]->sequence();

	// Fill in global mapping
	utility::vector1< std::string > global_to_pdb_chains;
	utility::vector1< std::string > global_to_pdb_segids;
	utility::vector1< int  > global_to_pdb_numbering; // Length is same as global sequence, values are PDB numbering
	core::sequence::get_conventional_chains_and_numbering( fasta_sequences, global_to_pdb_chains, global_to_pdb_numbering, global_to_pdb_segids );


	for ( Size fullmodel_idx = 1; fullmodel_idx <= conventional_numbering_.size(); fullmodel_idx++ ) {
		for ( Size global_idx = 1; global_idx <= global_to_pdb_numbering.size(); global_idx++ ) {
			// Do we need to also compare segids? what's a segid?
			if ( ( global_to_pdb_chains[ global_idx ] == conventional_chains_[ fullmodel_idx ] ) &&
					( global_to_pdb_numbering[ global_idx ] == conventional_numbering_[ fullmodel_idx ] ) ) {
				global_mapping_.push_back( global_idx );
				break;
			}
		}
	}
}

////////////////////////////////////////////////////
void
FullModelParameters::read_disulfides( std::string const & disulfide_file ) {

	if ( disulfide_file.size() == 0 ) return;
	utility::vector1< Size > disulfide_res_list;

	// could use core::io::raw_data::DisulfideFile
	// but, strangely, it does not seem to provide the
	// chain/residue renumbering functionality that it should.
	utility::io::izstream data( disulfide_file.c_str() );
	//Size count( 0 );
	while ( data.good() ) {
		std::string line;
		getline( data, line );
		if ( line[0] == '#' || line[0] == '\n' || line.size() == 0 ) continue;

		bool string_is_ok( false );
		std::tuple< utility::vector1< int >, utility::vector1< std::string >, utility::vector1< std::string > > resnum_and_chain = utility::get_resnum_and_chain_and_segid( line, string_is_ok );
		runtime_assert( string_is_ok );
		runtime_assert( std::get< 0 >( resnum_and_chain ).size() == 2 );
		runtime_assert( std::get< 1 >( resnum_and_chain ).size() == 2 );
		runtime_assert( std::get< 2 >( resnum_and_chain ).size() == 2 );
		//count++;
		for ( int k = 0; k <= 1; k++ ) {
			Size const full_model_number = conventional_to_full( std::get< 0 >( resnum_and_chain )[k], std::get< 1 >( resnum_and_chain )[k], std::get< 2 >( resnum_and_chain )[k] );
			runtime_assert( full_sequence_[ full_model_number - 1 ] == 'C' ); // better be cysteine.
			disulfide_res_list.push_back( full_model_number );
		}
	}
	data.close();

	set_parameter_as_res_list_in_pairs( DISULFIDE, disulfide_res_list );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// for disulfides, extra_min_jump_res, secstruct, etc. -- list of mutually exclusive pairs
void
FullModelParameters::set_parameter_as_res_list_in_pairs( FullModelParameterType const type, utility::vector1< Size > const & setting ) {
	runtime_assert( setting.size() % 2 == 0 ); // must be even
	utility::vector1< Size > parameter( size(), 0 );
	for ( Size n = 1; n <= setting.size()/2; n++ ) {
		for ( Size k = (2*n-1); k <= (2*n); k++ ) {
			runtime_assert( parameter[ setting[ k ] ] == 0 );
			parameter[ setting[ k ] ] =  n;
		}
	}

	set_parameter( type, parameter );
}


////////////////////////////////////////////////////////////////////////////////
// Create a sliced out version of this full_model_parameters and return it --
//  save the slice_res_list and 'parent' full_model_parameters inside the new object.
FullModelParametersOP
FullModelParameters::slice( utility::vector1< Size > const & slice_res ) const
{
	// following asserts may not actually be necessary -- perhaps at some
	// future time we'll want to have recursive slicing.
	runtime_assert( slice_res_list_.size() == 0 );
	runtime_assert( parent_full_model_parameters_ == nullptr );

	std::string full_sequence_new;
	utility::vector1< int >  conventional_numbering_new; // permits user to use numbering other than 1, 2, 3...
	utility::vector1< std::string > conventional_chains_new;    // permits user to use chains other than A, B, C, ..
	utility::vector1< std::string > conventional_segids_new;
	std::map< Size, std::string > non_standard_residue_map_new; // for DNA, non-natural protein/RNA, ligands, ions, etc.

	std::map< FullModelParameterType, utility::vector1< Size > > parameter_values_at_res_new;
	for ( Size q = 1; q < LAST_TYPE; q++ ) {
		auto type = static_cast< FullModelParameterType >( q );
		parameter_values_at_res_new[ type ] = utility::vector1< Size >();
	}

	// Needed for properly setting CUTPOINT_OPEN, which marks boundaries between chains:
	utility::vector1< Size > const chains = chains_in_full_model(); // derived from CUTPOINT_OPEN
	utility::vector1< Size > chains_new;

	for ( Size i = 1; i <= slice_res.size(); i++ ) {
		Size const & n = slice_res[ i ];
		full_sequence_new += full_sequence_[ n - 1 ];
		conventional_numbering_new.push_back( conventional_numbering_[ n ] );
		conventional_chains_new.push_back( conventional_chains_[ n ] );
		conventional_segids_new.push_back( conventional_segids_[ n ] );

		auto const it = non_standard_residue_map_.find( n );
		if ( it != non_standard_residue_map_.end() ) non_standard_residue_map_new[ i ] = it->second;

		chains_new.push_back( chains[ n ] );

		for ( Size q = 1; q < LAST_TYPE; q++ ) {
			auto type = static_cast< FullModelParameterType >( q );
			parameter_values_at_res_new[ type ].push_back( get_parameter( type )[ n ] );
		}
	}

	FullModelParametersOP full_model_parameters_new( new FullModelParameters( full_sequence_new ) );
	full_model_parameters_new->set_conventional_numbering( conventional_numbering_new );
	full_model_parameters_new->set_conventional_chains( conventional_chains_new );
	full_model_parameters_new->set_conventional_segids( conventional_segids_new );
	full_model_parameters_new->set_non_standard_residue_map( non_standard_residue_map_new );
	for ( Size q = 1; q < LAST_TYPE; q++ ) {
		auto type = static_cast< FullModelParameterType >( q );
		full_model_parameters_new->set_parameter( type, parameter_values_at_res_new[ type ] );
	}
	full_model_parameters_new->set_parameter_as_res_list( CUTPOINT_OPEN, get_cutpoint_open_from_chains( chains_new ) );
	full_model_parameters_new->set_cst_string( cst_string_ );

	// records information on parentage.
	full_model_parameters_new->set_slice_res_list( slice_res );
	full_model_parameters_new->set_parent_full_model_parameters( this->shared_from_this() );

	return full_model_parameters_new;
}

std::string
FullModelParameters::full_annotated_sequence() const
{
	std::string annotated_sequence = "";
	for ( Size n = 1; n <= full_sequence_.size(); n++ ) {
		annotated_sequence += full_sequence_[ n - 1 ];
		auto it =  non_standard_residue_map_.find( n );
		if ( it == non_standard_residue_map_.end() ) continue;
		annotated_sequence += "[" + it->second + "]";
	}
	return annotated_sequence;
}

} //full_model_info
} //pose
} //core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::full_model_info::FullModelParameters::save( Archive & arc ) const {
	arc( CEREAL_NVP( full_sequence_ ) ); // std::string
	arc( CEREAL_NVP( global_sequence_ ) ); // std::string
	arc( CEREAL_NVP( global_mapping_ ) ); // utility::vector1< Size >
	arc( CEREAL_NVP( conventional_numbering_ ) ); // utility::vector1<int>
	arc( CEREAL_NVP( conventional_chains_ ) ); // utility::vector1<char>
	arc( CEREAL_NVP( conventional_segids_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( non_standard_residue_map_ ) ); // std::map<Size, std::string>
	arc( CEREAL_NVP( cst_string_ ) ); // std::string
	arc( CEREAL_NVP( cst_set_ ) ); // core::scoring::constraints::ConstraintSetCOP
	arc( CEREAL_NVP( full_model_pose_for_constraints_ ) ); // pose::PoseCOP
	arc( CEREAL_NVP( parameter_values_at_res_ ) ); // std::map<FullModelParameterType, utility::vector1<Size> >
	arc( CEREAL_NVP( parameter_values_as_res_lists_ ) ); // std::map<FullModelParameterType, std::map<Size, utility::vector1<Size> > >
	arc( CEREAL_NVP( slice_res_list_ ) );
	arc( CEREAL_NVP( parent_full_model_parameters_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::full_model_info::FullModelParameters::load( Archive & arc ) {
	arc( full_sequence_ ); // std::string
	arc( global_sequence_ ); // std::string
	arc( global_mapping_ ); // utility::vector1< Size >
	arc( conventional_numbering_ ); // utility::vector1<int>
	arc( conventional_chains_ ); // utility::vector1<char>
	arc( conventional_segids_ ); // utility::vector1<std::string>
	arc( non_standard_residue_map_ ); // std::map<Size, std::string>
	arc( cst_string_ ); // std::string
	std::shared_ptr< core::scoring::constraints::ConstraintSet > local_cst_set;
	arc( local_cst_set ); // core::scoring::constraints::ConstraintSetCOP
	cst_set_ = local_cst_set; // copy the non-const pointer(s) into the const pointer(s)
	std::shared_ptr< core::pose::Pose > local_full_model_pose_for_constraints;
	arc( local_full_model_pose_for_constraints ); // pose::PoseCOP
	full_model_pose_for_constraints_ = local_full_model_pose_for_constraints; // copy the non-const pointer(s) into the const pointer(s)
	arc( parameter_values_at_res_ ); // std::map<FullModelParameterType, utility::vector1<Size> >
	arc( parameter_values_as_res_lists_ ); // std::map<FullModelParameterType, std::map<Size, utility::vector1<Size> > >
	arc( slice_res_list_ );
	std::shared_ptr< FullModelParameters > local_parent_full_model_parameters;
	arc( local_parent_full_model_parameters );
	parent_full_model_parameters_ = local_parent_full_model_parameters;
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::full_model_info::FullModelParameters );
CEREAL_REGISTER_TYPE( core::pose::full_model_info::FullModelParameters )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_full_model_info_FullModelParameters )
#endif // SERIALIZATION
