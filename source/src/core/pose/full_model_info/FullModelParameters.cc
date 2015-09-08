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
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>

static thread_local basic::Tracer TR( "core.pose.full_model_info.FullModelParameters" );

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
FullModelParameters::FullModelParameters( std::string const full_sequence ):
	full_sequence_( full_sequence )
{
	initialize_parameters( *this );
	for ( Size n = 1; n <= full_sequence.size(); n++ ) {
		conventional_numbering_.push_back( n );
		conventional_chains_.push_back( ' ' );
	}
}

//Constructorx
FullModelParameters::FullModelParameters( std::string const full_sequence,
	utility::vector1< Size > const & cutpoint_open_in_full_model,
	utility::vector1< Size > const & res_numbers_in_pose ):
	full_sequence_( full_sequence )
{
	initialize_parameters( *this );
	set_parameter( CUTPOINT_OPEN, convert_to_parameter_values_at_res( cutpoint_open_in_full_model ) );
	utility::vector1< Size > fixed_domain_map;
	for ( Size n = 1; n <= full_sequence.size(); n++ ) {
		if ( res_numbers_in_pose.has( n ) ) {
			fixed_domain_map.push_back( 1 );
		} else {
			fixed_domain_map.push_back( 0 );
		}
		conventional_numbering_.push_back( static_cast<int>( n ) );
		conventional_chains_.push_back( ' ' );
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
		res_list );

	initialize_parameters( *this );

	set_parameter( CUTPOINT_OPEN, convert_to_parameter_values_at_res( get_cutpoint_open_from_pdb_info( pose ) ) );

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
FullModelParameters::FullModelParameters( FullModelParameters const & src ) :
	ReferenceCount( src ),
	full_sequence_( src.full_sequence_ ),
	conventional_numbering_( src.conventional_numbering_ ),
	conventional_chains_( src.conventional_chains_ ),
	cst_string_( src.cst_string_ ),
	cst_set_( src.cst_set_ ),
	full_model_pose_for_constraints_( src.full_model_pose_for_constraints_ ),
	parameter_values_at_res_( src.parameter_values_at_res_ ),
	parameter_values_as_res_lists_( src.parameter_values_as_res_lists_ )
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
	for ( Size n = 1; n <= res_list.size(); n++ ) {
		runtime_assert ( res_list[n] >= 1 && res_list[n] <= size() );
		parameter_values_at_res[ res_list[n] ] = idx;
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
	for ( std::map< Size, utility::vector1< Size > >::const_iterator it = res_lists.begin(), end = res_lists.end();
			it != end; ++it ) {
		fill_parameter_values( parameter_values_at_res, it->first, it->second );
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
	for ( std::map< Size, utility::vector1< Size > >::const_iterator it = res_lists.begin(), end = res_lists.end();
			it != end; ++it ) {
		if ( it->first == 0 ) continue;
		if ( it->second.size() == 0 ) continue;
		runtime_assert( it->second.size() == 2 );
		Size const & res1_full =  it->second[1];
		Size const & res2_full =  it->second[2];
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
	utility::vector1< char > & conventional_chains,
	utility::vector1< Size > & res_list
) const {
	// should also be smart about not filling in n's between chains.
	// anyway. this is a quick hack for now.
	sequence = "";
	conventional_numbering.clear();
	conventional_chains.clear();

	utility::vector1< int > const pdb_res_list = get_res_num_from_pdb_info( pose );
	Size count( 0 );
	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		int const prev_res_num    = ( n > 1 && pose.chain( n-1 ) == pose.chain( n ) ) ? pdb_res_list[ n-1 ] : pdb_res_list[ n ];
		int const current_res_num = pdb_res_list[ n ];
		for ( int i = prev_res_num+1; i < current_res_num; i++ ) {
			sequence.push_back( 'n' );
			conventional_numbering.push_back( i );
			if ( pose.pdb_info() ) {
				conventional_chains.push_back( pose.pdb_info()->chain( n ) );
			} else {
				conventional_chains.push_back( ' ' );
			}
			count++;
		}
		sequence.push_back( pose.sequence()[ n-1 ] );
		conventional_numbering.push_back( current_res_num );
		if ( pose.pdb_info() ) {
			conventional_chains.push_back( pose.pdb_info()->chain( n ) );
		} else {
			conventional_chains.push_back( ' ' );
		}
		res_list.push_back( ++count );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelParameters::get_cutpoint_open_from_pdb_info( pose::Pose const & pose ) const {

	PDBInfoCOP pdb_info = pose.pdb_info();
	utility::vector1< Size > cutpoint_open;

	for ( Size n = 1; n < pose.total_residue(); n++ ) {

		if ( pdb_info &&  (pdb_info->chain( n ) != pdb_info->chain( n+1 )) ) {
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
FullModelParameters::conventional_to_full( std::pair< utility::vector1< int >, utility::vector1< char > > const & resnum_and_chain ) const {
	utility::vector1< Size > res_list_in_full_numbering;
	utility::vector1< int  > const & resnum = resnum_and_chain.first;
	utility::vector1< char > const & chain  = resnum_and_chain.second;
	runtime_assert( resnum.size() == chain.size() );
	for ( Size n = 1; n <= resnum.size(); n++ ) {
		res_list_in_full_numbering.push_back( conventional_to_full( resnum[n], chain[n] ) );
	}
	return res_list_in_full_numbering;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
FullModelParameters::conventional_to_full( int const res_num, char const chain ) const {
	using namespace ObjexxFCL;
	bool found_match( false );
	Size res_num_in_full_numbering( 0 );
	for ( Size n = 1; n <= conventional_numbering_.size() ; n++ ) {
		if ( res_num != conventional_numbering_[ n ] ) continue;
		if ( chain != ' ' && conventional_chains_.size() > 0 && conventional_chains_[ n ] != ' ' && chain != conventional_chains_[ n ] ) continue;
		if ( found_match ) utility_exit_with_message( "ambiguous res_num & chain "+ string_of(res_num)+" "+string_of(chain) );
		res_num_in_full_numbering = n;
		found_match = true;
	}
	if ( !found_match ) utility_exit_with_message( "Could not match residue number  " + string_of( res_num ) + " and chain " + string_of(chain) );
	return res_num_in_full_numbering;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
FullModelParameters::has_conventional_residue( int const res_num, char const chain ) const {
	for ( Size n = 1; n <= conventional_numbering_.size() ; n++ ) {
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
	for ( Size n = 1; n <= res_list.size(); n++ ) {
		res_list_in_conventional_numbering.push_back( full_to_conventional( res_list[n] ) );
	}
	return res_list_in_conventional_numbering;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::pair< utility::vector1< int >, utility::vector1< char > >
FullModelParameters::full_to_conventional_resnum_and_chain( utility::vector1< Size > const & res_list ) const {
	utility::vector1< int > conventional_numbering_in_res_list;
	utility::vector1< char > conventional_chains_in_res_list;
	for ( Size n = 1; n <= res_list.size(); n++ ) {
		Size const & res_num = res_list[ n ];
		runtime_assert( res_num >= 1 && res_num <= conventional_numbering_.size() );
		runtime_assert( res_num <= conventional_chains_.size() || conventional_chains_.size() == 0 );
		conventional_numbering_in_res_list.push_back( conventional_numbering_[ res_num ] );
		if ( conventional_chains_.size() > 0 ) {
			conventional_chains_in_res_list.push_back( conventional_chains_[ res_num ] );
		} else {
			conventional_chains_in_res_list.push_back( ' ' );
		}
	}
	return std::make_pair( conventional_numbering_in_res_list, conventional_chains_in_res_list );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int
FullModelParameters::full_to_conventional( Size const res_num ) const {
	runtime_assert( res_num >= 1 && res_num <= conventional_numbering_.size() );
	return conventional_numbering_[ res_num ];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::pair< int, char >
FullModelParameters::full_to_conventional_resnum_and_chain( Size const res_num ) const {
	runtime_assert( res_num >= 1 && res_num <= conventional_numbering_.size() );
	runtime_assert( res_num <= conventional_chains_.size() || conventional_chains_.size() == 0 );
	if ( conventional_chains_.size() > 0 ) {
		return std::make_pair( conventional_numbering_[ res_num ], conventional_chains_[ res_num ] );
	}
	return std::make_pair( conventional_numbering_[ res_num ], ' ' );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details  Nice one-line format summarizing everything in FullModelParameters.
std::ostream &
operator <<( std::ostream & os, FullModelParameters const & t )
{
	os << "FULL_MODEL_PARAMETERS";

	os << "  FULL_SEQUENCE " << t.full_sequence();
	if ( t.conventional_chains().size() > 0 ) {
		os << "  CONVENTIONAL_RES_CHAIN "  << make_tag_with_dashes( t.conventional_numbering(), t.conventional_chains(), ',' );
	} else {
		os << "  CONVENTIONAL_RES_CHAIN "  << make_tag_with_dashes( t.conventional_numbering(), ',' );
	}

	for ( Size n = 1; n < LAST_TYPE; n++ ) {
		FullModelParameterType type = static_cast< FullModelParameterType >( n );
		std::map< Size, utility::vector1< Size > > const & res_lists = t.get_parameter_as_res_lists( type );
		runtime_assert( res_lists.size() > 0 );

		bool has_domain_higher_than_one( false );
		for ( std::map< Size, utility::vector1< Size > >::const_iterator it = res_lists.begin(), end = res_lists.end();
				it != end; ++it ) {
			if ( it->first > 1 ) has_domain_higher_than_one = true;
		}

		std::ostringstream os_local;
		for ( std::map< Size, utility::vector1< Size > >::const_iterator it = res_lists.begin(), end = res_lists.end();
				it != end; ++it ) {
			if ( it->first == 0 )         continue; // don't bother with 0.
			if ( it->second.size() == 0 ) continue;
			os_local << ' ';
			if ( has_domain_higher_than_one ) os_local << it->first << ':'; // give index.
			os_local << make_tag_with_dashes( it->second, ',' );
		}
		if ( os_local.str().size() == 0 ) continue;
		os << "  " << to_string( type ) << os_local.str();
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
	std::pair< std::vector<int>, std::vector<char> > resnum_chain;
	runtime_assert ( !is.fail() && tag == "CONVENTIONAL_RES_CHAIN" );
	bool ok( true );
	while ( ok ) {
		is >> tag;
		resnum_chain = get_resnum_and_chain( tag, ok );
		if ( ok ) {
			t.conventional_numbering_ = resnum_chain.first;
			t.conventional_chains_    = resnum_chain.second;
		}
	}

	while ( !is.fail() && tag != "CONSTRAINTS" /*special*/ ) {
		FullModelParameterType type = full_model_parameter_type_from_string( tag );
		utility::vector1< Size > parameter_values_at_res( t.size(), 0 );

		bool ok( true );
		while ( ok )  {
			// the 'chain' will be any string before a ':'. Could be blank.
			is >> tag;
			utility::vector1< std::string > cols = string_split( tag, ':');
			runtime_assert( cols.size() == 1 || cols.size() == 2 );
			Size idx = ( cols.size() == 1 ) ? 1 : ObjexxFCL::int_of( cols[1] );
			std::vector< int > const resnum =  ObjexxFCL::ints_of( cols[ cols.size() ], ok );
			if ( ok ) {
				for ( Size q = 0; q < resnum.size(); q++ ) parameter_values_at_res[ resnum[q] ] = idx;
			}
			if ( is.fail() ) break;
		}
		t.set_parameter( type, parameter_values_at_res );
	}

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
FullModelParameters::read_cst_file( std::string const cst_file ) {
	cst_string_ = "";
	full_model_pose_for_constraints_ = 0;
	cst_set_ = 0;
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

	if ( full_model_pose_for_constraints_ != 0 && cst_set_ != 0 ) return;

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
	runtime_assert( cst_set_ != 0 ); // make sure that update_pose_and_cst_set_from_cst_string() has been called;
	return cst_set_;
}

////////////////////////////////////////////////////
Pose const &
FullModelParameters::full_model_pose_for_constraints() const {
	runtime_assert( full_model_pose_for_constraints_ != 0 ); // make sure that update_pose_and_cst_set_from_cst_string() has been called;
	return *full_model_pose_for_constraints_;
}

////////////////////////////////////////////////////
void
FullModelParameters::read_disulfides( std::string const disulfide_file ) {

	if ( disulfide_file.size() == 0 ) return;
	utility::vector1< Size > disulfide_res_list;

	// could use core::io::raw_data::DisulfideFile
	// but, strangely, it does not seem to provide the
	// chain/residue renumbering functionality that it should.
	utility::io::izstream data( disulfide_file.c_str() );
	Size count( 0 );
	while ( data.good() ) {
		std::string line;
		getline( data, line );
		if ( line[0] == '#' || line[0] == '\n' || line.size() == 0 ) continue;

		bool string_is_ok( false );
		std::pair< std::vector< int >, std::vector< char > > resnum_and_chain = utility::get_resnum_and_chain( line, string_is_ok );
		runtime_assert( string_is_ok );
		runtime_assert( resnum_and_chain.first.size()  == 2);
		runtime_assert( resnum_and_chain.second.size() == 2);
		count++;
		for ( int k = 0; k <= 1; k++ ) {
			Size const full_model_number = conventional_to_full( resnum_and_chain.first[k], resnum_and_chain.second[k] );
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

} //full_model_info
} //pose
} //core
