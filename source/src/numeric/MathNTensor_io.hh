// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


//////////////////////////////////////////////////////////////////////
///
/// @brief
/// File input/output for the MathNTensor class
///
/// @details
///
/// To avoid cryptic 'hardwiring' of binary file specification into
///  the code, look for details of # bins  as n_bins entry in JSON file.
///
/// File should be of the form '.bin.gz', with associated '.json' ASCII file.
///
/// @author Rhiju Das
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_numeric_MathNTensor_io_hh
#define INCLUDED_numeric_MathNTensor_io_hh

#include <numeric/MathNTensor.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>
#include <utility/tools/make_vector.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/json_spirit/json_spirit.h>
#include <utility/json_spirit/json_spirit_tools.hh>

namespace numeric {

template< class T, numeric::Size N >
bool
write_tensor_to_file_without_json( std::string const & filename,
	MathNTensor< T, N > const & tensor );


// @brief reads tensor from file. *requires* a .json file with n_bins & type.
template< class T, numeric::Size N >
void
read_tensor_from_file( std::string const & filename_input,
	MathNTensor< T, N > & tensor,
	utility::json_spirit::mObject & json )
{
	using namespace utility;
	using namespace utility::json_spirit;
	using namespace utility::file;
	std::string filename( filename_input ), filename_prefix( "" ), binary_file( "" );
	bool use_binary( true );
	if ( filename.find( ".txt" ) == filename.size() - 4 ||
			filename.find( ".txt.gz" ) == filename.size() - 7 ) {
		use_binary = false;
		filename_prefix = replace_in( replace_in( filename, ".txt.gz", "" ), ".txt", "" );
		binary_file = replace_in( filename, ".txt", ".bin" );
		if ( file::file_exists( binary_file ) ) {
			filename = binary_file;
			use_binary = true;
		}
	}
	if ( use_binary ) {
		runtime_assert( filename.find( ".bin" ) == filename.size() - 4 ||
			filename.find( ".bin.gz" ) == filename.size() - 7 );
		filename_prefix = replace_in( replace_in( filename, ".bin.gz", "" ), ".bin", "" );
	}
	if ( filename_prefix.size() == 0 ) utility_exit_with_message( "tensor filename must end in .bin, .bin.gz, .txt, or .txt.gz" );
	std::string const json_file( filename_prefix + ".json" );

	if ( !file::file_exists( json_file ) ) utility_exit_with_message( "Tensor_file should come with accompanying .json file" );

	utility::io::izstream json_stream( json_file );
	mValue json_data;
	read_or_throw( json_stream, json_data );
	json_stream.close();

	json = json_data.get_obj();
	mArray n_bins( get_mArray( json, "n_bins" ) );
	std::string const & type = get_string_or_empty( json, "type" );
	if ( ( type == "double" && !std::is_same<T,double>::value ) ||
			( type == "uint64" && !std::is_same<T,uint64_t>::value ) ) {
		utility_exit_with_message( filename+ " type in json " +type + " does not match requested tensor type." );
	}
	runtime_assert( N == n_bins.size() );
	Size tensor_size( 1 ), i( 0 );
	utility::fixedsizearray1< Size, N > nbinsarray;
	for ( auto n_bin : n_bins ) {
		nbinsarray[++i] = n_bin.get_int();
		tensor_size *= n_bin.get_int();
	}

	T * data =  new T[ tensor_size ];
	if ( use_binary ) {
		utility::io::izstream tensor_stream( filename, std::ios::in | std::ios::binary );
		tensor_stream.read( ( char* ) & data[0], tensor_size * sizeof( T )  );
		tensor_stream.close();
	} else {
		utility::io::izstream tensor_stream( filename, std::ios::in );
		for ( Size n = 0; n < tensor_size; n++ ) tensor_stream >> data[ n ];
		tensor_stream.close();
	}

	tensor = numeric::MathNTensor< double, 6 >( nbinsarray, data );

	if ( !use_binary ) {
		// try to create a binary variant -- should allow for fast i/o next time.
		// if user does not have permissions, should just return false.
		runtime_assert( binary_file.size() > 0 );
		bool success = write_tensor_to_file_without_json( binary_file, tensor );
		if ( success ) std::cout << "Created a binary version of your tensor for faster i/o next time: " << binary_file << std::endl;
	}
}

template< class T, numeric::Size N >
void
read_tensor_from_file( std::string const & filename,
	MathNTensor< T, N > & tensor ) {
	utility::json_spirit::mObject json;
	read_tensor_from_file( filename, tensor, json );
}

////////////////////////////////////////////////////////////////////////////////
template< class T, numeric::Size N >
bool
write_tensor_to_file_without_json( std::string const & filename,
	MathNTensor< T, N > const & tensor )
{
	utility::io::ozstream bin_out( filename, std::ios::out | std::ios::binary );
	if ( !bin_out.good() ) return false;
	bin_out.write( (const char*) & tensor.data( 0 ), tensor.size()*sizeof( T ) );
	bin_out.close();
	return true;
}

// @brief writes tensor to filename and an accompanying .json file with user-supplied json
template< class T, numeric::Size N >
bool
write_tensor_to_file( std::string const & filename,
	MathNTensor< T, N > const & tensor,
	utility::json_spirit::Value const & json_input )
{
	using namespace utility;
	using namespace utility::json_spirit;
	bool success = write_tensor_to_file_without_json( filename, tensor );
	if ( !success ) return false;

	// also output json information that goes with binary file.
	if ( filename.find( ".bin" ) == std::string::npos ) utility_exit_with_message( "Filename for tensor output should be *.bin or *.bin.gz" );
	std::string const json_file( replace_in( replace_in( filename, ".bin", ".json" ), ".json.gz", ".json" ) );
	utility::io::ozstream json_out( json_file, std::ios::out );
	if ( !json_out.good() ) return false;

	// originally was going to plop in n_bins & type but these json_spirit classes are hard to manipulate.
	Value json = json_input;

	json_out << write( json ) << std::endl;
	json_out.close();
	return true;
}

// @brief writes tensor to filename and an accompanying .json file with n_bins & type.
template< class T, numeric::Size N >
bool
write_tensor_to_file( std::string const & filename,
	MathNTensor< T, N > const & tensor )
{
	using namespace utility::json_spirit;
	using namespace utility::tools;
	std::vector< Value > n_bins;
	for ( auto const & v : tensor.n_bins() ) n_bins.push_back( Value(boost::uint64_t(v)) );

	std::string type;
	if ( std::is_same<T,double>::value )  type = "double";
	else if ( std::is_same<T,Size>::value ) type = "uint64";
	if ( type == "" ) utility_exit_with_message( "Did not have a string to go with the type of tensor, like double or unit64..." );

	return write_tensor_to_file( filename, tensor,
		make_vector( Pair( "n_bins", n_bins ), Pair( "type", type ) ));
}


}

#endif /* MATHNTENSOR_IO_HH_ */


