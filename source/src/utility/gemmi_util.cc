// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/gemmi_util.cc
/// @brief  Utilities for working with Gemmi CIF file data
///
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <utility/gemmi_util.hh>

#include <utility/json_utilities.hh>
#include <utility/exit.hh>
#include <utility/assert.hh>
#include <utility/stream_util.hh>

#include <numeric/types.hh>

#include <gemmi/cif.hpp>
#include <boost/endian/conversion.hpp>

#include <json.hpp>

#include <iostream>

using nlohmann::json;

namespace utility {

std::string
as_cif_value(json const & obj) {
	if ( obj.is_null() ) { return "?"; }
	if ( obj.is_boolean() ) { return obj.get<bool>() ? "YES" : "NO"; }
	if ( obj.is_string() ) { return gemmi::cif::quote(obj.get<std::string>()); }
	if ( obj.is_number_float() ) { return std::to_string(obj.get<double>()); }
	if ( obj.is_number_integer() ) { return std::to_string(obj.get<ssize_t>()); }
	if ( obj.is_array() ) {
		std::string s = obj[0].get<std::string>();
		for ( numeric::Size ii(1); ii < obj.size(); ++ii ) {
			s += ' ';
			s += obj[ii].get<std::string>();
		}
		return gemmi::cif::quote(s);
	}
	utility_exit_with_message("In ready mmJSON, cannot covert entry `" + obj.dump() + "` to cif-style entry.");
}

/// @brief Read an mmJSON file into a Gemmi Document object.
/// @details Implementation based on Gemmi's read_mmjson_insitu()/fill_document_from_sajson()
/// but reimplemented to keep the Gemmi dependency header-only.
gemmi::cif::Document
gemmi_load_mmjson(std::string const & contents_of_file, std::string const & filename) {
	gemmi::cif::Document doc;
	doc.source = filename;

	json obj = json::parse(contents_of_file);

	if ( ! obj.is_object() ) { utility_exit_with_message("mmJSON file " + filename + " is not a JSON (dict) object at top level." ); }
	if ( obj.empty() ) { utility_exit_with_message("mmJSON file " + filename + " does not contain valid contents" ); }
	// We don't have tracers at utility level ...
	//if ( obj.size() != 1 ) { TR.Warning << "mmJSON file " << filename << " contains more than one structure -- only loading the first." << std::endl; }

	std::string const & block_name = obj.begin().key();
	if ( block_name.substr(0,5) != "data_" ) { utility_exit_with_message("mmJSON file " + filename + " entry " + block_name + " does not start with `data_`, as required."); }
	doc.blocks.emplace_back( block_name.substr(5) );
	std::vector<gemmi::cif::Item>& items = doc.blocks[0].items;

	auto & top = obj.begin().value();
	if ( ! top.is_object() ) { utility_exit_with_message("mmJSON file " + filename + " entry " + block_name + " is not a JSON (dict) object." ); }

	for ( auto it = top.begin(); it != top.end(); ++it) {
		std::string category_name = "_" + it.key() + ".";
		auto const & category = it.value();
		if ( !category.is_object() || category.empty() ) {
			utility_exit_with_message("mmJSON file " + filename + " entry " + it.key() + " is malformed." );
		}
		numeric::Size cif_cols = category.size();
		numeric::Size cif_rows = category.begin().value().size();
		if ( cif_rows > 1 ) {
			items.emplace_back(gemmi::cif::LoopArg{});
			gemmi::cif::Loop & loop = items.back().loop;
			loop.tags.reserve(cif_cols);
			loop.values.resize(cif_cols * cif_rows);
		}
		numeric::Size jj = 0;
		for ( auto it2 = category.begin(); it2 != category.end(); ++it2, ++jj ) {
			std::string tag = category_name + it2.key();
			auto const & arr = it2.value();
			if ( !arr.is_array() ) { utility_exit_with_message("mmJSON file " + filename + " entry " + tag + " is malformed"); }
			if ( arr.size() != cif_rows ) { utility_exit_with_message("mmJSON file " + filename + " entry " + tag + " does not have the expected number of entries."); }
			if ( cif_rows == 1 ) {
				items.emplace_back(tag, as_cif_value(arr[0]));
			} else if ( cif_rows > 1 ) {
				gemmi::cif::Loop & loop = items.back().loop;
				loop.tags.emplace_back(tag);
				for ( numeric::Size kk = 0; kk < cif_rows; ++kk ) {
					loop.values[jj + kk*cif_cols] = as_cif_value(arr[kk]);
				}
			}
		}
	}

	return doc;
}

// Forward declare
json bcif_decode( json data, json const & encodings);

json
convert_by_source_type( numeric::Real value, numeric::Size src_type ) {
	// This is primarily signed/unsigned/float conversions
	if ( 1 <= src_type && src_type <= 3 ) { // signed int
		return ssize_t(value);
	} else if ( 4 <= src_type && src_type <= 6 ) { // unsigned int
		debug_assert( value >= 0 );
		return size_t(value);
	} else { // Float
		return value;
	}
}

// Utility function for bcif_decode
json
bcif_vector_decode( json const & data, json const & encoding ) {
	debug_assert( encoding.is_object() );

	std::string kind;
	if ( encoding.count("kind") == 0 ) {
		std::cout << encoding.dump() << std::endl;
	}
	extract_nonempty_string_from_json( encoding, "kind", kind );

	if ( kind == "FixedPoint" ) {
		// Since we're smashing everything through std:string eventually, the 32/64 bit difference isn't really needed.
		std::string srcType = encoding["srcType"].get<std::string>();
		if ( srcType != "Float32" && srcType != "Float64"  ) {
			utility_exit_with_message("Can't interpret " + srcType + " as FixedPoint source type.");
		}
		numeric::Real factor = encoding["factor"].get<numeric::Real>();
		std::vector< numeric::Real > vec;
		for ( auto & v: data ) {
			vec.push_back( v.get<numeric::Real>() / factor ); // Always real
		}
		return vec;
	} else if ( kind == "IntervalQuantization" ) {
		std::string srcType = encoding["srcType"].get<std::string>();
		if ( srcType != "Float32" && srcType != "Float64" ) {
			utility_exit_with_message("Can't interpret " + srcType + " as IntervalQuantization source type.");
		}
		numeric::Real min = encoding["min"].get<numeric::Real>();
		numeric::Real max = encoding["max"].get<numeric::Real>();
		numeric::Real numSteps = encoding["numSteps"].get<numeric::Real>(); // includes both fenceposts
		numeric::Real step = (max - min)/(numSteps-1);

		std::vector< numeric::Real > vec;
		for ( auto & v: data ) {
			vec.push_back( min + v.get<numeric::Size>() * step ); // always real
		}
		return vec;
	} else if ( kind == "RunLength" ) {
		numeric::Size src_type = encoding["srcType"].get<numeric::Size>();
		json vec = json::array();
		for ( numeric::Size ii(0); ii < data.size()-1; ii += 2 ) {
			std::int32_t v = data[ii].get< std::int32_t >();
			for ( numeric::Size jj(1); jj <= data[ii+1].get< numeric::Size >(); ++jj ) {
				vec.push_back( convert_by_source_type(v, src_type) );
			}
		}
		if ( vec.size() != encoding["srcSize"].get< numeric::Size >() ) {
			utility_exit_with_message("In RunLength encoding, decoded size of " + std::to_string(vec.size()) + " doesn't match source size of " + encoding["srcSize"].dump() );
		}
		return vec;
	} else if ( kind == "Delta" ) {
		numeric::Size src_type = encoding["srcType"].get<numeric::Size>();
		std::int32_t current = encoding["origin"].get< std::int32_t >();
		json vec = json::array();
		for ( auto & v: data ) {
			current += v.get< std::int32_t >();
			vec.push_back( convert_by_source_type(current, src_type) );
		}
		return vec;
	} else if ( kind == "IntegerPacking" ) {
		if ( encoding["isUnsigned"].get< bool >() ) {
			numeric::Size nbytes = encoding["byteCount"].get< numeric::Size >();
			std::uint32_t max_val;
			if ( nbytes == 1 ) {
				max_val = std::numeric_limits< std::uint8_t >::max();
			} else if ( nbytes == 2 ) {
				max_val = std::numeric_limits< std::uint16_t >::max();
			} else {
				utility_exit_with_message("Can't understand byteCount of " + std::to_string(nbytes) + " in IntegerPacking encoding");
			}
			bool was_max = false;
			std::vector< std::uint32_t > vec;
			for ( auto & v: data ) {
				std::uint32_t val = v.get< std::uint32_t >(); // Won't take the full range, but no harm in treating as such.
				if ( was_max ) {
					vec[ vec.size()-1 ] += val;
				} else {
					vec.push_back( val );
				}
				was_max = (val == max_val);
			}
			if ( vec.size() != encoding["srcSize"].get< numeric::Size >() ) {
				utility_exit_with_message("In IntegerPacking encoding, decoded size of " + std::to_string(vec.size()) + " doesn't match source size of " + encoding["srcSize"].dump() );
			}
			return vec;
		} else { // Signed IntegerPacking
			numeric::Size nbytes = encoding["byteCount"].get< numeric::Size >();
			std::int32_t max_val, min_val;
			if ( nbytes == 1 ) {
				max_val = std::numeric_limits< std::int8_t >::max();
				min_val = std::numeric_limits< std::int8_t >::min();
			} else if ( nbytes == 2 ) {
				max_val = std::numeric_limits< std::int16_t >::max();
				min_val = std::numeric_limits< std::int16_t >::min();
			} else {
				utility_exit_with_message("Can't understand byteCount of " + std::to_string(nbytes) + " in IntegerPacking encoding");
			}
			bool was_extreme = false;
			std::vector< std::int32_t > vec;
			for ( auto & v: data ) {
				std::int32_t val = v.get< std::int32_t >();
				if ( was_extreme ) {
					vec[ vec.size()-1 ] += val;
				} else {
					vec.push_back( val );
				}
				was_extreme = (val == max_val) || (val == min_val);
			}
			if ( vec.size() != encoding["srcSize"].get< numeric::Size >() ) {
				utility_exit_with_message("In IntegerPacking encoding, decoded size of " + std::to_string(vec.size()) + " doesn't match source size of " + encoding["srcSize"].dump() );
			}
			return vec;
		}
	} else {
		utility_exit_with_message("Cannot interpret BCIF encoding kind " + kind);
	}

	utility_exit_with_message("Error interpreting encoding");
}

json
bcif_byte_array_decode( json const & data, json const & encoding) {
	std::vector< unsigned char > byte_array;
	if ( data.is_string() ) {
		std::string string_data = data.get< std::string >();
		byte_array.assign( string_data.begin(), string_data.end() );
	} else if ( data.is_array() ) {
		for ( auto & v: data ) {
			byte_array.push_back( v.get< unsigned char >() );
		}
	} else if ( data.is_binary() ) {
		byte_array = data.get_binary();
	} else {
		utility_exit_with_message(std::string("Can't interpret ByteArray data of type ") + data.type_name() );
	}
	numeric::Size type = encoding["type"].get<numeric::Size>();
	// All values are little-ending coded.
	switch ( type ) {
	case 1: {// INT8
		std::vector< std::int8_t > vec;
		for ( numeric::Size offset(0); offset < byte_array.size(); offset += 1 ) {
			vec.push_back( boost::endian::load_little_s16( byte_array.data() + offset ) );
		}
		return vec;
	}
	case 2: {// INT16
		std::vector< std::int16_t > vec;
		for ( numeric::Size offset(0); offset < byte_array.size(); offset += 2 ) {
			vec.push_back( boost::endian::load_little_s16( byte_array.data() + offset ) );
		}
		return vec;
	}
	case 3: {// INT32
		std::vector< std::int32_t > vec;
		for ( numeric::Size offset(0); offset < byte_array.size(); offset += 4 ) {
			vec.push_back( boost::endian::load_little_s32( byte_array.data() + offset ) );
		}
		return vec;
	}
	case 4: {// UINT8
		std::vector< std::uint8_t > vec;
		for ( numeric::Size offset(0); offset < byte_array.size(); offset += 1 ) {
			vec.push_back( boost::endian::load_little_u16( byte_array.data() + offset ) );
		}
		return vec;
	}
	case 5: {// UNIT16
		std::vector< std::uint16_t > vec;
		for ( numeric::Size offset(0); offset < byte_array.size(); offset += 2 ) {
			vec.push_back( boost::endian::load_little_u16( byte_array.data() + offset ) );
		}
		return vec;
	}
	case 6: {// UNIT32
		std::vector< std::uint32_t > vec;
		for ( numeric::Size offset(0); offset < byte_array.size(); offset += 4 ) {
			vec.push_back( boost::endian::load_little_u32( byte_array.data() + offset ) );
		}
		return vec;
	}
	case 32: {// FLOAT32
		std::vector< float > vec;
		for ( numeric::Size offset(0); offset < byte_array.size(); offset += 4 ) {
			vec.push_back( boost::endian::endian_load<float, 4, boost::endian::order::little>( byte_array.data() + offset ) );
		}
		return vec;
	}
	case 33: {// FLOAT64
		std::vector< double > vec;
		for ( numeric::Size offset(0); offset < byte_array.size(); offset += 8 ) {
			vec.push_back( boost::endian::endian_load<double, 8, boost::endian::order::little>( byte_array.data() + offset ) );
		}
		return vec;
	}
	default:
		utility_exit_with_message("Can't interpret " + std::to_string(type) + " as BCIF encoding type");
	}
}

// Utility function for bcif_decode() -- handle StringArray case.
json
bcif_string_array_decode( json const & data, json const & encoding) {
	json decoded_data = bcif_decode( data, encoding["dataEncoding"] );
	json offsets = bcif_decode( encoding["offsets"], encoding["offsetEncoding"] );

	std::string stringData = encoding["stringData"].get< std::string >();


	std::vector< int > indexes; // int, as sometimes indices can be negative.
	if ( decoded_data.is_array() ) {
		for ( auto & v: decoded_data ) {
			indexes.push_back( v.get< int >() );
		}
	} else {
		utility_exit_with_message(std::string("Can't interpret StringArray data of type ") + data.type_name() );
	}

	std::vector< std::string > vec;
	for ( auto sindex: indexes ) {
		if ( sindex < 0 ) {
			vec.push_back("?"); // Theoretically will be masked later.
			continue;
		}
		numeric::Size index = sindex;
		if ( index+1 >= offsets.size() ) {
			std::cerr << "Bad index " << index << " in StringArray decoding. Maximum offset " << offsets.size() << std::endl;
			vec.push_back("?");
			continue;
		}
		numeric::Size begin = offsets[index].get< numeric::Size >();
		numeric::Size end = offsets[index+1].get< numeric::Size >();
		vec.push_back( stringData.substr( begin, end-begin ) );
	}
	return vec;
}

// Utility function for bcif_decode_data_obj() - decode a BCIF data/encoding pair, where the encoding pair can be an array.
// @details data is taken by copy, as we may pass it through a number of decoding steps.
json
bcif_decode( json data, json const & encodings) {
	if ( data.is_binary() ) {
		std::vector< std::uint8_t > binary_array = data.get_binary();
		data = binary_array; // Need this indirection to force proper type conversion.
	}

	debug_assert( encodings.is_array() );
	numeric::Size nencodings = encodings.size();

	for ( numeric::Size ee(1); ee <= nencodings; ++ee ) {
		auto const & encoding = encodings[ nencodings - ee ]; // Decode from the end backwards
		std::string encoding_kind;
		extract_nonempty_string_from_json( encoding, "kind", encoding_kind);

		// We split out the Byte & String arrays, as they have differences with the input data types.
		if ( encoding_kind == "ByteArray" ) {
			data = bcif_byte_array_decode( data, encoding );
		} else if ( encoding_kind == "StringArray" ) {
			data = bcif_string_array_decode( data, encoding );
		} else {
			data = bcif_vector_decode( data, encoding );
		}
	}

	return data;
}

// Utility function for gemmi_load_bcif() - decode a BCIF Data object.
std::vector< std::string >
bcif_decode_data_obj( json const & data_obj ) {
	debug_assert( data_obj.is_object() );
	debug_assert( data_obj.count("data") == 1 );

	json data = data_obj["data"];
	json encoding;
	extract_nonempty_array_from_json( data_obj, "encoding", encoding);

	json decoded_data = bcif_decode( data, encoding );

	std::vector< std::string > decoded;
	for ( auto & entry: decoded_data ) {
		decoded.push_back( as_cif_value(entry) );
	}

	return decoded;
}

gemmi::cif::Document
gemmi_load_bcif(std::string const & contents_of_file, std::string const & filename) {
	gemmi::cif::Document doc;
	doc.source = filename;

	json file = json::from_msgpack(contents_of_file);
	if ( ! file.is_object() ) { utility_exit_with_message("BCIF file " + filename + " is not a proper (dict) object at top level." ); }
	if ( file.count("dataBlocks") == 0 ) { utility_exit_with_message("BCIF file " + filename + " is missing the dataBlocks entry." ); }
	if ( ! file["dataBlocks"].is_array() || file["dataBlocks"].empty() ) { utility_exit_with_message("BCIF file " + filename + " does not have a properly formatted dataBlocks entry." ); }
	// We don't have tracers at utility level ...
	//if ( file["dataBlocks"].size() != 1 ) { TR.Warning << "BCIF file " << filename << " contains more than one structure -- only loading the first." << std::endl; }

	auto & dataBlock = file["dataBlocks"][0];
	doc.blocks.emplace_back( dataBlock["header"] );
	std::vector<gemmi::cif::Item>& items = doc.blocks[0].items;

	for ( auto & category: dataBlock["categories"] ) {
		std::string category_name = category["name"].get<std::string>() + ".";

		numeric::Size cif_cols = category["columns"].size();
		numeric::Size cif_rows = category["rowCount"];

		if ( cif_rows > 1 ) {
			items.emplace_back(gemmi::cif::LoopArg{});
			gemmi::cif::Loop & loop = items.back().loop;
			loop.tags.reserve(cif_cols);
			loop.values.resize(cif_cols * cif_rows);
		}

		for ( numeric::Size jj(0); jj < cif_cols; ++jj ) {
			auto & column = category["columns"][jj];
			std::string tag = category_name + column["name"].get<std::string>();

			std::vector< std::string > decoded = bcif_decode_data_obj( column["data"] );
			if ( decoded.size() != cif_rows ) {
				utility_exit_with_message("In BCIF file " + filename + " entry " + tag + " the number of decoded entries (" + std::to_string(decoded.size()) + ") does not match that expected (" + std::to_string(cif_rows) + ")");
			}
			if ( column.count("mask") && !column["mask"].is_null() ) {
					std::vector< std::string > mask = bcif_decode_data_obj( column["mask"] );
					if ( mask.size() > decoded.size() ) {
						utility_exit_with_message("In BCIF file " + filename + " entry " + tag + " the mask (" + std::to_string(mask.size()) + ") is larger than the data (" + std::to_string(decoded.size()) + ")" );
					}
					for ( numeric::Size mm(0); mm < mask.size(); ++mm ) {
						numeric::Size mask_val = std::stoi( mask[mm] );
						if ( mask_val == 0 ) {
							// do nothing
						} else if ( mask_val == 1 ) {
							decoded[mm] = ".";
						} else if ( mask_val == 2 ) {
							decoded[mm] = "?";
						} else {
							utility_exit_with_message("Can't interpret BCIF mask value of `" + std::to_string(mask_val) + "` in entry " + tag);
						}
					}
			}

			if ( cif_rows == 1 ) {
				items.emplace_back(tag, decoded[0]);
			} else if ( cif_rows > 1 ) {
				gemmi::cif::Loop & loop = items.back().loop;
				loop.tags.emplace_back(tag);
				for ( numeric::Size kk = 0; kk < cif_rows; ++kk ) {
					loop.values[jj + kk*cif_cols] = decoded[kk];
				}
			}
		}
	}

	return doc;
}


} // namespace utility
