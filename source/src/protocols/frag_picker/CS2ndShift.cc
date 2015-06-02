// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/CSTalosIO.cc
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

// unit headers
#include <protocols/frag_picker/CS2ndShift.hh>
// package headers

// disulfide compatibility
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/io/raw_data/DisulfideFile.hh>

#include <basic/database/open.hh>

// project headers
#include <basic/Tracer.hh>

// utility headers
#include <core/types.hh>

#include <string>
#include <map>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer tr( "protocols.frag_picker.SecondaryShiftCalculator" );

CS2ndShift::CS2ndShift(CSTalosIO & input_data, bool use_sslimit) {
	//Change to database files!
	std::map<char,std::map<std::string,Real> > rcprev(CS2ndShift::read_adjust_table("external/SPARTA+/tab/rcprev.tab") );
	std::map<char,std::map<std::string,Real> > rcnext(CS2ndShift::read_adjust_table("external/SPARTA+/tab/rcnext.tab") );
	std::map<char,std::map<std::string,Real> > rcadj(CS2ndShift::read_adjust_table("external/SPARTA+/tab/rcadj.tab") );
	std::map<char,std::map<std::string,Real> > randcoil(CS2ndShift::read_adjust_table("external/SPARTA+/tab/randcoil.tab") );

	std::map<char,std::map<std::string,std::pair< Real, Real > > > 	sslimit(CS2ndShift::read_sslimit_table("external/SPARTA+/tab/sslimit.tab") );

	//Recalculate Values:
	std::string const sequence( also_check_fix_disulf( input_data.get_sequence() ) );

	utility::vector1< std::pair< Size, std::string > > shift_types;

	shift_types.push_back(std::make_pair(1,"N"));//N
	shift_types.push_back(std::make_pair(2,"HA"));//HA but HA3 for Gly
	shift_types.push_back(std::make_pair(3,"C"));//C
	shift_types.push_back(std::make_pair(4,"CA"));//CA
	shift_types.push_back(std::make_pair(5,"CB"));//CB but HA2 for Gly
	shift_types.push_back(std::make_pair(6,"HN"));//HN

	//utility::vector1< utility::vector1< std::pair< Size, Real > > > secondary_shifts_;//( utility::vector1< Real > ( 0, 0 ) );


	for ( Size seqpos = 0; seqpos < sequence.length(); seqpos++ ) {
		//tr << "CALC_SECONDARY " << seqpos << " ";
		utility::vector1< std::pair< Size, Real > > res_shifts;

		for ( Size st_i = 1; st_i <= shift_types.size(); st_i++ ) {

			Real offset(0);

			std::pair<Size,std::string> shift_type(shift_types[st_i]);

			std::pair< Size, Real > shift;
			//static_cast< std::string > ( "HA" );
			static std::string const HA( "HA" );
			static std::string const CB( "CB" );
			static std::string const HA2( "HA2" );
			static std::string const HA3( "HA3" );

			Real shift_value(9999);
			bool has_shift(false);

			if ((( shift_type.second == HA ) || ( shift_type.second == CB)) && (sequence[seqpos] == 'G')) {

				offset = randcoil.find(sequence[seqpos])->second.find(HA)->second
					+rcadj.find(sequence[seqpos])->second.find(HA)->second;

				if (seqpos != 0) {
					offset = offset + rcprev.find(sequence[seqpos-1])->second.find(HA)->second;
				}
				if (seqpos != sequence.length()-1) {
					offset = offset + rcnext.find(sequence[seqpos+1])->second.find(HA)->second;
				}

				if ( input_data.has_atom(seqpos+1,HA2) && input_data.has_atom(seqpos+1,HA3) ) { // && (shift_type.second != CB) ) {
					//shift.first = shift_type.first;
					//shift.second = (input_data.get_shift(seqpos+1,HA2) + input_data.get_shift(seqpos+1,HA3))/2 - offset;

					shift_value = (input_data.get_shift(seqpos+1,HA2) + input_data.get_shift(seqpos+1,HA3))/2 - offset;
					has_shift = true;
					//res_shifts.push_back( shift );

					tr.Debug << "CALC_SECONDARY_G " << (seqpos+1) << " " << sequence[seqpos] << " " << input_data.get_shift(seqpos+1, HA2) + input_data.get_shift(seqpos+1,HA3) << " " << offset << " " << shift_type.second << " " << shift_value << std::endl;
				}

			} else {

				offset = randcoil.find(sequence[seqpos])->second.find(shift_type.second)->second
					+rcadj.find(sequence[seqpos])->second.find(shift_type.second)->second;

				if (seqpos != 0) {
					offset = offset + rcprev.find(sequence[seqpos-1])->second.find(shift_type.second)->second;
				}
				if (seqpos != sequence.length()-1) {
					offset = offset + rcnext.find(sequence[seqpos+1])->second.find(shift_type.second)->second;
				}

				if ( input_data.has_atom(seqpos+1,shift_type.second) ) {
					//shift.first = shift_type.first;

					//shift.second = input_data.get_shift(seqpos+1,shift_type.second) - offset;

					shift_value = input_data.get_shift(seqpos+1,shift_type.second) - offset;
					has_shift = true;
					//res_shifts.push_back( shift );

					tr.Debug << "CALC_SECONDARY_A " << (seqpos+1) << " " << sequence[seqpos] << " " << input_data.get_shift(seqpos+1, shift_type.second) << " " << offset << " " << shift_type.second << " " << shift_value << std::endl;
				}

			}

			if (has_shift == true) {
				if (sslimit.count(sequence[seqpos]) == 1) {
					std::map<char,std::map<std::string,std::pair< Real, Real > > >::const_iterator sslimit_at_seqpos( sslimit.find(sequence[seqpos]) );
					debug_assert( sslimit_at_seqpos != sslimit.end() );
					if (sslimit_at_seqpos->second.count(shift_type.second) == 1) {
						std::map<std::string,std::pair< Real, Real > >::const_iterator shift_type_itr( sslimit_at_seqpos->second.find(shift_type.second) );
						debug_assert( shift_type_itr != sslimit_at_seqpos->second.end() );
						Real min( shift_type_itr->second.first );
						Real max( shift_type_itr->second.second );

						// If use_sslimit == false, always accept
						if ( ((shift_value >= min) && ( shift_value <= max )) || ( !use_sslimit) ) {
							shift.first = shift_type.first;
							shift.second = shift_value;//(input_data.get_shift(seqpos+1,HA2) + input_data.get_shift(seqpos+1,HA3))/2 - offset;
							res_shifts.push_back( shift );

							tr.Debug << "USING_2ND_SHIFT: " << seqpos+1 << " " << sequence[seqpos] << "_aa" << " " << shift_type.second << " " << shift.second << std::endl;
						} else {
							if (shift_value < min) {
								tr.Debug << "SHIFT OUTLIER REMOVED: " << seqpos+1 << " " << sequence[seqpos] << " " << shift.second << " " << shift_value << " Limit: " << min << std::endl;
							} else {
								tr.Debug << "SHIFT OUTLIER REMOVED: " << seqpos+1 << " " << sequence[seqpos] << " " << shift.second << " " << shift_value << " Limit: " << max << std::endl;
							}
						}
					}
				}
			}//end: if (has_shift == true)

		}
		secondary_shifts_.push_back( res_shifts );
  }
}

std::map< char, std::map<std::string,Real> >
CS2ndShift::read_adjust_table(std::string const & file_name) {

	std::map< char, std::map<std::string,Real> > file_data_map;
	utility::vector1<std::string> column_names_;


	utility::io::izstream data(basic::database::full_name(file_name.c_str()));
  tr.Info << "read CS adjustment data from " << file_name << std::endl;
	if (!data)
		utility_exit_with_message("[ERROR] Unable to open talos file: "
															+ file_name);

	std::string line;
	std::string keyword;
	std::string entry;
	std::string junk;

	bool header_done = false;
	while (!header_done) {
		getline(data, line);
		std::istringstream line_stream(line);
		line_stream >> keyword;

		if (keyword == "VARS") {
			line_stream >> junk >> junk >> entry;

			while (!line_stream.eof()) {
				if (entry == "H") {
					entry = "HN";
				}
				column_names_.push_back(entry);
				line_stream >> entry;
			}

			if (entry == "H") {
				entry = "HN";
			}
			column_names_.push_back(entry);
		}

		if (keyword == "FORMAT") {
				header_done = true;
				getline(data,line);
		}
	}

	getline(data,line);
	while (!data.eof()) {
		std::istringstream line_stream(line);

		char aa;

		line_stream >> junk >> aa;

		std::map< std::string ,Real > linemap;
		for ( Size i = 1; i <= column_names_.size(); i++ ){
			Real offset(0.0);
			line_stream >> offset;
			linemap.insert(std::make_pair(column_names_[i], offset));
		}
		file_data_map.insert(std::make_pair(aa,linemap));

		getline(data,line);
	}

	return file_data_map;
}

std::map<char,std::map<std::string,std::pair< core::Real, core::Real > > >
CS2ndShift::read_sslimit_table(std::string const & file_name) {

	std::map<char,std::map<std::string,std::pair< core::Real, core::Real > > > file_data_map;

	utility::vector1<std::string> column_names_;


	utility::io::izstream data(basic::database::full_name(file_name.c_str()));
	tr.Info << "read CS sslimit data from " << file_name << std::endl;
	if (!data)
		utility_exit_with_message("[ERROR] Unable to open talos file: "
															+ file_name);

	std::string line;
	std::string keyword;
	std::string entry;
	std::string junk;

	std::string sub_entry;

	bool header_done = false;
	while (!header_done) {
		getline(data, line);
		std::istringstream line_stream(line);
		line_stream >> keyword;

		if (keyword == "VARS") {
			line_stream >> junk >> junk >> entry >> junk;

			while (!line_stream.eof()) {
				Size mid = entry.find_first_of('_');
				sub_entry = entry.substr(0,mid);
				column_names_.push_back(sub_entry);

				line_stream >> entry >> junk;
			}

			Size mid = entry.find_first_of('_');
			sub_entry = entry.substr(0,mid);
			column_names_.push_back(sub_entry);
		}

		if (keyword == "FORMAT") {
				header_done = true;
				getline(data,line);
		}
	}

	getline(data,line);
	while (!data.eof()) {
		std::istringstream line_stream(line);

		char aa;

		line_stream >> junk >> aa;

		std::map< std::string , std::pair< Real, Real> > linemap;
		for ( Size i = 1; i <= column_names_.size(); i++ ){
			Real min(0.0), max(0.0);
			line_stream >> min >> max;

			if ( (min < 1000) && (max < 1000) ) {
				tr.Debug << "READ_SSLIST " << column_names_[i] << " " << aa << " " << min << " " << max << std::endl;
				linemap.insert(std::make_pair(column_names_[i], std::make_pair(min,max)));
			}
		}
		file_data_map.insert(std::make_pair(aa,linemap));

		getline(data,line);
	}

	return file_data_map;
}

//TODO: move this with rest of disulfide code
std::string
CS2ndShift::also_check_fix_disulf( std::string instring ) {
	if (option[in::fix_disulf].user()) {

		core::io::raw_data::DisulfideFile ds_file( option[ in::fix_disulf ]() );

		utility::vector1< std::pair<Size,Size> > disulfides_in_file;

		ds_file.disulfides(disulfides_in_file);

		for ( Size i = 1; i <= disulfides_in_file.size(); ++i ) {

			Size l = disulfides_in_file[i].first;
			Size u = disulfides_in_file[i].second;

			if ( u <= l ) {
				utility_exit_with_message("[ERROR] Disulfide File Format: res2 must be > res1");
			}

			if ( !((instring[l-1] == 'C') || (instring[l-1] == 'c'))
					 || !((instring[u-1] == 'C') || (instring[u-1] == 'c')) ) {
				utility_exit_with_message("[ERROR] -fix_disulf residues do not map to cysteines in talos file");
			}

			instring[l-1] = 'c';
			instring[u-1] = 'c';
		}
	}

	return instring;
}

} // frag_picker
} // protocols
