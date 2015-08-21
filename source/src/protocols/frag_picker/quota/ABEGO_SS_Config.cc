// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/quota/ABEGO_SS_Config.hh
/// @brief Reads a config file for quota selector
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/quota/ABEGO_SS_Config.hh>
#include <protocols/frag_picker/quota/ABEGO_SS_Map.hh>

// utility headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>


// boost
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace quota {

static thread_local basic::Tracer trABEGO_SS_Config(
	"protocols.frag_picker.quota.ABEGO_SS_Config");

ABEGO_SS_Config::ABEGO_SS_Config(std::string & file_name) : source_file_name_(file_name) {

	utility::io::izstream data(file_name.c_str());
	trABEGO_SS_Config.Info << "reading SS-ABEGO quota config from " << file_name << std::endl;
	if ( !data ) {
		utility_exit_with_message("[ERROR] Unable to open quota config file: "
			+ file_name);
	}
	std::string line;
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sep(", :");
	while ( data ) {
		getline(data, line);
		if ( line.length() < 4 ) continue;
		if ( line[0]=='#' ) continue;
		tokenizer tokens(line, sep);
		utility::vector1<std::string> t;
		for ( tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter ) {
			t.push_back(*tok_iter);
		}
		trABEGO_SS_Config.Trace << "Parsing a line: "<<line<<"\n\t("<<t.size()<<" tokens)"<<std::endl;
		if ( line.find(':')!=line.npos ) {
			//Size pool_id;
			utility::vector1< std::pair<Size,Size> > bins;
			std::string pool_name;
			//pool_id = boost::lexical_cast<int>(t[1]);  set but never used ~Labonte
			pool_name = t[2];
			trABEGO_SS_Config.Trace << "a new pool defined: >"<<pool_name<<"<";
			for ( Size i=3; i<=t.size(); i+=2 ) {
				bins.push_back( std::pair<Size,Size>(ss_index(t[i][0]),abego_index(t[i+1][0])) );
				trABEGO_SS_Config.Trace << " "<<bins[bins.size()].first<<" : "<<bins[bins.size()].second;
			}
			trABEGO_SS_Config.Trace << std::endl;
			pool_defs_.push_back( bins );
			pool_names_.push_back( pool_name );
		} else {
			utility::vector1<Real> row;
			for ( Size i=5; i<=t.size(); i++ ) {
				row.push_back( boost::lexical_cast<double>(t[i]) );
			}
			bin_probs_.push_back(row);
		}
	}
	data.close();

	trABEGO_SS_Config.Debug << "Loded data from "<<source_file_name_<<
		" columns: "<<n_columns()<<" rows: " <<n_rows()<<std::endl;

	assert( pool_defs_.size() == pool_names_.size() ); // the number of pool names doesn't mach pool definitions
}

Real ABEGO_SS_Config::highest_probability(Size pos) {

	Real max = -100000.0;
	for ( Size i=1; i<=bin_probs_[pos].size(); i++ ) {
		if ( bin_probs_[pos][i] > max ) {
			max = bin_probs_[pos][i];
		}
	}
	return max;
}

Size ABEGO_SS_Config::most_probable_bin(Size pos) {

	Real max = -100000.0;
	Size i_max = 0;
	for ( Size i=1; i<=bin_probs_[pos].size(); i++ ) {
		if ( bin_probs_[pos][i] > max ) {
			max = bin_probs_[pos][i];
			i_max = i;
		}
	}

	return i_max;
}


} // quota
} // frag_picker
} // protocols

