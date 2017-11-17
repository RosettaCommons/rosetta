// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/antibody/grafting/regex_based_cdr_detection.hh
/// @brief RegEx based detection of CRS's
/// @author Sergey Lyskov
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)
/// @author Jeliazko Jeliazkov

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/regex_based_cdr_detection.hh>
#include <protocols/antibody/grafting/regex_manager.hh>
#include <protocols/antibody/grafting/exception.hh>

#include <basic/Tracer.hh>

#include <regex>
#include <string>


namespace protocols {
namespace antibody {
namespace grafting {

using std::string;

static basic::Tracer TR("protocols.antibody.grafting");

void RegEx_based_CDR_Detector::detect(AntibodySequence &A)
{
	*this << "RegEx_based_CDR_Detector run with arguments:\n  heavy: " + A.heavy.sequence + "\n  light: " + A.light.sequence + "\n\n";
	set("heavy", A.heavy.sequence);  set("light", A.light.sequence);

	detect_heavy_chain(A);
	detect_light_chain(A);

	*this << "RegEx_based_CDR_Detector results:\n";

	struct {
		string name;
		CDR_Bounds &cdr;
		string sequence;
	} J[] {
		{"h1", A.heavy.cdr1, A.h1_sequence()}, {"h2", A.heavy.cdr2, A.h2_sequence()}, {"h3", A.heavy.cdr3, A.h3_sequence()},
	    {"l1", A.light.cdr1, A.l1_sequence()}, {"l2", A.light.cdr2, A.l2_sequence()}, {"l3", A.light.cdr3, A.l3_sequence()},
	};

	for(auto &j : J) {
		string line = j.name + ": " + j.sequence;
		line += line.size() < 8 ? "\t\t\t" : (line.size() < 16 ? "\t\t" : "\t");
		*this << "\t" + line << j.cdr << "\n";
		utility::json_spirit::Object cdr;

		cdr.push_back( utility::json_spirit::Pair("sequence", j.sequence) );
		cdr.push_back( utility::json_spirit::Pair("begin",    int(j.cdr.begin) ) );
		cdr.push_back( utility::json_spirit::Pair("end",      int(j.cdr.end)) );
		set(j.name, cdr);
	}
	*this << "\n";
}


void RegEx_based_CDR_Detector::detect_heavy_chain(AntibodySequence &A)
{
	if( ! antibody_grafting_usable() ) {
		utility_exit_with_message("ERROR: Your compiler does not have full support for C++11 regex, and therefore can't support RegEx_based_CDR_Detector/antibody grafting.");
	}

	AntibodyChain & H(A.heavy);
	string const & heavy_chain_sequence(H.sequence);


	uint &H1_begin(H.cdr1.begin), &H1_end(H.cdr1.end), &H2_begin(H.cdr2.begin), &H2_end(H.cdr2.end), &H3_begin(H.cdr3.begin), &H3_end(H.cdr3.end);

	bool H1_detected(false), H3_detected(false);

	string heavy_first = ( heavy_chain_sequence.size() > 140) ? heavy_chain_sequence.substr(0, 70) : heavy_chain_sequence.substr(0, 60);

	// Identifying H1
	uint fr1_begin = 0;
	{
		string H1_pattern( RegExManager::get_instance()->H1_pattern() );
		std::regex re( H1_pattern, std::regex::extended ); 	std::smatch m;
		if( std::regex_search( heavy_first, m, re ) ) {
			string H1 = m.str().substr( 4, m.str().size() - 8 );  // H1 = res.group()[4:-4]

			H1_begin = heavy_chain_sequence.find(H1);
			H1_end   = H1_begin + H1.size();

			H1_detected = true;

			TR << "H1 detected: " << H1 << " (" << H1.size() << " residues at positions: " << H1_begin << " to " << H1_end << ")" << std::endl;

			if( H1_begin > 25 ) { fr1_begin = H1_begin - 25; }

		} else {
			std::cerr << "H1 detection: failed" << std::endl;
			throw CREATE_EXCEPTION(_AE_cdr_detection_failed_, "H1 detection: failed. Pattern:"+ H1_pattern + " query:" + heavy_chain_sequence);
		}
	}
	string heavy_second = ( heavy_chain_sequence.size() > 140 ) ? heavy_chain_sequence.substr(H1_end+33+15-1, 95+fr1_begin-1) : heavy_chain_sequence.substr(H1_end+33+15-1);

	// Identifying H3
	{
		string H3_pattern( RegExManager::get_instance()->H3_pattern() );
		std::regex re( H3_pattern, std::regex::extended ); 	std::smatch m;
		if( std::regex_search( heavy_second, m, re ) ) {
			string H3 = m.str().substr( 3, m.str().size() - 7 );  // H3 = res.group()[3:-4]

			H3_begin = heavy_chain_sequence.find(H3);
			H3_end   = H3_begin + H3.size();

			H3_detected = true;

			TR << "H3 detected: " << H3 << " (" << H3.size() << " residues at positions: " << H3_begin << " to " << H3_end << ")" << std::endl;

		} else {
			std::cerr << "H3 detection: failed" << std::endl;
			throw CREATE_EXCEPTION(_AE_cdr_detection_failed_, "H3 detection: failed. Pattern:"+ H3_pattern + " query:" + heavy_chain_sequence);
		}
	}

	if( H1_detected and H3_detected ) {
		H2_begin = H1_end + 15 - 1;
		H2_end   = H3_begin - 33 + 1;

		string H2( heavy_chain_sequence.substr(H2_begin, H2_end-H2_begin) );

		TR << "H2 detected: " << H2 << " (" << H2.size() << " residues at positions: " << H2_begin << " to " << H2_end << ")" << std::endl;

	} else {
		std::cerr << "H2 detection: failed" << std::endl;
		throw CREATE_EXCEPTION(_AE_cdr_detection_failed_, "H2 detection is based on the detection and position of H1 and H3!");
	}
}


void RegEx_based_CDR_Detector::detect_light_chain(AntibodySequence &A)
{
	if( ! antibody_grafting_usable() ) {
		utility_exit_with_message("ERROR: Your compiler does not have full support for C++11 regex, and therefore can't support RegEx_based_CDR_Detector/antibody grafting.");
	}

	AntibodyChain & L(A.light);
	string const & light_chain_sequence(L.sequence);


	uint &L1_begin(L.cdr1.begin), &L1_end(L.cdr1.end), &L2_begin(L.cdr2.begin), &L2_end(L.cdr2.end), &L3_begin(L.cdr3.begin), &L3_end(L.cdr3.end);

	bool L1_detected(false), L3_detected(false);

	// it's not clear to me why this is necessary ~BDW.  TODO: clarify origin
	std::string light_first = ( light_chain_sequence.size() > 130 )  ?  light_chain_sequence.substr( 0, 65 )  :  light_chain_sequence.substr( 0, 60 );

	try {

		// Identifying L1
		{
			string L1_pattern( RegExManager::get_instance()->L1_pattern() );
			std::regex re( L1_pattern, std::regex::extended );
			std::smatch m;

			if( std::regex_search( light_first, m, re ) ) {
				string L1 = m.str().substr( 1, m.str().size() - 4 );  // L1 = res.group()[1:-3]

				L1_begin = light_chain_sequence.find( L1 );
				L1_end   = L1_begin + L1.size();

				L1_detected = true;

				TR << "L1 detected: " << L1 << " (" << L1.size() << " residues at positions: " << L1_begin << " to " << L1_end << ")" << std::endl;

			} else {
				std::cerr << "L1 detected: false" << std::endl;
				throw CREATE_EXCEPTION(_AE_cdr_detection_failed_, "L1 detection: failed. Pattern:" + L1_pattern + " query:" + light_chain_sequence);
			}
		}

		string light_second =( light_chain_sequence.size() > 130 )  ?  light_chain_sequence.substr(L1_end+16+7-1, L1_end+16+7+80 - L1_end+16+7)  : light_chain_sequence.substr(L1_end+16+7-1);

		// Identifying L3
		{
			string L3_pattern( RegExManager::get_instance()->L3_pattern() );
			std::regex re( L3_pattern, std::regex::extended );
			std::smatch m;

			if( std::regex_search( light_second, m, re ) ) {
				string L3 = m.str().substr( 1, m.str().size() - 5 ); // L3 = res.group()[1:-4]

				L3_begin = light_chain_sequence.find( L3 );
				L3_end   = L3_begin + L3.size();

				L3_detected = true;

				TR << "L3 detected: " << L3 << " (" << L3.size() << " residues at positions: " << L3_begin << " to " << L3_end << ")" << std::endl;

			} else {
				std::cerr << "L3 detected: false" << std::endl;
				throw CREATE_EXCEPTION(_AE_cdr_detection_failed_, "L3 detection: failed. Pattern:"+ L3_pattern + " query:" + light_chain_sequence);
			}
		}

		// Identifying L2
		if( L1_detected and L3_detected ) {
			L2_begin = L1_end + 16 - 1;
			L2_end   = L2_begin + 7;

			string L2 = light_chain_sequence.substr(L2_begin, L2_end-L2_begin);  //L2 = light_chain[L2_begin:L2_begin+7]  # L2 is identified here. Current implementation can deal with only 7-resiue L2
			TR << "L2 detected: " << L2 << " (" << L2.size() << " residues at positions: " << L2_begin << " to " << L2_end << ")" << std::endl;

		} else {
			TR.Error << "L2 detection: failed" << std::endl; exit(1);
		}

	} catch( std::regex_error & e ) {
		std::cerr << e.what() << " " << e.code() << std::endl;
		if( e.code() == std::regex_constants::error_backref ) std::cerr << "error_backref" << std::endl;
		if( e.code() == std::regex_constants::error_brack ) std::cerr << "error_brack" << std::endl;
		exit(1);
	}

}

} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__
