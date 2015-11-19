// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/antibody/grafting/regex_based_cdr_detection.hh
/// @brief RegEx based detection of CRS's
/// @author Sergey Lyskov
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/regex_based_cdr_detection.hh>
#include <protocols/antibody/grafting/exception.hh>

#include <basic/Tracer.hh>


#include <regex>
#include <string>


namespace protocols {
namespace antibody {
namespace grafting {

using std::string;

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.grafting");

void RegEx_based_CDR_Detector::detect(AntibodySequence &A)
{
	detect_heavy_chain(A);
	detect_light_chain(A);
}


void RegEx_based_CDR_Detector::detect_heavy_chain(AntibodySequence &A)
{
    AntibodyChain & H(A.heavy);
	string const & heavy_chain_sequence(H.sequence);


	uint &H1_begin(H.cdr1.begin), &H1_end(H.cdr1.end), &H2_begin(H.cdr2.begin), &H2_end(H.cdr2.end), &H3_begin(H.cdr3.begin), &H3_end(H.cdr3.end);

	bool H1_detected(false), H3_detected(false);

	string H1_pattern{"C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C|G)(Q|K|H|E|L|R)"};
    string H3_pattern{"C[A-Z]{1,33}(W)(G|A|C)[A-Z]{1,2}(Q|S|G|R)"};


	// heavy_first = heavy_chain[:70] if len(heavy_chain) > 140 else heavy_chain[:60]
	string heavy_first = ( heavy_chain_sequence.size() > 140) ? heavy_chain_sequence.substr(0, 70) : heavy_chain_sequence.substr(0, 60);


    // ## H1
    // res = re.search(H1_pattern, heavy_first)
    // H1 = False
    // len_FR_H1 = 0
    // if res:
    //     H1 = res.group()[4:-4]
    //     H1_begin = heavy_chain.index(H1)
    //     H1_end = H1_begin + len(H1) - 1
    //     print "H1 detected: %s (%d residues at positions %d to %d)" % (H1, len(H1), H1_begin, H1_end)
    //     FR_H1 = heavy_chain[:H1_begin]
    //     if len(FR_H1) >  25:
    //         len_FR_H1 = len(FR_H1) - 25
    //         FR_H1 = heavy_chain[len_FR_H1:H1_begin]
    // else:
    //     print "H1 detected: False"

	// H1
	uint fr1_begin = 0;
	{
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
			throw _AE_cdr_detection_failed_("H1 detection: failed. Pattern:"+H1_pattern + " query:" + heavy_chain_sequence);
		}
	}

	// heavy_second = heavy_chain[H1_end+33+15:H1_end+33+15+95+len_FR_H1] if len(heavy_chain) > 140 else heavy_chain[H1_end+33+15:]
	string heavy_second = ( heavy_chain_sequence.size() > 140 ) ? heavy_chain_sequence.substr(H1_end+33+15-1, 95+fr1_begin-1) : heavy_chain_sequence.substr(H1_end+33+15-1);


    // ## H3
    // H3 = False
    // res = re.search(H3_pattern,heavy_second)
    // if res:
    //     H3 = res.group()[3:-4]
    //     H3_begin = heavy_chain.index(H3)
    //     H3_end = H3_begin + len(H3) - 1
    //     print "H3 detected: %s (%d residues at positions %d to %d)" % (H3, len(H3), H3_begin, H3_end)
    // else:
    //     print "H3 detected: False"


	{
		std::regex re( H3_pattern, std::regex::extended ); 	std::smatch m;
		if( std::regex_search( heavy_second, m, re ) ) {
			string H3 = m.str().substr( 3, m.str().size() - 7 );  // H3 = res.group()[3:-4]

			H3_begin = heavy_chain_sequence.find(H3);
			H3_end   = H3_begin + H3.size();

			H3_detected = true;

			TR << "H3 detected: " << H3 << " (" << H3.size() << " residues at positions: " << H3_begin << " to " << H3_end << ")" << std::endl;

		} else {
			std::cerr << "H3 detection: failed" << std::endl;
			throw _AE_cdr_detection_failed_("H3 detection: failed. Pattern:"+H3_pattern + " query:" + heavy_chain_sequence);
		}
	}

    // if H1 and H3:
    //     H2_begin = H1_end + 15
    //     H2_end = H3_begin - 33
    //     H2 = heavy_chain[H2_begin:H2_begin + H2_end-H2_begin+1]
    //     print "H2 detected: %s (%d residues at positions %d to %d)" % (H2, len(H2), H2_begin, H2_end)

	if( H1_detected and H3_detected ) {
		H2_begin = H1_end + 15 - 1;
		H2_end   = H3_begin - 33 + 1;

		string H2( heavy_chain_sequence.substr(H2_begin, H2_end-H2_begin) );

		TR << "H2 detected: " << H2 << " (" << H2.size() << " residues at positions: " << H2_begin << " to " << H2_end << ")" << std::endl;

	} else {
		std::cerr << "H2 detection: failed" << std::endl; exit(1);
	}
}


void RegEx_based_CDR_Detector::detect_light_chain(AntibodySequence &A)
{
    AntibodyChain & L(A.light);
	string const & light_chain_sequence(L.sequence);


	uint &L1_begin(L.cdr1.begin), &L1_end(L.cdr1.end), &L2_begin(L.cdr2.begin), &L2_end(L.cdr2.end), &L3_begin(L.cdr3.begin), &L3_end(L.cdr3.end);

	bool L1_detected(false), L3_detected(false);

	std::string L1_pattern{"C[A-Z]{1,17}(WYL|WLQ|WFQ|WYQ|WYH|WVQ|WVR|WWQ|WVK|WYR|WLL|WFL|WVF|WIQ|WYR|WNQ|WHL|WHQ|WYM|WYY)"};
	std::string L3_pattern{"C[A-Z]{1,15}(L|F|V|S)G[A-Z](G|Y)"};


	// it's not clear to me why this is necessary ~BDW.  TODO: clarify origin
	std::string light_first = ( light_chain_sequence.size() > 130 )  ?  light_chain_sequence.substr( 0, 65 )  :  light_chain_sequence.substr( 0, 60 );

	try {

		// Identifying L1
		/* Original Python code, remove after tests is passed:
			## L1
			len_FR_L1=0
			res = re.search(L1_pattern,light_first)
			FR_L1 = L1 = False
			if res:
				L1 = res.group()[1:-3]
				L1_begin = light_chain.index(L1)
				L1_end = L1_begin + len(L1) - 1
				print "L1 detected: %s (%d residues at positions %d to %d)" % (L1, len(L1), L1_begin, L1_end)
				FR_L1 = light_chain[:L1_begin]
				if len(FR_L1) >  24:
					len_FR_L1 = len(FR_L1) - 24
					FR_L1 = light_chain[len_FR_L1:L1_begin]
			else:
				print "L1 detected: False"
		 */
		{
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
				throw _AE_cdr_detection_failed_("L1 detection: failed. Pattern:"+L1_pattern + " query:" + light_chain_sequence);
			}
		}

		// Python: light_second = light_chain[L1_end+16+7:L1_end+16+7+80] if len(light_chain) > 130 else light_chain[L1_end+16+7:]
		string light_second =( light_chain_sequence.size() > 130 )  ?  light_chain_sequence.substr(L1_end+16+7-1, L1_end+16+7+80 - L1_end+16+7)  : light_chain_sequence.substr(L1_end+16+7-1); // light_chain[L1_end+16+7:]
		//TR << "light_second: " << light_second << std::endl;


		// Identifying L3
		/* Original Python code, remove after tests is passed:
			L3 = False
			res = re.search(L3_pattern,light_second)
			if res:
				L3 = res.group()[1:-4]
				L3_begin = light_chain.index(L3)
				L3_end = L3_begin + len(L3) - 1
				print "L3 detected: %s ( %d residues at positions %d to %d)" % (L3, len(L3), L3_begin, L3_end)
			else:
				print "L3 detected: False"
		 */
		{
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
				throw _AE_cdr_detection_failed_("L3 detection: failed. Pattern:"+L3_pattern + " query:" + light_chain_sequence);
			}
		}


		// Identifying L2
		/* Original Python code, remove after tests is passed:
		   if L1 and L3:
				L2_begin = L1_end + 16
				L2_end = L2_begin + 7 - 1

				L2 = light_chain[L2_begin:L2_begin+7]  # L2 is identified here. Current implementation can deal with only 7-resiue L2
				print "L2 detected: %s (%d residues at positions %d to %d)" % (L2, len(L2), L2_begin, L2_end)

				#FR_L1 = light_chain[:L1_begin]
				FR_L2 = light_chain[ L1_end + 1  :  L1_end + 1 + 15                    ]
				FR_L3 = light_chain[ L2_end + 1  :  L2_end + 1 + L3_begin - L2_end - 1 ]
				FR_L4 = light_chain[ L3_end + 1  :  L3_end + 1 + 12                    ]

				print "FR_L1: ", FR_L1
				print "FR_L2: ", FR_L2
				print "FR_L3: ", FR_L3
				print "FR_L4: ", FR_L4
				print "L segments: ",FR_L1,L1,FR_L2,L2,FR_L3,L3,FR_L4

				# Light chain sub-type classification. This is useful in the future. But currently this is not used.
				# ... skipped, see Google doc for details
				# FR classification by AHo. This might be useful in the future. But currently this is not used.
				# ... skipped, see Google doc for details
		 */
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



// struct ChainInfo {
// 	ChainInfo() :
// 		CDR1_begin( _CDR_undetected_ ), CDR1_stop( _CDR_undetected_ ),
// 		CDR2_begin( _CDR_undetected_ ), CDR2_stop( _CDR_undetected_ ),
// 		CDR3_begin( _CDR_undetected_ ), CDR3_stop( _CDR_undetected_ ),

// 		FR1_begin( _CDR_undetected_ ), FR1_stop( _CDR_undetected_ ),
// 		FR2_begin( _CDR_undetected_ ), FR2_stop( _CDR_undetected_ ),
// 		FR3_begin( _CDR_undetected_ ), FR3_stop( _CDR_undetected_ ),
// 		FR4_begin( _CDR_undetected_ ), FR4_stop( _CDR_undetected_ ) {}


// 	uint CDR1_begin, CDR1_stop,
// 		CDR2_begin, CDR2_stop,
// 		CDR3_begin, CDR3_stop,

// 		FR1_begin, FR1_stop,
// 		FR2_begin, FR2_stop,
// 		FR3_begin, FR3_stop,
// 		FR4_begin, FR4_stop;
// 	std::string sequence;
// };


} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__
