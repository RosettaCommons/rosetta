// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/grafting/chothia_numberer.cc
/// @brief Chothia numberer
/// @author Sergey Lyskov
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__


#include <protocols/antibody/grafting/chothia_numberer.hh>
#include <protocols/antibody/grafting/exception.hh>

#include <utility/string_util.hh>
#include <utility/stream_util.hh>

#include <sstream>
#include <regex>

namespace protocols {
namespace antibody {
namespace grafting {

using std::string;

_AE_unexpected_region_length_ unexpected_region_length_error(string const& region, int length, string const& msg)
{
	std::stringstream s;
	s << "ERROR: Unxpected length of " << region << " [length=" << length << "] " << msg;
	return _AE_unexpected_region_length_( s.str() );
}


_AE_unexpected_region_length_ numbering_region_lengths_error(string const& region, int length1, int length2)
{
	std::stringstream s;
	s << "ERROR: Length of numbering array for region " << region << " does not match does not match length of same region in AntibodyChain/AntibodyFramework! [" << length1 << "!=" << length2 <<"]";
	return _AE_unexpected_region_length_( s.str() );
}



/// @brief Check if numbering array length match content of AntibodyChain and AntibodyFramework
/// @brief throws _AE_unexpected_region_length_ if mismatch detected
void AntibodyChainNumbering::validate(AntibodyChain const &c, AntibodyFramework const &f, string const &chain_id)
{
	if( c.cdr1.size() != cdr1.size() ) throw numbering_region_lengths_error(chain_id+'1', c.cdr1.size(), cdr1.size());
	if( c.cdr2.size() != cdr2.size() ) throw numbering_region_lengths_error(chain_id+'2', c.cdr2.size(), cdr2.size());
	if( c.cdr3.size() != cdr3.size() ) throw numbering_region_lengths_error(chain_id+'3', c.cdr3.size(), cdr3.size());

	if( f.fr1.size() != fr1.size() ) throw numbering_region_lengths_error("fr_"+chain_id+'1', f.fr1.size(), fr1.size());
	if( f.fr2.size() != fr2.size() ) throw numbering_region_lengths_error("fr_"+chain_id+'2', f.fr2.size(), fr2.size());
	if( f.fr3.size() != fr3.size() ) throw numbering_region_lengths_error("fr_"+chain_id+'3', f.fr3.size(), fr3.size());
	if( f.fr4.size() != fr4.size() ) throw numbering_region_lengths_error("fr_"+chain_id+'4', f.fr4.size(), fr4.size());
}


AntibodyChainNumbering::NumberingVector AntibodyChainNumbering::all()
{
	NumberingVector n;
	n.insert(n.end(), fr1 .begin(), fr1 .end() );
	n.insert(n.end(), cdr1.begin(), cdr1.end() );
	n.insert(n.end(), fr2 .begin(), fr2 .end() );
	n.insert(n.end(), cdr2.begin(), cdr2.end() );
	n.insert(n.end(), fr3 .begin(), fr3 .end() );
	n.insert(n.end(), cdr3.begin(), cdr3.end() );
	n.insert(n.end(), fr4 .begin(), fr4 .end() );

	return n;
}


std::ostream & operator << (std::ostream & os, AntibodyChainNumbering const &f)
{
	struct { string name; AntibodyChainNumbering::NumberingVector const &n; } J[] {
		{" fr1", f.fr1}, {"cdr1", f.cdr1}, {" fr2", f.fr2}, {"cdr2", f.cdr2}, {" fr3", f.fr3}, {"cdr3", f.cdr3}, {" fr4", f.fr4},
	};

	os << "AntibodyChainNumbering {" << std::endl;

	for(auto &j : J) {
		os << "  " << j.name << " : " << j.n << std::endl;
	}
	os << '}';

	return os;
}

// std::ostream & operator << (std::ostream & os, AntibodyChainNumbering const &c)
// {
// 	using std::make_pair;
// 	std::pair <string, utility::vector0<string> > o[] {
// 		make_pair("fr1", c.fr1), make_pair("cdr1", c.cdr1),
// 		make_pair("fr2", c.fr2), make_pair("cdr2", c.cdr2),
// 		make_pair("fr3", c.fr3), make_pair("cdr3", c.cdr3), make_pair("fr4", c.fr4)
// 	};

// 	os << "AntibodyChainNumbering{";

// 	for(auto const &p : o) {
// 		os << p.first << ':' << p.second;
// 		if (p.first != "fr4") os << ", ";
// 	}
// 	os << '}';

// 	return os;
// }



std::ostream & operator << (std::ostream & os, AntibodyNumbering const &an)
{
	os << "Heavy " << an.heavy << std::endl;
	os << "Light " << an.light << std::endl;
	return os;
}



AntibodyNumbering Chothia_Numberer::number(AntibodySequence const &ab, AntibodyFramework const &heavy_fr, AntibodyFramework const &light_fr)
{
	AntibodyNumbering an;

	an.heavy = number_heavy_chain(ab, heavy_fr);
	an.light = number_light_chain(ab, light_fr);

	return an;
}


AntibodyChainNumbering::NumberingVector number_region(string const& name, uint length, std::vector<string> const & n /*numberings*/)
{
	std::vector<uint> possible_lengths;

	for(auto const &s : n) {
		auto nb = utility::string_split(s, ',' );
		if( nb.size() == length ) return nb;
		possible_lengths.push_back( nb.size() );
	}

	std::stringstream s;
	s << "ERROR: Unxpected length of " << name << " [length=" << length << "], length expected to be: " << possible_lengths << '!';
	throw _AE_unexpected_region_length_( s.str() );
}


AntibodyChainNumbering Chothia_Numberer::number_heavy_chain(AntibodySequence const &ab, AntibodyFramework const &f)
{
	ab.heavy.validate();

	//AntibodyFramework const f( ab.heavy_framework() );
	AntibodyChainNumbering n;

	n.fr1 = number_region("heavy-fr1", f.fr1.size(), {
			"10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25",
			"9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25",
			"8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25",
			"7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25",
			"6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25",
			"5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25",
			"4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25",
			"3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25",
			"2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25",
			"1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25",
	});

	n.fr2 = number_region("heavy-fr2", f.fr2.size(), { "36,37,38,39,40,41,42,43,44,45,46,47,48,49"} );

	n.fr3 = number_region("heavy-fr3", f.fr3.size(), {
			"66,67,68,69,70,71,72,73,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94",
			"66,67,68,69,70,71,72,73,74,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94",
			"66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94",
	});


	n.fr4 = number_region("heavy-fr4", f.fr4.size(), {
			"103,104,105,106,107,108,109,110,111",         // Case len=9.  Note: this is different from Python implementaion (it was silently failing it)
			"103,104,105,106,107,108,109,110,111,112",     // Case len=10. Note: this is different from Python implementaion (it was silently failing it)
			"103,104,105,106,107,108,109,110,111,112,113", // Case len=11. Note: this is different from Python implementaion (it was silently failing it)
			"103,104,105,106,107,108,109,110,111,112,113,114"
	});


	n.cdr1 = number_region("h1", ab.heavy.cdr1.size(), {
			"26,27,32,33,34,35",
			"26,27,28,32,33,34,35",
			"26,27,28,29,32,33,34,35",
			"26,27,28,29,30,32,33,34,35",
			"26,27,28,29,30,31,32,33,34,35",
			"26,27,28,29,30,31,31A,32,33,34,35",
			"26,27,28,29,30,31,31A,31B,32,33,34,35",
			"26,27,28,29,30,31,31A,31B,31C,32,33,34,35",
	});

	n.cdr2= number_region("h2", ab.heavy.cdr2.size(), {
			"50,51,52,57,58,59,60,61,62,63,64,65",
			"50,51,52,56,57,58,59,60,61,62,63,64,65",
			"50,51,52,55,56,57,58,59,60,61,62,63,64,65",
			"50,51,52,54,55,56,57,58,59,60,61,62,63,64,65",
			"50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65",
			"50,51,52,52A,53,54,55,56,57,58,59,60,61,62,63,64,65",
			"50,51,52,52A,52B,53,54,55,56,57,58,59,60,61,62,63,64,65",
			"50,51,52,52A,52B,52C,53,54,55,56,57,58,59,60,61,62,63,64,65",
			"50,51,52,52A,52B,52C,52D,53,54,55,56,57,58,59,60,61,62,63,64,65",
			"50,51,52,52A,52B,52C,52D,52E,53,54,55,56,57,58,59,60,61,62,63,64,65",
			"50,51,52,52A,52B,52C,52D,52E,52F,53,54,55,56,57,58,59,60,61,62,63,64,65",
	});

	n.cdr3= number_region("h3", ab.heavy.cdr3.size(), {
			"95,101,102",
			"95,96,101,102",
			"95,96,97,101,102",
			"95,96,97,98,101,102",
			"95,96,97,98,99,101,102",
			"95,96,97,98,99,100,101,102",
			"95,96,97,98,99,100,100A,101,102",
			"95,96,97,98,99,100,100A,100B,101,102",
			"95,96,97,98,99,100,100A,100B,100C,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,100X,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,100X,100Y,101,102",
			"95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,100X,100Y,100Z,101,102",
	});


	n.validate(ab.heavy, f, "h");
	return n;
}


AntibodyChainNumbering Chothia_Numberer::number_light_chain(AntibodySequence const &ab, AntibodyFramework const &f)
{
	ab.light.validate();

	//AntibodyFramework const f( ab.light_framework() );
	AntibodyChainNumbering n;

	std::smatch _;
	string pattern_a("[A-Z][QE][A-Z]{9}[A-Z][A-Z]{4}[LVIMF][A-Z]C");
	string pattern_b("[A-Z][QE][A-Z]{8}[A-Z][A-Z]{4}[LVIMF][A-Z]C");

	if( ! antibody_grafting_usable() ) {
		utility_exit_with_message("ERROR: Your compiler does not have full support for C++11 regex, and therefore can't support Chothia_Numberer/antibody grafting.");
	}
	std::regex re_a(pattern_a, std::regex::extended );
	std::regex re_b(pattern_b, std::regex::extended );
	if( std::regex_search( f.fr1, _, re_a ) ) {
		n.fr1 = number_region("light-fr1", f.fr1.size(), {
				"5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23",
				"4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23",
				"3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23",
				"2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23",
				"1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23",
				//"0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23",  // Note: this is different from Python implementaion (it was adjusteed L1 on the fly, deleting the first residue)
		});
	} else if ( std::regex_search( f.fr1, _, re_b ) ) {
		n.fr1 = number_region("light-fr1", f.fr1.size(), {
				"4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23",
				"3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23",
				"2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23",
				"1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23",
				//"0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23",  // Note: this is different from Python implementaion (it was adjusteed L1 on the fly, deleting the first residue)
				//"0,0A,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23",  // Note: this is different from Python implementaion (it was adjusteed L1 on the fly, deleting the first residue)
		});
	} else {
		throw _AE_unexpected_region_length_("ERROR: Current code could not assign Chothia numbering of FR_L1 in the query sequence!");
	}

	n.fr2 = number_region("light-fr2", f.fr2.size(), {"35,36,37,38,39,40,41,42,43,44,45,46,47,48,49"});

	n.fr3 = number_region("light-fr3", f.fr3.size(), {
			"57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88",
			"57,58,59,60,61,62,63,64,65,66,66A,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88",
			"57,58,59,60,61,62,63,64,65,66,66A,66B,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88",
	});

	n.fr4 = number_region("light-fr4", f.fr4.size(), {
			"98,99,100,101,102,103,104,105,106,107",      // Case len=10. Note: this is different from Python implementaion (it was silently failing it)
			"98,99,100,101,102,103,104,105,106,107,108",  // Case len=11. Note: this is different from Python implementaion (it was silently failing it)
			"98,99,100,101,102,103,104,105,106,107,108,109"
	});


	n.cdr1 = number_region("l1", ab.light.cdr1.size(), {
			"24,25,26,27,28,29,30,34",
			"24,25,26,27,28,29,30,33,34",
			"24,25,26,27,28,29,30,32,33,34",
			"24,25,26,27,28,29,30,31,32,33,34",
			"24,25,26,27,28,29,30,30A,31,32,33,34",
			"24,25,26,27,28,29,30,30A,30B,31,32,33,34",
			"24,25,26,27,28,29,30,30A,30B,30C,31,32,33,34",
			"24,25,26,27,28,29,30,30A,30B,30C,30D,31,32,33,34",
			"24,25,26,27,28,29,30,30A,30B,30C,30D,30E,31,32,33,34",
			"24,25,26,27,28,29,30,30A,30B,30C,30D,30E,30F,31,32,33,34",
	});

	n.cdr2 = number_region("l2", ab.light.cdr2.size(), {
			"50,51,52,53,54,55,56",
			"50,51,52,53,54,54A,55,56",
			"50,51,52,53,54,54A,54B,55,56",
			"50,51,52,53,54,54A,54B,54C,55,56",
			"50,51,52,53,54,54A,54B,54C,54D,55,56",
	});

	n.cdr3 = number_region("l3", ab.light.cdr3.size(), {
			"89,90,91,92,97",
			"89,90,91,92,93,97",
			"89,90,91,92,93,94,97",
			"89,90,91,92,93,94,95,97",
			"89,90,91,92,93,94,95,96,97",
			"89,90,91,92,93,94,95,95A,96,97",
			"89,90,91,92,93,94,95,95A,95B,96,97",
			"89,90,91,92,93,94,95,95A,95B,95C,96,97",
			"89,90,91,92,93,94,95,95A,95B,95C,95D,96,97",
			"89,90,91,92,93,94,95,95A,95B,95C,95D,95E,96,97",
			"89,90,91,92,93,94,95,95A,95B,95C,95D,95E,95F,96,97",
	});


	n.validate(ab.light, f, "l");
	return n;
}


} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__
