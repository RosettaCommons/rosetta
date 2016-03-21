// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/antibody/grafting/AntibodySequence.cc
/// @brief Helper classes to store parsed antibody sequence data
/// @author Sergey Lyskov
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__


#include <protocols/antibody/grafting/antibody_sequence.hh>

#include <protocols/antibody/grafting/exception.hh>



namespace protocols {
namespace antibody {
namespace grafting {

using std::string;


/// @brief /// Check if current value of begin/end could be valid for some sequence
bool CDR_Bounds::defined() const
{
	if( begin<_CDR_undefined_ && end<_CDR_undefined_ ) return true;
	else return false;
}


/// @brief Check if current value of begin/end if valid for given sequence. Always return true if CDR undefined.
bool CDR_Bounds::valid(string const &sequence) const
{
	if( defined() ) {
		if( begin <= end  &&  end < sequence.size() ) { return true; }
		else { return false; }
	} else {
		return true;
	}
}

/// @details validate CDR sequence
/// @throw _AE_cdr_undefined_ if reuested CDR regions is not defined.
/// @throw _AE_invalid_cdr_region_ if CDR sequence out of bound of Chain sequence
void CDR_Bounds::validate(string const &sequence) const
{
	if( !defined() ) throw _AE_cdr_undefined_();
	else {
		if( !valid(sequence) ) { throw _AE_invalid_cdr_region_(); }
	}
}


/// @details update fr1/4 sequences.
/// @throw _AE_invalid_cdr_region_ if CDR sequence out of bound of Chain sequence
void AntibodyFramework::update_sequences(std::string chain_sequence)
{
	for(auto i : {fr1_begin, fr2_begin, fr3_begin, fr3_begin, fr1_end, fr2_end, fr3_end, fr3_end} ) {
		if( i > chain_sequence.size() ) throw _AE_invalid_cdr_region_();
	}

	fr1 = chain_sequence.substr(fr1_begin, fr1_end-fr1_begin);
	fr2 = chain_sequence.substr(fr2_begin, fr2_end-fr2_begin);
	fr3 = chain_sequence.substr(fr3_begin, fr3_end-fr3_begin);
	fr4 = chain_sequence.substr(fr4_begin, fr4_end-fr4_begin);
}


/*! \details Calculate CDR sequence
 *  \throw _AE_cdr_undefined_ if reuested CDR regions is not defined.
 *  \throw _AE_invalid_cdr_region_ if CDR sequence out of bound of Chain sequence
 */
string cdr_sequence(string const & sequence, CDR_Bounds const &cdr)
{
	if( !cdr.defined() ) throw _AE_cdr_undefined_();
	else {
		if( cdr.valid(sequence) ) return sequence.substr(cdr.begin, cdr.end-cdr.begin);
		else throw _AE_invalid_cdr_region_();
	}
}

/*! \details generate CDR info object which hold CDR begin/end and loop sequence
 *  \throw _AE_cdr_undefined_ if reuested CDR regions is not defined.
 *  \throw _AE_invalid_cdr_region_ if CDR sequence out of bound of Chain sequence
 */
// CDR cdr_sequence(string const & sequence, CDR_Bounds const &cdr)
// {
// }

string AntibodyChain::cdr1_sequence() const { return cdr_sequence(sequence, cdr1); }
string AntibodyChain::cdr2_sequence() const { return cdr_sequence(sequence, cdr2); }
string AntibodyChain::cdr3_sequence() const { return cdr_sequence(sequence, cdr3); }


/*! \details  validate CDR sequence's
 *  @throw _AE_cdr_undefined_ if reuested CDR regions is not defined.
 *  @throw _AE_invalid_cdr_region_ if CDR sequence out of bound of Chain sequence
 */
void AntibodyChain::validate() const
{
	auto cdrs = {cdr1, cdr2, cdr3};

	for(auto const &c : cdrs ) c.validate(sequence);
}



AntibodyFramework AntibodySequence::heavy_framework() const
{
	AntibodyFramework f;

	heavy.validate();

	f.fr1_begin = 0;  f.fr1_end = heavy.cdr1.begin;
	// Moved to into scs_blast.cc:trim_framework function  if( f.fr1_end - f.fr1_begin > 24 ) f.fr1_begin = f.fr1_end - f.fr1_begin - 24;

	f.fr2_begin = heavy.cdr1.end;  f.fr2_end =  heavy.cdr1.end + heavy.cdr2.begin - heavy.cdr1.end;  // FR_H2 = heavy_chain[H1_end + 1: H1_end + 1 + H2_begin - H1_end - 1]
	f.fr3_begin = heavy.cdr2.end;  f.fr3_end =  heavy.cdr2.end + heavy.cdr3.begin - heavy.cdr2.end;  // FR_H3 = heavy_chain[H2_end + 1: H2_end + 1 + H3_begin - H2_end - 1]
	f.fr4_begin = heavy.cdr3.end;  f.fr4_end =  heavy.sequence.size();  // Moved to into scs_blast.cc:trim_framework function  f.fr4_begin = heavy.cdr3.end;  f.fr4_end =  heavy.cdr3.end + 12; // FR_H4 = heavy_chain[H3_end + 1: H3_end + 1 + 12] FR_H4 = heavy_chain[H3_end + 1: H3_end + 1 + 12]

	f.update_sequences(heavy.sequence);

	return f;
}


AntibodyFramework AntibodySequence::light_framework() const
{
	AntibodyFramework f;

	light.validate();

	f.fr1_begin = 0;  f.fr1_end = light.cdr1.begin;
	// Moved to into scs_blast.cc:trim_framework function 	if( f.fr1_end - f.fr1_begin > 23 ) f.fr1_begin = f.fr1_end - f.fr1_begin - 23;

	f.fr2_begin = light.cdr1.end;  f.fr2_end = light.cdr1.end + light.cdr2.begin - light.cdr1.end; // Moved to into scs_blast.cc:trim_framework function 	f.fr2_begin = light.cdr1.end;  f.fr2_end = light.cdr1.end + 15;  // FR_L2 = light_chain[ L1_end + 1  :  L1_end + 1 + 15]
	f.fr3_begin = light.cdr2.end;  f.fr3_end = light.cdr2.end + light.cdr3.begin - light.cdr2.end;  // FR_L3 = light_chain[ L2_end + 1  :  L2_end + 1 + L3_begin - L2_end - 1 ]
    f.fr4_begin = light.cdr3.end;  f.fr4_end = light.sequence.size();  // Moved to into scs_blast.cc:trim_framework function  f.fr4_begin = light.cdr3.end;  f.fr4_end = light.cdr3.end + 12;  // FR_L4 = light_chain[ L3_end + 1  :  L3_end + 1 + 12 ]

	f.update_sequences(light.sequence);

	return f;
}


string AntibodySequence::h1_sequence() const { return heavy.cdr1_sequence(); }
string AntibodySequence::h2_sequence() const { return heavy.cdr2_sequence(); }
string AntibodySequence::h3_sequence() const { return heavy.cdr3_sequence(); }

string AntibodySequence::l1_sequence() const { return light.cdr1_sequence(); }
string AntibodySequence::l2_sequence() const { return light.cdr2_sequence(); }
string AntibodySequence::l3_sequence() const { return light.cdr3_sequence(); }



/// @brief IO operator for debug and Python bindings
std::ostream & operator << (std::ostream & os, CDR_Bounds const &b)
{
	os << "CDR_Bounds{begin=" << b.begin << ", end=" << b.end << ", size=" << b.size();
	if( !b.defined() ) os << ", defined=" << b.defined();
	os << '}';
	return os;
}


/// @brief IO operator for debug and Python bindings
std::ostream & operator << (std::ostream & os, AntibodyFramework const &f)
{

	uint   begin[]    {f.fr1_begin, f.fr2_begin, f.fr3_begin, f.fr4_begin};
	uint   end[]      {f.fr1_end,   f.fr2_end,   f.fr3_end,   f.fr4_end};
	string sequence[] {f.fr1,       f.fr2,       f.fr3,       f.fr4};

	os << "AntibodyFramework{";

	for(int i=0; i<4; ++i) {
		os << "fr" << i+1 << '[' << begin[i] << ':' << end[i] << ", length=" << end[i]- begin[i] << ", " << sequence[i] << ']';
		if( i < 3 ) os << ", ";
	}
	os << '}';

	return os;
}


} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__
