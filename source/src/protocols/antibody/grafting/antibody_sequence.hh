// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/antibody/grafting/AntibodySequence.hh
/// @brief Helper classes to store parsed antibody sequence data
/// @author Sergey Lyskov
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

#ifndef INCLUDED_protocols_antibody_grafting_AntibodySequence_hh
#define INCLUDED_protocols_antibody_grafting_AntibodySequence_hh

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <core/types.hh>

#include <string>
#include <vector>

#include <iostream>
#include <limits>

namespace protocols {
namespace antibody {
namespace grafting {

typedef core::uint uint;


uint const _CDR_max_length_(std::numeric_limits< platform::Size >::max() / 10);
uint const _CDR_undefined_(_CDR_max_length_+1);
//uint const _CDR_undetected_(_CDR_max_length_+2);
uint const _FR_undefined_(_CDR_max_length_+1);


struct CDR_Bounds {
	CDR_Bounds() : begin(_CDR_undefined_), end(_CDR_undefined_) {}

	uint begin; //  @brief index of sequence for first element of region
	uint end;   /// @brief index of sequence for past-the-end of region


	/// @brief return size of region
	uint size() const { return end-begin; }


	/// @brief Check if current value of begin/end could be valid for some sequence
	bool defined() const;


	/// @brief Check if current value of begin/end if valid for given sequence. Always return true if CDR undefined.
	bool valid(std::string const &sequence) const;

	/// @brief validate CDR sequence
	/// @throw _AE_cdr_undefined_ if reuested CDR regions is not defined.
	/// @throw _AE_invalid_cdr_region_ if CDR sequence out of bound of Chain sequence
	void validate(std::string const &sequence) const;
};


struct AntibodyFramework {
	AntibodyFramework() : fr1_begin(_FR_undefined_), fr2_begin(_FR_undefined_), fr3_begin(_FR_undefined_), fr4_begin(_FR_undefined_),
						  fr1_end(_FR_undefined_),   fr2_end(_FR_undefined_),   fr3_end(_FR_undefined_),   fr4_end(_FR_undefined_)    {}

	uint fr1_begin, fr2_begin, fr3_begin, fr4_begin;
	uint fr1_end,   fr2_end,   fr3_end,   fr4_end;

	std::string fr1, fr2, fr3, fr4;

	/// @brief update fr1/4 sequences.
	/// @throw _AE_invalid_cdr_region_ if CDR sequence out of bound of Chain sequence
	void update_sequences(std::string chain_sequence);
};



/// @brief Hold information about heavy or light antibody chain.
///        This include:
///            - loop info (though CDR_Bounds cdr*)
///            - sequnce
struct AntibodyChain {
	typedef std::string string;

	AntibodyChain(string sequence_) : sequence(sequence_) {}

	CDR_Bounds cdr1;
	CDR_Bounds cdr2;
	CDR_Bounds cdr3;

	string sequence;

	/// @brief Calculate CDR sequence
	/// @throw _AE_cdr_undefined_ if reuested CDR regions is not defined.
	/// @throw _AE_invalid_cdr_region_ if CDR sequence out of bound of Chain sequence
	string cdr1_sequence() const; string  cdr2_sequence() const;  string cdr3_sequence() const;


	/// @brief validate CDR sequence's
	/// @throw _AE_cdr_undefined_ if reuested CDR regions is not defined.
	/// @throw _AE_invalid_cdr_region_ if CDR sequence out of bound of Chain sequence
	void validate() const;
};


/// Hold info about whole antibody sequence: heavy and light chains. This is basially a convenience class

struct AntibodySequence {
	typedef std::string string;

	AntibodySequence() : heavy(""), light("") {}
	AntibodySequence(string heavy_, string light_) : heavy(heavy_), light(light_) {}

	AntibodyChain heavy;
	AntibodyChain light;


	/// @brief Calculate Framework sequence's by using CDR's CDR_Bounds info
	/// @throw _AE_cdr_undefined_ if reuested CDR regions is not defined.
	/// @throw _AE_invalid_cdr_region_ if CDR sequence out of bound of Chain sequence
	AntibodyFramework heavy_framework() const;
	AntibodyFramework light_framework() const;


	/// @brief Calculate CDR sequence
	/// @throw _AE_cdr_undefined_ if reuested CDR regions is not defined.
	/// @throw _AE_invalid_cdr_region_ if CDR sequence out of bound of Chain sequence
	string h1_sequence() const;  string h2_sequence() const;  string h3_sequence() const;
	string l1_sequence() const;  string l2_sequence() const;  string l3_sequence() const;
};





/// @brief IO operator for debug and Python bindings
std::ostream & operator << (std::ostream & os, CDR_Bounds const &);
std::ostream & operator << (std::ostream & os, AntibodyFramework const &);



/*

// /// @brief READ-ONLY CDR data holder, this is convenience structure, all this info could be acquired but wit multiple calls
// struct CDR
// {
// 	CDR(uint begin, uint end, sequence)
// 	uint begin;
// 	uint end;
// 	std::string sequence;
// }

	// /// @brief Reatrive READ-ONLY object which hold CDR info (begin+end+cdr_sequence)
	// CDR h1() const;  CDR h2() const;  CDR h3() const;
	// CDR l1() const;  CDR l2() const;  CDR l3() const;


AbNumberer
StructuralComponentSelector
CDRDetector


enum NumberingScheme {
	_NS_Chothia_ = 0,
	_NS_AHo_,
	_NS_Aroop_ = 42,
	_NS_etc_,
	_NS_first_element_ = _NS_Chothia_,
	_NS_number_of_elements_ = _NS_etc_ + 1
}; // numbering schemes!


enum CDRs {
	H1 = 0,
	H2,
	H3,
	L1,
	L2,
	L3,
	_CDRs_first_element_ = H1,
	_CDRs_last_element_ = L3,
	_CDRs_number_of_elements_ = _CDRs_last_element_ + 1,
	_unknown_ = -1
};


enum AntibodyChainID {
	_AB_heavy_chain_ = 0,
	_AB_light_chain_,
	_AB_chains_first_element_ = _AB_heavy_chain_,
	_AB_chains_last_element_ = _AB_light_chain_,
	_AB_chains_number_of_elements_ = _AB_chains_last_element_ + 1,
	_AB_chain_unknown_ = -1
};


int const _Max_length_( INT_MAX / 10 );

struct ChainInfo {
	ChainInfo() :
		CDR1_begin( _Max_length_ ), CDR1_end( _Max_length_ ),
		CDR2_begin( _Max_length_ ), CDR2_end( _Max_length_ ),
		CDR3_begin( _Max_length_ ), CDR3_end( _Max_length_ ),

		FR1_begin( _Max_length_ ), FR1_end( _Max_length_ ),
		FR2_begin( _Max_length_ ), FR2_end( _Max_length_ ),
		FR3_begin( _Max_length_ ), FR3_end( _Max_length_ ),
		FR4_begin( _Max_length_ ), FR4_end( _Max_length_ ),

		id( _AB_chain_unknown_ ) {}


	int CDR1_begin, CDR1_end,
		CDR2_begin, CDR2_end,
		CDR3_begin, CDR3_end,

		FR1_begin, FR1_end,
		FR2_begin, FR2_end,
		FR3_begin, FR3_end,
		FR4_begin, FR4_end;


	std::string sequence;
	AntibodyChainID id;
};





class AntibodyChain {
public:
	AntibodyChain() : identity_( _AB_chain_unknown_ ), CDRs_(3), FR_(), sequence_("") {}
	AntibodyChain( AntibodyChainID id, std::string seq ) : identity_( id ), CDRs_(3), FR_(), sequence_( seq ) {}
	AntibodyChain( ChainInfo ci ) : identity_( ci.id ), CDRs_{ CDRSequence{ L1, ci.CDR1_begin, ci.CDR1_end }, CDRSequence{ L2, ci.CDR2_begin, ci.CDR2_end }, CDRSequence{ L3, ci.CDR3_begin, ci.CDR3_end } } {}


	AntibodyChainID identity() const { return identity_; }

	std::string sequence() const { return sequence_; }

private:


	/// constains sequence (as a string) plus information related to the sequence
	class CDRSequence {
	public:
		CDRSequence() : identity_( _unknown_ ), begin_( _Max_length_ ), end_( _Max_length_ ), template_{ "" }  {}
		CDRSequence( CDRs id, int begin, int end, std::string tmpl = "" ) : identity_( id ), begin_( begin ),
			end_( end ), template_( tmpl ) {}

		std::vector< std::string > compute_numbering( NumberingScheme ) const;
	private:
		CDRs identity_;
		int begin_;
		int end_;
		std::string template_; // maybe it would be useful to include the template we choose to graft for this loop
	};


	/// constains sequence (as a string) plus information related to the sequence
	class FrameworkSequence {
	public:
		FrameworkSequence() : sequence_chunks_(4), identity_( _unknown_ ) {}
		std::vector< std::string > compute_numbering( NumberingScheme ) const;
	private:
		std::vector< std::string > sequence_chunks_; // init to ""
		CDRs identity_;
		// std::string template_; // maybe it would be useful to include the template we choose to graft for this loop
	};


private:
	AntibodyChainID identity_;
	std::vector< CDRSequence > CDRs_;
	FrameworkSequence FR_;
	std::string sequence_;
};


/// Hold info about whole antibody sequence: heavy and light chains. This is basially a convenience class
class AntibodySequence {
public:
	AntibodySequence() : heavy_chain_(), light_chain_() {}

private:
	AntibodyChain heavy_chain_;
	AntibodyChain light_chain_;
};

*/



} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__

#endif // INCLUDED_protocols_antibody_grafting_AntibodySequence_hh
