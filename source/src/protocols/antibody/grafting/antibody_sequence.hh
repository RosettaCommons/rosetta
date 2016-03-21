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


uint const _CDR_max_length_(std::numeric_limits< platform::Size >::max() / 8);
uint const _CDR_undefined_(_CDR_max_length_+1);
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
						  fr1_end(_FR_undefined_),   fr2_end(_FR_undefined_),   fr3_end(_FR_undefined_),   fr4_end(_FR_undefined_) {}

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


} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__

#endif // INCLUDED_protocols_antibody_grafting_AntibodySequence_hh
