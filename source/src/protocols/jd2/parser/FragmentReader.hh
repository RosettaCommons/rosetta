// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/jd2/parser/FragmentReader.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

#ifndef INCLUDED_protocols_jd2_parser_FragmentReader_hh
#define INCLUDED_protocols_jd2_parser_FragmentReader_hh

// Unit Headers
#include <protocols/jd2/parser/FragmentReader.fwd.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>
// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>


// c++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {
namespace parser {

class FragmentReader : public utility::pointer::ReferenceCount {
public:


	typedef utility::pointer::ReferenceCount Parent;


public:

	typedef core::Size Size;
	typedef std::string String;
	typedef core::pose::Pose Pose;
	typedef core::fragment::FragSetOP FragSetOP;

	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;


public:


	/// @brief default constructor
	FragmentReader();

	/// @brief value constructor
	FragmentReader( TagCOP const & tag );

	/// @brief destructor
	virtual ~FragmentReader();

public:


	/// @brief main opeartion
	void apply( FragSetOP & fragset );


private:


	/// @brief parse tag
	void parse_tag( TagCOP const & tag );

	/// @brief set fragments just for helper function
	void set_fragments( Pose const & pose, FragSetOP const & fragset );


private:

	/// @brief way of reading fragments from pdbs, silent, fragfiles, or vall
	String read_type_;

	/// @brief input file name, pdbs, silent, fragfiles of vall
	String filename_;

	/// @brief length of fragments
	Size frag_size_;

	/// @brief the begin of sequence positions where fragments are stealed
	Size begin_;

	/// @brief the end of sequence positions where fragments are stealed
	Size end_;

	/// @brief number of stealing fragments if read_type_ is pdbs or silent
	Size steal_times_;

	/// @brief maximum number of fragments when a silent file is used
	Size nfrags_;

	/// @brief secondary structure assignment to pick fragment from vall
	String ss_;

	/// @brief amino acids to pick fragment from vall
	String aa_;

	/// @brief abego assignment to pick fragment from vall
	bool use_abego_;

	/// @brief blueprint
	protocols::jd2::parser::BluePrintOP blueprint_;

};

} //namespace parser
} //namespace jd2
} //namespace protocols

#endif
