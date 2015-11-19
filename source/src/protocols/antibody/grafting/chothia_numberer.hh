// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/grafting/chothia_numberer.hh
/// @brief Chothia numberer
/// @author Sergey Lyskov
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)



#ifndef INCLUDED_protocols_antibody_grafting_chothia_numberer_HH
#define INCLUDED_protocols_antibody_grafting_chothia_numberer_HH

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__


#include <protocols/antibody/grafting/antibody_sequence.hh>

#include <utility/vector0.hh>
#include <string>

namespace protocols {
namespace antibody {
namespace grafting {


struct AntibodyChainNumbering
{
	typedef utility::vector0<std::string> NumberingVector;

	NumberingVector fr1, cdr1, fr2, cdr2, fr3, cdr3, fr4;

	/// @brief Check if numbering array length match content of AntibodyChain and AntibodyFramework
	/// @brief throws _AE_unexpected_region_length_ if mismatch detected
	void validate(AntibodyChain const &, AntibodyFramework const &, std::string const & chain_id);

	NumberingVector all();
};


struct AntibodyNumbering
{
	AntibodyChainNumbering heavy, light;
};

/// @brief IO operator for debug and Python bindings
std::ostream & operator << (std::ostream & os, AntibodyChainNumbering const &);
std::ostream & operator << (std::ostream & os, AntibodyNumbering const &);



/// @brief Base class for antibody CDR detector. Sub-class it to implement particular detection methods
class Numberer {
public:
	virtual ~Numberer() {}

	/// @brief Number CDR's
	/// @throw _AE_unexpected_cdr_region_ (or it sub-class) if numbering is fails
	/// @throw _AE_cdr_undefined_ if reuested CDR regions is not defined.
	/// @throw _AE_invalid_cdr_region_ if CDR sequence out of bound of Chain sequence
	virtual AntibodyNumbering number(AntibodySequence const &, AntibodyFramework const &heavy_fr, AntibodyFramework const &light_fr) = 0;
};


/// @brief Use RegEx and antibody sequence information to detect CDR's
class Chothia_Numberer : Numberer {
public:

	/// @brief Detect CDR's
	/// @throw _AE_unexpected_cdr_region_ (sub-class of _AE_numbering_failed_) if numbering is fails
	/// @throw _AE_cdr_undefined_ if reuested CDR regions is not defined.
	/// @throw _AE_invalid_cdr_region_ if CDR sequence out of bound of Chain sequence
	virtual AntibodyNumbering number(AntibodySequence const &, AntibodyFramework const &heavy_fr, AntibodyFramework const &light_fr);

	AntibodyChainNumbering number_heavy_chain(AntibodySequence const &, AntibodyFramework const &heavy_fr);
	AntibodyChainNumbering number_light_chain(AntibodySequence const &, AntibodyFramework const &light_fr);
};


/// @brief IO operator for debug and Python bindings
std::ostream & operator << (std::ostream & os, AntibodyChainNumbering const &);



} // namespace grafting
} // namespace antibody
} // namespace protocols


#endif // __ANTIBODY_GRAFTING__


#endif // INCLUDED_protocols_antibody_grafting_chothia_numberer_HH
