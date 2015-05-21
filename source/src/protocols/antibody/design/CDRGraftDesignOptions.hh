// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/design/CDRGraftDesignOptions.hh
/// @brief GraftDesign options class and parser
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_design_CDRGraftDesignOptions_hh
#define INCLUDED_protocols_antibody_design_CDRGraftDesignOptions_hh

#include <protocols/antibody/design/CDRGraftDesignOptions.fwd.hh>

#include <protocols/antibody/AntibodyEnumManager.fwd.hh>
#include <protocols/antibody/design/AntibodyDesignEnum.hh>
#include <protocols/antibody/AntibodyEnum.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace antibody {
namespace design {

	typedef utility::vector1<CDRGraftDesignOptionsOP> AntibodyCDRGraftDesignOptions;

class CDRGraftDesignOptions : public utility::pointer::ReferenceCount {
public:

	CDRGraftDesignOptions();
	CDRGraftDesignOptions(CDRNameEnum cdr);
	CDRGraftDesignOptions(CDRGraftDesignOptions const & src);

	virtual ~CDRGraftDesignOptions();

public:

	/// @brief Set the CDR type for this instance of options.
	void
	set_cdr(CDRNameEnum cdr);

	/// @brief Should this CDR be designed?
	void
	design(bool design);

	/// @brief Should this CDR be designed?
	bool
	design() const {
		return design_;
	}

	/// @brief Set the weight we will use for choosing this particular CDR for design.
	void
	weight(core::Real weight);

	/// @brief Return the weight we will use for choosing this particular CDR for design.
	core::Real
	weight() const {
		return cdr_weight_;
	}


	/// @brief Options are: relax, minimize, repack, dualspace. Default is repack for time, but relax is the best
	void
	mintype(MinTypeEnum mintype);

	MinTypeEnum
	mintype() const {
		return mintype_;
	}

	/// @brief Minimize neighbor sidechains?
	void
	min_neighbor_sc(bool min_neighbor_sc);

	bool
	min_neighbor_sc() const {
		return min_neighbor_sc_;
	}

	/// @brief Minimize sidechains during minimization (relax/min/etc) (default True?
	void
	min_sc(bool min_sc);

	bool
	min_sc() const {
		return min_sc_;
	}

	/// @brief  Use rigid body optimization during minimization of the CDR (Minimize jump bt antigen and antibody)?
	/// @details Does nothing if repacking the cdr.  Default False.
	void
	min_rb(bool min_rb);

	bool
	min_rb() const {
		return min_rb_;
	}

	/// @brief Set other CDRs to minimize when this CDR is minimized
	void
	neighbor_min(utility::vector1<CDRNameEnum> neighbor_min);

	void
	neighbor_min_add(CDRNameEnum cdr);

	void
	neighbor_min_clear();

	utility::vector1<CDRNameEnum>
	neighbor_min() const {
		return neighbor_min_;
	}


public:

	CDRGraftDesignOptionsOP
	clone() const;

	/// @brief Set defaults - !!! these are not the defaults from the database default instruction file.
	/// Use the parser to parse those instructions if you wish.
	void
	set_defaults();

private:

	bool cdr_;
	bool design_;
	MinTypeEnum mintype_;
	bool min_neighbor_sc_;
	bool min_sc_;
	bool min_rb_;
	utility::vector1<CDRNameEnum> neighbor_min_;
	core::Real cdr_weight_;

};

class CDRGraftDesignOptionsParser : public utility::pointer::ReferenceCount {
public:

	CDRGraftDesignOptionsParser();
	// Undefined, commenting out to fix PyRosetta build  CDRGraftDesignOptionsParser(CDRGraftDesignOptions const & src);

	virtual ~CDRGraftDesignOptionsParser();

	CDRGraftDesignOptionsOP
	parse_options(CDRNameEnum cdr, std::string filename);

	/// @brief Parse default_instructions (mainly used for AbDesign) then parse user file
	CDRGraftDesignOptionsOP
	parse_default_and_user_options(CDRNameEnum cdr, std::string filename);


	//ALL CDRs

	utility::vector1<CDRGraftDesignOptionsOP>
	parse_options(std::string filename);

	utility::vector1<CDRGraftDesignOptionsOP>
	parse_default_and_user_options(std::string filename);

private:
	/// @brief Tries to find the path in either database, relative, or absolute, sets it.
	void
	check_path();

	void
	check_line_len(const utility::vector1<std::string> & lineSP, const core::Size len_check) const;

	void
	parse_cdr_option(std::string const mode, vector1<std::string> & lineSP);


	void
	parse_cdr_graft_option(std::string const adjective, utility::vector1<std::string> & lineSP) ;

	void
	parse_cdr_general_option(utility::vector1<std::string> & lineSP) ;


	void
	set_cdr_graft_mintype_options(std::string const mintype);

	void
	set_cdr_graft_neighbor_mintype_options(utility::vector1<std::string> & lineSP);


	void
	set_cdr_graft_general_option(std::string const option);


private:

	std::string instructions_path_;
	AntibodyEnumManagerOP ab_manager_;
	CDRGraftDesignOptionsOP cdr_options_;
	bool default_and_user_;

};

}
}
}


#endif	//INCLUDED_ CDRGraftDesignOptions.hh
