// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/enzdes/movers/AddorRemoveCsts.hh
/// @author Florian Richter (floric@u.washington.edu)

#ifndef INCLUDED_protocols_enzdes_AddorRemoveCsts_hh
#define INCLUDED_protocols_enzdes_AddorRemoveCsts_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace enzdes {

enum CstAction {
	VOID = 1,
	ADD_NEW,
	ADD_PREGENERATED,
	REMOVE
};


/// @brief A simple wrapper to get the functionality in EnzConstraintIO
/// into mover format
class AddOrRemoveMatchCsts : public protocols::moves::Mover {

public:  //Constructor / Destructor

	AddOrRemoveMatchCsts();

	AddOrRemoveMatchCsts( AddOrRemoveMatchCsts const & other );

	~AddOrRemoveMatchCsts();

public:

	protocols::moves::MoverOP clone() const;
	virtual std::string get_name() const;

	protocols::moves::MoverOP fresh_instance() const;

	void apply( core::pose::Pose & pose );

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

	// @brief Set the name of the constraint file. Overwrites the file name that might be read in during parse_my_tag.
	void cstfile( std::string const & setting );

	static
	toolbox::match_enzdes_util::EnzConstraintIOCOP
	get_const_EnzConstraintIO_for_cstfile( std::string cstfile = "" );

	void set_cst_action(CstAction action){ cst_action_=action; }
	void set_accept_blocks_missing_header( bool setting ){ accept_blocks_missing_header_ = setting; }
	void set_keep_covalent( bool setting ){ keep_covalent_ = setting; }


protected:

	toolbox::match_enzdes_util::EnzConstraintIOOP
	get_EnzConstraintIO_for_cstfile(
		std::string const & cstfile
	);

private:

	/// Save the contents of the constraint files that are read in for reuse.
	static std::map< std::string, toolbox::match_enzdes_util::EnzConstraintIOOP > cstfile_map_;

	std::string option_cstfile_; // Read from options system if no user-defined constraint file is given.
	std::string cstfile_; // May be set either by the parser or programmatically.
	CstAction cst_action_;
	bool keep_covalent_, accept_blocks_missing_header_, fail_on_constraints_missing_;

	///atm this scorefunction is only used if the user specifes a covalent ambiguous constraint
	///in which case the ambiguity is resolved at the time of newly adding the constraints, and the
	///covalent connection established according to the then best constraints.
	core::scoring::ScoreFunctionOP sfxn_;
};

} // enzdes
} // protocols


#endif /*INCLUDED_protocols_enzdes_AddorRemoveCsts_HH*/
