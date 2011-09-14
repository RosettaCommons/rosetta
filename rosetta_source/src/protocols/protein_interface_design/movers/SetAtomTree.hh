// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/SetAtomTree.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu), Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_SetAtomTree_hh
#define INCLUDED_protocols_protein_interface_design_movers_SetAtomTree_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief a mover that sets a user-defined atom tree
class SetAtomTree : public protocols::moves::Mover
{
public :
	SetAtomTree();
	SetAtomTree( core::kinematics::FoldTreeOP ft );
	virtual ~SetAtomTree();
	void fold_tree( core::kinematics::FoldTreeOP ft );
	core::kinematics::FoldTreeOP fold_tree() const;
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new SetAtomTree ); }
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	core::kinematics::FoldTreeOP create_atom_tree( core::pose::Pose const & pose, core::Size const host_chain, core::Size const resnum, core::Size const anchor_num_ = 0, std::string const connect_to = "", std::string connect_from = "" );//if connect_to or connect_from = "" optimal_connection_point is invoked
private :
	bool docking_ft_; //dflt false
	core::Size jump_; //dflt true
	std::string resnum_, connect_to_, anchor_res_, connect_from_; //as parsed
	core::Size host_chain_; //dflt 2
};

} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_SetAtomTree_HH*/
