// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/AddChainBreak.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu), Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_AddChainBreak_hh
#define INCLUDED_protocols_protein_interface_design_movers_AddChainBreak_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Residue.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief a mover that sets a chainbreak in a specified position
class AddChainBreak : public protocols::moves::Mover
{
public :
	AddChainBreak();
	virtual ~AddChainBreak();
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new AddChainBreak ); }
	void parse_my_tag( utility::tag::TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	void resnum( std::string const r ){ resnum_=r; }
	std::string resnum() const{ return resnum_;}
	void change_foldtree( bool const c ){ change_foldtree_ = c; }
	bool change_foldtree() const{ return change_foldtree_; }
	bool find_automatically() const{ return find_automatically_; }
	void find_automatically( bool const b ){ find_automatically_ = b; }
	core::Real automatic_distance_cutoff() const{ return automatic_distance_cutoff_; }
	void automatic_distance_cutoff( core::Real const a ){ automatic_distance_cutoff_ = a; }
private :
	std::string resnum_;
	bool change_foldtree_; //dflt true; should we add a jump around the chainbreak?
	bool find_automatically_; // dflt false; allow the mover to decide where the cutpoint is based on the distance between subsequent C and N atoms?
	core::Real automatic_distance_cutoff_; //dflt 2.5; very large (probably 1.5 is enough), but this a very primitive method for finding breaks!
};

} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_AddChainBreak_HH*/
