// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/SaveAndRetrieveSidechains.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_SaveAndRetrieveSidechains_hh
#define INCLUDED_protocols_protein_interface_design_movers_SaveAndRetrieveSidechains_hh
#include <protocols/protein_interface_design/movers/SaveAndRetrieveSidechains.fwd.hh>


#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/datacache/DataMapObj.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>



namespace protocols {
namespace protein_interface_design {
namespace movers {


/// @brief saves a pose and reapplies its sequence and rotamers at a later stage.
/// The constructor saves the initial pose, and then any calls to apply replace the residues on the input pose with
/// that saved pose. Notice, that only ALA positions will be replaced, so this is meant to work strictly along with
/// BuildAlaPose moves. This way, if in the design process an interface residue is designed, that will not be reverted
/// to w/t
class SaveAndRetrieveSidechains : public simple_moves::DesignRepackMover
{
public:
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
public:
	SaveAndRetrieveSidechains();
	SaveAndRetrieveSidechains(
		core::pose::Pose const & pose,
		bool const allsc=false,
		bool const ensure_variant_matching=false,
		core::Size const jumpid=1
	);
	virtual ~SaveAndRetrieveSidechains();
	bool allsc() const { return allsc_; }
	void allsc( bool const allsc ) { allsc_ = allsc; }
	void apply( Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new SaveAndRetrieveSidechains ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	bool two_step() const{ return two_step_; }
	void two_step( bool const b ) { two_step_ = b; }

	bool multi_use() const{ return multi_use_; }
	void multi_use( bool const b ) { multi_use_ = b; }
private:
	PoseOP init_pose_;
	bool allsc_, ensure_variant_matching_, two_step_,multi_use_; // two_step: dflt false; on first apply, record sidechains, on second apply, enforce them.
	core::Size jumpid_;
	utility::pointer::owning_ptr< basic::datacache::DataMapObj< bool > > first_apply_; // internal for two_step_
};

} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_SaveAndRetrieveSidechains_HH*/
