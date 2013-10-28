// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SewingDesigner.hh
///
/// @brief A Mover that uses loophash to find fragments that can
/// bridge a gap with minimal modifications to the original pose.
/// @author Tim Jacobs


#ifndef INCLUDED_devel_sewing_SewingDesigner_HH
#define INCLUDED_devel_sewing_SewingDesigner_HH

//Core
#include <core/pose/Pose.hh>

//Protocols
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>
#include <devel/loop_creation/LoopCreationMover.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/features/FeaturesReporter.fwd.hh>

//Utility
#include <utility/vector1.hh>

//Devel
#include <devel/sewing/util.hh>

//C++
#include <set>

namespace devel {
namespace sewing {
	
struct SewingDesignerData
{
	core::pose::Pose pose;
	protocols::loops::Loops loops;
	NativeRotamersMap native_helix_residues;
	std::map<core::Size, char> native_loop_residues;
};

class SewingDesigner : public protocols::moves::Mover {
	
public:

	SewingDesigner();
	
	SewingDesigner(
		NativeRotamersMap native_helix_residues
	);
	
	SewingDesigner(
		devel::loop_creation::LoopCreationMoverOP loop_creation_mover,
		NativeRotamersMap native_helix_residues,
		utility::vector1<core::Size> loop_anchors
	);
	
	protocols::moves::MoverOP
	clone() const;

	protocols::moves::MoverOP
	fresh_instance() const;
	
	std::string
	get_name() const;

	void
	init();

	void
	apply( Pose & pose );
	
	void
	detect_loop_anchors(
		Pose const & pose
	);
	
	void
	enforce_repeat(
		core::pose::Pose & pose,
		protocols::loops::Loops & loops,
		NativeRotamersMap & nat_ro_map
	);
	
	void
	design(
		core::pose::Pose & pose,
		NativeRotamersMap const & native_helix_residues
	);
	
	void
	design_neighborhood(
		core::pose::Pose & pose,
		NativeRotamersMap const & nat_ro_map,
		std::set<core::Size> const & neighbor_residues
	);
	
	void
	record_statistics(
		core::pose::Pose & pose,
		NativeRotamersMap const & native_helix_residues,
		std::map<core::Size, char> const & native_loop_residues
	);
	
	void
	update_anchors(
		utility::vector1<core::Size> & loop_anchors,
		protocols::loops::Loop const & new_loop,
		core::Size index_of_new_loop
	);
	
	void
	update_native_helix_residues(
		NativeRotamersMap & native_helix_residues,
		protocols::loops::Loop const & new_loop
	);
	
	void
	rearrange_pose(
		Pose & pose,
		NativeRotamersMap & native_helix_residues,
		utility::vector1<core::Size> & loop_anchors
	);
	
	void
	parse_my_tag(
		TagCOP const tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const &movers,
		Pose const & pose
	);
	
private:
	core::scoring::ScoreFunctionOP scorefxn_;
	protocols::simple_moves::PackRotamersMoverOP repack_mover_;
	
	devel::loop_creation::LoopCreationMoverOP loop_creation_mover_;
	
	NativeRotamersMap native_helix_residues_;
	utility::vector1<core::Size> loop_anchors_;
	
	core::Size asym_size_;
	
//	core::Size num_helices_in_repeat_;
};
	
} // sewing
} // devel

#endif
