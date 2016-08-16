// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file IRCollection.hh
/// @brief class declaration for collections of inverse rotamers
/// @author havranek


#ifndef INCLUDED_protocols_motifs_IRCollection_hh
#define INCLUDED_protocols_motifs_IRCollection_hh

#include <protocols/loops/Loops.fwd.hh>
#include <protocols/motifs/Motif.fwd.hh>
#include <protocols/motifs/MotifLibrary.fwd.hh>
#include <protocols/motifs/IRCollection.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>

#include <map>

//Auto Headers
namespace protocols {
namespace motifs {

class IRCollection{

public:
	IRCollection();
	IRCollection( core::pose::Pose & pose, MotifLibrary & motifs, utility::vector1< core::Size > const & build_sites );

	core::Size nirotamers() const;

	void find_closest_backbone( core::pose::Pose & pose, protocols::loops::LoopsOP const flexible_positions,
		utility::vector1< core::Size > & closest_pos, utility::vector1< core::Real > & closest_rmsd );

	void incorporate_motifs( core::pose::Pose & pose, protocols::loops::LoopsOP const flexible_positions, utility::vector1< core::Size > & trim_positions );

	void try_for_more( core::pose::Pose & pose,
		protocols::loops::LoopsOP const flexible_positions,
		std::map< core::Size, MotifCOP > setpos,
		std::map< core::Size, core::conformation::ResidueCOP > setpos_ir,
		std::map< core::Size, bool > setpos_forward_info,
		core::Size start_depth
	);

	bool
	successful_loop_closure(
		core::pose::Pose & pose,
		protocols::loops::LoopsOP flexible_regions,
		std::map< core::Size, MotifCOP > & setpos,
		std::map< core::Size, core::conformation::ResidueCOP > & setpos_ir,
		std::map< core::Size, bool > & setpos_forward_info,
		core::conformation::Residue const & this_rotamer,
		core::Size const this_pos,
		MotifCOP this_motif,
		bool const this_forward_info
	);

	std::string
	make_motif_filename(
		std::map< core::Size, MotifCOP > & setpos,
		std::map< core::Size, bool > & setpos_forward_info,
		core::pose::Pose & pose
	);

	core::Size unique_id() { return unique_id_; };
	void reset_unique_id() { unique_id_ = 0; };
	void increment_unique_id() { ++unique_id_; };

private:
	core::Size unique_id_;
	MotifCOPs motif_source_;
	utility::vector1< bool > motif_forward_;
	utility::vector1< core::Size > target_positions_;
	utility::vector1< core::pack::rotamer_set::RotamerSetOP > rotamer_sets_;
};

}
}

#endif // INCLUDED_protocols_motifs_IRCollection
