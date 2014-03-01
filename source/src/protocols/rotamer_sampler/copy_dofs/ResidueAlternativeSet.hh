// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/copy_dofs/ResidueAlternativeSet.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_copy_dofs_ResidueAlternativeSet_HH
#define INCLUDED_protocols_rotamer_sampler_copy_dofs_ResidueAlternativeSet_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueAlternativeSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <map>

namespace protocols {
namespace rotamer_sampler {
namespace copy_dofs {

	class ResidueAlternativeSet: public utility::pointer::ReferenceCount {

	private:
		// Can't use a default constructor
		ResidueAlternativeSet();

	public:

		//constructor
		ResidueAlternativeSet(  utility::vector1< core::pose::PoseOP > const & pose_list,
														std::map< Size, Size > const & res_map,
														Size const representative_seqpos);

		//constructor
		ResidueAlternativeSet(  utility::vector1< core::pose::PoseOP > const & pose_list,
														Size const representative_seqpos);

		//destructor
		~ResidueAlternativeSet();

	public:

		std::map < Size, Size > res_map() const{ return res_map_; }
		utility::vector1< core::pose::PoseOP > pose_list() const;
		core::pose::PoseOP pose( Size const n ) const;
		core::Size representative_seqpos() const{ return representative_seqpos_; }

	private:

		utility::vector1< core::pose::PoseOP > pose_list_;
		std::map< Size, Size > res_map_;
		Size representative_seqpos_;

	};

} //copy_dofs
} //rotamer_sampler
} //protocols

#endif
