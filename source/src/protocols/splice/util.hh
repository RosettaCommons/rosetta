// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Utility functions useful in Splice mover.
/// @author Gideon Lapidoth (glapidoth@gmail.com)


#ifndef INCLUDED_protocols_splice_util_hh
#define INCLUDED_protocols_splice_util_hh

#include <protocols/splice/Splice.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <utility/vector1.hh>
#include <protocols/splice/SampleRotamersFromPDB.fwd.hh>


// C++ headers
#include <string>
#include <set>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>
namespace protocols {
namespace splice {
core::Size nearest_to_entry_stop_on_pose (core::pose::Pose const & pose,core::pose::Pose const & template_pose,core::Size residue,std::string tail_segment, std::string protein_family,core::Size chain,std::string segment);
void modify_dbase_with_compatible_backbones(utility::vector1< ResidueBBDofs > torsion_database, utility::vector1< core::Size > & dbase,utility::vector1 <std::string> intersect_pdbs);
void load_pdb_segments_from_pose_comments(core::pose::Pose const & pose, std::map< std::string/*which segment (L1,L2...)*/, std::string/*pdb name*/ > & pdb_segments);
//@brief concatenate 2 res matrix objects (see SampleRotamerFromPDB)
core::Size report_coordinate_constraints(core::pose::Pose const & pose);
bool calculate_rmsd(core::pose::Pose  & pose, core::pose::Pose const & source_pose,core::Size total_residue_new,
	core::Size nearest_to_from, core::Size from_res, core::Real rms_A, core::Real rms_B=-1  );
void min_seg(core::pose::Pose & pose, ResidueBBDofs dofs,bool debug, core::Size from_res, std::string tail_segment, core::Size cut_site, core::Size cut_vl_vh_after_llc,ResidueBBDofs tail_dofs,core::scoring::ScoreFunctionOP scorefxn,std::string segment_type);
utility::vector1< numeric::xyzVector< core::Real > > coords( core::pose::Pose const & pose, utility::vector1< core::Size > const positions );
utility::vector1<core::Size> find_residues_on_chain1_inside_interface(core::pose::Pose const & pose, core::Size chainNum);
std::string parse_pdb_code(std::string pdb_file_name);
void report_designable_packable_residues(core::pack::task::TaskFactoryOP const tf, core::pose::Pose pose);
void fix_chain_break_residue(core::pose::Pose & pose,SpliceManager splicemanager );
void report_des_pack_res(core::pose::Pose const & pose,core::pack::task::TaskFactoryOP tf);
} // splice
} // protocols


#endif /*INCLUDED_protocols_splice_util_hh*/
