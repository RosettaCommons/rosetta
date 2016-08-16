// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/brunette/minimalCstRelaxUtil.hh
/// @author TJ Brunette (tjbrunette@gmail.com)

#ifndef MINIMAL_CST_RELAX_UTIL
#define MINIMAL_CST_RELAX_UTIL

#include <core/types.hh>
#include <set>
#include <map>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <utility/vector1.hh>


namespace devel {
namespace cstEnergyBalance {

using core::pose::Pose;
using core::Size;
using core::Real;
using std::string;
using namespace core::sequence;

/// @brief Gets the x,y,z center of mass of a pose
numeric::xyzVector< core::Real > get_centerOfMass(const core::pose::Pose& pose );
/// @brief Generates a set of coordinate constraints that correspond to the
/// CA that are being constrained.
core::scoring::constraints::ConstraintSetOP convert_caAtomsToConstrain_to_coordCsts(std::set< core::Size > caAtomsToConstrain,const core::pose::Pose& pose);

/// @brief Generates a list of residues to constrain based on the center of the 3 residues clossest to the center of mass of contiguous regions
std::set<Size> get_coreDistDeviationResiduesToConstrain(const Real distDeviationThresh, Pose& relaxed_pose, const Pose& unmodified_pose);

/// @brief Generates a list of core residues to constrain. The residues chosen are the number of residues to constrain closest to the center of mass
std::set<core::Size> get_coreResiduesToConstrain(core::Size numbResiduesToConstrain,const core::pose::Pose& pose);

/// @brief Generates a list of residues to constrain based on movement between relaxed pose and an unmodified pose.
std::set<core::Size> get_distDeviationResiduesToConstrain(const core::Real distDeviationThresh,core::Size gapBtwConstrainedResidues, core::pose::Pose& relaxed_pose, const core::pose::Pose& unmodified_pose);

/// @brief This util uses the above utilities to create a list of which CA residues to constrain. *Pose is not const because gaps of more than 3.5 angstroms are changed to be jumps
std::set< core::Size> get_residuesToConstrain(const core::Size coordinate_cst_gap, const Real gdt_thresh, core::pose::Pose& pose);

/// @brief Same as above but with defaults for gap size and gdt_thresh
std::set<core::Size> get_residuesToConstrain(core::pose::Pose& pose);

/// @brief Outputs coordinate contraints
void output_coordCsts(const std::set< core::Size > & caAtomsToConstrain,std::ostream & out, core::pose::Pose& pose);

/// @brief Outputs coordinate contraints
void output_coordCsts(const std::set< Size > & caAtomsToConstrain,std::ostream & out, Pose& pose, const SequenceAlignment aln,string query_sequence, bool only_res_out);

/// @brief input coordinate constraints from file in standard coordinate constraint format. I assume everything is CA constraints
std::set< core::Size > input_coordCsts(const std::string inputFileName);

/// @brief outputs Fasta with virtual atoms attached.
void output_fastaWVirtual(const std::string & fastaFileName, std::ostream & out);

/// @brief inputs the sequence alignments and maps them to the appropriate
/// pdbid
std::map< std::string,SequenceAlignment> input_alignmentsMapped(bool mapToPdbid);

/// @brief inputs the sequence alignments and keeps them ordered the way they
/// were in the alingment.filt file.
utility::vector1<SequenceAlignment> input_alignments();
/// @brief gets the first and last residues from an alignment ignoring residues that are gapped alnIdx can be 1 or 2 depending if you want the target or template.
void get_terminal_aln_res(const SequenceAlignment aln, const Size alnIdx, Size & firstRes, Size & lastRes);
} //namespace cstEnergyBalance
}//namespace devel
#endif //MINIMAL_CST_RELAX_UTIL
