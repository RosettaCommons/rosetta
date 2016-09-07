// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/enzdes/SecondaryMatchProtocol.hh
///
/// @brief
/// @author Florian Richter


#ifndef INCLUDED_protocols_enzdes_SecondaryMatchProtocol_hh
#define INCLUDED_protocols_enzdes_SecondaryMatchProtocol_hh


#include <protocols/enzdes/EnzdesBaseProtocol.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace enzdes {

class SecondaryMatchProtocol;
typedef utility::pointer::shared_ptr< SecondaryMatchProtocol > SecondaryMatchProtocolOP;
typedef utility::pointer::shared_ptr< SecondaryMatchProtocol const > SecondaryMatchProtocolCOP;
typedef utility::pointer::weak_ptr< SecondaryMatchProtocol const > SecondaryMatchProtocolCAP;

class PoseFoundResiduesCombination;

typedef utility::pointer::shared_ptr< PoseFoundResiduesCombination > PoseFoundResiduesCombinationOP;


class SecondaryMatchProtocol : public protocols::enzdes::EnzdesBaseProtocol
{

public:

	SecondaryMatchProtocol();
	~SecondaryMatchProtocol() override;

	void apply( core::pose::Pose & start_pose) override;

	std::string get_name() const override;

	core::Size
	residues_compatible(
		core::conformation::ResidueCOP res1,
		core::conformation::ResidueCOP res2
	) const;


	bool
	do_matching(
		core::pose::Pose & start_pose
	);

	void
	set_trial_positions(
		utility::vector1< core::Size > const & trial_pos ){
		trial_positions_ = trial_pos; }


protected:

	void
	add_enz_cst_interaction_to_pose(
		core::pose::Pose & pose,
		toolbox::match_enzdes_util::EnzConstraintParametersCOP params,
		toolbox::match_enzdes_util::EnzCstTemplateResCOP missing_template,
		toolbox::match_enzdes_util::EnzCstTemplateResCOP present_template,
		toolbox::match_enzdes_util::EnzConstraintIOCOP cstio
	);

	void
	find_all_allowed_positions(
		core::pose::Pose const & pose
	);

	bool
	generate_and_dump_pose_found_residues_combinations( core::pose::PoseCOP ref_poseCOP );


	bool
	restype_possible_at_position(
		core::pose::Pose const & pose,
		core::chemical::ResidueTypeCOP restype,
		core::conformation::ResidueCOP target_residue,
		core::Size const trial_pos
	);

	void
	determine_found_residues_compatibility( core::pose::PoseCOP ref_poseCOP );


private:
	utility::vector1< utility::vector1< core::conformation::ResidueOP > > found_resis_;

	std::map< core::conformation::ResidueCOP, std::map< core::conformation::ResidueCOP, core::Size > > found_res_compatibility_;
	bool found_res_compatibility_determined_;

	utility::vector1< core::Size >trial_positions_;
	core::scoring::ScoreFunctionOP reduced_scofx_;
	core::Size cut1_, cut2_, cut3_, cut4_;
	utility::vector1< toolbox::match_enzdes_util::EnzConstraintParametersCOP > match_params_;

}; //class SecondaryMatchProtocol


/// @brief helper class to process and output the different found poses
class PoseFoundResiduesCombination : public utility::pointer::ReferenceCount
{
public:
	~PoseFoundResiduesCombination() override;
	PoseFoundResiduesCombination(
		core::pose::PoseCOP ref_pose_in,
		SecondaryMatchProtocolCAP seqmatch_in
	);

	void
	add_residue( core::conformation::ResidueOP res_in );

	bool
	construct_and_dump_outpose( utility::vector1< toolbox::match_enzdes_util::EnzConstraintParametersCOP > match_params );

private:

	//the reference pose for this class
	core::pose::PoseCOP ref_pose_;

	//the residues that need to be put into the ref pose
	utility::vector1< core::conformation::ResidueCOP > combine_resis_;

	//the secondary match protocol that a particular instance comes from
	SecondaryMatchProtocolCAP secmatch_prot_;

}; //PoseFoundResidueCombination


} //namespace enzdes
} //namespace protocols


#endif // INCLUDED_protocols_enzdes_EnzdesFixBBProtocol_HH
