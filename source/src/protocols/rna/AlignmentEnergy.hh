// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/AlignmentEnergy.hh
/// @brief An RMSD-type energy to a reference pose, complete with derivatives hacked out of coordinate constraints
/// @author Andy Watkins (amw579@nyu.edu)

#ifndef INCLUDED_protocols_rna_AlignmentEnergy_hh
#define INCLUDED_protocols_rna_AlignmentEnergy_hh

// Unit headers
#include <protocols/rna/AlignmentEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>

// Project Headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <protocols/stepwise/modeler/align/StepWisePoseAligner.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <utility/vector1.hh>
#include <core/types.hh>

namespace protocols {
namespace rna {

///@brief An RMSD-type energy to a reference pose, complete with derivatives hacked out of coordinate constraints
class AlignmentEnergy : public core::scoring::methods::WholeStructureEnergy {

	typedef core::scoring::methods::WholeStructureEnergy parent;

public:

	//AlignmentEnergy();

	AlignmentEnergy(
		core::scoring::methods::EnergyMethodOptions const & options //core::pose::PoseOP const & align_pose,
		//utility::vector1< Size > const & moving_res_list
	);

	// copy constructor (not needed unless you need deep copies)
	//AlignmentEnergy( AlignmentEnergy const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~AlignmentEnergy() {};

	/// @brief Clone: create a copy of this object, and return an owning pointer
	/// to the copy.
	virtual core::scoring::methods::EnergyMethodOP clone() const override;

	/// @brief Indicate required setup steps for scoring
	virtual
	void setup_for_scoring( core::pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	/// @brief Is the score context dependent or context independent?
	virtual void indicate_required_context_graphs( utility::vector1< bool > &context_graphs_required ) const override;

	/// @brief Indicates the current version of this score term
	virtual core::Size version() const override;

	/// @brief Actually calculate the total energy
	/// @details Called by the scoring machinery.
	virtual void finalize_total_energy( core::pose::Pose & pose, core::scoring::ScoreFunction const &, core::scoring::EnergyMap & totals ) const override;

	virtual
	void
	eval_atom_derivative(
		core::id::AtomID const & id,
		core::pose::Pose const & pose,
		core::kinematics::DomainMap const & domain_map,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & emap,
		core::Vector & F1,
		core::Vector & F2
	) const override;

	void align_pose( core::pose::PoseOP const & align_pose );

	protocols::stepwise::modeler::align::StepWisePoseAlignerOP const &
	pose_aligner() const { return pose_aligner_; }

	void rmsd_screen( core::Real const setting );

private:

	core::pose::PoseOP align_pose_;
	utility::vector1< Size > moving_res_list_;
	core::scoring::func::FuncOP func_;
	protocols::stepwise::modeler::align::StepWisePoseAlignerOP pose_aligner_;
	core::Real rmsd_screen_;

};

} //protocols
} //rna

#endif //protocols/rna_AlignmentEnergy_hh
