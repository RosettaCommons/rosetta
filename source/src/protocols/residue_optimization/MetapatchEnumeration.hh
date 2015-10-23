// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/residue_optimization/MetapatchEnumeration.hh
/// @brief a class that enumerates all metapatched variants
/// @author Andy Watkins (amw579@nyu.edu)

#ifndef INCLUDED_protocols_residue_optimization_MetapatchEnumeration_hh
#define INCLUDED_protocols_residue_optimization_MetapatchEnumeration_hh

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/Metapatch.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/pointer/owning_ptr.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <utility/vector1.functions.hh>

// Utility Headers
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

// C++ headers
#include <string>
#include <sstream>

namespace protocols {
namespace residue_optimization {

using namespace core;
using namespace kinematics;
using namespace chemical;
using namespace scoring;
using namespace pose;

class MetapatchEnumeration {
public:
	MetapatchEnumeration( Pose pose, MoveMapOP mm, ScoreFunctionOP score_fxn, Size mode, bool pack, Real threshold, utility::vector1< std::string > metapatch_names, Real taboo = 5 ) :
		pose_( pose ),
		mm_( mm ),
		sampling_score_fxn_( score_fxn ),
		evaluation_score_fxn_( score_fxn ),
		mode_( mode ),
		pack_( pack ),
		threshold_( threshold ),
		metapatch_names_( metapatch_names ),
		taboo_( taboo )
	{}

	utility::vector1< Real > final_scores() { return final_scores_; }
	utility::vector1< std::string > summary_lines() { return summary_lines_; }

	void generate_metapatched_variants( Size resi );

	Real score() { return score_; }

private:

	void generate_derived_types( Size resi, Size tp, utility::vector1< std::string > & types_considered );
	Real binding( Pose pose );
	void initial_sampling( Pose & pose );
	void final_sampling( Pose & mut_pose, Size resi );
	bool tabooed( utility::vector1< std::string > const & patch_names );

	utility::vector1< Real > final_scores_;
	utility::vector1< std::string > summary_lines_;

	Pose pose_;
	MoveMapOP mm_;
	ScoreFunctionOP sampling_score_fxn_;
	ScoreFunctionOP evaluation_score_fxn_;
	Size mode_;
	bool pack_;
	Real threshold_;
	utility::vector1< std::string > metapatch_names_;

	Real score_;
	// Score past which a mutation is no longer investigated
	Real taboo_;
	utility::vector1< utility::vector1< std::string > > tabooed_;

};

}
}

#endif
