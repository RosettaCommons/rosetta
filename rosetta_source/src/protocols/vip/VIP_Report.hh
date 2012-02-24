// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/vip/VIP_Report.hh
/// @brief class for reports from vip mover implementation


#ifndef INCLUDED_protocols_vip_VIP_Report_HH
#define INCLUDED_protocols_vip_VIP_Report_HH


#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/AddCavitiesMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/Energies.hh>


#include <core/init.hh>
#include <core/types.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/packstat/types.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/packstat/packing_score_params.hh>
#include <core/scoring/packstat/AtomRadiusMap.hh>
#include <core/scoring/packstat/SimplePDB.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <core/conformation/Residue.hh>
#include <numeric/all.hh>
#include <utility/all.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/MoverContainer.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
//#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cp.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <utility/options/keys/OptionKey.hh>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <utility/io/izstream.hh>
#include <ObjexxFCL/format.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace vip {

class VIP_Report
{

	public:
                VIP_Report();
                virtual ~VIP_Report();

//		void define_report_file();
//		void close_report_file();
		void get_GOE_repack_report( core::pose::Pose & goe_native, utility::vector1<core::Real> & goe_repack_e, utility::vector1<core::conformation::ResidueOP> & goe_repack_res, utility::vector1<core::Size> & goe_repack_pos );
		void get_GOE_relaxed_report( core::pose::Pose & goe_native, utility::vector1<core::Real> & goe_repack_e, utility::vector1<core::conformation::ResidueOP> & goe_repack_res, utility::vector1<core::Size> & goe_repack_pos );
		void get_GOE_packstat_report(   core::pose::Pose & goe_native, utility::vector1<core::pose::Pose> & goe_relax );

//		friend class VIP_Mover;
};

}}
#endif

