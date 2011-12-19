#ifndef INCLUDED_protocols_vip_VIP_Mover_HH
#define INCLUDED_protocols_vip_VIP_Mover_HH

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
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
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cp.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <utility/options/keys/OptionKey.hh>
#include <protocols/simple_moves/AddCavitiesMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <utility/io/izstream.hh>
#include <ObjexxFCL/format.hh>
#include <utility/string_util.hh>

#include <protocols/vip/VIP_Report.hh>

namespace protocols {
namespace vip {

class VIP_Mover
{

                core::pose::Pose initial_pose;
		core::pose::Pose cavity_pose;
		core::pose::Pose final_pose;
		utility::vector1<core::pose::Pose> temp_poses;
		utility::vector1<core::pose::Pose> favorable_poses;
		core::Size number_cavities;
		utility::vector1<core::Size> cavity_balls;
		utility::vector1<std::string> favorable_mutations;
		utility::vector1<core::Size> void_neighbors;
		utility::vector1<core::Size> void_mutatables;
		core::Real final_energy;

	public:
		VIP_Mover();
		VIP_Mover(
                        core::pose::Pose,
                        core::pose::Pose,
                        core::pose::Pose,
                        utility::vector1<core::pose::Pose>,
                        utility::vector1<core::pose::Pose>,
                        core::Size,
                        utility::vector1<core::Size>,
                        utility::vector1<std::string>,
                        utility::vector1<core::Size>,
                        utility::vector1<core::Size>,
			core::Real);
		virtual ~VIP_Mover();

		void set_initial_pose( core::pose::Pose );
		void minimize_conformation();
		void compute_number_cavities();
		void get_cavity_positions();
		void apply_holes();
		void dump_pdb_to_file( core::pose::Pose, std::string );
		void get_neighbors();
		void try_point_mutants();
		void relax_favorable_poses();
		void cull_mutatable_residues();
		void sort_fill_energies();
		void print_favorable_mutations();
		void sort_relaxed_poses();
		void print_pack_report();
		void print_relax_report();
		void nook_finder();
		void cranny_packer();

		core::Real get_final_energy(){
			return final_energy;}
		core::pose::Pose get_final_pose(){
			return final_pose;}
		void apply();


};

}}
#endif
