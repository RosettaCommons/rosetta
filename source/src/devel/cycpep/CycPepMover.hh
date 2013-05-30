/*
 * CycPepProtocol.h
 *
 *  Created on: Nov 18, 2009
 *      Author: Lior Zimmerman
 */

#ifndef CYCPEPPROTOCOL_H_
#define CYCPEPPROTOCOL_H_
#include <numeric/constants.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>
#include <core/conformation/util.hh>
#include <sstream>
#include <core/graph/Graph.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/init/init.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/util.hh>//option.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.fwd.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <core/kinematics/FoldTree.hh>
#include <numeric/random/random.hh>
#include <protocols/loops/Loops.hh>
#include <core/chemical/VariantType.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/util/disulfide_util.hh>
#include <devel/init.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <utility/vector1.hh>
#include <basic/basic.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/disulfides/FullatomDisulfidePotential.hh>
using namespace core;
namespace protocols {
namespace CycPepMover {
class CycPepMover;
typedef utility::pointer::owning_ptr< CycPepMover > CycPepMoverOP;
typedef utility::pointer::owning_ptr< CycPepMover const > CycPepMoverCOP;
class CycPepMover : public moves::Mover {
public:
	CycPepMover();
	virtual ~CycPepMover();
	void minimize(pose::Pose& workpose);
	void packRotamers(pose::Pose& workpose);
	void apply(pose::Pose& workpose);
	virtual std::string get_name() const;
private:
	core::scoring::constraints::ConstraintSetOP prolineConstraint(pose::Pose& workpose);
	core::scoring::constraints::ConstraintSetOP IleConstraint(pose::Pose& workpose);
	void printEnergies(pose::Pose& workpose);
	void rotateUntilCys(pose::Pose& workpose, Size untilCys);
	void updateSSAtoms(pose::Pose& workpose, Size upNum, utility::vector1_int& vec, Size pepsize);
	Real scoreNoConstraint(pose::Pose& workpose);
	core::scoring::constraints::ConstraintSetOP createDihedralConstraint(pose::Pose& workpose);
	void modelSSLoop(Size startCys, Size endCys, pose::Pose& workpose);
	protocols::comparative_modeling::LoopRelaxMover _loop_relax_mover;
	core::scoring::ScoreFunctionOP _scorefxn;
	pack::task::PackerTaskOP _packTask;
	//prevents repacking of disulfides bonds
	pack::task::operation::NoRepackDisulfides _noRepackDisulf;
	//prevents the repacking of N-C terminal (repacking of these residues disconnects the structures)
	pack::task::operation::PreventRepacking _preventer;
	core::pack::task::TaskFactoryOP _packfactory;
};
}
}
#endif /* CYCPEPPROTOCOL_H_ */
