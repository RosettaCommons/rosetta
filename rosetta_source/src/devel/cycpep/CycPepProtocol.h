/*
 * CycPepProtocol.h
 *
 *  Created on: Nov 18, 2009
 *      Author: liorz06
 */

#ifndef CYCPEPPROTOCOL_H_
#define CYCPEPPROTOCOL_H_
#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/loops/LoopRelaxMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/PackRotamersMover.fwd.hh>
#include <core/pack/task/operation/TaskOperations.fwd.hh>

#include <devel/FlexPepDocking/FlexPepDockingFlags.fwd.hh>
#include <devel/FlexPepDocking/FlexPepDockingPoseMetrics.hh>


//Auto Headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
namespace protocols {
namespace CycPep {
class CycPepProtocol : public moves::Mover {
public:
	CycPepProtocol();
	virtual ~CycPepProtocol();
};
}
}
#endif /* CYCPEPPROTOCOL_H_ */
