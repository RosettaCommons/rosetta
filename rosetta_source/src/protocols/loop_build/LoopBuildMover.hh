/*
 * LoopBuildMover.hh
 *
 *  Created on: Jun 9, 2012
 *      Author: nannemdp
 */

#ifndef INCLUDED_protocols_loop_build_LoopBuildMover_hh
#define INCLUDED_protocols_loop_build_LoopBuildMover_hh

#include <core/pose/Pose.hh>

#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <protocols/loop_build/LoopBuildMover.fwd.hh>
#include <protocols/moves/Mover.hh>

namespace protocols{
namespace loop_build{


class LoopBuildMover : public moves::Mover {
public:

	LoopBuildMover(protocols::comparative_modeling::LoopRelaxMover loop_relax_mover);

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

private:
	protocols::comparative_modeling::LoopRelaxMover loop_relax_mover_;

private:
	void setup_loop_definition();
};


} // namespace loop_build
} // namespace protocols

#endif /* INCLUDED_protocols_loop_build_LoopBuildMover_hh */
