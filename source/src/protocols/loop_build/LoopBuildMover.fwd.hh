/*
 * LoopBuilderMover.fwd.hh
 *
 *  Created on: Jun 9, 2012
 *      Author: nannemdp
 */

#ifndef INCLUDED_protocols_loop_build_LoopBuildMover_fwd_hh
#define INCLUDED_protocols_loop_build_LoopBuildMover_fwd_hh
#include <utility/pointer/owning_ptr.hh>


namespace protocols{
namespace loop_build{

class LoopBuildMover;


typedef utility::pointer::shared_ptr< LoopBuildMover > LoopBuildMoverOP;
typedef utility::pointer::shared_ptr< LoopBuildMover const > LoopBuildMoverCOP;



}//loop_build
}//protocols

#endif /* LOOPBUILDERMOVER_FWD_HH_ */
