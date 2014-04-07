////////////////////////////////////////////////////////////////////////////////
//     pathrover.hh: module for computing pathways between two conformations
//                   based on various criteria, either for kinetic or sampling
//                   applications
//
//     Typical Usage with initialization from misc:
//          PathRover pw(); // TODO: direct initialization from PDB
//          pw.run();
//          fp.dump_info(); // TODO: implement this?
//
//     Originally Created in R++: 29/11/2007 (Barak Raveh & Angela Enosh)
//     Transformed to mini: 27/10/2009 (Barak Raveh)
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_devel_path_rover_pathrover_hh
#define INCLUDED_devel_path_rover_pathrover_hh

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include "pathrover_algorithms.hh"
#include "pathrover_DOFs_manager.hh"
#include "pathrover_parameters.hh"
#include "pathrover_partial_data_manager.hh"

# #include <ObjexxFCL/ObjexxFCL.hh>
# #include "loop_class.h"
# #include "score_data.h"

namespace protocols{
  namespace pathrover{

    // wrapping function for invoke pathrover protocol for generating conformation pathrover
    // invoked from job_distributor::run_jobs() or something similar
    void pathrover_generator_main(bool &failed);

    /* Main class for generating pathways of conformations */
    class PathRover {
      private:
      ////////////////////
      // class properties:
      ////////////////////
      core::pose::PoseOP  _src_pose, _trg_pose, _trg_pose_partial,_src_pose_origin, _trg_pose_origin; // pose of the strucutres
      core::scoring::ScoreFunctionOP _score_weight_map; // current weight map
      Params_handler _params_h;
      Partial_Data _partial_data;
      bool _full_atom;
      bool _is_fail; // flag for failures in class operation

      // flags of fold tree:
      // ===================
      DOFs_manager _dofs_manager;

      // flags for building pose:
      // =======================
      bool ideal_pos;
      bool coords_init;
      bool full_atom;	 // flag for building pose - fullatom vs. centroid mode

    public:

      // ctr - loads pose from misc + loop list
      PathRover();

      // start generating pathways
      void run();

      // check fail flag
      bool fail()
      { return _is_fail; }

      Partial_Data& get_partial_data() {return _partial_data;}

      DOFs_manager& get_dofs_manager()
      { return _dofs_manager; }

      /* Getters & Setters: */

      core::pose::PoseOP& get_src_pose(){return _src_pose;}

      core::pose::PoseOP& get_trg_pose(){return _trg_pose;}


      core::pose::PoseOP& get_trg_partial_data_pose(){return _trg_pose_partial;}

      core::scoring::ScoreFunction& get_score_weight_map(){return _score_weight_map;} // current weight map

      // set fail flag to false
      void resetFail()
      { _is_fail = false; }

      Params_handler& get_params_h(){ return _params_h; }

      bool is_full_atom(){ return _full_atom; }

    private:

      // create a pose, a fold tree and a DOFs manager
      // according to parameters
      void initialize_from_params(Params_handler& ph);

    }; // class PathRover

    /* function wrapping the macro */
    // NOTE: this is just a small utility stuff to wrap the macro "toupper",
    // so that it can be passed as function ptr to std::transform().
    // TODO: move somewhere else, to a utility header of some form
    char toupper_wrapper (const char c);

    // auxilary class for storing chain residue range
    struct Chain_boundaries{
      int first, last;
      bool valid;
      // ctrs
      Chain_boundaries() :
	valid(false) {}
      Chain_boundaries(int _first, int _last) :
	first(_first), last(_last), valid(true) {}
    };



} // namespace pathways
} // namespace protocols

#endif
