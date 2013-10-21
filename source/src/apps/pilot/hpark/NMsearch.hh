// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Pose
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/io/silent/SilentStruct.hh>

// Scoring
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>

// Normal Mode
#include <protocols/normalmode/NormalModeRelaxMover.fwd.hh>
#include <protocols/normalmode/NormalModeRelaxMover.hh>

// LoopHash
#include <protocols/loophash/LoopHashSampler.fwd.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LocalInserter.fwd.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <core/io/silent/SilentStruct.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

// Relax, min
#include <core/optimization/AtomTreeMinimizer.fwd.hh>
#include <core/optimization/CartesianMinimizer.fwd.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

//
#include <protocols/simple_moves/SwitchResidueTypeSetMover.fwd.hh>

using namespace core;

namespace myspace{
class Evaluator{

public:

Evaluator();
~Evaluator();

Evaluator( pose::Pose const &pose,
					 pose::Pose const &native );

Real apply( pose::Pose &pose ){ return scoring::CA_gdtmm( pose, native_, resmap_ ); }
Real apply( pose::Pose &pose1, pose::Pose &pose2 ) { return scoring::CA_gdtmm( pose1, pose2 ); }

void store_gdt0( pose::Pose &pose ){ gdt0 = apply( pose ); }

public:
  Real gdt0;

private:
  pose::Pose native_;
  std::map< Size, Size > resmap_;
};

////////////////////
class Scheduler{

public:

Scheduler(){}
~Scheduler(){}

Scheduler( pose::Pose &pose,
					 pose::Pose const &native,
					 std::string const mode,
					 Size const maxiter = 100 );

void
run_NM( pose::Pose &pose, std::string const mode );

void
run_NM_evalmodes( pose::Pose &pose, 
									Size const maxmode );

void set_maxiter( Size const value ){ maxiter_ = value; }

void run_LH( pose::Pose &pose, std::string const mode );

void run_LHeval( pose::Pose &pose );

void run_combine( pose::Pose &pose, bool const report_full );

void run_md( pose::Pose &pose );

private:

void 
set_default_parameters();

void
set_pertscale( Real const scale );

pose::Pose
NM_linesearch( pose::Pose const &pose, 
							 Real &gdtmin, Real &pert_best,
							 std::vector< io::silent::SilentStructOP > &ss
							 );

void 
control_schedule( pose::Pose const &pose,
									pose::Pose &pose_ref, 
									Real const gdt_to_ref );

void
setup_LHsampler( pose::Pose const &pose,
								 Size const max_struct = 20 );
 
pose::Pose
loophash_search( pose::Pose const &pose,
								 Real &gdt, Size const ires,
								 std::vector< io::silent::SilentStructOP > &ss
								 );

Size
pick_hashing_res( pose::Pose const &pose,
									std::string const mode,
									Size &Lcount );

void
put_coordinate_constraints( pose::Pose & pose );

void 
accept_and_report( pose::Pose &pose_tmp,
									 pose::Pose &pose,
									 Size const iter,
									 Size const arg1,
									 Real const arg2,
									 Real const gdt,
									 Real const gdt_to_ref,
									 bool const is_lh,
									 Real const pdb_dump_cut = 0.005 );

private:
  myspace::Evaluator evaluator_;
  protocols::loophash::LoopHashSamplerOP LHsampler_;
  protocols::normalmode::NormalModeRelaxMoverOP NMmover_;

  // On-the-fly values
  Real gdtmax_;

  // Parameters
  Size iter_updated_last_;
  Size iter_added_last_;
  Size maxmode_;
	Size maxiter_;
  Real max_scale_;
  Real min_scale_;
	Size loopcut_;
  Real gdt_refresh_cut_;
  
  // LH options
  Size looplen_;

  // Pert
  utility::vector1< Real > pert_scales_;
};

} //myspace

