// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file   apps/public/design/zinc1_homodimer_design.cc
/// @brief  Symmetric design of a zinc-mediated homodimer.  One zinc at the interface.  Gridsearch (or not) of perturbations to the metalsite geometry using rollmoves about axes that pass through zinc.  Metalsite constraints are in place to allow backbone minimization of the designs.
/// @author Bryan Der

#include <devel/init.hh>
#include <protocols/metal_interface/MetalSiteResidue.hh>
#include <protocols/metal_interface/ZincSiteFinder.hh>
#include <protocols/metal_interface/AddZincSiteConstraints.hh>

#include <protocols/anchored_design/InterfaceAnalyzerMover.hh>
#include <protocols/rigid/RollMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceOperation.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jobdist/Jobs.hh>

#include <core/conformation/Conformation.hh>
//#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Symmetry Headers
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh> //create symmetric homodimer from input monomer via symmetry:symmetry_definition option
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymRotamerTrialsMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/conformation/symmetry/util.hh>
//Constraint Headers
#include <core/scoring/func/Func.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/id/AtomID.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <numeric/xyz.io.hh> //print a point
#include <numeric/xyz.functions.hh>

#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

#include <string>
#include <sstream>

typedef core::pose::Pose Pose;
typedef numeric::xyzVector<core::Real> point;
typedef point axis;

using namespace core;
using basic::Warning;

static basic::Tracer TR("apps.public.design.zinc1_homodimer_design");

basic::options::IntegerOptionKey const lowres_symmetric_design_cycles( "lowres_symmetric_design_cycles" );
basic::options::IntegerOptionKey const highres_symmetric_design_cycles( "highres_symmetric_design_cycles" );
basic::options::BooleanOptionKey const gridsearch_rollmove_angles( "gridsearch_rollmove_angles" );
basic::options::RealOptionKey const fav_nat_bonus("fav_nat_bonus");
basic::options::IntegerOptionKey const nstruct_iterations( "nstruct_iterations" );


///@brief
class zinc1_homodimer_design : public protocols::moves::Mover {
public:
  zinc1_homodimer_design()
  {
  }
  zinc1_homodimer_design( Size lowres_symmetric_design_cycles, Size highres_symmetric_design_cycles, bool gridsearch_rollmove_angles)
    : lowres_symmetric_design_cycles_(lowres_symmetric_design_cycles), highres_symmetric_design_cycles_(highres_symmetric_design_cycles), gridsearch_rollmove_angles_(gridsearch_rollmove_angles)
  {
  }
  virtual ~zinc1_homodimer_design(){};

  virtual void
  apply( Pose & pose ){

    utility::file::FileName pdbfilename( pose.pdb_info()->name() );
    pdbname_base_ = pdbfilename.base();

    setup( pose );

    if(gridsearch_rollmove_angles_) {
      gridsearch_rollmove_angles( pose ); // calls design_symmetric_homodimer_interface
    }

    else {
      Pose const designable_start_pose( pose );

      for (Size n(1); n <= basic::options::option[nstruct_iterations]; n++) {
				std::stringstream ss;
				ss << n;
				std::string n_string;
				if(n < 10) { n_string = "000" + ss.str(); }
				else if (n < 100) { n_string = "00" + ss.str(); }
				else if (n < 1000) { n_string = "0" + ss.str(); }
				else { n_string = ss.str(); }


				design_symmetric_homodimer_interface( pose );

				std::string nstruct_dump_name = pdbname_base_ + "_" + n_string + ".pdb";
				pose.dump_pdb(nstruct_dump_name);

				TR << "FINISHED " << nstruct_dump_name << std::endl;

				pose = designable_start_pose;

      }
    }
    return;
  }


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////// Setup  ///////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  virtual void
  setup ( Pose & pose ) {
    //pose.dump_pdb("initial_mono_pose.pdb");

    protocols::simple_moves::symmetry::SetupForSymmetryMoverOP make_monomeric_input_pose_symmetrical = new protocols::simple_moves::symmetry::SetupForSymmetryMover(); // according to symm definition file included as an option
    make_monomeric_input_pose_symmetrical->apply( pose );

    //pose.dump_pdb("initial_symm_pose.pdb");

    setup_output_tag();
    setup_metalsite( pose );
    setup_metalsite_constraints( pose );
    setup_taskfactory( pose );
    setup_scorefunctions();
    setup_movemap_bb_chi( pose );
    setup_movers();
    setup_rollmoving( pose );

  }


  //setup
  virtual void
  setup_output_tag() {

    using protocols::jd2::JobDistributor;
    protocols::jd2::JobOP job_me( JobDistributor::get_instance()->current_job() );
    std::string me( JobDistributor::get_instance()->job_outputter()->output_name( job_me ) );
    output_tag_ = me; //output_tag_base_ + ".pdb";
  }


  //setup
  virtual void
  setup_metalsite( Pose const & pose ) {
    TR << "Parsing metalsite... " << std::endl;
    protocols::metal_interface::ZincSiteFinderOP find_zinc = new protocols::metal_interface::ZincSiteFinder();
    msr_ = find_zinc->find_zinc_site(pose);
    TR << "metalsite 1. " << msr_[1]->get_seqpos() << std::endl;
    TR << "metalsite 2. " << msr_[2]->get_seqpos() << std::endl;
    TR << "metalsite 3. " << msr_[3]->get_seqpos() << std::endl;
    TR << "metalsite 4. " << msr_[4]->get_seqpos() << std::endl;
    TR << "metalsite 5. " << msr_[5]->get_seqpos() << std::endl;
  }

  //setup
  virtual void
  setup_metalsite_constraints( Pose & pose ) {
    TR << "Adding metalsite constraints..." << std::endl;
    protocols::metal_interface::AddZincSiteConstraintsOP constraints_adder = new protocols::metal_interface::AddZincSiteConstraints( msr_ );
    //constraints_adder->add_constraints( pose );

    //Favor native residue
    using namespace core::scoring::constraints;
    Size const nres( pose.total_residue() );
    Real const favor_native_res_bonus( basic::options::option[fav_nat_bonus].value() );
    for ( Size i=1; i<= nres;  ++i ) {
      pose.add_constraint( new ResidueTypeConstraint( pose, i,  favor_native_res_bonus) );
    }

    Pose copy( pose );
    starting_fullatom_pose_ = copy;
  }



  //setup
  virtual void
  setup_taskfactory( Pose & pose ) {
    using namespace core::pack::task;
    using namespace basic::options;

    TR << "Generating task factory..." << std::endl;
    TaskFactoryOP task_factory = new TaskFactory();
    task_factory->push_back(new operation::InitializeFromCommandline());
    if ( option[ OptionKeys::packing::resfile ].user() ) { // probably will not be used
      TR << "Reading resfile..." << std::endl;
      task_factory->push_back( new operation::ReadResfile );
    }
    TR << "Restricting to interface..." << std::endl;
    task_factory->push_back(new protocols::toolbox::task_operations::RestrictToInterfaceOperation(1, 2));

    operation::PreventRepackingOP prevent_repack = new operation::PreventRepacking();
    TR << "Preventing repacking of residues ";
    for(core::Size i(1); i <= 4 /*msr_.size() - 1*/; ++i) {
      TR << msr_[i]->get_seqpos() << " ";
      prevent_repack->include_residue( msr_[i]->get_seqpos() );
    }
    TR << std::endl;
    task_factory->push_back(prevent_repack);

    TR << "Packer Task: " << *(task_factory->create_task_and_apply_taskoperations(pose));

    taskfactory_ = task_factory;
  }

  //setup
  virtual void
  setup_scorefunctions() {
    TR << "Generating scorefunction..." << std::endl;
    using namespace core::scoring;

    ScoreFunctionOP softrep_scorefxn = ScoreFunctionFactory::create_score_function( SOFT_REP_DESIGN_WTS ); // SoftRep, will become SymmetricScoreFunction when symmetry definition file is included as an option
    softrep_scorefxn->set_weight( atom_pair_constraint, 4.0 ); // 4 distances
    softrep_scorefxn->set_weight( angle_constraint, 2.0 );     // 10 angles (6 tetr + 4)
    softrep_scorefxn->set_weight( dihedral_constraint, 8.0 );  // 1 dihedral per His
    softrep_sym_scorefxn_ = new core::scoring::symmetry::SymmetricScoreFunction( *softrep_scorefxn );
    TR << "Softrep Scorefunction: " << *softrep_sym_scorefxn_ << std::endl;

    ScoreFunctionOP scorefxn = get_score_function(); // SCORE12, will become SymmetricScoreFunction when symmetry definition file is included as an option
		//scorefxn->set_weight( atom_pair_constraint, 4.0 ); // 4 distances
    //scorefxn->set_weight( angle_constraint, 2.0 );     // 10 angles (6 tetr + 4)
    //scorefxn->set_weight( dihedral_constraint, 8.0 );  // 1 dihedral per His
    scorefxn->set_weight( res_type_constraint, basic::options::option[fav_nat_bonus].value() );
    sym_scorefxn_ = new core::scoring::symmetry::SymmetricScoreFunction( *scorefxn );
    TR << "Score12 Scorefunction: " << *sym_scorefxn_ << std::endl;
    no_constraints_scorefxn_for_ddG_calc_ = new core::scoring::symmetry::SymmetricScoreFunction( *scorefxn );

    metal_scorefxn_ = new ScoreFunction;
    metal_scorefxn_->set_weight( atom_pair_constraint, 4.0 ); // 4 distances
    metal_scorefxn_->set_weight( angle_constraint, 2.0 );     // 10 angles (6 tetr + 4)
    metal_scorefxn_->set_weight( dihedral_constraint, 8.0 );  // 1 dihedral per His

    pair_scorefxn_ = new ScoreFunction;
    pair_scorefxn_->set_weight( atom_pair_constraint, 4.0 ); // 4 distances

    angle_scorefxn_ = new ScoreFunction;
    angle_scorefxn_->set_weight( angle_constraint, 2.0 );     // 10 angles (6 tetr + 4)

    dihed_scorefxn_ = new ScoreFunction;
    dihed_scorefxn_->set_weight( dihedral_constraint, 8.0 );  // 1 dihedral per His

  }


  //setup
  virtual void
  setup_movemap_bb_chi( Pose & pose ) {
    TR << "Making movemap..." << std::endl;
    sc_move_map_ = new core::kinematics::MoveMap;
    sc_move_map_->set_chi(true);
    sc_move_map_->set_bb(false);

    sc_bb_move_map_ = new core::kinematics::MoveMap;
    sc_bb_move_map_->set_chi(true);
    sc_bb_move_map_->set_bb(true);

    Size zn = msr_[5]->get_seqpos();

    sc_move_map_->set_chi(zn, false);
    sc_bb_move_map_->set_chi(zn, false);


    //     using namespace core::pose::metrics;
    //     core::util::MetricValue< std::set< Size > > mv_interface_set;
    //     pose.metric("SymmInterfaceNeighborDefinition", "interface_residues", mv_interface_set);
    //     std::set< Size > const interface_set(mv_interface_set.value());

    //     TR << "Interface min_movemap PyMol selection: sel resi ";
    //     for (std::set<Size>::const_iterator it = interface_set.begin(); it != interface_set.end(); ++it) {
    //       TR << *it << "+";
    //       sc_move_map_->set_chi( *it, true );
    //       //      sc_move_map_->set_bb(*it, true);
    //     }

  }



  //setup
  virtual void
  setup_movers() {
    using namespace protocols::simple_moves::symmetry;

    TR << "Generating sym rottrials mover..." << std::endl;
    sym_rottrials_mover_ = new SymRotamerTrialsMover;
    sym_rottrials_mover_->task_factory( taskfactory_ );
    sym_rottrials_mover_->score_function( softrep_sym_scorefxn_ );

    TR << "Generating sym pack mover..." << std::endl;
    sym_pack_mover_ = new SymPackRotamersMover;
    sym_pack_mover_->task_factory( taskfactory_ );
    sym_pack_mover_->score_function( sym_scorefxn_ );

    TR << "Generating softrep and score12 minmovers..." << std::endl;
    softrep_min_mover_ = new SymMinMover( sc_move_map_, softrep_sym_scorefxn_ /*softrep + constraints*/, "dfpmin_armijo", 0.01, true );
    sc_min_mover_ = new SymMinMover( sc_move_map_, sym_scorefxn_ /*score12 + constraints*/, "dfpmin_armijo", 0.01, true );
    sc_bb_min_mover_ = new SymMinMover( sc_bb_move_map_, sym_scorefxn_ /*score12 + constraints*/, "dfpmin_armijo", 0.01, true );

    interface_analyzer_ = new protocols::anchored_design::InterfaceAnalyzerMover( 2, false, centroid_scorefxn_for_ddG_calc_ );

  }


  //setup
  virtual void
  setup_rollmoving( Pose const & pose ) {
    //	TR << "define_rollmove_axes" << std::endl;
    //p1a and p1b are the atom xyz's of liganding residues on the target protein
    point p1a = pose.residue(msr_[1]->get_seqpos()).atom(msr_[1]->get_ligand_atom_name()).xyz(); // hard-coded for now
    point p1b = pose.residue(msr_[2]->get_seqpos()).atom(msr_[2]->get_ligand_atom_name()).xyz();
    point zinc = pose.residue(msr_[5]->get_seqpos()).atom(1).xyz();
    Size const chain_begin = pose.conformation().chain_begin(2);
    Size const chain_end =  pose.conformation().chain_end(2);
    //	TR << "Chain_Begin " << chain_begin << " Chain_End " << chain_end << std::endl;

    //create Z axis
    axis const zaxis = cross_product( zinc-p1a, zinc-p1b );
    zmover_ = new protocols::rigid::RollMover( chain_begin, chain_end, -10, 10, zaxis, zinc );
    //create Y axis
    core::Angle const create_y_rot( angle_of(zinc-p1a, zinc-p1b) );
    numeric::xyzMatrix< core::Real > z_rotation_matrix( numeric::rotation_matrix( zaxis, create_y_rot ) );
    axis const yaxis = ( z_rotation_matrix*(zinc - p1a) );
    ymover_ = new protocols::rigid::RollMover( chain_begin, chain_end, -10, 10, yaxis, zinc );
    //create X axis
    axis const xaxis = cross_product(zaxis, yaxis);
    xmover_ = new protocols::rigid::RollMover( chain_begin, chain_end, -10, 10, xaxis, zinc );
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////// End of Setup  ////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////





  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////// Protocol  ////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //protocol
  virtual void
  gridsearch_rollmove_angles( Pose & pose ) {

		interface_analyzer_->set_use_centroid_dG( true );


    Real x, y, z;
    Real ddG_centroid;
    for( z = -10.0; z <= 10; z=z+10 ) {
      zmover_->set_min_max_angles( z, z );

      for( y = -10.0; y <= 10; y=y+10 ) {
				ymover_->set_min_max_angles( y, y );

				for( x = -10.0; x <= 10; x=x+10 ) {

					std::stringstream zstr;
					zstr << z;
					std::stringstream ystr;
					ystr << y;
					std::stringstream xstr;
					xstr << x;

					xmover_->set_min_max_angles( x, x );
					TR << "Rollmoving to " << z << " " << y << " " << x << " " << std::endl;
					zmover_->apply( pose );
					ymover_->apply( pose );
					xmover_->apply( pose );

					interface_analyzer_->apply( pose );
					if( interface_analyzer_->get_centroid_dG() < 10.0 /*cutoff*/ ) {
						TR << "DESIGNABLE" << std::endl;
						design_symmetric_homodimer_interface( pose );
						interface_analyzer_->apply( pose );
						if(interface_analyzer_->get_separated_interface_energy() < -10.0) {
							std::string name = output_tag_ + "_" + zstr.str() + "_" + ystr.str() + "_" + xstr.str() + ".pdb";
							pose.dump_scored_pdb( name, *no_constraints_scorefxn_for_ddG_calc_ );
						}
					}
					else{ TR << "NOT DESIGNABLE" << std::endl; }

					pose = starting_fullatom_pose_;
				}
      }
    }
  }


  //protocol
  virtual void
  design_symmetric_homodimer_interface( Pose & pose ) {

    protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( pose , *sym_scorefxn_ , 0.6 );

    for( Size i(1); i <= lowres_symmetric_design_cycles_; ++i ) {
      TR << "Lowres design cycle " << i << " out of " << lowres_symmetric_design_cycles_ << std::endl;

      sym_rottrials_mover_->apply( pose );
      softrep_min_mover_->apply( pose );
      mc->boltzmann( pose );
    }

    for( Size i(1); i <= highres_symmetric_design_cycles_; ++i ) {
      TR << "Highres design cycle " << i << " out of " << highres_symmetric_design_cycles_ << std::endl;

      sym_pack_mover_->apply( pose );
      TR << "Score after packing " << sym_scorefxn_->score( pose ) << std::endl;
      sc_min_mover_->apply( pose );
      TR << "Score after packing and minimization" << sym_scorefxn_->score( pose ) << std::endl;
      mc->boltzmann( pose );
    }

    mc->recover_low( pose );

    //     TR << "Minimizing sc bb" << std::endl;
    //     TR << "Score before minimization " << sym_scorefxn_->score( pose ) << std::endl;
    //     sc_bb_min_mover_->apply( pose );
    //     TR << "Score after minimization 1 " << sym_scorefxn_->score( pose ) << std::endl;
    //     sc_bb_min_mover_->apply( pose );
    //     TR << "Score after minimization 2 " << sym_scorefxn_->score( pose ) << std::endl;
    //     sc_bb_min_mover_->apply( pose );
    //     TR << "Score after minimization 3"  << sym_scorefxn_->score( pose ) << std::endl;


    return;
  }



  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////// End of Protocol  /////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////



	virtual
	std::string
	get_name() const { return "zinc1_homodimer_design"; }




private:
  //options
  Size lowres_symmetric_design_cycles_;
  Size highres_symmetric_design_cycles_;
  bool gridsearch_rollmove_angles_;

  Pose starting_monomer_pose_;
  Pose centroid_; // centroid pose is rollmoved and scored during gridsearch
  Pose /*const*/ starting_centroid_pose_;
  Pose /*const*/ starting_fullatom_pose_;

  utility::vector1< protocols::metal_interface::MetalSiteResidueOP > msr_;
  core::pack::task::TaskFactoryCOP taskfactory_;
  core::scoring::symmetry::SymmetricScoreFunctionOP softrep_sym_scorefxn_;
  core::scoring::symmetry::SymmetricScoreFunctionOP sym_scorefxn_;
  core::scoring::symmetry::SymmetricScoreFunctionOP no_constraints_scorefxn_for_ddG_calc_;
  core::scoring::symmetry::SymmetricScoreFunctionOP centroid_scorefxn_for_ddG_calc_;

  core::scoring::ScoreFunctionOP metal_scorefxn_;
  core::scoring::ScoreFunctionOP pair_scorefxn_;
  core::scoring::ScoreFunctionOP angle_scorefxn_;
  core::scoring::ScoreFunctionOP dihed_scorefxn_;

  core::kinematics::FoldTree interface_scorable_tree_;
  core::kinematics::MoveMapOP sc_move_map_;
  core::kinematics::MoveMapOP sc_bb_move_map_;

  protocols::simple_moves::symmetry::SymPackRotamersMoverOP sym_pack_mover_;
  protocols::simple_moves::symmetry::SymRotamerTrialsMoverOP sym_rottrials_mover_;
  protocols::simple_moves::symmetry::SymMinMoverOP softrep_min_mover_;
  protocols::simple_moves::symmetry::SymMinMoverOP sc_min_mover_;
  protocols::simple_moves::symmetry::SymMinMoverOP sc_bb_min_mover_;
  protocols::moves::MonteCarloOP mc_;
  protocols::anchored_design::InterfaceAnalyzerMoverOP interface_analyzer_;

  protocols::rigid::RollMoverOP zmover_;
  protocols::rigid::RollMoverOP ymover_;
  protocols::rigid::RollMoverOP xmover_;


  //int nstruct_count_;
  Pose best_pose_;
  Real best_ddG_;
  std::string output_tag_;
  std::string output_tag_base_;
  std::string pdbname_;
  std::string pdbname_base_;



};

typedef utility::pointer::owning_ptr< zinc1_homodimer_design > zinc1_homodimer_designOP;


int
main( int argc, char* argv[] )
{
	try {
  using basic::options::option;
  option.add( lowres_symmetric_design_cycles, "lowres_symmetric_design_cycles" ).def(10);
  option.add( highres_symmetric_design_cycles, "highres_symmetric_design_cycles" ).def(10);
  option.add( gridsearch_rollmove_angles, "gridsearch_rollmove_angles" ).def(false);
  option.add( fav_nat_bonus, "fav_nat_bonus" ).def(1.5);
  option.add( nstruct_iterations, "nstruct_iterations" ).def(1);

  devel::init(argc, argv);
  protocols::jd2::JobDistributor::get_instance()->go(new zinc1_homodimer_design( option[lowres_symmetric_design_cycles].value(), option[highres_symmetric_design_cycles].value(), option[gridsearch_rollmove_angles].value() ));
  TR << "************************d**o**n**e**************************************" << std::endl;


  } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }


  return 0;
}

