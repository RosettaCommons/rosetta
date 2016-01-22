// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/bder/HbondZinc.cc
/// @brief
/// @details
/// @author Bryan Der

#include <devel/init.hh>
#include <protocols/metal_interface/MatchGrafter.hh>
#include <protocols/metal_interface/ZincSiteFinder.hh>
#include <protocols/metal_interface/MetalSiteResidue.hh>
#include <protocols/metal_interface/AddZincSiteConstraints.hh>

#include <core/graph/Graph.hh>
#include <core/pack/packer_neighbors.hh>

#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/backrub/BackrubMover.hh>

#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>

#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/MetricValue.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <utility/file/FileName.hh>
#include <utility/vector1.hh>

#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <numeric/conversions.hh> //degrees-radians

//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static basic::Tracer TR("apps.pilot.bder.HbondZinc");

typedef std::set< core::Size > SetSize;

using namespace core;

basic::options::FileOptionKey const scaffold_protein( "scaffold_protein" );


/// @brief
class HbondZinc : public protocols::moves::Mover {
public:
  HbondZinc()
  {
  }
  virtual ~HbondZinc(){};


  	virtual
  	void
  	apply( core::pose::Pose & pose ){

		TR << "HbondZinc protocol" << std::endl;
	  
		find_zinc(pose); //pose is match
		graft_match(pose); //pose is now grafted scaffold
		repack(pose);
		constrain_zinc(pose);
		find_neighbors(pose);
		build_rotamers(pose); //nested with hbond constraint
		
    	return;
  	}


	virtual
	void
	find_zinc(core::pose::Pose & pose) {
		
		protocols::metal_interface::ZincSiteFinderOP zinc_finder = new protocols::metal_interface::ZincSiteFinder();
		msr_ = zinc_finder->find_zinc_site(pose);
		
		for(Size i(2); i<=msr_.size(); ++i) {
			Size pdbres = pose.pdb_info()->number( msr_[i]->get_seqpos() );
			msr_[i]->set_seqpos(pdbres);
		}

	}
	
	
	
	virtual
	void
	graft_match(core::pose::Pose & pose) {
		pose::Pose scaffold;
		core::import_pose::pose_from_pdb( scaffold, basic::options::option[scaffold_protein].value() );
		
		TR << "Scaffold: " << scaffold.pdb_info()->name() << std::endl;
		
		utility::vector1< Size > match_residues;
		for( Size i(1); i <= pose.total_residue() - 1; ++i ) {
			match_residues.push_back( pose.pdb_info()->number(i) );
		}
		
		protocols::metal_interface::MatchGrafterOP match_grafter = new protocols::metal_interface::MatchGrafter;
		pose::Pose grafted_pose = match_grafter->graft( pose, scaffold );
		
		pose = grafted_pose;
		
	}
	
	
	virtual
	void
	repack(core::pose::Pose & pose) {
		
		using namespace core::pack::task;
		using namespace basic::options;
		
		TR << "initialize task factory" << std::endl;
		basic::options::option[basic::options::OptionKeys::packing::use_input_sc](true);
		TaskFactoryOP task_factory = new TaskFactory();
		task_factory->push_back(new operation::InitializeFromCommandline()); //ex1, ex1, minimize sidechains, use_input_sc
		
		core::pack::task::operation::RestrictToRepackingOP restrict_to_repack( new core::pack::task::operation::RestrictToRepacking() );
		task_factory->push_back( restrict_to_repack );
		
		operation::PreventRepackingOP prevent_repack = new operation::PreventRepacking();
		for(Size i(2); i <= msr_.size(); ++i) {
			prevent_repack->include_residue(msr_[i]->get_seqpos());
		}
		task_factory->push_back(prevent_repack);
		
		using namespace core::scoring;
		ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH );
		
		protocols::simple_moves::PackRotamersMoverOP packrot_mover = new protocols::simple_moves::PackRotamersMover;
		packrot_mover->score_function( scorefxn );
		packrot_mover->task_factory( task_factory );
		
		packrot_mover->apply( pose );
	}
	
	
	virtual
	void
	constrain_zinc(core::pose::Pose & pose) {
		protocols::metal_interface::AddZincSiteConstraintsOP constrain_zinc = new protocols::metal_interface::AddZincSiteConstraints( msr_ );
		constrain_zinc->add_constraints(pose);
	}
	

	virtual
	void
	find_neighbors(core::pose::Pose & pose) {
		using namespace pose::metrics;
		using namespace pose::metrics::simple_calculators;
		
		primary_shell_neighbors_.clear();
		
		
		for(Size i(2); i <= msr_.size(); ++i) {
			
			std::string const calc_stem("nbr_calc_");
			std::stringstream calcname;
			calcname << calc_stem << msr_[i]->get_seqpos();
			
			TR << "Registering NeighborsByDistance Calculator " << calcname.str() << std::endl;
			if( !CalculatorFactory::Instance().check_calculator_exists( calcname.str() ) ){
				CalculatorFactory::Instance().register_calculator( calcname.str(), new protocols::toolbox::pose_metric_calculators::NeighborsByDistanceCalculator( msr_[i]->get_seqpos() ) );
			}
			
			basic::MetricValue< SetSize > this_neighbors;
			pose.metric( calcname.str(), "neighbors", this_neighbors );
			
			
			TR << "These are the neighbor residues: " << std::endl;
			utility::vector1<Size> neighbors_list;
			for( std::set<Size>::const_iterator it(this_neighbors.value().begin()), end(this_neighbors.value().end()); it != end; ++it){
				TR << *it << "," << std::endl;
				bool is_zinc_coordinating_residue(false);
				for(Size j(1); j<=msr_.size(); ++j) {  //neighbor cannot be a zinc coordinating residue
					if(*it == msr_[j]->get_seqpos()) {
						is_zinc_coordinating_residue = true;
					}
				}
				if(!is_zinc_coordinating_residue) {
					neighbors_list.push_back(*it);
				}
			}
			TR << std::endl;
			primary_shell_neighbors_.push_back( neighbors_list ); //there will be three primary shell residues
			
		}
		return;
	}
	
	
	virtual
	void
	build_rotamers(core::pose::Pose & pose) {
		
		core::scoring::ScoreFunction dummy_sfxn;
		dummy_sfxn( pose );
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
		
		
		utility::vector1< utility::vector1<bool> > allow_single_aas;
		for(Size i(1); i<=20; ++i) {
			utility::vector1<bool> allow_single_aa(20, false);
			allow_single_aa[i] = true;
			allow_single_aas.push_back(allow_single_aa);
		}
		
		utility::vector1<Size> polar_aas;
		polar_aas.push_back(3);  //D
		polar_aas.push_back(4);  //E
		polar_aas.push_back(12); //N
		polar_aas.push_back(14); //Q
		polar_aas.push_back(16); //S
		polar_aas.push_back(17); //T

		
		core::pose::Pose const starting_pose( pose );
		starting_pose_ = starting_pose;
		
		for(Size i(1); i <= primary_shell_neighbors_.size(); ++i) {
			
			TR << "Shell " << i << std::endl;


			for(Size j(1); j <= primary_shell_neighbors_[i].size()-1; ++j) {
				Size this_res = primary_shell_neighbors_[i][j];
				
				if(this_res - 2 < 1 || this_res + 2 > pose.total_residue()) {
					continue;
				}
				
				TR << "we are on residue " << i << " " << this_res << std::endl;
				if(!pose.residue(this_res).is_protein()) {
					continue;
				}
				
				for(Size polar_aa_index(1); polar_aa_index <= polar_aas.size(); ++polar_aa_index) {
					
					Size aa_index = polar_aas[polar_aa_index];
					core::pack::task::TaskFactoryOP task_factory = new core::pack::task::TaskFactory();
					core::pack::task::operation::PreventRepackingOP prevent_repack = new core::pack::task::operation::PreventRepacking();
					core::pack::task::operation::RestrictAbsentCanonicalAASOP restrict_absent = new core::pack::task::operation::RestrictAbsentCanonicalAAS();
					restrict_absent->include_residue(0);
					//restrict_absent->keep_aas("S");
					restrict_absent->keep_aas(allow_single_aas[aa_index]);
					
					for(Size res(1); res <= pose.total_residue(); ++res) {
						if(res != this_res) {
							prevent_repack->include_residue( res );
						}
					
					}
					task_factory->push_back( restrict_absent );
					task_factory->push_back( prevent_repack );

					
					core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
					task->nonconst_residue_task(this_res).restrict_absent_canonical_aas( allow_single_aas[aa_index] );
					
					protocols::simple_moves::PackRotamersMoverOP pack_rotamers = new protocols::simple_moves::PackRotamersMover();
					pack_rotamers->score_function(scorefxn);
					pack_rotamers->task_factory(task_factory);
					
					pack_rotamers->apply( pose );
									
					core::pack::rotamer_set::RotamerSetFactory rsf;
					core::pack::rotamer_set::RotamerSetOP rotset( rsf.create_rotamer_set( pose.residue( this_res ) ) );
					rotset->set_resid( this_res );
					core::graph::GraphOP dummy_png = core::pack::create_packer_graph( pose, dummy_sfxn, task );
					rotset->build_rotamers( pose, dummy_sfxn, *task, dummy_png );
									
					for(Size krot(1); krot <= rotset->num_rotamers(); ++krot) {
						//TR << "ik_lys_ctp_asp" << rsd1 << " " << rsd2 << " " << krot << std::endl;
						pose.set_chi(1,this_res,rotset->rotamer(krot)->chi(1));
						pose.set_chi(2,this_res,rotset->rotamer(krot)->chi(2));
						//pose.set_chi(3,this_res,rotset->rotamer(krot)->chi(3));
						//pose.set_chi(4,this_res,rotset->rotamer(krot)->chi(4));

						TR << "Rotamer: " << pose.residue(this_res).name3() << " resi " << this_res << " rotamer " << krot << " " << pose.residue(this_res).chi(1) << " " << pose.residue(this_res).chi(2) << std::endl;
						
						//std::stringstream output_name;
						//output_name << "out_" << this_res << "_" << aa_index << "_" << krot << "_" << pose.pdb_info()->name();
						//pose.dump_pdb(output_name.str());
						
						TR << "Calling hbond_constraint function with " << this_res << " " << aa_index << " " << i << " " << krot << std::endl;
						hbond_constraint(pose, this_res, aa_index, i, krot);
					}
					
					pose = starting_pose_;
				}
			}
		}
	}
	
	
	virtual
	void
	hbond_constraint(core::pose::Pose & pose, Size hbond_res, Size aa_index, Size primary_shell_index, Size krot) {
		
		//Asp(3)  is OD1 or OD2
		//Glu(4)  is OE1 or OE1
		//Asn(12) is OD1
		//Gln(14) is OE1
		//Ser(16) is OG
		//Thr(17) is OG1
		
		core::pose::Pose const pose_nohbondcst_nomin( pose );
		
		//get atom index of backside His
		id::AtomID backside_id;
		id::AtomID hydrogen_id;

		std::string ligand_atom_name = msr_[primary_shell_index+1]->get_ligand_atom_name();
		Size primary_shell_resnum = msr_[primary_shell_index+1]->get_seqpos();
		
		if(ligand_atom_name == " NE2") { //histidine
			backside_id = core::id::AtomID( pose.residue(primary_shell_resnum).atom_index(" ND1"), primary_shell_resnum );
			hydrogen_id = core::id::AtomID( pose.residue(primary_shell_resnum).atom_index(" HD1"), primary_shell_resnum );
		}
		else if(ligand_atom_name == " ND1") { //histidine
			backside_id = core::id::AtomID( pose.residue(primary_shell_resnum).atom_index(" NE2"), primary_shell_resnum );
			hydrogen_id = core::id::AtomID( pose.residue(primary_shell_resnum).atom_index(" HE2"), primary_shell_resnum );
		}
		
		
		//get atom index of polar atom for hbond
		id::AtomID hbond_id;
		if(aa_index == 3) {
			hbond_id = id::AtomID( pose.residue(hbond_res).atom_index(" OD1"), hbond_res );
		}
		if(aa_index == 4) {
			hbond_id = id::AtomID( pose.residue(hbond_res).atom_index(" OE1"), hbond_res );
		}
		if(aa_index == 12) {
			hbond_id = id::AtomID( pose.residue(hbond_res).atom_index(" OD1"), hbond_res );
		}
		if(aa_index == 14) {
			hbond_id = id::AtomID( pose.residue(hbond_res).atom_index(" OE1"), hbond_res );
		}
		if(aa_index == 16) {
			hbond_id = id::AtomID( pose.residue(hbond_res).atom_index("OG"), hbond_res );
		}
		if(aa_index == 17) {
			hbond_id = id::AtomID( pose.residue(hbond_res).atom_index("OG1"), hbond_res );
		}
		
		
		point backside_atom_xyz = pose.residue(backside_id.rsd()).atom(backside_id.atomno()).xyz();
		point hbond_atom_xyz = pose.residue(hbond_id.rsd()).atom(hbond_id.atomno()).xyz();
		
		Real dist = backside_atom_xyz.distance(hbond_atom_xyz);
		
		TR << "Distance " << backside_id << " " << hbond_id << " " << dist << std::endl;
		
		//Constrain distance to 2.6, constrain D-H::O angle to 180
		using namespace scoring;
		using namespace scoring::constraints;
		using namespace conformation;
		using namespace id;
		using namespace basic::options;
		using numeric::conversions::radians;
		using numeric::conversions::degrees;
		
		FuncOP const hbond_dist_func( new HarmonicFunc( 2.6, 0.2 ) );
		FuncOP const hbond_angle_func( new HarmonicFunc( 180, 20 ) );
		
		core::scoring::constraints::AtomPairConstraintCOP hbond_dist_cst = new AtomPairConstraint( backside_id, hbond_id, hbond_dist_func );
		core::scoring::constraints::AngleConstraintCOP hbond_angle_cst = new AngleConstraint( backside_id, hydrogen_id, hbond_id, hbond_angle_func );
		
		TR << "Adding hbond constraints" << std::endl;
		pose.add_constraint( hbond_dist_cst );
		pose.add_constraint( hbond_angle_cst );
		
		
		scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
		scorefxn->set_weight( atom_pair_constraint, 1.0 );
		scorefxn->set_weight( angle_constraint, 1.0 );
		scorefxn->set_weight( dihedral_constraint, 1.0 );
		
		//kinematics::MoveMap mm; //I should really just allow residues i-4 to i+4 to minimize
		//mm.set_bb  ( true );
		//mm.set_chi ( true );
		//mm.set_jump( true );
		//mm.set( core::id::THETA, basic::options::option[ OptionKeys::relax::minimize_bond_angles ]() ); //will this be true for all residues in the protein?
		//mm.set( core::id::D, basic::options::option[ OptionKeys::relax::minimize_bond_lengths ]() );
		//core::optimization::MinimizerOptions cmin_options( "dfpmin", 0.00001, true, /*debug_derivs*/false, /*debug_verbose*/false );
		//core::optimization::CartesianMinimizer cartesian_minimizer;
		//cartesian_minimizer.run( pose, mm, *scorefxn, cmin_options );
		
		kinematics::MoveMapOP mm = new kinematics::MoveMap(); //I should really just allow residues i-4 to i+4 to minimize
		mm->set_bb  ( false );
		mm->set_chi ( hbond_res, true ); //will only minimize one sidechain
		mm->set_chi ( msr_[2], true );
		mm->set_chi ( msr_[3], true );
		mm->set_chi ( msr_[4], true );
		mm->set_jump(true); //to move the zinc atom
		
		
		//mm->set_jump( true );
		protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm, scorefxn, "dfpmin_armijo", 0.01, true );
		//TR << "Minimizing" << std::endl;
		
		
		//backrub
		//#include <basic/options/keys/backrub.OptionKeys.gen.hh> // we might not need this after all...
		protocols::backrub::BackrubMoverOP backrubmover = new protocols::backrub::BackrubMover();
		const core::Size backrub_iterations = 10;
		backrubmover->clear_segments();
		
		utility::vector1<Size> pivot_residues;
		pivot_residues.push_back(hbond_res - 2);
		pivot_residues.push_back(hbond_res - 1);
		pivot_residues.push_back(hbond_res);
		pivot_residues.push_back(hbond_res + 1);
		pivot_residues.push_back(hbond_res + 2);
		backrubmover->set_pivot_residues(pivot_residues);
		
		utility::vector1<std::string> pivot_atoms;
		pivot_atoms.push_back("N");
		pivot_atoms.push_back("CA");
		pivot_atoms.push_back("C");
		
		backrubmover->set_pivot_atoms(pivot_atoms);
		//backrubmover->set_input_pose(pose_nohbondcst_nomin);
		
		for(Size bi(1); bi <= 20; ++bi) {
			
			TR << "Backrubbing " << bi << std::endl;
			backrubmover->apply(pose);
			TR << "Minimizing" << std::endl;
			min_mover->apply(pose);
			
			point new_backside_atom_xyz = pose.residue(backside_id.rsd()).atom(backside_id.atomno()).xyz();
			point new_hbond_atom_xyz = pose.residue(hbond_id.rsd()).atom(hbond_id.atomno()).xyz();
			point new_hydrogen_xyz = pose.residue(hydrogen_id.rsd()).atom(hydrogen_id.atomno()).xyz();
			
			//Real hbond_dist_score = hbond_dist_cst->score(new_backside_atom_xyz, new_hbond_atom_xyz);
			//Real hbond_angle_score = hbond_angle_cst->score(new_backside_atom_xyz, new_hydrogen_xyz, new_hbond_atom_xyz);
			
			Real current_hbond_dist = new_backside_atom_xyz.distance(new_hbond_atom_xyz);
			Real current_hbond_angle = degrees(angle_of( new_backside_atom_xyz, new_hydrogen_xyz, new_hbond_atom_xyz ));
			
			
			TR << "Evaluating new hbond.  Dist: " << current_hbond_dist << ".  Angle: " << current_hbond_angle << std::endl;
			
			if(current_hbond_dist < 3.0 && current_hbond_angle > 120) {
				std::stringstream output_name;
				output_name << "out_" << hbond_res << "_" << aa_index << "_" << krot << "_" << bi << "_" << pose.pdb_info()->name();
				pose.dump_pdb(output_name.str());
			}
			
			
			//movemap for minimization (what is the scorefunction?  hbond only? hbond constraints and zinc constraints only?)
			
			
			
			pose = pose_nohbondcst_nomin; //will keep zinc constraints but lose hbond constraints, will lose cartesian minimization
		}

		
	}
	
	
	virtual
	std::string
	get_name() const { return "HbondZinc"; }



private:

	core::pose::Pose starting_pose_;
	utility::vector1< protocols::metal_interface::MetalSiteResidueOP > msr_;
	utility::vector1< utility::vector1<Size> > primary_shell_neighbors_;
	
	// basic::MetricValue< id::AtomID_Map< Real > > atom_sasa_;
	// basic::MetricValue< id::AtomID_Map< Size > > atom_hbonds_;

};

typedef utility::pointer::owning_ptr< HbondZinc > HbondZincOP;

int main( int argc, char* argv[] )
{
	try {
	using basic::options::option;
	option.add( scaffold_protein, "scaffold for grafting" ).def("1YZM.pdb");

  	devel::init(argc, argv);
  	protocols::jd2::JobDistributor::get_instance()->go(new HbondZinc);
		


  	TR << "************************d**o**n**e**************************************" << std::endl;

	} catch (utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

  return 0;
}

