// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/bder/HbondZincBackbone.cc
/// @brief
/// @details
/// @author Bryan Der

#include <devel/init.hh>
#include <protocols/metal_interface/MatchGrafter.hh>
#include <protocols/metal_interface/ZincSiteFinder.hh>
#include <protocols/metal_interface/MetalSiteResidue.hh>
#include <protocols/metal_interface/AddZincSiteConstraints.hh>

#include <utility/graph/Graph.hh>
#include <core/pack/packer_neighbors.hh>


#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

#include <core/types.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/backrub/BackrubMover.hh>

#include <protocols/simple_pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/simple_pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/pose_metric_calculators/NeighborsByDistanceCalculator.hh>

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
static basic::Tracer TR("apps.pilot.bder.HbondZincBackbone");

typedef std::set< core::Size > SetSize;

using namespace core;

basic::options::FileOptionKey const scaffold_protein( "scaffold_protein" );


/// @brief
class HbondZincBackbone : public protocols::moves::Mover {
public:
	HbondZincBackbone()
	{
	}
	virtual ~HbondZincBackbone(){};


	virtual
	void
	apply( core::pose::Pose & pose ){

		TR << "HbondZincBackbone protocol" << std::endl;

		find_zinc(pose); //pose is match
		match_name_ = pose.pdb_info()->name();

		graft_match(pose); //pose is now grafted scaffold
		scaffold_name_ = pose.pdb_info()->name();

		repack(pose);
		constrain_zinc(pose);
		find_neighbors(pose);
		carbonyl_hbond_search(pose);

		return;
	}


	virtual
	void
	find_zinc(core::pose::Pose & pose) {

		protocols::metal_interface::ZincSiteFinderOP zinc_finder = new protocols::metal_interface::ZincSiteFinder();
		msr_ = zinc_finder->find_zinc_site(pose);

		for ( Size i(2); i<=msr_.size(); ++i ) {
			Size pdbres = pose.pdb_info()->number( msr_[i]->get_seqpos() );
			msr_[i]->set_seqpos(pdbres);

			TR << "MSR " << i << " " << pdbres << std::endl;
		}

	}



	virtual
	void
	graft_match(core::pose::Pose & pose) {
		pose::Pose scaffold;
		core::import_pose::pose_from_file( scaffold, basic::options::option[scaffold_protein].value() , core::import_pose::PDB_file);

		TR << "Scaffold: " << scaffold.pdb_info()->name() << std::endl;

		utility::vector1< Size > match_residues;
		for ( Size i(1); i <= pose.size() - 1; ++i ) {
			match_residues.push_back( pose.pdb_info()->number(i) );
		}

		protocols::metal_interface::MatchGrafterOP match_grafter = new protocols::metal_interface::MatchGrafter;
		pose::Pose grafted_pose = match_grafter->graft( pose, scaffold );

		pose = grafted_pose;

		Size znres = pose.size();
		msr_[1]->set_seqpos(znres);
		TR << "MSR " << 1 << " " << znres << std::endl;

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
		for ( Size i(2); i <= msr_.size(); ++i ) {
			prevent_repack->include_residue(msr_[i]->get_seqpos());
		}
		task_factory->push_back(prevent_repack);

		using namespace core::scoring;
		ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH );

		protocols::minimization_packing::PackRotamersMoverOP packrot_mover = new protocols::minimization_packing::PackRotamersMover;
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


		for ( Size i(2); i <= msr_.size(); ++i ) {

			std::string const calc_stem("nbr_calc_");
			std::stringstream calcname;
			calcname << calc_stem << msr_[i]->get_seqpos();

			//TR << "Registering NeighborsByDistance Calculator " << calcname.str() << std::endl;
			if ( !CalculatorFactory::Instance().check_calculator_exists( calcname.str() ) ) {
				CalculatorFactory::Instance().register_calculator( calcname.str(), new protocols::pose_metric_calculators::NeighborsByDistanceCalculator( msr_[i]->get_seqpos() ) );
			}

			basic::MetricValue< SetSize > this_neighbors;
			pose.metric( calcname.str(), "neighbors", this_neighbors );

			//TR << "These are the neighbor residues: " << std::endl;
			utility::vector1<Size> neighbors_list;
			for ( std::set<Size>::const_iterator it(this_neighbors.value().begin()), end(this_neighbors.value().end()); it != end; ++it ) {
				//TR << *it << "," << std::endl;
				bool is_zinc_coordinating_residue(false);
				for ( Size j(1); j<=msr_.size(); ++j ) {  //neighbor cannot be a zinc coordinating residue
					if ( *it == msr_[j]->get_seqpos() ) {
						is_zinc_coordinating_residue = true;
					}
				}
				if ( !is_zinc_coordinating_residue ) {
					neighbors_list.push_back(*it);
				}
			}
			//TR << std::endl;
			primary_shell_neighbors_.push_back( neighbors_list ); //there will be three primary shell residues

		}
		return;
	}


	virtual
	void
	carbonyl_hbond_search(core::pose::Pose & pose) {

		core::pose::Pose const pose_nohbondcst_nomin( pose );

		for ( Size his_res(1); his_res <= 3; ++his_res ) {

			//get atom index of backside His
			id::AtomID backside_id;
			id::AtomID hydrogen_id;

			std::string ligand_atom_name = msr_[his_res+1]->get_ligand_atom_name();
			Size primary_shell_resnum = msr_[his_res+1]->get_seqpos();

			if ( ligand_atom_name == " NE2" ) { //histidine
				backside_id = core::id::AtomID( pose.residue(primary_shell_resnum).atom_index(" ND1"), primary_shell_resnum );
				hydrogen_id = core::id::AtomID( pose.residue(primary_shell_resnum).atom_index(" HD1"), primary_shell_resnum );
			} else if ( ligand_atom_name == " ND1" ) { //histidine
				backside_id = core::id::AtomID( pose.residue(primary_shell_resnum).atom_index(" NE2"), primary_shell_resnum );
				hydrogen_id = core::id::AtomID( pose.residue(primary_shell_resnum).atom_index(" HE2"), primary_shell_resnum );
			}


			for ( Size neighbor_res_index(1); neighbor_res_index <= primary_shell_neighbors_[his_res].size(); ++neighbor_res_index ) {
				Size neighbor_res = primary_shell_neighbors_[his_res][neighbor_res_index];
				//TR << "neighbor res " << neighbor_res << std::endl;


				//get atom index of polar atom for hbond
				id::AtomID carbonyl_id( pose.residue(neighbor_res).atom_index("O"), neighbor_res );

				point backside_xyz = pose.residue(backside_id.rsd()).atom(backside_id.atomno()).xyz();
				point hydrogen_xyz = pose.residue(hydrogen_id.rsd()).atom(hydrogen_id.atomno()).xyz();
				point carbonyl_xyz = pose.residue(carbonyl_id.rsd()).atom(carbonyl_id.atomno()).xyz();

				Real distance = backside_xyz.distance(carbonyl_xyz);
				Real distance_hyd = hydrogen_xyz.distance(carbonyl_xyz);

				if ( distance < 4.0 && distance_hyd < 3.5 ) {

					TR << scaffold_name_ << " " << match_name_ << " His_index " << his_res << " His_resnum " << msr_[his_res+1]->get_seqpos() << " neighbor_res " << neighbor_res << " distance " << distance << " Hbond_energy " << has_backbone_hbond(pose, carbonyl_id) << std::endl;

				}



				//TR << pose.pdb_info()->name() << " His_index " << his_res << " His_resnum " << msr_[his_res+1]->get_seqpos() << " neighbor_res " << neighbor_res << " distance " << distance << std::endl;

			}
		}
	}


	virtual
	Real
	has_backbone_hbond(core::pose::Pose & pose, id::AtomID carbonyl_id) {

		core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::STANDARD_WTS, core::scoring::SCORE12_PATCH );
		scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn->energy_method_options() );
		energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
		scorefxn->set_energy_method_options( energymethodoptions );
		scorefxn->score(pose);

		Size res = carbonyl_id.rsd();

		Real bb_hbond_energy = pose.energies().residue_total_energies(res)[scoring::hbond_sr_bb] + pose.energies().residue_total_energies(res)[scoring::hbond_lr_bb]; //compiler warning

		return bb_hbond_energy;

	}



	virtual
	std::string
	get_name() const { return "HbondZincBackbone"; }



private:

	core::pose::Pose starting_pose_;
	utility::vector1< protocols::metal_interface::MetalSiteResidueOP > msr_;
	utility::vector1< utility::vector1<Size> > primary_shell_neighbors_;


	std::string match_name_;
	std::string scaffold_name_;
	// basic::MetricValue< id::AtomID_Map< Real > > atom_sasa_;
	// basic::MetricValue< id::AtomID_Map< Size > > atom_hbonds_;

};

typedef utility::pointer::owning_ptr< HbondZincBackbone > HbondZincBackboneOP;

int main( int argc, char* argv[] )
{
	try {
		using basic::options::option;
		option.add( scaffold_protein, "scaffold for grafting" ).def("1YZM.pdb");

		devel::init(argc, argv);
		protocols::jd2::JobDistributor::get_instance()->go(new HbondZincBackbone);



		TR << "************************d**o**n**e**************************************" << std::endl;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}

