// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/bder/HbondZinc.cc
/// @brief
/// @details
/// @author Bryan Der

#include <devel/init.hh>
#include <protocols/metal_interface/MatchGrafter.hh>
#include <protocols/metal_interface/ZincSiteFinder.hh>
#include <protocols/metal_interface/MetalSiteResidue.hh>

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
#include <core/types.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
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
#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>


//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static THREAD_LOCAL basic::Tracer TR("apps.pilot.bder.HbondZinc");

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
		find_neighbors(pose);
		build_rotamers(pose);

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

	}


	virtual
	void
	repack(core::pose::Pose & pose) {

		using namespace core::pack::task;
		using namespace basic::options;

		TR << "initialize task factory" << std::endl;
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

		protocols::simple_moves::PackRotamersMoverOP packrot_mover = new protocols::simple_moves::PackRotamersMover;
		packrot_mover->score_function( scorefxn );
		packrot_mover->task_factory( task_factory );

		packrot_mover->apply( pose );
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

			TR << "Registering NeighborsByDistance Calculator " << calcname.str() << std::endl;
			if ( !CalculatorFactory::Instance().check_calculator_exists( calcname.str() ) ) {
				CalculatorFactory::Instance().register_calculator( calcname.str(), new protocols::toolbox::pose_metric_calculators::NeighborsByDistanceCalculator( msr_[i]->get_seqpos() ) );
			}

			basic::MetricValue< SetSize > this_neighbors;
			pose.metric( calcname.str(), "neighbors", this_neighbors );


			TR << "These are the neighbor residues: " << std::endl;
			utility::vector1<Size> neighbors_list;
			for ( std::set<Size>::const_iterator it(this_neighbors.value().begin()), end(this_neighbors.value().end()); it != end; ++it ) {
				TR << *it << "," << std::endl;
				neighbors_list.push_back(*it);
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


		utility::vector1<bool> allow_Ser(20, false);
		allow_Ser[16] = true;

		core::pose::Pose const starting_pose( pose );

		for ( Size i(1); i <= primary_shell_neighbors_.size(); ++i ) {

			TR << "Shell " << i << std::endl;


			for ( Size j(1); j <= primary_shell_neighbors_[i].size()-1; ++j ) {
				Size this_res = primary_shell_neighbors_[i][j];

				TR << "we are on residue " << i << " " << this_res << std::endl;
				if ( !pose.residue(this_res).is_protein() ) {
					continue;
				}

				core::pack::task::TaskFactoryOP task_factory = new core::pack::task::TaskFactory();
				core::pack::task::operation::PreventRepackingOP prevent_repack = new core::pack::task::operation::PreventRepacking();
				core::pack::task::operation::RestrictAbsentCanonicalAASOP restrict_absent = new core::pack::task::operation::RestrictAbsentCanonicalAAS();
				restrict_absent->include_residue(0);
				restrict_absent->keep_aas("S");

				for ( Size res(1); res <= pose.size(); ++res ) {
					if ( res != this_res ) {
						prevent_repack->include_residue( res );
					}

				}
				task_factory->push_back( restrict_absent );
				task_factory->push_back( prevent_repack );


				core::pack::task::PackerTaskOP Ser_task = core::pack::task::TaskFactory::create_packer_task( pose );
				Ser_task->nonconst_residue_task(this_res).restrict_absent_canonical_aas( allow_Ser );

				protocols::simple_moves::PackRotamersMoverOP pack_rotamers = new protocols::simple_moves::PackRotamersMover();
				pack_rotamers->score_function(scorefxn);
				pack_rotamers->task_factory(task_factory);

				pack_rotamers->apply( pose );

				core::pack::rotamer_set::RotamerSetFactory rsf;
				core::pack::rotamer_set::RotamerSetOP rotset( rsf.create_rotamer_set( pose.residue( this_res ) ) );
				rotset->set_resid( this_res );
				utility::graph::GraphOP dummy_png = core::pack::create_packer_graph( pose, dummy_sfxn, Ser_task );
				rotset->build_rotamers( pose, dummy_sfxn, *Ser_task, dummy_png );

				for ( Size krot(1); krot <= rotset->num_rotamers(); ++krot ) {
					//TR << "ik_lys_ctp_asp" << rsd1 << " " << rsd2 << " " << krot << std::endl;
					pose.set_chi(1,this_res,rotset->rotamer(krot)->chi(1));
					pose.set_chi(2,this_res,rotset->rotamer(krot)->chi(2));
					//pose.set_chi(3,this_res,rotset->rotamer(krot)->chi(3));
					//pose.set_chi(4,this_res,rotset->rotamer(krot)->chi(4));

					TR << "Serine: " << pose.residue(this_res).name3() << " resi " << this_res << " rotamer " << krot << " " << pose.residue(this_res).chi(1) << " " << pose.residue(this_res).chi(2) << std::endl;

					std::stringstream output_name;
					output_name << "out_" << this_res << "_" << krot << "_" << pose.pdb_info()->name();
					pose.dump_pdb(output_name.str());
				}

				pose = starting_pose;
			}
		}



	}


	virtual
	std::string
	get_name() const { return "HbondZinc"; }



private:

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

