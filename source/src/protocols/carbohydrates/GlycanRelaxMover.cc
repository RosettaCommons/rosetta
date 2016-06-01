// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/carbohydrates/GlycanRelaxMover.cc
/// @brief Main mover for Glycan Relax, which optimizes glycans in a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) and Jason W. Labonte (JWLabonte@jhu.edu)

#include <protocols/carbohydrates/GlycanRelaxMover.hh>
#include <protocols/carbohydrates/GlycanRelaxMoverCreator.hh>
#include <protocols/carbohydrates/LinkageConformerMover.hh>
#include <protocols/simple_moves/BBDihedralSamplerMover.hh>




#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/selection.hh>
#include <core/pose/util.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.hh>
#include <core/id/types.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/bb_sampler/SugarBBSampler.hh>
#include <protocols/simple_moves/bb_sampler/SmallBBSampler.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <numeric/random/random.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/carbohydrates.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.carbohydrates.GlycanRelaxMover" );


namespace protocols {
namespace carbohydrates {
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::bb_sampler;
using namespace core::pack::task;
using namespace basic::options;
using namespace core::select::residue_selector;

GlycanRelaxMover::GlycanRelaxMover():
	protocols::moves::Mover( "GlycanRelaxMover" ),
	full_movemap_(/*NULL*/),
	glycan_movemap_(/*NULL*/),
	tf_(/*NULL*/),
	mc_(/*NULL*/),
	scorefxn_(/* NULL */),
	linkage_mover_(/* NULL */)
{
	set_defaults();
}

GlycanRelaxMover::GlycanRelaxMover(
	core::kinematics::MoveMapCOP mm,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::Size rounds):
	protocols::moves::Mover("GlycanRelaxMover"),
	glycan_movemap_(/*NULL*/),
	tf_(/*NULL*/),
	mc_(/* NULL */),
	scorefxn_(scorefxn),
	linkage_mover_(/*NULL*/)
{
	full_movemap_ = mm->clone();
	set_defaults();
	set_rounds(rounds);
}

GlycanRelaxMover::~GlycanRelaxMover(){}

GlycanRelaxMover::GlycanRelaxMover( GlycanRelaxMover const & src ):
	protocols::moves::Mover( src ),
	full_movemap_(src.full_movemap_),
	glycan_movemap_(src.glycan_movemap_),
	tf_(src.tf_),
	mc_(src.mc_),
	scorefxn_(src.scorefxn_),
	linkage_mover_(src.linkage_mover_),
	weighted_random_mover_(src.weighted_random_mover_),
	min_mover_(src.min_mover_),
	packer_(src.packer_),
	rounds_(src.rounds_),
	kt_(src.kt_),
	accept_log_(src.accept_log_),
	test_(src.test_),
	final_min_(src.final_min_),
	pack_glycans_(src.pack_glycans_),
	random_start_(src.random_start_),
	sugar_bb_start_(src.sugar_bb_start_),
	total_glycan_residues_(src.total_glycan_residues_),
	pymol_movie_(src.pymol_movie_),
	ref_pose_name_(src.ref_pose_name_),
	use_branches_( src.use_branches_),
	parsed_positions_(src.parsed_positions_)

{
}

void
GlycanRelaxMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & pose)
{
	kt_ = tag->getOption< core::Real >("kt", kt_);
	rounds_ = tag->getOption< core::Size >("rounds", rounds_);

	pack_glycans_ = tag->getOption< bool >("pack_glycans", pack_glycans_);
	final_min_ = tag->getOption< bool >("final_min", final_min_);

	pymol_movie_ = tag->getOption< bool >("pymol_movie", pymol_movie_);


	//Movemap
	if ( protocols::rosetta_scripts::has_branch(tag, "MoveMap") ) {
		full_movemap_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap() );

		//protocols::rosetta_scripts::add_movemaps_to_datamap(tag, pose, data, false);
		protocols::rosetta_scripts::parse_movemap( tag, pose, full_movemap_, datamap, false );
	}

	//Scorefunction.
	if ( tag->hasOption("scorefxn") ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, datamap );
	}

	ref_pose_name_ = tag->getOption< std::string >( "ref_pose_name", ref_pose_name_ );



	if ( tag->hasOption("branch") ) {
		parsed_positions_.push_back( tag->getOption< std::string >("branch"));
		use_branches_ = true;
	} else if ( tag->hasOption("branches") ) {
		parsed_positions_ = utility::string_split_multi_delim( tag->getOption< std::string >("branches"), ",'`~+*&|;. ");
		use_branches_ = true;
	}

	random_start_ = tag->getOption< bool >( "random_start", random_start_);
	sugar_bb_start_ = tag->getOption< bool >("sugar_bb_start", sugar_bb_start_);
	
	if (tag->hasOption("task_operations") ){
		TaskFactoryOP tf( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
		set_taskfactory( tf );
	}
	
	pack_distance_ = tag->getOption< core::Real >("pack_distance", pack_distance_);

}

void
GlycanRelaxMover::set_defaults(){
	rounds_ = 75; //Means absolutely nothing right now - Actually set from cmd line.
	test_ = false;

	pack_glycans_ = false;
	final_min_ = true;

	set_cmd_line_defaults();

	total_glycan_residues_ = 0;
	ref_pose_name_ = "";
	pack_distance_ = 6.0;
}


void
GlycanRelaxMover::set_cmd_line_defaults(){

	rounds_ = option [OptionKeys::carbohydrates::glycan_relax::glycan_relax_rounds]();
	test_ = option [OptionKeys::carbohydrates::glycan_relax::glycan_relax_test]();
	pack_glycans_ = option[ OptionKeys::carbohydrates::glycan_relax::pack_glycans]();
	final_min_ = option[ OptionKeys::carbohydrates::glycan_relax::final_min_glycans]();
	pymol_movie_ = option[ OptionKeys::carbohydrates::glycan_relax::glycan_relax_movie]();
	kt_ = option[ OptionKeys::carbohydrates::glycan_relax::glycan_relax_kt]();
	random_start_ = option[ OptionKeys::carbohydrates::glycan_relax::glycan_relax_random_start]();
	sugar_bb_start_ = option[ OptionKeys::carbohydrates::glycan_relax::glycan_relax_sugar_bb_start]();


}

protocols::moves::MoverOP
GlycanRelaxMover::clone() const{
	return protocols::moves::MoverOP( new GlycanRelaxMover( *this ) );
}

/*
GlycanRelaxMover & GlycanRelaxMoveroperator=( GlycanRelaxMover const & src){
return GlycanRelaxMover( src );
}
*/


moves::MoverOP
GlycanRelaxMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new GlycanRelaxMover );
}

std::string
GlycanRelaxMover::get_name() const {
	return "GlycanRelaxMover";
}

void
GlycanRelaxMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

std::ostream &operator<< (std::ostream &os, GlycanRelaxMover const &mover)
{
	mover.show(os);
	return os;
}

void
GlycanRelaxMover::set_movemap(core::kinematics::MoveMapCOP movemap){
	full_movemap_ = movemap->clone();
}

void
GlycanRelaxMover::set_taskfactory(core::pack::task::TaskFactoryCOP tf){
	tf_ = tf;
}

void
GlycanRelaxMover::setup_default_task_factory(utility::vector1< bool > const & glycan_positions ){
	using namespace core::pack::task::operation;
	
	
	TaskFactoryOP tf = TaskFactoryOP( new TaskFactory());
	tf->push_back(InitializeFromCommandlineOP( new InitializeFromCommandline));
	
	//If a resfile is provided, we just use that and get out.
	if (option[ OptionKeys::packing::resfile ].user() ){
		tf->push_back( ReadResfileOP( new ReadResfile()) );
	}
	else {
	
		//Select all the NON-glycan and neighbors and then turn them off.
		NeighborhoodResidueSelectorOP neighbor_selector = NeighborhoodResidueSelectorOP( new NeighborhoodResidueSelector(glycan_positions, pack_distance_, true /* include focus */));
		PreventRepackingRLTOP prevent_repacking = PreventRepackingRLTOP( new PreventRepackingRLT());
		
		OperateOnResidueSubsetOP subset_op = OperateOnResidueSubsetOP( new OperateOnResidueSubset( prevent_repacking, neighbor_selector, true /* flip */));
		tf->push_back( subset_op );
		tf->push_back( RestrictToRepackingOP( new RestrictToRepacking()));
	}
	
	set_taskfactory(tf);
}




void
GlycanRelaxMover::set_scorefunction(core::scoring::ScoreFunctionCOP scorefxn){
	scorefxn_ = scorefxn;
}

void
GlycanRelaxMover::set_kt(core::Real kt){
	kt_ = kt;
}

void
GlycanRelaxMover::set_rounds(core::Size rounds){
	rounds_ = rounds;
}


void
GlycanRelaxMover::init_objects(core::pose::Pose & pose ){


	using namespace core::pose::carbohydrates;
	using namespace protocols::moves;
	using namespace protocols::simple_moves;
	using namespace core::kinematics;

	TR << "initializing objects " << std::endl;
	total_glycan_residues_ = 0;
	
	utility::vector1< bool > glycan_positions( pose.total_residue(), false);
	

	//Create Scorefunction if needed.
	if ( ! scorefxn_ ) {
		scorefxn_ = core::scoring::get_score_function();
	}

	//Create a MonteCarlo object
	scorefxn_->score(pose);
	mc_ = MonteCarloOP( new MonteCarlo(pose, *scorefxn_, kt_) );
	TR << "mc" <<std::endl;

	//Setup Movemaps
	glycan_movemap_ = MoveMapOP( new MoveMap );

	if ( ! full_movemap_ ) {
		MoveMapOP mm = MoveMapOP( new MoveMap);
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			mm->set_bb( i, true);
			mm->set_chi( i, true);
		}
		full_movemap_ = mm;
	}

	///////////  Setup Glycan Movemap ////////////////
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		//TR << i << std::endl;
		core::Size resnum = i;
		if ( ref_pose_name_ != "" ) {
			resnum = pose.corresponding_residue_in_current( i, ref_pose_name_ );
		}

		if ( full_movemap_->get_bb( resnum ) && pose.residue( resnum ).is_carbohydrate() ) {
			total_glycan_residues_+=1;
			glycan_movemap_->set_bb( resnum , true );
			glycan_movemap_->set_chi(resnum , true );
			glycan_positions[ resnum ] = true;

		} else {
			glycan_movemap_->set_bb( resnum , false );
			glycan_movemap_->set_chi(resnum , false ); //Not using Movemap to pack here - must deal with that another time!
		}
	}

	///////////  Expand Branches and Parsed Residues ////////////////
	/// Used for RosettaScripts.
	/// In code, you can do this manually using get_carbohydrate_residues_upstream in pose/carbo../util
	///   Use this for creating a Movemap.
	///
	///
	if ( parsed_positions_.size() > 0 ) {
		for ( core::Size i = 1; i <= parsed_positions_.size(); ++ i ) {
			core::Size resnum = core::pose::parse_resnum( parsed_positions_[ i ], pose);

			if ( ref_pose_name_ != "" ) {
				resnum = pose.corresponding_residue_in_current( resnum, ref_pose_name_ );
			}

			if ( use_branches_ ) {
				std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > > res_and_tips;

				res_and_tips = get_carbohydrate_residues_upstream( pose, resnum );
				utility::vector1< core::Size > branching_residues = res_and_tips.first;
				branching_residues.push_back( resnum );
				for ( core::Size x = 1; x <= branching_residues.size(); ++x ) {
					core::Size branching_resnum = branching_residues[ x ];

					if ( pose.residue( branching_resnum ).is_carbohydrate() ) {
						glycan_movemap_->set_bb( branching_resnum , true );
						glycan_movemap_->set_chi(branching_resnum , true );
						total_glycan_residues_+=1;
						glycan_positions[ branching_resnum ] = true;
					}
				}
			}
		}
	}

	if ( total_glycan_residues_ == 0 ) {
		utility_exit_with_message( " No glycan residues in pose.  Cannot continue. ");
	}

	//////////  Setup Phi/Psi/N-Omega movemaps. /////////////////////////
	MoveMapOP sugar_bb_movemap = MoveMapOP( new MoveMap() );
	MoveMapOP glycan_dih_movemap = MoveMapOP( new MoveMap() );

	SmallBBSamplerOP random_sampler = SmallBBSamplerOP( new SmallBBSampler( 360.0 ) ); // +/- 180 degrees
	SugarBBSamplerOP random_sugar_sampler = SugarBBSamplerOP( new SugarBBSampler( ) );


	core::Size max_glycan_dihedrals = 2;
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {

		//Turn off if not carbohydrates or if N terminal carbohydrate not attached to anything.
		if ( ! glycan_movemap_->get_bb( i )  || find_seqpos_of_saccharides_parent_residue( pose.residue( i ) ) == 0 ) {
			sugar_bb_movemap->set_bb( i, false );
			glycan_dih_movemap->set_bb( i, false );
			continue;
		}

		//Need Psi Movemap (Non-Exocyclic) and Omega for only Exocyclic
		core::Size n_dihedrals = get_n_glycosidic_torsions_in_res( pose, i );

		sugar_bb_movemap->set_bb(i , 1, true); //Turn on for Phi
		sugar_bb_movemap->set_bb(i, 2, true);

		//Turn off psi for residue if it has any omegas.  We don't have data for this.
		if ( n_dihedrals > 2 ) {
			sugar_bb_movemap->set_bb( i, 2, false );
		}

		if ( n_dihedrals > max_glycan_dihedrals ) {
			max_glycan_dihedrals = n_dihedrals;
		}

		//Turn on only dihedrals for which these residues actually have
		for ( core::Size torsion_id = 1; torsion_id <= n_dihedrals; ++torsion_id ) {

			glycan_dih_movemap->set_bb( i, torsion_id, true );

			//Randomize starting structure if set.
			if ( random_start_ ) {
				random_sampler->set_torsion_type( torsion_id);
				random_sampler->set_torsion_to_pose(pose, i);
				
			} else if ( sugar_bb_start_ ) {

				//Continue if we don't have sugar bb data.  Its ok.
				try {
					random_sugar_sampler->set_torsion_type( torsion_id );
					random_sugar_sampler->set_torsion_to_pose(pose, i);
					
				} catch ( utility::excn::EXCN_Base& excn ) {
					continue;
				}

			}
		}

	}

	TR << "Modeling " << total_glycan_residues_ << " glycan residues" << std::endl;
	//pose.dump_pdb("post_random.pdb");

	////////////////// Mover Setup //////////////
	//Create Movers that will be part of our sequence mover.
	SugarBBSamplerOP phi_sugar_sampler = SugarBBSamplerOP( new SugarBBSampler( core::id::phi_dihedral ) );
	SugarBBSamplerOP psi_sugar_sampler = SugarBBSamplerOP( new SugarBBSampler( core::id::psi_dihedral ) );


	BBDihedralSamplerMoverOP sugar_sampler_mover = BBDihedralSamplerMoverOP( new BBDihedralSamplerMover );
	BBDihedralSamplerMoverOP small_sampler_mover = BBDihedralSamplerMoverOP( new BBDihedralSamplerMover );

	//Tolerance of .001 can significantly decrease energies, but it takes ~ 25% longer.  This tolerance should be optimized here.
	min_mover_ = MinMoverOP( new MinMover( glycan_movemap_->clone(), scorefxn_, "dfpmin_armijo_nonmonotone", 0.01, false /* use_nblist*/ ) );
	linkage_mover_ = LinkageConformerMoverOP( new LinkageConformerMover( glycan_movemap_ ));

	//Setup Sugar Samplers
	sugar_sampler_mover->add_sampler( phi_sugar_sampler );
	sugar_sampler_mover->add_sampler( psi_sugar_sampler );
	sugar_sampler_mover->set_movemap( sugar_bb_movemap );

	////////// Sequence Mover Setup //////////////
	//Settings for linkage conformer mover here!
	weighted_random_mover_ = RandomMoverOP(new RandomMover);
	weighted_random_mover_->add_mover(sugar_sampler_mover, .40);
	weighted_random_mover_->add_mover(linkage_mover_, .20);
	weighted_random_mover_->add_mover(min_mover_, .05);


	//Setup Small Sampler
	BBDihedralSamplerMoverOP glycan_small_mover = BBDihedralSamplerMoverOP( new BBDihedralSamplerMover() );
	BBDihedralSamplerMoverOP glycan_medium_mover= BBDihedralSamplerMoverOP( new BBDihedralSamplerMover() );
	BBDihedralSamplerMoverOP glycan_large_mover = BBDihedralSamplerMoverOP( new BBDihedralSamplerMover() );


	glycan_small_mover->set_movemap(  glycan_dih_movemap );
	glycan_medium_mover->set_movemap( glycan_dih_movemap );
	glycan_large_mover->set_movemap(  glycan_dih_movemap );


	for ( core::Size i =1; i <= max_glycan_dihedrals; ++i ) {
		core::id::MainchainTorsionType dih_type = static_cast< core::id::MainchainTorsionType >( i );

		SmallBBSamplerOP small_sampler = SmallBBSamplerOP( new SmallBBSampler( dih_type, 30 ) ); // +/- 15 degrees
		SmallBBSamplerOP medium_sampler= SmallBBSamplerOP( new SmallBBSampler( dih_type, 90 ) ); // +/- 45 degrees
		SmallBBSamplerOP large_sampler = SmallBBSamplerOP( new SmallBBSampler( dih_type, 180) ); // +/- 90 degrees


		glycan_small_mover->add_sampler( small_sampler );
		glycan_medium_mover->add_sampler( medium_sampler );
		glycan_large_mover->add_sampler( large_sampler );
	}

	weighted_random_mover_->add_mover( glycan_small_mover, 0.17142857142857143 );
	weighted_random_mover_->add_mover( glycan_medium_mover, 0.08571428571428572 );
	weighted_random_mover_->add_mover( glycan_large_mover, 0.04285714285714286 );
	
	if (! tf_ ){
		setup_default_task_factory(glycan_positions);
	}
	
	//Setup pack rotamers mover, add it.
	packer_ = PackRotamersMoverOP( new PackRotamersMover() );
	packer_->score_function(scorefxn_);
	packer_->task_factory(tf_);
	
	weighted_random_mover_->add_mover( packer_, .05);

}

void
GlycanRelaxMover::apply( core::pose::Pose& pose ){
	using namespace core::kinematics;
	using namespace core::chemical::carbohydrates;
	using namespace core::scoring;
	using namespace protocols::moves;
	using utility::to_string;

	accept_log_.clear();

	init_objects( pose );

	PyMolMover pmm_accepts = PyMolMover();
	PyMolMover pmm_trials  = PyMolMover();
	pmm_accepts.keep_history( true );
	pmm_trials.keep_history( true );
	pmm_accepts.set_PyMol_model_name( "accepts_"+ pmm_accepts.get_PyMol_model_name( pose ));
	pmm_trials.set_PyMol_model_name(  "trials_" + pmm_trials.get_PyMol_model_name( pose ));

	TR << "Initialized" << std::endl;
	utility::vector1< core::Size > bb_residues;

	core::Real energy = scorefxn_->score( pose );
	TR << "starting energy: "<< energy << std::endl;

	core::Size total_rounds = total_glycan_residues_ * rounds_;
	TR << "Total Rounds = "<< total_rounds << " ( " << total_glycan_residues_ << " glycans * " << rounds_ << " )"<<std::endl;

	bool accepted = false;
	mc_->set_last_accepted_pose(pose);
	for ( core::Size round = 1; round <= total_rounds; ++round ) {
		TR << "Round: "<< round << std::endl;
		weighted_random_mover_->apply(pose);
		if ( weighted_random_mover_->get_last_move_status() == protocols::moves::MS_SUCCESS ) {
			core::Real energy = scorefxn_->score(pose);
			if ( TR.Debug.visible() ) {
				TR.Debug << "energy post MC: "<< energy << std::endl;
			}

			if ( pymol_movie_ ) {
				pmm_trials.apply( pose );
			}

			accepted = mc_->boltzmann( pose );

			if ( pymol_movie_ && accepted ) {
				pmm_accepts.apply( pose );
			}

			std::string out = to_string( round )+" "+to_string( energy )+" "+to_string( accepted );
			accept_log_.push_back( out );

			out = "FINAL "+ to_string( round )+" "+to_string( scorefxn_->score( pose ));
			accept_log_.push_back( out );
		} else {
			TR << "Last mover failed.  Continueing!" << std::endl;
			continue;
		}
	}

	mc_->recover_low( pose );

	energy = scorefxn_->score( pose );
	TR << "energy final: "<< energy << std::endl;

	if ( final_min_ ) {
	
		
		min_mover_->apply( pose );
		packer_->apply( pose );
		min_mover_->apply( pose );
		packer_->apply( pose);
		min_mover_->apply( pose );
		
		
		energy = scorefxn_->score( pose );
		TR << "energy final post min: "<< energy << std::endl;
	}

	energy = scorefxn_->score( pose );
	TR << "energy final: "<< energy << std::endl;

	std::string out = "FINAL END "+to_string( scorefxn_->score(pose) );
	accept_log_.push_back(out);
	
	//Add Accept log to pose, have it written:
	for ( core::Size i = 1; i <= accept_log_.size(); ++i ) {
		core::pose::add_comment(pose, "ACCEPT LOG "+utility::to_string(i), accept_log_[i]);
	}
	
	mc_->show_counters();

}

/////////////// Creator ///////////////

protocols::moves::MoverOP
GlycanRelaxMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new GlycanRelaxMover );
}

std::string
GlycanRelaxMoverCreator::keyname() const {
	return GlycanRelaxMoverCreator::mover_name();
}

std::string
GlycanRelaxMoverCreator::mover_name(){
	return "GlycanRelaxMover";
}

} //protocols
} //carbohydrates


