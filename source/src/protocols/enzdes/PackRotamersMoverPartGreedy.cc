/// @file
/// @brief Partly greedy pack rotamers mover. given a set of target residues, choose their neighbors greedily and then
/// @brief pack with usual simulated annealing maintaining these neighbor choices
/// @author Sagar Khare (khares@uw.edu)

// Unit headers
#include <protocols/enzdes/PackRotamersMoverPartGreedy.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/enzdes/PackRotamersMoverPartGreedyCreator.hh>

#include <protocols/moves/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/IGEdgeReweightContainer.hh>
#include <protocols/toolbox/IGEdgeReweighters.hh>
#include <protocols/toolbox/pose_manipulation.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/scoring/Energies.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>
#include <protocols/enzdes/enzdes_util.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <core/scoring/EnergyGraph.hh>

//#include <core/pack/task/operation/TaskOperation.fwd.hh>
//#include <utility/vector0.hh>
#include <utility/vector1.hh>
//#include <utility/string_util.hh>
// Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <utility/sort_predicates.hh>


namespace protocols {
namespace enzdes {

using namespace core;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;

using basic::Warning;
using basic::t_warning;
static basic::Tracer TR("protocols.enzdes.PackRotamersMoverPartGreedy");
static numeric::random::RandomGenerator trg_RG(10805); // <- Magic number, do not change it!!!

std::string
PackRotamersMoverPartGreedyCreator::keyname() const
{
  return PackRotamersMoverPartGreedyCreator::mover_name();
}

protocols::moves::MoverOP
PackRotamersMoverPartGreedyCreator::create_mover() const
{
  return new PackRotamersMoverPartGreedy;
}

PackRotamersMoverPartGreedy::PackRotamersMoverPartGreedy(
		ScoreFunctionOP scorefxn,
		PackerTaskOP task,
		utility::vector1 <core::Size> target_residues
	) :

	protocols::moves::Mover("PackRotamersMoverPartGreedy"),
	scorefxn_repack_ (scorefxn),
	scorefxn_minimize_ (scorefxn),
	task_ (task),
	target_residues_ (target_residues)
{}

PackRotamersMoverPartGreedy::PackRotamersMoverPartGreedy() : 
	protocols::moves::Mover("PackRotamersMoverPartGreedy"),
	scorefxn_repack_(0),
	scorefxn_minimize_(0),
	task_factory_(0),
	target_residues_(0),
	threshold_(0),
	n_best_(0)
{}

std::string
PackRotamersMoverPartGreedyCreator::mover_name()
{
        return "PackRotamersMoverPartGreedy";
}


PackRotamersMoverPartGreedy::~PackRotamersMoverPartGreedy(){}

void
PackRotamersMoverPartGreedy::parse_my_tag(
	TagPtr const tag,
	protocols::moves::DataMap & datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & pose
)
{
	//task operations
	if( tag->hasOption("task_operations") ) task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
  	else task_factory_ = NULL;
	//Scorefunctions
	scorefxn_repack_ = protocols::rosetta_scripts::parse_score_function( tag, "scorefxn_repack", datamap )->clone();
	scorefxn_repack_greedy_ = protocols::rosetta_scripts::parse_score_function( tag, "scorefxn_repack_greedy", datamap )->clone();
	scorefxn_minimize_ = protocols::rosetta_scripts::parse_score_function( tag, "scorefxn_minimize", datamap )->clone();
	//target residues for greedy opt around
	 if( tag->hasOption("target_residues") ) {
      		target_residues_ = core::pose::get_resnum_list(tag, "target_residues",pose);
	}
	use_cstids_ = false;
	if( tag->hasOption("target_cstids") ) {
		cstid_list_ = tag->getOption<std::string>("target_cstids", "" );
		use_cstids_ = true;
    	}
	//distance threshold for neighbor check in greedy opt
	 threshold_  =  tag->getOption<core::Real>("distance_threshold", 8.0 );
	//if choose targets based on N-best residues
	if( tag->hasOption("choose_best_n") ) {
		n_best_ = tag->getOption<core::Size>("choose_best_n", 0 );
	} 
}

std::string
PackRotamersMoverPartGreedy::get_name() const {
  return PackRotamersMoverPartGreedyCreator::mover_name();
}

protocols::moves::MoverOP
PackRotamersMoverPartGreedy::clone() const
{
  return new PackRotamersMoverPartGreedy( *this );
}

protocols::moves::MoverOP
PackRotamersMoverPartGreedy::fresh_instance() const
{
  return new PackRotamersMoverPartGreedy;
}


void
PackRotamersMoverPartGreedy::apply( Pose & pose )
{
  //create task
  using namespace core::pack::task;
  core::pose::Pose greedy_pose = pose;
  TR<<"Creating packer task based on specified task operations..."<< std::endl;
  if ( task_factory_ !=0 ){
  	task_factory_->push_back( new core::pack::task::operation::InitializeFromCommandline );
  	task_ = task_factory_->create_task_and_apply_taskoperations( greedy_pose );
  }
  else{
	core::pack::task::TaskFactoryOP tf = new core::pack::task::TaskFactory();
        task_ = tf->create_task_and_apply_taskoperations( greedy_pose );
  }
  (*scorefxn_repack_)(greedy_pose);
  // make sure we have some target residues
  if (use_cstids_){ // enzdes csts if specified
	utility::vector1< core::Size > trg_res;
	protocols::enzdes::enzutil::get_resnum_from_cstid_list( cstid_list_, greedy_pose, trg_res );
	target_residues_.insert( target_residues_.begin(), trg_res.begin(), trg_res.end() );
        std::unique( target_residues_.begin(), target_residues_.end() );

   }
  if (n_best_>0){
	utility::vector1< core::Size > n_best_res = choose_n_best( greedy_pose, n_best_ );
	target_residues_.insert( target_residues_.begin(), n_best_res.begin(), n_best_res.end() );
        std::unique( target_residues_.begin(), target_residues_.end() );
  }
  runtime_assert(target_residues_.size()>0);
  //randomly shuffle targets
  random_permutation( target_residues_ , trg_RG );
  //Make sure target residues are held fixed
   for( utility::vector1< core::Size >::const_iterator pos_it = target_residues_.begin(); pos_it != target_residues_.end(); ++pos_it ){
	task_->nonconst_residue_task( *pos_it ).prevent_repacking();
   }
  //Convert to Ala
  utility::vector1< core::Size > toAla;
  for (core::Size ii=1;ii<=greedy_pose.total_residue();++ii){
	if (task_->being_designed(ii)) toAla.push_back( ii );
  }
  protocols::toolbox::pose_manipulation::construct_poly_ala_pose( greedy_pose, toAla, true, true, true );
  //Then, greedy opt around user specified target residues, this will modify the task and pose
  greedy_around( greedy_pose, target_residues_, task_, scorefxn_repack_greedy_);
  //greedy_pose.dump_pdb("after_greedy.pdb");
  //Now just replace the target_residues and the greedily optimized positions in original pose with greedy pose, set task to reflect these
  for( utility::vector1< core::Size >::const_iterator pos_it = target_residues_.begin(); pos_it != target_residues_.end(); ++pos_it ){
	pose.replace_residue( *pos_it, greedy_pose.residue(*pos_it), true);
   }
  for( utility::vector1< core::Size >::const_iterator res=restrict_to_repacking_.begin(); res!=restrict_to_repacking_.end(); ++res ){
	if (!(greedy_pose.residue(*res).name3() == "ALA")) pose.replace_residue( *res, greedy_pose.residue(*res), true);
  }
  //old task is no longer valid. create a new one.
  PackerTaskOP task = task_factory_->create_task_and_apply_taskoperations( pose );
  TR<<" Restrcting the following greedily found residues to repack: ";
  for( utility::vector1< core::Size >::const_iterator res=restrict_to_repacking_.begin(); res!=restrict_to_repacking_.end(); ++res )
    {
	 if (!(greedy_pose.residue(*res).name3() == "ALA")) task->nonconst_residue_task( *res ).restrict_to_repacking();
	TR<< *res << "+";
    }
  TR<< std::endl;
  task->initialize_from_command_line().or_include_current( true );
  //Next hand over these to regular PackRotamers mover and we're done!
  protocols::simple_moves::PackRotamersMoverOP pack =  new protocols::simple_moves::PackRotamersMover( scorefxn_repack_ , task );
  pack->apply( pose );
  
}

void
PackRotamersMoverPartGreedy::greedy_around(
	core::pose::Pose & pose,
	utility::vector1<core::Size > const & target_residues,
	core::pack::task::PackerTaskOP task,
	core::scoring::ScoreFunctionCOP scorefxn
)
{	
	using namespace core::pack;
	using namespace core::pack::task;

	core::pose::Pose cur_pose = pose;
	
	
	for( utility::vector1< core::Size >::const_iterator pos_it = target_residues.begin(); pos_it != target_residues.end(); ++pos_it ){
		TR<<"Considering target residue "<< cur_pose.residue( *pos_it ).name3() << *pos_it <<": "<<std::endl;	
		utility::vector1< core::Size > neighbors  = compute_designable_neighbors( *pos_it, task, cur_pose );
		utility::vector1< core::Size > cur_neighbors = neighbors;
		utility::vector1< bool > allow_minimization (pose.total_residue(), false);
		for (utility::vector1< core::Size>::const_iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter){
			allow_minimization[ *iter ] = true;
		}
		allow_minimization[*pos_it] = true;
		// Iterate neighbors.size() times
		TR << "Total number of neighbors are "<<neighbors.size()<<std::endl;
		for (core::Size ii=1;ii<=neighbors.size();++ii){
	 		protocols::toolbox::pose_manipulation::construct_poly_ala_pose( cur_pose, cur_neighbors, false, false, true );	
		//	cur_pose.dump_pdb("polyA_"+utility::to_string(ii)+".pdb");
			TR << "At iteration "<< ii << ", current neighbor size is "<< cur_neighbors.size()<<std::endl;	
			core::Size best_neigh (0);
			//calculate energy in polyAla state before choosing
			core::Real best_energy (cur_pose.energies().total_energies()[core::scoring::total_score]);
			TR <<" Poly Ala energy is :"<< best_energy<<std::endl;
			core::pose::Pose inner_pose = cur_pose; 
			// Iterate over all neighbors at which no choice has been made so far and choose best in the context of current state
		 	for( utility::vector1< core::Size >::const_iterator neigh_it = cur_neighbors.begin(); neigh_it != cur_neighbors.end(); ++neigh_it ){
				core::pose::Pose working_pose = inner_pose;
				// make a task
				core::pack::task::TaskFactoryOP mut_res = new core::pack::task::TaskFactory();
				core::pack::task::PackerTaskOP mutate_residue = mut_res->create_task_and_apply_taskoperations( working_pose );
			        mutate_residue->initialize_from_command_line().or_include_current( true );
				// upweight interactions of this position to those within the cluster (to favor intracluster interactions)
				utility::vector1< core::Size > upweight_1;
				upweight_1.push_back( *neigh_it );
				utility::vector1< core::Size > upweight_2 = neighbors;
				upweight_2.push_back( *pos_it );
				//upweight interactions between this position and neighbors
				IGEdgeReweighterOP ig_up = new protocols::toolbox::ResidueGroupIGEdgeUpweighter( 2.5, upweight_1 , upweight_2 );
				mutate_residue->set_IGEdgeReweights()->add_reweighter( ig_up );
				// fix everything else except this position
				utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, true );
				
				for (core::Size i=1; i<=pose.total_residue(); ++i){
				if (i== *neigh_it) mutate_residue->nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas );
				else mutate_residue->nonconst_residue_task( i ).prevent_repacking();
				}
				// Now run design on one position
				protocols::simple_moves::PackRotamersMoverOP prm =  new protocols::simple_moves::PackRotamersMover( scorefxn, mutate_residue, 1 );
				prm->apply( working_pose );
				//Minimize the sidechains of all designed residues so far
				core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap();
				movemap->set_jump( false );
  				movemap->set_bb( false );
				movemap->set_chi( allow_minimization );
				protocols::simple_moves::MinMoverOP minmover = new protocols::simple_moves::MinMover( movemap, scorefxn_minimize_, "dfpmin_armijo_nonmonotone_atol", 0.02, true /*use_nblist*/);
				minmover->apply(working_pose);

				//Check energy
				core::Real this_energy = working_pose.energies().total_energies()[core::scoring::total_score];
				//TR << "Energy for best residue at position "<< *neigh_it << "is "<< this_energy << std::endl;
				// see if this energy is best so far. if so replace this neighbor in cur_pose 
				if ( (this_energy<best_energy) || (best_neigh==0) ) { 
					//TR<< "Best neighbor at position "<<*neigh_it <<" is "<< working_pose.residue( *neigh_it ).name3() << " with energy "<< this_energy<<std::endl;
					cur_pose = working_pose;
					best_neigh = *neigh_it;
					best_energy = this_energy;
				}
			} // all currently unchosen neighbors so far
			//Now that a choice has been made at the position best_neigh, inform cur_neighbors and task to take this position out from future design 
			update_task_and_neighbors( best_neigh, task, cur_neighbors );
		} //iterate neighbors.size() times around each target residue
	} //all target residues
	pose = cur_pose;
} // greedy around

utility::vector1< core::Size > 
PackRotamersMoverPartGreedy::compute_designable_neighbors(
core::Size const & position,
core::pack::task::PackerTaskCOP task,
core::pose::Pose const & pose
)
{
	utility::vector1< core::Size > neighbors;	
	TR << "neighbors of position "<< position << " are: ";
	core::conformation::Residue const res_pos (pose.residue( position ));
	for (core::Size ii=1;ii<=pose.total_residue();++ii){
		if ( pose.residue( ii ).is_protein() && task->being_designed( ii ) ){
			core::conformation::Residue const res_ii (pose.residue( ii ));
			core::Real distance = res_ii.xyz( res_ii.nbr_atom() ).distance( res_pos.xyz(res_pos.nbr_atom()) );
			if (distance<=threshold_)  {
				neighbors.push_back( ii );
				TR <<ii<<" + "; 
			}
		}
	}
	TR<<std::endl;
//	TR<< "In function total neighbors are"<< neighbors.size() <<std::endl;
	return neighbors;
}

void
PackRotamersMoverPartGreedy::update_task_and_neighbors(
core::Size const & best_neigh,
core::pack::task::PackerTaskOP task,
utility::vector1< core::Size > & current_neighbors
)
{
	task->nonconst_residue_task( best_neigh ).restrict_to_repacking();
	restrict_to_repacking_.push_back( best_neigh );
	utility::vector1< core::Size > temp_neighbors;
	for( utility::vector1< core::Size >::const_iterator neigh_it = current_neighbors.begin(); neigh_it != current_neighbors.end(); ++neigh_it ){
		if (*neigh_it != best_neigh) temp_neighbors.push_back(*neigh_it);
		else TR<<" Updated current neighbors by removing "<<*neigh_it << " from it."<<std::endl; 
	}
	current_neighbors = temp_neighbors;


}

void
PackRotamersMoverPartGreedy::task_factory( core::pack::task::TaskFactoryOP p ) {
  task_factory_ = p;
}

void 
PackRotamersMoverPartGreedy::task( core::pack::task::PackerTaskOP task ){
	task_ = task;
}

void 
PackRotamersMoverPartGreedy::target_residues( utility::vector1< core::Size > & trg_res ){
	target_residues_ = trg_res;
}

utility::vector1<core::Size>  
PackRotamersMoverPartGreedy::choose_n_best( core::pose::Pose const & pose , core::Size const & n_best ){

	using namespace core::graph;
	using namespace core::scoring;
	std::list< std::pair< core::Size, core::Real > > residue_energies;
	utility::vector1<core::Size> chosen_residues;
	core::pose::Pose nonconst_pose( pose );
	(*scorefxn_repack_)(nonconst_pose);
 	EnergyMap const curr_weights = nonconst_pose.energies().weights();
	core::Size res (1);
	for (core::Size ii=1;ii<=pose.total_residue();++ii){
  		if (pose.residue( ii ).is_ligand()) res = ii;
	}
	//Fill residue_energies of interface residues by traversing energy graph
	for( EdgeListConstIterator egraph_it = nonconst_pose.energies().energy_graph().get_node( res )->const_edge_list_begin(); egraph_it != nonconst_pose.energies().energy_graph().get_node( res )->const_edge_list_end(); ++egraph_it){
		core::Size const int_resi = (*egraph_it)->get_other_ind( res );
		EnergyEdge const * Eedge = static_cast< EnergyEdge const * > (*egraph_it);
		core::Real intE = Eedge->dot( curr_weights );
		residue_energies.push_back( std::make_pair( int_resi, intE));	
  	}//for each egraph_it

	//Sort and return n_best
	residue_energies.sort( utility::SortSecond< core::Size, core::Real >() );
	core::Size counter(1);
	for( std::list< std::pair<core::Size, core::Real> >::iterator it = residue_energies.begin();counter<=n_best; ++it){
		chosen_residues.push_back( it->first );
		TR << "Chose residue "<<it->first<<" with interface energy "<<it->second <<std::endl; 
		counter++;
	}
	return chosen_residues;
}
/*core::PackerEnergy PackRotamersMoverPartGreedy::run_with_ig( Pose & pose, utility::vector0< int > rot_to_pack, InteractionGraphBaseOP ig) const
{
	return pack_rotamers_run( pose, task(), rotamer_sets(), ig, rot_to_pack );
}
*/

}//moves
}//protocols
