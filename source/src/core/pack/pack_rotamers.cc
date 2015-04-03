// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/pack_rotamers.cc
/// @brief  pack rotamers module
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/pack_rotamers.hh>

// Package Headers
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/task/IGEdgeReweightContainer.hh>

#include <core/pack/annealer/AnnealerFactory.hh>
#include <core/pack/annealer/SimAnnealerBase.hh>
#include <core/pack/interaction_graph/InteractionGraphFactory.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>


// Project Headers
#include <core/conformation/Residue.hh>
#include <core/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// util
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>

// option key includes

#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


using namespace ObjexxFCL;

namespace core {
namespace pack {

using core::conformation::symmetry::SymmetryInfoCOP;
using core::conformation::symmetry::SymmetricConformation;

static thread_local basic::Tracer tt( "core.pack.pack_rotamers", basic::t_info );

// @details Wraps the two very distinct and separate stages of rotamer packing, which are factored so that they may be called asynchronously.  Use this wrapper as a base model for higher-level packing routines (such as pack_rotamers_loop)
void
pack_rotamers(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task
)
{
	using namespace interaction_graph;
	using namespace rotamer_set;

	//fpd safety check for symmetry
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		symmetric_pack_rotamers( pose, scfxn, task );
		return;
	}

	PROF_START( basic::PACK_ROTAMERS );

	pack_scorefxn_pose_handshake( pose, scfxn);

	//replace this with RotSetsFactory
	rotamer_set::RotamerSetsOP rotsets( new rotamer_set::RotamerSets() );
	InteractionGraphBaseOP ig = NULL;

	pack_rotamers_setup( pose, scfxn, task, rotsets, ig ); //TODO: make compatible with D-amino acids
	pack_rotamers_run( pose, task, rotsets, ig ); //TODO: make compatible with D-amino acids

	// rescore here to make the state of the Energies good.
	scfxn( pose );

	PROF_STOP ( basic::PACK_ROTAMERS );
}

// @details run the FixbbSimAnnealer multiple times using the same InteractionGraph, storing the results
void
pack_rotamers_loop(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	Size const nloop
)
{
	utility::vector1< std::pair< Real, std::string > > results;
	pack_rotamers_loop( pose, scfxn, task, nloop, results );
}

// @details run the FixbbSimAnnealer multiple times using the same InteractionGraph, storing the results
void
pack_rotamers_loop(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	Size const nloop,
	utility::vector1< std::pair< Real, std::string > > & results
)
{
	utility::vector1< pose::PoseOP > pose_list;
	pack_rotamers_loop( pose, scfxn, task, nloop, results, pose_list );
}

// @details run the FixbbSimAnnealer multiple times using the same InteractionGraph, storing the results
void
pack_rotamers_loop(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	Size const nloop,
	utility::vector1< std::pair< Real, std::string > > & results,
	utility::vector1< pose::PoseOP > & pose_list
)
{
	using namespace ObjexxFCL::format;
	using namespace interaction_graph;
	using namespace rotamer_set;

	rotamer_set::RotamerSetsOP rotsets( new rotamer_set::RotamerSets() );
	InteractionGraphBaseOP ig = NULL;
	pack_rotamers_setup( pose, scfxn, task, rotsets, ig );

	setup_IG_res_res_weights( pose, task, rotsets, ig );

	Real best_bestenergy( 0.0 );
	pose::Pose best_pose;
	best_pose = pose;

	for ( Size run(1); run <= nloop; ++run ) {

		Real bestenergy( pack_rotamers_run( pose, task, rotsets, ig ) );

		Real const final_score( scfxn( pose ) );
		// show the resulting sequence
		std::string final_seq;
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( task->design_residue(i) ) final_seq+= pose.residue(i).name1();
		}
		if ( final_seq.size() == 0 ) final_seq = "-";

		tt << "packloop: " << I(4,run) << " rescore-deltaE: " << F(9,3,final_score - bestenergy ) <<
			" seq: " << final_seq << " simannealerE: " << bestenergy << " rescoreE: " << final_score << std::endl;

		results.push_back( std::make_pair( final_score, pose.sequence() ) );
		pose_list.push_back( core::pose::PoseOP( new pose::Pose( pose ) ) ); //saves a copy.
		if ( run == 1 || bestenergy < best_bestenergy ) {
			best_pose = pose;
			best_bestenergy = bestenergy;
		}
	}
	pose = best_pose;
}

// @details get rotamers, compute energies
void
pack_rotamers_setup(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	rotamer_set::RotamerSetsOP rotsets,
	interaction_graph::InteractionGraphBaseOP & ig
)
{
	using namespace interaction_graph;

	pack_scorefxn_pose_handshake( pose, scfxn);

	pose.update_residue_neighbors();

	scfxn.setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );

	graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scfxn, task );

	rotsets->set_task( task );

	rotsets->build_rotamers( pose, scfxn, packer_neighbor_graph );
	rotsets->prepare_sets_for_packing( pose, scfxn );

	if ( basic::options::option[ basic::options::OptionKeys::packing::dump_rotamer_sets ] ) { // hacking
		static int counter(0);
		++counter;
		std::string const filename( "rotset"+lead_zero_string_of( counter,4 )+".pdb" );
		tt << "dump rotsets: " << filename << std::endl;
		rotsets->dump_pdb( pose, filename );
	}

	ig = InteractionGraphFactory::create_interaction_graph( *task, *rotsets, pose, scfxn );

	tt << "built " << rotsets->nrotamers() << " rotamers at "
		 << rotsets->nmoltenres() << " positions." << std::endl;

	//for ( Size i=1; i<= pose.total_residue(); ++i ) {
	//	if ( task->design_residue(i) ) tt << "designing at position: " << i << std::endl;
	//}

	PROF_START( basic::GET_ENERGIES );
	rotsets->compute_energies( pose, scfxn, packer_neighbor_graph, ig );
	PROF_STOP( basic::GET_ENERGIES );

	tt << "IG: " << ig->getTotalMemoryUsage() << " bytes" << std::endl;

}


// PyRosetta compatible version
interaction_graph::InteractionGraphBaseOP
pack_rotamers_setup(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	rotamer_set::RotamerSetsOP rotsets)
{
	interaction_graph::InteractionGraphBaseOP ig;
	pack_rotamers_setup(pose, scfxn, task, rotsets, ig);
	return ig;
}


/// @brief upweights certain edges in the interaction graph if this is specified in the task
void
setup_IG_res_res_weights(
	pose::Pose const & pose,
	task::PackerTaskCOP task,
	rotamer_set::RotamerSetsCOP rotsets,
	interaction_graph::InteractionGraphBaseOP ig
)
{
	task::IGEdgeReweightContainerCOP edge_reweights = task->IGEdgeReweights();

	if( edge_reweights ){

		for( Size ii = 1; ii<= rotsets->nmoltenres(); ++ii){
			Size const res1id = rotsets->moltenres_2_resid( ii );

			for( ig->reset_edge_list_iterator_for_node( ii ); !ig->edge_list_iterator_at_end(); ig->increment_edge_list_iterator() ){

				interaction_graph::EdgeBase const & edge( ig->get_edge() );
				Size const other_node = edge.get_other_ind( ii );
				if ( other_node < ii ) continue; // only deal with upper edges
				Size const res2id = rotsets->moltenres_2_resid( other_node );
				ig->set_edge_weight( ii, other_node, edge_reweights->res_res_weight( pose, *task, res1id, res2id ) );
			}
		}
	} // if( edge_reweights )

}//setup_IG_res_res_weights


// @details as simple as possible -- runs simulated annealing and then places
// the optimal rotamers onto the backbone of the input pose.
Real
pack_rotamers_run(
	pose::Pose & pose,
	task::PackerTaskCOP task,
	rotamer_set::FixbbRotamerSetsCOP rotsets,
	interaction_graph::InteractionGraphBaseOP ig,
	utility::vector0< int > rot_to_pack // defaults to an empty vector (no effect)
)
{
	using namespace ObjexxFCL::format;

	FArray1D_int bestrotamer_at_seqpos( pose.total_residue() );
	core::PackerEnergy bestenergy( 0.0 );

	pack_rotamers_run( pose, task, rotsets, ig, rot_to_pack, bestrotamer_at_seqpos, bestenergy );

	// place new rotamers on input pose
	for ( uint ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
		uint iiresid = rotsets->moltenres_2_resid( ii );
		uint iibestrot = rotsets->rotid_on_moltenresidue( bestrotamer_at_seqpos( iiresid ) );
		conformation::ResidueCOP bestrot( rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( iibestrot ) );

//		for( Size i = 1 ; i <= bestrot->natoms() ; ++i ) {
//			std::cout << "bestrot " << bestrot->atom_name(i) <<
//					" x " << bestrot->xyz(i)[0] <<
//					" y " << bestrot->xyz(i)[1] <<
//					" z " << bestrot->xyz(i)[2] << std::endl;
//		}

		conformation::ResidueOP newresidue( bestrot->create_residue() );
		pose.replace_residue ( iiresid, *newresidue, false );
	}
	return bestenergy;
}

/// @brief Runs simulated annealing and returns the
void
pack_rotamers_run(
	pose::Pose const & pose,
	task::PackerTaskCOP task,
	rotamer_set::FixbbRotamerSetsCOP rotsets,
	interaction_graph::InteractionGraphBaseOP ig,
	utility::vector0< int > rot_to_pack,
	ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy
)
{
	using namespace annealer;

	bool start_with_current = false;
	FArray1D_int current_rot_index( pose.total_residue(), 0 );
	bool calc_rot_freq = false;
	FArray1D< core::PackerEnergy > rot_freq( ig->get_num_total_states(), 0.0 );

	// too many parameters for annealer's constructor! should replace with a task only;
	// the annealer should then provide read access to the data its collected after
	// annealing has completed. some data it won't bother collecting because the task
	// did not instruct it to.

	/// Parameters passed by reference in task's constructor to which it writes at the
	/// completion of sim annealing.

	SimAnnealerBaseOP annealer = AnnealerFactory::create_annealer(
		task, rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig,
		rotsets, current_rot_index, calc_rot_freq, rot_freq );

	// Following additional initialization would be cleaner if above suggestion
	// regarding the removal of "too many parameters" is implemented properly.
	if ( task->low_temp()  > 0.0 ) annealer->set_lowtemp(  task->low_temp()  );
	if ( task->high_temp() > 0.0 ) annealer->set_hightemp( task->high_temp() );
	annealer->set_disallow_quench( task->disallow_quench() );

	PROF_START( basic::SIMANNEALING );
	annealer->run();
	PROF_STOP( basic::SIMANNEALING );
}

void
symmetric_pack_rotamers(
  pose::Pose & pose,
  scoring::ScoreFunction const & scfxn,
  task::PackerTaskCOP non_symmetric_task
)
{
  using namespace interaction_graph;
  using namespace rotamer_set;
  PROF_START( basic::PACK_ROTAMERS );

  task::PackerTaskCOP task = make_new_symmetric_PackerTask_by_requested_method(pose,non_symmetric_task);

  pack_scorefxn_pose_handshake( pose, scfxn);

  //replace this with RotSetsFactory
  rotamer_set::symmetry::SymmetricRotamerSetsOP rotsets( new rotamer_set::symmetry::SymmetricRotamerSets() );
  InteractionGraphBaseOP ig = NULL;

  symmetric_pack_rotamers_setup( pose, scfxn, task, rotsets, ig );
  symmetric_pack_rotamers_run( pose, task, rotsets, ig );

	// rescore here to make the state of the Energies good.
  scfxn( pose );

  PROF_STOP ( basic::PACK_ROTAMERS );
}

void
symmetric_pack_rotamers_setup(
  pose::Pose & pose,
  scoring::ScoreFunction const & scfxn,
  task::PackerTaskCOP task,
  rotamer_set::symmetry::SymmetricRotamerSetsOP rotsets,
  interaction_graph::InteractionGraphBaseOP & ig
)
{
  using namespace interaction_graph;

  pack_scorefxn_pose_handshake( pose, scfxn);

  pose.update_residue_neighbors();
  scfxn.setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );

  graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scfxn, task );

  //rotsets->set_symmetrical_task( task, pose );
  rotsets->set_task( task );
  rotsets->build_rotamers( pose, scfxn, packer_neighbor_graph );
  rotsets->prepare_sets_for_packing( pose, scfxn );

  if ( basic::options::option[ basic::options::OptionKeys::packing::dump_rotamer_sets ] ) { // hacking
    static int counter(0);
    ++counter;
    std::string const filename( "rotset"+lead_zero_string_of( counter,4 )+".pdb" );
    tt << "dump rotsets: " << filename << std::endl;
    rotsets->dump_pdb( pose, filename );
  }

  ig = InteractionGraphFactory::create_interaction_graph( *task, *rotsets, pose, scfxn );

  tt << "built " << rotsets->nrotamers() << " rotamers at "
     << rotsets->nmoltenres() << " positions." << std::endl;

  //for ( Size i=1; i<= pose.total_residue(); ++i ) {
  //  if ( task->design_residue(i) ) tt << "designing at position: " << i << std::endl;
  //}

  PROF_START( basic::GET_ENERGIES );
  rotsets->compute_energies( pose, scfxn, packer_neighbor_graph, ig );
  PROF_STOP( basic::GET_ENERGIES );

  tt << "IG: " << ig->getTotalMemoryUsage() << " bytes" << std::endl;

}

// PyRosetta compatible version
interaction_graph::InteractionGraphBaseOP
symmetric_pack_rotamers_setup(
  pose::Pose & pose,
  scoring::ScoreFunction const & scfxn,
  task::PackerTaskCOP task,
  rotamer_set::symmetry::SymmetricRotamerSetsOP rotsets
)
{
	interaction_graph::InteractionGraphBaseOP ig;
	symmetric_pack_rotamers_setup(pose, scfxn, task, rotsets, ig);
	return ig;
}


// @details as simple as possible
Real
symmetric_pack_rotamers_run(
	pose::Pose & pose,
	task::PackerTaskCOP task,
	rotamer_set::symmetry::SymmetricRotamerSetsCOP rotsets,
	interaction_graph::InteractionGraphBaseOP ig,
	utility::vector0< int > rot_to_pack // defaults to an empty vector (no effect)
)
{
	using namespace ObjexxFCL::format;
	using namespace annealer;

	// too many parameters for annealer's constructor! should replace with a task only;
	// the annealer should then provide read access to the data its collected after
	// annealing has completed. some data it won't bother collecting because the task
	// did not instruct it to.

	// Symmetry info
		SymmetricConformation const & SymmConf (
    dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
  	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );


	/// Parameters passed by reference in task's constructor to which it writes at the
	/// completion of sim annealing.
	FArray1D_int bestrotamer_at_seqpos( pose.total_residue() );
	core::PackerEnergy bestenergy;
	bool start_with_current = false;
	FArray1D_int current_rot_index( pose.total_residue(), 0 );
	bool calc_rot_freq = false;
	FArray1D< core::PackerEnergy > rot_freq( ig->get_num_total_states(), 0.0 );

	SimAnnealerBaseOP annealer = AnnealerFactory::create_annealer(
		task, rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig,
		rotsets, current_rot_index, calc_rot_freq, rot_freq );

//	SimAnnealerBaseOP annealer = new annealer::symmetry::SymFixbbSimAnnealer(
//      rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig,
//      rotsets, current_rot_index, calc_rot_freq, rot_freq, symm_info );

	// Following additional initialization would be cleaner if above suggestion
	// regarding the removal of "too many parameters" is implemented properly.
	if ( task->low_temp()  > 0.0 ) annealer->set_lowtemp(  task->low_temp()  );
	if ( task->high_temp() > 0.0 ) annealer->set_hightemp( task->high_temp() );
	annealer->set_disallow_quench( task->disallow_quench() );

	PROF_START( basic::SIMANNEALING );
	annealer->run();
	PROF_STOP( basic::SIMANNEALING );

	// place new rotamers on input pose
	for ( uint ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
		uint iiresid = rotsets->moltenres_2_resid( ii );
		uint iibestrot = rotsets->rotid_on_moltenresidue( bestrotamer_at_seqpos( iiresid ) );
		conformation::ResidueCOP bestrot( rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( iibestrot ) );

//		for( Size i = 1 ; i <= bestrot->natoms() ; ++i ) {
//			std::cout << "bestrot " << bestrot->atom_name(i) <<
//					" x " << bestrot->xyz(i)[0] <<
//					" y " << bestrot->xyz(i)[1] <<
//					" z " << bestrot->xyz(i)[2] << std::endl;
//		}

		conformation::ResidueOP newresidue( bestrot->create_residue() );
		pose.replace_residue ( iiresid, *newresidue, false );

	//fpd replace residue is symmetric now
	//for ( std::vector< Size>::const_iterator
	//      clone     = symm_info.bb_clones( iiresid ).begin(),
	//      clone_end = symm_info.bb_clones( iiresid ).end();
	//      clone != clone_end; ++clone ){
	//	conformation::ResidueOP sym_rsd = newresidue->clone();
	// 	sym_rsd->orient_onto_residue(pose.residue( *clone) );
	// 	pose.replace_residue ( *clone, *sym_rsd, false );
	// }
	}
	return bestenergy;
}


} // namespace pack
} // namespace core
