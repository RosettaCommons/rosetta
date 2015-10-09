// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pack_daemon/PackDaemon.cc
/// @brief  Implementation for class PackDaemon
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

/// MPI Headers
#ifdef USEMPI
#include <mpi.h>
#endif

// Unit headers
#include <protocols/pack_daemon/PackDaemon.hh>

// Package headers
#include <protocols/pack_daemon/EntityCorrespondence.hh>
#include <protocols/pack_daemon/util.hh>

// Project headers
#include <core/chemical/ChemicalManager.fwd.hh>

#ifdef USEMPI
#include <core/io/pdb/pose_io.hh>
#endif
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/annealer/FASTERAnnealer.hh>
#include <core/pack/interaction_graph/DensePDInteractionGraph.hh>
#include <core/pack/interaction_graph/DoubleDensePDInteractionGraph.hh>
#include <core/pack/interaction_graph/FASTERInteractionGraph.hh>
#include <core/pack/interaction_graph/FixedBBInteractionGraph.hh>
//#include <core/pack/interaction_graph/LazyInteractionGraph.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSubsets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>

#include <protocols/multistate_design/MultiStatePacker.hh>

#include <utility/exit.hh>
#include <utility/integer_mapping.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#ifdef USEMPI
#include <utility/mpi_util.hh>
#endif


namespace protocols {
namespace pack_daemon {

static THREAD_LOCAL basic::Tracer TR( "protocols.pack_daemon.PackDaemon" );

PackDaemon::PackDaemon() :
	include_background_energies_( true ),
	background_energies_( 0.0 ),
	setup_complete_( false ),
	best_assignment_valid_( false )
{
	best_assignment_.second = 12345;
	last_assignment_.second = 12345;
}

PackDaemon::~PackDaemon() {}

// Initialize the PackDaemon with the appropriate data before
// calling setup().
void PackDaemon::set_pose_and_task( Pose const & pose, PackerTask const & task )
{
	pose_ = PoseOP( new Pose( pose ) );
	task_ = task.clone();
	for ( Size ii = 1; ii <= task.total_residue(); ++ii ) {
		task_->nonconst_residue_task( ii ).and_extrachi_cutoff( 1 );
	}
	task_->or_double_lazy_ig( true );
	task_->set_bump_check( false ); /// bump check must be disabled.
	setup_complete_ = false;
}

void PackDaemon::set_score_function( ScoreFunction const & sfxn )
{
	score_function_ = sfxn.clone();
	setup_complete_ = false;
}

void PackDaemon::set_entity_correspondence( EntityCorrespondence const &  ec )
{
	if ( ! pose_ ) {
		utility_exit_with_message( "PackDaemon::set_entity_correspondence may only be called after the set_pose_and_task" );
	}
	if ( ec.num_residues() != pose_->total_residue() ) {
		utility_exit_with_message( "Num residue disgreement between input EntityCorrespondence and existing pose_" );
	}
	if ( ! correspondence_ ) {
		correspondence_ = EntityCorrespondenceOP( new EntityCorrespondence( ec ) );
		setup_complete_ = false;
	}
	correspondence_->set_pose( pose_ ); // discard the pose that's already being pointed to by the ec
}

void PackDaemon::set_dlig_nmeg_limit( Size setting )
{
	TR << "Setting dlig nmeg limit" << std::endl;
	task_->decrease_double_lazy_ig_memlimit( 1024 * 1024 * setting );
	setup_complete_ = false;
}

void PackDaemon::set_include_background_energies( bool setting )
{
	include_background_energies_ = setting;
	setup_complete_ = false;
}

void PackDaemon::setup()
{
	if ( !pose_ || !task_ ) {
		utility_exit_with_message( "PackDaemon::set_pose_and_task must be called before setup() can be" );
	}
	if ( !score_function_ ) {
		utility_exit_with_message( "PackDaemon::set_score_function must be called before setup() can be" );
	}
	if ( !correspondence_ ) {
		utility_exit_with_message( "PackDaemon::set_entity_correspondence must be called before setup() can be" );
	}
	TR << "PackDaemon::setup()" << std::endl;

	rot_sets_ = RotamerSetsOP( new RotamerSets );
	core::pack::interaction_graph::AnnealableGraphBaseOP ig;
	core::pack::pack_rotamers_setup( *pose_, *score_function_, task_, rot_sets_, ig );

	ig_ = utility::pointer::dynamic_pointer_cast< core::pack::interaction_graph::FixedBBInteractionGraph > ( ig );
	if ( ! ig_ )  {
		throw utility::excn::EXCN_Msg_Exception( "Interaction graph returned by pack_rotamers_setup is not a"
			" fixed-backbone two-body interaction graph.  Cannot continue" );
	}

	best_assignment_.first.resize( rot_sets_->nmoltenres() );
	last_assignment_.first.resize( rot_sets_->nmoltenres() );

	//repacker_ = new BasicSimAnnealerRepacker( pose_, task_, ig_, rot_sets_ );
	//repacker_ = new DenseIGRepacker( pose_, task_, ig_, rot_sets_ );
	//repacker_ = new DoubleDenseIGRepacker( pose_, task_, ig_, rot_sets_ );
	//FASTER_IG_Repacker * temp;

	/* temp = new FASTER_IG_Repacker( pose_, task_, ig_, rot_sets_ );
	temp->set_num_sa( 8 );
	repacker2_ = temp;

	temp = new FASTER_IG_Repacker( pose_, task_, ig_, rot_sets_ );
	temp->set_num_sa( 4 );
	repacker3_ = temp;

	temp = new FASTER_IG_Repacker( pose_, task_, ig_, rot_sets_ );
	temp->set_num_sa( 2 );
	repacker4_ = temp;

	temp = new FASTER_IG_Repacker( pose_, task_, ig_, rot_sets_ );
	temp->set_sa_scale( 0.025 );
	temp->set_num_sa( 8 );
	repacker5_ = temp;

	temp = new FASTER_IG_Repacker( pose_, task_, ig_, rot_sets_ );
	temp->set_sa_scale( 0.025 );
	temp->set_num_sa( 4 );
	repacker6_ = temp;

	temp = new FASTER_IG_Repacker( pose_, task_, ig_, rot_sets_ );
	temp->set_sa_scale( 0.025 );
	temp->set_num_sa( 2 );
	repacker7_ = temp;

	temp = new FASTER_IG_Repacker( pose_, task_, ig_, rot_sets_ );
	temp->set_ciBR_only( true );
	repacker_ = temp;*/

	DenseIGRepackerOP mca_repacker( new DenseIGRepacker( pose_, task_, ig_, rot_sets_ ) );
	mca_repacker->set_MCA();
	repacker_ = mca_repacker;

	if ( include_background_energies_ ) calculate_background_energies();

	setup_complete_ = true;
}

/// @brief Repack the structure with the Entity
/// This function proceeds in two steps: it creates a list of
/// rotamer indices to be used during the repacking, and then
/// it uses that list to repack the rotamers.  The first step
/// is taken care of by the select_rotamer_subset method.
PackDaemon::Real
PackDaemon::compute_energy_for_assignment( Entity const & entity )
{
	if ( ! last_entity_ ) {
		last_entity_ = entity.clone();
	} else {
		(*last_entity_) = entity;
	}

	utility::vector0< int > rot_to_pack( select_rotamer_subset( entity ) );

	//clock_t starttime, stoptime;

	//starttime = clock();
	last_assignment_ = repacker_->repack( rot_to_pack );
	last_assignment_.second += background_energies_; // add in the background energies to give the total energy for the pose.

	//stoptime = clock();

	/*TR << "repacker_ " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " " << last_assignment_.second << std::endl;

	starttime = clock();
	last_assignment_ = repacker2_->repack( rot_to_pack );
	stoptime = clock();

	TR << "repacker2_ " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " " << last_assignment_.second << std::endl;

	starttime = clock();
	last_assignment_ = repacker3_->repack( rot_to_pack );
	stoptime = clock();

	TR << "repacker3_ " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " " << last_assignment_.second << std::endl;

	starttime = clock();
	last_assignment_ = repacker4_->repack( rot_to_pack );
	stoptime = clock();

	TR << "repacker4_ " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " " << last_assignment_.second << std::endl;

	starttime = clock();
	last_assignment_ = repacker5_->repack( rot_to_pack );
	stoptime = clock();

	TR << "repacker5_ " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " " << last_assignment_.second << std::endl;

	starttime = clock();
	last_assignment_ = repacker6_->repack( rot_to_pack );
	stoptime = clock();

	TR << "repacker6_ " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " " << last_assignment_.second << std::endl;

	starttime = clock();
	last_assignment_ = repacker7_->repack( rot_to_pack );
	stoptime = clock();

	TR << "repacker7_ " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " " << last_assignment_.second << std::endl;

	starttime = clock();
	last_assignment_ = repacker8_->repack( rot_to_pack );
	stoptime = clock();

	TR << "repacker8_ " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " " << last_assignment_.second << std::endl;*/

	if ( ! best_assignment_valid_ || last_assignment_.second < best_assignment_.second ) {
		best_assignment_ = last_assignment_;
		best_assignment_valid_ = true;
		if ( ! best_entity_ ) {
			best_entity_ = entity.clone();
		} else {
			(*best_entity_) = entity;
		}
	}

	return last_assignment_.second;
}

utility::vector0< int >
PackDaemon::select_rotamer_subset( Entity const & entity ) const
{
	using namespace core::pack::rotamer_set;
	using namespace protocols::multistate_design;

	utility::vector0< Size > rotamer_subset;
	for ( Size ii = 1; ii <= rot_sets_->nmoltenres(); ++ii ) {
		Size ii_resid = rot_sets_->moltenres_2_resid( ii );
		Size ii_entity_id = correspondence_->entity_for_residue( ii_resid );
		Size ii_offset = rot_sets_->nrotamer_offset_for_moltenres( ii );
		RotamerSetCOP ii_rotamers = rot_sets_->rotamer_set_for_moltenresidue( ii );
		//std::cout << "Entity id for " << ii_resid << " " << ii_entity_id << std::endl;
		Size ii_rotamers_appended( 0 );
		if (  ii_entity_id == 0 ) {
			// include all the rotamers for this residue
			for ( Size jj = 1, ii_nrots = ii_rotamers->num_rotamers(); jj <= ii_nrots; ++jj ) {
				rotamer_subset.push_back( ii_offset + jj );
				++ii_rotamers_appended;
			}
			if ( ii_rotamers_appended == 0 ) {
				throw utility::excn::EXCN_Msg_Exception( "Failed to find any rotamers for residue "
					+ utility::to_string( ii_resid ) + "." );
			}
		} else {
			PosTypeCOP pt_ptr( utility::pointer::dynamic_pointer_cast< protocols::multistate_design::PosType const > ( entity.traits()[ ii_entity_id ] ));
			if ( ! pt_ptr ) {
				utility_exit_with_message( "Failed to downcast entity trait[" +
					utility::to_string( ii ) + "] to PosTye; trait name: " +
					entity.traits()[ ii ]->name() );
			}
			core::chemical::AA entity_aa( pt_ptr->type() );
			// include only the rotamers from the amino acid specified in this residue's entity element
			for ( Size jj = 1, ii_nrots = ii_rotamers->num_rotamers(); jj <= ii_nrots; ++jj ) {
				if ( ii_rotamers->rotamer( jj )->aa() == entity_aa ) {
					rotamer_subset.push_back( ii_offset + jj );
					++ii_rotamers_appended;
				}
			}

			if ( ii_rotamers_appended == 0 ) {
				throw utility::excn::EXCN_Msg_Exception( "Failed to find any rotamers for residue "
					+ utility::to_string( ii_resid ) + " corresponding to entity "
					+ utility::to_string( ii_entity_id ) + " when looking for "
					+ utility::to_string( entity_aa ) + " rotamers." );
			}
		}
	}
	return rotamer_subset;
}

void PackDaemon::mark_last_entity_as_important()
{

	EntityElements traits_copy( last_entity_->traits().size() );
	for ( Size ii = 1; ii <= traits_copy.size(); ++ii ) {
		traits_copy[ ii ] = last_entity_->traits()[ ii ]->clone();
	}

	prev_state_hash_[ traits_copy ] = last_assignment_;

}

void PackDaemon::mark_entity_as_unimportant( Entity const & ent )
{
	prev_state_hash_.erase( ent.traits() );
}


PackDaemon::PoseCOP                 PackDaemon::pose() const { return pose_; }
PackDaemon::ScoreFunctionCOP        PackDaemon::score_function() const { return score_function_; }
PackDaemon::PackerTaskCOP           PackDaemon::task() const { return task_; }
EntityCorrespondenceCOP             PackDaemon::correspondence() const { return correspondence_; }
PackDaemon::FixedBBInteractionGraphCOP PackDaemon::ig() const { return ig_; }
PackDaemon::RotamerSetsCOP          PackDaemon::rot_sets() const { return rot_sets_; }

PackDaemon::RotamerAssignmentAndEnergy const &
PackDaemon::best_assignment() const { return best_assignment_; }

PackDaemon::RotamerAssignmentAndEnergy const &
PackDaemon::last_assignment() const { return last_assignment_; }

PackDaemon::PoseOP PackDaemon::recreate_pose_for_entity( Entity const & ent ) const
{
	EntityToRotamerHash::const_iterator iter = prev_state_hash_.find( ent.traits() );

	if ( iter == prev_state_hash_.end() ) {
		TR << "Failed to find entity in entity-hash: " << ent << std::endl;
		print_entity_history();
		/// Failed to find this entity -- that's a problem!
		return 0;
	}
	RotamerAssignmentAndEnergy const & rotamer_assignment( iter->second );
	PoseOP return_pose( new Pose( *pose_ ) );
	for ( Size ii = 1; ii <= rot_sets_->nmoltenres(); ++ii ) {
		Size iiresid = rot_sets_->moltenres_2_resid( ii );
		Size iibestrot = rot_sets_->rotid_on_moltenresidue( rotamer_assignment.first[ ii ] );

		return_pose->replace_residue(
			iiresid,
			*rot_sets_->rotamer_set_for_moltenresidue( ii )->rotamer( iibestrot ),
			false
		);
	}
	return return_pose;
}

void PackDaemon::assign_last_rotamers_to_pose( Pose & pose ) const
{
	for ( Size ii = 1; ii <= rot_sets_->nmoltenres(); ++ii ) {
		Size iiresid = rot_sets_->moltenres_2_resid( ii );
		Size iibestrot = rot_sets_->rotid_on_moltenresidue( last_assignment_.first[ ii ] );

		pose.replace_residue(
			iiresid,
			*rot_sets_->rotamer_set_for_moltenresidue( ii )->rotamer( iibestrot ),
			false
		);
	}

}


void PackDaemon::print_entity_history() const
{
	TR << "Entity History:\n";
	for ( EntityToRotamerHash::const_iterator hashiter = prev_state_hash_.begin(),
			hashiter_end = prev_state_hash_.end(); hashiter != hashiter_end; ++hashiter ) {
		TR << "Stored state:";
		for ( Size ii = 1; ii <= hashiter->first.size(); ++ii ) {
			TR << " " << hashiter->first[ ii ]->to_string();
		}
		TR << " fitness: " << hashiter->second.second << " rotamers: ";
		for ( Size ii = 1; ii <= hashiter->second.first.size(); ++ii ) {
			TR << " " << hashiter->second.first[ ii ];
		}
		TR << "\n";
	}
	TR << std::endl;
}

void PackDaemon::calculate_background_energies()
{
	using namespace core;
	using namespace core::scoring;
	using namespace core::graph;

	background_energies_ = 0;
	Energies const & energies( pose_->energies() );
	EnergyGraph const & energy_graph( energies.energy_graph() );
	for ( Size ii = 1; ii <= pose_->total_residue(); ++ii ) {
		if ( task_->being_packed( ii ) ) continue;
		background_energies_ += energies.onebody_energies( ii ).dot( score_function_->weights() );
		for ( Graph::EdgeListConstIter
				iru  = energy_graph.get_node(ii)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(ii)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge const & edge( static_cast< EnergyEdge const & > (**iru) );
			Size const jj = edge.get_second_node_ind();
			if ( task_->being_packed( jj ) ) continue;
			EnergyMap emap = edge.fill_energy_map();
			background_energies_ += score_function_->weights().dot( emap );
		}

	}
}

DaemonSet::DaemonSet() :
	num_entities_( 0 ),
	task_factory_( core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory ) ),
	include_background_energies_( true ),
	limit_dlig_mem_usage_( false ),
	dlig_nmeg_limit_( 0 ),
	ndaemons_( 0 ),
	n_npd_properties_( 0 )
{}

DaemonSet::~DaemonSet()
{}

/// @details  The entity resfile is slightly different from a regular resfile.
/// Its first line should consist of the number of residues that are entity; the
/// rest of the resfile lists those residues numbering from 1.
void DaemonSet::set_entity_resfile(
	std::string const & resfile
)
{
	utility::io::izstream infile( resfile );
	if ( ! infile ) {
		utility_exit_with_message( "Failed to open entity resfile: " + resfile );
	}
	set_entity_resfile( infile, resfile );
}


void DaemonSet::set_entity_resfile(
	std::istream & resfile,
	std::string const & resfile_name
)
{

	if ( ! daemons_.empty() ) {
		utility_exit_with_message( "Entity resfile must be set before any daemons may be added and all daemons must use the same correspondece file" );
	}

	create_entity_resfile_contents( resfile, resfile_name, entity_resfile_, entity_task_, num_entities_ );

}

void DaemonSet::set_score_function( ScoreFunction const & sfxn )
{
	score_function_ = sfxn.clone();
}

void DaemonSet::set_task_factory( core::pack::task::TaskFactoryOP factory )
{
	task_factory_ = factory;
}

void DaemonSet::set_include_background_energies( bool setting )
{
	include_background_energies_ = setting;
}

void DaemonSet::set_dlig_nmeg_limit( Size setting )
{
	if ( setting != 0 ) {
		limit_dlig_mem_usage_ = true;
		dlig_nmeg_limit_ = setting;
	} else {
		limit_dlig_mem_usage_ = false;
	}
}

void DaemonSet::add_npdpro_calculator_creator(
	NPDPropCalculatorCreatorOP creator
)
{
	npd_calculator_creators_[ creator->calculator_name() ] = creator;
}


/// @brief Each daemon is associated with an index representing its position in
/// some master list somewhere.  The DaemonSet is responsible for keeping this index.
void
DaemonSet::add_pack_daemon(
	Size daemon_index,
	std::string const & pdb_name,
	std::string const & correspondence_file_name,
	std::string const & secondary_resfile_name
)
{
	/// Read in the pdb
	Pose pose;
	core::import_pose::pose_from_pdb( pose, pdb_name );

	// Open the correspondence file
	utility::io::izstream correspondence_file( correspondence_file_name );
	if ( ! correspondence_file ) {
		utility_exit_with_message( "Could not open correspondence file named: " + correspondence_file_name );
	}

	// Open the resfile
	utility::io::izstream secondary_resfile( secondary_resfile_name );
	if ( ! secondary_resfile ) {
		utility_exit_with_message( "Could not open secondary resfile named: " + secondary_resfile_name );
	}
	add_pack_daemon( daemon_index, pdb_name, pose, correspondence_file_name, correspondence_file, secondary_resfile_name, secondary_resfile );
}

void
DaemonSet::add_pack_daemon(
	Size daemon_index,
	std::string const & pose_file_name,
	Pose const & pose,
	std::string const & correspondence_file_filename,
	std::istream & correspondence_file,
	std::string const & secondary_refile_file_filename,
	std::istream & secondary_resfile
)
{
	TR << "Adding daemon: " << daemon_index << " " << pose_file_name << " " << correspondence_file_filename << " " << secondary_refile_file_filename << std::endl;

	using namespace core;
	using namespace core::pack::task;

	if ( ! score_function_ ) {
		utility_exit_with_message( "ScoreFunction must be set before DaemonSet::add_pack_daemon may be called" );
	}
	if ( num_entities_ == 0 ) {
		utility_exit_with_message( "Entity resfile must be set before DaemonSet::add_pack_daemon may be called" );
	}

	/// Read the correspondence file
	EntityCorrespondenceOP ec( new EntityCorrespondence );
	ec->set_pose( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( pose ) ) ));
	ec->set_num_entities( num_entities_ );
	ec->initialize_from_correspondence_file( correspondence_file );

	/// Setup the task
	PackerTaskOP task = task_factory_->create_task_and_apply_taskoperations( pose );
	ResfileContentsOP secondary_resfile_contents( new ResfileContents( pose, secondary_resfile ) );

	initialize_task_from_entity_resfile_and_secondary_resfile( pose, ec, *entity_resfile_, *secondary_resfile_contents, task );
	task->initialize_from_command_line();

	/// At this point, the packer task is initialized and ready to be handed to the
	/// PackDaemon.
	///
	/// Take a moment and make sure that the task is in fact compatible with
	/// the entity resfile: Make sure the amino-acid contents of the residues
	/// with non-zero entity-id have the same allowable residues as the
	/// entity packer task.  If there is a discrepancy, then the program should
	/// be halted now.
	std::string error_message;
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		Size ii_entity_id = ec->entity_for_residue( ii );
		if ( ii_entity_id == 0 ) {
			continue;
		}

		core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter
			final_task_allowed = task->residue_task( ii ).allowed_residue_types_begin(),
			entity_task_allowed = entity_task_->residue_task( ii_entity_id ).allowed_residue_types_begin(),
			final_task_allowed_end = task->residue_task( ii ).allowed_residue_types_end(),
			entity_task_allowed_end = entity_task_->residue_task( ii_entity_id ).allowed_residue_types_end();

		utility::vector1< bool > final_aas( core::chemical::num_canonical_aas );
		utility::vector1< bool > entity_aas(  core::chemical::num_canonical_aas );

		/// NO SUPPORT YET FOR NONCANONICAL AMINO ACIDS

		bool bad( false );
		while ( final_task_allowed != final_task_allowed_end ) {
			if ( (*final_task_allowed)->aa() <= core::chemical::num_canonical_aas ) {
				final_aas[ (*final_task_allowed)->aa() ] = true;
			} else {
				bad = true;
				break;
			}
			++final_task_allowed;
		}
		while ( entity_task_allowed != entity_task_allowed_end ) {
			if ( (*entity_task_allowed)->aa() <= core::chemical::num_canonical_aas ) {
				entity_aas[ (*entity_task_allowed)->aa() ] = true;
			} else {
				bad = true;
				break;
			}
			++entity_task_allowed;
		}

		if ( ! bad ) {
			for ( Size ii = 1; ii <= core::chemical::num_canonical_aas; ++ii ) {
				if ( entity_aas[ ii ] != final_aas[ ii ] ) {
					bad = true;
					break;
				}
			}
		}
		if ( bad ) {
			std::string ii_error;
			ii_error = "Discrepancy between allowed residue types in the entity resfile and "
				"the final resfile for residue " + utility::to_string( ii ) + " in structure " +
				pose_file_name + " which corresponds to entity " + utility::to_string( ii_entity_id ) + "\n";

			ii_error += "Original residue type: " + pose.residue_type(ii).name() + "\n";

			ii_error += "Final packer task allows residue types\n";
			final_task_allowed = task->residue_task( ii ).allowed_residue_types_begin();
			final_task_allowed_end = task->residue_task( ii ).allowed_residue_types_end();
			while ( final_task_allowed != final_task_allowed_end ) {
				ii_error += utility::to_string( (*final_task_allowed)->aa()) + " (" + (*final_task_allowed)->name() + ")\n";
				++final_task_allowed;
			}

			ii_error += "Shared packer task allows residue types\n";
			entity_task_allowed = entity_task_->residue_task( ii_entity_id ).allowed_residue_types_begin();
			entity_task_allowed_end = entity_task_->residue_task( ii_entity_id ).allowed_residue_types_end();
			while ( entity_task_allowed != entity_task_allowed_end ) {
				ii_error += utility::to_string((*entity_task_allowed)->aa()) + "\n";
				++entity_task_allowed;
			}
			error_message += ii_error;
		}
	}
	if ( error_message != "" ) {
		//std::cerr << error_message << std::endl;
		throw utility::excn::EXCN_Msg_Exception( error_message );
	}

	// Create the new daemon
	PackDaemonOP daemon( new PackDaemon );
	daemon->set_pose_and_task( pose, *task );
	daemon->set_score_function( *score_function_ );
	daemon->set_entity_correspondence( *ec );
	daemon->set_include_background_energies( include_background_energies_ );
	if ( limit_dlig_mem_usage_ ) {
		daemon->set_dlig_nmeg_limit( dlig_nmeg_limit_ );
	}

	daemons_.push_back( std::make_pair( daemon_index, daemon ));
	daemon_poses_.push_back( core::pose::PoseOP( new core::pose::Pose( pose ) ));
	daemon_tasks_.push_back( task );
	npd_calcs_for_poses_.resize( daemon_tasks_.size() );
	++ndaemons_;
}

void
DaemonSet::add_npd_property_calculator_for_state(
	Size daemon_index,
	std::string const & npd_property,
	Size npd_index
)
{
	if ( npd_calculator_creators_.find( npd_property ) == npd_calculator_creators_.end() ) {
		std::string msg;
		if ( npd_calculator_creators_.empty() ) {
			msg = "Requested non-pairwise-decomposable (NPD) property '" + npd_property +"' but there are"
				" no NPD Property Calculator Creators that have been registered with the DaemonSet";
		} else {
			msg = "Requested non-pairwise-decomposable (NPD) property '" + npd_property +"' but no such NPD"
				" calculator creator has been registered with the DaemonSet.  Available calculators include\n";
			for ( std::map< std::string, NPDPropCalculatorCreatorOP >::const_iterator
					npditer = npd_calculator_creators_.begin(),
					npditer_end = npd_calculator_creators_.end();
					npditer != npditer_end; ++npditer ) {
				msg += "   " + npditer->first + "\n";
			}
		}
		throw utility::excn::EXCN_Msg_Exception( msg );
	}

	Size which_daemon = 0;
	for ( Size ii = 1; ii <= daemons_.size(); ++ii ) {
		if ( daemons_[ ii ].first == daemon_index ) {
			which_daemon = ii;
			break;
		}
	}
	if ( which_daemon == 0 ) {
		throw utility::excn::EXCN_Msg_Exception( "Internal error: could not locate requested daemon " +
			utility::to_string( daemon_index ) + " while trying to add NPD Property Calculator for that daemon." );
	}

	/// if we've made it this far, create, initialize and add the creator.
	NPDPropCalculatorOP calculator = npd_calculator_creators_[ npd_property ]->new_calculator();
	calculator->setup( *daemon_poses_[ which_daemon ], *daemon_tasks_[ which_daemon ] );
	npd_calcs_for_poses_[ which_daemon ].push_back( std::make_pair( npd_index, calculator ));
	++n_npd_properties_;
}


void DaemonSet::setup_daemons()
{
	for ( DaemonListIter iter = daemons_.begin(), iter_end = daemons_.end(); iter != iter_end; ++iter ) {
		iter->second->setup();
	}
}


core::Size DaemonSet::ndaemons() const { return ndaemons_; }

core::Size DaemonSet::n_npd_properties() const { return n_npd_properties_; }

DaemonSet::StateEsAndNPDs
DaemonSet::compute_energy_for_assignment( Entity const & entity )
{
#ifdef APL_MEASURE_MSD_LOAD_BALANCE
	std::clock_t starttime = clock();
#endif


	SizeRealPairs daemon_scores;
	for ( DaemonListIter iter = daemons_.begin(), iter_end = daemons_.end();
			iter != iter_end; ++iter ) {
		core::Real energy = iter->second->compute_energy_for_assignment( entity );
		daemon_scores.push_back( std::make_pair( iter->first, energy ) );
		//std::cout << "Computed energy " << energy << " for state " << iter->first << std::endl;
	}
#ifdef APL_MEASURE_MSD_LOAD_BALANCE
	std::clock_t midtime = clock();
#endif
	SizeRealPairs npd_properties;
	if ( n_npd_properties_ != 0 ) {
		for ( Size ii = 1; ii <= npd_calcs_for_poses_.size(); ++ii ) {
			if ( npd_calcs_for_poses_[ ii ].empty() ) continue;
			daemons_[ ii ].second->assign_last_rotamers_to_pose( *daemon_poses_[ ii ] );
			for ( std::list< NPDIndAndCalc >::const_iterator
					npditer = npd_calcs_for_poses_[ ii ].begin(),
					npditer_end = npd_calcs_for_poses_[ ii ].end();
					npditer != npditer_end; ++npditer ) {
				core::Real const npd_property = npditer->second->calculate( *daemon_poses_[ ii ] );
				npd_properties.push_back( std::make_pair( npditer->first, npd_property ));
				//std::cout << "Computed non-pairwise decomposable property " << npd_property << " for npdindex " << npditer->first << std::endl;
			}
		}
	}
#ifdef APL_MEASURE_MSD_LOAD_BALANCE
	std::clock_t stoptime = clock();
	packing_runtime_ = ((double) midtime  - starttime ) / CLOCKS_PER_SEC;
	npd_runtime_     = ((double) stoptime - midtime   ) / CLOCKS_PER_SEC;
#endif
	return std::make_pair( daemon_scores, npd_properties );
}

DaemonSet::ConstDaemonList
DaemonSet::daemons() const
{
	ConstDaemonList return_daemons;
	for ( DaemonListIter iter = daemons_.begin(), iter_end = daemons_.end();
			iter != iter_end; ++iter ) {
		std::pair< core::Size, PackDaemonCOP > new_element;
		new_element.first = iter->first;
		new_element.second = iter->second;
		return_daemons.push_back( new_element );
	}
	return return_daemons;
}

void DaemonSet::mark_last_entity_as_important()
{
	for ( DaemonListIter iter = daemons_.begin(), iter_end = daemons_.end(); iter != iter_end; ++iter ) {
		iter->second->mark_last_entity_as_important();
	}
}

void DaemonSet::mark_entity_as_unimportant( Entity const & ent )
{
	for ( DaemonListIter iter = daemons_.begin(), iter_end = daemons_.end(); iter != iter_end; ++iter ) {
		iter->second->mark_entity_as_unimportant( ent );
	}

}

std::list< std::pair< core::Size, core::pose::PoseOP > >
DaemonSet::retrieve_relevant_poses_for_entity(
	Entity const & ent,
	DaemonIndices const & daemon_indices
) const
{
	std::list< std::pair< Size, PoseOP > > return_list;
	for ( DaemonListIter iter = daemons_.begin(), iter_end = daemons_.end(); iter != iter_end; ++iter ) {
		bool generate_pose_for_daemon( false );
		for ( DaemonIndices::const_iterator
				index_iter = daemon_indices.begin(),
				index_iter_end = daemon_indices.end();
				index_iter != index_iter_end; ++index_iter ) {
			if ( iter->first == *index_iter ) {
				generate_pose_for_daemon = true;
			}
		}
		if ( ! generate_pose_for_daemon ) continue;
		return_list.push_back( std::make_pair( iter->first, iter->second->recreate_pose_for_entity( ent )) );
	}
	return return_list;
}

core::pack::task::PackerTaskOP       DaemonSet::entity_task() const { return entity_task_; }
core::pack::task::ResfileContentsCOP DaemonSet::entity_resfile() const { return entity_resfile_; }

void DaemonSet::activate_daemon_mode()
{
#ifdef USEMPI
	while ( true ) {
		/// Wait for a signal from node 0, and then process that request
		int signal = utility::receive_integer_from_node( 0 );
		DaemonSetMessage message( static_cast< DaemonSetMessage > (signal) );
		bool leave_main_while_loop( false );
		switch ( message ) {
			case error_message: {
				graceful_exit();
				break;
			}
			case success_message: {
				// noop?  Why are we recieving a success signal?
				break;
			}
			case add_daemon: {
				process_add_daemon_message();
				break;
			}
			case evaluate_entity: {
				process_state_energy_evaluations_for_entity();
				break;
			}
			case keep_rotamer_assignment_for_last_entity: {
				mark_last_entity_as_important();
				break;
			}
			case discard_old_entity: {
				process_discard_entity_message();
				break;
			}
			case geneate_pose_from_old_state: {
				process_pose_request_for_entity();
				break;
			}
			case spin_down: {
				leave_main_while_loop = true;
				break;
			}
			default: {
				utility::send_integer_to_node( 0, error_message );
				graceful_exit();
			}
		}
		if ( leave_main_while_loop ) break;
	}
#else
	utility_exit_with_message( "MPI-related function requested of class DaemonSet in a non-MPI build" );
#endif
}


void DaemonSet::process_add_daemon_message()
{
#ifdef USEMPI
	int n_daemons = utility::receive_integer_from_node( 0 );
	utility::vector1< int > daemon_indices( n_daemons );

	utility::vector1< std::string > pdb_names( n_daemons );
	utility::vector1< std::string > pdbs( n_daemons );

	utility::vector1< std::string > correspondence_file_names( n_daemons );
	utility::vector1< std::string > correspondence_files( n_daemons );

	utility::vector1< std::string > secondary_resfile_names( n_daemons );
	utility::vector1< std::string > secondary_resfiles( n_daemons );

	utility::vector1< std::list< std::pair< Size, std::string > > > npd_properties( n_daemons );

	for ( int ii = 1; ii <= n_daemons; ++ii ) {
		daemon_indices[ ii ]            = utility::receive_integer_from_node( 0 );
		pdb_names[ ii ]                 = utility::receive_string_from_node( 0 );
		pdbs[ ii ]                      = utility::receive_string_from_node( 0 );
		correspondence_file_names[ ii ] = utility::receive_string_from_node( 0 );
		correspondence_files[ ii ]      = utility::receive_string_from_node( 0 );
		secondary_resfile_names[ ii ]   = utility::receive_string_from_node( 0 );
		secondary_resfiles[ ii ]        = utility::receive_string_from_node( 0 );
		int n_npd_properties_for_state  = utility::receive_integer_from_node( 0 );
		for ( int jj = 1; jj <= n_npd_properties_for_state; ++jj ) {
			Size npd_property_id         = utility::receive_integer_from_node( 0 );
			std::string property         = utility::receive_string_from_node( 0 );
			npd_properties[ ii ].push_back( std::make_pair( npd_property_id, property ) );
		}
	}

	try {
		for ( int ii = 1; ii <= n_daemons; ++ii ) {
			Pose pose;
			core::import_pose::pose_from_pdbstring( pose, pdbs[ ii ], pdb_names[ ii ] );
			std::istringstream correspondence_file( correspondence_files[ ii ] );
			std::istringstream secondary_resfile( secondary_resfiles[ ii ] );
			add_pack_daemon(
				daemon_indices[ ii ],
				pdb_names[ ii ], pose,
				correspondence_file_names[ ii ], correspondence_file,
				secondary_resfile_names[ ii ], secondary_resfile );
			for ( std::list< std::pair< Size, std::string > >::const_iterator
					npditer     = npd_properties[ ii ].begin(),
					npditer_end = npd_properties[ ii ].end();
					npditer != npditer_end; ++npditer ) {
				add_npd_property_calculator_for_state( daemon_indices[ ii ], npditer->second, npditer->first );
			}

		}
	} catch ( utility::excn::EXCN_Base & excn )  {
		// send the error message to node 0 and exit gracefully.
		std::string message( excn.msg() );
		utility::send_integer_to_node( 0, error_message );
		utility::send_string_to_node( 0, message );
		graceful_exit();
	}

	/// If we got here, then we added all daemons successfully.  Let node 0 know that.
	/// It's there waiting to hear from us.
	utility::send_integer_to_node( 0, success_message );

	/// Now, ask, did all the other nodes successfully get added?
	int others_successful = utility::receive_integer_from_node( 0 );
	if ( others_successful == success_message ) {
		// Green light: go and precompute all rotamer pair energies
		setup_daemons();
	} else {
		graceful_exit();
	}
#endif
}


void DaemonSet::process_state_energy_evaluations_for_entity()
{
#ifdef USEMPI

#ifdef APL_MEASURE_MSD_LOAD_BALANCE
	std::clock_t starttime = clock();
#endif

	EntityOP entity = recieve_entity();
	StateEsAndNPDs entity_energies =
		compute_energy_for_assignment( *entity );

#ifdef APL_MEASURE_MSD_LOAD_BALANCE
	std::clock_t stoptime = clock();
#endif

	assert( ndaemons() == entity_energies.first.size() );
	assert( n_npd_properties_ == entity_energies.second.size() );

	int n_daemons( ndaemons() );
	utility::send_integer_to_node( 0, n_daemons );
	for ( std::list< std::pair< core::Size, core::Real > >::const_iterator
			iter = entity_energies.first.begin(),
			iter_end = entity_energies.first.end();
			iter != iter_end; ++iter ) {
		utility::send_integer_to_node( 0, iter->first );
		utility::send_double_to_node( 0, iter->second );
	}
	utility::send_integer_to_node( 0, n_npd_properties_ );
	for ( std::list< std::pair< core::Size, core::Real > >::const_iterator
			iter = entity_energies.second.begin(),
			iter_end = entity_energies.second.end();
			iter != iter_end; ++iter ) {
		utility::send_integer_to_node( 0, iter->first );
		utility::send_double_to_node( 0, iter->second );
	}

#ifdef APL_MEASURE_MSD_LOAD_BALANCE
	core::Real evaltime = ((double) stoptime - starttime) / CLOCKS_PER_SEC;
	//TR << "Node " << utility::mpi_rank() << " " << evaltime << " seconds to evaluate" << std::endl;
	utility::send_double_to_node( 0, evaltime );
	utility::send_double_to_node( 0, packing_runtime_ );
	utility::send_double_to_node( 0, npd_runtime_ );
#endif

	//TR << "Finished sending state energies" << std::endl;


#endif
}


void DaemonSet::process_discard_entity_message()
{
#ifdef USEMPI
	EntityOP entity = recieve_entity();
	mark_entity_as_unimportant( *entity );
#endif
}


void DaemonSet::process_pose_request_for_entity()
{
#ifdef USEMPI
	EntityOP entity = recieve_entity();
	DaemonIndices indices = recieve_daemon_inds_requiring_pose_creation();

	std::list< std::pair< Size, PoseOP > > poses = retrieve_relevant_poses_for_entity( *entity, indices );
	utility::send_integer_to_node( 0, poses.size() );
	for ( std::list< std::pair< Size, PoseOP > >::const_iterator
			iter = poses.begin(), iter_end = poses.end(); iter != iter_end; ++iter ) {
		std::ostringstream oss;
		core::io::pdb::dump_pdb( *( iter->second ), oss );
		utility::send_integer_to_node( 0, iter->first );
		utility::send_string_to_node(  0, oss.str() );
	}
#endif
}


DaemonSet::EntityOP
DaemonSet::recieve_entity() const
{
#ifdef USEMPI
	std::string entity_string = utility::receive_string_from_node( 0 );
	EntityOP ent( new Entity( entity_string ) );
	return ent;
#endif
	return 0;
}


DaemonSet::DaemonIndices
DaemonSet::recieve_daemon_inds_requiring_pose_creation() const
{
	DaemonIndices indices;
#ifdef USEMPI
	utility::vector1< int > int_indices = utility::receive_integers_from_node( 0 );
	indices = int_indices; // cast from unsigned to signed integers
#endif
	return indices;
}

void
DaemonSet::graceful_exit() const
{
#ifdef USEMPI
	MPI_Finalize();
#endif
	utility_exit();
}

NPDPropCalculator::NPDPropCalculator() {}
NPDPropCalculator::~NPDPropCalculator() {}

void
NPDPropCalculator::setup(
	core::pose::Pose const &,
	core::pack::task::PackerTask const &
)
{}


NPDPropCalculatorCreator::NPDPropCalculatorCreator() {}
NPDPropCalculatorCreator::~NPDPropCalculatorCreator() {}

QuickRepacker::QuickRepacker(
	PoseOP                     pose,
	PackerTaskOP               task,
	FixedBBInteractionGraphOP  ig,
	RotamerSetsOP              rot_sets
) :
	pose_( pose ),
	task_( task ),
	ig_( ig ),
	rot_sets_( rot_sets )
{}

QuickRepacker::~QuickRepacker() {}

QuickRepacker::PoseOP                    QuickRepacker::pose() { return pose_; }
QuickRepacker::PackerTaskOP              QuickRepacker::task() { return task_; }
QuickRepacker::FixedBBInteractionGraphOP QuickRepacker::ig()   { return ig_; }
QuickRepacker::RotamerSetsOP             QuickRepacker::rot_sets() { return rot_sets_; }
void QuickRepacker::task( PackerTaskOP setting ) { task_ = setting; }


BasicSimAnnealerRepacker::BasicSimAnnealerRepacker(
	PoseOP                 pose,
	PackerTaskOP           task,
	FixedBBInteractionGraphOP ig,
	RotamerSetsOP          rotsets
) :
	parent( pose, task, ig, rotsets )
{
	ig->prepare_for_simulated_annealing();
}

BasicSimAnnealerRepacker::~BasicSimAnnealerRepacker() {}

BasicSimAnnealerRepacker::RotamerAssignmentAndEnergy
BasicSimAnnealerRepacker::repack( utility::vector0< int > const & rot_to_pack )
{
	ObjexxFCL::FArray1D< int > rotamer_assignment( pose()->total_residue() );
	core::PackerEnergy   rotamer_energy;

	core::pack::pack_rotamers_run(
		*pose(), task(), rot_sets(), ig(), rot_to_pack,
		rotamer_assignment, rotamer_energy );

	RotamerAssignmentAndEnergy result;
	result.first.resize( rot_sets()->nmoltenres() );

	for ( Size ii = 1; ii <= rot_sets()->nmoltenres(); ++ii ) {
		result.first[ ii ] = rotamer_assignment( rot_sets()->moltenres_2_resid( ii ) );
	}
	result.second = rotamer_energy;
	return result;
}


RotamerSubsetRepacker::RotamerSubsetRepacker(
	PoseOP                 pose,
	PackerTaskOP           task,
	FixedBBInteractionGraphOP ig,
	RotamerSetsOP          rotsets
) :
	parent( pose, task, ig, rotsets )
{}

RotamerSubsetRepacker::~RotamerSubsetRepacker() {}

RotamerSubsetRepacker::RotamerSubsetsOP
RotamerSubsetRepacker::create_rotamer_subsets_from_rot_to_pack(
	utility::vector0< int > const & rot_to_pack
)
{
	return RotamerSubsetRepacker::RotamerSubsetsOP( new core::pack::rotamer_set::RotamerSubsets( *rot_sets(), rot_to_pack ) );
}


/// DENSE IG REPACKER

DenseIGRepacker::DenseIGRepacker(
	PoseOP                     pose,
	PackerTaskOP               task,
	FixedBBInteractionGraphOP  ig,
	RotamerSetsOP              rotsets
) :
	parent( pose, task, ig, rotsets )
{
	ig->prepare_for_simulated_annealing();
}

DenseIGRepacker::~DenseIGRepacker() {}

void DenseIGRepacker::set_MCA() {
	PackerTaskOP copy_task = task()->clone();
	copy_task->or_multi_cool_annealer( true ); // cannot be undone
	task( copy_task );
}

DenseIGRepacker::RotamerAssignmentAndEnergy
DenseIGRepacker::repack( utility::vector0< int > const & rot_to_pack )
{
	using namespace core::pack::rotamer_set;
	using namespace core::pack::interaction_graph;
	utility::vector0< int > local_rot_to_pack( rot_to_pack );
	/// SHOULD be in sorted order already, but make sure!
	std::sort( local_rot_to_pack.begin(), local_rot_to_pack.end() );

	utility::subset_mapping rotamer_subset_map( rot_sets()->nrotamers() );
	rotamer_subset_map.reserve_destination_size( local_rot_to_pack.size() );
	for ( Size ii = 0; ii < local_rot_to_pack.size(); ++ii ) {
		rotamer_subset_map.set_next_correspondence( local_rot_to_pack[ ii ] );
	}

	RotamerSubsetsOP rsubset = create_rotamer_subsets_from_rot_to_pack( local_rot_to_pack );
	DensePDInteractionGraphOP dense_ig =
		create_dense_pdig_from_rot_to_pack( local_rot_to_pack, rsubset );

	//core::Size nmoltenres = rot_sets()->nmoltenres();
	ObjexxFCL::FArray1D< int > rotamer_assignment( pose()->total_residue() );
	core::PackerEnergy   rotamer_energy;
	utility::vector0< int > all_rots( local_rot_to_pack.size() );
	for ( Size ii = 0; ii < all_rots.size(); ++ii ) all_rots[ ii ] = ii+1;

	core::pack::pack_rotamers_run(
		*pose(), task(), rsubset, dense_ig, all_rots,
		rotamer_assignment, rotamer_energy );

	RotamerAssignmentAndEnergy result;
	result.first.resize( rot_sets()->nmoltenres() );

	for ( Size ii = 1; ii <= rot_sets()->nmoltenres(); ++ii ) {
		result.first[ ii ] = rotamer_subset_map.d2s( rotamer_assignment( rot_sets()->moltenres_2_resid( ii ) ));
	}

	/*TR << "Assignment:";
	for ( Size ii = 1; ii <= rot_sets()->nmoltenres(); ++ii ) {
	TR << " " << rotamer_assignment( rot_sets()->moltenres_2_resid( ii ));
	}
	TR << std::endl;*/
	result.second = rotamer_energy;
	return result;
}

DenseIGRepacker::DensePDInteractionGraphOP
DenseIGRepacker::create_dense_pdig_from_rot_to_pack(
	utility::vector0< int > const & rot_to_pack,
	RotamerSubsetsOP rot_subsets
)
{
	using namespace core::pack::interaction_graph;

	core::Size nmoltenres = rot_sets()->nmoltenres();
	DensePDInteractionGraphOP dense_ig( new DensePDInteractionGraph( nmoltenres ) );
	dense_ig->initialize( *rot_subsets );

	/// If the rotamers are in a contiguous block and all have the same "aa" according to
	/// the PDInteractionGraph, then we can rapidly transfer the energies from the PDIG
	/// to the new dense interaction graph using the new get_aa_submatrix_energies_for_edge
	/// function; otherwise, we have to fall back on the old
	utility::vector1< int > aaind_for_moltres( nmoltenres, -1 );

	for ( Size ii = 1; ii <= rot_to_pack.size(); ++ii ) {
		Size const ii_moltenres = rot_subsets->moltenres_for_rotamer( ii );
		Size const ii_local_id  = rot_subsets->rotid_on_moltenresidue( ii );
		Size const ii_old_rotid = rot_sets()->rotid_on_moltenresidue(rot_to_pack[ ii - 1 ]);

		if ( aaind_for_moltres[ ii_moltenres ] == -1 ) {
			aaind_for_moltres[ ii_moltenres ] = ig()->aatype_for_node_state( ii_moltenres, ii_old_rotid );
		} else if ( aaind_for_moltres[ ii_moltenres ] != 0 ) {
			if ( aaind_for_moltres[ ii_moltenres ] != ig()->aatype_for_node_state( ii_moltenres, ii_old_rotid ) ) {
				aaind_for_moltres[ ii_moltenres ] = 0;
			}
		}

		dense_ig->add_to_nodes_one_body_energy(
			ii_moltenres,
			ii_local_id,
			ig()->get_one_body_energy_for_node_state( ii_moltenres, ii_old_rotid ));
		//std::cout << "Setting rotamer 1-body-energy " << ii_moltenres << " rot# " << ii_local_id << " (";
		//std::cout << rot_to_pack[ ii - 1 ] << ") = ";
		//std::cout << ig()->get_one_body_energy_for_node_state( ii_moltenres, ii_old_rotid) << std::endl;
	}

	for ( Size ii = 1; ii <= nmoltenres; ++ii ) {
		for ( Size jj = ii+1; jj <= nmoltenres; ++jj ) {
			if ( ig()->get_edge_exists( ii, jj ) ) {
				if ( aaind_for_moltres[ ii ] > 0 && aaind_for_moltres[ jj ] > 0 ) {
					if ( ! ig()->get_sparse_aa_info_for_edge( ii, jj,aaind_for_moltres[ ii ], aaind_for_moltres[ jj ] ) ) {
						/// NO non-zero energies for this edge
						continue;
					}
				}
				dense_ig->add_edge( ii, jj );
				bool fast_edge_energy_transfer_successfull = false;
				if ( aaind_for_moltres[ ii ] > 0 && aaind_for_moltres[ jj ] > 0 ) {
					ObjexxFCL::FArray2D< core::PackerEnergy > edge_energies =
						ig()->get_aa_submatrix_energies_for_edge( ii, jj, aaind_for_moltres[ ii ], aaind_for_moltres[ jj ]  );
					if ( edge_energies.size1() == rot_subsets->nrotamers_for_moltenres( jj ) &&
							edge_energies.size2() == rot_subsets->nrotamers_for_moltenres( ii ) ) {
						// Safe to transfer the edge energy table to the DenseIG
						fast_edge_energy_transfer_successfull = true;
						dense_ig->swap_edge_energies( ii, jj, edge_energies );
					} else {
						TR << "Surprisingly, edge " << ii << " " << jj << " did not have the right sized"
							<< " edge energy table from the PDIG:" << rot_subsets->nrotamers_for_moltenres( ii )
							<< " " << edge_energies.size2()  << " " << rot_subsets->nrotamers_for_moltenres( jj )
							<< " " << edge_energies.size1() << std::endl;
					}
				}

				if ( ! fast_edge_energy_transfer_successfull ) {
					for ( Size kk = 1, kk_end = rot_subsets->nrotamers_for_moltenres( ii ); kk <= kk_end; ++kk ) {
						Size const kk_old = rot_sets()->rotid_on_moltenresidue(
							rot_to_pack[ rot_subsets->moltenres_rotid_2_rotid( ii, kk ) - 1 ] );
						for ( Size ll = 1, ll_end = rot_subsets->nrotamers_for_moltenres( jj ); ll <= ll_end; ++ll ) {
							Size const ll_old = rot_sets()->rotid_on_moltenresidue(
								rot_to_pack[ rot_subsets->moltenres_rotid_2_rotid( jj, ll ) - 1 ] );
							dense_ig->set_two_body_energy_for_edge( ii, jj, kk, ll,
								ig()->get_two_body_energy_for_edge( ii, jj, kk_old, ll_old ) );
							//std::cout << "Setting energy for edge " << ii << " " << jj;
							//std::cout << " rotamer pair " << kk << " (" << kk_old << ") and " << ll << " (" << ll_old;
							//std::cout << " ) to " << ig()->get_two_body_energy_for_edge( ii, jj, kk_old, ll_old ) << std::endl;
						}
					}
				}
			}
		}
	}
	return dense_ig;

}

//// DOUBLE DENSE IG REPACKER


DoubleDenseIGRepacker::DoubleDenseIGRepacker(
	PoseOP                                    pose,
	PackerTaskOP                              task,
	FixedBBInteractionGraphOP ig,
	RotamerSetsOP                             rotsets
) :
	parent( pose, task, ig, rotsets )
{
	ig->prepare_for_simulated_annealing();
}

DoubleDenseIGRepacker::~DoubleDenseIGRepacker() {}

DoubleDenseIGRepacker::RotamerAssignmentAndEnergy
DoubleDenseIGRepacker::repack( utility::vector0< int > const & rot_to_pack )
{
	using namespace core::pack::rotamer_set;
	using namespace core::pack::interaction_graph;
	utility::vector0< int > local_rot_to_pack( rot_to_pack );
	/// SHOULD be in sorted order already, but make sure!
	std::sort( local_rot_to_pack.begin(), local_rot_to_pack.end() );

	utility::subset_mapping rotamer_subset_map( rot_sets()->nrotamers() );
	rotamer_subset_map.reserve_destination_size( local_rot_to_pack.size() );
	for ( Size ii = 0; ii < local_rot_to_pack.size(); ++ii ) {
		rotamer_subset_map.set_next_correspondence( local_rot_to_pack[ ii ] );
	}

	RotamerSubsetsOP rsubset = create_rotamer_subsets_from_rot_to_pack( local_rot_to_pack );
	DoubleDensePDInteractionGraphOP dense_ig =
		create_dense_pdig_from_rot_to_pack( local_rot_to_pack, rsubset );

	//core::Size nmoltenres = rot_sets()->nmoltenres();
	ObjexxFCL::FArray1D< int > rotamer_assignment( pose()->total_residue() );
	core::PackerEnergy   rotamer_energy;
	utility::vector0< int > all_rots( local_rot_to_pack.size() );
	for ( Size ii = 0; ii < all_rots.size(); ++ii ) all_rots[ ii ] = ii+1;

	core::pack::pack_rotamers_run(
		*pose(), task(), rsubset, dense_ig, all_rots,
		rotamer_assignment, rotamer_energy );

	RotamerAssignmentAndEnergy result;
	result.first.resize( rot_sets()->nmoltenres() );

	for ( Size ii = 1; ii <= rot_sets()->nmoltenres(); ++ii ) {
		result.first[ ii ] = rotamer_subset_map.d2s( rotamer_assignment( rot_sets()->moltenres_2_resid( ii ) ));
	}

	/*TR << "Assignment:";
	for ( Size ii = 1; ii <= rot_sets()->nmoltenres(); ++ii ) {
	TR << " " << rotamer_assignment( rot_sets()->moltenres_2_resid( ii ));
	}
	TR << std::endl;*/
	result.second = rotamer_energy;
	return result;
}


DoubleDenseIGRepacker::DoubleDensePDInteractionGraphOP
DoubleDenseIGRepacker::create_dense_pdig_from_rot_to_pack(
	utility::vector0< int > const & rot_to_pack,
	RotamerSubsetsOP rot_subsets
)
{
	using namespace core::pack::interaction_graph;

	core::Size nmoltenres = rot_sets()->nmoltenres();
	DoubleDensePDInteractionGraphOP dense_ig( new DoubleDensePDInteractionGraph( nmoltenres ) );
	dense_ig->initialize( *rot_subsets );
	for ( Size ii = 1; ii <= rot_to_pack.size(); ++ii ) {
		Size const ii_moltenres = rot_subsets->moltenres_for_rotamer( ii );
		Size const ii_local_id  = rot_subsets->rotid_on_moltenresidue( ii );
		Size const ii_old_rotid = rot_sets()->rotid_on_moltenresidue(rot_to_pack[ ii - 1 ]);
		dense_ig->add_to_nodes_one_body_energy(
			ii_moltenres,
			ii_local_id,
			ig()->get_one_body_energy_for_node_state( ii_moltenres, ii_old_rotid ));
		//std::cout << "Setting rotamer 1-body-energy " << ii_moltenres << " rot# " << ii_local_id << " (";
		//std::cout << rot_to_pack[ ii - 1 ] << ") = ";
		//std::cout << ig()->get_one_body_energy_for_node_state( ii_moltenres, ii_old_rotid) << std::endl;
	}

	for ( Size ii = 1; ii <= nmoltenres; ++ii ) {
		for ( Size jj = ii+1; jj <= nmoltenres; ++jj ) {
			if ( ig()->get_edge_exists( ii, jj ) ) {
				dense_ig->add_edge( ii, jj );
				for ( Size kk = 1, kk_end = rot_subsets->nrotamers_for_moltenres( ii ); kk <= kk_end; ++kk ) {
					Size const kk_old = rot_sets()->rotid_on_moltenresidue(
						rot_to_pack[ rot_subsets->moltenres_rotid_2_rotid( ii, kk ) - 1 ] );
					for ( Size ll = 1, ll_end = rot_subsets->nrotamers_for_moltenres( jj ); ll <= ll_end; ++ll ) {
						Size const ll_old = rot_sets()->rotid_on_moltenresidue(
							rot_to_pack[ rot_subsets->moltenres_rotid_2_rotid( jj, ll ) - 1 ] );
						dense_ig->set_two_body_energy_for_edge( ii, jj, kk, ll,
							ig()->get_two_body_energy_for_edge( ii, jj, kk_old, ll_old ) );
						//std::cout << "Setting energy for edge " << ii << " " << jj;
						//std::cout << " rotamer pair " << kk << " (" << kk_old << ") and " << ll << " (" << ll_old;
						//std::cout << " ) to " << ig()->get_two_body_energy_for_edge( ii, jj, kk_old, ll_old ) << std::endl;
					}
				}
			}
		}
	}
	return dense_ig;

}

FASTER_IG_Repacker::FASTER_IG_Repacker(
	PoseOP                 pose,
	PackerTaskOP           task,
	FixedBBInteractionGraphOP ig,
	RotamerSetsOP          rotsets
) :
	parent( pose, task, ig, rotsets ),
	num_sa_( 5 ),
	sa_scale_( 0.05 ),
	ciBR_only_( false )
{
	ig->prepare_for_simulated_annealing();
}

FASTER_IG_Repacker::~FASTER_IG_Repacker() {}

FASTER_IG_Repacker::RotamerAssignmentAndEnergy
FASTER_IG_Repacker::repack( utility::vector0< int > const & rot_to_pack )
{
	using namespace core::pack::rotamer_set;
	using namespace core::pack::interaction_graph;
	utility::vector0< int > local_rot_to_pack( rot_to_pack );
	/// SHOULD be in sorted order already, but make sure!
	std::sort( local_rot_to_pack.begin(), local_rot_to_pack.end() );

	utility::subset_mapping rotamer_subset_map( rot_sets()->nrotamers() );
	rotamer_subset_map.reserve_destination_size( local_rot_to_pack.size() );
	for ( Size ii = 0; ii < local_rot_to_pack.size(); ++ii ) {
		rotamer_subset_map.set_next_correspondence( local_rot_to_pack[ ii ] );
	}

	RotamerSubsetsOP rsubset = create_rotamer_subsets_from_rot_to_pack( local_rot_to_pack );
	FASTERInteractionGraphOP faster_ig =
		create_faster_ig_from_rot_to_pack( local_rot_to_pack, rsubset );

	//core::Size nmoltenres = rot_sets()->nmoltenres();
	ObjexxFCL::FArray1D< int > rotamer_assignment( pose()->total_residue() );
	core::PackerEnergy   rotamer_energy;
	utility::vector0< int > all_rots( local_rot_to_pack.size() );
	for ( Size ii = 0; ii < all_rots.size(); ++ii ) all_rots[ ii ] = ii+1;

	{ // scope
		using namespace ObjexxFCL;
		using namespace core::pack;
		using namespace core::pack::annealer;

		bool start_with_current = false;
		FArray1D_int current_rot_index( pose()->total_residue(), 0 );
		bool calc_rot_freq = false;
		FArray1D< core::PackerEnergy > rot_freq( faster_ig->get_num_total_states(), 0.0 );

		FASTERAnnealer fa( rotamer_assignment, rotamer_energy, start_with_current,
			faster_ig, rsubset, current_rot_index, calc_rot_freq, rot_freq );

		fa.set_ciBR_only( ciBR_only_ );
		fa.set_num_sa_trajectories( num_sa_ );
		fa.set_sa_length_scale( sa_scale_ );

		fa.run();
	}


	RotamerAssignmentAndEnergy result;
	result.first.resize( rot_sets()->nmoltenres() );

	for ( Size ii = 1; ii <= rot_sets()->nmoltenres(); ++ii ) {
		result.first[ ii ] = rotamer_subset_map.d2s( rotamer_assignment( rot_sets()->moltenres_2_resid( ii ) ));
	}

	/*TR << "Assignment:";
	for ( Size ii = 1; ii <= rot_sets()->nmoltenres(); ++ii ) {
	TR << " " << rotamer_assignment( rot_sets()->moltenres_2_resid( ii ));
	}
	TR << std::endl;*/
	result.second = rotamer_energy;
	return result;

}

FASTER_IG_Repacker::FASTERInteractionGraphOP
FASTER_IG_Repacker::create_faster_ig_from_rot_to_pack(
	utility::vector0< int > const & rot_to_pack,
	RotamerSubsetsOP rot_subsets
)
{
	using namespace core::pack::interaction_graph;

	core::Size nmoltenres = rot_sets()->nmoltenres();
	FASTERInteractionGraphOP faster_ig( new FASTERInteractionGraph( nmoltenres ) );
	faster_ig->initialize( *rot_subsets );


	/// If the rotamers are in a contiguous block and all have the same "aa" according to
	/// the PDInteractionGraph, then we can rapidly transfer the energies from the PDIG
	/// to the new dense interaction graph using the new get_aa_submatrix_energies_for_edge
	/// function; otherwise, we have to fall back on the old
	utility::vector1< int > aaind_for_moltres( nmoltenres, -1 );

	for ( Size ii = 1; ii <= rot_to_pack.size(); ++ii ) {
		Size const ii_moltenres = rot_subsets->moltenres_for_rotamer( ii );
		Size const ii_local_id  = rot_subsets->rotid_on_moltenresidue( ii );
		Size const ii_old_rotid = rot_sets()->rotid_on_moltenresidue(rot_to_pack[ ii - 1 ]);

		if ( aaind_for_moltres[ ii_moltenres ] == -1 ) {
			aaind_for_moltres[ ii_moltenres ] = ig()->aatype_for_node_state( ii_moltenres, ii_old_rotid );
		} else if ( aaind_for_moltres[ ii_moltenres ] != 0 ) {
			if ( aaind_for_moltres[ ii_moltenres ] != ig()->aatype_for_node_state( ii_moltenres, ii_old_rotid ) ) {
				aaind_for_moltres[ ii_moltenres ] = 0;
			}
		}

		faster_ig->add_to_nodes_one_body_energy(
			ii_moltenres,
			ii_local_id,
			ig()->get_one_body_energy_for_node_state( ii_moltenres, ii_old_rotid ));
		//std::cout << "Setting rotamer 1-body-energy " << ii_moltenres << " rot# " << ii_local_id << " (";
		//std::cout << rot_to_pack[ ii - 1 ] << ") = ";
		//std::cout << ig()->get_one_body_energy_for_node_state( ii_moltenres, ii_old_rotid) << std::endl;
	}

	for ( Size ii = 1; ii <= nmoltenres; ++ii ) {
		for ( Size jj = ii+1; jj <= nmoltenres; ++jj ) {
			if ( ig()->get_edge_exists( ii, jj ) ) {
				if ( aaind_for_moltres[ ii ] > 0 && aaind_for_moltres[ jj ] > 0 ) {
					if ( ! ig()->get_sparse_aa_info_for_edge( ii, jj, aaind_for_moltres[ ii ], aaind_for_moltres[ jj ] ) ) {
						/// NO non-zero energies for this edge
						continue;
					}
				}
				faster_ig->add_edge( ii, jj );
				bool fast_edge_energy_transfer_successfull = false;
				if ( aaind_for_moltres[ ii ] > 0 && aaind_for_moltres[ jj ] > 0 ) {
					ObjexxFCL::FArray2D< core::PackerEnergy > edge_energies =
						ig()->get_aa_submatrix_energies_for_edge( ii, jj, aaind_for_moltres[ ii ], aaind_for_moltres[ jj ]  );
					if ( edge_energies.size1() == rot_subsets->nrotamers_for_moltenres( jj ) &&
							edge_energies.size2() == rot_subsets->nrotamers_for_moltenres( ii ) ) {
						// Safe to transfer the edge energy table to the DenseIG
						fast_edge_energy_transfer_successfull = true;
						faster_ig->swap_edge_energies( ii, jj, edge_energies );
					} else {
						TR << "Surprisingly, edge " << ii << " " << jj << " did not have the right sized"
							<< " edge energy table from the PDIG:" << rot_subsets->nrotamers_for_moltenres( ii )
							<< " " << edge_energies.size2()  << " " << rot_subsets->nrotamers_for_moltenres( jj )
							<< " " << edge_energies.size1() << std::endl;
					}
				}

				if ( ! fast_edge_energy_transfer_successfull ) {
					for ( Size kk = 1, kk_end = rot_subsets->nrotamers_for_moltenres( ii ); kk <= kk_end; ++kk ) {
						Size const kk_old = rot_sets()->rotid_on_moltenresidue(
							rot_to_pack[ rot_subsets->moltenres_rotid_2_rotid( ii, kk ) - 1 ] );
						for ( Size ll = 1, ll_end = rot_subsets->nrotamers_for_moltenres( jj ); ll <= ll_end; ++ll ) {
							Size const ll_old = rot_sets()->rotid_on_moltenresidue(
								rot_to_pack[ rot_subsets->moltenres_rotid_2_rotid( jj, ll ) - 1 ] );
							faster_ig->set_two_body_energy_for_edge( ii, jj, kk, ll,
								ig()->get_two_body_energy_for_edge( ii, jj, kk_old, ll_old ) );
							//std::cout << "Setting energy for edge " << ii << " " << jj;
							//std::cout << " rotamer pair " << kk << " (" << kk_old << ") and " << ll << " (" << ll_old;
							//std::cout << " ) to " << ig()->get_two_body_energy_for_edge( ii, jj, kk_old, ll_old ) << std::endl;
						}
					}
				}
			}
		}
	}
	return faster_ig;

}

void FASTER_IG_Repacker::set_num_sa( int setting ) { num_sa_ = setting; }
void FASTER_IG_Repacker::set_sa_scale( core::Real setting ) { sa_scale_ = setting; }
void FASTER_IG_Repacker::set_ciBR_only( bool setting ) { ciBR_only_ = setting; }


}
}
