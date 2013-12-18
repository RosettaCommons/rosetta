// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pack_daemon/PackDaemon.hh
/// @brief  declaration for class PackDaemon
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_pack_daemon_PackDaemon_hh
#define INCLUDED_protocols_pack_daemon_PackDaemon_hh

// Unit headers
#include <protocols/pack_daemon/PackDaemon.fwd.hh>

// Package headers
#include <protocols/pack_daemon/EntityCorrespondence.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/interaction_graph/DensePDInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/DoubleDensePDInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/FixedBBInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/FASTERInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/RotamerSubsets.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/ResfileReader.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/genetic_algorithm/Entity.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL headers
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>

// C++ headers
// AUTO-REMOVED #include <utility>
#include <list>

// Boost headers
#include <boost/unordered_map.hpp>

#include <utility/vector1.hh>

#ifdef WIN32
	#include <core/pack/task/PackerTask.hh>
#endif


namespace protocols {
namespace pack_daemon {

class PackDaemon : public utility::pointer::ReferenceCount {
public:
	typedef core::Real Real;
	typedef core::Size Size;

	typedef core::pose::Pose                Pose;
	typedef core::pose::PoseOP              PoseOP;
	typedef core::pose::PoseCOP             PoseCOP;
	typedef core::scoring::ScoreFunction    ScoreFunction;
	typedef core::scoring::ScoreFunctionOP  ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::pack::task::PackerTask    PackerTask;
	typedef core::pack::task::PackerTaskOP  PackerTaskOP;
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;

	typedef core::pack::interaction_graph::InteractionGraphBaseOP  InteractionGraphBaseOP;
	typedef core::pack::interaction_graph::InteractionGraphBaseCOP InteractionGraphBaseCOP;
	typedef core::pack::interaction_graph::FixedBBInteractionGraphOP  FixedBBInteractionGraphOP;
	typedef core::pack::interaction_graph::FixedBBInteractionGraphCOP FixedBBInteractionGraphCOP;
	typedef core::pack::rotamer_set::RotamerSets    RotamerSets;
	typedef core::pack::rotamer_set::RotamerSetsOP  RotamerSetsOP;
	typedef core::pack::rotamer_set::RotamerSetsCOP RotamerSetsCOP;

	typedef protocols::genetic_algorithm::Entity           Entity;
	typedef protocols::genetic_algorithm::EntityOP         EntityOP;
	typedef protocols::genetic_algorithm::EntityElement    EntityElement;
	typedef protocols::genetic_algorithm::EntityElements   EntityElements;
	typedef protocols::genetic_algorithm::EntityElementOP  EntityElementOP;
	typedef protocols::genetic_algorithm::Vec1Hash         Vec1Hash;

	typedef utility::vector1< int >              RotamerAssignment;
	typedef std::pair< RotamerAssignment, Real > RotamerAssignmentAndEnergy;
	typedef protocols::genetic_algorithm::EntityElementsEqual  EntElemEq;
	typedef boost::unordered_map< EntityElements, RotamerAssignmentAndEnergy, Vec1Hash, EntElemEq > EntityToRotamerHash;

public:
	PackDaemon();
	virtual ~PackDaemon();

	// Initialize the PackDaemon with the appropriate data before
	// calling setup().
	void set_pose_and_task( Pose const &, PackerTask const & );
	void set_score_function( ScoreFunction const & );
	void set_entity_correspondence( EntityCorrespondence const &  );
	void set_include_background_energies( bool setting );
	/// @brief restrict the amount of memory spent on storing Rotamer Pair Energies in the
	/// DoubleLazyInteractionGraph;
	void set_dlig_nmeg_limit( Size setting );

	void setup(); // must be called before compute_energy_for_assignment;

	/// @brief Repack the structure with the Entity
	/// This function proceeds in two steps: it creates a list of
	/// rotamer indices to be used during the repacking, and then
	/// it uses that list to repack the rotamers.  The first step
	/// is taken care of by the select_rotamer_subset method.
	Real compute_energy_for_assignment( Entity const & );

	utility::vector0< int >
	select_rotamer_subset( Entity const & ) const;

	void mark_last_entity_as_important();
	void mark_entity_as_unimportant( Entity const & );

	PoseCOP                 pose() const;
	ScoreFunctionCOP        score_function() const;
	PackerTaskCOP           task() const;
	EntityCorrespondenceCOP correspondence() const;
	FixedBBInteractionGraphCOP ig() const;
	RotamerSetsCOP          rot_sets() const;

	RotamerAssignmentAndEnergy const & best_assignment() const;
	RotamerAssignmentAndEnergy const & last_assignment() const;

	PoseOP recreate_pose_for_entity( Entity const & ) const;
	void assign_last_rotamers_to_pose( Pose & pose ) const;

	void print_entity_history() const;

private:
	void calculate_background_energies();

private:
	// Primary data -- set from outside
	PoseOP                 pose_;
	ScoreFunctionOP        score_function_;
	PackerTaskOP           task_;
	EntityCorrespondenceOP correspondence_;

	/// Should one-body energies for background residues and
	/// two body energies for pairs of background residues
	/// be included in the total energy for the state after repacking,
	/// along with the energies of those residues which have been repacked?
	bool                   include_background_energies_;
	Real                   background_energies_;

	// Derived data -- set internally
	bool                   setup_complete_;
	FixedBBInteractionGraphOP ig_;
	RotamerSetsOP          rot_sets_;

	QuickRepackerOP        repacker_;
	/*QuickRepackerOP        repacker2_;
	QuickRepackerOP        repacker3_;
	QuickRepackerOP        repacker4_;
	QuickRepackerOP        repacker5_;
	QuickRepackerOP        repacker6_;
	QuickRepackerOP        repacker7_;
	QuickRepackerOP        repacker8_;*/

	bool                       best_assignment_valid_;
	RotamerAssignmentAndEnergy best_assignment_;
	EntityOP                   best_entity_;
	RotamerAssignmentAndEnergy last_assignment_;
	EntityOP                   last_entity_;

	EntityToRotamerHash prev_state_hash_;

};

enum DaemonSetMessage {
	error_message = 0,
	success_message,
	add_daemon,
	evaluate_entity,
	keep_rotamer_assignment_for_last_entity,
	discard_old_entity,
	geneate_pose_from_old_state,
	spin_down,
	n_daemon_set_messages
};

class DaemonSet : public utility::pointer::ReferenceCount
{
public:
	typedef core::pose::Pose                       Pose;
	typedef core::pose::PoseOP                     PoseOP;
	typedef core::scoring::ScoreFunction           ScoreFunction;
	typedef core::scoring::ScoreFunctionOP         ScoreFunctionOP;
	typedef core::pack::task::PackerTask           PackerTask;
	typedef core::pack::task::PackerTaskOP         PackerTaskOP;
	typedef core::pack::task::ResfileContentsOP    ResfileContentsOP;
	typedef protocols::genetic_algorithm::Entity   Entity;
	typedef protocols::genetic_algorithm::EntityOP EntityOP;

	typedef utility::vector1< core::Size >                             DaemonIndices;
	typedef utility::vector1< std::pair< core::Size, PackDaemonCOP > > ConstDaemonList;
	typedef ConstDaemonList::const_iterator                            ConstDaemonListIter;
	typedef utility::vector1< std::pair< core::Size, PackDaemonOP > >  DaemonList;
	typedef DaemonList::const_iterator                                 DaemonListIter;

	typedef std::pair< core::Size, core::Real >          SizeRealPair;
	typedef std::list< SizeRealPair >                    SizeRealPairs;
	typedef std::pair< SizeRealPairs, SizeRealPairs >    StateEsAndNPDs;

	typedef std::pair< core::Size, NPDPropCalculatorOP > NPDIndAndCalc;
public:
	DaemonSet();
	virtual ~DaemonSet();

	void set_entity_resfile( std::string const & resfile );
	void set_entity_resfile( std::istream & resfile, std::string const & resfile_name );

	void set_score_function( ScoreFunction const & );
	void set_task_factory( core::pack::task::TaskFactoryOP factory );

	void set_include_background_energies( bool setting );

	/// @brief restrict the amount of memory spent on storing Rotamer Pair Energies in the
	/// DoubleLazyInteractionGraph;
	void set_dlig_nmeg_limit( Size setting );

	void add_npdpro_calculator_creator( NPDPropCalculatorCreatorOP );

	/// @brief Each daemon is associated with an index representing its position in
	/// some master list somewhere.  The DaemonSet is responsible for keeping this index.
	void add_pack_daemon(
		Size daemon_index,
		std::string const & pdb_name,
		std::string const & correspondence_file_name,
		std::string const & secondary_resfile
	);

	void
	add_pack_daemon(
		Size daemon_index,
		std::string const & pose_file_name,
		Pose const & pose,
		std::string const & correspondence_file_filename,
		std::istream & correspondence_file,
		std::string const & secondary_refile_file_filename,
		std::istream & secondary_resfile
	);

	void
	add_npd_property_calculator_for_state(
		Size daemon_index,
		std::string const & npd_property,
		Size npd_index
	);

	/// @brief call daemon->setup() on all daemons, which will trigger the
	/// precomputation of all rotamer pair energies.
	void setup_daemons();

	core::Size ndaemons() const;
	core::Size n_npd_properties() const;

	/// @brief Compute the state energies and, for those states requiring non-pairwise-decomposable-properties,
	/// the  non-pairwise decomposable properties as well.  Return them as a pair of lists.
	StateEsAndNPDs
	compute_energy_for_assignment( Entity const & entity );

	ConstDaemonList daemons() const;

	void mark_last_entity_as_important();
	void mark_entity_as_unimportant( Entity const & );

	std::list< std::pair< Size, PoseOP > >
	retrieve_relevant_poses_for_entity( Entity const &, DaemonIndices const & ) const;

	core::pack::task::PackerTaskOP       entity_task() const;
	core::pack::task::ResfileContentsCOP entity_resfile() const;

#ifdef APL_MEASURE_MSD_LOAD_BALANCE
	core::Real last_packing_runtime() const { return packing_runtime_; }
	core::Real last_npd_runtime() const { return npd_runtime_; }
#endif


public:
	/// MPI related methods

	/// @brief receive and respond to remote repacking requests from node 0 until
	/// a spin-down signal is broadcast.  This function does not return until
	/// all MPI requests have completed.
	void activate_daemon_mode();


private:
	/// Private MPI related methods

	/// @brief When we recieve an add_daemons message, recieve all the
	/// input files that go along with the set of daemons that should be added,
	/// and then proceed to read these files.  Send back a shut-down signal
	/// and a message if any of the input files cannot be properly read.  If all
	/// goes well, send back the ok signal.
	void process_add_daemon_message();

	/// @brief When we receive a evaluate_entity signal, calculate
	/// the energies for that entity for each PackDaemon and send those
	/// energies back to the master node.
	void process_state_energy_evaluations_for_entity();

	/// @brief When we receive a discard_entity signal, pass that
	/// message on to the Daemons
	void process_discard_entity_message();

	/// @brief When we receive a pose_request signal, have the desired PackDaemons
	/// create the poses for that entity, and, ship them as PDB strings
	/// back to node 0.
	void process_pose_request_for_entity();

	/// @brief Receive an entity string from node 0 and create an Entity object from that string
	EntityOP recieve_entity() const;

	/// @brief Accept a list of daemon indices which should return a pose based on
	/// their state in the presence of a particular entity.
	DaemonIndices recieve_daemon_inds_requiring_pose_creation() const;

	/// @brief Call MPI_Finalize and exit.
	void graceful_exit() const;

private:
	ScoreFunctionOP                                      score_function_;
	core::Size                                           num_entities_;
	core::pack::task::TaskFactoryOP                      task_factory_;

	bool                                                 include_background_energies_;

	bool                                                 limit_dlig_mem_usage_;
	Size                                                 dlig_nmeg_limit_;


	core::pack::task::PackerTaskOP                       entity_task_;
	core::pack::task::ResfileContentsOP                  entity_resfile_;
	DaemonList                                           daemons_;
	core::Size                                           ndaemons_;

	std::map< std::string, NPDPropCalculatorCreatorOP >  npd_calculator_creators_;

	utility::vector1< PoseOP >                           daemon_poses_;
	utility::vector1< core::pack::task::PackerTaskOP >   daemon_tasks_;
	utility::vector1< std::list< NPDIndAndCalc > >       npd_calcs_for_poses_;
	Size                                                 n_npd_properties_;

#ifdef APL_MEASURE_MSD_LOAD_BALANCE
	core::Real packing_runtime_;
	core::Real npd_runtime_;
#endif

};

////////////////////////////////////////////////
/////////      NPDPropCalculator      //////////
////////////////////////////////////////////////


class NPDPropCalculator : public utility::pointer::ReferenceCount
{
public:
	NPDPropCalculator();
	virtual ~NPDPropCalculator();

	virtual
	void
	setup(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask const & task
	); // no-op default implementation


	virtual
	core::Real
	calculate( core::pose::Pose const & p ) = 0;
};

class NPDPropCalculatorCreator : public utility::pointer::ReferenceCount
{
public:
	NPDPropCalculatorCreator();
	virtual ~NPDPropCalculatorCreator();

	virtual
	std::string
	calculator_name() const = 0;

	virtual
	NPDPropCalculatorOP
	new_calculator() const = 0;
};

////////////////////////////////////////////////
/////////       QuickRepacker         //////////
////////////////////////////////////////////////

class QuickRepacker : public utility::pointer::ReferenceCount
{
public:
	typedef utility::pointer::ReferenceCount                         parent;
	typedef PackDaemon::RotamerAssignmentAndEnergy                   RotamerAssignmentAndEnergy;
	typedef core::pose::PoseOP                                       PoseOP;
	typedef core::pack::task::PackerTaskOP                           PackerTaskOP;
	typedef core::pack::interaction_graph::FixedBBInteractionGraphOP FixedBBInteractionGraphOP;
	typedef core::pack::rotamer_set::RotamerSetsOP                   RotamerSetsOP;

public:
	QuickRepacker(
		PoseOP                 pose,
		PackerTaskOP           task,
		FixedBBInteractionGraphOP ig,
		RotamerSetsOP          rot_sets
	);

	virtual ~QuickRepacker();

	virtual
	RotamerAssignmentAndEnergy
	repack( utility::vector0< int > const & rot_to_pack ) = 0;

protected:
	/// Read access to derived classes
	PoseOP                    pose();
	PackerTaskOP              task();
	FixedBBInteractionGraphOP ig();
	RotamerSetsOP             rot_sets();

	void task( PackerTaskOP ); // allow derived classes to reset their tasks
private:
	PoseOP                    pose_;
	PackerTaskOP              task_;
	FixedBBInteractionGraphOP ig_;
	RotamerSetsOP             rot_sets_;
};

class BasicSimAnnealerRepacker : public QuickRepacker
{
public:
	typedef QuickRepacker parent;

public:
	BasicSimAnnealerRepacker(
		PoseOP                    pose,
		PackerTaskOP              task,
		FixedBBInteractionGraphOP ig,
		RotamerSetsOP             rot_sets
	);

	virtual ~BasicSimAnnealerRepacker();

	virtual
	RotamerAssignmentAndEnergy
	repack( utility::vector0< int > const & rot_to_pack );

private:

	PoseOP                 pose_;
	PackerTaskOP           task_;
	FixedBBInteractionGraphOP ig_;
	RotamerSetsOP          rot_sets_;

};

class RotamerSubsetRepacker : public QuickRepacker
{
public:
	typedef core::pack::rotamer_set::RotamerSubsetsOP RotamerSubsetsOP;
	typedef QuickRepacker parent;
public:
	RotamerSubsetRepacker(
		PoseOP                 pose,
		PackerTaskOP           task,
		FixedBBInteractionGraphOP ig,
		RotamerSetsOP          rot_sets
	);

	~RotamerSubsetRepacker();

public:

	RotamerSubsetsOP
	create_rotamer_subsets_from_rot_to_pack(
		utility::vector0< int > const & rot_to_pack
	);

};

class DenseIGRepacker : public RotamerSubsetRepacker
{
public:
	typedef RotamerSubsetRepacker parent;
	typedef core::pack::interaction_graph::DensePDInteractionGraphOP DensePDInteractionGraphOP;

public:
	DenseIGRepacker(
		PoseOP                    pose,
		PackerTaskOP              task,
		FixedBBInteractionGraphOP ig,
		RotamerSetsOP             rotsets
	);

	virtual ~DenseIGRepacker();
	void set_MCA(); // once done, cannot be undone

	virtual
	RotamerAssignmentAndEnergy
	repack( utility::vector0< int > const & rot_to_pack );

	DensePDInteractionGraphOP
	create_dense_pdig_from_rot_to_pack(
		utility::vector0< int > const & rot_to_pack,
		RotamerSubsetsOP rot_subsets
	);

};

class DoubleDenseIGRepacker : public RotamerSubsetRepacker
{
public:
	typedef RotamerSubsetRepacker parent;
	typedef core::pack::interaction_graph::DoubleDensePDInteractionGraphOP DoubleDensePDInteractionGraphOP;

public:
	DoubleDenseIGRepacker(
		PoseOP                    pose,
		PackerTaskOP              task,
		FixedBBInteractionGraphOP ig,
		RotamerSetsOP             rotsets
	);

	virtual ~DoubleDenseIGRepacker();

	virtual
	RotamerAssignmentAndEnergy
	repack( utility::vector0< int > const & rot_to_pack );

	DoubleDensePDInteractionGraphOP
	create_dense_pdig_from_rot_to_pack(
		utility::vector0< int > const & rot_to_pack,
		RotamerSubsetsOP rot_subsets
	);

};

class FASTER_IG_Repacker : public RotamerSubsetRepacker
{
public:
	typedef RotamerSubsetRepacker parent;
	typedef core::pack::interaction_graph::FASTERInteractionGraphOP FASTERInteractionGraphOP;

public:
	FASTER_IG_Repacker(
		PoseOP                    pose,
		PackerTaskOP              task,
		FixedBBInteractionGraphOP ig,
		RotamerSetsOP             rotsets
	);

	virtual ~FASTER_IG_Repacker();

	virtual
	RotamerAssignmentAndEnergy
	repack( utility::vector0< int > const & rot_to_pack );

	FASTERInteractionGraphOP
	create_faster_ig_from_rot_to_pack(
		utility::vector0< int > const & rot_to_pack,
		RotamerSubsetsOP rot_subsets
	);

	void set_num_sa( int setting );
	void set_sa_scale( core::Real setting );
	void set_ciBR_only( bool setting );

private:
	int num_sa_; // how many rounds of SA followed by sPR should the FASTERAnnealer perform?
	core::Real sa_scale_;
	bool ciBR_only_; // only perform ciBR?

};

}
}

#endif
