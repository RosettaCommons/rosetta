// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/hbnet/HBNet.hh
/// @brief base class HBNet; to detect and design h-bond networks
/// @author Scott Boyken (sboyken@gmail.com)

#ifndef INCLUDED_protocols_hbnet_HBNet_hh
#define INCLUDED_protocols_hbnet_HBNet_hh

#include <protocols/hbnet/HBNet.fwd.hh>
#include <protocols/hbnet/NetworkState.fwd.hh>

// BASIC INCLUDES
#include <basic/datacache/DataMap.fwd.hh>

// UTILITY INCLUDES
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <utility/graph/Graph.fwd.hh>
#include <utility/pointer/memory.hh>

// CORE INCLUDES
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/hbonds/graph/HBondGraph.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pack/rotamer_set/RotamerLinks.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>

// PROTOCOL INCLUDES
#include <protocols/moves/Mover.hh>

#include <set> // MANUAL IWYU
#include <core/scoring/hbonds/graph/HBondGraph.hh> // AUTO IWYU For HBondGraph, HBondEdge (ptr only), HBondNode

#ifdef PYROSETTA
#include <core/scoring/hbonds/HBondSet.hh>
#endif

namespace protocols {
namespace hbnet {

/////////////////////////////////////////////////////////////////////////////////////////////////
//                        HBondResStruct and HBondResStructOP/COP                            //
/////////////////////////////////////////////////////////////////////////////////////////////////

///@brief struct that represents minimal info for residue in an h-bond network
struct HBondResStruct : public utility::VirtualBase {
	core::Size resnum;
	platform::uint rot_index;
	char aa;
	std::string chainid;
	bool is_protein;
	bool is_solvent;
	bool is_ligand;

	HBondResStruct(){}
	HBondResStruct( core::Size const res, platform::uint const rot, char const a, std::string const & c, bool const prot, bool const solv, bool const lig ) :
		resnum(res),
		rot_index(rot), //needed for quick look-ups in ig_ and rotamer_sets_
		aa(a),
		chainid(c),
		is_protein(prot),
		is_solvent(solv),
		is_ligand(lig)
	{}

	bool operator<( HBondResStruct const & a) const
	{
		return resnum < a.resnum;
	}
	bool operator==( HBondResStruct const & a ) const
	{
		return resnum == a.resnum;
	}
};

//compare function for sorting vectors of HBondResStructOP's, by Res # and AA type
struct compare_hbond_residues {
	bool operator()( HBondResStructCOP const & a, HBondResStructCOP const & b) const{
		if ( a->resnum == b->resnum ) {
			return a->aa < b->aa;
		} else {
			return a->resnum < b->resnum;
		}
	}
};

//compare function for sorting vectors of HBondResStructCOP's, by Res #
struct compare_hbond_resnums {
	bool operator()( HBondResStructCOP const & a, HBondResStructCOP const & b) const{
		return a->resnum < b->resnum;
	}
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//                          HBondNetStruct and HBondNetStructOP/COP                          //
/////////////////////////////////////////////////////////////////////////////////////////////////

// TO-DO: HBondNetStruct is big enough now, should probably make it it's own class, HBondNetwork,
//   with pointer types HBondNetworkOP, COP, AP, CAP, and with default operator<
// TO-DO: add graph connectivity and placeholders for bridging waters

///@brief struct that contains info needed for hbond networks
struct HBondNetStruct : public utility::VirtualBase {
	bool is_native;
	bool is_extended;
	bool term_w_bb;                                     //network terminates with sc_bb h-bond
	bool term_w_start;                                  //network cycles to starting res (e.g. ligand)
	bool term_w_cycle;                                  //network has a cycle
	bool scored;                                        //has network been scored?
	bool sort_first_by_tot_unsat;                       //if true, will sort first by total #unsat polar atoms
	bool sort_by_percent_hbond_capacity;                          //if true, will sort next by % percent_hbond_capacity
	bool cst_file_written;
	bool network_pdb_written;
	bool pml_file_written;
	std::string outstring;
	core::Size id;
	core::Size total_hbonds;                            //total h-bonds in the network
	core::Size total_polar_atoms;
	core::Size num_intermolecular_hbs;                  //number of interface h-bonds
	core::Size num_unsat_Hpol;                               // # unsatisfied polar atoms in the entire network (excluding the ligand if ligand_)
	core::Size num_heavy_unsat;                         // how many heavy atoms unsatisfied in network
	core::Size lig_num_unsatisfied;                     //for ligand_, # unsatisfied polar atoms on the ligand
	core::Size num_core_residues;
	core::Size num_boundary_residues;
	core::Real percent_hbond_capacity;                            // % of polar atoms that participate in h-bonds
	core::Real score;                                   //energy score of the network
	utility::vector1< HBondResStructCOP > residues;      //list of residues in the network
	utility::vector1< HBondResStructCOP > asymm_residues; //asymmetric residue list to accurately score networks in symmetric cases (find networks that span entire interface)
	utility::vector1<core::id::AtomID> unsat_Hpols;         //buried unsatisfied polar H's
	utility::vector1<core::id::AtomID> unsat_accs;          //buried unsatisfied acceptors
	utility::vector1< core::scoring::hbonds::HBondCOP > hbond_vec; //all h-bonds in the pose using Rosetta's HBond Object's
	core::scoring::hbonds::HBondSetOP hbond_set;        //HBondSet of all hbonds in the network; HBondSet has useful machinery for searching and unsat checks
	std::vector< platform::uint > lig_state_list;       //rotamer states of other ligand rots compatible with this network
	std::vector< core::Size > net_indices;                    //network_vector_ indices of other h-bond networks that are compatible with this one
	//std::set< core::Size > net_sets;                    //network_vector_ indices of other h-bond networks that are compatible with this one
	utility::vector1<core::conformation::ResidueCOP> rotamers;
	//utility::vector1< core::conformation::ResidueCOP > waterrots;
	//utility::vector1< std::pair< core::id::AtomID, core::conformation::ResidueOP > > waterrots;

	HBondNetStruct() :
		is_native(false),
		is_extended(false),
		term_w_bb(false),
		term_w_start(false),
		term_w_cycle(false),
		scored(false),
		sort_first_by_tot_unsat(true),
		sort_by_percent_hbond_capacity(true),
		cst_file_written(false),
		network_pdb_written(false),
		pml_file_written(false),
		outstring(""),
		id(0),
		total_hbonds(0),
		total_polar_atoms(0),
		num_intermolecular_hbs(0),
		num_unsat_Hpol(0),
		num_heavy_unsat(0),
		lig_num_unsatisfied(0),
		num_core_residues(0),
		num_boundary_residues(0),
		percent_hbond_capacity(0.0),
		score(0.0),
		residues(0),
		asymm_residues(0),
		unsat_Hpols(0),
		unsat_accs(0),
		hbond_vec(0),
		hbond_set(nullptr),
		lig_state_list(0),
		net_indices(0),
		rotamers(0)
		//waterrots(0)
		//network()
	{};

	//copy constructor
	HBondNetStruct( HBondNetStruct const & hbns ) :
		VirtualBase( hbns ),
		is_native(hbns.is_native),
		is_extended(hbns.is_extended),
		term_w_bb(hbns.term_w_bb),
		term_w_start(hbns.term_w_start),
		term_w_cycle(hbns.term_w_cycle),
		scored(hbns.scored),
		sort_first_by_tot_unsat(hbns.sort_first_by_tot_unsat),
		sort_by_percent_hbond_capacity(hbns.sort_by_percent_hbond_capacity),
		cst_file_written(hbns.cst_file_written),
		network_pdb_written(hbns.network_pdb_written),
		pml_file_written(hbns.pml_file_written),
		outstring(hbns.outstring),
		id(0),
		total_hbonds(hbns.total_hbonds),
		total_polar_atoms(hbns.total_polar_atoms),
		num_intermolecular_hbs(hbns.num_intermolecular_hbs),
		num_unsat_Hpol(hbns.num_unsat_Hpol),
		num_heavy_unsat(hbns.num_heavy_unsat),
		lig_num_unsatisfied(hbns.lig_num_unsatisfied),
		num_core_residues(hbns.num_core_residues),
		num_boundary_residues(hbns.num_boundary_residues),
		percent_hbond_capacity(hbns.percent_hbond_capacity),
		score(hbns.score),
		residues(hbns.residues),
		asymm_residues(hbns.asymm_residues),
		unsat_Hpols(hbns.unsat_Hpols),
		unsat_accs(hbns.unsat_accs),
		hbond_vec(hbns.hbond_vec),
		hbond_set(hbns.hbond_set),
		lig_state_list(hbns.lig_state_list),
		net_indices(hbns.net_indices),
		rotamers(hbns.rotamers)
		//waterrots(hbns.waterrots)
		//network()
	{};

	/// @brief Move constructor (mainly for copy/swap idiom)
	/// We're just moving the structure from one address to the other, so a default bitwise copy is appropriate.
	HBondNetStruct( HBondNetStruct && ) = default;

	/// @brief Move assignment operator (mainly for copy/swap idiom)
	HBondNetStruct & operator =( HBondNetStruct && ) = default;

	/// @brief Assigment operator, needed as copy constructor is specified
	HBondNetStruct & operator =( HBondNetStruct const & rhs) {
		// Use the copy/swap idiom
		HBondNetStruct cpy( rhs );
		std::swap( *this, cpy );
		return *this;
	}


	bool operator<( HBondNetStruct const & a ) const
	{
		if ( sort_first_by_tot_unsat ) {
			if ( !( lig_num_unsatisfied == a.lig_num_unsatisfied ) ) { //only happens in ligand case
				return lig_num_unsatisfied < a.lig_num_unsatisfied; // want lowest number of unsats
			} else if ( num_unsat_Hpol == a.num_unsat_Hpol ) {
				if ( sort_by_percent_hbond_capacity && percent_hbond_capacity != a.percent_hbond_capacity ) {
					return percent_hbond_capacity > a.percent_hbond_capacity; // want higher percent_hbond_capacity
				} else {
					return score < a.score;
				}
			} else {
				return num_unsat_Hpol < a.num_unsat_Hpol; // want lowest number of unsats
			}
		} else if ( sort_by_percent_hbond_capacity && percent_hbond_capacity != a.percent_hbond_capacity ) {
			return percent_hbond_capacity > a.percent_hbond_capacity;
		} else {
			return score < a.score; // want lowest score
		}
	}
};

//compare function for sorting vectors of HBondNetStructOP's
// This is needed to compare vectors of OP's; if the vector was of the structrs themselves then
//  the overloading within the struct would work without calling this compare function from std::sort()
// make "inline bool operator()"?? this comparison happens many many times
struct compare_net_vec {
	bool operator()( HBondNetStructOP const & a, HBondNetStructOP const & b) const
	{
		// if exent_existing_networks_=true then use this logic; else is_extends will be false for all networks
		//  network are only flagged as extended at the end for outputing the networks and poses
		if ( ( a->is_extended && !b->is_extended ) || ( !a->is_extended && b->is_extended ) ) {
			return ( a->is_extended && !b->is_extended );
		} else if ( a->sort_first_by_tot_unsat || b->sort_first_by_tot_unsat ) {
			if ( !( a->lig_num_unsatisfied == b->lig_num_unsatisfied ) ) { //only happens in ligand case
				return a->lig_num_unsatisfied < b->lig_num_unsatisfied;
			} else if ( a->num_unsat_Hpol == b->num_unsat_Hpol ) {
				if ( ( a->sort_by_percent_hbond_capacity || b->sort_by_percent_hbond_capacity ) && a->percent_hbond_capacity != b->percent_hbond_capacity ) {
					return a->percent_hbond_capacity > b->percent_hbond_capacity;
				} else {
					return a->score < b->score;
				}
			} else {
				return a->num_unsat_Hpol < b->num_unsat_Hpol;
			}
		} else if ( ( a->sort_by_percent_hbond_capacity || b->sort_by_percent_hbond_capacity ) && a->percent_hbond_capacity != b->percent_hbond_capacity ) {
			return a->percent_hbond_capacity > b->percent_hbond_capacity;
		} else {
			return a->score < b->score;
		}
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////


//MONTE CARLO STRUCTS///////////////////////////////////////////////////////////////////////////////////

struct compare_by_x {
	inline bool operator() ( std::pair< unsigned int, numeric::xyzVector< float > > const & lhs, std::pair< unsigned int, numeric::xyzVector< float > > const & rhs ) const {
		//sorts by x first, then y, then z
		if ( lhs.second.x() != rhs.second.x() ) {
			return lhs.second.x() < rhs.second.x();
		}

		//tiebreaker 1
		if ( lhs.second.y() != rhs.second.y() ) {
			return lhs.second.y() < rhs.second.y();
		}

		//tiebreaker 2
		return lhs.second.z() < rhs.second.z();
	}
};

struct DecoratedNetworkState {
	DecoratedNetworkState(
		core::Real sat,
		unsigned int n_unsat_H,
		NetworkState const * state
	) :
		saturation( sat ),
		n_unsat_Hpols( n_unsat_H ),
		network_state( state )
	{}

	core::Real saturation;
	unsigned int n_unsat_Hpols;
	NetworkState const * network_state;

	bool operator < ( DecoratedNetworkState const & other ) const {
		if ( n_unsat_Hpols == other.n_unsat_Hpols ) {
			return saturation < other.saturation;
		} else {
			return n_unsat_Hpols < other.n_unsat_Hpols;
		}
	}
};


////////////////////////////////////////////////////////////////////////////////////////////////////////


static core::Real const MIN_HB_E_CUTOFF = { -0.1 };
//static core::Real const HB_DIST_CUTOFF = { 3.0 };
static core::Real const SC_RMSD_CUTOFF = { 0.5 };
//static core::Real const DEFAULT_HB_THRESHOLD = { -0.5 };


class HBNet : public moves::Mover
{
public:

	// constructors
	HBNet();
	HBNet( std::string const name );

	HBNet(
		core::scoring::ScoreFunctionCOP scorefxn,
		core::Size max_unsat,
		core::Size min_network_size = 3,
		core::Real hb_threshold = -0.75,
		core::Size max_network_size = 15,
		std::string des_residues = "STRKHYWNQDE",
		bool find_native = false,
		bool only_native = false,
		bool keep_existing=false,
		bool extend_existing=false,
		bool only_extend=false
	);

	//HBNet( HBNet const & other );

	// destructor
	~HBNet() override;

	//virtuals derived from Mover class
	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;
	void apply( core::pose::Pose & pose ) override;

	//optional virtuals that can be derived from HBNet to control behavior
	virtual void setup_packer_task_and_starting_residues( core::pose::Pose const & pose );
	virtual void trim_additional_rotamers( core::pose::Pose & ){}
	virtual void search_IG_for_networks( core::pose::Pose & pose );
	virtual void prepare_output();
	///@brief initial criteria for screening that just reuqires HBondNetStruct (no scoring or placement of rotamers on the pose)
	virtual bool network_meets_initial_criteria( HBondNetStruct const & ){ return true; }
	///@brief final criteria that reuqires network rotamers placed on pose
	virtual bool network_meets_final_criteria( core::pose::Pose const &, HBondNetStruct & ){ return true; }
	virtual bool state_is_starting_aa_type( core::Size const, core::Size const ){ return true; }
	virtual bool pair_meets_starting_criteria( core::Size const, core::Size const, core::Size const, core::Size const ){ return true; }
	core::Real upweight_starting_twobody_energy(){ return upweight_twobody_; }
	virtual core::Real scale_twobody_energy( core::Real input_twobody_energy, char, char ){ return input_twobody_energy; }
	virtual std::string print_additional_info_for_net( HBondNetStruct &, core::pose::Pose const & ){ return ""; }
	virtual std::string print_additional_headers(){ return ""; }
	virtual core::Size ligand(){ return 0; }

	// setters and getters for storing/accessing private vars
	bool symmetric() const { return symmetric_; }
	void symmetric ( bool const is_symmetric_pose ){ symmetric_ = is_symmetric_pose; }
	core::conformation::symmetry::SymmetryInfoCOP get_symm_info() const { return symm_info_; }
	void set_symm_info( core::conformation::symmetry::SymmetryInfoCOP const & symminfo ){ symm_info_ = symminfo; }
	bool multi_component() const { return multi_component_; }
	void multi_component( bool const multi_component ){ multi_component_ = multi_component; }
	bool verbose() const { return verbose_; }
	void verbose( bool const set_verbose ){ verbose_ = set_verbose; }
	bool pdb_numbering() const { return use_pdb_numbering_; }
	void pdb_numbering( bool const pdb_numbering ){ use_pdb_numbering_ = pdb_numbering; }
	std::set< core::Size > get_start_res_vec(){ return start_res_vec_; }
	void add_start_res( core::Size res ){ start_res_vec_.insert( res ); };
	void set_start_res_vec( std::set< core::Size > const start_resnums ){ start_res_vec_ = start_resnums; }
	void set_start_resnums( std::set< core::Size > const start_resnums ){ start_res_vec_ = start_resnums; }
	core::pack::task::PackerTaskOP get_task() const { return task_; }
	void set_task( core::pack::task::PackerTaskOP const & task ){ task_ = task; }
	inline void task_factory( core::pack::task::TaskFactoryOP const & task_factory ){ task_factory_ = task_factory; }
	core::pack::task::TaskFactoryOP task_factory() const{ return task_factory_; }
	void set_core_residues( core::select::residue_selector::ResidueSubset core_residues ){ core_residues_ = core_residues; }
	core::select::residue_selector::ResidueSubset get_core_residues(){ return core_residues_; }
	core::select::residue_selector::ResidueSelectorOP get_core_selector(){ return core_selector_; }
	void set_boundary_residues( core::select::residue_selector::ResidueSubset boundary_residues ){ boundary_residues_ = boundary_residues; }
	core::select::residue_selector::ResidueSubset get_boundary_residues(){ return boundary_residues_; }
	core::select::residue_selector::ResidueSelectorOP get_boundary_selector(){ return boundary_selector_; }
	void update_core_and_boundary_residues( core::pose::Pose const & pose );
	void set_upweight_starting_twobody( core::Real upweight_factor ){ upweight_twobody_ = upweight_factor; }
	///@brief return all hbond networks, sorted with best networks at the front, THIS SHOULD RETURN A POINTER (OP)
	std::vector< HBondNetStructOP > & get_net_vec() // don't want COP here because of sorting / compare; use const_iterator
	{
		std::sort( network_vector_.begin(), network_vector_.end(), compare_net_vec() );
		return network_vector_;
	}
	std::vector< HBondNetStructOP > & get_native_vec() // don't want COP here because of sorting / compare; use const_iterator
	{
		std::sort( native_networks_.begin(), native_networks_.end(), compare_net_vec() );
		return native_networks_;
	}

	std::vector< std::set< core::Size > > & get_output_vector()
	{
		return output_vector_;
	}

	// if call this, be sure to check that id exists and nullptr was not returned
	HBondNetStructOP get_network_by_id( core::Size id )
	{
		for ( auto & net : network_vector_ ) {
			if ( net->id == id ) return net;
		}
		return nullptr; // if id does not exist
	}

	inline core::pose::Pose const & get_orig_pose() { return *orig_pose_; }

	inline void set_orig_pose( core::pose::Pose & pose ){
		orig_pose_ = utility::pointer::make_shared< core::pose::Pose >( pose );
	}

	inline core::pose::Pose const & get_ala_pose() { return *ala_pose_; }

	inline void set_ala_pose( core::pose::Pose & pose ){
		ala_pose_ = utility::pointer::make_shared< core::pose::Pose >( pose );
	}

	char get_aa_for_state( core::Size const res, core::Size const rot ) const;
	bool res_is_boundary( core::Size const res ) const;
	bool res_is_core( core::Size const res ) const;

	void set_min_networks_size( core::Size min ){ min_network_size_ = min; }
	void set_max_networks_size( core::Size max ){ max_network_size_ = max; }
	void set_max_unsat( core::Size max ){ max_unsat_Hpol_ = max; }
	core::Size get_min_networks_size(){ return min_network_size_; }
	core::Size get_max_networks_size(){ return max_network_size_; }
	core::Size get_max_unsat(){ return max_unsat_Hpol_; }
	void set_find_native( bool native ){ find_native_ = native; }
	bool find_native(){ return find_native_; }
	void set_find_only_native( bool only )
	{
		only_native_ = only;
		if ( only_native_ == true ) {
			find_native_ = true;
		}
	}
	bool only_native(){ return only_native_; }
	void set_keep_existing_networks( bool keep_existing_networks ){ keep_existing_networks_ = keep_existing_networks; }
	bool get_keep_existing_networks(){ return keep_existing_networks_; }
	void set_extend_existing_networks( bool extend_existing_networks ){ extend_existing_networks_ = extend_existing_networks; }
	bool get_extend_existing_networks(){ return extend_existing_networks_; }
	void set_only_extend_existing( bool only_extend_existing )
	{
		only_extend_existing_ = only_extend_existing;
		if ( only_extend_existing_ == true ) {
			extend_existing_networks_ = true;
		}
	}
	bool get_only_extend_existing(){ return only_extend_existing_; }


	utility::graph::GraphOP get_packer_graph(){ return packer_neighbor_graph_; }

	core::Size num_core_res( HBondNetStruct const & network );
	core::Size num_boundary_res( HBondNetStruct const & network );
	void select_best_networks();

	//return all hbond networks, sorted with best networks at the front, DEEP COPY
	std::vector< HBondNetStructOP >
	get_networks()
	{
		std::sort( network_vector_.begin(), network_vector_.end(), compare_net_vec() );
		return network_vector_;
	}

	inline void set_score_function( core::scoring::ScoreFunctionCOP const & sf ) {
		runtime_assert( sf != nullptr );
		scorefxn_ = sf->clone();
	}

	///@brief checks if two h-bond networks clash; returns true if they do clash
	bool net_clash(HBondNetStruct const & i, HBondNetStruct const & j);

	void set_symmetry( core::pose::Pose & pose );

	core::Size get_ind_res( core::pose::Pose const & pose, core::Size const res_i);

	bool quick_and_dirty_network_has_heavy_atom_unsat( core::pose::Pose const & pose, HBondNetStruct const & network );
	bool quick_and_dirty_heavy_atom_is_unsat( core::pose::Pose const & pose, core::id::AtomID const at_id );
	bool atom_hbonds_to_bridging_water( core::pose::Pose const & pose, core::id::AtomID const at_id );
	void find_unsats( core::pose::Pose const & pose, HBondNetStruct & i );

	//bool atom_is_buried( core::pose::Pose const & pose, core::id::AtomID id );

	///@brief places the rotamers of the provided h-bond network onto the provided pose
	utility::vector1< core::Size > place_rots_on_pose( core::pose::Pose & pose, HBondNetStruct & residues, bool use_pose=false );

	///@brief return the number of rotamers in a network that are identical in seq or rot to the original input pose:
	core::Size get_num_native_rot(core::pose::Pose & pose, utility::vector1< HBondResStructCOP > const & residues, core::Real sc_rmsd_cut=0.25, bool super=true);
	core::Size get_num_native_seq(core::pose::Pose & pose, utility::vector1< HBondResStructCOP > const & residues);

	// bool
	// water_clashes(
	//  Pose const & pose,
	//  utility::graph::GraphOP packer_neighbor_graph,
	//  core::Size const anchor_i,
	//  core::Size const anchor_j,
	//  core::Vector const water_O
	// );
	//
	// //    bool
	// //    water_clashes(
	// //                  core::pose::Pose const & pose,
	// //                  core::Vector const water_O,
	// //                  core::Real const clash_dist_cut=1.5
	// //                  );
	// //
	// bool
	// water_clashes(
	//  core::pose::Pose const & pose,
	//  core::conformation::Residue const & water
	// );

	bool water_clashes( core::pose::Pose const & pose, core::Vector const water_O );
	bool water_oxygen_clashes_with_residue( core::Vector const water_oxygen, core::Size const resnum, int const rot_state );
	bool water_oxygen_clashes_with_residue( core::Vector const water_oxygen, core::conformation::Residue const & res );

	///@brief used by the job distributor and MultiplePoseMover (MPM) to return poses OTF as requested (without need to store in memory)
	core::pose::PoseOP get_additional_output() override;

	std::string get_file_name( core::Size id, std::string prefix, std::string extension );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	attributes_for_hbnet( utility::tag::AttributeList& );

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void
	provide_citation_info( basic::citation_manager::CitationCollectionList & ) const override;

public: //Monte Carlo Protocol Public methods
	///@brief turn on or off the monte carlo branching protocol
	inline void set_monte_carlo_branch( bool setting ){
		monte_carlo_branch_ = setting;
	}

	///@brief set number of monte carlo runs to be divided over all the seed hbonds. A single monte carlo run appears to take roughly 1 ms (very loose estimate).
	inline void set_total_num_mc_runs( core::Size setting ){
		total_num_mc_runs_ = setting;
	}

	///@brief Maybe you only want to branch from strong hbonds. If this value is -1.2, then only hbonds with a strength of -1.2 or lower will be branched from.
	inline void set_monte_carlo_seed_threshold( core::Real setting ){
		monte_carlo_seed_threshold_ = setting;
	}

	///@brief Only branch from hbonds where at least one residue is buried. Effectively, this results in only finding networks that have at least one buried residue.
	inline void set_monte_carlo_seed_must_be_buried( bool setting ){
		monte_carlo_seed_must_be_buried_ = setting;
	}

	///@brief Only branch from hbonds where both residues are buried. This results in only finding networks that have at least one buried hbond but this does not prevent having additional exposed hbonds.
	inline void set_monte_carlo_seed_must_be_fully_buried( bool setting ){
		monte_carlo_seed_must_be_fully_buried_ = setting;
	}

	// functions to be accessed by classes that inherit from HBNet
protected:

	///@brief sets up rotamer_sets and makes preparations for populating the IG based on task_ops
	void setup( core::pose::Pose const & pose );

	void get_native_networks( core::pose::Pose const & pose );

	///@brief run common steps universal to all HBNet, e.g. finalize rotamer sets and scoring then populate and search IG for networks.
	void run( core::pose::Pose & pose );

	///@brief remove unwanted rotamers from rotamer_sets_ before IG is populated and searched.
	///  calls virtual trim_additional_rotamers() that can be overriden by specialized HBNet classes to remove additional rotmaers
	void trim_rotamers( core::pose::Pose & pose );

	///@brief Do a recursive traversal of an EnergyGraph of a static pose;
	///will find and store all native h-bond networks, regardless of score/parameters.
	void traverse_native(core::pose::Pose const & pose, core::Real const hb_threshold );
	void rec_trav_native(core::pose::Pose const & pose, core::Size new_res, core::Size prev_res,
		utility::vector1< HBondResStructCOP > residues, core::Real const hb_threshold );

	///@brief Recursively traverse the Interaction Graph ig_, and find all possible h-bond networks,
	/// given the parameters that have been set.
	void traverse_IG( core::Real const hb_threshold );
	void recursive_traverse( int const new_node_ind, int const newstate, core::Size const newres, core::Size const prevres, utility::vector1< HBondResStructCOP > residues,
		core::Size network_rec_count, core::Real init_sc, core::Real const hb_threshold, bool const second_search=false );

	///@brief A monte carlo alternative to traverse_IG() and recursive_traverse(). Basically just populates hbond_graph_
	void MC_traverse_IG( );

	///@breif for efficiency, makes sure that an equivalent network has not already been stored
	bool network_already_stored( utility::vector1< HBondResStructCOP > & residues, utility::vector1< HBondResStructCOP > & i_residues );

	// store the h-bond network; writes it to a HBondNetStruct and pushed HBondNetStruct to the back of network_vector_
	void store_network(utility::vector1< HBondResStructCOP > residues,  core::Real init_score=0.0,
		bool term_w_start=false, bool term_w_cycle=false, bool score_now=false, bool native=false );

	///@brief void score_networks( bool minimize=false);
	void score_network_on_pose( core::pose::Pose & pose, HBondNetStruct & i );
	void minimize_network( core::pose::Pose & pose, HBondNetStruct & network, bool residues_already_placed=true );

	///@brief check if a residue clashes with other reisdues in the network before adding it to the h-bond network.  true = clashes
	bool check_clash(utility::vector1< HBondResStructCOP > const & residues, platform::uint my_node_ind,
		core::Size mystate, core::Size myres, core::Real & init_score, bool & cycle);

	///@brief used by net_clash( HBondNetStruct & i, HBondNetStruct & j ), should not be called externally!
	bool net_clash(utility::vector1< HBondResStructCOP > const & residues_i, utility::vector1< HBondResStructCOP > const & residues_j);

	///@brief monte carlo protocol needs to use a different net_clash call because its IG is not populated
	bool monte_carlo_net_clash( utility::vector1< HBondResStructCOP > const & residues_i, utility::vector1< HBondResStructCOP > const & residues_j ) const;

	// I commented this out because it does not seem to be declared anywhere - Jack Maguire
	//bool rotamer_state_compatible( HBondResStructCOP i, HBondResStructCOP j );

	///@brief This function finds all networks that share a common rotamer, and then for ones that do not clash, branches them together in
	///   all possible combination to find all possible networks.
	void branch_overlapping_networks();
	///@brief used by branch_overlapping() to efficiently search for all combinations of compatible networks that can be merged
	void rec_set_intersection(std::vector< core::Size > add_index_vec, std::vector< core::Size > next_index_list, core::Size pos);
	//void alternative_branch_overlapping_networks();
	//void recursive_branch_networks( HBondNetStruct const & original_network, HBondNetStruct const & new_network, HBondNetStructOP output_network, std::set< core::Size > & sets_of_networks_already_added );
	//used by branch_overlapping_networks() to efficiently search for all combinations of compatible networks that can be merged

	///@brief monte carlo alternative to branch_overlapping_networks(). Requires MC_traverse_IG() to be called ahead of time. Theoretically O( total_num_mc_runs_ )
	void MC_branch_overlapping_networks();


	///@brief Merges 2 networks (i,j) into a single h-bond network, new_network.
	void merge_2_branched_networks( HBondNetStruct const & i, HBondNetStruct const & j, HBondNetStructOP new_network );
	void merge_2_branched_networks(utility::vector1< HBondResStructCOP > const & residues1, utility::vector1< HBondResStructCOP > const & residues2, utility::vector1< HBondResStructCOP > & new_residues);

	bool networks_unique( HBondNetStruct const & i, HBondNetStruct const & j, bool no_surface=true );
	bool networks_identical_aa_sequence( HBondNetStruct const & i, HBondNetStruct const & j );
	bool residues_identical( utility::vector1< HBondResStructCOP > & residues1, utility::vector1< HBondResStructCOP > & residues2 );
	bool residues_not_unique( utility::vector1< HBondResStructCOP > & residues1, utility::vector1< HBondResStructCOP > & residues2 );

	///@brief Returns true if a network is a subset of another; it's symmetric, i.e. returns true if i is subset of j, or if j is a subset of i
	bool is_sub_residues( utility::vector1< HBondResStructCOP > & residues1, utility::vector1< HBondResStructCOP > & residues2 );
	bool is_sub_residues( utility::vector1< HBondResStructCOP > & residues1, utility::vector1< HBondResStructCOP > & residues2, bool & branch, bool true_if_identical=true );
	//NOTE cannot use const & residues here because need to sort them and use std:: functions to compare

	///@brief Sorts all networks by score, removes replicate networks and networks > upper_score_limit;
	/// also removes subnetworks if store_subnetworks=false.
	void remove_replicate_networks( core::Size same_max=1 );

	//void pack_with_net_rotamers( core::pose::Pose & pose, core::scoring::ScoreFunctionOP sfxn, std::list<core::conformation::ResidueCOP> & final_res_list, core::Size max_unsat=1 );

	///@brief Writes h-bond networks out to TR, and stores them in output_net_vec_ if finalize=true
	/// so they can be returned by get_additional_output().
	void output_networks( bool finalize=false );

	// for setting up task_operations
	core::pack::task::PackerTaskOP create_ptask( core::pose::Pose const & pose, bool initialize_from_commandline=false );
	bool task_is_valid( core::pose::Pose const & pose ) const;

	bool store_network_scores_in_pose() const {
		return store_network_scores_in_pose_;
	}

	void store_network_scores_in_pose( bool setting ) {
		store_network_scores_in_pose_ = setting;
	}

	// // writes out enzdes .cst constraint file for the h-bond network
	//    void write_cst_files( core::pose::Pose & in_pose, core::pose::Pose & out_pose, HBondNetStructOP o,
	//  std::string cst_fname, bool use_enzdes_cst, bool write_pymol_file=false, bool dump_resfile=false );
	//
	//    void write_pymol_files( std::string name );

	void write_network_pdb( HBondNetStructOP & p );
	void set_constraints( core::pose::Pose & pose, core::scoring::constraints::ConstraintSet & constraints, HBondNetStruct & p, bool write_cst_file=false );


	// ///@brief turns on .cst contraints for the network
	// void turn_on_enzdes_cst( Pose & pose, std::string cst_fname );

protected://Monte Carlo Protocol Protected Methods

	inline
	core::Size chain_num_for_resid( core::Size resid ) const {
		return orig_pose_->chain( resid );
	}

	std::string chain_for_moltenres( core::Size mres ) const;

	inline void set_monte_carlo_data_to_default();

	///@brief This transfers rotamer sets data to hbond_graph_. Need to call this before traversing the IG.
	void initialize_hbond_graph();

	void register_clash( core::Size rotamerA, core::Size rotamerB );

	///@brief The edge_can_yield_monte_carlo_seed() methods determine if a hydrogen bond can be used as a seed for monte carlo branching. Examples of when this might fail include when the hbond does not cross the interface (for HBNetStapleInterface) or when the hbond is not buried to the level of the user's specifications
	virtual bool edge_can_yield_monte_carlo_seed( core::Size resid1, core::Size resid2 ) const;
	virtual bool edge_can_yield_monte_carlo_seed( core::pack::interaction_graph::EdgeBase const * ) const;

	///@brief get_next_node() randomly decides how the current_state will grow during the monte carlo branching protocol
	core::scoring::hbonds::graph::HBondNode * get_next_node( NetworkState & current_state );

	///@brief returns true if this hbond can not create a network
	bool monte_carlo_seed_is_dead_end( core::scoring::hbonds::graph::HBondEdge const * monte_carlo_seed );

	///@brief called at the end of the monte carlo branching protocol. Adds all of the monte carlo networks to the data elements used by the traditional HBNet protocol
	void append_to_network_vector( utility::vector1< DecoratedNetworkState > const & designed_networks );

	void add_residue_to_network_state( NetworkState & current_state, core::scoring::hbonds::graph::HBondNode * node_being_added ) const;

	///@brief quick and dirty satisfaction check
	bool network_state_is_satisfied( NetworkState & current_state ) const;

	///@brief evaluates state and stores the score in state.score
	void score_network_state( NetworkState & state ) const;

	bool register_new_network(
		NetworkState & current_network_state,
		utility::vector0< std::list< utility::vector1< core::Size > > > & hash_table,//hash of elements in current_network_state. hash key is based on sum of rotamer_ids
		std::list< NetworkState > & designed_networks
	) const;

	///@brief returns false if there is a clash between the node_being_added and any node currently in the current_state
	bool node_is_compatible( NetworkState const & current_state, core::scoring::hbonds::graph::HBondNode const * node_being_added );

	///@brief This is only being used for debug purposes right now. Makes sure that there are no possible nodes that can be added to current_state
	bool network_state_is_done_growing( NetworkState const & current_state, core::scoring::hbonds::graph::HBondGraphCOP hbond_graph );

	///@brief TODO
	core::Real estimate_saturation( NetworkState const & ) const;

private:
	core::pose::Pose & nonconst_get_orig_pose() { return *orig_pose_; }

private:
	//bool use_enzdes_cst_;
	bool write_network_pdbs_;                               //write .pdb's of the h-bond networks on Poly-Ala for easy visualization
	bool write_cst_files_;
	bool output_poly_ala_background_;
	bool store_network_scores_in_pose_;
	bool find_native_;
	bool only_native_;                                      //only search for native networks in a static pose
	bool keep_existing_networks_;
	bool extend_existing_networks_;
	bool only_extend_existing_;
	bool verbose_;                                          //print details of mover actions
	bool symmetric_;
	bool multi_component_;
	bool show_task_;
	bool minimize_;
	bool allow_no_hbnets_;
	bool tyr_hydroxyls_must_donate_;
	bool hydroxyls_must_donate_;
	bool use_pdb_numbering_;
	bool no_heavy_unsats_allowed_;
	bool keep_start_selector_rotamers_fixed_;
	//core::Real bw_binary_cut_;
	core::Size min_network_size_;
	core::Size max_network_size_;
	core::Size min_unique_networks_;
	core::Size min_core_res_;
	core::Size min_boundary_res_;
	core::Size max_unsat_Hpol_;
	core::Size max_lig_unsat_;
	core::Size max_rep_;                                    //maximum replicates allowed; default is 1 (no replicate networks)
	core::Size max_replicates_before_branch_;
	core::Size max_replicates_before_unsat_check_;
	std::string const allaas_;
	std::string const hbond_default_;
	std::string const hbond_disallow_default_;
	std::string hbond_disallow_;
	std::string des_residues_;
	std::set< core::Size > start_res_vec_;                  //starting residues for IG traversal; all interface residues if interface_
	std::vector< HBondNetStructOP > network_vector_;        //VECTOR OF ALL STORED H-BOND NETWORKS
	std::vector< HBondNetStructOP > native_networks_;
	std::vector< std::vector< core::Size > > merged_vecs_;
	std::vector< std::set< core::Size > > output_vector_;
	//core::Real pore_radius_;                                //for SASA calculations
	std::map<core::Size,std::pair<core::Size,core::Size> > chain_bounds_; // chain number to residue numbers
	//core::Real atom_burial_cutoff_;
	core::Real hydrogen_bond_threshold_;                    //2-body cutoff; if < means we found h-bond
	core::Real onebody_hb_threshold_;                       //1-body cutoff for symmetric one-residue networks
	core::Real charge_charge_rep_cutoff_;
	core::PackerEnergy clash_threshold_;                    // cutoff; if > then we count as a clash
	core::Real upper_score_limit_;
	core::Real min_percent_hbond_capacity_;
	core::pose::PoseOP ala_pose_;                           //OP to Poly-Ala design/repack shell pose (keeps PRO/GLY/CYS)
	core::pose::PoseOP orig_pose_;                          //OP to original pose
	core::pack::task::PackerTaskOP task_;
	core::pack::task::TaskFactoryOP task_factory_;
	core::scoring::ScoreFunctionOP init_scorefxn_;               //for IG traversal
	core::scoring::ScoreFunctionOP scorefxn_;               //for scoring h-bond networks
	core::scoring::ScoreFunctionOP bw_sfxn_;
	core::pack::rotamer_set::RotamerSetsOP rotamer_sets_;
	core::pack::interaction_graph::InteractionGraphBaseOP ig_;
	core::conformation::symmetry::SymmetryInfoCOP symm_info_;
	core::pack::rotamer_set::RotamerLinksCOP rotamer_links_;
	bool store_subnetworks_;                       //store all subnetworks of detected networks as independent networks to be considered
	bool secondary_search_;
	core::Real secondary_threshold_;
	core::Real upweight_twobody_;
	core::select::residue_selector::ResidueSelectorOP start_selector_;
	core::select::residue_selector::ResidueSelectorOP core_selector_;
	core::select::residue_selector::ResidueSelectorOP boundary_selector_;
	core::select::residue_selector::ResidueSubset core_residues_;
	core::select::residue_selector::ResidueSubset boundary_residues_;
	core::select::residue_selector::ResidueSubset input_hbnet_info_residues_;
	utility::graph::GraphOP packer_neighbor_graph_;
	core::scoring::hbonds::HBondDatabaseCOP hb_database_;

	//MONTE CARLO DATA:
	bool monte_carlo_branch_;
	core::Size total_num_mc_runs_;
	core::Size max_mc_nets_;

	core::scoring::hbonds::graph::HBondGraphOP hbond_graph_;

	//"init state" == "seed hbond"
	utility::vector1< core::scoring::hbonds::graph::HBondEdge const * > monte_carlo_seeds_;
	core::Real monte_carlo_seed_threshold_;//twobody energy threshold
	bool monte_carlo_seed_must_be_buried_;//only branch from hbonds where both residues are in the core
	bool monte_carlo_seed_must_be_fully_buried_;//only branch from hbonds where both residues are in the core

};
// end HBNet

inline void HBNet::set_monte_carlo_data_to_default(){
	monte_carlo_branch_ = false;
	total_num_mc_runs_ = 10000;
	monte_carlo_seeds_.reserve( 1024 );
	monte_carlo_seed_threshold_ = 0;

	monte_carlo_seed_must_be_buried_ = false;
	monte_carlo_seed_must_be_fully_buried_ = false;
	max_mc_nets_ = 0;
}

inline void HBNet::register_clash( core::Size rotamerA, core::Size rotamerB ){
	hbond_graph_->get_node( rotamerA )->register_clash( core::scoring::hbonds::graph::NodeIDSize( rotamerB ) );
	hbond_graph_->get_node( rotamerB )->register_clash( core::scoring::hbonds::graph::NodeIDSize( rotamerA ) );
}


} // hbnet
} // protocols

#endif
