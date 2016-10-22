// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/hbnet/hbnet_util.hh
/// @brief utility functions used by many of the HBNet movers and filters
/// @author Scott Boyken (sboyken@gmail.com)


#ifndef INCLUDED_protocols_hbnet_HBNet_util_hh
#define INCLUDED_protocols_hbnet_HBNet_util_hh

#include <protocols/hbnet/HBNet.hh>
#include <protocols/hbnet/HBNet.fwd.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/graph/Graph.hh>

#include <core/types.hh>
//#include <core/graph/Graph.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/id/AtomID.hh>
#include <ObjexxFCL/FArray2D.hh>

// Unit headers

//Utility Headers

// C++ Headers

namespace protocols {
namespace hbnet {

using namespace core;
using namespace core::pose;

std::string print_list_to_string( utility::vector1< HBondResStructCOP > const & residues, bool chainid=true, bool term_w_start=false,
	bool term_w_cycle=false, bool term_w_bb=false );
std::string print_list_to_string( Pose const & pose, utility::vector1< HBondResStructCOP > const & residues, bool chainid=true, bool term_w_start=false,
	bool term_w_cycle=false, bool term_w_bb=false, bool use_pdb_numbering=true );

//std::string print_network( hbond_net_struct & i, Size net_id, bool chainid=true );
std::string print_network( hbond_net_struct & i, bool chainid=true );
std::string print_network_w_pdb_numbering( core::pose::Pose const & pose, hbond_net_struct const & i, bool chainid /* true */ );
//void print_network( hbond_net_struct & i, Size net_id, std::ostream & out, bool chainid=true );

std::string print_headers();

utility::vector1< HBondResStructCOP >::const_iterator find_hbond_res_struct( utility::vector1< HBondResStructCOP > const & residues, Size resnum );

//Size get_intermolecular_hbonds( hbond_net_struct & i );

//Size get_total_hbonds( hbond_net_struct & i );

///@brief return all of the hbonds in the network
//utility::vector1<core::scoring::hbonds::HBondCOP> get_hbond_atom_pairs( utility::vector1< HBondResStructCOP > const & residues, core::pose::Pose & pose, bool bb_exlcusion=false, core::Real hb_e_cutoff=HB_E_CUTOFF );
//void get_hbond_atom_pairs( hbond_net_struct & network, core::pose::Pose & pose, bool bb_exlcusion=false, core::Real hb_e_cutoff=HB_E_CUTOFF );
void get_hbond_atom_pairs( hbond_net_struct & network, core::pose::Pose & pose, bool bb_exlcusion=false );

bool hbond_exists_in_vector( utility::vector1<core::scoring::hbonds::HBondCOP> const & hbond_vec, core::scoring::hbonds::HBondCOP h2 );

//utility::vector1<core::scoring::hbonds::HBondCOP> get_hbond_atom_pairs( utility::vector1< HBondResStructCOP > const & residues, core::pose::Pose & pose, core::Real hb_e_cutoff=HB_E_CUTOFF, bool bb_sc=true );
//utility::vector1<core::scoring::hbonds::HBondCOP> get_hbond_atom_pairs( utility::vector1< HBondResStructCOP > const & residues, core::pose::Pose & pose, core::graph::GraphOP packer_neighbor_graph, core::Real hb_e_cutoff=HB_E_CUTOFF, bool bb_sc=true );

//utility::vector1<core::scoring::hbonds::HBondCOP> get_hbond_atom_pairs( utility::vector1< HBondResStructCOP > const & residues, core::pose::Pose & pose, bool bb_sc, core::scoring::hbonds::HBondDatabaseCOP & hb_database );

void add_reslabels_to_pose( Pose & pose, hbond_net_struct & i, std::string label="HBNet" );

//ObjexxFCL::FArray2D_int get_path_dists( core::pose::Pose & pose, hbond_net_struct & i );

Size get_num_protein_sc_sc_hbonds( core::pose::Pose & pose, hbond_net_struct & i );

Size get_num_edges_for_res( Size const res, ObjexxFCL::FArray2D_int & path_dists );

//void make_network_assymmetric( core::pose::Pose & pose, hbond_net_struct & i );

//void remove_res_from_hbond_vec( hbond_net_struct & i, Size resnum );

//bool is_donorH_satisfied( core::pose::Pose const & pose, core::graph::GraphOP packer_neighbor_graph, Size resnum, Size datm_ind, Size hatm_ind,
//	core::scoring::hbonds::HBondDatabaseCOP & hb_database );

//bool is_donor_satisfied( core::pose::Pose const & pose, core::graph::GraphOP packer_neighbor_graph, Size resnum, Size datm_ind,
//	core::scoring::hbonds::HBondDatabaseCOP & hb_database );

//bool is_acc_satisfied( core::pose::Pose const & pose, core::graph::GraphOP packer_neighbor_graph, Size resnum, Size aatm_ind,
//	core::scoring::hbonds::HBondDatabaseCOP & hb_database );

//bool is_donorH_satisfied( core::pose::Pose const & pose, Size resnum, Size datm_ind, Size hatm_ind,
//	core::scoring::hbonds::HBondDatabaseCOP & hb_database );

//bool is_donor_satisfied( core::pose::Pose const & pose, Size resnum, Size datm_ind,
//	core::scoring::hbonds::HBondDatabaseCOP & hb_database );

//bool is_acc_satisfied( core::pose::Pose const & pose, Size resnum, Size aatm_ind,
//	core::scoring::hbonds::HBondDatabaseCOP & hb_database );

//void update_sasa( core::pose::Pose & pose, core::id::AtomID_Map< core::Real > & atom_sasa, Real pore_radius );

void hbnet_symm_one_body_energies(
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSetOP rotset_op,
	core::scoring::ScoreFunction const & sf,
	core::pack::task::PackerTask const & task,
	utility::graph::GraphCOP packer_neighbor_graph,
	utility::vector1< core::PackerEnergy > & energies
);
void hbnet_one_body_energies(
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSetOP rotset_op,
	core::scoring::ScoreFunction const & sf,
	utility::vector1< core::PackerEnergy > & energies
);

bool network_contains_aa( char aa_one_letter, hbond_net_struct const & i );
bool network_contains_aa( char aa_one_letter, utility::vector1< HBondResStructCOP > const & residues );
bool his_tyr_connectivity( core::pose::Pose & pose, hbond_net_struct & i );

} //hbnet
} //protocols


#endif
