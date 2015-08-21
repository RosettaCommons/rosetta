// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Interface - information about the interface between to partners
/// @brief contains the following information:
///  calculate the interface between the two (or other?) partners
///  set the packer task to only pack the interface
/// @author Monica Berrondo


#include <protocols/scoring/Interface.hh>

// Rosetta Headers
#include <core/conformation/Residue.hh>


#include <core/pack/task/PackerTask.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/EnergyGraph.hh>


#include <basic/options/option.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ Headers

//Utility Headers


#include <numeric/trig.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <basic/Tracer.hh>

// symmetry includes
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/symmetry/SymmetricEnergies.hh>
#include <core/scoring/symmetry/SymmetricEnergies.fwd.hh>
#include <core/pose/symmetry/util.hh>


// option key includes

#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>


using basic::T;
using basic::Error;
using basic::Warning;
using namespace ObjexxFCL;

static thread_local basic::Tracer TR( "core.conformation.Interface" );

namespace protocols {
namespace scoring {

static core::Size max_interchain_sites ( 2 );

////////////////////////////////////////////////////////////////////////////////////
///
/// @brief base for calculating the interface
/// @details
///    decide which type of an interface calculation to use
///
/// @author Monica Berrondo October 19 2007
/////////////////////////////////////////////////////////////////////////////////
void
Interface::calculate( core::pose::Pose const & pose )
{
	/// create a dummy array to initialize all the members of the partner
	/// and is_interface array to false
	using namespace core;
	FArray1D_bool false_array ( pose.total_residue(), false );
	partner_ = false_array;
	is_interface_ = false_array;
	pair_list_.resize( max_interchain_sites );
	contact_list_.resize( pose.total_residue() );
	kinematics::FoldTree const & fold_tree ( pose.fold_tree() );

	// Check if we are symmetric. We don't need the FoldTree partitioning for this
	// so we check for symmetry here
	if ( pose::symmetry::is_symmetric( pose.energies() ) ) {
		symmetric_protein_calculate( pose );
		return;
	}

	// partner is the same as is_upstream of Ian's code in LigandDockProtocol
	if ( !use_input_partners_ ) {
		fold_tree.partition_by_jump( jump_number_, partner_ );
	}

	//for ( Size i=1; i<=pose.total_residue(); ++i ) {
	// assuming if it is not a polymer residue, it must be a ligand
	Size upstream_jump_res, downstream_jump_res;
	upstream_jump_res = fold_tree.upstream_jump_residue( jump_number_ );
	downstream_jump_res = fold_tree.downstream_jump_residue( jump_number_ );
	if ( pose.residue( upstream_jump_res ).is_ligand()  || pose.residue( downstream_jump_res ).is_ligand()  ) {
		TR.Debug << "One of the residues is a ligand, calculating the interface between ligand and protein" << std::endl;
		ligand_calculate( pose );
		return;
		// check if some sort of nucleic acid
	} else if ( pose.residue( upstream_jump_res ).is_NA() || pose.residue( downstream_jump_res ).is_NA() ) {
		TR.Debug << "One of the residues is a nucleic acid, calculating interface of nucleic acids" << std::endl;
		NA_calculate( pose );  //set to protein_calculate(pose)
		return;
	}
	//}
	// if it gets through the entire protein and all the residues were of type protein
	// then use protein-protein interface calculation
	TR.Debug << "Calculating protein-protein interface" << std::endl;
	protein_calculate( pose );
}

////////////////////////////////////////////////////////////////////////////////////
///
/// @brief calculate the protein-protien interface
/// @details
///    decide which type of an interface calculation to use
///    calculate the residues that are at the interface
///    this uses the CAlpha class (at least for now) which
///    gets the distances between c-alpha atoms
///    Returns a vector of bools, true if at interface, false otherwise
///
///    This uses partition_by_jump to determine which residues belong
///    to each partner those on one side of the jump are set to 0,
///    the others are set to 1
///    A residue is at the interface if it is within 8A radius of the
///    residue in question on the other partner
///
/// @references pose_docking_calc_interface from pose_docking.cc
///
/// @author Monica Berrondo June 14 2007
///
/////////////////////////////////////////////////////////////////////////////////
void
Interface::protein_calculate( core::pose::Pose const & pose )
{
	using namespace core;
	using namespace conformation;

	// update_residue_neighbors should be called before calling the energy_graph method
	// however, the pose is const and can't be changed; needs to be fixed!
	// pose.update_residue_neighbors();

	core::scoring::EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	std::vector< int>::iterator new_end_pos;

	for ( Size i=1; i<=Size(energy_graph.num_nodes()); ++i ) {
		for ( graph::Graph::EdgeListConstIter
				iru = energy_graph.get_node(i)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			core::scoring::EnergyEdge const * edge( static_cast< core::scoring::EnergyEdge const *> (*iru) );
			Size const j( edge->get_second_node_ind() );
			if ( partner_(i) == partner_(j) ) continue;
			// if ( is_interface_(i) && is_interface_(j) ) continue;
			Real const cendist = edge->square_distance();
			if ( cendist < distance_squared_ ) {
				//    TR << "interface edge: " << i << ' ' << j << ' ' << cendist << std::endl;
				is_interface_(i) = is_interface_(j) = true;
				pair_list_[1].push_back(i);
				pair_list_[2].push_back(j);

				contact_list_[i].push_back(j);
				contact_list_[j].push_back(i);
			} ///< set pack residue to true if within a predetermined distance (default 8A)
		} ///< for j
	} ///< for i

	/// this is really ugly, there has to be a better way of doing this
	std::sort( pair_list_[1].begin(), pair_list_[1].end() );
	new_end_pos = std::unique( pair_list_[1].begin(), pair_list_[1].end() );
	pair_list_[1].erase( new_end_pos, pair_list_[1].end() );
	std::sort( pair_list_[2].begin(), pair_list_[2].end() );
	new_end_pos = std::unique( pair_list_[2].begin(), pair_list_[2].end() );
	pair_list_[2].erase( new_end_pos, pair_list_[2].end() );
}

// The Rosetta++ criterion for which sidechains repack/minimize in docking:
//  ligand heavy atom within paircutoff(aa,GLY)+1 of aa's CB
// See
//  docking_minimize.cc   docking_MCM_pack_side_chains()
//  docking_movement.cc   docking_repack()
//  docking_scoring.cc    docking_interf_residues()
//  ligand.cc             detect_ligand_interface[_res]()
//                        hetero_atom_amino_acid_distance()
// Ian's criterion to approximate this in Mini:
//  ligand heavy atom within rsd.nbr_radius()+6 of rsd.nbr_atom()
// 6A is an eyeballed magic number to get ~ agreement w/ Rosetta++ paircutoffs+1

void
Interface::NA_calculate( core::pose::Pose const & pose ){ protein_calculate( pose );}

void
Interface::ligand_calculate(
	core::pose::Pose const & pose
)
{

	using namespace core;
	for ( Size i = 1, i_end = pose.total_residue(); i <= i_end; ++i ) {
		// all residues on ligand side can move
		if ( ! partner_(i) ) {
			is_interface_(i) = true;
			continue;
		}
		// on protein side, have to do distance check
		conformation::Residue const & prot_rsd = pose.residue(i);
		for ( Size j = 1, j_end = pose.total_residue(); j <= j_end; ++j ) {
			if ( partner_(j) ) continue; // compare against only ligand residues
			conformation::Residue const & lig_rsd = pose.residue(j);
			for ( Size k = 1, k_end = lig_rsd.nheavyatoms(); k <= k_end; ++k ) {
				double dist2 = lig_rsd.xyz(k).distance_squared( prot_rsd.xyz(prot_rsd.nbr_atom()) );
				double cutoff = prot_rsd.nbr_radius() + 6.0;
				if ( dist2 <= cutoff * cutoff ) {
					is_interface_(i) = true;

					contact_list_[i].push_back(j);
					contact_list_[j].push_back(i);
					goto END_LIGRES_LOOP; // C++ lacks multi-level break  :(
				}
			}
		}
		END_LIGRES_LOOP: ; // compiler needs ; as a no-op before end of loop
	}
}

/// @brief find the nearest residue at the interface to a given residue
/// @author Monica Berrondo November 18, 2010
core::Size
Interface::closest_interface_residue( core::pose::Pose const & pose, Size src_rsd, core::Real & distance )
{
	using namespace core;
	using namespace conformation;

	Size ret_rsd (0);
	Real min_distance (1000000.0);

	for ( Size i=1; i<=pose.total_residue(); ++i ) {
		if ( partner_(src_rsd) == partner_(i) ) continue;
		Real const cendist = pose.residue(i).nbr_atom_xyz().distance(pose.residue(src_rsd).nbr_atom_xyz()) ;
		if ( cendist < min_distance ) {
			ret_rsd = i;
			min_distance = cendist;
		}
	} ///< for j
	distance = min_distance;
	return ret_rsd;
}

////////////////////////////////////////////////////////////////////////////////////
/// @brief print out the interface information
/// @author Monica Berrondo November 07 2007
////////////////////////////////////////////////////////////////////////////////////
void
Interface::print( core::pose::Pose const & pose )
{
	show(pose);
}

void
Interface::show( core::pose::Pose const & pose )
{
	TR << "Interface residues:" << std::endl;

	std::string selection;
	for ( Size i=1; i<=max_interchain_sites; ++i ) {
		TR << "Site " << i << std::endl;
		for ( Size j=1; j<= pair_list_[i].size(); ++j ) {
			core::conformation::Residue const & rsd = pose.residue( pair_list_[i][j] );
			TR << "     " << rsd.aa() << ' ' << rsd.seqpos() << std::endl;
			selection += string_of(rsd.seqpos()) + '+';
		}
		if ( pair_list_[i].size() > 0 ) {
			TR << "     sele interface" << i << ", resi " << selection << std::endl;
		}
		if ( i==1 ) selection.clear();
	}
}

////////////////////////////////////////////////////////////////////////////////////
// need to check to make sure that the interface residues are actually interface
// residues
///
/// @brief sets up which residues are to be packed
/// @details
///
/// @references pose_docking_repack from pose_docking.cc
///
/// @author Monica Berrondo June 14 2007
///
/////////////////////////////////////////////////////////////////////////////////
void
Interface::set_pack(
	core::pose::Pose const & pose,
	core::pack::task::PackerTaskOP task
)
{
	///TODO Get rid of dependency on option system
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// initially set pack_residue to false for all residues
	// apl -- this logic is now inverted to start from a task of "repack everything" and to then
	// produce a task of "repack only a few things"
	task->restrict_to_repacking();
	for ( Size ii=1; ii<=pose.total_residue(); ++ii ) {
		// apl if the residue is set to false (as not being an interface residue), set pack to false
		// Disable packing completely for ligands, not supported yet.
		if ( !is_interface_(ii) || pose.residue(ii).is_ligand() ) {
			task->nonconst_residue_task( ii ).prevent_repacking();
		}
	}

	Size cutpoint ( pose.fold_tree().cutpoint_by_jump( jump_number_ ) );
	// sc - fixed residue selection for norepack1 and norepack2 options to be compatible with docking foldtree
	if ( option[ docking::norepack1 ]() ) {
		for ( Size ii = 1 ; ii <= cutpoint; ++ii ) {
			task->nonconst_residue_task( ii ).prevent_repacking();
		}
	}
	if ( option[ docking::norepack2 ]() ) {
		for ( Size ii = cutpoint ; ii <= pose.total_residue(); ++ii ) {
			task->nonconst_residue_task( ii ).prevent_repacking();
		}
	}
	// set to true for now
	// TODO fix this
	task->or_include_current( true );
}

/// @details
///  Function to determine whether two residues are a "pair" for docking-type
///  scoring calculations, such as vdw and pair across an interface
bool
Interface::is_pair(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2
)
{
	/// this is only necessary for Atom and Residue
	using namespace core;
	using namespace core::conformation;

	bool is_pair = false;

	for ( Size i = 1; i <= contact_list_[rsd1.seqpos()].size(); i++ ) {
		if ( rsd2.seqpos() == contact_list_[rsd1.seqpos()][i] ) is_pair = true;
	}

	return is_pair;
	//
	// Real const cendist ( rsd1.nbr_atom_xyz().distance_squared( rsd2.nbr_atom_xyz() ) );
	// if ( partner_(rsd1.seqpos()) == partner_(rsd2.seqpos()) ) return false;
	// if ( is_interface_(rsd1.seqpos()) && is_interface_(rsd2.seqpos())
	//  && (cendist < distance_ * distance_) )
	//  return true;
	// else
	//  return false;
}

/////////////////////////////////////////////////////////////////////////////////
///
/// @brief calculates the center of mass of interface residues
/// @details
///    calculate the center of mass of the interface
///    loops over all residues at the interface and gets the xyz coordinates for
///    the c-alpha atom of that residue
///
/// @references pose_docking_calc_interface from pose_docking.cc
///
/// @author Monica Berrondo June 14 2007
///
/////////////////////////////////////////////////////////////////////////////////
// check to make sure this is working correctly
core::Vector
Interface::center (
	core::pose::Pose const & pose
)
{
	int count ( 0 );
	std::vector < bool > interface ( pose.total_residue(), false );
	using namespace core;
	Vector center ( 0.0 );

	// first, calculate the residues that are at the interface
	// this should be already calculated?
	//  calculate( pose );
	for ( Size i=1; i<=pose.total_residue(); ++i ) {
		if ( interface[i] ) {
			// the c-alpha atom, this should probably be more generalized so that it can work
			// with a ligand, surface, or dna.
			// in that case it wouldn't be the ca atom, it would be something like the backbone
			// phospate atom from dna, etc.
			// jec generalized to be residue neighbor atom. Not the same as Ca, but much more general.
			Vector const nbr_pos( pose.residue( i ).nbr_atom_xyz() );
			center += nbr_pos;
			count++;
		}
	}
	center /= (count);
	return center;
}

bool Interface::is_interface( core::conformation::Residue const & rsd ) const { return is_interface_( rsd.seqpos() ); }

bool Interface::is_interface( Size const position ) const { return is_interface_(position); }
void Interface::distance( Real const distance_in ) { distance_squared_ = distance_in * distance_in; }
void Interface::jump( Size const jump_num ) { jump_number_ = jump_num; }

//returns the number of residues at the interface
core::Size Interface::interface_nres()
{
	Size nres = 0;
	for ( Size i=1; i<=is_interface_.size(); i++ ) {
		if ( is_interface_(i) ) nres++;
	}
	return nres;
}

// symmetric interfaces

void
Interface::symmetric_protein_calculate( core::pose::Pose const & pose )
{
	using namespace core;
	using namespace conformation;
	using namespace conformation::symmetry;

	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
	core::scoring::symmetry::SymmetricEnergies const & energies
		( dynamic_cast< core::scoring::symmetry::SymmetricEnergies const & > ( pose.energies() ) );
	core::scoring::EnergyGraph const & energy_graph( energies.energy_graph() );

	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( symm_info->bb_is_independent(i) ) {
			partner_(i) = true;
		}
	}
	std::vector< int>::iterator new_end_pos;

	for ( Size i=1; i<=Size(energy_graph.num_nodes()); ++i ) {
		for ( graph::Graph::EdgeListConstIter
				iru = energy_graph.get_node(i)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			core::scoring::EnergyEdge const * edge( static_cast< core::scoring::EnergyEdge const *> (*iru) );
			Size const j( edge->get_second_node_ind() );
			bool symm_add;
			if ( basic::options::option[basic::options::OptionKeys::matdes::num_subs_building_block].user() ) {
				Size num_subs = basic::options::option[basic::options::OptionKeys::matdes::num_subs_building_block]();
				symm_add = (symm_info->subunit_index(i) <= num_subs && symm_info->subunit_index(j) > num_subs) || (symm_info->subunit_index(j) <= num_subs && symm_info->subunit_index(i) > num_subs); // NK
			} else {
				symm_add = ( ( (symm_info->bb_is_independent(i) && !symm_info->bb_is_independent(j)) ) );//||
				//(symm_info->bb_is_independent(i) && !symm_info->bb_is_independent(j)) ) );
			}
			if ( !symm_add ) continue;
			Size i_sym = i;
			Size j_sym = j;
			//if ( is_interface_(i_sym) && is_interface_(j_sym) ) continue; //commented out by Sid
			Real const cendist = edge->square_distance();
			if ( cendist < distance_squared_ ) {
				// TR << "interface edge: " << i << ' ' << j << ' ' << cendist << std::endl;
				is_interface_(i_sym) = is_interface_(j_sym) = true;
				pair_list_[1].push_back(i_sym);
				pair_list_[2].push_back(j_sym);

				contact_list_[i_sym].push_back(j_sym);
				contact_list_[j_sym].push_back(i_sym);

			} ///< set pack residue to true if within a predetermined distance (default 8A)
		} ///< for j
	} ///< for i

	/// this is really ugly, there has to be a better way of doing this
	std::sort( pair_list_[1].begin(), pair_list_[1].end() );
	new_end_pos = std::unique( pair_list_[1].begin(), pair_list_[1].end() );
	pair_list_[1].erase( new_end_pos, pair_list_[1].end() );
	std::sort( pair_list_[2].begin(), pair_list_[2].end() );
	new_end_pos = std::unique( pair_list_[2].begin(), pair_list_[2].end() );
	pair_list_[2].erase( new_end_pos, pair_list_[2].end() );
}

void
Interface::set_symmetric_pack(
	core::pose::Pose const & pose,
	core::pack::task::PackerTaskOP task
)
{
	set_pack( pose, task );
	using namespace core::conformation::symmetry;

	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( !symm_info->chi_is_independent(i) ) {
			task->nonconst_residue_task( i ).prevent_repacking();
		}
	}
}

} // namespace scoring
} // namespace protocols
