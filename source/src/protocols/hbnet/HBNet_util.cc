// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/hbnet/hbnet_util.cc
/// @brief utility functions used by many of the HBNet movers and filters
/// @author Scott Boyken (sboyken@gmail.com)

//Headers
#include <protocols/hbnet/HBNet_util.hh>
#include <protocols/hbnet/HBNet.hh>
//#include <protocols/hbnet/HBNet.fwd.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/graph/Graph.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/symmetry/util.hh>
//#include <core/io/Remarks.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
//#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/sasa.hh>
//#include <core/graph/Graph.hh>
#include <core/id/AtomID.hh>
#include <ObjexxFCL/FArray2D.hh>


namespace protocols {
namespace hbnet {

using namespace core;
using namespace core::pose;
using namespace core::scoring::hbonds;

std::string
print_list_to_string( utility::vector1< HBondResStructCOP > const & residues, bool chainid/* true */, bool term_w_start/*=false*/,
	bool term_w_cycle/*=false*/, bool term_w_bb/*=false*/ )
{
	std::stringstream ret_str;
	Size count(0);
	for ( const auto & residue : residues ) {
		if ( count > 0 ) {
			ret_str << ",";
		}
		if ( chainid ) {
			ret_str << residue->chainid << "_";
		}
		ret_str << residue->aa  << "_" << residue->resnum << ",";
		count++;
	}
	if ( term_w_bb ) {
		ret_str << "backbone, ";
	}
	if ( term_w_cycle ) {
		ret_str << "cycle_";
		if ( term_w_start ) {
			ret_str << "start";
		}
	}
	return ret_str.str();
}

std::string
print_list_to_string( Pose const & pose, utility::vector1< HBondResStructCOP > const & residues, bool chainid, bool term_w_start/*=false*/,
	bool term_w_cycle/*=false*/, bool term_w_bb/*=false*/, bool use_pdb_numbering/*=true*/ )
{
	Size total( pose.total_residue() );
	//need to fix for water
	//        if ( core::pose::symmetry::is_symmetric( *pose ) ){
	//            core::conformation::symmetry::SymmetricConformation const & SymmConf(dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose->conformation()));
	//            core::conformation::symmetry::SymmetryInfoCOP symm_info = SymmConf.Symmetry_Info();
	//            total = symm_info->num_independent_residues();
	//        }
	std::stringstream ret_str;
	Size count(0);
	for ( const auto & residue : residues ) {
		if ( count > 0 ) {
			ret_str << ",";
		}
		if ( chainid ) {
			ret_str << residue->chainid << "_";
		}
		ret_str << residue->aa << "_";
		if ( use_pdb_numbering && ( residue->resnum <= total ) && pose.pdb_info() ) {
			ret_str << pose.pdb_info()->number(residue->resnum); //ASSUMING RESNUMS DON"T CHANGE...COULD BE TRICKY WITH WATER
		} else {
			ret_str << residue->resnum;
		}
		count++;
	}
	if ( term_w_bb ) {
		ret_str << "backbone, ";
	}
	if ( term_w_cycle ) {
		ret_str << "cycle_";
		if ( term_w_start ) {
			ret_str << "start";
		}
	}
	return ret_str.str();
}

//BETTER TO PASS REFERENCES THAT OP's HERE SINCE WE ARE OUTSIDE OF ANY CLASS
std::string
print_network( hbond_net_struct & i, bool chainid /* true */ )
{
	utility::vector1< HBondResStructCOP > const residues( (i.asymm_residues.empty()) ? i.residues : i.asymm_residues );
	//utility::vector1< HBondResStructCOP > const residues( i.residues );
	std::string net_prefix("");
	if ( i.is_native ) net_prefix = "native";
	else if ( i.is_extended ) net_prefix = "extended";
	std::stringstream output;
	output << net_prefix << "network_" << i.id << "\t" << print_list_to_string( residues, chainid, i.term_w_start, i.term_w_cycle, i.term_w_bb ) << "\t"<< residues.size() << "\t" << i.score << "\t" << i.total_hbonds << "\t" << i.connectivity << "\t" << i.num_unsat << "\t";
	return output.str();
}

std::string
print_network_w_pdb_numbering( Pose const & pose, hbond_net_struct const & i, bool chainid /* true */ )
{
	utility::vector1< HBondResStructCOP > const residues( (i.asymm_residues.empty()) ? i.residues : i.asymm_residues );
	//utility::vector1< HBondResStructCOP > const & residues( (i.residues ) );
	std::string net_prefix("");
	if ( i.is_native ) net_prefix = "native";
	else if ( i.is_extended ) net_prefix = "extended";
	std::stringstream output;
	output << net_prefix << "network_" << i.id << "\t" << print_list_to_string( pose, residues, chainid, i.term_w_start, i.term_w_cycle, i.term_w_bb ) << "\t"<< residues.size() << "\t" << i.score << "\t" << i.total_hbonds << "\t" << i.connectivity << "\t" << i.num_unsat << "\t";
	return output.str();
}

std::string
print_headers()
{
	return "HBNet_rank \t residues \t size \t score \t num_hbonds \t connectivity \t num_unsat \t";
}

utility::vector1< HBondResStructCOP >::const_iterator
find_hbond_res_struct( utility::vector1< HBondResStructCOP > const & residues, Size resnum )
{
	auto r = residues.begin();
	for ( ; r != residues.end(); ++r ) {
		if ( resnum == (*r)->resnum ) {
			return r;
		}
	}
	return r;
}

//utility::vector1< HBondCOP >
void
//get_hbond_atom_pairs( hbond_net_struct & network, Pose & pose, bool bb_exclusion /* false */, Real hb_e_cutoff /* HB_E_CUTOFF */ )
//get_hbond_atom_pairs( utility::vector1< HBondResStructCOP > const & residues, Pose & pose, core::utility::graph::GraphOP packer_neighbor_graph, Real hb_e_cutoff /* HB_E_CUTOFF */, bool bb_sc /* true */ )
get_hbond_atom_pairs( hbond_net_struct & network, Pose & pose, bool bb_exclusion /* false */ )
{
	//std::cout << "inside get_hbond_atom_pairs" << std::endl;
	//utility::vector1 < std::pair< id::AtomID, id::AtomID > > hbond_atom_pairs(0);
	utility::vector1< HBondCOP > hbond_vec(0);

	HBondSet temp_hbond_set;
	core::scoring::hbonds::HBondOptions new_options( temp_hbond_set.hbond_options() );
	new_options.use_hb_env_dep(false);
	new_options.bb_donor_acceptor_check(bb_exclusion); // don't use bb exclusion logic when penalizing unsatisfied -- ideally would only ecluse N-H donors and not exclude C=O with only 1 h-bond
	HBondSetOP full_hbond_set( new HBondSet(new_options) );

	//setting this to false will calculate all hbonds in pose
	//full_hbond_set.setup_for_residue_pair_energies( pose, false/*calculate_derivative*/, false/*backbone_only*/ );

	core::scoring::hbonds::fill_hbond_set( pose, false /* deriv */, *full_hbond_set, true /* exclude bb-bb */, false /* exclude bb-sc */, false /* exclude sc-bb */, false /* exclude sc-sc */ );

	for ( utility::vector1< HBondResStructCOP >::const_iterator i = network.residues.begin(); i != network.residues.end(); ++i ) {
		Size res_i((*i)->resnum);

		//std::cout << "getting hbonds for residue " << res_i;
		utility::vector1< HBondCOP > hbonds_for_res_i( full_hbond_set->residue_hbonds(res_i, false) );
		for ( auto & ih : hbonds_for_res_i ) {
			//std::cout << ": arsd = " << (*ih)->acc_res() << "; drsd = " << (*ih)->don_res();
			// NEED TO FIX; may need to uncomment
			//if ( ((*ih)->acc_atm_is_protein_backbone() && (*ih)->acc_res() == res_i) || ((*ih)->don_hatm_is_protein_backbone() && (*ih)->don_res() == res_i) ) {
			// continue;
			//}
			//   if ( (*ih)->energy() <= hb_e_cutoff && !( hbond_exists_in_vector( hbond_vec, *ih ) ) ) { //unweighted h-bond energy
			//    hbond_vec.push_back( *ih );
			//std::cout << "stored" << std::endl;
			if ( !( hbond_exists_in_vector( hbond_vec, ih )) ) { //unweighted h-bond energy
				hbond_vec.push_back( ih );
			}
		}
	}
	// //remove duplicates (doing one_way both ways is needed but results in some duplicates)
	// utility::vector1<HBondCOP>::iterator h1 = hbond_vec.begin();
	// for ( ; h1 != hbond_vec.end(); ) {
	//  for ( utility::vector1<HBondCOP>::iterator h2 = h1+1; h2 != hbond_vec.end(); ) {
	//   if ( (*h1)->acc_res() == (*h2)->acc_res() && (*h1)->don_res() == (*h2)->don_res() && (*h1)->acc_atm() == (*h2)->acc_atm() && (*h1)->don_hatm() == (*h2)->don_hatm() ) {
	//    h2 = hbond_vec.erase(h2);
	//   } else {
	//    ++h2;
	//   }
	//  }
	//  ++h1;
	// }
	network.hbond_vec = hbond_vec;
	//network.hbond_set = HBondSetOP(full_hbond_set->get_self_ptr());
	network.hbond_set = full_hbond_set;

	//return hbond_vec;
	//return new_hbond_set;
}

bool
hbond_exists_in_vector( utility::vector1<HBondCOP> const & hbond_vec, HBondCOP h2 )
{
	for ( const auto & h1 : hbond_vec ) {
		if ( h1->acc_res() == h2->acc_res() && h1->don_res() == h2->don_res() && h1->acc_atm() == h2->acc_atm() && h1->don_hatm() == h2->don_hatm() ) {
			return true;
		}
	}
	return false;
}

void
add_reslabels_to_pose( Pose & pose, hbond_net_struct & i, std::string label /* "HBNet" */ )
{
	if ( !( pose.pdb_info() ) ) {
		pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( pose, true ) ) );
	}
	for ( utility::vector1< HBondResStructCOP >::const_iterator res = i.residues.begin(); res != i.residues.end(); ++res ) {
		pose.pdb_info()->add_reslabel( (*res)->resnum, label );
	}
}

Size
get_num_protein_sc_sc_hbonds( Pose & pose, hbond_net_struct & i )
{
	if ( i.hbond_vec.empty() ) {
		return 0;
	}

	Size num_protein_sc_sc_hbonds(0);
	for ( utility::vector1<HBondCOP>::const_iterator h = i.hbond_vec.begin(); h != i.hbond_vec.end(); ++h ) {
		if ( pose.residue((*h)->acc_res()).is_protein() && pose.residue((*h)->don_res()).is_protein() ) {
			num_protein_sc_sc_hbonds++;
		}
	}
	return num_protein_sc_sc_hbonds;
}

Size
get_num_edges_for_res( Size const res, ObjexxFCL::FArray2D_int & path_dists )
{
	//for ( utility::vector1< HBondResStructCOP >::iterator r = i.residues.begin(); r != i.residues.end(); ++r )
	Size num_edges(0);
	int r( static_cast<int>(res));
	int size( static_cast<int>(path_dists.size2()));
	for ( int s = 1; s <= size; ++s ) { // assume path_dists is symmetrical and dim1 = dim2
		if ( path_dists( r, s ) == 1 ) {
			num_edges++;
		}
	}
	return num_edges;
}

void
hbnet_one_body_energies(
	pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSetOP rotset_op,
	core::scoring::ScoreFunction const & sf,
	utility::vector1< core::PackerEnergy > & energies
	//bool add_background_energies
)
{
	std::fill( energies.begin(), energies.end(), core::PackerEnergy( 0.0 ) );

	int const nrotamers = (int)(rotset_op->num_rotamers());
	//Size const theresid = rotset_op->resid();

	for ( int ii = 1; ii <= nrotamers; ++ii ) {
		core::scoring::EnergyMap emap;
		core::conformation::ResidueCOP temprot = rotset_op->rotamer(ii);
		//sf.eval_ci_1b( *temprot, pose, emap );
		sf.eval_cd_1b( *temprot, pose, emap );
		energies[ ii ] += static_cast< core::PackerEnergy > (sf.weights().dot( emap ));
	}
}

//Yes, some of this is highly similar to SymmetricRotamerSet_::compute_one_body_energies(), but it is necessarily different:
// here we ONLY want to store self 2-body interactions between independent residues and its symm_clones, as well as certain cd_1b's
void
hbnet_symm_one_body_energies(
	pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSetOP rotset_op,
	core::scoring::ScoreFunction const & sf,
	core::pack::task::PackerTask const & task,
	utility::graph::GraphCOP packer_neighbor_graph,
	utility::vector1< core::PackerEnergy > & energies
)
{
	if ( !( core::pose::symmetry::is_symmetric( pose ) ) ) {
		return;
	}
	core::conformation::symmetry::SymmetricConformation const & SymmConf(dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation()));
	core::conformation::symmetry::SymmetryInfoCOP symm_info = SymmConf.Symmetry_Info();

	std::fill( energies.begin(), energies.end(), core::PackerEnergy( 0.0 ) );
	utility::vector1< core::PackerEnergy > temp_energies = energies;

	int const nrotamers = (int)(rotset_op->num_rotamers());
	Size const theresid = rotset_op->resid();

	for ( int ii = 1; ii <= nrotamers; ++ii ) {
		core::scoring::EnergyMap emap;
		core::conformation::ResidueCOP temprot = rotset_op->rotamer(ii);
		sf.eval_cd_1b( *temprot, pose, emap );
		energies[ ii ] += static_cast< core::PackerEnergy > (sf.weights().dot( emap ))*(symm_info->score_multiply_factor());
	}

	// define a factory
	core::pack::rotamer_set::RotamerSetFactory rsf;

	//Detect self 2-body interactions and store in 1-body
	for ( utility::graph::Graph::EdgeListConstIter
			ir  = packer_neighbor_graph->get_node( theresid )->const_edge_list_begin(),
			ire = packer_neighbor_graph->get_node( theresid )->const_edge_list_end();
			ir != ire; ++ir ) {
		int const neighbor_id( (*ir)->get_other_ind( theresid ) );
		if ( task.pack_residue( neighbor_id ) ) continue;
		//core::conformation::Residue const & neighbor( pose.residue( neighbor_id ) );
		platform::uint const symm_clone( symm_info->chi_follows( neighbor_id ) );
		if ( symm_clone != 0 && task.pack_residue( symm_clone )
				&& symm_clone != theresid )  continue;
		// Residue IS interating with itself
		if ( symm_clone == theresid ) {
			// We have a self interaction. Self 2-body interactions between ind res and its symm_clones stored in 1-body e
			for ( int jj = 1; jj <= nrotamers; ++jj ) {
				// make a new rotamer set that is going to be translated to the neighbor interation residue
				conformation::ResidueOP sym_rsd( rotset_op->rotamer( jj )->clone() );
				core::pack::rotamer_set::RotamerSetOP one_rotamer_set = rsf.create_rotamer_set( *sym_rsd );
				one_rotamer_set->set_resid( theresid );
				one_rotamer_set->add_rotamer( *sym_rsd );
				// place rotamer set at neighbor position
				//RotamerSetOP sym_rotset(core::pack::rotamer_set::symmetry::SymmetricRotamerSet_::orient_rotamer_set_to_symmetric_partner( pose, sym_rsd, neighbor_id, one_rotamer_set ) );
				core::pack::rotamer_set::RotamerSetOP sym_rotset = rsf.create_rotamer_set( *sym_rsd );
				for ( auto rot = one_rotamer_set->begin(), rot_end = one_rotamer_set->end(); rot != rot_end; ++rot ) {
					core::conformation::Residue target_rsd( *sym_rsd );
					target_rsd.orient_onto_residue( pose.residue( neighbor_id ) );
					target_rsd.copy_residue_connections_from( pose.residue( neighbor_id ) );
					sym_rotset->set_resid( neighbor_id );
					sym_rotset->add_rotamer( target_rsd );
				}

				sf.prepare_rotamers_for_packing( pose, *one_rotamer_set );
				sf.prepare_rotamers_for_packing( pose, *sym_rotset );
				// make a temporary core::PackerEnergy object with on rotamer in it
				ObjexxFCL::FArray2D< core::PackerEnergy > temp_e( 1, 1, 0.0 );
				// evaluate the energy for this rotamer-rotamer interaction
				sf.evaluate_rotamer_pair_energies(*one_rotamer_set, *sym_rotset, pose, temp_e );
				// add the energy of this interaction. Mulitply with the number of subunits in the system
				energies[ jj ] += temp_e[ 0 ]*(symm_info->score_multiply( theresid, neighbor_id ));
			}
		}
	}
}

bool
network_contains_aa( char aa_one_letter, hbond_net_struct const & i )
{
	return network_contains_aa( aa_one_letter, i.residues );
}

bool
network_contains_aa( char aa_one_letter, utility::vector1< HBondResStructCOP > const & residues )
{
	for ( const auto & residue : residues ) {
		if ( residue->aa == aa_one_letter ) {
			return true;
		}
	}
	return false;
}

bool
his_tyr_connectivity( Pose & pose, hbond_net_struct & i )
{
	//HBondSetOP hbond_set(new HBondSet(get_hbond_atom_pairs(pose, i.residues, false)) );
	//utility::vector1<HBondOP> hbond_vec(0);

	if ( i.hbond_vec.empty() ) {
		//i.hbond_vec = get_hbond_atom_pairs( i.residues, pose );
		get_hbond_atom_pairs( i, pose );
	}

	//runtime_assert( !(i.hbond_vec.empty() ) );
	if ( i.hbond_vec.empty() ) {
		return false;
	}

	bool found(false);


	//new way: His must accept hydrogen from Tyr
	for ( utility::vector1<HBondCOP>::const_iterator h = i.hbond_vec.begin(); h != i.hbond_vec.end(); ++h ) {
		Size arsd((*h)->acc_res());
		Size drsd((*h)->don_res());
		char a_aa = pose.residue( arsd ).name1();
		char d_aa = pose.residue( drsd ).name1();
		if ( a_aa == 'H' && !((*h)->acc_atm_is_backbone()) && d_aa == 'Y' && !((*h)->don_hatm_is_backbone()) ) {
			found = true;
		}
	}
	return found;
}

} //hbnet
} //protocols
