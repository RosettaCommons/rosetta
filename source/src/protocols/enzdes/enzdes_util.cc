// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/enzdes/enzdes_util.cc
/// @brief a bunch of utility functions used in enzdes
/// @author Florian Richter, floric@u.washington.edu


// Unit headers
#include <protocols/enzdes/enzdes_util.hh>

#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintParameters.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>
#include <protocols/enzdes/ModifyStoredLigandRBConfsMovers.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <utility/graph/Graph.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/CacheableString.hh>
#include <core/pack/task/TaskFactory.hh> //task shit
#include <core/pack/task/PackerTask.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/IGEdgeReweightContainer.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/string_util.hh>
#include <boost/foreach.hpp>

#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/id/AtomID.hh>

namespace protocols {
namespace enzdes {
namespace enzutil {

static THREAD_LOCAL basic::Tracer tr( "protocols.enzdes.enzdes_util" );

bool
is_catalytic_seqpos(
	core::pose::Pose const & pose,
	core::Size const seqpos
)
{
	toolbox::match_enzdes_util::EnzdesCacheableObserverCOP enz_obs( toolbox::match_enzdes_util::get_enzdes_observer( pose ) );
	if ( !enz_obs ) return false;
	toolbox::match_enzdes_util::EnzdesCstCacheCOP cst_cache( enz_obs->cst_cache() );
	if ( !cst_cache ) return false;
	for ( core::Size i = 1; i <= cst_cache->ncsts(); ++i ) {
		if ( cst_cache->param_cache( i )->contains_position( seqpos ) ) return true;
	}
	return false;
}


utility::vector1< core::Size >
catalytic_res( core::pose::Pose const & pose)
{

	using namespace core;
	utility::vector1< Size > to_return;

	toolbox::match_enzdes_util::EnzdesCacheableObserverCOP enz_obs( toolbox::match_enzdes_util::get_enzdes_observer( pose ) );
	if ( enz_obs ) {
		toolbox::match_enzdes_util::EnzdesCstCacheCOP cst_cache( enz_obs->cst_cache() );
		if ( cst_cache ) to_return = cst_cache->ordered_constrained_positions( pose );
	}

	if ( to_return.size() == 0 ) {
		for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
			if ( pose.residue_type( i ).is_ligand() ) to_return.push_back( i );
		}
	}
	return to_return;
}


core::Real
sum_constraint_scoreterms(
	core::pose::Pose const & pose,
	int which_res
)
{

	using namespace core::scoring;

	EnergyMap const & all_weights = pose.energies().weights();
	EnergyMap const & scores( which_res == -1 ? pose.energies().total_energies() : pose.energies().residue_total_energies( which_res ) );

	//if( which_res == -1){ //means we want to know stuff for the whole pose
	//  scores = pose.energies().total_energies();
	//}
	//else{ scores = pose.energies().residue_total_energies( which_res ); }

	return scores[ coordinate_constraint ] * all_weights[ coordinate_constraint ] + scores[atom_pair_constraint] * all_weights[ atom_pair_constraint] +
		scores[ angle_constraint ] * all_weights[ angle_constraint ] + scores[ dihedral_constraint ] * all_weights[ dihedral_constraint ];

}


/// @brief convenience function to read in a pose for enzdes
/// this function basically just calls pose_from_pdb, but it
/// will create an enzdes cacheable observer in the pose,
/// also, 1, if a multimodel pdb is requested, if there are different
/// ligand positions, these will be saved in the enzdes cacheable observer
/// and 2, if option additional_packing_ligand_rb_confs is active, random
/// ligand rigid body confs will be generated and added to the cacheable
/// observer
void
read_pose_from_file(
	core::pose::Pose & pose,
	std::string const & filename
)
{
	utility::vector1< core::pose::Pose > input_poses;
	core::import_pose::pose_from_file( input_poses, filename );
	runtime_assert( input_poses.size() > 0 );
	pose = input_poses[1];
	toolbox::match_enzdes_util::EnzdesCacheableObserverOP enz_obs( toolbox::match_enzdes_util::get_enzdes_observer( pose ) );

	//now the basic pose is ready. let's see if there are any additional ligands
	//in multimodel pdbs
	if ( input_poses.size() > 1 ) {
		tr << "PDB " << filename << " is a multimodel pdb containing " << input_poses.size() << " models." << std::endl;
		utility::vector1< utility::vector1< core::conformation::ResidueOP > > add_residue_confs( pose.size() );

		for ( core::Size i = 2; i <= input_poses.size(); ++i ) {
			for ( core::Size j = 1; j <= input_poses[i].size(); ++j ) {
				core::Size pdbres( input_poses[i].pdb_info()->number( j ) );
				char pdbchain( input_poses[i].pdb_info()->chain( j ) );
				char pdbicode( input_poses[i].pdb_info()->icode( j ) );
				core::Size orig_seqpos( pose.pdb_info()->pdb2pose( pdbchain, pdbres, pdbicode ) );
				if ( orig_seqpos > 0 ) {
					core::conformation::ResidueOP newres( input_poses[i].residue(j).clone() );
					newres->seqpos( orig_seqpos );
					newres->chain( pose.residue( orig_seqpos ).chain() );
					add_residue_confs[ orig_seqpos ].push_back( newres );
				}
			}
		}
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( add_residue_confs[i].size() > 0 ) {
				if ( pose.residue(i).is_ligand() ) {
					enz_obs->set_rigid_body_confs_for_lig( i, add_residue_confs[i] );
					tr << "For ligand at seqpos " << i << ", " << add_residue_confs[i].size() << " additional rigid body locations were read in." << std::endl;
				}
			}
		}
	} //if there was more than one pose

	//2.
	if ( basic::options::option[basic::options::OptionKeys::enzdes::additional_packing_ligand_rb_confs] > 0 ) {
		GenerateStoredRBConfs rb_generator( basic::options::option[basic::options::OptionKeys::enzdes::additional_packing_ligand_rb_confs], false );
		rb_generator.apply( pose );
	}
}

void
make_continuous_true_regions_in_bool_vector(
	utility::vector1< bool > & the_vector,
	core::Size const min_number_continuous_trues
)
{

	if (  min_number_continuous_trues > the_vector.size() ) {
		utility_exit_with_message("ridiculous. continuous region is requested to be longer than the actual vector. go play in traffic.\n");
	}

	utility::vector1< std::pair< core::Size, core::Size > > continuous_regions;

	bool in_cont_region = false;
	core::Size last_cont_begin(0);

	for ( core::Size i = 1; i <= the_vector.size(); ++i ) {

		if ( ( ! in_cont_region ) && (the_vector[ i ]==true) ) {
			last_cont_begin = i;
			in_cont_region = true;
			//tr << "cont region begin at " << i << ", ";
		}

		if ( in_cont_region &&  (the_vector[ i ]==false ) ) {
			continuous_regions.push_back( std::pair< core::Size, core::Size > ( last_cont_begin, (i - 1) ));
			in_cont_region = false;
			//tr << "cont_region end at " << (i-1) << " with last_cont_begin "<< last_cont_begin << std::endl;
		}
	}
	//in case the last residue in the vector is set to true
	if ( in_cont_region ) continuous_regions.push_back( std::pair< core::Size, core::Size > (last_cont_begin, the_vector.size() ) );

	if ( continuous_regions.size() == 0 ) {
		utility_exit_with_message("The passed in vector doesn't have a single element set to true.\n");
	}

	for ( core::Size j = 1; j<= continuous_regions.size(); ++j ) {
		std::pair< core::Size, core::Size > & cur_region = continuous_regions[j];

		//tr << "processing cont region between " << cur_region.first << " and " << cur_region.second << std::endl;
		if ( cur_region.second - cur_region.first + 1 < min_number_continuous_trues ) {

			//we need to do something
			core::Size to_fill = min_number_continuous_trues - ( cur_region.second - cur_region.first + 1 );

			core::Size to_fill_each_side = (to_fill / 2) + 1;

			core::Size left_gap(0), right_gap(0);
			core::Size prev_region_first(1), next_region_second( the_vector.size () );

			if ( j == 1 ) left_gap = cur_region.first - 1;
			else {
				left_gap = cur_region.first - continuous_regions[ j - 1 ].second - 1;
				prev_region_first = continuous_regions[ j - 1 ].first;
			}

			if ( j == continuous_regions.size() ) right_gap = the_vector.size() - cur_region.second;
			else {
				right_gap = continuous_regions[ j + 1 ].first - cur_region.second - 1;
				next_region_second = continuous_regions[ j + 1 ].second;
			}

			//tr << "left gap is " << left_gap << " right is " << right_gap << std::endl;
			core::Size cur_region_tmp_first = cur_region.first - to_fill_each_side;
			core::Size cur_region_tmp_second = cur_region.second + to_fill_each_side;

			for ( core::Size k = 1; k <= to_fill_each_side; ++k ) {

				if ( k <= left_gap ) {
					the_vector[ cur_region.first - k ] = true;
					//tr << (cur_region.first - k) << "set true , ";

					if ( k == left_gap )  {
						if ( j == 1 ) cur_region_tmp_first = 1;

						else cur_region_tmp_first = continuous_regions[ j -1 ].second + 1;


						if ( ((cur_region.second + k - 1 ) - prev_region_first + 1) >= min_number_continuous_trues ) {
							cur_region_tmp_second = cur_region.second + k - 1;
							//tr << "breaking because of left gap ";
							break;
						}
					}
				} else if ( ((cur_region.second + k - 1 ) - prev_region_first + 1) >= min_number_continuous_trues ) {
					cur_region_tmp_second = cur_region.second + k - 1;
					break;
				}

				if ( k <= right_gap ) {
					the_vector[ cur_region.second + k ] = true;
					//tr << (cur_region.second + k) << "set true , ";

					if ( k == right_gap ) {
						if ( j == continuous_regions.size() ) cur_region_tmp_second = the_vector.size();
						else cur_region_tmp_second = continuous_regions[ j + 1 ].first - 1;

						if ( (next_region_second - (cur_region.first - k ) + 1 ) >= min_number_continuous_trues ) {
							cur_region_tmp_first =  (cur_region.first - k );
							//tr << "breaking because of right gap ";
							break;
						}
					}
				} else if (  (next_region_second - (cur_region.first - k ) + 1 ) >= min_number_continuous_trues ) {
					cur_region_tmp_first =  (cur_region.first - k );
					break;
				}
			}

			cur_region.first = cur_region_tmp_first;
			cur_region.second = cur_region_tmp_second;
			//tr << std::endl;
		} //if we need to fill up this region

	}
	//tr << " making vector continuous over." << std::endl;

} // make_contiuous_true_regions


core::pack::task::PackerTaskOP
recreate_task(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask const & orig_task
)
{

	using namespace core::pack::task;

	if ( orig_task.total_residue() != pose.size() ) utility_exit_with_message("old task and pose don't have same number of residues.");

	PackerTaskOP mod_task = TaskFactory::create_packer_task( pose );
	mod_task->initialize_from_command_line();

	for ( core::Size i = 1; i <= pose.size(); ++i ) {

		//first, we need to copy the rotamer and rotamerset operations
		for ( auto rot_it = orig_task.residue_task(i).rotamer_operations().begin(); rot_it != orig_task.residue_task(i).rotamer_operations().end(); ++rot_it ) {
			mod_task->nonconst_residue_task( i ).append_rotamer_operation( *rot_it );
		}
		for ( auto rotset_it = orig_task.residue_task(i).rotamer_set_operation_begin(); rotset_it != orig_task.residue_task(i).rotamer_set_operation_end(); ++rotset_it ) {
			mod_task->nonconst_residue_task( i ).append_rotamerset_operation( *rotset_it );
		}

		if ( !orig_task.residue_task( i ).being_packed() ) mod_task->nonconst_residue_task(i).prevent_repacking();
		else if ( !orig_task.residue_task( i ).being_designed() ) mod_task->nonconst_residue_task(i).restrict_to_repacking();
		else {
			utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );

			for ( auto res_it = orig_task.residue_task( i ).allowed_residue_types_begin(); res_it != orig_task.residue_task( i ).allowed_residue_types_end(); ++res_it ) {

				keep_aas[ (*res_it)->aa() ] = true;
			}


			//keep_aas[ core::chemical::aa_cys ] = false;
			mod_task->nonconst_residue_task(i).restrict_absent_canonical_aas( keep_aas );
		}
	}

	if ( orig_task.IGEdgeReweights() ) {
		for ( auto it = orig_task.IGEdgeReweights()->reweighters_begin(); it != orig_task.IGEdgeReweights()->reweighters_end(); ++it ) {
			mod_task->set_IGEdgeReweights()->add_reweighter( *it );
		}
	}

	return mod_task;
} //recreate_task

toolbox::match_enzdes_util::EnzConstraintIOCOP
get_enzcst_io( core::pose::Pose const & pose )
{
	toolbox::match_enzdes_util::EnzdesCacheableObserverCOP enz_obs( toolbox::match_enzdes_util::get_enzdes_observer( pose ) );
	if ( !enz_obs ) return nullptr;
	toolbox::match_enzdes_util::EnzdesCstCacheCOP enzcache( enz_obs->cst_cache() );
	if ( !enzcache ) return nullptr;
	return enzcache->enzcst_io();
}

void
remove_remark_header_for_geomcst(
	core::pose::Pose & pose,
	core::Size geomcst
)
{

	core::io::Remarks remarks( pose.pdb_info()->remarks() );
	core::io::Remarks newremarks;
	for ( core::Size i = 0; i < remarks.size(); ++i ) {

		std::string chainA(""), chainB(""), resA(""), resB("");
		core::Size cst_block(0), exgeom_id( 0 );
		int pdbposA(0), pdbposB(0);
		if ( toolbox::match_enzdes_util::split_up_remark_line( remarks[i].value, chainA, resA, pdbposA, chainB, resB, pdbposB, cst_block, exgeom_id ) ) {
			if ( cst_block == geomcst ) continue;
		}
		newremarks.push_back( remarks[i] );
	}
	pose.pdb_info()->remarks( newremarks );
}

void
create_remark_headers_from_cstcache(
	core::pose::Pose & pose
)
{
	core::io::Remarks remarks( pose.pdb_info()->remarks() );
	core::io::Remarks newremarks;

	toolbox::match_enzdes_util::EnzdesCstCacheCOP cstcache( toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache() );

	if ( cstcache ) {

		for ( core::Size i(1); i <= cstcache->ncsts(); ++i ) {
			toolbox::match_enzdes_util::EnzdesCstParamCache const & param_cache( *(cstcache->param_cache( i )) );
			if ( param_cache.missing_in_pose() ) continue;
			std::string chainA(""), chainB(""), resA(""), resB("");
			int pdbposA(0), pdbposB(0), seqposA(0), seqposB(0);

			if ( param_cache.template_res_cache_size() != 2 ) utility_exit_with_message("Can't build a remark line from a param cache that does not have 2 template res saved.");
			if ( (param_cache.template_res_cache(1)->seqpos_map_size() != 1 ) || (param_cache.template_res_cache(2)->seqpos_map_size() != 1 ) ) utility_exit_with_message("Can't build a remark line from a templateres cahce that does not have 1 position saved.");
			seqposA = param_cache.template_res_cache(1)->seqpos_map_begin()->first;
			seqposB = param_cache.template_res_cache(2)->seqpos_map_begin()->first;

			chainA = pose.pdb_info()->chain( seqposA );
			chainB = pose.pdb_info()->chain( seqposB );
			pdbposA = pose.pdb_info()->number( seqposA );
			pdbposB = pose.pdb_info()->number( seqposB );
			resA = pose.residue( seqposA ).name3();
			resB = pose.residue( seqposB ).name3();

			core::io::RemarkInfo ri;
			ri.num = 666;
			ri.value = toolbox::match_enzdes_util::assemble_remark_line( chainA, resA, pdbposA, chainB, resB, pdbposB, i, 1 );
			//std::cout << "adding " << ri.value << " to header... " << std::endl;
			newremarks.push_back( ri );
		}
	}

	//we also have to make sure that no old remark headers remain
	for ( core::Size i = 0; i < remarks.size(); ++i ) {

		std::string chainA(""), chainB(""), resA(""), resB("");
		core::Size cst_block(0), exgeom_id( 0 );
		int pdbposA(0), pdbposB(0);

		if ( !toolbox::match_enzdes_util::split_up_remark_line( remarks[i].value, chainA, resA, pdbposA, chainB, resB, pdbposB, cst_block, exgeom_id ) ) newremarks.push_back( remarks[i] );
	}
	pose.pdb_info()->remarks( newremarks );
}

/// @detail function not implemented very slick at the moment, need to find way to compile boost regex library :((
/// @detail to extract the pdb code from the pose tag without a regular expression module, some explicit functions
/// @detail have been written below
std::string
get_pdb_code_from_pose_tag( core::pose::Pose const & pose ){
	using namespace core::pose::datacache;

	std::string outtag = pose.data().get_const_ptr< basic::datacache::CacheableString >( CacheableDataType::JOBDIST_OUTPUT_TAG )->str();

	utility::vector1< std::string > pdb_matches;

	//std::cerr <<"string that's supposed to contain a pdb code is " << outtag << ", found to contain the following pdb tags: " << std::endl;

	utility::vector1< std::string > tagparts = utility::string_split( outtag, '_' );

	for ( utility::vector1< std::string >::const_iterator it = tagparts.begin(); it != tagparts.end(); ++it ) {

		if ( it->size() != 4 ) continue;
		std::string cand_str = *it;
		//ok, no boost regex, so clumsy implementation to look for pdb code in string
		if ( is_digit( &cand_str[0]) ) {
			if ( is_digit( &cand_str[1]) && is_digit( &cand_str[2]) && is_digit ( &cand_str[3] ) ) continue;

			if ( (( is_uppercase_letter( & cand_str[1] )|| is_digit( & cand_str[1] ) )
					&& ( is_uppercase_letter( & cand_str[2] )|| is_digit( & cand_str[2] ) )
					&& ( is_uppercase_letter( & cand_str[3] )|| is_digit( & cand_str[3] ) ) )
					||(( is_lowercase_letter( & cand_str[1] )|| is_digit( & cand_str[1] ) )
					&& ( is_lowercase_letter( & cand_str[2] )|| is_digit( & cand_str[2] ) )
					&& ( is_lowercase_letter( & cand_str[3] )|| is_digit( & cand_str[3] ) ) )
					) {

				//std::cerr << "yeah, found putative pdb code " << cand_str << std::endl;
				pdb_matches.push_back( cand_str );
			}
		}
	}
	/*

	//assemble regular expression to match pdb codes:
	//4 char string, first one is a digit, remaining 3 are digits or letters
	boost::regex pdb_re("\d([a-z]{3}|[A-Z]{3})");

	boost::regex pdb_re("\d((\d|[a-z]){3}) | \d((\d|[A-Z]){3})");
	boost::cmatch pdb_rematches;

	if( boost::regex_match(outtag.c_str(), pdb_rematches, pdb_re) ){

	for (core::Size i = 1; i < pdb_rematches.size(); i++){
	pdb_matches.push_back( std::string (pdb_rematches[i].first, pdb_rematches[i].second ) );
	}
	}

	*/

	//for( core::Size i = 1; i <= pdb_matches.size(); ++i) std::cerr << pdb_matches[i] << std::endl;

	if ( pdb_matches.size() == 0 ) {
		std::cerr << "protocols/enzdes/enzdes_util: WARNING: string " << outtag << "does not seem to contain a pdb code, returning N/A. " << std::endl;
		pdb_matches.push_back( "N/A" );
	}

	if ( pdb_matches.size() > 1 ) {
		tr << "WARNING WARNING: in tag " << outtag << ", more than 1 pdbcode like pattern has been identified. assuming the first one (" << pdb_matches[1] << ") is the correct one." << std::endl;
	}

	return pdb_matches[1];

	// for( std::vector< std::string >::const_iterator it = tagparts.begin(); it != tagparts.end(); ++it){
	// if( it->size() != 4 ) continue;
	//}

}

bool
is_digit( char * cha )
{

	//std::cerr << "comparing " << cha[0] << " to digits. " << std::endl;

	if ( cha[0] == '0' ) return true;
	else if ( cha[0] == '1' ) return true;
	else if ( cha[0] == '2' ) return true;
	else if ( cha[0] == '3' ) return true;
	else if ( cha[0] == '4' ) return true;
	else if ( cha[0] == '5' ) return true;
	else if ( cha[0] == '6' ) return true;
	else if ( cha[0] == '7' ) return true;
	else if ( cha[0] == '8' ) return true;
	else if ( cha[0] == '9' ) return true;

	return false;
}


bool
is_uppercase_letter( char * cha)
{

	if ( cha[0] == 'A' ) return true;
	else if ( cha[0] == 'B' ) return true;
	else if ( cha[0] == 'C' ) return true;
	else if ( cha[0] == 'D' ) return true;
	else if ( cha[0] == 'E' ) return true;
	else if ( cha[0] == 'F' ) return true;
	else if ( cha[0] == 'G' ) return true;
	else if ( cha[0] == 'H' ) return true;
	else if ( cha[0] == 'I' ) return true;
	else if ( cha[0] == 'J' ) return true;
	else if ( cha[0] == 'K' ) return true;
	else if ( cha[0] == 'L' ) return true;
	else if ( cha[0] == 'M' ) return true;
	else if ( cha[0] == 'N' ) return true;
	else if ( cha[0] == 'O' ) return true;
	else if ( cha[0] == 'P' ) return true;
	else if ( cha[0] == 'Q' ) return true;
	else if ( cha[0] == 'R' ) return true;
	else if ( cha[0] == 'S' ) return true;
	else if ( cha[0] == 'T' ) return true;
	else if ( cha[0] == 'U' ) return true;
	else if ( cha[0] == 'V' ) return true;
	else if ( cha[0] == 'W' ) return true;
	else if ( cha[0] == 'X' ) return true;
	else if ( cha[0] == 'Y' ) return true;
	else if ( cha[0] == 'Z' ) return true;

	return false;
}


bool
is_lowercase_letter( char * cha)
{

	if ( cha[0] == 'a' ) return true;
	else if ( cha[0] == 'b' ) return true;
	else if ( cha[0] == 'c' ) return true;
	else if ( cha[0] == 'd' ) return true;
	else if ( cha[0] == 'e' ) return true;
	else if ( cha[0] == 'f' ) return true;
	else if ( cha[0] == 'g' ) return true;
	else if ( cha[0] == 'h' ) return true;
	else if ( cha[0] == 'i' ) return true;
	else if ( cha[0] == 'j' ) return true;
	else if ( cha[0] == 'k' ) return true;
	else if ( cha[0] == 'l' ) return true;
	else if ( cha[0] == 'm' ) return true;
	else if ( cha[0] == 'n' ) return true;
	else if ( cha[0] == 'o' ) return true;
	else if ( cha[0] == 'p' ) return true;
	else if ( cha[0] == 'q' ) return true;
	else if ( cha[0] == 'r' ) return true;
	else if ( cha[0] == 's' ) return true;
	else if ( cha[0] == 't' ) return true;
	else if ( cha[0] == 'u' ) return true;
	else if ( cha[0] == 'v' ) return true;
	else if ( cha[0] == 'w' ) return true;
	else if ( cha[0] == 'x' ) return true;
	else if ( cha[0] == 'y' ) return true;
	else if ( cha[0] == 'z' ) return true;

	return false;
}
void
disable_constraint_scoreterms(core::scoring::ScoreFunctionOP scorefxn){

	scorefxn->set_weight(core::scoring::coordinate_constraint, 0.0 );
	scorefxn->set_weight(core::scoring::atom_pair_constraint, 0.0 );
	scorefxn->set_weight(core::scoring::angle_constraint, 0.0 );
	scorefxn->set_weight(core::scoring::dihedral_constraint, 0.0 );
	scorefxn->set_weight(core::scoring::res_type_constraint, 0.0 );

}
void
enable_constraint_scoreterms(core::scoring::ScoreFunctionOP scorefxn){

	if ( scorefxn->has_zero_weight( core::scoring::coordinate_constraint ) ) scorefxn->set_weight(core::scoring::coordinate_constraint, 1 );
	if ( scorefxn->has_zero_weight( core::scoring::atom_pair_constraint ) )  scorefxn->set_weight(core::scoring::atom_pair_constraint, 1 );
	if ( scorefxn->has_zero_weight( core::scoring::angle_constraint ) )      scorefxn->set_weight(core::scoring::angle_constraint, 1 );
	if ( scorefxn->has_zero_weight( core::scoring::dihedral_constraint ) )   scorefxn->set_weight(core::scoring::dihedral_constraint, 1 );
	if ( basic::options::option[basic::options::OptionKeys::enzdes::favor_native_res].user() || basic::options::option[ basic::options::OptionKeys::in::file::pssm ].user() ) {
		if ( scorefxn->has_zero_weight( core::scoring::res_type_constraint ) ) scorefxn->set_weight(core::scoring::res_type_constraint, 1 );
	}
}
bool
is_scofx_cstfied(core::scoring::ScoreFunctionCOP scorefxn){
	if ( scorefxn->has_zero_weight( core::scoring::coordinate_constraint ) && scorefxn->has_zero_weight( core::scoring::atom_pair_constraint ) && scorefxn->has_zero_weight( core::scoring::angle_constraint ) && scorefxn->has_zero_weight( core::scoring::dihedral_constraint ) ) return false;
	else return true;
}

void
scorefxn_update_from_options(core::scoring::ScoreFunctionOP scorefxn){

	if ( basic::options::option[ basic::options::OptionKeys::docking::ligand::old_estat ].user() ) {
		core::scoring::methods::EnergyMethodOptions options_repack( scorefxn->energy_method_options() );
		options_repack.exclude_protein_protein_fa_elec( basic::options::option[basic::options::OptionKeys::docking::ligand::old_estat ]  );
		scorefxn->set_energy_method_options( options_repack );
	}
}

void
remove_all_enzdes_constraints( core::pose::Pose & pose )
{
	protocols::toolbox::match_enzdes_util::EnzdesCstCacheOP cst_cache( toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache() );
	if ( !cst_cache ) return;
	cst_cache->enzcst_io()->remove_constraints_from_pose( pose, false, false );
	toolbox::match_enzdes_util::get_enzdes_observer( pose )->set_cst_cache( nullptr );
}

void
get_resnum_from_cstid_list( std::string const& cstidlist, core::pose::Pose const& pose, utility::vector1<core::Size>& resnums )
{
	resnums.clear();
	utility::vector1< std::string > const cstids ( utility::string_split( cstidlist , ',' ) );
	BOOST_FOREACH ( std::string const cstid, cstids ) {
		if ( cstid=="" ) continue;
		core::Size const resnum (get_resnum_from_cstid( cstid, pose) );
		runtime_assert( resnum>0 && resnum <=pose.size() );
		resnums.push_back( resnum );
	}
	auto last = std::unique(resnums.begin(), resnums.end());
	resnums.erase( last, resnums.end() );
	//tr <<" In util function size of resnums is "<< resnums.size() << std::endl;
}
/// @brief Extracts residue number from cstid string
/// @detail Expects cstid string to be of format [0-9]+[A-Z], where the number is constraint number and trailing letter is temlplate id i.e. either A or B
///     - see enzdes style REMARKs format used to specify match type constraint between two residues (template A and template B)
core::Size
get_resnum_from_cstid( std::string const& cstid, core::pose::Pose const &pose)
{
	char const templ(cstid[ cstid.length()-1] );
	runtime_assert(templ == 'A' || templ == 'B');
	core::Size template_num;
	if ( templ == 'A' ) template_num = 1;
	else template_num = 2;

	std::stringstream ss( cstid.substr( 0, cstid.length() - 1) );
	core::Size cstnum;
	ss >> cstnum;
	//tr << "Cstid " <<cstid <<" parsed as template:"<< templ <<" and cstnum: "<< cstnum<<std::endl;
	protocols::toolbox::match_enzdes_util::EnzdesCstCacheCOP cst_cache( toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache() );
	runtime_assert( cst_cache != nullptr );
	runtime_assert( cstnum <= cst_cache->ncsts());
	bool found (false);
	for ( core::Size seqpos =1; seqpos <=pose.size(); ++seqpos ) {
		if ( cst_cache->param_cache(cstnum)->template_res_cache( template_num )->contains_position( seqpos ) ) {
			//   tr <<"Cstid "<<cstid<<" corresponds to "<< seqpos <<std::endl;
			found = true;
			return seqpos;
		}
	}
	if ( !found ) utility_exit_with_message("Could not parse " + cstid + "to resnum");
	return (0); // should not get here
}

} //enzutil
} //enzdes
} //protocols
