// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file IO-functionality for enzyme design constraints
/// @brief
/// @author Florian Richter, floric@u.washington.edu

// Unit headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintParameters.hh>

// Package headers
#include <protocols/toolbox/match_enzdes_util/EnzCstTemplateRes.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh> //reading remarks
#include <core/pose/PDBPoseMap.hh> //for PDB-info-to-resid functionality
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/scoring/ScoreFunction.hh> //scoring ambiguous constraints
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/Constraints.hh>

#include <basic/options/option.hh> //options

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// C++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

// Basic Headers
#include <basic/Tracer.hh>

// option key includes
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.toolbox.match_enzdes_util.EnzConstraintIO" );

#ifdef    SERIALIZATION
// Project serialization headers
#include <core/chemical/ResidueTypeSet.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

/// @ brief constructor for EnzConstraintIO class, builds up function types
EnzConstraintIO::EnzConstraintIO (core::chemical::ResidueTypeSetCOP src_restype_set) {
	restype_set_ = src_restype_set;
	//favor_native_constraints_.clear();
	mcfi_lists_.clear();
	cst_pairs_.clear();
}

EnzConstraintIO::~EnzConstraintIO() {}

/// @brief reads the enzyme cstfile and for each block of residue residue constraints, creates a new
/// @brief instance of EnzConstraintParameters
void
EnzConstraintIO::read_enzyme_cstfile( std::string fname ) {

	utility::io::izstream data( fname.c_str() );
	std::istringstream line_stream;

	if ( !data ) {
		std::cerr << "ERROR:: Unable to open constraints file: "
			<< fname << std::endl;
		std::exit( 1 );
	}
	tr.Info << "read enzyme constraints from " << fname << " ...";

	EnzConstraintParametersOP parameters;
	toolbox::match_enzdes_util::MatchConstraintFileInfoListOP mcfil;
	cst_pairs_.clear();
	mcfi_lists_.clear();

	std::string line, key("");
	core::Size counted_blocks(0);

	bool in_variable_block( false );

	while ( !data.eof() ) {
		key = "";
		getline(data,line);
		line_stream.clear();
		line_stream.str(line);
		line_stream >> key;

		if ( key == "VARIABLE_CST::BEGIN" ) {
			counted_blocks++; // APL MOD HERE
			in_variable_block = true;
			mcfil = toolbox::match_enzdes_util::MatchConstraintFileInfoListOP( new toolbox::match_enzdes_util::MatchConstraintFileInfoList( restype_set_ ) );
		}

		if ( key == "VARIABLE_CST::END" ) {

			if ( !in_variable_block ) {
				utility_exit_with_message("Error when reading cstfile. Stray VARIABLE_CST::END tag in file.");
			}
			in_variable_block = false;
			parameters = EnzConstraintParametersOP( new EnzConstraintParameters() );
			parameters->init(counted_blocks, restype_set_, get_self_weak_ptr());
			//parameters->set_mcfi_list( mcfil );
			cst_pairs_.push_back( parameters );
			mcfi_lists_.push_back( mcfil );
		}

		if ( key == "CST::BEGIN" ) {


			if ( !in_variable_block ) {
				mcfil = toolbox::match_enzdes_util::MatchConstraintFileInfoListOP( new toolbox::match_enzdes_util::MatchConstraintFileInfoList( restype_set_ ) );
				counted_blocks++;
			}


			if ( mcfil->read_data( data ) ) {

				if ( !in_variable_block ) {
					parameters = EnzConstraintParametersOP( new EnzConstraintParameters );
					parameters->init(counted_blocks, restype_set_, get_self_weak_ptr());
					//parameters->set_mcfi_list( mcfil );
					cst_pairs_.push_back( parameters );
					mcfi_lists_.push_back( mcfil );
				}
			} else {
				utility_exit_with_message("Undefined error when reading cstfile. Something is wrong with the format (no CST::END tag maybe? ).\n");
			}
		}

	} // file reading

	if ( in_variable_block ) utility_exit_with_message("Error when reading cstfile. VARIABLE_CST::BEGIN tag without corresponding VARIABLE_CST::END tag found.");

	tr.Info << " done, " << cst_pairs_.size() << " cst blocks were read." << std::endl;

	//if we're doing matching or building inverse rotamer trees later on,
	//this information is necessary
	this->determine_target_downstream_res();

} //funtion read enzyme cst


toolbox::match_enzdes_util::MatchConstraintFileInfoListCOP
EnzConstraintIO::mcfi_list( core::Size block ) const
{

	runtime_assert( block <= mcfi_lists_.size() );

	return mcfi_lists_[ block ];
}


/// @brief reads the residue numbers that the constraints will be applied to.
void
EnzConstraintIO::process_pdb_header(
	core::pose::Pose & pose,
	bool accept_missing_blocks)
{


	core::pose::PDBPoseMap PDB_map( pose.pdb_info()->pdb2pose() );
	EnzdesCstCacheOP cst_cache = protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache();

	std::set< Size > found_cst_blocks;

	//std::istringstream line_stream;
	//std::map< Size, std::pair< utility::vector1<Size >, utility::vector1< Size > > > cst_block_to_residues;

	//basic::datacache::BasicDataCache const & pose_cache = pose.data();
	core::pose::PDBInfoCOP pose_pdbinfo = pose.pdb_info();
	core::io::Remarks const & pose_remarks = pose_pdbinfo->remarks();
	std::string /*line, */buffer(""), tag("");
	Size cst_block(0), counted_blocks(0);

	//std::cerr << "There are " << pose_remarks->size() << " remark lines." << std::endl;
	for ( std::vector< core::io::RemarkInfo >::const_iterator remark_it = pose_remarks.begin(), end = pose_remarks.end(); remark_it != end; ++remark_it ) {

		//line_stream.clear();
		//line_stream.str( remark_it->value );

		//std::cerr << "this remark string is: " << remark_it->value << std::endl;
		//std::cerr << "and the number is: " << remark_it->num << std::endl;

		//line_stream >> buffer >> tag;
		std::string remark_line( remark_it->value), resA_type(""), resB_type("");
		int resA_num(0), resB_num(0);
		Size pose_resnumA(0), pose_resnumB(0), ex_geom_id(0);
		std::string resA_chain(""),resB_chain("");
		//if( tag == "TEMPLATE"){
		if ( split_up_remark_line( remark_line, resA_chain, resA_type, resA_num, resB_chain, resB_type, resB_num, cst_block, ex_geom_id ) ) {

			if ( cst_block > cst_pairs_.size() ) utility_exit_with_message("The cst_block given in line:\n"+remark_line+"\n is larger than the number of blocks in the constraint file.");

			if ( ex_geom_id > mcfi_lists_[ cst_block ]->num_mcfis() ) utility_exit_with_message("The external geometry ID specified for cst block "+utility::to_string( cst_block)+" is larger than the number of sub-blocks for that block in the cst file.");

			counted_blocks++;
			cst_pairs_[ cst_block ]->set_mcfi( mcfi_lists_[ cst_block ]->mcfi( ex_geom_id ) );

			//Size resA_pose_num(0), resB_pose_num(0);

			//line_stream >> resA_chain >> resA_type >> resA_num;
			//line_stream >> buffer >> buffer >> resB_chain >> resB_type >> resB_num >> cst_block;
			//if( resA_type.size() == 2 ) resA_type = " " + resA_type;
			//if( resB_type.size() == 2 ) resB_type = " " + resB_type;
			//note: if the chain is '_', this means the pose doesn't have a chain info.
			//we'll set the chain to ' ' to prevent a crash
			if ( resA_chain[0] == '_' ) resA_chain[0] = ' ';
			if ( resB_chain[0] == '_' ) resB_chain[0] = ' ';

			found_cst_blocks.insert( cst_block );
			//debug output
			//tr.Info << "remark debug " << buffer << tag << resA_chain << resA_type << resA_num;
			//tr.Info << resB_chain << resB_type << resB_num << cst_block << std::endl;
			if ( resA_num != 0 ) {
				pose_resnumA = PDB_map.find(resA_chain[0], resA_num);

				if ( pose_resnumA == 0 ) utility_exit_with_message( "residue at chain "+resA_chain+" position "+utility::to_string( resA_num )+" not found in pose.");

			}

			if ( resB_num != 0 ) {
				pose_resnumB = PDB_map.find(resB_chain[0], resB_num);

				if ( pose_resnumB == 0 ) utility_exit_with_message( "residue at chain "+resB_chain+" position "+utility::to_string( resB_num )+" not found in pose.");

			}

			// first do sanity checks for format and consistency of REMARKs with actual atom data in the pdb
			if ( cst_block == 0 || cst_block > cst_pairs_.size() ) {
				std::cerr << "Error: catalytic map in pdb file and information in cst file don't match. Either there is no correctly formatted info given in the REMARK block, or there are more constraint REMARKS than blocks in the .cst file." << std::endl;
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			}
			if ( (resA_chain.size() != 1) || (resB_chain.size() !=1) ) {
				std::cerr << "Error: format in pdb file header is wrong, missing chains of catalytic residues. Information readfor resA_chain is " << resA_chain << ", for resB_chain is " << resB_chain  << std::endl;
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			}
			if ( ( (resA_num != 0 ) && (pose.residue(pose_resnumA).name3() != resA_type ) ) ||
					( (resB_num != 0 ) && (pose.residue(pose_resnumB).name3() != resB_type) ) ) {
				std::cerr << "Error: residue names/positions in catalytic header map in pdb file don't match actual protein residues:" << std::endl;
				std::cerr << "Error: residue " << pose_resnumA << " ( " << resA_chain << " " << resA_num << " ) should be " << resA_type << ", is " << ( pose_resnumA?(pose.residue(pose_resnumA).name3()):"<autofind>" ) << std::endl;
				std::cerr << "Error: residue " << pose_resnumB << " ( " << resB_chain << " " << resB_num << " ) should be " << resB_type << ", is " << ( pose_resnumB?(pose.residue(pose_resnumB).name3()):"<autofind>" ) << std::endl;
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			}

			//then get the correct residues in the pose
			//special feature: if res(A/B)_num is 0, all residues with 3 letter name res(A/B)_type will be added


			if ( resA_num != 0 ) {
				cst_cache->param_cache( cst_block )->template_res_cache( 1 )->set_position_in_pose( pose_resnumA );
			}
			if ( resB_num != 0 ) {
				cst_cache->param_cache( cst_block )->template_res_cache( 2 )->set_position_in_pose( pose_resnumB );
			}
			if ( resA_num == 0 || resB_num == 0 ) {

				bool resA_missing( true ), resB_missing( true );
				if ( resA_num != 0 ) resA_missing = false;
				if ( resB_num != 0 ) resB_missing = false;

				for ( Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
					if ( resA_num == 0 && pose.residue(i).name3() == resA_type ) {
						cst_cache->param_cache( cst_block )->template_res_cache( 1 )->add_position_in_pose( i );
						resA_missing = false;
					}
					if ( resB_num == 0 && pose.residue(i).name3() == resB_type ) {
						cst_cache->param_cache( cst_block )->template_res_cache( 2 )->add_position_in_pose( i );
						resB_missing = false;
					}
				}
				//make sure that even though the respos was declared 0 in the header,
				//at least one of the demanded type was found
				if ( resA_missing ) utility_exit_with_message("Residue with name "+resA_type+" declared in header not found in pose.");
				if ( resB_missing ) utility_exit_with_message("Residue with name "+resB_type+" declared in header not found in pose.");

			} //if one of the residue positions is given as 0 in the header

			//resA_pose_num = pose.pdb_info()->pdb2pose(resA_num, resA_chain);
			//resB_pose_num = pose.pdb_info()->pdb2pose(resB_num, resB_chain);

		}//remark line processing done


	}// end file reading

	//start to process missing blocks: try to find position of residues based on information in cst file
	if ( counted_blocks != cst_pairs_.size() ) {

		for ( Size i = 1; i <= cst_pairs_.size(); ++i ) {

			if ( found_cst_blocks.find( i ) == found_cst_blocks.end() ) {
				//assuming it's mcfi 1 for now
				cst_pairs_[ i ]->set_mcfi( mcfi_lists_[ i ]->mcfi( 1 ) );
				bool resA_missing = !cst_pairs_[i]->nonconst_resA()->find_in_pose_if_missing_from_header( pose );
				bool resB_missing = !cst_pairs_[i]->nonconst_resB()->find_in_pose_if_missing_from_header( pose );

				if ( resA_missing || resB_missing ) {

					//make sure at least one of the residues is in the pose
					if ( !accept_missing_blocks ) {
						std::cerr << "Error: catalytic map in pdb file and information in cst file don't match, unequal number of constraints. should be " << cst_pairs_.size() << ", is " << counted_blocks << std::endl;
						utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
					}
				}
			}
		}
	} //if counted blocks != cst pairs size

}// end of process pdb header info function


/// @brief prepares the class for reading in data from a pose/pdb
//void
//EnzConstraintIO::clear_pose_specific_data()
//{
// favor_native_constraints_.clear();
//}//clear pdb specific data function


void
EnzConstraintIO::generate_pose_specific_data(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scofx) const
{
	for ( core::Size block = 1; block <= cst_pairs_.size(); ++block ) {
		generate_pose_specific_data_for_block( pose, scofx, block );
	}
}// end check consistency function


void
EnzConstraintIO::generate_pose_specific_data_for_block(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scofx,
	core::Size cst_block) const
{

	if ( cst_pairs_[ cst_block]->missing_in_pose(pose) ) {
		tr.Info << "Block " << cst_block << " wasn't found in the pose, so no constraints will be generated." << std::endl;
	} else {
		tr.Info << "checking cst data consistency for block " << cst_block << "... ";
		cst_pairs_[ cst_block]->generate_pose_specific_data( pose, scofx );

		tr.Info << " done" << std::endl;
	}
}

void
EnzConstraintIO::add_constraints_to_pose(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scofx,
	bool accept_blocks_missing_header
)
{
	//tmp hack
	EnzdesCstCacheOP cst_cache = protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache();
	if ( !cst_cache ) protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->set_cst_cache( toolbox::match_enzdes_util::EnzdesCstCacheOP( new EnzdesCstCache( get_self_ptr(), cst_pairs_.size() ) ) );
	//tmp hack over

	//in case this function gets called twice, we remove constraints from the pose
	remove_constraints_from_pose( pose, false, false );

	//new
	protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->set_cst_cache( toolbox::match_enzdes_util::EnzdesCstCacheOP( new EnzdesCstCache( get_self_ptr(), cst_pairs_.size() ) ) );

	process_pdb_header( pose, accept_blocks_missing_header );

	if ( basic::options::option[basic::options::OptionKeys::enzdes::enz_debug] ) show_cst_definitions();

	tr.Info << "Generating constraints for pose... " << std::endl;

	for ( core::Size block = 1; block <= cst_pairs_.size(); ++block ) {

		add_constraints_to_pose_for_block_without_clearing_and_header_processing( pose, scofx, block);

		tr.Info << "Cst Block " << block << "done... " << std::endl;
	}

	//tr.Info << std::endl << "All constraints generated and added. " << std::endl;
} //add_constraints_to_pose functino


void
EnzConstraintIO::add_constraints_to_pose_for_block_without_clearing_and_header_processing(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scofx,
	core::Size cst_block
) const
{

	generate_pose_specific_data_for_block( pose, scofx, cst_block);

	//debug shit
	//utility::vector1< core::scoring::constraints::ConstraintCOP > pcst = cst_pairs_[ cst_block - 1]->active_pose_constraints();
	//core::Size num_to_add = pcst.size();
	//std::cerr << "about to add " << num_to_add << " constraints for block " << cst_block << " to pose." << std::endl;
	//debug shit over

	EnzdesCstParamCacheOP param_cache = protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_block );
	using namespace core::scoring::constraints;
	ConstraintCOPs tmp = pose.add_constraints( param_cache->active_pose_constraints() );
	param_cache->set_active_pose_constraints( tmp );

} // add_constraints_to_pose_for_block_without_clearing_and_header_processing


void
EnzConstraintIO::remove_constraints_from_pose(
	core::pose::Pose & pose,
	bool const keep_covalent,
	bool const fail_on_constraints_missing
) const
{
	EnzdesCstCacheOP cst_cache = protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache();
	if ( !cst_cache ) return;
	//std::cerr << "removing constraints from pose" << std::endl;
	//need to do this in decreasing fashion in case several constraints are covalent
	for ( core::Size block = cst_pairs_.size(); block >= 1; --block ) {
		if ( !cst_cache->param_cache( block )->missing_in_pose() && (! (keep_covalent && cst_pairs_[ block ]->is_covalent() )) ) {
			remove_constraints_from_pose_for_block( pose, block, fail_on_constraints_missing);
		}
	}
} //remove constraints from pose function


void
EnzConstraintIO::remove_constraints_from_pose_for_block(
	core::pose::Pose & pose,
	core::Size cst_block,
	bool const fail_on_constraints_missing
) const
{

	bool constraints_found(false);
	EnzdesCstParamCacheOP param_cache = protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_block );
	if ( pose.remove_constraints( param_cache->active_pose_constraints(), true )  ) {

		constraints_found = true;
	} else if ( fail_on_constraints_missing ) {

		utility::vector1< core::scoring::constraints::ConstraintCOP > pcst = param_cache->active_pose_constraints();
		std::cerr << "trying to remove the following " << pcst.size() << " constraints for cst block " << cst_block << "... " << std::endl;

		for ( utility::vector1< core::scoring::constraints::ConstraintCOP >::const_iterator cst_it = pcst.begin();
				cst_it != pcst.end(); ++cst_it ) {
			(*cst_it)->show( std::cerr );
		}
		utility_exit_with_message("Error: an enzdes constraint that should be in the pose got lost along the way.\n");

	}

	if ( constraints_found && cst_pairs_[cst_block]->is_covalent() ) {
		cst_pairs_[cst_block]->remove_covalent_connections_from_pose( pose );
	}

} //remove_constraints_from_pose_for_block

void
EnzConstraintIO::remove_position_from_template_res_for_block(
	core::pose::Pose & pose,
	core::Size pos,
	core::Size cst_block ) const
{
	protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_block )->remove_seqpos_from_template_res( pos );
}

void
EnzConstraintIO::remove_position_from_template_res(
	core::pose::Pose & pose,
	core::Size pos ) const
{
	EnzdesCstCacheCOP cst_cache( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache() );
	if ( !cst_cache ) return;

	for ( core::Size i =1; i <= cst_pairs_.size(); ++i ) {
		if ( cst_cache->param_cache( i )->contains_position( pos ) ) remove_position_from_template_res_for_block( pose, pos, i );
	}
}


void
EnzConstraintIO::add_pregenerated_constraints_to_pose(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scofx ) const
{

	using namespace core::scoring::constraints;
	//some precautions to avoid segfaults
	EnzdesCstCacheOP cst_cache( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache() );
	if ( !cst_cache ) {
		tr << "Notice: trying to add pregenerated enzdes constraints even though no constraints have been generated, function will have no effect... " << std::endl;
		return;
	}
	//std::cerr << "adding pregenerated constraints" << std::endl;
	utility::vector1< ConstraintCOP > all_pose_constraints = pose.constraint_set()->get_all_constraints();

	for ( core::Size i = 1; i<= cst_pairs_.size(); ++i ) {
		EnzdesCstParamCacheOP param_cache( cst_cache->param_cache( i ) );

		if ( cst_pairs_[i]->missing_in_pose(pose) ) continue;

		if ( ( !  cst_pairs_[i]->is_empty() ) && ( param_cache->active_pose_constraints().size() == 0 )  ) {
			utility_exit_with_message("trying to add pregenerated constraints to the pose even though they haven't been generated yet.");
		}

		bool covalent_kept(false);

		utility::vector1< ConstraintCOP > & cur_active_constraints = param_cache->active_pose_constraints();

		for ( utility::vector1< ConstraintCOP >::iterator cst_it = cur_active_constraints.begin();
				cst_it != cur_active_constraints.end(); ++cst_it ) {

			utility::vector1< ConstraintCOP >::iterator cst_find = find(all_pose_constraints.begin(), all_pose_constraints.end(), *cst_it);

			if ( cst_find != all_pose_constraints.end() ) {
				if ( ! cst_pairs_[i]->is_covalent() ) {
					tr << "WARNING: tried to add an enzdes constraint that's already in the pose. Something's a bit unclean somewhere." << std::endl;
				} else covalent_kept = true;
			} else {
				//std::cerr << "There are a total of " << cur_active_constraints.size() << " constraints in this parameter.\n" << std::endl;
				//std::cerr << "showing definition for a constraint with " << (*cst_it)->natoms() << " atoms... ";
				//(*cst_it)->show( std::cerr );
				if ( !  cst_pairs_[i]->is_covalent() ) *cst_it = pose.add_constraint( *cst_it );
			}

		}
		if ( !covalent_kept &&  cst_pairs_[i]->is_covalent() ) {
			add_constraints_to_pose_for_block_without_clearing_and_header_processing( pose, scofx, i );
		}

	} //loop over cst_pairs

} //add pregenerated constraints function


bool
EnzConstraintIO::contains_position(
	core::pose::Pose const & pose,
	core::Size const seqpos ) const
{
	protocols::toolbox::match_enzdes_util::EnzdesCacheableObserverCOP enz_obs( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose ) );
	if ( !enz_obs ) return false; // get_enzdes_observer() can return NULL for const pose
	EnzdesCstCacheCOP cst_cache( enz_obs->cst_cache() );
	if ( !cst_cache ) return false;
	for ( core::Size i = 1; i <= cst_cache->ncsts(); ++i ) {
		if ( cst_cache->param_cache( i )->contains_position( seqpos ) ) return true;
	}
	return false;
}

bool
EnzConstraintIO::is_backbone_only_cst(
	core::pose::Pose const & pose,
	core::Size const seqpos ) const
{

	bool to_return(false);
	protocols::toolbox::match_enzdes_util::EnzdesCacheableObserverCOP enz_obs( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose ) );
	if ( !enz_obs ) return false; // get_enzdes_observer() can return NULL for const pose
	EnzdesCstCacheCOP cst_cache( enz_obs->cst_cache() );
	if ( !cst_cache ) return false;

	for ( core::Size i = 1; i<= cst_pairs_.size(); ++i ) {

		if ( cst_cache->param_cache(i)->template_res_cache( 1 )->contains_position( seqpos ) ) {
			if ( cst_pairs_[i]->resA()->is_backbone() ) to_return = true;
			else return false;
		}

		if ( cst_cache->param_cache(i)->template_res_cache( 2 )->contains_position( seqpos ) ) {
			if ( cst_pairs_[i]->resB()->is_backbone() ) to_return = true;
			else return false;
		}
	}
	return to_return;
}

void
EnzConstraintIO::update_pdb_remarks_for_backbone_params(
	core::pose::Pose & pose ) const
{

	for ( utility::vector1< EnzConstraintParametersOP >::const_iterator it = cst_pairs_.begin(); it != cst_pairs_.end(); ++it ) {

		if ( !(*it)->resA()->is_backbone() && !(*it)->resB()->is_backbone() ) continue;
		if ( (*it)->missing_in_pose(pose) ) continue;

		if ( !(*it)->update_pdb_remarks( pose ) ) utility_exit_with_message("Error when trying to update pdb remarks.");
	}
}

utility::vector1< std::string >
EnzConstraintIO::allowed_res_name3_at_position(
	core::pose::Pose const & pose,
	core::Size const seqpos ) const
{

	utility::vector1< std::string > to_return;
	if ( !protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache() ) return to_return;

	std::set< std::string > found;

	for ( utility::vector1< EnzConstraintParametersOP >::const_iterator it = cst_pairs_.begin(); it != cst_pairs_.end(); ++it ) {

		std::set< std::string > res_this_param = (*it)->allowed_res_name3_at_position( pose, seqpos );
		if ( res_this_param.empty()/*size() == 0*/ ) continue;
		//we need to make two checks:
		//1. are there already residues in the found set? if so, dont include anything,
		//but remove all that are not part of the second set
		//2. if not, the found set will become res_this_param
		if ( found.empty()/*size() == 0*/ ) found = res_this_param;
		else {
			for ( std::set< std::string >::iterator set_it = found.begin(); set_it != found.end();  ) {
				if ( res_this_param.find( *set_it ) == res_this_param.end() ) {
					std::set< std::string >::iterator to_erase = set_it;
					++set_it;
					found.erase( to_erase );
				} else ++set_it;
			}
		}
	} //iterator over all params
	for ( std::set< std::string >::iterator set_it = found.begin(); set_it != found.end();  ++set_it ) {
		to_return.push_back( *set_it );
	}
	return to_return;
}


void
EnzConstraintIO::show_cst_definitions() const
{
	if ( cst_pairs_.size() == 0 ) {
		tr.Info << "No constraints have been read in." << std::endl;
	} else {
		tr.Info << cst_pairs_.size() << " constraint blocks have been read in: " << std::endl;
		for ( utility::vector1< EnzConstraintParametersOP >::const_iterator it = cst_pairs_.begin(); it != cst_pairs_.end(); ++it ) { (*(*it)).show_definitions(); }
	}
}

void
EnzConstraintIO::remap_resid( core::id::SequenceMapping const & smap )
{

	for ( utility::vector1< EnzConstraintParametersOP >::iterator it = cst_pairs_.begin(); it != cst_pairs_.end(); ++it ) {
		(*it)->remap_resid( smap );
	}
}

void
EnzConstraintIO::set_position_for_missing_res_in_parameter_block(
	core::pose::Pose & pose,
	core::Size cst_block,
	core::Size respos
) const
{

	runtime_assert( cst_block <= cst_pairs_.size() );

	protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_block )->set_position_for_missing_res( respos );
}

void
EnzConstraintIO::clear_active_pose_constraints_for_block(
	core::pose::Pose & pose,
	core::Size cst_block
) const{
	runtime_assert( cst_block <= cst_pairs_.size() );
	protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache()->param_cache( cst_block)->clear_active_pose_constraints();
}


void
EnzConstraintIO::set_external_position_for_resA_in_parameter_block(
	core::Size cst_block,
	core::Size respos
)
{
	runtime_assert( cst_block <= cst_pairs_.size() );
	cst_pairs_[ cst_block ]->set_external_position_for_resA( respos );
}

void
EnzConstraintIO::set_external_position_for_resB_in_parameter_block(
	core::Size cst_block,
	core::Size respos
)
{
	runtime_assert( cst_block <= cst_pairs_.size() );
	cst_pairs_[ cst_block ]->set_external_position_for_resB( respos );
}

// migrate in second step
/*
void
EnzConstraintIO::setup_favor_native_constraints(
core::pose::Pose & pose,
core::pack::task::PackerTaskCOP task,
core::pose::Pose const & native_pose
)
{
using namespace basic::options;

if( option[OptionKeys::enzdes::favor_native_res].user() ) {
using namespace core::scoring::constraints;

core::Real bonus = option[OptionKeys::enzdes::favor_native_res].value();

tr.Info << "favor_native_res: adding a bonus of " << bonus << " for native residues to pose." << std::endl;

//safety check first
if( favor_native_constraints_.size() != 0 ){
tr.Info << "Warning: when setting up favor native constraints, there might already be some previously generated favor_native constraints in the pose, trying to remove these first." << std::endl;
remove_favor_native_constraints( pose );

}

favor_native_constraints_.clear();
for( core::Size i = 1; i <= pose.size(); ++i){

if( task->design_residue(i) ){

ConstraintOP resconstraint = new ResidueTypeConstraint( native_pose, i, bonus );
favor_native_constraints_.push_back( resconstraint );

}
}
favor_native_constraints_ = pose.add_constraints( favor_native_constraints_ );
} else if (option[ OptionKeys::in::file::pssm ].user() ) {
//multiple sequence aligniment (adapted from MSA app)

using namespace core;
using namespace scoring;
using namespace constraints;
using namespace sequence;

using namespace protocols;
using namespace constraints_additional;

tr << " Starting MSA design " << std::endl;

// register SequenceProfileConstraint with the ConstraintFactory so that it can be constructed from a constraint file
//ConstraintIO::get_cst_factory().add_type(
//new core::scoring::constraints::SequenceProfileConstraint( Size(), utility::vector1< id::AtomID >(), NULL ) );

// add constraints to bias design toward a sequence profile
SequenceProfileOP profile = new SequenceProfile;
utility::file::FileName filename( option[ OptionKeys::in::file::pssm ]().front() );

profile->read_from_file( filename );
profile->convert_profile_to_probs( 1 ); // was previously implicit in read_from_file()

tr << *profile << std::endl;

for ( Size seqpos(1), end( pose.size() ); seqpos <= end; ++seqpos ) {
// add individual profile constraint for each residue position
// because of the underlying constraint implementation, this enures that the constraint is
// a context-independent 1-body energy, or (intra)residue constraint
pose.add_constraint( new core::scoring::constraints::SequenceProfileConstraint( pose, seqpos, profile ) );
}
}// else if ( option[ OptionKeys::constraints::in::file::pssm ].user() ){
//}

}//setup_favor_native_constraints

//migrate in second step
void
EnzConstraintIO::remove_favor_native_constraints(
core::pose::Pose & pose
)
{
if( !( pose.remove_constraints( favor_native_constraints_ ) ) ){
tr.Info << "Warning: some of the favor native constraints that were previously added to the pose are not there anymore, something's a little unclean somewhere." << std::endl;
}
favor_native_constraints_.clear();
}
*/

EnzConstraintParametersCOP
EnzConstraintIO::enz_cst_params( core::Size block ) const
{
	return cst_pairs_[ block];
}

//add pose to interface
utility::vector1< EnzConstraintParametersCOP >
EnzConstraintIO::enz_cst_params_missing_in_pose( core::pose::Pose const & pose ) const
{
	utility::vector1< EnzConstraintParametersCOP > to_return;

	for ( Size i = 1; i <= cst_pairs_.size(); ++i ) {
		if ( cst_pairs_[i]->missing_in_pose( pose ) ) {
			to_return.push_back( cst_pairs_[i] );

			//runtime_assert( cst_pairs_[i]->active_pose_constraints().size() == 0 ); //old
		}
	}

	return to_return;

}

utility::vector1< core::Size >
EnzConstraintIO::ordered_constrained_positions( core::pose::Pose const & pose ) const
{
	using namespace core;
	//utility::vector1< Size > found_protein_positions;
	//utility::vector1< Size > found_lig_positions;
	protocols::toolbox::match_enzdes_util::EnzdesCacheableObserverCOP enz_obs( protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose ) );
	if ( !enz_obs ) return utility::vector1< Size > (); // get_enzdes_observer() can return NULL for const pose
	EnzdesCstCacheCOP cst_cache( enz_obs->cst_cache() );
	if ( !cst_cache ) return utility::vector1< Size > ();

	return cst_cache->ordered_constrained_positions( pose );

	/*
	for( core::Size i = 1; i <= cst_pairs_.size(); ++i ){

	for(std::map< Size, EnzCstTemplateResAtomsOP >::const_iterator resA_it = cst_cache->param_cache(i)->template_res_cache( 1 )->seqpos_map_begin(), resA_end = cst_cache->param_cache(i)->template_res_cache( 1 )->seqpos_map_end();
	resA_it !=  resA_end; ++resA_it ){
	if( pose.residue_type( resA_it->first ).is_ligand() ) {
	if( find( found_lig_positions.begin(), found_lig_positions.end(), resA_it->first ) == found_lig_positions.end() ) {
	found_lig_positions.push_back( resA_it->first );
	}
	}
	else found_protein_positions.push_back( resA_it->first );
	}

	for(std::map< Size, EnzCstTemplateResAtomsOP >::const_iterator resB_it = cst_cache->param_cache(i)->template_res_cache( 2 )->seqpos_map_begin(), resB_end = cst_cache->param_cache(i)->template_res_cache( 2 )->seqpos_map_end();
	resB_it !=  resB_end; ++resB_it ){
	if( pose.residue_type( resB_it->first ).is_ligand() ) {
	if( find( found_lig_positions.begin(), found_lig_positions.end(), resB_it->first ) == found_lig_positions.end() ) {
	found_lig_positions.push_back( resB_it->first );
	}
	}
	else found_protein_positions.push_back( resB_it->first );
	}

	}// loop over EnzConstraintParams

	for(  utility::vector1< Size >::const_iterator lig_it = found_lig_positions.begin(); lig_it != found_lig_positions.end(); ++lig_it ){
	found_protein_positions.push_back( *lig_it );
	}

	return found_protein_positions;
	*/
}//ordered catalytic positions

core::Size
EnzConstraintIO::mcfi_lists_size() const
{
	return mcfi_lists_.size();
}

void
EnzConstraintIO::determine_target_downstream_res()
{
	target_downstream_res_.clear();

	for ( core::Size i =1; i <= mcfi_lists_.size(); ++i ) {

		std::map< std::string, utility::vector1< std::string > > const &
			alg_info( (*mcfi_lists_[i]).mcfi(1)->algorithm_inputs() );

		if ( alg_info.find( "match" ) == alg_info.end() ) target_downstream_res_.push_back( std::pair<core::Size, core::Size >(1,1) );
		else {
			utility::vector1< std::string > const & info( alg_info.find( "match" )->second );
			std::pair<core::Size, core::Size > this_target(1,1);
			for ( core::Size ll = 1; ll <= info.size(); ++ll ) {
				std::string llstr = info[ ll ];
				std::istringstream llstream( llstr );
				std::string first, second;
				llstream >> first >> second;
				if ( first == "SECONDARY_MATCH:" && second == "UPSTREAM_CST" ) {
					this_target.second = 2;
					llstream >> this_target.first;
					break;
				}
			}
			target_downstream_res_.push_back( this_target );
		}
	}
}

}
} //enzdes
} //protocols



#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
protocols::toolbox::match_enzdes_util::EnzConstraintIO::EnzConstraintIO() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::EnzConstraintIO::save( Archive & arc ) const {
	arc( CEREAL_NVP( cst_pairs_ ) ); // utility::vector1<EnzConstraintParametersOP>
	arc( CEREAL_NVP( mcfi_lists_ ) ); // utility::vector1<toolbox::match_enzdes_util::MatchConstraintFileInfoListOP>
	arc( CEREAL_NVP( target_downstream_res_ ) ); // utility::vector1<std::pair<core::Size, core::Size> >
	core::chemical::serialize_residue_type_set( arc, restype_set_ ); // core::chemical::ResidueTypeSetCAP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::EnzConstraintIO::load( Archive & arc ) {
	arc( cst_pairs_ ); // utility::vector1<EnzConstraintParametersOP>
	arc( mcfi_lists_ ); // utility::vector1<toolbox::match_enzdes_util::MatchConstraintFileInfoListOP>
	arc( target_downstream_res_ ); // utility::vector1<std::pair<core::Size, core::Size> >
	core::chemical::deserialize_residue_type_set( arc, restype_set_ ); // core::chemical::ResidueTypeSetCAP
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::match_enzdes_util::EnzConstraintIO );
CEREAL_REGISTER_TYPE( protocols::toolbox::match_enzdes_util::EnzConstraintIO )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_match_enzdes_util_EnzConstraintIO )
#endif // SERIALIZATION
