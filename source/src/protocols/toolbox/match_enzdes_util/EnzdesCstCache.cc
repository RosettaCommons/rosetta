// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file .cc file for enzdes cacheable observer
/// @brief
/// @author Florian Richter, floric@u.washington.edu, september 09


//unit headers
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>

//package headers
#include <protocols/toolbox/match_enzdes_util/EnzCstTemplateRes.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintParameters.hh> //for CovalentConnectionReplaceInfo

//project headers
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/id/SequenceMapping.hh>
#include <core/scoring/constraints/Constraint.hh>


//utility headers
#include <utility/string_util.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

EnzdesCstCache::EnzdesCstCache(
	EnzConstraintIOCOP enzcst_io,
	core::Size num_cst_blocks
) : enzcst_io_(enzcst_io)
{
	param_cache_.clear();
	for ( core::Size i = 1; i <= num_cst_blocks; ++i ) {
		param_cache_.push_back( EnzdesCstParamCacheOP( new EnzdesCstParamCache() ) );
	}
}

EnzdesCstCache::EnzdesCstCache( EnzdesCstCache const & other )
: ReferenceCount( other ), enzcst_io_(other.enzcst_io_)
{
	param_cache_.clear();
	for ( utility::vector1< EnzdesCstParamCacheOP >::const_iterator param_it( other.param_cache_.begin()), param_end(other.param_cache_.end()); param_it != param_end; ++param_it ) {
		param_cache_.push_back( EnzdesCstParamCacheOP( new EnzdesCstParamCache( **param_it ) ) );
	}
}

EnzdesCstCache::~EnzdesCstCache() {}

EnzConstraintIOCOP
EnzdesCstCache::enzcst_io() const
{
	return enzcst_io_;
}

EnzdesCstParamCacheOP
EnzdesCstCache::param_cache(
	core::Size cst_block )
{
	runtime_assert( cst_block <= param_cache_.size() );
	return param_cache_[ cst_block ];
}

EnzdesCstParamCacheCOP
EnzdesCstCache::param_cache(
	core::Size cst_block ) const
{
	if ( cst_block > param_cache_.size() ) {
		std::cerr << "param_cache_.size:" << param_cache_.size() << " cst_block:" << cst_block << std::endl;
		std::cerr << "cst_block should be smaller or equal to param_cache_.size()" << std::endl;
		runtime_assert( cst_block <= param_cache_.size() );
	}
	return param_cache_[ cst_block ];
}

void
EnzdesCstCache::set_param_cache(
	core::Size cst_block,
	EnzdesCstParamCacheOP param_cache )
{
	runtime_assert( cst_block <= param_cache_.size() );
	param_cache_[ cst_block ] = param_cache;
}

bool
EnzdesCstCache::contains_position(
	core::Size seqpos) const
{
	for ( core::Size ii = 1; ii <= this->ncsts(); ++ii ) {
		if ( this->param_cache( ii )->contains_position( seqpos ) ) return true;
	}
	return false;
}

utility::vector1< core::Size >
EnzdesCstCache::ordered_constrained_positions( core::pose::Pose const & pose ) const
{

	using namespace core;
	utility::vector1< Size > found_protein_positions;
	utility::vector1< Size > found_lig_positions;

	for ( Size i =1; i<= param_cache_.size(); ++i ) {

		for ( std::map< Size, EnzCstTemplateResAtomsOP >::const_iterator resA_it = param_cache_[i]->template_res_cache( 1 )->seqpos_map_begin(), resA_end = param_cache_[i]->template_res_cache( 1 )->seqpos_map_end(); resA_it !=  resA_end; ++resA_it ) {
			if ( pose.residue_type( resA_it->first ).is_ligand() ) {
				if ( find( found_lig_positions.begin(), found_lig_positions.end(), resA_it->first ) == found_lig_positions.end() ) {
					found_lig_positions.push_back( resA_it->first );
				}
			} else found_protein_positions.push_back( resA_it->first );
		}

		for ( std::map< Size, EnzCstTemplateResAtomsOP >::const_iterator resB_it = param_cache_[i]->template_res_cache( 2 )->seqpos_map_begin(), resB_end = param_cache_[i]->template_res_cache( 2 )->seqpos_map_end(); resB_it !=  resB_end; ++resB_it ) {
			if ( pose.residue_type( resB_it->first ).is_ligand() ) {
				if ( find( found_lig_positions.begin(), found_lig_positions.end(), resB_it->first ) == found_lig_positions.end() ) {
					found_lig_positions.push_back( resB_it->first );
				}
			} else found_protein_positions.push_back( resB_it->first );
		}
	} //loop over params

	for (  utility::vector1< Size >::const_iterator lig_it = found_lig_positions.begin(); lig_it != found_lig_positions.end(); ++lig_it ) {
		found_protein_positions.push_back( *lig_it );
	}
	return found_protein_positions;
} //ordered constrained positions

void
EnzdesCstCache::remap_resid(
	core::id::SequenceMapping const & smap )
{
	for ( utility::vector1< EnzdesCstParamCacheOP >::iterator param_it = param_cache_.begin(), param_end = param_cache_.end(); param_it != param_end; ++param_it ) {
		(*param_it)->remap_resid( smap );
	}
}

/// @brief Default constructor
/// right now two instances of cached template res get created
/// probably not ideal in general case, but is always the case
/// in the current implementation of enzdes, and this is the only
/// thing using this at the moment
EnzdesCstParamCache::EnzdesCstParamCache()
{
	template_res_cache_.push_back( EnzCstTemplateResCacheOP( new EnzCstTemplateResCache() ) );
	template_res_cache_.push_back( EnzCstTemplateResCacheOP( new EnzCstTemplateResCache() ) );
}

EnzdesCstParamCache::EnzdesCstParamCache( EnzdesCstParamCache const & other )
: ReferenceCount( other ),
	active_pose_constraints_( other.active_pose_constraints_ ),
	covalent_connections_(other.covalent_connections_ )
{
	template_res_cache_.clear();
	for ( core::Size i = 1; i <= other.template_res_cache_.size(); ++i ) {
		template_res_cache_.push_back( EnzCstTemplateResCacheOP( new EnzCstTemplateResCache( *(other.template_res_cache_[i] ) ) ) );
	}
}


EnzdesCstParamCache::~EnzdesCstParamCache(){}


void
EnzdesCstParamCache::set_active_pose_constraints(
	utility::vector1< core::scoring::constraints::ConstraintCOP > const & constraints
){
	active_pose_constraints_ = constraints;
}

void
EnzdesCstParamCache::clear_active_pose_constraints()
{
	active_pose_constraints_.clear();
}

bool
EnzdesCstParamCache::missing_in_pose() const
{
	for ( utility::vector1< EnzCstTemplateResCacheOP >::const_iterator template_res_it = template_res_cache_.begin(), template_res_end = template_res_cache_.end(); template_res_it != template_res_end; ++template_res_it ) {
		if ( (*template_res_it)->not_in_pose() ) return true;
	}
	return false;
}

void
EnzdesCstParamCache::remove_seqpos_from_template_res( core::Size seqpos ){

	core::Size find_count(0);
	for ( utility::vector1< EnzCstTemplateResCacheOP >::iterator template_res_it = template_res_cache_.begin(), template_res_end = template_res_cache_.end(); template_res_it != template_res_end; ++template_res_it ) {
		if ( (*template_res_it)->remove_seqpos( seqpos ) ) find_count++;
	}

	if ( find_count > 1 ) {
		utility_exit_with_message("Error: Several template residues were apparently at pose position "+utility::to_string(seqpos )+". This shouldn't happen...\n");
	}

	if ( find_count == 0 ) {
		utility_exit_with_message("Error: No template residues were found at pose position "+utility::to_string(seqpos )+". Something's unclean somewhere...\n");
	}
} //remove_seqpos_from_template_res

void
EnzdesCstParamCache::set_position_for_missing_res( core::Size seqpos ){

	core::Size find_count(0);
	for ( utility::vector1< EnzCstTemplateResCacheOP >::iterator template_res_it = template_res_cache_.begin(), template_res_end = template_res_cache_.end(); template_res_it != template_res_end; ++template_res_it ) {

		if ( (*template_res_it)->not_in_pose() ) {
			(*template_res_it)->set_position_in_pose( seqpos );
			find_count++;
		}
	}
	if ( find_count > 1 ) utility_exit_with_message("Error: Several template residues are missing in the pose. This shouldn't happen...\n");
	if ( find_count == 0 ) utility_exit_with_message("Error: no template residue is missing in the pose, this shouldn't have happened... \n");
}


bool
EnzdesCstParamCache::contains_position( core::Size seqpos ) const
{
	for ( utility::vector1< EnzCstTemplateResCacheOP >::const_iterator template_res_it = template_res_cache_.begin(), template_res_end = template_res_cache_.end(); template_res_it != template_res_end; ++template_res_it ) {
		if ( (*template_res_it)->contains_position( seqpos ) ) return true;
	}
	return false;
}


/// @details have to remap all the data that's cached
/// 1. the information about the template res
/// 2. the actual constraints
/// 3. the information about covalent connections
void
EnzdesCstParamCache::remap_resid( core::id::SequenceMapping const & smap )
{
	// 1.
	for ( utility::vector1< EnzCstTemplateResCacheOP >::iterator template_res_it = template_res_cache_.begin(), template_res_end = template_res_cache_.end(); template_res_it != template_res_end; ++template_res_it ) {
		(*template_res_it)->remap_resid( smap );
	}

	//2.
	for ( utility::vector1< core::scoring::constraints::ConstraintCOP >::iterator cst_it = active_pose_constraints_.begin();
			cst_it != active_pose_constraints_.end(); ++cst_it ) {
		*cst_it = (*cst_it)->remap_resid( smap );
		if ( ! (*cst_it) ) utility_exit_with_message("Remapping of catalytic constraints failed");
	}

	//3.
	for ( utility::vector1< CovalentConnectionReplaceInfoCOP >::iterator con_it = covalent_connections_.begin(), con_end = covalent_connections_.end(); con_it != con_end; ++con_it ) {
		CovalentConnectionReplaceInfoOP newcov( new CovalentConnectionReplaceInfo( **con_it ) );
		newcov->remap_resid( smap );
		(*con_it) = newcov;
	}
}


EnzCstTemplateResCache::EnzCstTemplateResCache() :
	ReferenceCount(),
	not_in_pose_(true), pose_data_uptodate_(false)
{
	seqpos_map_.clear();
}

EnzCstTemplateResCache::EnzCstTemplateResCache( EnzCstTemplateResCache const & other ) :
	ReferenceCount( other ),
	not_in_pose_(other.not_in_pose_),
	pose_data_uptodate_(other.pose_data_uptodate_)
{
	seqpos_map_.clear();
	for ( SeqposTemplateAtomsMap::const_iterator map_it(other.seqpos_map_.begin()), map_end(other.seqpos_map_.end()); map_it != map_end; ++map_it ) {
		seqpos_map_.insert( std::pair< core::Size, EnzCstTemplateResAtomsOP >( map_it->first, EnzCstTemplateResAtomsOP( new EnzCstTemplateResAtoms( *(map_it->second) ) ) ) );
	}
}

EnzCstTemplateResCache::~EnzCstTemplateResCache(){}

void
EnzCstTemplateResCache::set_position_in_pose( core::Size seqpos ){
	seqpos_map_.clear();
	seqpos_map_.insert( std::pair<core::Size,EnzCstTemplateResAtomsOP>(seqpos, EnzCstTemplateResAtomsOP( new EnzCstTemplateResAtoms() ) ) );
	not_in_pose_ = false;
	pose_data_uptodate_ = false;
}

void
EnzCstTemplateResCache::add_position_in_pose( core::Size seqpos ){
	seqpos_map_.insert( std::pair<core::Size,EnzCstTemplateResAtomsOP>(seqpos, EnzCstTemplateResAtomsOP( new EnzCstTemplateResAtoms() ) ) );
	not_in_pose_ = false;
	pose_data_uptodate_ = false;
}

bool
EnzCstTemplateResCache::remove_seqpos( core::Size seqpos ){

	bool return_val = seqpos_map_.erase( seqpos );
	if ( seqpos_map_.size() == 0  ) {
		not_in_pose_ = true;
		pose_data_uptodate_ = false;
	}
	return return_val;
}

void
EnzCstTemplateResCache::remap_resid( core::id::SequenceMapping const & smap ){

	// we have to remap the maps, little complicated to change the indexes of a map
	utility::vector1< std::pair< Size, EnzCstTemplateResAtomsOP > > temp_vec;
	for ( SeqposTemplateAtomsMap::const_iterator map_it(seqpos_map_.begin()), map_end(seqpos_map_.end()); map_it != map_end; ++map_it ) {
		temp_vec.push_back( *map_it );
	}
	seqpos_map_.clear();

	for ( utility::vector1< std::pair< Size, EnzCstTemplateResAtomsOP > >::iterator vec_it = temp_vec.begin(), vec_end = temp_vec.end(); vec_it != vec_end; ++vec_it ) {
		vec_it->second->remap_resid( smap );
		core::Size newpos( smap[ vec_it->first  ] );
		if ( newpos == 0 ) utility_exit_with_message("A catalytic residue is apparently missing from the pose");
		seqpos_map_.insert( std::pair< Size, EnzCstTemplateResAtomsOP >( newpos, vec_it->second ) );
	}
}

}
} // enzdes
} //protocols
