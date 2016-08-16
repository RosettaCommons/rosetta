// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file .hh file for enzdes cacheable observer
/// @brief
/// @author Florian Richter, floric@u.washington.edu

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_EnzdesCstCache_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_EnzdesCstCache_hh

//unit headers
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.fwd.hh>

//package headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>

#ifdef WIN32
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintParameters.hh>
#endif

//project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {


/// @brief a simple class to store the pose specific
/// enzdes constraint information in the actual pose
class EnzdesCstCache : public utility::pointer::ReferenceCount {

public: //constructor/destructor/

	EnzdesCstCache(
		EnzConstraintIOCOP enz_cstio,
		core::Size num_cst_blocks
	);

	EnzdesCstCache( EnzdesCstCache const & other );

	virtual ~EnzdesCstCache();

public: //accessors

	EnzConstraintIOCOP
	enzcst_io() const;

	EnzdesCstParamCacheOP
	param_cache(
		core::Size cst_block );

	EnzdesCstParamCacheCOP
	param_cache(
		core::Size cst_block ) const;

	core::Size
	ncsts() const { return param_cache_.size(); }

	void
	set_param_cache(
		core::Size cst_block,
		EnzdesCstParamCacheOP param_cache
	);

	bool
	contains_position(
		core::Size seqpos) const;

	/// @brief returns all seqpos present in the cache in order of the params
	/// which they are in.
	///note: ligands are explicityly put at the end of the vector
	utility::vector1< core::Size >
	ordered_constrained_positions( core::pose::Pose const & pose ) const;

	/// @brief remapping sequence positions
	void
	remap_resid( core::id::SequenceMapping const & smap );

private:

	utility::vector1< EnzdesCstParamCacheOP > param_cache_;
	EnzConstraintIOCOP enzcst_io_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	EnzdesCstCache();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


class EnzdesCstParamCache : public utility::pointer::ReferenceCount {

public: //typedefs

	friend class EnzConstraintParameters;

public: //constructor/destructor

	EnzdesCstParamCache();

	EnzdesCstParamCache( EnzdesCstParamCache const & other );

	virtual ~EnzdesCstParamCache();

public: //accessors

	EnzCstTemplateResCacheOP
	template_res_cache( core::Size index ) {
		return template_res_cache_[index]; }

	EnzCstTemplateResCacheCOP
	template_res_cache( core::Size index ) const {
		return template_res_cache_[index]; }

	core::Size
	template_res_cache_size() const {
		return template_res_cache_.size(); }

	utility::vector1< core::scoring::constraints::ConstraintCOP > const &
	active_pose_constraints() const{
		return active_pose_constraints_; }

	utility::vector1< core::scoring::constraints::ConstraintCOP > &
	active_pose_constraints() {
		return active_pose_constraints_; }

	void
	set_active_pose_constraints(
		utility::vector1< core::scoring::constraints::ConstraintCOP > const & constraints
	);

	void
	clear_active_pose_constraints();

	//utility::vector1< CovalentConnectionReplaceInfo > &
	//covalent_connections();

	utility::vector1< CovalentConnectionReplaceInfoCOP > const &
	covalent_connections() const{
		return covalent_connections_; }

	bool
	missing_in_pose() const;

	void
	remove_seqpos_from_template_res( core::Size seqpos );

	void
	set_position_for_missing_res( core::Size seqpos );

	bool
	contains_position( core::Size seqpos ) const;

	/// @brief remapping sequence positions
	void
	remap_resid( core::id::SequenceMapping const & smap );

private:

	utility::vector1< EnzCstTemplateResCacheOP > template_res_cache_;
	utility::vector1< core::scoring::constraints::ConstraintCOP > active_pose_constraints_;
	utility::vector1< CovalentConnectionReplaceInfoCOP > covalent_connections_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class EnzCstTemplateResCache : public utility::pointer::ReferenceCount {

public: //typedefs

	friend class EnzCstTemplateRes;

	typedef std::map< core::Size, EnzCstTemplateResAtomsOP > SeqposTemplateAtomsMap;

public: //constructor/destructor

	EnzCstTemplateResCache();

	EnzCstTemplateResCache( EnzCstTemplateResCache const & other );

	virtual ~EnzCstTemplateResCache();

public:  //accessors

	std::map< Size, EnzCstTemplateResAtomsOP >::const_iterator
	seqpos_map_begin() const{
		return seqpos_map_.begin();
	}

	std::map< Size, EnzCstTemplateResAtomsOP >::const_iterator
	seqpos_map_end() const{
		return seqpos_map_.end();
	}

	core::Size
	seqpos_map_size() const{
		return seqpos_map_.size(); }

	bool
	contains_position( core::Size const seqpos ) const {
		return seqpos_map_.find( seqpos ) != seqpos_map_.end();
	}

	void
	set_position_in_pose( core::Size seqpos );

	void
	add_position_in_pose( core::Size seqpos );

	bool
	remove_seqpos( core::Size seqpos );

	/// @brief remapping sequence positions
	void
	remap_resid( core::id::SequenceMapping const & smap );

	bool
	not_in_pose() const {
		return not_in_pose_; }

private:

	SeqposTemplateAtomsMap seqpos_map_;
	bool not_in_pose_;
	bool pose_data_uptodate_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}
} //toolbox
} //protocols


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_match_enzdes_util_EnzdesCstCache )
#endif // SERIALIZATION


#endif
