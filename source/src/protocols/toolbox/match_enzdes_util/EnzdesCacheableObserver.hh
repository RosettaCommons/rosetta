// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file .hh file for enzdes cacheable observer
/// @brief
/// @author Florian Richter, floric@u.washington.edu

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_EnzdesCacheableObserver_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_EnzdesCacheableObserver_hh

//unit headers
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.fwd.hh>
#include <core/pose/datacache/CacheableObserver.hh>

//package headers
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesSeqRecoveryCache.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesLoopsFile.fwd.hh>

//project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>

//utility headers
#include <utility/signals/Link.hh>

#include <map>

#include <core/scoring/constraints/Constraint.fwd.hh>
#include <utility/vector1.hh>

#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh>
#endif


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {


/// @brief convenience function to get a cacheable observer from a pose
EnzdesCacheableObserverOP
get_enzdes_observer(
	core::pose::Pose & pose );

EnzdesCacheableObserverCOP
get_enzdes_observer(
	core::pose::Pose const & pose );

class EnzdesCacheableObserver : public core::pose::datacache::CacheableObserver {

public: //typedefs

	//typedef utility::signals::Link Link;

public: //constructor/destructor/

	EnzdesCacheableObserver();

	EnzdesCacheableObserver( EnzdesCacheableObserver const & other );

	~EnzdesCacheableObserver();

	virtual
	core::pose::datacache::CacheableObserverOP
	clone();

	virtual
	core::pose::datacache::CacheableObserverOP
	create();

public: //observer interface

	virtual
	bool is_attached() const;

protected: //observer interface

	virtual
	void attach_impl( core::pose::Pose & pose );

	virtual
	void detach_impl();

	void
	on_length_change( core::conformation::signals::LengthEvent const & event );

public: //enzdes specific stuff

	void
	set_cst_cache(
		toolbox::match_enzdes_util::EnzdesCstCacheOP cst_cache
	);

	toolbox::match_enzdes_util::EnzdesCstCacheOP
	cst_cache();

	toolbox::match_enzdes_util::EnzdesCstCacheCOP
	cst_cache() const;

	std::map< core::Size, utility::vector1< core::conformation::ResidueCOP > > const &
	lig_rigid_body_confs() const;

	void
	set_rigid_body_confs_for_lig(
		core::Size seqpos,
		utility::vector1< core::conformation::ResidueCOP > const & rg_confs
	);

	void
	erase_rigid_body_confs_for_lig(
		core::Size seqpos );

	void set_seq_recovery_cache( EnzdesSeqRecoveryCacheOP seq_recovery_cache);

	EnzdesSeqRecoveryCacheOP get_seq_recovery_cache();
	EnzdesSeqRecoveryCacheCOP get_seq_recovery_cache() const;

	void
	set_enzdes_loops_file(
		EnzdesLoopsFileCOP loopfile_in
	);

	EnzdesLoopsFileCOP
	enzdes_loops_file() const;

	void
	setup_favor_native_constraints(
		core::pose::Pose & pose,
		core::pack::task::PackerTaskCOP task,
		core::pose::Pose const & native_pose
	);

	void
	remove_favor_native_constraints(
		core::pose::Pose & pose
	);


private:

	//enzdes constraint storage
	toolbox::match_enzdes_util::EnzdesCstCacheOP cst_cache_;
	//contains wt seq and keeps track of the designed residues
	EnzdesSeqRecoveryCacheOP seq_recovery_cache_;

	//favor native constraints
	core::scoring::constraints::ConstraintCOPs favor_native_constraints_;

	//evtl loops file
	EnzdesLoopsFileCOP enz_loops_file_;

	// storage for different rigid body conformations for ligands
	//e.g. from matching or docking
	std::map< core::Size, utility::vector1< core::conformation::ResidueCOP > > lig_rigid_body_confs_;

	utility::signals::Link length_event_link_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //match_enzdes_util
} //toolbox
} //protocols


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_match_enzdes_util_EnzdesCacheableObserver )
#endif // SERIALIZATION


#endif
