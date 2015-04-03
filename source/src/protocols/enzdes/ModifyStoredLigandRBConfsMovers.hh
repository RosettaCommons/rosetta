// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file .hh file for movers that mess around with additional ligand rigid body conformations
/// stored in the enzdes cacheable observer
/// @brief
/// @author Florian Richter, floric@u.washington.edu, oct 09

#ifndef INCLUDED_protocols_enzdes_ModifyStoredLigandRBConfsMovers_hh
#define INCLUDED_protocols_enzdes_ModifyStoredLigandRBConfsMovers_hh

//unit headers
#include <protocols/enzdes/ModifyStoredLigandRBConfsMovers.fwd.hh>
#include <protocols/moves/Mover.hh>
//#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.fwd.hh>
//#include <core/pose/datacache/CacheableObserver.hh>

//package headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace enzdes {


/// @pure virtual base class that has some common functionality
class ModifyStoredRBConfs : public protocols::moves::Mover {

public:
	typedef core::Size Size;
	typedef core::Real Real;
	typedef protocols::moves::Mover parent;
	typedef utility::vector1< utility::vector1< core::conformation::ResidueCOP > > RBConfLists;

public:

	ModifyStoredRBConfs( std::string const name );
	~ModifyStoredRBConfs();

	virtual
	void
	apply( core::pose::Pose & pose ) = 0;

	virtual std::string get_name() const;

	/// @brief helper function to exchange coordinates in a pose
	/// this function doesn't actually exchange the residues,
	/// but sets all the coords of pose.residue( rescoords.seqpos() )
	/// to what is in rescoords and vice versa
	void
	swap_coordinates_in_pose(
		core::pose::Pose & pose,
		core::conformation::Residue & rescoords
	) const;

	/// @brief calculates the closest orient atom (i.e. atoms used in rotamer placement
	/// msd between every conf and the confs coming before it in the pose or in the vector
	utility::vector1< Real >
	closest_orient_atoms_msd(
		core::pose::Pose const & pose,
		utility::vector1< core::conformation::ResidueCOP > const & confs
	) const;

protected:

	/// @brief this function is the function that all derived classes
	/// should call to get access to pose stored additional rb confs.
	/// right now it only returns the confs in the proper format,
	/// but in the future some logic might be added as to for which seqpos
	/// to return the stored confs
	RBConfLists
	get_rigid_body_confs( core::pose::Pose const & pose ) const;

	/// @brief same as for above function. all derived classes
	/// should use this function to communicate new/changed rigid
	/// body confs back to the pose/enzdes cacheable observer
	void
	set_rigid_body_confs( RBConfLists rbs, core::pose::Pose & pose ) const;

	/// @brief see comments for set_rigid_body_confs
	void
	set_rigid_body_confs_for_seqpos(
		core::Size seqpos,
		utility::vector1< core::conformation::ResidueCOP > & confs,
		core::pose::Pose & pose
	) const;

};

/// @brief generates random rbconfs until a total of num_total_rbconfs_
/// are present in the cacheable observer. The diversifier is used to
/// ensure that all newly generated confs are different.
/// note: no scorefunction used
class GenerateStoredRBConfs : public ModifyStoredRBConfs {

public:
	typedef ModifyStoredRBConfs parent;

public:

	GenerateStoredRBConfs(
		Size num_total_rbconfs,
		bool include_metals );

	~GenerateStoredRBConfs();

	void
	apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:

	Size num_total_rbconfs_;
	bool include_metals_;

};

/// @brief for every ligand that has additional stored rb conformations
/// in the enzdes cacheable observer, the one currently in the pose gets
/// exchanged with a random stored one
class ApplyRandomStoredRBConf : public ModifyStoredRBConfs {

public:
	typedef ModifyStoredRBConfs parent;

public:

	ApplyRandomStoredRBConf();
	~ApplyRandomStoredRBConf();

	void
	apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

};

/// @brief rotamer trials/minimizes each of the ligand conformations stored in the enzdes
/// cacheable observer and overwrites them with the minimized position.
/// if the new, minimized position is too close to one that already exists,
///  a random small perturbation is applied to ensure diversity
/// note: only jump minimization
class MinimizeStoredRBConfs : public ModifyStoredRBConfs {

public:
	typedef ModifyStoredRBConfs parent;

public:

	MinimizeStoredRBConfs( core::scoring::ScoreFunctionCOP sfxn );
	~MinimizeStoredRBConfs();

	void
	apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	void
	rb_minimize_all_confs(
		core::pose::Pose const & pose,
		utility::vector1< core::conformation::ResidueCOP > & confs
	) const;

private:

	core::scoring::ScoreFunctionCOP sfxn_;
	Real min_rms_;

};


/// @brief uses a docking mover to diversiy the stored confs
/// until they're all min_rms_ away from each other
/// note: no scorefunction used here
class DiversifyStoredRBConfs : public ModifyStoredRBConfs {

public:
	typedef ModifyStoredRBConfs parent;

public:

	DiversifyStoredRBConfs(
		Real min_rms );

	~DiversifyStoredRBConfs();

	void
	apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	void
	diversify_all_confs(
		core::pose::Pose const & pose,
		utility::vector1< core::conformation::ResidueCOP > & confs
	) const;

private:
	Real min_rms_;
	Size max_trials_;

};


} // enzdes
} //protocols


#endif
