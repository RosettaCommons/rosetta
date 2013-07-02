// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/MembraneTopology.hh
/// @brief  MembraneTopology
/// @author Bjorn Wallner
/*
This code reads in the spanfile and lipo file that is provided by the user.
Residues in the membrane are stored in private data which is accessed by
cacheable data from the pose.
*/

#ifndef INCLUDED_core_scoring_MembraneTopology_hh
#define INCLUDED_core_scoring_MembraneTopology_hh

#include <core/types.hh>

// Unit headers
#include <core/scoring/MembraneTopology.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/CacheableData.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>


namespace core {
namespace scoring {



class MembraneTopology : public basic::datacache::CacheableData {

public:
	MembraneTopology(): LipoDefined_(false),init_(false),beta_barrel_(false),N_term_inside_(false),initialized_(false) {};

	MembraneTopology( MembraneTopology const & src);

	basic::datacache::CacheableDataOP
	clone() const
	{
		return new MembraneTopology( *this );
	}

//pba
  bool
  initialized() const
  {
    return initialized_;
  }

  bool &
  initialized()
  {
    return initialized_;
  }

	Size
	tmhelix() const {
		return total_tmhelix_;
	}

	Size
	span_begin(Size n) const {
		return span_(n,1);
	}

	Size
	span_end(Size n) const {
		return span_(n,2);
	}

	Size
	helix_id(Size n) const {
		return helix_id_(n);
	}

	Size
	relative_tmh_ori(Size const n, Size const m) const {
		return(relative_tmh_ori(n,m));
	}

	Real
	depth(Size const seqpos) {
		return depth_[ seqpos ];
	}

	utility::vector1< core::Real > &
	depth() {
		return depth_;
	}
	void
	initialize(std::string const & spanfile);

	//void
	//initialize(pose::Pose & pose, std::string const & spanfile);
	///

	void
	print() const;

	void
	shift_span(Size shift);

 	//void attach_to_pose(pose::Pose & pose);


	void
	get_subset( utility::vector1< Size > & TMH_list , MembraneTopology & src);

	bool
	allow_scoring(Size const seqpos) const
	{
		if(seqpos > allow_scoring_.size())
		{
			throw utility::excn::EXCN_RangeError(
					"tried to get Membrane score for residue " +
					utility::to_string(seqpos) +
					" but spanfile only specifies " +
					utility::to_string(allow_scoring_.size())+
					" residues. Check your spanfile.");
			return false; //so the compiler is ok with it
		}
		return allow_scoring_[seqpos];
	}
/*
	bool &
	allow_scoring(Size const seqpos)
	{
		return allow_scoring_[seqpos];
	}
*/
	Size
	tmh_inserted() const
	{
		return tmh_inserted_;
	}

	void
	reset_tmh_insert()
	{
		tmh_inserted_=0;
	}
	void
	reset_allowed_scoring()
	{
		tmh_inserted_=0;
		for(Size i=1;i<=allow_tmh_scoring_.size();++i)
		{
			allow_tmh_scoring_[i]=false;
		}
		for(Size i=1;i<=allow_scoring_.size();++i)
		{
			allow_scoring_[i]=false;
		}

	}
	void
	set_tmh_inserted(Size tmh_inserted)
	{
		tmh_inserted_=tmh_inserted;
	}

	bool
	tmregion(Size const pos) const
	{
		return tmregion_[pos];
	}

	bool
	allow_tmh_scoring(Size const tmh) const
	{
		return allow_tmh_scoring_[tmh];
	}

	void
	set_allow_tmh_scoring(Size const tmh,bool setting)
	{
		allow_tmh_scoring_[tmh]=setting;
	}
	void
	set_allow_scoring(Size const pos, bool setting)
	{
		allow_scoring_[pos]=setting;
	}

	Real
	LipidExposure(Size const n) const {
		return LipidExposure_[n];
	}
	Real
	LipidBurial(Size const n) const {
		return LipidBurial_[n];
	}

	bool
	LipoDefined() const {
		return LipoDefined_;
	}
	/*
	bool &
	allow_tmh_scoring(Size const tmh)
	{
		return allow_tmh_scoring_[tmh];
	}
	bool &
	allow_scoring(Size const pos)
	{
		return allow_scoring_[pos];
	}
	 */
 protected:



/*	Real const cen_dist_cutoff2;


	CenListInfo const & cenlist_from_pose( pose::Pose const & ) const;
	CenListInfo & nonconst_cenlist_from_pose( pose::Pose & ) const;
*/
private:


private: // data

	ObjexxFCL::FArray1D< Size > helix_id_;
	ObjexxFCL::FArray2D< Size > span_;
	ObjexxFCL::FArray2D< Size > full_span_;
	ObjexxFCL::FArray2D< Size > relative_tmh_ori_;
	Size total_tmhelix_;
	utility::vector1< core::Real > depth_; // this is just to speed up the cross talk between scoring fxns and pose.metrics.
	utility::vector1< core::Real > LipidExposure_;
	utility::vector1< core::Real > LipidBurial_;
	bool LipoDefined_;
	bool init_;
	bool beta_barrel_;
	bool N_term_inside_;
  bool initialized_; //pba
	utility::vector1< bool > tmregion_; //stores if the residue is in the TM or not
	utility::vector1< bool > allow_scoring_;
	utility::vector1< bool > allow_tmh_scoring_;
	Size tmh_inserted_;


};

//extern MembraneTopology &  MembraneTopology_from_pose(core::pose:Pose const pose);
//pbadebug
MembraneTopology const & MembraneTopology_from_pose( pose::Pose const & pose );
MembraneTopology & nonconst_MembraneTopology_from_pose( pose::Pose & pose );


} // ns scoring
} // ns core

#endif
