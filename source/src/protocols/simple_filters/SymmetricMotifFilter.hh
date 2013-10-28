// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simplefilters/SymmetricMotifFilter.hh
/// @brief position-independent RMS filter evaluating how close a set of interfaces is to symmetric
/// @author Frank DiMaio

#ifndef INCLUDED_protocols_simple_filters_SymmetricMotifFilter_hh
#define INCLUDED_protocols_simple_filters_SymmetricMotifFilter_hh


#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <list>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>


namespace protocols {
namespace simple_filters {

struct Quat { core::Real x,y,z,w; };
void R2quat( numeric::xyzMatrix< core::Real > R, Quat &Q );
void quat2R( Quat &Q ,numeric::xyzMatrix< core::Real > R );
core::Real RMSwrapper( utility::vector1 <numeric::xyzVector< core::Real > > chainA,
                       utility::vector1 <numeric::xyzVector< core::Real > > chainB,
                       numeric::xyzMatrix< core::Real > &R, 
                      numeric::xyzVector< core::Real > &preT, numeric::xyzVector< core::Real > &postT);

//
class SymmetricMotifFilter : public protocols::filters::Filter
{
public:
	SymmetricMotifFilter();
	SymmetricMotifFilter(
		utility::vector1<core::pose::PoseOP> reference_motifs,
		std::string symm_type_in="D2");

	void set_defaults();

	bool apply( core::pose::Pose const & pose ) const;

	void add_motif( core::pose::PoseOP motif );
	void set_symm( std::string symm_type_in ) { symm_type_ = symm_type_in; };

	void set_thresholds( core::Real angle_thresh_in, core::Real trans_thresh_in, core::Real rmsd_thresh_in, core::Size /*clash_thresh_in*/ ) {
		angle_thresh_ = angle_thresh_in;
		trans_thresh_ = trans_thresh_in;
		rmsd_thresh_ = rmsd_thresh_in;
	}

	void set_weights( core::Real angle_thresh_in, core::Real trans_thresh_in, core::Real rmsd_thresh_in ) { 
		angle_thresh_ = angle_thresh_in;
		trans_thresh_ = trans_thresh_in;
		rmsd_thresh_ = rmsd_thresh_in;
	}

	// call after all motifs are added; verifies motifs are of correct symmetry & precomputes all transformations
	void process_motifs();

	protocols::filters::FilterOP clone() const;
	protocols::filters::FilterOP fresh_instance() const { return new SymmetricMotifFilter(); }

	void report( std::ostream & out, core::pose::Pose const & pose ) const;

	// compute angle and rms offsets
	bool compute( core::pose::Pose const & pose, core::Real &best_score, std::string &motifhit ) const;

	// symmetry-specific variants
	bool compute_d2( core::pose::Pose const & pose, core::Real &best_score, std::string &motifhit ) const;
	core::Real score_d2( core::Real rms, core::Real angle, core::Real trans, core::Size clash ) const {
		if (rms<=rmsd_thresh_ && angle<=angle_thresh_ && clash<=clash_thresh_) { 
			return(
				angle_wt_ * (angle-angle_thresh_)
				+ rmsd_wt_ * (rms-rmsd_thresh_)
				+ trans_wt_ * (trans-trans_thresh_)
				+ clash_wt_ * (clash-((int)clash_thresh_)) );
		} else {
			return(999);  // very large
		}
	}


	core::Real report_sm( core::pose::Pose const & pose ) const;

	virtual ~SymmetricMotifFilter();

	void parse_my_tag(
		utility::tag::TagCOP const tag,
		basic::datacache::DataMap & data_map,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & reference_pose );

private:
	// symmetry group and motifs
	std::string symm_type_;
	utility::vector1<core::pose::PoseOP> ref_motifs_;
	core::Size nsegs_;

	// filter parameters
	//  * cutoffs for pass/fail
	core::Real angle_thresh_, trans_thresh_,  rmsd_thresh_;
	core::Size clash_thresh_;
	//  * weights for report_sm
	core::Real angle_wt_, trans_wt_,  rmsd_wt_, clash_wt_;

	// force the motifs to a particular location
	utility::vector1< int > forced_pos_;

	// motif coords & transformations
	utility::vector1< Quat > Qs;
	utility::vector1< numeric::xyzMatrix< core::Real > > Rdimers;
	utility::vector1< numeric::xyzVector< core::Real > > delta_coms;
	utility::vector1< numeric::xyzVector< core::Real > > symm_axes;
	utility::vector1< core::Real  > symm_orders;
	utility::vector1< utility::vector1< numeric::xyzVector< core::Real > > > cas_chainA, cas_chainB;
	utility::vector1< utility::vector1< core::Size > > motif_cuts;
};

} // filters
}

#endif //INCLUDED_protocols_protein_interface_design_filters_SymmetricMotifFilter_HH_

