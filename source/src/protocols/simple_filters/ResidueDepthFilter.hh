// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/simple_filter/SidechainDepthFilter.hh
/// @brief Filter for measuring minimum distance from any sidechain atom in a reside
//  to the closest water molecule
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_simple_filters_ResidueDepthFilter_hh
#define INCLUDED_protocols_simple_filters_ResidueDepthFilter_hh

#include <protocols/simple_filters/ResidueDepthFilter.fwd.hh>

#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>
#include <numeric/xyzVector.hh>

#ifdef WIN32
#include <utility/tag/Tag.hh>
#endif

#include <utility/vector1.hh>
#include <math.h>

namespace protocols {
namespace simple_filters {

using namespace core;

///////////////////////////////////////////////////////////////////////////
// Data structures
struct ResidueDepthData
{
	Vector CAcrd;
	Vector CENcrd;
	core::Size ipdb, ires, n8;
	std::string aa;
	core::Real depth;
	utility::vector1< core::Size > neigh_ress;
};

class ResidueDepthFrag : public utility::pointer::ReferenceCount
{
public:
	ResidueDepthFrag();
	~ResidueDepthFrag();
	utility::vector1< Vector > get_CENcrd() const;
	utility::vector1< Vector > get_CAcrd() const;

	Vector get_CENcrd( core::Size const i ) const { return rdds_[i]->CENcrd; }
	Vector get_CAcrd( core::Size const i ) const { return rdds_[i]->CAcrd; }

	std::string aa() const { return rdds_[5]->aa; }
	core::Size ipdb() const { return rdds_[5]->ipdb; }
	core::Size ires() const { return rdds_[5]->ires; }

	void set_rdd( ResidueDepthDataCOP rdr, core::Size const i ){ rdds_[i] = ResidueDepthDataCOP( rdr ); }

private:
	utility::vector1< ResidueDepthDataCOP > rdds_;
};

///////////////////////////////////////////////////////////////////////////
// Depth calculator
class ResidueDepthCalculator
{
public:
	ResidueDepthCalculator(){}
	ResidueDepthCalculator( core::pose::Pose const &pose );
	~ResidueDepthCalculator();

	// this is main
	utility::vector1< core::Real >
	estimate_sidechain_depth( core::pose::Pose const &pose ) const;

	// setters
	void report_crd( bool const value ){ report_crd_ = value; }
	void niter( core::Size const value ){ niter_ = value; }
	void set_dcut1( core::Real const value ){ dcut1_ = value; }
	void set_dcut2( core::Real const value ){ dcut2_ = value; }

	// helpers
	core::Real get_scdepth_avrg( core::Size const ires ) const
	{ if ( ires > sc_depth_avrg_.size() ) return -1.0;
		return sc_depth_avrg_[ires]; }

	Size nres() const { return nres_; };
	utility::vector1< core::Real > get_scdepth_avrg() const { return sc_depth_avrg_; }
	utility::vector1< core::Real > get_scdepth_sdev() const { return sc_depth_sdev_; }
	utility::vector1< core::Real > get_scdepth_fvar() const { return sc_depth_fvar_; }

private:

	void
	initialize( core::pose::Pose const &pose );

	utility::vector1< Vector >
	read_unit_waterbox( Vector &boxwidth) const;

	void
	get_pose_crd_and_index( core::pose::Pose const &pose,
		utility::vector1< Vector > &protein_crd,
		utility::vector1< core::Size > &coarse_index,
		utility::vector1< core::Size > &res_id
	) const;

	utility::vector1< Vector >
	get_coarse_crd( utility::vector1< Vector > const &protein_crd,
		utility::vector1< core::Size > const &coarse_index
	) const;

	void
	duplicate_waterbox( utility::vector1< Vector > const &unit_waterbox_crd,
		Vector const boxwidth,
		Vector const mincrds, Vector const maxcrds ) const;
	Vector
	bring_to_origin( utility::vector1< Vector > &protein_crd,
		Vector &maxcrds,
		Vector &mincrds
	) const;
	void
	append_unitbox( utility::vector1< Vector > const &unitbox,
		utility::vector1< Vector > &waterbox_dupl,
		Vector const boxwidth,
		int const i, int const j, int const k ) const;

	void
	pert_protein( utility::vector1< Vector > &protein_crd ) const;

	utility::vector1< Vector >
	quat2U( Vector const & quat3, core::Real const &quatw ) const;

	utility::vector1< bool >
	get_exclusion_index( utility::vector1< Vector > const & protein_crd,
		utility::vector1< Vector > const & protein_coarse_crd ) const;

	utility::vector1< core::Real >
	get_scdepth( utility::vector1< bool > const &excluded_wat,
		utility::vector1< Vector > const & protein_crd,
		utility::vector1< core::Size > const & res_id
	) const;

	bool
	stack_and_getaverage(
		utility::vector1< utility::vector1< core::Real > > &sc_depth_stack,
		utility::vector1< core::Real > const &sc_depth,
		core::Size const niter
	) const;

private:

	core::Size niter_;
	core::Size nres_;
	std::string waterbox_file_;
	Real dcut1_, dcut2_;

	bool use_bb_, use_sc_;
	bool report_crd_;

	mutable utility::vector1< Vector > waterbox_;
	mutable utility::vector1< core::Real > sc_depth_avrg_;
	mutable utility::vector1< core::Real > sc_depth_sdev_;
	mutable utility::vector1< core::Real > sc_depth_fvar_;

}; // ResidueDepthCalculator

///////////////////////////////////////////////////////////////////////////
class ResidueDepthFilter : public filters::Filter
{
public:
	ResidueDepthFilter() : filters::Filter( "ResidueDepth" ) {}
	~ResidueDepthFilter();

	//ResidueDepthFilter( ResidueDepthFilter const &init );
	ResidueDepthFilter( core::pose::Pose const &pose );

	bool apply( core::pose::Pose const & ) const;

	filters::FilterOP clone() const {
		return filters::FilterOP( new ResidueDepthFilter( *this ) );
	}

	filters::FilterOP fresh_instance() const{
		return filters::FilterOP( new ResidueDepthFilter() );
	}

	// this is the main
	core::Real get_SDE_score( core::pose::Pose const &pose );

	// wrapper for residue depth calculation
	utility::vector1< core::Real >
	get_residue_depth( core::pose::Pose const &pose ) const;

	void
	parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

	// setters
	void set_report_crd( bool const value ){ RDC_.report_crd( value ); }
	void set_niter( core::Size const value ){ RDC_.niter( value ); }

	// helpers
	utility::vector1< core::Real > get_scdepth_avrg() const { return RDC_.get_scdepth_avrg(); }
	utility::vector1< core::Real > get_scdepth_sdev() const { return RDC_.get_scdepth_sdev(); }
	utility::vector1< core::Real > get_scdepth_fvar() const { return RDC_.get_scdepth_fvar(); }


private:

	/*
	bool
	mycomp( const std::pair< core::Real, core::Size >& lhs,
	const std::pair< core::Real, core::Size >& rhs )
	{ return lhs.first < rhs.first; }
	*/

	void initialize( core::pose::Pose const &pose );

	void read_db( std::string const infile );
	void read_GUIP_matrix( std::string const infile );

	void fill_GUIP();

	utility::vector1< core::Size >
	get_n8( core::pose::Pose const &pose ) const;

	utility::vector1< core::Size >
	search_close_frags( utility::vector1< Vector > const &frag_crd,
		core::Size const n8,
		core::Size const npick
	) const;

	core::Real
	get_residue_similarity( core::conformation::Residue const& rsd,
		utility::vector1< Vector > const &frag_crd,
		core::Size const n8,
		utility::vector1< ResidueDepthDataCOP > const context_rdd,
		utility::vector1< core::Size > const &close_ids,
		core::Size const n ) const;

	core::Real
	get_simscore( core::chemical::AA const aa1,
		core::chemical::AA const aa2 ) const;

	core::Real
	compare_by_superposition( utility::vector1< ResidueDepthDataCOP > const & context_rdd,
		utility::vector1< Vector > const & frag_crd,
		utility::vector1< ResidueDepthDataCOP > const & context_rdd_ref,
		utility::vector1< Vector > const & frag_crd_ref
	) const;

	utility::vector1< ResidueDepthDataCOP >
	make_context( core::pose::Pose const &pose,
		core::Size const ires ) const;

	utility::vector1< ResidueDepthDataCOP >
	make_context_ref( core::Size const ipdb,
		core::Size const ires ) const;
	void
	fill_neighs( ResidueDepthData &rdd,
		utility::vector1< ResidueDepthData > const & pdb_rdds ) const;

private:
	std::string dbfile_, GUIP_matrix_file_;

	ResidueDepthCalculator RDC_;
	core::Size ncandidate1_, ncandidate2_;
	std::string similarity_mode_;

	core::Real gamma_;

	// pdbid/resno
	utility::vector1< utility::vector1< ResidueDepthDataCOP > > db_;
	utility::vector1< std::string > pdbid_;

	// index by 12-res
	std::map< core::Size, utility::vector1<ResidueDepthFrag> > frag_db_;

	utility::vector1< utility::vector1< core::Real > > GUIP_matrix_;

	// For filtering
	utility::vector1< bool > evalres_;
	core::Real maxdist_, mindist_;


}; //ResidueDepthFilter

} // namespace simple_filters
} // namespace protocols
#endif
