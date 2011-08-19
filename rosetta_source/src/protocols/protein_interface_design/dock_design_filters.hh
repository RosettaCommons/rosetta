// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/protein_interface_design/dock_design_filters.hh
/// @brief definition of filter classes for iterations of docking/design.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_dock_design_filters_hh
#define INCLUDED_protocols_protein_interface_design_dock_design_filters_hh


// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <utility/exit.hh>

//Auto Headers
#include <utility/options/keys/BooleanOptionKey.hh>


// C++ headers

// Unit headers
//#include <protocols/moves/DataMap.hh>

namespace protocols {
namespace protein_interface_design {

using protocols::filters::Filter;
using protocols::filters::FilterOP;
using protocols::filters::Filters_map;

class ResiduesInInterfaceFilter : public Filter
{
public:
	ResiduesInInterfaceFilter() : Filter( "ResInInterface" ) {}
	ResiduesInInterfaceFilter( core::Size const residues_in_interface_threshold, core::Size const rb_jump ) : Filter( "ResInInterface" ) {
		residues_in_interface_threshold_ = residues_in_interface_threshold;
		rb_jump_ = rb_jump;
	}
	bool apply( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return new ResiduesInInterfaceFilter( *this );
	}
	FilterOP fresh_instance() const{
		return new ResiduesInInterfaceFilter();
	}
  ResiduesInInterfaceFilter( ResiduesInInterfaceFilter const & init ) : 
	//utility::pointer::ReferenceCount(),
	Filter( init ) {
    residues_in_interface_threshold_ = init.residues_in_interface_threshold_;
		rb_jump_ = init.rb_jump_;
  };
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Size compute( core::pose::Pose const & pose ) const;
	virtual ~ResiduesInInterfaceFilter();
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Size residues_in_interface_threshold_, rb_jump_;
};

class AlaScan : public Filter
{
public :
	AlaScan() : Filter( "AlaScan" ) {}
	AlaScan( bool const chain1, bool const chain2, core::Size const repeats, core::Real const dist, core::scoring::ScoreFunctionCOP scorefxn, core::Size const jump, bool const symmetry );
	bool apply( core::pose::Pose const & ) const{ return true; }
	FilterOP clone() const {
		return new AlaScan( *this );
	}
	FilterOP fresh_instance() const{
		return new AlaScan();
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	void report_symmetry( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & ) const { return (0.0); };
	void chain1( bool const c1 ){ chain1_ = c1; }
	void chain2( bool const c2 ){ chain2_ = c2; }
	void dist( core::Real const d ){ distance_threshold_ = d; }
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	core::Real ddG_for_single_residue( core::pose::Pose const & pose, core::Size const resi ) const;
	bool chain1() const{ return chain1_; }
	bool chain2() const{ return chain2_; }
	core::Size repeats() const { return repeats_; }
	void repeats( core::Size const r ) { repeats_ = r; }
	core::Size jump() const { return jump_; }
	void jump( core::Size const j ) { jump_ = j; }
	core::Real dist() const{ return distance_threshold_; }
	virtual ~AlaScan();
	void repack( bool const repack );
	bool repack() const;
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	bool chain1_, chain2_;
	core::Size repeats_;
	core::Real distance_threshold_;
	core::Size jump_;
	bool symmetry_;
	core::scoring::ScoreFunctionOP scorefxn_;
	bool repack_; //dflt true; do you want to repack the partners in the bound and unbound states?
};

class ScoreTypeFilter : public Filter
{
public:
	ScoreTypeFilter() : Filter( "ScoreType" ) {}
	ScoreTypeFilter( core::scoring::ScoreFunctionCOP scorefxn, core::scoring::ScoreType const score_type, core::Real const score_type_threshold );
	bool apply( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return new ScoreTypeFilter( *this );
	}
	FilterOP fresh_instance() const{
		return new ScoreTypeFilter();
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const &pose ) const;
	virtual ~ScoreTypeFilter();
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Real score_type_threshold_;
	core::scoring::ScoreType score_type_;
	core::scoring::ScoreFunctionOP scorefxn_;
};

class InterfaceSasaFilter : public Filter
{
public:
	InterfaceSasaFilter();
	InterfaceSasaFilter( core::Real const lower_threshold, bool const hydrophobic=false, bool const polar=false );
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const; // which residue numbers are neighbors
	FilterOP clone() const;
	FilterOP fresh_instance() const;

	virtual ~InterfaceSasaFilter();
	void jump( core::Size const jump );
	core::Size jump() const;
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Real lower_threshold_;
	bool hydrophobic_, polar_; /// count only hydrophobics? polars?
	core::Size jump_; // dflt 1; across which jump to compute sasa
};

class NeighborTypeFilter : public Filter
{
public:
	NeighborTypeFilter() : Filter( "NeighborType" ){};
	NeighborTypeFilter( core::Size const target_residue, utility::vector1< bool > const residue_types, core::Real const distance_threshold ) :
		Filter( "NeighborType" ) {
		target_residue_ = target_residue; residue_types_ = residue_types; distance_threshold_ = distance_threshold;
	}
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	std::vector< core::Size > compute( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return new NeighborTypeFilter( *this );
	}
	FilterOP fresh_instance() const{
		return new NeighborTypeFilter();
	}
	void clear() { residue_types_.clear(); }
	virtual ~NeighborTypeFilter();
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Size target_residue_;
	utility::vector1< bool > residue_types_;
	core::Real distance_threshold_;
};

class ResidueBurialFilter : public Filter
{
public:
	ResidueBurialFilter() : Filter( "ResidueBurial"  ) {}
	ResidueBurialFilter( core::Size const target_residue, core::Size const neighbors, core::Real const distance_threshold ) :
		Filter( "ResidueBurial" ), target_residue_( target_residue ), neighbors_( neighbors ), distance_threshold_( distance_threshold ) {}
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Size compute( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return new ResidueBurialFilter( *this );
	}
	FilterOP fresh_instance() const{
		return new ResidueBurialFilter();
	}

	virtual ~ResidueBurialFilter();
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Size target_residue_;
	core::Size neighbors_;
	core::Real distance_threshold_;
};

class ResidueDistanceFilter : public Filter
{
public:
	ResidueDistanceFilter() : Filter( "ResidueDistance"  ) {}
	ResidueDistanceFilter( core::Size const res1, core::Size const res2, core::Real const distance_threshold ) :
		Filter( "ResidueDistance" ), res1_( res1 ), res2_( res2 ), distance_threshold_( distance_threshold ) {}
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return new ResidueDistanceFilter( *this );
	}
	FilterOP fresh_instance() const{
		return new ResidueDistanceFilter();
	}

	virtual ~ResidueDistanceFilter();
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Size res1_, res2_;
	core::Real distance_threshold_;
};

class DdgFilter : public Filter
{
public:
	DdgFilter();
	DdgFilter( core::Real const ddg_threshold, core::scoring::ScoreFunctionCOP scorefxn, core::Size const rb_jump=1, core::Size const repeats=1, bool const symmetry=false );
	bool apply( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return new DdgFilter( *this );
	}
	FilterOP fresh_instance() const{
		return new DdgFilter();
	}

	void repack( bool const repack );
	bool repack() const;
	void repeats( core::Size const repeats );
	core::Size repeats() const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~DdgFilter();
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Real ddg_threshold_; //dflt -15
	core::scoring::ScoreFunctionOP scorefxn_; //dflt NULL/score12 in cstrctr/rosettascripts
	core::Size rb_jump_; // dflt 1
	core::Size repeats_;//average of how many repeats? defaults to 1
	bool symmetry_; //dflt false
	bool repack_; //dflt true; Do you want to repack in the bound and unbound states (ddG) or merely compute the dG
};

/// @brief returns true if the number of hbonding partners to a particular residue exceeds a certain value
/// This filter is useful in conjunction with DesignMinimizeHbonds class
class HbondsToResidueFilter : public Filter
{
public:
	typedef core::Real Real;
	typedef core::Size Size;
public :
	HbondsToResidueFilter() : Filter( "HbondsToResidue" ) {}
	HbondsToResidueFilter( Size const resnum, Size const partners, Real const energy_cutoff=-0.5,
						   bool const backbone=false, bool const sidechain=true ) : Filter( "HbondsToResidue" ) {
		resnum_ = resnum; partners_ = partners; energy_cutoff_ = energy_cutoff; backbone_ = backbone;
		sidechain_ = sidechain;
		runtime_assert( backbone_ || sidechain_ );
		runtime_assert( partners_ );
		runtime_assert( energy_cutoff_ <= 0 );
	}
	bool apply( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return new HbondsToResidueFilter( *this );
	}
	FilterOP fresh_instance() const{
		return new HbondsToResidueFilter();
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Size compute( core::pose::Pose const & pose ) const;
	virtual ~HbondsToResidueFilter();
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	Size resnum_, partners_;
	Real energy_cutoff_;
	bool backbone_, sidechain_;
};

class EnergyPerResidueFilter : public Filter
{
public:
EnergyPerResidueFilter() : Filter( "EnergyPerResidue" ) {}

			EnergyPerResidueFilter( core::Size const resnum, core::scoring::ScoreFunctionCOP scorefxn,
								   core::scoring::ScoreType const score_type, core::Real const threshold,
								   bool const whole_interface = false, core::Size const rb_jump = 1,
								   core::Real const interface_distance_cutoff =  8.0 );

	EnergyPerResidueFilter( EnergyPerResidueFilter const &init );
	bool apply( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return new EnergyPerResidueFilter( *this );
	}
	FilterOP fresh_instance() const{
		return new EnergyPerResidueFilter();
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~EnergyPerResidueFilter();
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Size resnum_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreType score_type_;
	core::Real threshold_;
	bool whole_interface_;
	core::Size rb_jump_;
	core::Real interface_distance_cutoff_;

};


/// @brief filters based on an upper bound # of buried unsatisfied polar residues
class BuriedUnsatHbondFilter : public Filter
{
public:
	BuriedUnsatHbondFilter() : Filter( "BuriedUnsatHbonds" ) {}
	BuriedUnsatHbondFilter( core::Size const upper_threshold, core::Size const jump_num );
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return new BuriedUnsatHbondFilter( *this );
	}
	FilterOP fresh_instance() const{
		return new BuriedUnsatHbondFilter();
	}

	virtual ~BuriedUnsatHbondFilter();
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Size upper_threshold_;
	core::Size jump_num_;
};


class TerminusDistanceFilter : public Filter
{
public:
	TerminusDistanceFilter() : Filter( "TerminusDistance" ) {}
	//TerminusDistanceFilter( core::Size const distance, core::Size const jump_num );
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return new TerminusDistanceFilter( *this );
	}
	FilterOP fresh_instance() const{
		return new TerminusDistanceFilter();
	}

	virtual ~TerminusDistanceFilter();
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Size distance_;
	core::Size jump_num_;
};


} // protein_interface_design
} // devel


#endif /*INCLUDED_DOCK_DESIGN_FILTERS_H_*/

