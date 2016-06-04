// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/enzdes/EnzFilters.hh
/// @brief definition of filter classes for enzdes in rosetta_scripts framework
/// @author Sagar Khare (khares@uw.edu)
/// @author Roland A Pache

#ifndef INCLUDED_protocols_enzdes_EnzFilters_hh
#define INCLUDED_protocols_enzdes_EnzFilters_hh

// unit headers
#include <protocols/enzdes/EnzFilters.fwd.hh>

// package headers
#include <protocols/enzdes/DesignVsNativeComparison.hh>
#include <protocols/match/downstream/LigandConformerBuilder.fwd.hh> //ResidueConformerFilter
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>

// Project Headers
#include <core/chemical/ResidueType.fwd.hh> //ResidueConformerFilter
#include <core/io/silent/SilentEnergy.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/filters/Filter.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace enzdes {


class LigDSasaFilter : public protocols::filters::Filter {

public:

	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef protocols::filters::Filters_map Filters_map;

public:
	LigDSasaFilter() : Filter( "LigDSasa" ) {}
	LigDSasaFilter( core::Real const lower_threshold, core::Real const upper_threshold );
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return FilterOP( new LigDSasaFilter( *this ) );
	}
	FilterOP fresh_instance() const{
		return FilterOP( new LigDSasaFilter() );
	}

	virtual ~LigDSasaFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Real lower_threshold_, upper_threshold_;
};

class LigBurialFilter : public protocols::filters::Filter
{
public:

	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef protocols::filters::Filters_map Filters_map;

public:
	LigBurialFilter() : Filter( "LigBurial"  ) {}
	LigBurialFilter( core::Size lig_id, core::Size const numneighbors, core::Real const distance_threshold ) :
		Filter( "LigBurial" ), lig_id_( lig_id), neighbors_( numneighbors ), distance_threshold_( distance_threshold ) {}
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Size compute( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return FilterOP( new LigBurialFilter( *this ) );
	}
	FilterOP fresh_instance() const{
		return FilterOP( new LigBurialFilter() );
	}

	virtual ~LigBurialFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Size lig_id_;
	core::Size neighbors_;
	core::Real distance_threshold_;
};


class LigInterfaceEnergyFilter : public protocols::filters::Filter
{
public:

	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef protocols::filters::Filters_map Filters_map;

public:
	LigInterfaceEnergyFilter() : Filter( "LigInterfaceEnergy" ) {}

	LigInterfaceEnergyFilter( core::scoring::ScoreFunctionOP scorefxn,
		core::Real const threshold, bool const include_cstE = false, core::Size const rb_jump = 1,
		core::Real const interface_distance_cutoff =  8.0 );

	LigInterfaceEnergyFilter( LigInterfaceEnergyFilter const &init );
	bool apply( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return FilterOP( new LigInterfaceEnergyFilter( *this ) );
	}
	FilterOP fresh_instance() const{
		return FilterOP( new LigInterfaceEnergyFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	core::Real constraint_energy( core::pose::Pose const & pose, int which_res ) const;

	virtual ~LigInterfaceEnergyFilter();

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &pose );


private:
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Real threshold_;
	bool include_cstE_;
	core::Size rb_jump_;
	core::Real interface_distance_cutoff_;

};

class EnzScoreFilter : public protocols::filters::Filter
{
public:

	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef protocols::filters::Filters_map Filters_map;

public:
	EnzScoreFilter() : Filter( "EnzScore" ) {}

	EnzScoreFilter(  core::Size const resnum, std::string const & cstid, core::scoring::ScoreFunctionOP scorefxn,
		core::scoring::ScoreType const score_type, core::Real const threshold, bool const whole_pose, bool const is_cstE
	);

	EnzScoreFilter( EnzScoreFilter const &init );
	bool apply( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return FilterOP( new EnzScoreFilter( *this ) );
	}
	FilterOP fresh_instance() const{
		return FilterOP( new EnzScoreFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	//core::Size get_resid_from_cstid( core::pose::Pose const & pose, core::Size const & cstid) const;
	virtual ~EnzScoreFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &pose );
private:
	core::Size resnum_;
	std::string cstid_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreType score_type_;
	core::Real threshold_;
	bool whole_pose_;
	bool is_cstE_;
};

class DiffAtomSasaFilter : public protocols::filters::Filter
{
public:

	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef protocols::filters::Filters_map Filters_map;

public:
	DiffAtomSasaFilter() : Filter( "DiffAtomBurial" ) {}

	DiffAtomSasaFilter( core::Size resid1, core::Size resid2, std::string atomname1, std::string atomane2, std::string sample_type );
	//DiffAtomSasaFilter( DiffAtomSasaFilter const &init );
	bool apply( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return FilterOP( new DiffAtomSasaFilter( *this ) );
	}
	FilterOP fresh_instance() const{
		return FilterOP( new DiffAtomSasaFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	bool compute( core::pose::Pose const & pose ) const;
	virtual ~DiffAtomSasaFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &pose );
private:
	core::Size resid1_, resid2_;
	std::string aname1_, aname2_, sample_type_;
};

class RepackWithoutLigandFilter : public protocols::filters::Filter
{
public:

	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef protocols::filters::Filters_map Filters_map;

public:
	RepackWithoutLigandFilter() : Filter( "RepackWithoutLigand" ) {}

	RepackWithoutLigandFilter( core::scoring::ScoreFunctionOP scorefxn, core::Real rms_thresh, core::Real energy_thresh, utility::vector1< core::Size > rms_target_res  );
	bool apply( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return FilterOP( new RepackWithoutLigandFilter( *this ) );
	}
	FilterOP fresh_instance() const{
		return FilterOP( new RepackWithoutLigandFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~RepackWithoutLigandFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &pose );
	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn){ scorefxn_ = scorefxn; }
	void set_cstid_list( std::string setting){ cstid_list_ = setting; }
private:
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Real rms_threshold_, energy_threshold_;
	bool calc_dE_, calc_rms_, use_cstids_, rms_all_rpked_;
	std::string cstid_list_;
	utility::vector1< core::Size > target_res_;
};

/// @brief tiny helper struct for EnzdesScoreFileFilter
struct ValueEvaluator
{
	enum CompareMode {
		SMALLER = 0,
		LARGER,
		EQUALS
	};

	ValueEvaluator( CompareMode const mode, core::Real const cutoff );
	~ValueEvaluator();

	bool
	value_passes( core::Real const value ) const;

	CompareMode mode_;
	core::Real cutoff_;
};


/// brief Class that calculates silent energies for an enzdes type scorefile
/// these silent energies can then be written into a scorefile.
/// class is derived from Filter because eventually it should be possible
/// to have this class read in a requirement file and return false if
/// any of the silent energies don't have the required value
class EnzdesScorefileFilter : public protocols::filters::Filter
{

	//constructor / destructor
public:

	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef protocols::filters::Filters_map Filters_map;

	EnzdesScorefileFilter();
	EnzdesScorefileFilter( EnzdesScorefileFilter const & other);
	~EnzdesScorefileFilter();

	//filter interface
public:
	bool apply( core::pose::Pose const & pose ) const;

	FilterOP clone() const {
		return FilterOP( new EnzdesScorefileFilter( *this ) );
	}

	FilterOP fresh_instance() const{
		return FilterOP( new EnzdesScorefileFilter() );
	}

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &pose );

	// specific stuff
public:

	void
	set_cstio(
		toolbox::match_enzdes_util::EnzConstraintIOCOP enzcst_io );

	void
	examine_pose(
		core::pose::Pose const & pose
	) const;

	utility::vector1< core::io::silent::SilentEnergy > const & silent_Es() const{
		return silent_Es_;
	}

	core::pose::PoseOP rnl_pose();

	/// @brief clear rnl pose to save some memory
	void
	clear_rnl_pose();


private:

	void
	initialize_value_evaluators_from_file( std::string const & filename );

	void
	compute_metrics_for_residue_subset(
		std::string sub_name,
		bool separate_out_constraints,
		core::pose::Pose const & calc_pose,
		utility::vector1< core::Size > const & res_subset ) const;

	void
	setup_pose_metric_calculators ( core::pose::Pose const & pose, bool separate_out_constraints ) const;

	bool no_packstat_calc_, native_comparison_;
	bool repack_no_lig_, keep_rnl_pose_;
	mutable core::pose::PoseOP rnl_pose_;
	core::scoring::ScoreFunctionOP sfxn_;
	toolbox::match_enzdes_util::EnzConstraintIOCOP enzcst_io_;

	mutable std::map< Size, utility::vector1< std::pair< std::string, std::string > > > residue_calculators_;
	mutable utility::vector1< std::pair< std::string, std::string > > native_compare_calculators_;

	DesignVsNativeComparisonOP native_comp_;
	mutable utility::vector1< core::io::silent::SilentEnergy > silent_Es_;
	utility::vector1< std::string > relevant_scoreterms_;
	mutable utility::vector1< std::pair< core::Size, core::Size > > spec_segments_;

	//filter necessary stuff
	std::string reqfile_name_;
	static std::map< std::string, std::map< std::string, ValueEvaluator > > evaluator_map_;


};

/// @brief filter that figures out which rotamer of a given
///rotamer lib is in the pose at apply time, and can be used
///to filter on it. supposed to be used for ligands,
///and for now only tested for them, but
///should also work with any other residue.
///can be used for example in specificity redesign, if one
///wants to divide up a bunch of designs according to the
///orientation in which they bind the ligand
class ResidueConformerFilter : public protocols::filters::Filter
{

public:

	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef protocols::filters::Filters_map Filters_map;

	ResidueConformerFilter();
	ResidueConformerFilter( ResidueConformerFilter const & other);
	~ResidueConformerFilter();

	bool apply( core::pose::Pose const & pose ) const;

	FilterOP clone() const {
		return FilterOP( new ResidueConformerFilter( *this ) );
	}

	FilterOP fresh_instance() const{
		return FilterOP( new ResidueConformerFilter() );
	}

	virtual void report( std::ostream &, core::pose::Pose const & pose ) const;

	virtual core::Real report_sm( core::pose::Pose const & pose ) const;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &pose );

	core::Size
	get_current_conformer( core::pose::Pose const & pose ) const;

private:

	void initialize_internal_data();

private:

	//seqpos and residue type the filter acts on
	core::chemical::ResidueType const * restype_;
	core::Size seqpos_;

	core::Size desired_conformer_; //the conformer that we'd like to have

	core::Real max_rms_;  //rms difference for two conformers to count as the same
	utility::vector1< core::Size > relevant_atom_indices_; //which atoms to consider


	//the actual work is being done by the ligconformer builder class,
	//since this already splits things up by conformers
	match::downstream::LigandConformerBuilderCOP lig_conformer_builder_;

};


} // enzdes
} // protocols


#endif /*INCLUDED_protocols_enzdes_EnzFilters_HH*/

