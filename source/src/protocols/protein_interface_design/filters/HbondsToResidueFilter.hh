// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/HbondsToResidueFilter.hh
/// @brief definition of filter classes for iterations of docking/design.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)
/// @author Refactored considerably by Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory.

#ifndef INCLUDED_protocols_protein_interface_design_filters_HbondsToResidueFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_HbondsToResidueFilter_hh

#include <protocols/protein_interface_design/filters/HbondsToResidueFilter.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <utility/exit.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

#include <utility/vector1.hh>


// C++ headers

// Unit headers
//#include <basic/datacache/DataMap.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

using protocols::filters::Filter;
using protocols::filters::FilterOP;
using protocols::filters::Filters_map;

/// @brief returns true if the number of hbonding partners to a particular residue exceeds a certain value
/// This filter is useful in conjunction with DesignMinimizeHbonds class
class HbondsToResidueFilter : public Filter
{
public:
	typedef core::Real Real;
	typedef core::Size Size;
public :

	/// @brief Default constructor.
	///
	HbondsToResidueFilter();

	/// @brief Constructor
	///
	HbondsToResidueFilter(
		Size const resnum,
		Size const partners,
		Real const &energy_cutoff=-0.5,
		bool const backbone=false,
		bool const sidechain=true,
		bool const bb_bb=true,
		bool const from_other_chains=true,
		bool const from_same_chain=true
	);

	/// @brief Copy constructor
	///
	HbondsToResidueFilter( HbondsToResidueFilter const &src );


	bool apply( core::pose::Pose const & pose ) const;
	FilterOP clone() const {
		return FilterOP( new HbondsToResidueFilter( *this ) );
	}
	FilterOP fresh_instance() const{
		return FilterOP( new HbondsToResidueFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;

	/// @brief Actually compute the number of hydrogen bonds to the target residue.
	///
	core::Size compute( core::pose::Pose const & pose, core::Size const resnum_rosetta ) const;
	virtual ~HbondsToResidueFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

	/// @brief Set the minimum number of H-bond partners that this residue must have for the filter to pass.
	///
	inline void set_partners( core::Size const val) {
		partners_=val;
		return;
	}

	/// @brief Get the minimum number of H-bond partners that this residue must have for the filter to pass.
	///
	inline core::Size partners() const { return partners_; }

	/// @brief Set the threshhold for the hbond score term at which two residues are counted as being hydrogen bonded.
	///
	inline void set_energy_cutoff( core::Real const val) {
		runtime_assert_string_msg( val<=0, "Error in HbondsToResidueFilter::set_energy_cutoff(): The energy cutoff must be less than or equal to zero." );
		energy_cutoff_=val;
		return;
	}

	/// @brief Get the threshhold for the hbond score term at which two residues are counted as being hydrogen bonded.
	///
	inline core::Real energy_cutoff() const { return energy_cutoff_; }

	/// @brief Set whether to include backbone hydrogen bonds.
	/// @details I'm not sure that this is implemented properly.  (VKM -- 6 July 2015).
	inline void set_backbone( bool const val ) { backbone_ = val; return; }

	/// @brief Set whether to include backbone hydrogen bonds.
	/// @details I'm not sure that this is implemented properly.  (VKM -- 6 July 2015).
	inline void set_sidechain( bool const val ) { sidechain_ = val; return; }

	/// @brief Set whether to include backbone-backbone hydrogen bonds.
	/// @details I'm not sure that this is implemented properly.  (VKM -- 6 July 2015).
	inline void set_bb_bb( bool const val ) { bb_bb_ = val; return; }

	/// @brief Get whether to include backbone hydrogen bonds.
	/// @details I'm not sure that this is implemented properly.  (VKM -- 6 July 2015).
	inline bool backbone( ) const { return backbone_; }

	/// @brief Get whether to include backbone hydrogen bonds.
	/// @details I'm not sure that this is implemented properly.  (VKM -- 6 July 2015).
	inline bool sidechain( ) const { return sidechain_; }

	/// @brief Get whether to include backbone-backbone hydrogen bonds.
	/// @details I'm not sure that this is implemented properly.  (VKM -- 6 July 2015).
	inline bool bb_bb( ) const { return bb_bb_; }

	/// @brief Set the residue number (as a string to be parsed at apply time).
	///
	inline void set_resnum( std::string const &input ) { resnum_=input; return; }

	/// @brief Set the residue number (as an integer -- Rosetta numbering).
	inline void set_resnum( core::Size const val ) {
		runtime_assert(val>0);
		std::stringstream ss("");
		ss << val;
		set_resnum(ss.str());
		return;
	}

	/// @brief Get the residue number (a string to be parsed at apply time).
	///
	inline std::string resnum() const { return resnum_; }

	/// @brief Set whether hydrogen bonds from other chains should be counted.
	///
	inline void set_from_other_chains( bool const val ) { from_other_chains_=val; return; }

	/// @brief Get whether hydrogen bonds from other chains should be counted.
	///
	inline bool from_other_chains() const { return from_other_chains_; }

	/// @brief Set whether hydrogen bonds from the same chain should be counted.
	///
	inline void set_from_same_chain( bool const val ) { from_same_chain_=val; return; }

	/// @brief Get whether hydrogen bonds from the same chain should be counted.
	///
	inline bool from_same_chain() const { return from_same_chain_; }

	/// @brief Set the scorefunction to use for hbond calculation.
	///
	void set_scorefxn( core::scoring::ScoreFunctionCOP sfxn_in);

	/// @brief Set the ResidueSelector to use.
	/// @details Only hydrogen bonds between this residue and the residues selected by the ResidueSelector will be counted,
	/// if a ResidueSelector is provided.
	void set_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

private:

	/// @brief The current residue, stored as a string to be parsed at apply time.
	/// @details This could be a PDB number (e.g. 32A), a Rosetta number (e.g. 32), or a reference pose number (e.g. refpose(snapshot1,34).
	std::string resnum_;

	/// @brief The minimum number of H-bond partners that this residue must have for the filter to pass.
	///
	Size partners_;

	/// @brief The threshhold for the hbond score term at which two residues are counted as being hydrogen bonded.
	///
	Real energy_cutoff_;

	/// @brief Include backbone hydrogen bonds?
	/// @details I'm not sure that this is implemented properly.  (VKM -- 6 July 2015).
	bool backbone_;

	/// @brief Include sidechain hydrogen bonds?
	/// @details I'm not sure that this is implemented properly.  (VKM -- 6 July 2015).
	bool sidechain_;

	/// @brief Include backbone-backbone hydrogen bonds?
	/// @details I'm pretty sure that this is not implemented properly.  (VKM -- 6 July 2015).
	bool bb_bb_;

	/// @brief If true, hydrogen bonds from other chains will be counted.  True by default.
	///
	bool from_other_chains_;

	/// @brief If true, hydrogen bonds from the same chain will be counted.  True by default.
	///
	bool from_same_chain_;

	/// @brief Owning pointer to the scorefunction to use.
	///
	core::scoring::ScoreFunctionOP sfxn_;

	/// @brief Owning pointer to a ResidueSelector, optionally used to select the residues to count.
	/// @details Only hydrogen bonds between this residue and the residues selected by the ResidueSelector will be counted,
	/// if a ResidueSelector is provided.
	core::select::residue_selector::ResidueSelectorCOP selector_;

};

}
} // protein_interface_design
} // devel


#endif /*INCLUDED_DOCK_DESIGN_FILTERS_H_*/
