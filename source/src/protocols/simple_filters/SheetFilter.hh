// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/filters/SheetFilters.hh
/// @brief header file for SheetSheetFilter class.
/// @details
/// @author Robert Vernon
/// @author James Thompson

#ifndef INCLUDED_protocols_simple_filters_SheetFilter_hh
#define INCLUDED_protocols_simple_filters_SheetFilter_hh

// Unit Headers
#include <protocols/simple_filters/AbinitioBaseFilter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/options/keys/OptionKeys.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>

// Utility headers

//// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {

class SheetFilter : public AbinitioBaseFilter {
public:

	/// c-tor and d-tor
	SheetFilter() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
	}
	~SheetFilter() override = default;

	filters::FilterOP clone() const override {
		return filters::FilterOP( new SheetFilter( *this ) ); }

	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new SheetFilter() ); }


	/// @brief Returns true if the given pose passes the filter, false otherwise.
	
	bool apply( core::pose::Pose const & pose ) const override;

	static int const max_nstr = 50;

	std::string name() const override {
		return "SheetFilter";
	}

private:

	// data
	/// @brief initialized_ is true if init() has been called, false otherwise.
	bool initialized_;
	mutable float handedness_score_;

	// methods

	void ingo_diamers(
		int const & nres,
		ObjexxFCL::FArray1A_int strnm,
		int natm,
		ObjexxFCL::FArray2D_float const & atmps,
		ObjexxFCL::FArray1D_int & rnm,
		ObjexxFCL::FArray1D_int & indC,
		ObjexxFCL::FArray1D_int & indN,
		ObjexxFCL::FArray2A_float strdm,
		ObjexxFCL::FArray1A_int inddm
	) const;

	void
	ingo_find_dir(
		int shnm,
		int hm,
		int & nstr,
		ObjexxFCL::FArray1A_int slct,
		ObjexxFCL::FArray1A_int order,
		ObjexxFCL::FArray1A_int strlbl,
		ObjexxFCL::FArray2A_float strdr,
		// int & proper, // Not used
		ObjexxFCL::FArray1A_int directions
	) const;

	void
	ingo_find_ord(
		int shnm,
		int hm,
		int & nstr,
		ObjexxFCL::FArray1A_int strlbl,
		ObjexxFCL::FArray2A_float dstrmin,
		float ngbhct,
		ObjexxFCL::FArray1A_int order,
		ObjexxFCL::FArray1A_int sequence,
		ObjexxFCL::FArray1A_int slct,
		int & rubbish
	) const;

	void
	ingo_hand(
		int const k,
		int const hm,
		ObjexxFCL::FArray1A_int stpppt,
		ObjexxFCL::FArray1A_int strtpt,
		int const str1,
		int const str2,
		ObjexxFCL::FArray1A_int order,
		ObjexxFCL::FArray2A_float strdr,
		int const nstr,
		int const nres,
		ObjexxFCL::FArray1A_int sequence,
		ObjexxFCL::FArray1A_int directions,
		ObjexxFCL::FArray1A_int scstr,
		ObjexxFCL::FArray2A_float lctn,
		bool const use_whole_helix,
		int & rubbish
	) const;

	void
	ingo_ident_sheets(
		int &nstr,
		ObjexxFCL::FArray2A_float dstrmin,
		float ngbhct,
		ObjexxFCL::FArray2A_int strprs,
		int & nsht,
		ObjexxFCL::FArray1A_int strsht,
		ObjexxFCL::FArray1A_int strlbl,
		int & rubbish
	) const;

	void
	ingo_lnl(
		int hm,
		ObjexxFCL::FArray1A_int order,
		int & nloc,
		int & nnloc
	) const;

	void
	ingo_locations(
		int const & nres,
		int natm,
		ObjexxFCL::FArray1D_int & indCA,
		ObjexxFCL::FArray2D_float const & atmps,
		ObjexxFCL::FArray2A_float lctn
	) const;

	void
	ingo_number_of_strands(
		int const & nres,
		ObjexxFCL::FArray1A_int scstr,
		int & nstr
	) const;

	void
	ingo_proper_sheets(
		int & nstr,
		int & nsht,
		float ngbhct,
		ObjexxFCL::FArray2A_float dstrmin,
		ObjexxFCL::FArray2A_float strdr,
		ObjexxFCL::FArray1A_int strsht,
		ObjexxFCL::FArray1A_int strlbl,
		int & rubbish,
		float & maxdist,
		float & mindotprodabs,
		ObjexxFCL::FArray1A_int strtpt,
		ObjexxFCL::FArray1A_int stpppt,
		ObjexxFCL::FArray2A_int locdsm,
		ObjexxFCL::FArray2A_float lctn,
		int const & nres
	) const;

	void
	ingo_sheet_stuff(
		int const & nres,
		ObjexxFCL::FArray1A_int scstr,
		int natm,
		ObjexxFCL::FArray2D_float const & atmps,
		ObjexxFCL::FArray1D_int & indN,
		ObjexxFCL::FArray1D_int & indCA,
		ObjexxFCL::FArray1D_int & indC,
		ObjexxFCL::FArray1D_int & rnm,
		float ngbhct,
		int & nstr,
		int & nsht,
		float & maxdt,
		float & mindp,
		int & nloc,
		int & nnloc,
		int & rubbish
	) const;

	void
	ingo_start_stop(
		int const & nres,
		ObjexxFCL::FArray1A_int scstr,
		int & nstr,
		ObjexxFCL::FArray1A_int strnm,
		ObjexxFCL::FArray1A_int strtpt,
		ObjexxFCL::FArray1A_int stpppt
	) const;

	void
	ingo_strand_dirs(
		int const & nres,
		int & nstr,
		ObjexxFCL::FArray1A_int strtpt,
		ObjexxFCL::FArray1A_int stpppt,
		ObjexxFCL::FArray2A_float lctn,
		ObjexxFCL::FArray2A_float strdr
	) const;

	void
	ingo_strand_dists_min(
		int const & nres,
		int & nstr,
		ObjexxFCL::FArray1A_int inddm,
		ObjexxFCL::FArray1A_int strnm,
		ObjexxFCL::FArray2A_float strdm,
		ObjexxFCL::FArray2A_float dstrmin,
		ObjexxFCL::FArray2A_int locdsm
	) const;

	void ingo_clean_ss(
		int const & nres,
		ObjexxFCL::FArray1A_int scstr
	) const;
};

} // filters
} // protocols

#endif
