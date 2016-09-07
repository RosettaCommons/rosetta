// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file CrossPeakList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_PeakFileFormat_xpk_hh
#define INCLUDED_protocols_noesy_assign_PeakFileFormat_xpk_hh


// Unit Headers
#include <protocols/noesy_assign/PeakFileFormat.hh>

// Package Headers

// Project Headers
#include <core/types.hh>
//#include <core/id/NamedAtomID.fwd.hh>
//#include <core/chemical/AA.hh>

// Utility headers
// #include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

//Auto Headers
// #include <numeric/numeric.functions.hh>
// #include <basic/prof.hh>
//#include <basic/Tracer.hh>
// #include <basic/options/option.hh>
// #include <basic/options/keys/abinitio.OptionKeys.gen.hh>
// #include <basic/options/keys/run.OptionKeys.gen.hh>
//#include <basic/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers

namespace protocols {
namespace noesy_assign {

class PeakFileFormat_xpk : public PeakFileFormat {
public:
	PeakFileFormat_xpk() {};
	//  PeakFileFormat_xpk( ResonanceListOP const& );
	//  virtual ~PeakFileFormat_xpk();

	void set_format_from_peak( CrossPeak const& ) override;
	void write_peak( std::ostream&, core::Size ct, CrossPeak const& ) const override;
	//  virtual void write_resonances( std::ostream&, CrossPeak const& ) const;
	//  virtual void write_strength( std::ostream&, CrossPeak const& ) const;
	//  virtual void write_assignments( std::ostream&, CrossPeak const&, std::string const& first_line_end ) const;
	void write_assignment( std::ostream&, PeakAssignment const& ) const override;
	void write_assignment_indent( std::ostream&, CrossPeak const& ) const override;
	void write_assignment_stats( std::ostream&, PeakAssignment& ) const override {}; //don't write these
	void write_nil_assignment( std::ostream& ) const override;
	//   virtual void read_resonances( std::istream&, CrossPeak& ) const;
	//   virtual void read_assignments( std::istream& is, std::istream& rest_line, CrossPeak& ) const;
	//   virtual void read_strength( std::istream&, CrossPeak& ) const;

	//   virtual CrossPeakOP read_peak( std::istream& ) const;
	//   virtual void read_header( std::istream& );
	//  virtual void write_header( std::ostream& );

	//   virtual void output_diagnosis( std::ostream& ) const;

	//   virtual void set_format_from_peak( CrossPeak const& );
	void write_header( std::ostream& ) override;
	//   virtual bool compatible_with_current_format( CrossPeak const& ) const;

	//static void register_options();
};

}
}

#endif
