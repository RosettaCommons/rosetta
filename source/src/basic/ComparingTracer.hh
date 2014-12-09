// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/basic/ComparingTracer.hh
/// @brief  version of Tracer calss for to tracking down the source of an instability by compare
///         output with file.
/// @author Sergey Lyskov
///

#ifndef INCLUDED_basic_ComparingTracer_hh
#define INCLUDED_basic_ComparingTracer_hh

#include <basic/Tracer.hh>

#include <fstream>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/vector1_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <ostream>
#include <sstream>
#include <basic/Tracer.fwd.hh>



namespace basic {

class ComparingTracer : public basic::otstream
{
public:
	ComparingTracer(std::string const & file_name);
	virtual ~ComparingTracer();

protected:
	/// @brief overload member function.
	virtual void t_flush(std::string const & s);

private:
	ComparingTracer(ComparingTracer const & );

	//std::string  file_name_;
	std::fstream file_;

	/// @brief mode of operation:
	///                  if true: create model file, do not compare
	///                    false: read model file, and compare output against it, generating 'assert false'
	///                    if difference encountered.
	bool dry_run_;
};



} // namespace basic

#endif // INCLUDED_basic_comparing_tracer_hh

