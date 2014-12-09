// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TracerToFile.hh
/// @brief Class for a tracer that writes all output to a file.
/// @author Matt O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_basic_TracerToFile_hh
#define INCLUDED_basic_TracerToFile_hh



#include <basic/Tracer.hh>

#include <utility/pointer/owning_ptr.hh>

#include <fstream>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/vector1_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <ostream>
#include <sstream>
#include <basic/Tracer.fwd.hh>


namespace basic {


class TracerToFile : public basic::otstream
{
public:
  TracerToFile( std::string const & file_name );
  virtual ~TracerToFile();


protected:
  void t_flush( std::string const & s );

private:
  std::ofstream file_;
};



} // namespace basic


#endif // INCLUDED_basic_TracerToFile_hh
