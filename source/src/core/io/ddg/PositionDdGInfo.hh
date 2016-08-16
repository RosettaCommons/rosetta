// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_core_io_ddg_PositionDdGInfo_hh
#define INCLUDED_core_io_ddg_PositionDdGInfo_hh

// Unit header
#include <core/io/ddg/PositionDdGInfo.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <fstream>

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <map>


namespace core {
namespace io {
namespace PositionDdGInfo {

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
using namespace scoring;

// the following functions added by flo may '11 to read in ddg prediction output files
// so that this info can be used in other protocols

/// @brief small helper class  that stores the ddGs for mutations
/// at a given position. camel case gets weird when trying to write
/// words containing ddG...
class PositionDdGInfo : public utility::pointer::ReferenceCount {

public:

	PositionDdGInfo(
		core::Size seqpos,
		core::chemical::AA wt_aa
	);

	virtual ~PositionDdGInfo();

	void
	add_mutation_ddG(
		core::chemical::AA aa,
		core::Real ddG
	);

	core::Size
	seqpos() const {
		return seqpos_; }

	core::chemical::AA
	wt_aa() const {
		return wt_aa_; }

	std::map< core::chemical::AA, core::Real > const &
	mutation_ddGs() const {
		return mutation_ddGs_; }


private:
	core::Size seqpos_;
	core::chemical::AA wt_aa_;
	std::map< core::chemical::AA, core::Real > mutation_ddGs_;

};


/// @brief function that reads in a ddg predictions out file
/// and returns the info in it as a map of PositionDdGInfo
const std::map< core::Size, PositionDdGInfoOP >
read_ddg_predictions_file( std::string filename );


}//namespace ddG
}//namespace io
} //namespace core

#endif
