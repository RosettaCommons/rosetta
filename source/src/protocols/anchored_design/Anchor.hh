// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/AnchoredDesign/Anchor.hh
/// @brief header for Anchor
/// @author Steven Lewis

#ifndef INCLUDED_protocols_anchored_design_Anchor_hh
#define INCLUDED_protocols_anchored_design_Anchor_hh

// Unit Headers
#include <protocols/anchored_design/Anchor.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace anchored_design {

/// @brief Anchor class provides for the "anchor" part of a scaffold in anchored interface design
/// @details Anchor class implements the code to read in an anchor file, stores the definition of the anchor, and provides access methods for setting move maps.  It reads in PDB-keyed information, but converts it to pose resid information and internally stores only the latter.
class Anchor : public utility::pointer::ReferenceCount
{

public:
	/// @brief default ctor uses option system to set anchorfile; can't set start and end until a pose is available
	Anchor();

	/// @brief input constructor for arbitrary anchor, if you know resids; ignores anchorfile datum
	Anchor( core::Size const start, core::Size const end );

	/// @brief input constructor from pose and anchorfile
	Anchor( core::pose::Pose const & pose, std::string const & filename );

	/// @brief input constructor from pose, assumes anchorfile
	Anchor( core::pose::Pose const & pose );

	virtual ~Anchor();

	/// @brief copy ctor
	Anchor( Anchor const & rhs );

	/// @brief assignment operator
	Anchor & operator=( Anchor const & rhs );

	/// @brief method to read an anchor file. Pose is necessary to reference against
	void read_anchorfile(core::pose::Pose const & pose, std::string const & filename);

	/// @brief method to read an anchor file. Pose is necessary to reference against.  Checks internal data for filename.
	void read_anchorfile(core::pose::Pose const & pose);

	/// @brief returns start of anchor (as pose resid)
	inline core::Size start() const { return start_; }
	/// @brief returns end of anchor (as pose resid)
	inline core::Size end() const { return end_; }

	//option system replacement
	/// @brief setter for filename for anchorfile
	void set_filename(std::string const & filename);
	/// @brief getter for filename for anchorfile
	std::string const & get_filename() const;
	/// @brief read option system into internal data
	void read_options();

private:

	//Anchor class converts pdb info to pose resid and stores here
	/// @brief start in pose resid terms
	core::Size start_;
	/// @brief end in pose resid terms
	core::Size end_;
	/// @brief filename source for anchorfile
	std::string anchorfile_;

}; //class Anchor

}//AnchoredDesign
}//protocols

#endif //INCLUDED_protocols_AnchoredDesign_Anchor_HH
