// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_protocols_cryst_cryst_movers_creator_hh
#define INCLUDED_protocols_cryst_cryst_movers_creator_hh

// Project headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace cryst {

class ReportGradientsMoverCreator : public moves::MoverCreator {
public:
	moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	static  std::string mover_name();
};

class SetCrystWeightMoverCreator : public moves::MoverCreator {
public:
	moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	static  std::string mover_name();
};

class RecomputeDensityMapMoverCreator : public moves::MoverCreator {
public:
	moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	static  std::string mover_name();
};

class LoadDensityMapMoverCreator : public moves::MoverCreator {
public:
	moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	static  std::string mover_name();
};

class FitBfactorsMoverCreator : public moves::MoverCreator {
public:
	moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	static  std::string mover_name();
};

class UpdateSolventMoverCreator : public moves::MoverCreator {
public:
	moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	static  std::string mover_name();
};

class TagPoseWithRefinementStatsMoverCreator : public moves::MoverCreator {
public:
	moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	static  std::string mover_name();
};

class SetRefinementOptionsMoverCreator : public moves::MoverCreator {
public:
	moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	static  std::string mover_name();
};

}
}

#endif
