// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/indexed_structure_store/Datatypes.hh
/// @brief indexed_structure_store primitive datatypes.
/// @author Alex Ford <fordas@uw.edu>

#ifndef INCLUDED_protocols_indexed_structure_store_Datatypes_hh
#define INCLUDED_protocols_indexed_structure_store_Datatypes_hh

#include <stdint.h>
#include <math.h>
#include <array>
#include <limits>

#include <protocols/indexed_structure_store/Datatypes.fwd.hh>

namespace protocols
{
namespace indexed_structure_store
{

struct StructureEntry
{
	uint32_t id;
	char name[32];

	static bool comp_id(StructureEntry a, StructureEntry b)
	{
		return a.id < b.id;
	}
};

struct ResidueBackboneEntry
{
	float phi = 0;
	float psi = 0;
	float omega = 0;
};
static_assert(sizeof(ResidueBackboneEntry) % sizeof(float) == 0, "ResidueBackbone must be even multiple of float size.");

struct ResidueSidechainEntry
{
	float chi1 = std::numeric_limits<float>::quiet_NaN();
	float chi2 = std::numeric_limits<float>::quiet_NaN();
	float chi3 = std::numeric_limits<float>::quiet_NaN();
	float chi4 = std::numeric_limits<float>::quiet_NaN();
	char aa = 'x';
};
static_assert(sizeof(ResidueSidechainEntry) % sizeof(float) == 0, "ResidueSidechain must be even multiple of float size.");

struct ResidueOrientEntry
{
	std::array<float, 3> N = {{0, 0, 0}};
	std::array<float, 3> C = {{0, 0, 0}};
	std::array<float, 3> CA = {{0, 0, 0}};
	std::array<float, 3> O = {{0, 0,0}};
};
static_assert(sizeof(ResidueOrientEntry) % sizeof(float) == 0, "ResidueOrient must be even multiple of float size.");

struct ResidueEntry
{
	uint32_t structure_id = 0;
	uint32_t residue_id = 0;

	ResidueBackboneEntry bb;
	ResidueSidechainEntry sc;
	ResidueOrientEntry orient;

	bool chain_ending = false;

	static bool comp_id(ResidueEntry a, ResidueEntry b)
	{
		return (a.structure_id < b.structure_id) | ((a.structure_id == b.structure_id) & (a.residue_id < b.residue_id));
	}
};

static_assert(sizeof(ResidueEntry) % sizeof(float) == 0, "ResidueEntry must be even multiple of float size.");

}
}

#endif
