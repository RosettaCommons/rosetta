// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ABPSWrapper.fwd.hh
/// @brief  APBSWrapper class declaration
/// @author Sachko Honda (honda@apl.washington.edu)

#ifndef INCLUDED_core_scoring_APBSWrapper_HH
#define INCLUDED_core_scoring_APBSWrapper_HH

#include <core/scoring/APBSWrapper.fwd.hh>
#include <utility/VirtualBase.hh>

// Utility headers
#include <numeric/xyzVector.hh>

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <vector>
#include <string>
#include <map>

namespace core {
namespace scoring {

///-------------------------------------------------------------------------------------
/// APBS wrapper
//--------------------------------------------------------------------------------------
class APBSWrapper : public utility::VirtualBase {

	PQROP pqr;
	APBSConfigOP config;
	APBSResultOP result;
public:
	APBSWrapper(pose::Pose const & pose,
		id::AtomID_Map<bool> const & charged_atoms,
		id::AtomID_Map<bool> const & present_atoms,
		int dbg, bool calcenergy);
	~APBSWrapper() override;
	APBSResultCOP exec();
private:
	// Count the number of non-virtual atoms
	int count_atoms( id::AtomID_Map<bool> const & present_atoms ) const;
};
///-------------------------------------------------------------------------------------
/// PQR
//--------------------------------------------------------------------------------------
class PQR : public utility::VirtualBase {
public:
	PQR( pose::Pose const &pose,
		int natoms,
		id::AtomID_Map<bool> const & charged_atoms,
		id::AtomID_Map<bool> const & present_atoms);
	~PQR() override;
	inline int get_natoms() { return natoms_; }
	static std::string const chains;
	int natoms_;
	std::vector<core::Real> x, y, z, charge, radius;
};
///-------------------------------------------------------------------------------------
/// APBSResult
//--------------------------------------------------------------------------------------
class  APBSResult  : public utility::VirtualBase {
public:
	APBSResult(int nsims, int natoms, int const grid_dimes[3],
		int calcforce, int calcenergy,
		int write_pot, int write_charge, int write_smol,
		int write_kappa, int write_diel, int write_atompot ) ;
	~APBSResult() override;

	std::vector<core::Real> esEnergy;
	std::vector<core::Real> npEnergy;
	std::vector<core::Real> dx;
	std::vector<core::Real> dy;
	std::vector<core::Real> dz;
	std::vector<core::Real> qfx;
	std::vector<core::Real> qfy;
	std::vector<core::Real> qfz;
	std::vector<core::Real> ibx;
	std::vector<core::Real> iby;
	std::vector<core::Real> ibz;
	std::vector<core::Real> npx;
	std::vector<core::Real> npy;
	std::vector<core::Real> npz;
	std::vector<core::Real> dbx;
	std::vector<core::Real> dby;
	std::vector<core::Real> dbz;
	int nwrites;
	core::Real grid_meta[13];
	std::vector< std::vector<core::Real> > grid_data;
};
///-------------------------------------------------------------------------------------
/// APBSConfig
//--------------------------------------------------------------------------------------
class APBSConfig : public utility::VirtualBase{

public:

	///---------------------------------------------------------
	/// I_PARAM
	///---------------------------------------------------------
	class I_PARAM : public utility::VirtualBase
	{
		mutable int array[25];
	public:
		int sim_type; // 0=mg-manual, 1=mg-auto, 2=mg-parallel
		int nlev;
		int grid_centering_mode; // 0=use ccener[0-2], 1=mol 1
		int coarse_centering_mode; // 0=use ccener[0-2], 1=mol 1
		int fine_centering_mode;// 0=use ccener[0-2], 1=mol 1
		int chgm; // 0=sp10, 1=sp12, 2=sp14
		int pbe_mode; // 0=lpbe, 1=npbe, 2=lrpbe, 3=nrpbe, 4=smpbe
		int bcfl; //0=zero, 1=sdh, 2=mdh, 4=focus
		int srfm; //0=mol, 1=smol, 2=sp12, 3=sp14
		int calcforce; // 0=no, 1=yes
		int calcenergy; // 0=no, 1=yes
		int write_pot;
		int write_charge;
		int write_smol;
		int write_kappa;
		int write_diel;
		int write_atompot;
		int use_pot;
		int iparam19;
		int apol_calcforce;
		int apol_calcenergy;
		int nions;
		int use_charge;
		int use_kappa;
		int use_diel;

		I_PARAM();
		~I_PARAM() override;
		int * raw_array();
	};

	///---------------------------------------------------------
	/// R_PARAM
	///---------------------------------------------------------
	class R_PARAM : public utility::VirtualBase
	{
		mutable core::Real array[9];
	public:
		core::Real pdie, nsdie, srad, swin, temp, sdens, gamma, smpbe_vol, smpbe_size;
		R_PARAM();
		~R_PARAM() override;
		core::Real * raw_array();
	};

	APBSConfig(pose::Pose const & pose, int natoms, int dbg, bool calcenergy, id::AtomID_Map<bool> const & present_atoms);
	~APBSConfig() override;

	// APBS debug level
	int dbg;

	// Number of simulations
	int nsims;

	// Number of atoms
	int natoms;

	// Grid generation parameters:
	core::Real cfac;
	core::Real fadd;
	core::Real space;

	// APBS config parameters
	core::Real grid[3];  // For manual & energy
	int dime[3];
	int pdime[3];  // Non 1 for parallel
	core::Real glen[3]; // manual only
	core::Real center[3]; // manual only
	core::Real cglen[3];
	core::Real fglen[3];
	core::Real ccenter[3];
	core::Real fcenter[3];

	core::Real ionq[4], ionc[4], ionr[4];
	core::Real ofrac; // parallel only
	I_PARAM i_param;
	R_PARAM r_param;

	bool calcenergy;  // calculate energy?

}; // end of APBSConfig

} // scoring
} // core
#endif
