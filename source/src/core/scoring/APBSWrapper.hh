// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ABPSWrapper.fwd.hh
/// @brief  APBSWrapper class declaration
/// @author Sachko Honda (honda@apl.washington.edu)

#ifndef INCLUDED_core_scoring_APBSWrapper_HH
#define INCLUDED_core_scoring_APBSWrapper_HH

#include <core/scoring/APBSWrapper.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// Utility headers
#include <numeric/xyzVector.hh>

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <vector>
#include <string>
#include <map>

namespace core {
namespace scoring {

///-------------------------------------------------------------------------------------
/// APBS wrapper
//--------------------------------------------------------------------------------------
class APBSWrapper : public utility::pointer::ReferenceCount {

	PQROP pqr;
	APBSConfigOP config;
	APBSResultOP result;
public:
	APBSWrapper(pose::Pose const & pose,
		std::map<std::string, bool> const & charged_residues,
		int dbg, bool calcenergy);
	virtual ~APBSWrapper();
	APBSResultCOP exec();
private:
	// Count the number of non-virtual atoms
	int count_atoms( pose::Pose const & pose ) const;
};
///-------------------------------------------------------------------------------------
/// PQR
//--------------------------------------------------------------------------------------
class PQR : public utility::pointer::ReferenceCount {
public:
	PQR(pose::Pose const &pose, int natoms, std::map<std::string, bool> const & charged_residues);
	virtual ~PQR();
	inline int get_natoms() { return natoms_; }
	static std::string const chains;
	int natoms_;
	std::vector<double> x, y, z, charge, radius;
};
///-------------------------------------------------------------------------------------
/// APBSResult
//--------------------------------------------------------------------------------------
class  APBSResult  : public utility::pointer::ReferenceCount {
public:
	APBSResult(int nsims, int natoms, int grid_dimes[3],
		int calcforce, int calcenergy,
		int write_pot, int write_charge, int write_smol,
		int write_kappa, int write_diel, int write_atompot ) ;
	virtual ~APBSResult();

	std::vector<double> esEnergy;
	std::vector<double> npEnergy;
	std::vector<double> dx;
	std::vector<double> dy;
	std::vector<double> dz;
	std::vector<double> qfx;
	std::vector<double> qfy;
	std::vector<double> qfz;
	std::vector<double> ibx;
	std::vector<double> iby;
	std::vector<double> ibz;
	std::vector<double> npx;
	std::vector<double> npy;
	std::vector<double> npz;
	std::vector<double> dbx;
	std::vector<double> dby;
	std::vector<double> dbz;
	int nwrites;
	double grid_meta[13];
	std::vector< std::vector<double> > grid_data;
};
///-------------------------------------------------------------------------------------
/// APBSConfig
//--------------------------------------------------------------------------------------
class APBSConfig : public utility::pointer::ReferenceCount{

public:

	///---------------------------------------------------------
	/// I_PARAM
	///---------------------------------------------------------
	class I_PARAM : public utility::pointer::ReferenceCount
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
		~I_PARAM();
		int * raw_array();
	};

	///---------------------------------------------------------
	/// R_PARAM
	///---------------------------------------------------------
	class R_PARAM : public utility::pointer::ReferenceCount
	{
		mutable double array[9];
	public:
		double pdie, nsdie, srad, swin, temp, sdens, gamma, smpbe_vol, smpbe_size;
		R_PARAM();
		~R_PARAM();
		double * raw_array();
	};

	APBSConfig(pose::Pose const & pose, int natoms, int dbg, bool calcenergy);
	virtual ~APBSConfig();

	// APBS debug level
	int dbg;

	// Number of simulations
	int nsims;

	// Number of atoms
	int natoms;

	// Grid generation parameters:
	double cfac;
	double fadd;
	double space;

	// APBS config parameters
	double grid[3];  // For manual & energy
	int dime[3];
	int pdime[3];  // Non 1 for parallel
	double glen[3]; // manual only
	double center[3]; // manual only
	double cglen[3];
	double fglen[3];
	double ccenter[3];
	double fcenter[3];

	double ionq[4], ionc[4], ionr[4];
	double ofrac; // parallel only
	I_PARAM i_param;
	R_PARAM r_param;

	bool calcenergy;  // calculate energy?

}; // end of APBSConfig

} // scoring
} // core
#endif
