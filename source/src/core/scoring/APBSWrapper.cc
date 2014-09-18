// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ABPSWrapper.fwd.hh
/// @brief  APBSWrapper class definition
/// @author Sachko Honda (honda@apl.washington.edu)

#include <core/scoring/APBSWrapper.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// Numeric Headers

// Utility Headers
#include <basic/Tracer.hh>
#include <string>
#include <algorithm>
#include <numeric>
#include <vector>

#ifdef LINK_APBS_LIB

#include <apbs/routines.h>
#include <apbs_driver.h>

#endif
static thread_local basic::Tracer TR( "core/scoring/APBSWrapper" );

namespace core {
namespace scoring{
std::string const PQR::chains( " ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890" );
APBSWrapper::~APBSWrapper()
{
}

APBSWrapper::APBSWrapper(core::pose::Pose const & pose,
	      std::map<std::string, bool> const & charged_residues,
	      int dbg,
	      bool calcenergy)
{

		int natoms = count_atoms(pose);
		pqr = new PQR(pose, natoms, charged_residues);
    TR << "PQR data is prepared." << std::endl;
		config = new APBSConfig(pose, natoms, dbg, calcenergy);
		TR << "APBS config is prepared." << std::endl;
		result = new APBSResult(config->nsims, config->natoms, config->dime, 
														config->i_param.calcforce,
														config->i_param.calcenergy,
														config->i_param.write_pot,
														config->i_param.write_charge,
														config->i_param.write_smol,
														config->i_param.write_kappa,
														config->i_param.write_diel,
														config->i_param.write_atompot);
		TR << "APBS result data structure is prepared." << std::endl;
}
APBSResultCOP
APBSWrapper::exec() {

	int ret = 0;

#ifdef LINK_APBS_LIB
	const int nwrites = result->nwrites;
	std::vector< double * > raw_grid_data;
	for(int i=0; i<nwrites; i++ ) {
    raw_grid_data.push_back(result->grid_data[i].data());
	}

	ret = apbsdrv_(&pqr->natoms_,
									pqr->x.data(),
									pqr->y.data(),
									pqr->z.data(),
									pqr->radius.data(),
									pqr->charge.data(),
									config->r_param.raw_array(),
									config->i_param.raw_array(),
									config->grid,
									config->dime,
									config->pdime,
									config->glen,
									config->center,
									config->cglen,
									config->fglen,
									config->ccenter,
									config->fcenter,
									&config->ofrac,
									&config->dbg,
									config->ionq,
									config->ionc,
									config->ionr,
									result->esEnergy.data(),
									result->npEnergy.data(),
									result->dx.data(),
									result->dy.data(),
									result->dz.data(),
									result->qfx.data(),
									result->qfy.data(),
									result->qfz.data(),
									result->ibx.data(),
									result->iby.data(),
									result->ibz.data(),
									result->npx.data(),
									result->npy.data(),
									result->npz.data(),
									result->dbx.data(),
									result->dby.data(),
									result->dbz.data(),
									result->grid_meta,
									raw_grid_data.data() );
#endif

	if( ret != 0 ) {
		return NULL;
	}
	return result;
}

int APBSWrapper::count_atoms( core::pose::Pose const & pose ) const 
{
	int nres = pose.total_residue();
	int cntAtoms=0;
	for ( int i=1; i<= nres; ++i ) {
		conformation::ResidueCAP rsd = &pose.residue(i);
		for ( Size j=1; j<= rsd->natoms(); ++j ) {
			if ( rsd->atom_type(j).is_virtual() ) continue;
			++cntAtoms;
		}
	}
	return cntAtoms;
}

PQR::PQR(core::pose::Pose const & pose, int natoms,  
				 std::map<std::string, bool> const & charged_residues)
	: natoms_(natoms)
{
	int nres = pose.total_residue();
	int cntAtoms = 0;
 
	for ( int i=1; i<= nres; ++i ) {
		conformation::ResidueCAP rsd = &pose.residue(i);
		bool residue_charged = const_cast<std::map<std::string,bool>&>(charged_residues)[rsd->type().name()];
		for ( Size j=1; j<=rsd->natoms(); ++j ) {
			conformation::Atom const & atom( rsd->atom(j) );
			
			//skip outputing virtual atom unless specified.
			//fixed so that the last atom in atom type set can be 
			//something other than a virtual atom --steven combs
			if ( rsd->atom_type(j).is_virtual() ) continue;
			
			runtime_assert( rsd->chain() < chains.size() ); // silly restriction
			
			x.push_back(atom.xyz()(1));
			y.push_back(atom.xyz()(2));
			z.push_back(atom.xyz()(3));
			charge.push_back(residue_charged? pose.residue_type(i).atom(j).charge() : 0.);
			radius.push_back(rsd->atom_type(j).lj_radius());
			++cntAtoms;
		}
	}
}
PQR::~PQR(){}

APBSResult::APBSResult(int nsims, int natoms, int grid_dimes[3],
											 int calcforce, int calcenergy,
											 int write_pot, int write_charge, int write_smol,
											 int write_kappa, int write_diel, int write_atompot ) 
	: esEnergy(nsims),
		npEnergy(nsims)
{
	// allocate if force is to be calculated
  if(calcforce) {
		// allocate memory for these natoms-dimensioned data
		for( int i=0; i<natoms; i++ ) {
			dx.push_back(0);
			dy.push_back(0);
			dz.push_back(0);
		  qfx.push_back(0);
			qfy.push_back(0); 
			qfz.push_back(0);
			ibx.push_back(0);
			npx.push_back(0);
			dbx.push_back(0);
		}
	}
	if(calcenergy) {
	}
	nwrites = write_pot
		+ write_charge 
		+ write_smol
		+ write_kappa
		+ write_diel*3 
		+ write_atompot;
	grid_meta[0] = nwrites;
	for(int i=1; i<13; i++) grid_meta[i] = 0.;

	for(int i=0; i<nwrites; i++) {
		std::vector<double> i_data( grid_dimes[0] * grid_dimes[1] * grid_dimes[2], 0.);
		grid_data.push_back( i_data );
	}
}
APBSResult::~APBSResult(){}

APBSConfig::I_PARAM::I_PARAM() :
	sim_type(1),
	nlev(4),
	grid_centering_mode(0),
	coarse_centering_mode(0),
	fine_centering_mode(0),
	chgm(1),
	pbe_mode(1),
	bcfl(1),
	srfm(1),
	calcforce(0),
	calcenergy(0),
	write_pot(1),
  write_charge(0),
	write_smol(0),
	write_kappa(0),
	write_diel(0),
	write_atompot(0),
	use_pot(0),
	iparam19(0),
	apol_calcforce(0),
	apol_calcenergy(0),
	nions(4),
	use_charge(0),
	use_kappa(0),
	use_diel(0)
{ 
}
APBSConfig::I_PARAM::~I_PARAM() {}

int * APBSConfig::I_PARAM::raw_array()
{
	array[0] = sim_type;
	array[1] = nlev;
	array[2] = grid_centering_mode;
	array[3] = coarse_centering_mode;
	array[4] = fine_centering_mode;
	array[5] = chgm;
	array[6] = pbe_mode;
	array[7] = bcfl;
	array[8] = srfm;
	array[9] = calcforce;
	array[10] = calcenergy;
	array[11] = write_pot;
	array[12] = write_charge;
	array[13] = write_smol;
	array[14] = write_kappa;
	array[15] = write_diel;
	array[16] = write_atompot;
	array[17] = use_pot;
	array[18] = iparam19;
	array[19] = apol_calcforce;
	array[20] = apol_calcenergy;
	array[21] = nions;
	array[22] = use_charge;
	array[23] = use_kappa;
	array[24] = use_diel;
	return array;
}

APBSConfig::R_PARAM::R_PARAM() :
	pdie(4.),
	nsdie(80.),
	srad(1.4),
	swin(0.3),
	temp(310),
	sdens(10),
	gamma(0.105)
{
}
APBSConfig::R_PARAM::~R_PARAM() {}

double * APBSConfig::R_PARAM::raw_array()
{
	array[0] = pdie;
	array[1] = nsdie;
	array[2] = srad;
	array[3] = swin;
	array[4] = temp;
	array[5] = sdens;
	array[6] = gamma;
	array[7] = smpbe_vol;
	array[8] = smpbe_size;
	return array;
}

APBSConfig::APBSConfig(core::pose::Pose const & pose, int natomsIn, int dbgIn, bool calcenergyIn) 
	:
	dbg(dbgIn),
	nsims(1),
	natoms(natomsIn),
	cfac(1.7),
	fadd(20),
	space(0.5),
	ofrac(0.),
	calcenergy(calcenergyIn)
{
	ionq[0] = 1; ionq[1] = -1, ionq[2] = 2; ionq[3] = -2;
	ionc[0] = .15; ionc[1] = .15, ionc[2] = 0; ionc[3] = 0;
	ionr[0] = 2; ionr[1] = 2, ionr[2] = 2; ionr[3] = 2;
	dime[0] = dime[1] = dime[2] = 1;

	double min_r[] = {9999,9999,9999};
  double max_r[] = {-9999,-9999,-9999};
    // Find the min & max coords within the moleculer system to define the grid.
	for (Size ires=1; ires<=pose.total_residue(); ++ires) {
      for (Size iatom=1; iatom<=pose.residue(ires).natoms(); ++iatom) {
				for (int i=0; i<3; ++i) {
					min_r[i] = std::min(pose.residue(ires).xyz(iatom)[i], min_r[i]);
					max_r[i] = std::max(pose.residue(ires).xyz(iatom)[i], max_r[i]);
				}
      }
    }

    double length[3];
		for(int i=0; i<3; i++ ) {
			length[i] = max_r[i]- min_r[i];       // grid widths
			ccenter[i] = (min_r[i] + max_r[i])/2.;  // grid center coords
			fcenter[i] = ccenter[i];
			fglen[i] = length[i] + fadd;
			cglen[i] = length[i] * cfac;
			if (cglen[i] <= fglen[i]) cglen[i] = fglen[i] + 1;
			dime[i] = static_cast<int>(fglen[i] / space + 1);
		}
}
APBSConfig::~APBSConfig()
{}


}
}
