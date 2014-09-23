// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/DipolarCouplingEnergy.cc
/// @brief  DC energy - Orientation dependent chemical shift
/// @author Lei Shi


//Unit headers
#include <core/scoring/methods/DipolarCouplingEnergy.hh>
#include <core/scoring/methods/DipolarCouplingEnergyCreator.hh>
#include <core/scoring/DipolarCoupling.hh>
#include <core/scoring/DipolarCoupling.fwd.hh>
#include <core/scoring/ScoreType.hh>
//Package headers

#include <core/conformation/Residue.hh>
//#include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
#include <core/pose/Pose.hh>
//#include <core/pose/datacache/CacheableDataType.hh>

//numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/xyzVector.hh>

#include <core/id/NamedAtomID.hh>

//utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>

//Objexx headers
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
// AUTO-REMOVED #include <ObjexxFCL/Fmath.hh>

#include <utility/io/ozstream.hh> //for dump_weights

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
//#include <basic/options/keys/dc.OptionKeys.gen.hh>

//C++ headers
#include <iostream>

//Auto Headers
#include <core/scoring/EnergyMap.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/format.hh>


static thread_local basic::Tracer tr( "core.scoring.DipolarCoupling" );

namespace core {
namespace scoring {
namespace methods {

using namespace ObjexxFCL::format;

/// @details This must return a fresh instance of the DipolarCouplingEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
DipolarCouplingEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new DipolarCouplingEnergy;
}

ScoreTypes
DipolarCouplingEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( dc );
	return sts;
}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
DipolarCouplingEnergy::DipolarCouplingEnergy() :
	parent( methods::EnergyMethodCreatorOP( new DipolarCouplingEnergyCreator ) )
{}

//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
EnergyMethodOP
DipolarCouplingEnergy::clone() const
{
  return new DipolarCouplingEnergy();
}

void DipolarCouplingEnergy::setup_for_scoring(
  pose::Pose & pose,
  ScoreFunction const &
) const
{
	dc_score_ = eval_dc( pose );
}

void DipolarCouplingEnergy::finalize_total_energy(
  pose::Pose &,
  ScoreFunction const &,
  EnergyMap & totals
) const
{
	totals[ dc ] = dc_score_;
}

void DipolarCouplingEnergy::setup_for_minimizing(
  pose::Pose & pose,
  ScoreFunction const &,
	kinematics::MinimizerMapBase const &
) const
{
	DipolarCoupling const& dc_data( * retrieve_DC_from_pose( pose ) );
	DipolarCoupling::DC_lines const& All_DC_lines( dc_data.get_DC_data() );
	DipolarCoupling::DC_lines::const_iterator it;
	Size ct = 0;
	for( it = All_DC_lines.begin(); it != All_DC_lines.end(); ++it) {
		id::AtomID atom1( pose.residue(it->res1()).atom_index(it->atom1()), it->res1());
		id::AtomID atom2( pose.residue(it->res2()).atom_index(it->atom2()), it->res2());
		//tr.Trace << "method: it->res1(): " << it->res1() << " it->atom1() " << it->atom1() << std::endl;
		//tr.Trace << "method: it->res2(): " << it->res2() << " it->atom2() " << it->atom2() << std::endl;
		//tr.Trace << "insert in atom-map " << atom1 << " " << atom2 << " " << atom3 << std::endl;
		++ct;
		utility::vector1< core::Size > atm1_map = atom2dc_map_.get( atom1 );
		utility::vector1< core::Size > atm2_map = atom2dc_map_.get( atom2 );
		atm1_map.push_back( ct );
		atm2_map.push_back( ct );
		atom2dc_map_.set( atom1, atm1_map );
		atom2dc_map_.set( atom2, atm2_map );
	}
}

//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
DipolarCoupling &
DipolarCouplingEnergy::dc_from_pose(
	pose::Pose & pose
) const
{
 	DipolarCouplingOP dc_info( retrieve_DC_from_pose( pose ) );
	if ( !dc_info ) {
		dc_info = new DipolarCoupling;
		store_DC_in_pose( dc_info, pose );
	}
	return *dc_info;
}

//////////////////////////////////////////////////////
//@brief main computation routine for DC energy... everything is happening here right now.
// this has to be spread out over different routines to make this energy yield derivatives
//////////////////////////////////////////////////////
Real DipolarCouplingEnergy::eval_dc(
  pose::Pose & pose
) const
{

	DipolarCoupling& dc_data( dc_from_pose( pose ) );
	Real score = dc_data.compute_dcscore( pose );
	return score;
}

void
DipolarCouplingEnergy::eval_atom_derivative(
			 id::AtomID const & aid,
			 pose::Pose const & pose,
    	 kinematics::DomainMap const &,
	   	 ScoreFunction const &,
		   EnergyMap const & score_weights,
		   Vector & F1,
		   Vector & F2
) const {

	if ( !atom2dc_map_.has( aid ) ) return; //damn this "has" isn't correct at all
	utility::vector1< Size > const dc_nrs( atom2dc_map_[ aid ] );
//	tr.Trace << " aid " << aid << std::endl;

	if ( dc_nrs.size() == 0 ) {
		//		tr.Trace << "no DC entry for " << aid << " skipping.. "<< std::endl;
		return;
	}

	Vector fij(0,0,0);

	for (core::Size ii=1; ii<=dc_nrs.size(); ++ii) {
		core::Size dc_nr = dc_nrs[ ii ];
		DipolarCoupling const& dc_cache( *retrieve_DC_from_pose( pose ) );
		utility::vector1< core::scoring::DC > All_DC_lines( dc_cache.get_DC_data() );
		runtime_assert( dc_nr <= All_DC_lines.size() );
		DC const& dc_data( All_DC_lines[ dc_nr ] );
		conformation::Residue const& rsd1( pose.residue( dc_data.res1() ) );
		conformation::Residue const& rsd2( pose.residue( dc_data.res2() ) );

		if ( aid.rsd() == dc_data.res1() && utility::trimmed_compare( rsd1.atom_name( aid.atomno() ), dc_data.atom1() ) ) {
//	    tr.Trace << "aid.rsd(): " << aid.rsd() << " rsd1.atom_name( aid.atomno() ) "<< rsd1.atom_name( aid.atomno() ) << " dc_data.atom1() " << dc_data.atom1() <<  std::endl;
			fij += dc_data.f1ij();
		} else if ( aid.rsd() == dc_data.res2() && utility::trimmed_compare( rsd2.atom_name( aid.atomno() ), dc_data.atom2() ) ){
//	    tr.Trace << "aid.rsd(): " << aid.rsd() << " rsd2.atom_name( aid.atomno() ) "<< rsd2.atom_name( aid.atomno() ) << " dc_data.atom2() " << dc_data.atom2() << std::endl;
			fij += dc_data.f2ij();
		} else return;

	}

//	tr.Trace << "fij[0]: " << fij[0]<< " fij[1]: " << fij[1]<< " fij[2]: " << fij[2]<< std::endl;
//	tr.Trace << "torsion gradient: " << aid << std::endl;
//	tr.Trace << "score_weights[ dc ]: " << score_weights[ dc ] << std::endl;
  //thanks to Will Sheffler:
	numeric::xyzVector<core::Real> atom_x = pose.xyz(aid);
	numeric::xyzVector<core::Real> const f2( fij );
	numeric::xyzVector<core::Real> const atom_y = atom_x - f2;   // a	"fake" atom in the direcion of the gradient
	numeric::xyzVector<core::Real> const f1( atom_x.cross( atom_y ) );

	F1 += score_weights[ dc ] * f1;
	F2 += score_weights[ dc ] * f2;

}

core::Size
DipolarCouplingEnergy::version() const
{
	return 1; // Initial versioning
}

} // methods
} // scoring
} // core
