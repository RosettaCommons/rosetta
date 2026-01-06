// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/energy_methods/CCS_IMMSEnergy.cc
/// @author SM Bargeen Alam Turzo (turzo.1@osu.edu)

#include <core/energy_methods/CCS_IMMSEnergy.hh>
#include <core/energy_methods/CCS_IMMSEnergyCreator.hh>


#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/pose/Pose.hh>
#include <basic/prof.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>

#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/vector1.hh>

// Numeric header
#include <numeric/geometry/projection_area.hh>
#include <numeric/random/random_xyz.hh>
#include <numeric/statistics/functions.hh>
#include <numeric/NumericTraits.hh>

#include <basic/Tracer.hh>

namespace core {
namespace energy_methods {

static basic::Tracer tr( "core.energy_methods.CCS_IMMSEnergy" );
#define HEAVY_ATOM_EFFECTIVE_RADIUS 1.91 // effective vdw radius for heavy atoms
#define HYDROGEN_ATOM_EFFECTIVE_RADIUS 1.20 // effective vdw radius for hydrogen

/// @details This must return a fresh instance of the CCS_IMMSEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
CCS_IMMSEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const & ) const
{
	return utility::pointer::make_shared< CCS_IMMSEnergy >();
}

core::scoring::ScoreTypes
CCS_IMMSEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back(ccs_imms);
	return sts;
}

//IM complex score
core::scoring::methods::EnergyMethodOP
CCS_IMMSComplexEnergyCreator::create_energy_method(core::scoring::methods::EnergyMethodOptions const &) const
{
    return utility::pointer::make_shared<CCS_IMMSComplexEnergy>();
}

core::scoring::ScoreTypes
CCS_IMMSComplexEnergyCreator::score_types_for_method() const {
    using namespace core::scoring;
    ScoreTypes sts;
    sts.push_back(ccs_imms_complex);
    return sts;
}

//IMMS with Cryo EM Energy
core::scoring::methods::EnergyMethodOP
CCS_IMMS_with_CryoEMEnergyCreator::create_energy_method(core::scoring::methods::EnergyMethodOptions const &) const
{
    return utility::pointer::make_shared<CCS_IMMS_with_CryoEMEnergy>();
}

core::scoring::ScoreTypes
CCS_IMMS_with_CryoEMEnergyCreator::score_types_for_method() const {
    using namespace core::scoring;
    ScoreTypes sts;
    sts.push_back(ccs_imms_cryoem);
    return sts;
}
// constructor
CCS_IMMSEnergy::CCS_IMMSEnergy() :
	core::scoring::methods::WholeStructureEnergy( utility::pointer::make_shared< CCS_IMMSEnergyCreator >() )
{
	runtime_assert_string_msg( basic::options::option[ basic::options::OptionKeys::score::ccs_exp ].user(), "Error in CCS_IMMSEnergy constructor: the use must specify an experimental value using the -score:ccs_exp flag." );
	ccs_exp_  = basic::options::option[ basic::options::OptionKeys::score::ccs_exp ]();
	nrot_     = basic::options::option[ basic::options::OptionKeys::score::ccs_nrots ]();
	prad_     = basic::options::option[ basic::options::OptionKeys::score::ccs_prad ]();
}

CCS_IMMSComplexEnergy::CCS_IMMSComplexEnergy() :
    core::scoring::methods::WholeStructureEnergy(utility::pointer::make_shared<CCS_IMMSComplexEnergyCreator>())
{
    runtime_assert_string_msg(basic::options::option[basic::options::OptionKeys::score::ccs_exp].user(), 
        "Error in CCS_IMMSComplexEnergy constructor: the user must specify an experimental value using the -score:ccs_exp flag.");
    ccs_exp_ = basic::options::option[basic::options::OptionKeys::score::ccs_exp]();
    nrot_ = basic::options::option[basic::options::OptionKeys::score::ccs_nrots]();
    prad_ = basic::options::option[basic::options::OptionKeys::score::ccs_prad]();
    // Always use multimer for this energy
    use_multimer_ = true;
    
    // Using hardcoded values for lb_, ub_, and weight_
    // These values match the defaults defined in the header file
    // lb_ = 10;  // Default value 
    // ub_ = 100; // Default value
    // weight_ = 100; // Default value
}

CCS_IMMS_with_CryoEMEnergy::CCS_IMMS_with_CryoEMEnergy() :
    core::scoring::methods::WholeStructureEnergy(utility::pointer::make_shared<CCS_IMMS_with_CryoEMEnergyCreator>())
{
    runtime_assert_string_msg(basic::options::option[basic::options::OptionKeys::score::ccs_exp].user(), 
        "Error in CCS_IMMS_with_CryoEMEnergy constructor: the user must specify an experimental value using the -score:ccs_exp flag.");
    ccs_exp_ = basic::options::option[basic::options::OptionKeys::score::ccs_exp]();
	nrot_ = basic::options::option[basic::options::OptionKeys::score::ccs_nrots]();
    prad_ = basic::options::option[basic::options::OptionKeys::score::ccs_prad]();
    
    // Never use multimer for this energy
    use_multimer_ = false;
    
    // Using hardcoded values for lb_, ub_, and weight_
    // These values match the defaults defined in the header file
    // lb_ = 10;  // Default value 
    // ub_ = 100; // Default value
    // weight_ = 100; // Default value
}


/// clone
core::scoring::methods::EnergyMethodOP
CCS_IMMSEnergy::clone() const {
	return utility::pointer::make_shared< CCS_IMMSEnergy >( *this );
}

core::scoring::methods::EnergyMethodOP
CCS_IMMSComplexEnergy::clone() const {
    return utility::pointer::make_shared<CCS_IMMSComplexEnergy>(*this);
}

core::scoring::methods::EnergyMethodOP
CCS_IMMS_with_CryoEMEnergy::clone() const {
    return utility::pointer::make_shared<CCS_IMMS_with_CryoEMEnergy>(*this);
}
/* Up to close everything is working */
/////////////////////////////////////////////////////////////////////////////
// collision cross section calculation
/////////////////////////////////////////////////////////////////////////////

core::Real
parcs_ccs(core::pose::Pose &mypose , core::Size const nrot, core::Real const prad, bool use_multimer) // Function for CCS calculation
{
	core::Size tota_natoms = mypose.total_atoms();
	utility::vector1 < core::Real > eff_radius;
	eff_radius.reserve(tota_natoms);
	core::Real accumulator = 0.0;
	core::Size counter = 0;
	for ( core::Size rt = 1; rt <= nrot; ++rt ) { // Do everything below for every rotations
		numeric::xyzMatrix< core::Real > M = numeric::random::random_rotation(); // Generate a Random Matrix.
		numeric::xyzVector<core::Real> zero_vec; // Since we are not translating, get zero vector
		zero_vec.zero(); // Fill it with zeroes
		mypose.apply_transform_Rx_plus_v(M, zero_vec); // Apply the random rotations.

		// Vectors to store x, y, z and eff_radius of protein
		utility::vector1 < core::Real > xcoord;
		utility::vector1 < core::Real > ycoord;
		utility::vector1 < core::Real > zcoord;
		xcoord.reserve(tota_natoms);
		ycoord.reserve(tota_natoms);
		zcoord.reserve(tota_natoms);
		for ( core::Size res = 1; res <= mypose.size(); ++res ) {
			for ( core::Size atom = 1; atom <= mypose.residue(res).natoms(); ++atom ) {
				if ( !mypose.residue(res).atom_type(atom).is_virtual() ) { // Avoid virtual atoms
					xcoord.push_back(mypose.residue(res).xyz(atom).x());
					ycoord.push_back(mypose.residue(res).xyz(atom).y());
					zcoord.push_back(mypose.residue(res).xyz(atom).z());
					if ( rt==1 ) {
						if ( !mypose.residue(res).atom_type(atom).is_hydrogen() ) { // If atom is not hydrogen store 1.91 as eff radius
							eff_radius.push_back(HEAVY_ATOM_EFFECTIVE_RADIUS);
						} else eff_radius.push_back(HYDROGEN_ATOM_EFFECTIVE_RADIUS);  // else store 1.20 as atom radius
					}
				}
			}
		}
		accumulator += numeric::geometry::projection_area(xcoord,ycoord,eff_radius,prad);
		accumulator += numeric::geometry::projection_area(xcoord,zcoord,eff_radius,prad);
		accumulator += numeric::geometry::projection_area(ycoord,zcoord,eff_radius,prad);
		counter += 3;
	}
	core::Real ccs = ( accumulator /  static_cast<core::Real>(counter) ); //Compute the average as the sum over the number of elements, and return it.
	if (use_multimer){
		ccs = 0.8 * std::pow(ccs, 1.04);
	}
	return ccs;
}




/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

core::Real
CCS_IMMSEnergy::calc_IMMS_score(const core::Real &CCS_pred, const core::Real &CCS_exp)
const {
	core::Real diff = std::abs(CCS_pred - CCS_exp);
	core::Real ub = 100;
	core::Real lb = 10;
	core::Real CCS_score;
	if ( diff <=lb ) {
		CCS_score = 0;
	} else if ( diff < ub ) {
		core::Real b = -(diff - ub) / (ub-lb);
		core::Real b2= b*b;
		core::Real b3= b2*b;
		CCS_score = 100*( 2*b3 -3*b2 +1 );
	} else {
		CCS_score = 100.0;
	}
	return CCS_score;
}

core::Real
CCS_IMMSComplexEnergy::calc_IMMS_score(const core::Real &CCS_pred, const core::Real &CCS_exp) const {
    core::Real diff = std::abs(CCS_pred - CCS_exp);
	core::Real ub = 1350;
	core::Real lb = 60;
    core::Real CCS_score;
    
    if (diff <= lb_) {
        CCS_score = 0;
    } else if (diff < ub_) {
        core::Real b = -(diff - ub_) / (ub_ - lb_);
        core::Real b2 = b*b;
        core::Real b3 = b2*b;
        CCS_score = 100 * (2*b3 - 3*b2 + 1);
    } else {
        CCS_score = 100;
    }
    
    return CCS_score;
}
core::Real
CCS_IMMS_with_CryoEMEnergy::calc_IMMS_score(const core::Real &CCS_pred, const core::Real &CCS_exp) const {
    core::Real diff = std::abs(CCS_pred - CCS_exp);
	core::Real ub = 125;
	core::Real lb = 25;
    core::Real CCS_score;
    
    if (diff <= lb_) {
        CCS_score = 0;
    } else if (diff < ub_) {
        core::Real b = -(diff - ub_) / (ub_ - lb_);
        core::Real b2 = b*b;
        core::Real b3 = b2*b;
        CCS_score = 100 * (2*b3 - 3*b2 + 1);
    } else {
        CCS_score = 125;
    }
    
    return CCS_score;
}

void
CCS_IMMSEnergy::finalize_total_energy(core::pose::Pose &mypose, core::scoring::ScoreFunction const &, core::scoring::EnergyMap & emap) const
{
	core::Real pred_ccs = parcs_ccs(mypose, nrot_, prad_);
	core::Real ccs_score = calc_IMMS_score(pred_ccs, ccs_exp_);
	emap[ core::scoring::ccs_imms] = ccs_score;
}

void
CCS_IMMSComplexEnergy::finalize_total_energy(core::pose::Pose &mypose, core::scoring::ScoreFunction const &, core::scoring::EnergyMap & emap) const
{
	core::Real pred_ccs = parcs_ccs(mypose, nrot_, prad_, use_multimer_);
	core::Real ccs_score = calc_IMMS_score(pred_ccs, ccs_exp_);
	emap[ core::scoring::ccs_imms_complex] = ccs_score;
}

void
CCS_IMMS_with_CryoEMEnergy::finalize_total_energy(core::pose::Pose &mypose, core::scoring::ScoreFunction const &, core::scoring::EnergyMap & emap) const
{
	core::Real pred_ccs = parcs_ccs(mypose, nrot_, prad_, use_multimer_);
	core::Real ccs_score = calc_IMMS_score(pred_ccs, ccs_exp_);
	emap[ core::scoring::ccs_imms_cryoem] = ccs_score;
}

void 
CCS_IMMSEnergy::indicate_required_context_graphs(utility::vector1<bool>&) const 
{
	// No context graphs required
}

void 
CCS_IMMSComplexEnergy::indicate_required_context_graphs(utility::vector1<bool>&) const 
{
	// No context graphs required
}

void 
CCS_IMMS_with_CryoEMEnergy::indicate_required_context_graphs(utility::vector1<bool>&) const 
{
	// No context graphs required
}

core::Size 
CCS_IMMSEnergy::version() const 
{ 
	return 1; 
}

core::Size 
CCS_IMMSComplexEnergy::version() const 
{ 
	return 1; 
}

core::Size 
CCS_IMMS_with_CryoEMEnergy::version() const 
{ 
	return 1; 
}
}
}
