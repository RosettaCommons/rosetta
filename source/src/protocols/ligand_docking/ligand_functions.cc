// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/ligand_functions.cc
///
/// @brief
/// @author Ian W. Davis


#include <protocols/ligand_docking/ligand_functions.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <basic/Tracer.hh>

#include <protocols/ligand_docking/grid_functions.hh>

#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>

#include <core/pose/Pose.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

static thread_local basic::Tracer TR( "protocols.ligand_docking.ligand_functions" );

/// @brief Helper function.
core::scoring::constraints::ConstraintOP
torsion_constraints_from_mean_sd(
	core::Size rsd_no,
	core::Size chino,
	core::chemical::ResidueType const & rsd_type,
	utility::vector1< std::pair< core::Real, core::Real > > const & mean_sd_degrees
)
{
	using namespace core;
	using namespace core::scoring::constraints;
	using core::chemical::AtomIndices;
	using core::id::AtomID;

	TR << "Torsion constraint for Rsd " << rsd_no << " chi " << chino << " =";
	ConstraintCOPs csts;
	for(Size j = 1; j <= mean_sd_degrees.size(); ++j) {
		TR << " " << mean_sd_degrees[j].first;
		TR << " std dev " << mean_sd_degrees[j].second;
		// input is in degrees, but dihedral constraints deal in radians
		Real const chi_radians = numeric::conversions::radians( mean_sd_degrees[j].first );
		Real const stddev_radians = numeric::conversions::radians( mean_sd_degrees[j].second );
		core::scoring::func::FuncOP restr_func( new core::scoring::func::CircularHarmonicFunc( chi_radians, stddev_radians ) );
		AtomIndices chi_idx = rsd_type.chi_atoms(chino);
		ConstraintCOP constraint( new DihedralConstraint(
			AtomID(chi_idx[1], rsd_no),
			AtomID(chi_idx[2], rsd_no),
			AtomID(chi_idx[3], rsd_no),
			AtomID(chi_idx[4], rsd_no),
			restr_func
		) );
		csts.push_back( constraint );
	}
	TR << std::endl;
	ConstraintOP cst( new AmbiguousConstraint( csts ) );
	return cst;
}


core::scoring::constraints::ConstraintOP
torsion_constraints_from_rotamers(
	core::Size rsd_no,
	core::Size chino,
	utility::vector1< core::conformation::ResidueCOP > const & rsds,
	core::Real stddev_degrees
)
{
	using namespace core;
	chemical::ResidueType const & rsd_type = rsds[1]->type();

	Real const tol_d = stddev_degrees / 10.0; // within this range, considered to be the same minimum
	utility::vector1< std::pair< Real, Real> > minima_d;
	for(Size i = 1; i <= rsds.size(); ++i) {
		//runtime_assert( rsds[i].type().name() == rsdtype.name() );
		Real const chi_d = rsds[i]->chi(chino);
		bool found = false;
		for(Size j = 1; j <= minima_d.size(); ++j) {
			Real const min = minima_d[j].first;
			if(std::abs(min - numeric::nearest_angle_degrees(chi_d, min)) < tol_d) {
				found = true;
				break;
			}
		}
		if( !found ) {
			minima_d.push_back( std::make_pair(chi_d, stddev_degrees) );
		}
	}

	return torsion_constraints_from_mean_sd(rsd_no, chino, rsd_type, minima_d);
}


core::scoring::constraints::ConstraintOP
torsion_constraints_from_chi_rotamers(
	core::Size rsd_no,
	core::Size chino,
	core::chemical::ResidueType const & rsdtype
)
{
	return torsion_constraints_from_mean_sd( rsd_no, chino, rsdtype, rsdtype.chi_rotamers(chino) );
}


void
get_ligand_torsion_constraints(
	core::pose::Pose & pose,
	core::Size rsd_no,
	core::Real stddev_degrees,
	utility::vector1< core::scoring::constraints::ConstraintOP > & csts_out,
	bool const constrain_all_torsions_equally
)
{
	using namespace core;
	using namespace basic::options;
	using core::chemical::ResidueType;

	ResidueType const & rsdtype = pose.residue_type(rsd_no);
	for(Size i = 1; i <= rsdtype.nchi(); ++i) {
		bool has_diversity(false);
		if( rsdtype.chi_rotamers(i).size() == 0 ) {
				utility::vector1< core::conformation::ResidueOP > rotamers;
				rotamers_for_trials(pose, rsd_no, rotamers);
				if( rotamers.empty() || option[ OptionKeys::packing::use_input_sc ]() ) {
					rotamers.push_back( pose.residue(rsd_no).clone() );
					if ( (rotamers.size()==1) && constrain_all_torsions_equally) has_diversity=true;
				}
				for (Size j=2;j<=rotamers.size();++j) {
					Real const chi_d1 = rotamers[j]->chi(i);
					Real const chi_d2 = rotamers[j-1]->chi(i);
					if (std::abs(chi_d2 - chi_d1)> stddev_degrees ){
					//	TR<<"Debug: found diversity for CHI " << i << "values "<< chi_d2 << " and " << chi_d1 << std::endl;
						has_diversity=true;
						break;
					}
				}
			 	if (has_diversity || constrain_all_torsions_equally) csts_out.push_back( torsion_constraints_from_rotamers(rsd_no, i, rotamers, stddev_degrees) );
			 	else csts_out.push_back( torsion_constraints_from_rotamers(rsd_no, i, rotamers, 0.1/*make sure non-diverse torsions are highly constrained*/) );
		} else {
			csts_out.push_back( torsion_constraints_from_chi_rotamers(rsd_no, i, rsdtype) );
		}
	}
}


void
constrain_ligand_torsions(
	core::pose::Pose & pose,
	core::Real stddev_degrees,
	bool constrain_all_torsions_equally
)
{
	using namespace core;
	using namespace core::scoring::constraints;
	ConstraintSetOP new_constraint_set = pose.constraint_set()->clone();
	for(Size rsdno = 1; rsdno <= pose.total_residue(); ++rsdno) {
		if( pose.residue_type(rsdno).is_polymer() ) continue;
		utility::vector1< ConstraintOP > csts;
		get_ligand_torsion_constraints(pose, rsdno, stddev_degrees, csts, constrain_all_torsions_equally);
		for(Size cstno = 1; cstno <= csts.size(); ++cstno) {
			//csts[cstno]->show(TR);
			new_constraint_set->add_constraint( csts[cstno] );
		}
	}
	pose.constraint_set( new_constraint_set );
}

utility::vector1< core::Size >
get_ligand_seqpos(
	core::pose::Pose const & pose
)
{
	utility::vector1< core::Size > to_return;
	for( core::Size i =1; i <= pose.total_residue(); ++i){
		if( pose.residue_type(i).is_ligand() ) to_return.push_back( i );
	}
	return to_return;
}


} // namespace ligand_docking
} // namespace protocols
