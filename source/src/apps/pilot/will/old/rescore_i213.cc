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


#include <basic/database/open.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
//#include <core/conformation/symmetry/SymDof.hh>
//#include <core/conformation/symmetry/SymmetricConformation.hh>
//#include <core/conformation/symmetry/SymmetryInfo.hh>
//#include <core/conformation/symmetry/util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/packing/HolesParams.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <fstream>
#include <iostream>
#include <math.h>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <protocols/docking/util.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/Mover.hh>
//#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
//#include <protocols/minimization_packing/symmetry/SymMinMover.hh>
#include <protocols/viewer/viewers.hh>
#include <sstream>
#include <string>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <protocols/protein_interface_design/movers/ddG.hh>
#include <protocols/pose_metric_calculators/RotamerBoltzCalculator.hh>
#include <protocols/calc_taskop_filters/RotamerBoltzmannWeight.hh>
#include <core/pack/task/ResfileReader.hh>
//Auto Headers

static basic::Tracer TR( "rescore_i213" );

using core::Size;
using core::Real;
using core::pose::Pose;
using std::string;
using utility::vector1;
typedef vector1<Size> Sizes;

// Pose needs to be scored before this will work.
void
new_sc(core::pose::Pose &pose, Real& int_area, Real& sc) {
	using namespace core;
	Size nres_monomer = pose.size()/2;
	core::scoring::sc::ShapeComplementarityCalculator scc;
	scc.Init();
	for ( Size i=1; i<=nres_monomer; ++i ) scc.AddResidue(0, pose.residue(i));
	for ( Size i=1; i<=nres_monomer; ++i ) scc.AddResidue(1, pose.residue(i+nres_monomer));
	if ( scc.Calc() ) {
		sc = scc.GetResults().sc;
		int_area = scc.GetResults().surface[2].trimmedArea;
	}
}

utility::vector1<Real>
sidechain_sasa(Pose const & pose, Real probe_radius) {
	using core::id::AtomID;
	utility::vector1<Real> rsd_sasa(pose.size(),0.0);
	core::id::AtomID_Map<Real> atom_sasa;
	core::id::AtomID_Map<bool> atom_mask;
	core::pose::initialize_atomid_map(atom_sasa,pose,0.0);
	core::pose::initialize_atomid_map(atom_mask,pose,false);
	for ( Size i = 1; i <= pose.size(); i++ ) {
		for ( Size j = 1; j <= pose.residue(i).nheavyatoms(); j++ ) {
			atom_mask[AtomID(j,i)] = true;
		}
	}
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, false, atom_mask );
	utility::vector1<Real> sc_sasa(pose.size(),0.0);
	for ( Size i = 1; i <= pose.size(); i++ ) {
		// Use CA as the side chain for Glys
		if ( pose.residue(i).name3()=="GLY" ) sc_sasa[i] += atom_sasa[AtomID(2,i)];
		for ( Size j = 5; j <= pose.residue(i).nheavyatoms(); j++ ) {
			sc_sasa[i] += atom_sasa[AtomID(j,i)];
		}
	}
	return sc_sasa;
}


Real
average_degree (Pose const &pose, vector1<Size> mutalyze_pos, Real distance_threshold=10.0)
{

	Size nres_monomer = pose.size()/2;
	Size count_neighbors( 0 );

	for ( Size i = 1; i <= mutalyze_pos.size(); ++i ) {
		Size ires = mutalyze_pos[i];
		core::conformation::Residue const resi( pose.conformation().residue( ires ) );
		Size resi_neighbors( 0 );
		for ( Size jres = 1; jres <= nres_monomer; ++jres ) {
			core::conformation::Residue const resj( pose.residue( jres ) );
			Real const distance( resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) );
			if ( distance <= distance_threshold ) {
				++count_neighbors;
				++resi_neighbors;
			}
		}
		//TR << "avg_deg of " << resi.name3() << ires << " = " << resi_neighbors << std::endl;
	}

	return( (Real) count_neighbors / mutalyze_pos.size() );

}


vector1<Size> get_des_pos(core::pose::Pose & pose_for_design) {
	// Find out which positions are near the inter-subunit interfaces
	// These will be further screened below, then passed to design() and minimize()
	Real const contact_dist = 10.0;
	Real const contact_dist_sq = contact_dist * contact_dist;

	core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();

	////////////////////////////////////
	vector1<bool> indy_resis(pose_for_design.size(),false);
	for ( Size i = 1; i <= pose_for_design.size()/2; ++i ) indy_resis[i] = true;
	vector1<bool> subunit_index(pose_for_design.size());
	for ( Size i = 1; i <= pose_for_design.size(); ++i ) subunit_index[i] = (6*(i-1))/pose_for_design.size();
	vector1<Size> intra_subs; intra_subs.push_back(1); intra_subs.push_back(2); intra_subs.push_back(3);


	Sizes interface_pos;
	for ( Size ir=1; ir<= pose_for_design.size()/2; ir++ ) {
		std::string atom_i = "";
		if ( pose_for_design.residue(ir).name3() == "GLY" ) {
			atom_i = "CA";
		} else {
			atom_i = "CB";
		}
		for ( Size jr=1+pose_for_design.size()/2; jr<= pose_for_design.size(); jr++ ) {
			std::string atom_j = "";
			if ( pose_for_design.residue(jr).name3() == "GLY" ) {
				atom_j = "CA";
			} else {
				atom_j = "CB";
			}
			if ( pose_for_design.residue(ir).xyz(atom_i).distance_squared(pose_for_design.residue(jr).xyz(atom_j)) <= contact_dist_sq ) {
				interface_pos.push_back(ir);
				break;
			}
		}
	}

	return interface_pos;

}


core::io::silent::SilentFileData sfd;

void rescore(core::pose::Pose & pose, string tag) {
	using namespace basic::options;

	vector1<Size> design_pos = get_des_pos(pose);
	vector1<Size> intra_subs; intra_subs.push_back(1); intra_subs.push_back(2); intra_subs.push_back(3);

	// Calculate the surface area and surface complementarity for the interface
	Real int_area = 0; Real sc = 0;
	new_sc(pose, int_area, sc);
	Real avg_deg = average_degree(pose, design_pos);

	Real avg_nchi = 0.0;
	for ( Size i = 1; i <= design_pos.size(); ++i ) {
		avg_nchi += pose.residue(design_pos[i]).nchi();
	}
	avg_nchi /= design_pos.size();

	// Create a scorefile struct, add custom metrics to it
	core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
	ss_out->fill_struct(pose,tag);
	ss_out->add_energy("avg_deg", avg_deg);
	ss_out->add_energy("sc_int_area", int_area);
	ss_out->add_energy("sc", sc);
	ss_out->add_energy("avg_nchi", avg_nchi);

	// // Write the scorefile
	sfd.write_silent_struct( *ss_out, option[OptionKeys::out::file::o]() + "/" + option[ OptionKeys::out::file::silent ]() );


}


int
main (int argc, char *argv[])
{

	try {

		devel::init(argc,argv);
		using namespace basic::options;

		for ( Size ifile = 1; ifile <= option[OptionKeys::in::file::s]().size(); ++ifile ) {
			string fname = option[OptionKeys::in::file::s]()[ifile];
			Pose p;
			core::import_pose::pose_from_file(p,fname, core::import_pose::PDB_file);
			rescore(p,fname);

		}


	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


