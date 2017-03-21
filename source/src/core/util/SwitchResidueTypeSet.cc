// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/util/SwitchResidueTypeSet.cc
/// @brief Functions for switching the residue type set of a pose
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit Headers
#include <core/util/SwitchResidueTypeSet.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/chemical/residue_support.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

#include <core/util/disulfide_util.hh>

// Basic Header
#include <basic/Tracer.hh>

// Option Headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// Utility header
#include <utility/vector1.hh>


namespace core {
namespace util {

static THREAD_LOCAL basic::Tracer TR( "core.util.switchresiduetypeset" );

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details the function allows a pose to use a different residue_type_set to represent all its residues,
///such as from fullatom residues to centroid residues, or vice versa. During the switch, corresponding atoms
///will be copied. Redundant atoms will be removed (in case from fullatom to centroid) and missing atoms will be
///built by ideal geometry (in the case from centroid to fullatom).
void
switch_to_residue_type_set(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & type_set,
	bool allow_sloppy_match,
	bool switch_protein_res_only,
	bool keep_energies
)
{
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace std;
	using utility::vector1;

	//SML 04/06/09
	//Energies object is not properly "aware" of typeset changes, and can attempt to score your pose with an incompatible
	//scorefunction if you go FA->CEN (or vice versa) and access the Energies without rescoring.
	//So, we'll eject the Energies to be safe!
	if ( !keep_energies ) { // but there are some specific cases when we actually want to keep the energies
		pose.energies().clear();
	}

	//fpd only do this over the asymmetric unit if the pose is symmetric
	symmetry::SymmetryInfoCOP symm_info=nullptr;
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		symmetry::SymmetricConformation const & symm_conf ( dynamic_cast< symmetry::SymmetricConformation const & > ( pose.conformation() ) );
		symm_info = symm_conf.Symmetry_Info();
	}

	if ( type_set.mode() == chemical::CENTROID_ROT_t ) {
		switch_to_centroid_rot_set( pose, symm_info, type_set, allow_sloppy_match );
		return;
	}

	bool warned(false); // Have we warned about already-in-this-typeset yet?

	// There's cases where we can create a new ResidueType when switching - cache these to reuse.
	// These need to be ResidueTypeOPs because we're later sticking them into the Pose (all in one go).
	std::map< std::string, core::chemical::ResidueTypeOP > generated_restypes;

	// loop each position and find new type that matches from the new type set
	for ( core::Size i=1; i<= pose.size(); ++i ) {
		if ( !pose.residue( i ).is_protein() && switch_protein_res_only ) continue;
		if ( symm_info && !symm_info->bb_is_independent(i) ) continue;

		core::conformation::Residue const & rsd( pose.residue(i) );
		// in future we may have a conformation using mixed type set, so check this by residue
		if ( rsd.type().mode() == type_set.mode() ) {
			if ( ! warned ) {
				TR.Warning << "When switching to a " << type_set.mode() << " ResidueTypeSet: "
					<< " Pose already contains " << type_set.mode() << " ResidueTypes." << std::endl;
				warned = true;
			}
			TR.Debug << "core::util::switch_to_residue_type_set: residue " << i << " already in a " <<
				type_set.mode() << " ResidueTypeSet" << std::endl;
			continue;
		}

		// get all residue types with same AA

		core::chemical::ResidueTypeCOP new_rsd_type( nullptr );

		if ( ( rsd.aa() == aa_unk ) || ( rsd.name().substr(0,5) == "HIS_D" ) ) {
			// ligand or metal ions are all defined as "UNK" AA, so check a rsdtype with same name
			// for HIS_D tautomer, we want to keep its tautomer state
			if ( type_set.has_name( rsd.name() ) ) {
				new_rsd_type = type_set.name_mapOP( rsd.name() );
			}
		} else  {
			// for a normal AA/DNA/RNA residue, now look for a rsdtype with same variants
			new_rsd_type = type_set.get_representative_type_name3( rsd.name().substr(0,3), rsd.type().variant_types() );
			if ( !new_rsd_type && allow_sloppy_match ) {
				TR.Warning << "Did not find perfect match for residue: "  << rsd.name()
					<< "at position " << i << ". Trying to find acceptable match. " << std::endl;
				core::chemical::ResidueTypeCOP new_rsd_type2( type_set.get_representative_type_name3( rsd.name().substr(0,3) ) );
				if ( new_rsd_type2 && rsd.type().name3()  == new_rsd_type2->name3() ) { // Would the name3's ever not match?
					TR.Warning << "Found an acceptable match: " << rsd.type().name() << " --> " << new_rsd_type2->name() << std::endl;
					new_rsd_type = new_rsd_type2;
				}
			}
		}

		if ( ! new_rsd_type ) {
			if ( generated_restypes.count( rsd.name() ) == 1 ) {
				new_rsd_type = generated_restypes[ rsd.name() ];
			} else {
				core::chemical::ResidueTypeOP gen_restype( generate_replacement_restype( rsd.type(), type_set.mode() ) );
				if ( gen_restype ) {
					TR.Warning << "When switching to " << type_set.mode() << " mode, a normal replacement type for residue " << rsd.name() << " can't be found." << std::endl;
					TR.Warning << "    an autogenerated replacement type is being used instead." << std::endl;
					generated_restypes[ rsd.name() ] = gen_restype;
					new_rsd_type = gen_restype;
				}
			}
		}

		if ( ! new_rsd_type ) {
			TR.Error << pose.sequence() << std::endl;
			TR.Error  << "can not find a residue type that matches the residue " << rsd.name()
				<< " at position " << i << std::endl;
			utility_exit_with_message( "core::util::switch_to_residue_type_set fails\n" );
		}

		core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type, rsd, pose.conformation() ) );

		// switch to corresponding residue type in the new set.
		if ( !rsd.is_protein() ) {
			// rethink this logic, phil
			TR.Debug << "trying to preserve existing coords for non-protein residue: " << rsd.seqpos() << ' ' << rsd.name() << std::endl;
			core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( rsd, *new_rsd, pose.conformation() );
		}
		pose.replace_residue( i, *new_rsd, false );
	}

	// Put the newly-generated ResidueTypes, if any, into the Pose.
	if ( generated_restypes.size() > 0 ) {
		TR.Debug << "Adding " << generated_restypes.size() << " newly generated ResidueTypes to the pose's " << type_set.mode() << " ResidueTypeSet." << std::endl;
		core::chemical::PoseResidueTypeSetOP mod_rts( pose.conformation().modifiable_residue_type_set_for_conf( type_set.mode() ) );
		for ( std::pair< std::string, core::chemical::ResidueTypeOP > const & str_rt_pair : generated_restypes ) {
			if ( ! mod_rts->has_name( str_rt_pair.second->name() ) ) {
				mod_rts->add_unpatchable_residue_type( str_rt_pair.second ); // Unpatchable, as this residue type might already be patched.
			} else {
				TR << "Pose already has a " << str_rt_pair.second->name() << " residue type of mode " << mod_rts->mode() << ": skipping addition." << std::endl;
			}
		}
		pose.conformation().reset_residue_type_set_for_conf(mod_rts);
	}

	// After a CEN->FA transition, rebuild the disulfides
	if ( type_set.mode() == chemical::FULL_ATOM_t && basic::options::option[ basic::options::OptionKeys::run::rebuild_disulf ]() ) {
		rebuild_fa_disulfides(pose);
	}
}

void
switch_to_residue_type_set(
	core::pose::Pose & pose,
	core::chemical::TypeSetMode type_set_mode,
	bool allow_sloppy_match,
	bool switch_protein_res_only,
	bool keep_energies
) {
	core::chemical::ResidueTypeSetCOP type_set( pose.residue_type_set_for_pose( type_set_mode ) );
	switch_to_residue_type_set( pose, *type_set, allow_sloppy_match, switch_protein_res_only, keep_energies );
}

void
switch_to_residue_type_set(
	core::pose::Pose & pose,
	std::string const & type_set_name,
	bool allow_sloppy_match,
	bool switch_protein_res_only,
	bool keep_energies
) {
	core::chemical::ResidueTypeSetCOP type_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( type_set_name ) );
	switch_to_residue_type_set( pose, *type_set, allow_sloppy_match, switch_protein_res_only, keep_energies );
}

void
switch_to_centroid_rot_set(
	core::pose::Pose & pose,
	core::conformation::symmetry::SymmetryInfoCOP symm_info,
	core::chemical::ResidueTypeSet const & rsd_set,
	bool /* allow_sloppy_match */
) {
	using namespace core::chemical;
	using namespace core::conformation;
	using utility::vector1;

	debug_assert( rsd_set.mode() == chemical::CENTROID_ROT_t );

	bool warned(false); // Have we warned about already-in-this-typeset yet?

	//loop for each residue
	for ( core::Size i=1; i<= pose.total_residue(); ++i ) {
		if ( symm_info && !symm_info->bb_is_independent(i) ) continue;

		//get the residue
		Residue const & rsd( pose.residue(i) );
		if ( (rsd.aa()==aa_unk) || (rsd.aa()==aa_vrt) ) continue; //skip invalid residue

		TypeSetMode existing_rsd_typeset_mode( rsd.type().mode() );
		//check current restype
		if ( existing_rsd_typeset_mode == rsd_set.mode() ) {
			if ( ! warned ) {
				TR.Warning << "When switching to a " << rsd_set.mode() << " ResidueTypeSet: "
					<< " Pose already contains " << existing_rsd_typeset_mode << " ResidueTypes." << std::endl;
				warned = true;
			}
			TR.Debug << "core::util::switch_to_residue_type_set: residue " << i << " already in a " <<
				existing_rsd_typeset_mode << " ResidueTypeSet" << std::endl;
			continue;
		}

		//get temperature
		core::Real maxB=0.0;
		if ( pose.pdb_info() ) {
			for ( Size k=rsd.first_sidechain_atom(); k<=rsd.nheavyatoms(); ++k ) {
				if ( rsd.is_virtual(k) ) continue;
				Real B = pose.pdb_info()->temperature( i, k );
				maxB = std::max( B, maxB );
			}
			if ( maxB==0.0 ) { maxB = pose.pdb_info()->temperature( i, rsd.atom_index("CA") );}

			//TR.Debug << "maxB=" << maxB << std::endl;
		}

		//gen new residue
		core::conformation::ResidueOP new_rsd( 0 );
		if ( (rsd.aa()==aa_unk) || (rsd.aa()==aa_vrt) ) {
			//skip
			continue;
		} else if ( rsd.name().substr(0,5)=="HIS_D" ) {
			//else if ( rsd.name().substr(0,3)=="CYD" ) {
			// core::chemical::ResidueTypeCOPs const & rsd_types( rsd_set->name3_map( "CYS" ) );
			// for (core::Size j=1; j<=rsd_types.size(); ++j ) {
			//  core::chemical::ResidueType const & new_rsd_type( *rsd_types[j] );
			//  if ( new_rsd_type.name3()=="CYS" ) {
			//   new_rsd = core::conformation::ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );
			//   break;
			//  }
			// }
			//}
			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set.get_representative_type_name3( rsd.name3() ) );
			if ( new_rsd_type && new_rsd_type->name3() == "HIS" ) { // Is it ever the case that it wouldn't be "HIS"?
				new_rsd = core::conformation::ResidueFactory::create_residue( *new_rsd_type, rsd, pose.conformation() );
			}
		} else  if ( rsd.is_terminus() ) {
			//get the terminal type (maybe no need, but to consist with reading a cenrot pdb)
			//TR.Debug << "TER" << std::endl;
			core::chemical::ResidueType const & new_rsd_type( rsd_set.name_map(rsd.name()) );
			new_rsd = core::conformation::ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );
		} else {
			//just find the standard aa restype
			//TR.Debug << "looking for " << rsd.name() << std::endl;
			core::chemical::ResidueType const & new_rsd_type( rsd_set.name_map(rsd.name()) );
			new_rsd = core::conformation::ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );
		}

		if ( ! new_rsd ) {
			TR.Warning << "Did not find perfect match for residue: "  << rsd.name()
				<< " at position " << i << ". Trying to find acceptable match. " << std::endl;

			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set.get_representative_type_name3( rsd.name3() ) );
			if ( new_rsd_type && new_rsd_type->name3() == rsd.name3() ) { // Is it ever the case the name3's wouldn't match?
				new_rsd = core::conformation::ResidueFactory::create_residue( *new_rsd_type, rsd, pose.conformation() );
			}
			if ( new_rsd ) {
				TR.Warning << "Found an acceptable match: " << rsd.type().name() << " --> " << new_rsd->name() << std::endl;
			} else {
				//bug here?
				utility_exit_with_message( "switch to "+core::chemical::string_from_type_set_mode(rsd_set.mode())+" type set fails\n" );
			}
		}

		//find the centroid postion based on fa sidechain
		PointPosition cenrotxyz(0,0,0);
		PointPosition cbxyz(0,0,0);
		if ( existing_rsd_typeset_mode == chemical::FULL_ATOM_t ) {
			std::map<std::string, Real> masslst;
			masslst["C"]=12.0107;
			masslst["O"]=15.9994;
			masslst["S"]=32.066;
			masslst["N"]=14.00674;
			masslst["H"]=1.00794;

			Real mass = 0.0;
			if ( rsd.name3()=="GLY" ) {
				cenrotxyz = rsd.atom("CA").xyz();
			} else if ( rsd.name3()=="ALA" ) {
				cenrotxyz = rsd.atom("CB").xyz();
			} else {
				//std::cout<<rsd.name()<<std::endl;
				for ( Size na=rsd.type().first_sidechain_atom();
						na<=rsd.type().nheavyatoms();
						na++ ) {
					if ( rsd.atom_name(na)==" CB " ) continue;
					std::string elem = rsd.atom_name(na).substr(1,1);
					cenrotxyz += (rsd.atoms()[na].xyz()*masslst[elem]);
					mass += masslst[elem];
					// TR.Debug << "|" << rsd.atom_name(na) << "|"
					// << " (" << masslst[elem] << ") "
					// << rsd.atoms()[na].xyz().x() << ","
					// << rsd.atoms()[na].xyz().y() << ","
					// << rsd.atoms()[na].xyz().z() << std::endl;
				}
				cenrotxyz = cenrotxyz/mass;
			}
		} else if ( existing_rsd_typeset_mode == chemical::CENTROID_t && rsd.name()!="VRT" ) {
			//keep the cen position
			cenrotxyz = rsd.atom("CEN").xyz();

			// preserve CB
			if ( !(rsd.name3()=="GLY") ) {
				cbxyz = rsd.atom("CB").xyz();
			}
		}

		//replace
		if ( ! new_rsd ) {
			std::cerr << pose.sequence() << std::endl;
			std::cerr  << "can not find a residue type that matches the residue " << rsd.name()
				<< " at position " << i << std::endl;
			utility_exit_with_message( "switch_to_cenrot_residue_type_set fails\n" );
			//continue;
		}

		//replace it
		pose.replace_residue( i, *new_rsd, false );
		if ( existing_rsd_typeset_mode == chemical::FULL_ATOM_t ) {
			//set centroid_rot xyz
			pose.set_xyz(id::AtomID(pose.residue(i).atom_index("CEN"), i), cenrotxyz);
		} else if ( existing_rsd_typeset_mode == chemical::CENTROID_t && new_rsd->name()!="VRT" ) {
			// preserve CB
			if ( !(new_rsd->name3()=="GLY") ) {
				pose.set_xyz(id::AtomID(pose.residue(i).atom_index("CB"), i), cbxyz);
			}

			pose.set_xyz(id::AtomID(pose.residue(i).atom_index("CEN"), i), cenrotxyz);
		}

		//set temperature
		if ( pose.pdb_info() ) pose.pdb_info()->temperature(i, pose.residue(i).nbr_atom(), maxB);
	}

	return;
}

core::chemical::ResidueTypeOP
generate_replacement_restype(
	core::chemical::ResidueType const & restype,
	core::chemical::TypeSetMode mode)
{
	switch( restype.mode() ) {
	case core::chemical::FULL_ATOM_t :
		switch( mode ) {
		case core::chemical::CENTROID_t :
			// From FULL_ATOM to CENTROID
			return core::chemical::make_centroid( restype );
		default :
			return nullptr;
		}
	default :
		return nullptr;
	}
}

void
rebuild_fa_disulfides(
	core::pose::Pose & pose
) {
	debug_assert( pose.is_fullatom() );

	utility::vector1< std::pair<core::Size,core::Size> > disulfides;
	core::conformation::disulfide_bonds(pose.conformation(), disulfides);

	if ( disulfides.size() > 0 ) {
		// Setup Packer & Minimizer
		core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
		task->initialize_from_command_line().or_include_current( true );
		task->restrict_to_repacking();

		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb( false );

		// Set up each residue individually
		for ( core::Size i(1); i <= pose.size(); ++i ) {
			core::conformation::Residue const& res(pose.residue(i));
			if ( !res.is_protein() ) {
				continue;
			}

			// Determine if i is part of disulfides
			bool is_disulf = false;
			for ( utility::vector1< std::pair<core::Size, core::Size> >::const_iterator
					disulf(disulfides.begin()), end_disulf(disulfides.end());
					disulf != end_disulf; ++disulf ) {
				if ( i == disulf->first || i == disulf->second ) {
					is_disulf = true;
					break;
				}
			}

			if ( is_disulf ) {
				// repack & minimize disulfides
				mm->set_chi(i, true);
			} else {
				// Other residues are unchanged
				task->nonconst_residue_task(i).prevent_repacking();
			}
		}

		// Rebuild disulfides
		core::util::rebuild_disulfide(pose, disulfides, task, nullptr, mm, nullptr);
	}
}

} // util
} // core
