// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/toolbox/match_enzdes_util/util_functions.cc
/// @brief bunch of utility functions
/// @author Florian Richter, floric@u.washington.edu
/// @modified Tom Linsky, tlinsky@uw.edu

//unit headers
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>

//package headers
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>

// project headers
#include <basic/Tracer.hh>

#include <core/types.hh>

#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/BackboneStubConstraint.hh>
#include <core/scoring/constraints/MultiConstraint.hh>


//utility headers
#include <utility/string_util.hh>

//stl headers
#include <map>

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

static basic::Tracer tr( "protocols.toolbox.match_enzdes_util.util_functions" );

void
replace_residue_keeping_all_atom_positions(
	core::pose::Pose & pose,
	core::conformation::Residue new_res,
	core::Size res_pos
)
{

	//have to set the position of the new res to their old values, so we gotta save them now
	std::map< std::string, core::PointPosition > atom_name_to_xyz;

	for ( core::Size at_ct = 1; at_ct <= pose.residue(res_pos).natoms(); at_ct++ ) {
		// If new_res lacks this atom, check if it has a spare virt.
		atom_name_to_xyz.insert(  std::pair< std::string, core::PointPosition > (pose.residue(res_pos).atom_name(at_ct), pose.residue(res_pos).xyz( at_ct ) ) );
	}

	//replacing the residue
	pose.replace_residue( res_pos, new_res, true);

	//and resetting the xyz positions

	// We may have new virtual atoms (or hydrogens) we need to place.
	core::id::AtomID_Mask missing(false);

	for ( core::Size at_ct = 1; at_ct <= pose.residue(res_pos).natoms(); at_ct++ ) {

		auto xyz_map_it = atom_name_to_xyz.find( pose.residue(res_pos).atom_name(at_ct) );


		if ( xyz_map_it != atom_name_to_xyz.end() ) {
			pose.set_xyz( core::id::AtomID(at_ct, res_pos), xyz_map_it->second );
		} else if ( pose.residue(res_pos).is_virtual(at_ct) || pose.residue(res_pos).atom_is_hydrogen(at_ct) ) {
			missing.set( core::id::AtomID(at_ct, res_pos), true );
		} else {
			std::cerr << "ERROR: when trying to make dsflkj constraint covalent, atom " << pose.residue(res_pos).atom_name(at_ct) << " was not found for residue " << pose.residue(res_pos).name3() << " at position " << res_pos << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
	}

	if ( ! missing.empty() ) {
		pose.conformation().fill_missing_atoms( missing );
	}

} //replace_residues_keeping_positions


/// @details forwarding function, allows stuff to call this functionality
///          without caller having to specify root
core::scoring::constraints::AmbiguousConstraintCOP
constrain_pose_res_to_invrots(
	std::list< core::conformation::ResidueCOP > const & invrots,
	utility::vector1< core::Size > const & seqpos,
	core::pose::Pose const & pose,
	core::scoring::func::FuncOP constraint_func
)
{
	core::id::AtomID fixed_pt( pose.atom_tree().root()->atom_id() );
	return constrain_pose_res_to_invrots( invrots, seqpos, pose, fixed_pt, constraint_func );
}


core::scoring::constraints::AmbiguousConstraintCOP
constrain_pose_res_to_invrots(
	std::list< core::conformation::ResidueCOP > const & invrots,
	utility::vector1< core::Size > const & seqpos,
	core::pose::Pose const & pose,
	core::id::AtomID const & fixed_pt,
	core::scoring::func::FuncOP constraint_func
)
{
	using namespace core::scoring::constraints;

	if ( !constraint_func ) constraint_func = core::scoring::func::FuncOP( new BoundFunc( 0, 0.05, 0.4, "invrot") );
	//see the comment in protocols/ligand_docking/LigandBaseProtocol.cc::restrain_protein_Calphas
	//core::id::AtomID fixed_pt( pose.atom_tree().root()->atom_id() );
	//tr << "Hackack fixed_pt was passed in to be AtomID " << fixed_pt << std::endl;

	utility::vector1< ConstraintCOP > all_res_invrot_csts;
	core::Size totrescount(0);

	for ( core::Size i =1; i <= seqpos.size(); ++i ) {

		core::conformation::ResidueCOP cur_remodel_res( pose.residue( seqpos[i] ).get_self_ptr() );
		if ( cur_remodel_res->name3() == "GLY" ) continue;
		totrescount++;

		core::id::AtomID rem_CA( cur_remodel_res->type().atom_index("CA"), seqpos[i] );
		core::id::AtomID rem_CB( cur_remodel_res->type().atom_index("CB"), seqpos[i] );
		core::id::AtomID rem_N( cur_remodel_res->type().atom_index("N"), seqpos[i] );

		for ( auto const & invrot : invrots ) {

			utility::vector1< ConstraintCOP > cur_res_invrot_csts;
			cur_res_invrot_csts.push_back( core::scoring::constraints::ConstraintOP( new BackboneStubConstraint( pose, seqpos[i], fixed_pt, *invrot, -20.0, 0.8) ) );

			//old style: coordinate constraints for all atoms, backbone stub csts
			// might be working better
			cur_res_invrot_csts.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( rem_CA, fixed_pt, invrot->xyz("CA"), constraint_func ) ) );
			cur_res_invrot_csts.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( rem_CB, fixed_pt, invrot->xyz("CB"), constraint_func ) ) );
			cur_res_invrot_csts.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint( rem_N, fixed_pt, invrot->xyz("N"), constraint_func ) ) );

			all_res_invrot_csts.push_back( core::scoring::constraints::ConstraintOP( new MultiConstraint( cur_res_invrot_csts ) ) );
		}// loop over invrots
	}//loop over seqpos

	tr << "Created a total of " << all_res_invrot_csts.size() << " constraints between " << invrots.size() << " inverse rotamers and " << totrescount << " residues." << std::endl;

	return core::scoring::constraints::AmbiguousConstraintCOP( core::scoring::constraints::AmbiguousConstraintOP( new AmbiguousConstraint( all_res_invrot_csts ) ) );

} //constrain_pose_res_to_invrots


core::conformation::ResidueCOP
cst_residue_in_pose(
	core::pose::Pose const & pose,
	core::Size geomcst,
	core::Size geomcst_template_res
)
{

	EnzdesCacheableObserverCOP enz_ob( get_enzdes_observer( pose ) );
	if ( !enz_ob ) return nullptr;

	EnzCstTemplateResCacheCOP res_cache( enz_ob->cst_cache()->param_cache( geomcst )->template_res_cache( geomcst_template_res) );
	runtime_assert( res_cache != nullptr );
	if ( res_cache->not_in_pose() ) return nullptr;

	runtime_assert( res_cache->seqpos_map_size() == 1 );

	return core::conformation::ResidueCOP( core::conformation::ResidueOP( new core::conformation::Residue( pose.residue( res_cache->seqpos_map_begin()->first ) ) ) );
}


std::string
assemble_remark_line(
	std::string const & chainA,
	std::string const & resA,
	int seqposA,
	std::string const & chainB,
	std::string const & resB,
	int seqposB,
	core::Size cst_block,
	core::Size ex_geom_id
)
{

	std::string posA = utility::pad_left( seqposA, 4 );

	std::string posB = utility::pad_left( seqposB, 4 );

	return "MATCH TEMPLATE "+ chainA +" "+ resA +" "+ posA +  " MATCH MOTIF "+ chainB + " " + resB + " "+posB + "  " + utility::to_string( cst_block ) + "  " + utility::to_string( ex_geom_id );

} //assemble remark line function


bool
split_up_remark_line(
	std::string line,
	std::string & chainA,
	std::string & resA,
	int & seqposA,
	std::string & chainB,
	std::string & resB,
	int & seqposB,
	core::Size & cst_block,
	core::Size & ex_geom_id
)
{

	std::istringstream line_stream;
	std::string buffer(""), tag("");

	line_stream.clear();
	line_stream.str( line );

	line_stream >> buffer >> tag;
	if ( tag == "TEMPLATE" ) {
		line_stream >> chainA >> resA >> seqposA >> buffer >> buffer;
		line_stream >> chainB >> resB >> seqposB >> cst_block;
		if ( resA.size() == 2 ) resA = " " + resA;
		if ( resB.size() == 2 ) resB = " " + resB;

		if ( line_stream.bad() ) {
			tr << "ERROR when trying to split up pdb remark line. Not all fields seem to have been specified." << std::endl;
			return false;
		}

		line_stream >> ex_geom_id;
		if ( !line_stream.good() ) ex_geom_id = 1;

		return true;
	}

	return false;
}  //split up remark line function

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// @brief finds the first non-ligand residue in the pose  (should be the N-terminus)
core::Size
get_first_protein_residue( core::pose::Pose const & pose )
{
	// This will not work properly on structures with multiple chains, for now

	// briefly, loop through the list of residues in the pose, from 1 to N
	// if the pose is not in the list of ligands, it is the N-terminal residue
	for ( core::Size i = 1; i <= pose.size(); i++ ) {
		if ( pose.residue( i ).is_protein() ) {
			return i;
		}
	}

	tr.Fatal << "No non-ligand residues were detected!!" << std::endl;
	utility_exit_with_message("No non-ligand residues were detected in pose.");
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// @brief finds the first non-ligand residue in the pose  (should be the N-terminus)
core::Size
get_last_protein_residue( core::pose::Pose const & pose )
{
	// This should fail on structures with multiple chains, for now

	// briefly, loop through the list of residues in the pose, from 1 to N
	// if the pose is not in the list of ligands, it is the N-terminal residue
	for ( core::Size i = pose.size(); i >= 1; i-- ) {
		if ( pose.residue( i ).is_protein() ) {
			return i;
		}
	}

	tr.Fatal << "No non-ligand residues were detected!!" << std::endl;
	utility_exit_with_message("No non-ligand residues were detected in pose.");
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< core::chemical::ResidueTypeCOP >
sort_residue_type_pointers_by_name( utility::vector1< core::chemical::ResidueTypeCOP > const & restype_temp_set )
{
	utility::vector1< std::pair< std::string, core::chemical::ResidueTypeCOP > > restype_temp_set_with_name;
	for ( core::Size n = 1; n <= restype_temp_set.size(); n++ ) restype_temp_set_with_name.push_back( std::make_pair( restype_temp_set[ n ]->name(), restype_temp_set[ n ] ) );
	std::sort( restype_temp_set_with_name.begin(), restype_temp_set_with_name.end() );
	// for ( Size n = 1; n <= restype_temp_set.size(); n++ ) tr << n << " original order: " << restype_temp_set[ n ]->name() << " new order " <<  (restype_temp_set_with_name[ n ].second)->name() << std::endl;

	//finally put all the restypes into the storage vector
	utility::vector1< core::chemical::ResidueTypeCOP > restype_set_sorted;
	for ( auto & set_it : restype_temp_set_with_name ) {
		restype_set_sorted.push_back( set_it.second );
	}

	return restype_set_sorted;
}


}  // match_enzdes_util
} // toolbox
} //protocols
