// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/conformation/symmetry/util.hh
/// @brief utility functions for handling of symmetric conformations
/// @author Ingemar Andre

// Unit headers
#include <core/conformation/symmetry/util.hh>

// Package headers
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/MirrorSymmetricConformation.hh>
#include <core/conformation/symmetry/VirtualCoordinate.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/util.hh>

#include <core/kinematics/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/kinematics/Edge.hh>
#include <core/id/AtomID.hh>

// Basic Headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>


namespace core {
namespace conformation {
namespace symmetry {



//This string is duplicated in utility/tag/XMLSchemaGeneration.cc; DO NOT MODIFY!
static std::string const chr_chains("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcba");

std::string
get_chain_id_string(){
	return chr_chains;
}

static basic::Tracer TR( "core.conformation.symmetry.util" );


/// @details A very specific helper function that counts the number of mirror ops from root->each subunit
/// returns:
///   a bool vector 'mirrored_subs' that defines which subunits are mirrored w.r.t. master
///   a bool pair vector 'mirrored_jumps' that defines which jumps are mirrored w.r.t. master (upstream/downstream)
void
calculate_inverting_virtuals(
	core::kinematics::FoldTree const & ft,
	Conformation const & conf,  // unsymmetrized conf
	SymmetryInfo const &symm_info,
	utility::vector1<bool> &mirrored_subs,
	utility::vector1< std::pair<bool,bool> > &inverted_jumps )
{
	// find all current_res -> other res jumps
	for ( auto const & it : ft ) {
		if ( it.is_polymer() ) continue;

		int const upstream (it.start());
		int const downstream (it.stop());

		// if downstream res is non-vrt, mark the subunit array
		if (
				downstream <= (int)symm_info.num_total_residues_without_pseudo() &&
				upstream > (int)symm_info.num_total_residues_without_pseudo()
				) {
			// note: this logic neets to be changed for two-comp with different inversion
			//       within the ASU
			mirrored_subs[ symm_info.subunit_index( downstream ) ] =
				conf.residue_type(upstream).is_inverted_virtual_residue();
			inverted_jumps[ it.label() ] = std::make_pair<bool,bool>(
				conf.residue_type(upstream).is_inverted_virtual_residue(),
				conf.residue_type(upstream).is_inverted_virtual_residue() );
		} else {
			inverted_jumps[ it.label() ] = std::make_pair<bool,bool> (
				conf.residue_type(upstream).is_inverted_virtual_residue(),
				conf.residue_type(downstream).is_inverted_virtual_residue() );
		}
	}

	// now mark all intra-res jumps
	for ( auto const & it : ft ) {
		if ( it.is_polymer() ) continue;

		int const upstream (it.start());
		int const downstream (it.stop());

		// if downstream res is non-vrt, mark the subunit array
		if (
				downstream <= (int)symm_info.num_total_residues_without_pseudo() &&
				upstream <= (int)symm_info.num_total_residues_without_pseudo()
				) {
			inverted_jumps[ it.label() ] = std::make_pair<bool,bool>(
				mirrored_subs[ symm_info.subunit_index(upstream) ],
				mirrored_subs[ symm_info.subunit_index(downstream) ] );
		}
	}
}



Size fold_tree_entry_point(core::kinematics::FoldTree const & ft, Size lb_resi=0, Size ub_resi=0) {
	int resi = -1;
	if ( lb_resi==0 ) lb_resi = 1;
	if ( ub_resi==0 ) ub_resi = ft.nres();
	for ( int i = 1; i <= (int)ft.num_jump(); ++i ) {
		Size up = ft.upstream_jump_residue(i);
		Size dn = ft.downstream_jump_residue(i);
		if ( lb_resi <= dn && dn <= ub_resi && ( up < lb_resi || ub_resi < up ) ) resi = dn; // dn in bounds, up out
	}
	if ( (int)lb_resi <= resi && resi <= (int)ub_resi ) return resi; // if resi good, return it
	resi = ft.root();
	if ( (int)lb_resi <= resi && resi <= (int)ub_resi ) return resi; // try fold_tree root
	utility_exit_with_message("can'd get fold_tree root in range!!!");
	return 0; // this will never happen
}

// sheffler
/// @details get residue index
Size process_residue_request(conformation::Conformation const & src_conf, std::string input, Size lb_resi=0, Size ub_resi=0) {
	int resi = -1;
	if ( lb_resi==0 ) lb_resi = 1;
	if ( ub_resi==0 ) ub_resi = src_conf.size();
	if ( input=="FT"    || input=="KEEP_FOLDTREE_ANCHOR" ) resi = fold_tree_entry_point( src_conf.fold_tree(), lb_resi, ub_resi );
	if ( input=="COM"   || input=="CENTER_OF_MASS"      ) resi = residue_center_of_mass( src_conf, lb_resi, ub_resi );
	if ( input=="FIRST" || input=="FIRST_RESIDUE"       ) resi = lb_resi;
	if ( input=="LAST"  || input=="LAST_RESIDUE"        ) resi = ub_resi;
	if ( input=="" ) resi = 0;
	if ( resi == -1 ) resi = utility::string2int(input);
	if ( resi < 0 ) utility_exit_with_message("error in process_residue_request: '"+input+"'");
	// TR << "process_residue_request: '" << input << "' to " << resi << std::endl;
	return resi;
}

bool is_jump_intracomponent(std::map<char,std::pair<Size,Size> > chain2range, Size up, Size dn) {
	for ( std::map<char,std::pair<Size,Size> >::const_iterator it = chain2range.begin(); it != chain2range.end(); ++it ) {
		Size const lb = it->second.first, ub = it->second.second;
		if ( lb <= up && up <= ub && lb <= dn && dn <= ub ) return true;
	}
	return false;
}

char which_component(std::map<char,std::pair<Size,Size> > chain2range, Size resi) {
	for ( std::map<char,std::pair<Size,Size> >::const_iterator it = chain2range.begin(); it != chain2range.end(); ++it ) {
		Size const lb = it->second.first, ub = it->second.second;
		if ( lb <= resi && resi <= ub ) return it->first;
	}
	utility_exit_with_message("resi not in any component!");
	return (char)NULL; // this will never happen
}

core::kinematics::FoldTree
get_component_contiguous_foldtree(
	core::kinematics::FoldTree const & f_orig,
	std::map<char,std::pair<Size,Size> > const & /*chain2range*/ // maybe use later
) {
	if ( f_orig.num_cutpoint() != f_orig.num_jump() ) utility_exit_with_message("FoldTree ISANITY!!!!!!!");
	ObjexxFCL::FArray1D< Size > cuts( f_orig.num_cutpoint() );
	ObjexxFCL::FArray2D< Size > jumps( 2, f_orig.num_jump() );
	// utility::vector1<std::string> upatom(f_orig.num_jump() );
	utility::vector1<std::string> dnatom(f_orig.num_jump() );

	utility::vector1< Size > cutpoints(f_orig.num_cutpoint());
	for ( Size i = 1; i <= f_orig.num_cutpoint(); ++i ) {
		cutpoints[i] = f_orig.cutpoint(i);
	}
	std::sort(cutpoints.begin(),cutpoints.end());
	cutpoints.push_back(f_orig.nres());

	// simple way for now
	for ( Size i = 1; i < cutpoints.size(); ++i ) {
		cuts(i) = cutpoints[i];
		jumps(1,i) = cutpoints[i];
		jumps(2,i) = fold_tree_entry_point(f_orig, cutpoints[i]+1,cutpoints[i+1]);
		// upatom[i] = f_orig.upstream_atom(i);
		dnatom[i] = f_orig.downstream_atom(i);
	}
	core::kinematics::FoldTree f_contig;
	f_contig.tree_from_jumps_and_cuts(f_orig.nres(),f_orig.num_jump(),jumps,cuts);
	for ( Size i = 1; i < cutpoints.size(); ++i ) {
		if ( dnatom[i]!="" ) f_contig.set_jump_atoms(i,"N",dnatom[i]);
	}
	f_contig.reorder(1);

	return f_contig;
}


// sheffler
/// @details get a mapping of chain chars to resi ranges
std::map<char,std::pair<Size,Size> >
get_chain2range( Conformation const & src_conf, std::map< core::Size, char > src_conf2pdb_chain ) {
	// TR << src_conf.num_chains() << endl;
	// for(Size i = 1; i <= src_conf.num_chains(); ++i){
	//  TR << "src_conf chain " << i << " " << src_conf.chain_begin(i) << "-" << src_conf.chain_end  (i) << endl;
	// }
	// TR << "PDBINFO CHAIN MAP" << endl;
	// for(map< core::Size, char >::const_iterator i = src_conf2pdb_chain.begin(); i != src_conf2pdb_chain.end(); ++i) {
	//  TR << i->first << " " << i->second << endl;
	// }
	// TR << "END PDBINFO CHAIN MAP" << endl;
	std::map<char,std::pair<Size,Size> > crange;
	for ( Size i = 1; i <= src_conf.num_chains(); ++i ) {
		char chain = src_conf2pdb_chain.count(i)==0 ? chr_chains[(i-1)%chr_chains.size()] : src_conf2pdb_chain[i];
		if ( crange.count(chain)==0 ) crange[chain] = std::make_pair(99999999,0);
		crange[chain].first  = std::min( src_conf.chain_begin(i), crange[chain].first  );
		crange[chain].second = std::max( src_conf.chain_end  (i), crange[chain].second );
	}
	// sanity check
	for ( std::map<char,std::pair<Size,Size> >::const_iterator i = crange.begin(); i != crange.end(); ++i ) {
		// TR << "get_chain2range " << i->first << " " << i->second.first << "-" << i->second.second << std::endl;
		if ( i->second.first > i->second.second ) utility_exit_with_message("THIS SHOULD NEVER HAPPEN!");
		for ( std::map<char,std::pair<Size,Size> >::const_iterator j = crange.begin(); j != crange.end(); ++j ) {
			if ( i->second.first==j->second.first && i->second.second==j->second.second ) continue; // same
			if ( i->second.first > j->second.first && i->second.first > j->second.second ) continue; // disjoint, i < j
			if ( j->second.first > i->second.first && j->second.first > i->second.second ) continue; // disjoint, j > i
			utility_exit_with_message("[ERROR] overlapping chain ranges in pose!!!!!!!!");
		}
	}
	return crange;
}

/// @details  Attempt to detect whether a conformation is symmetric
bool
is_symmetric( conformation::Conformation const & conf )
{
	return ( dynamic_cast< conformation::symmetry::SymmetricConformation const * >( &conf ) );
}

/// @details  Attempt to detect whether a conformation is mirror symmetric
bool
is_mirror_symmetric( conformation::Conformation const & conf )
{
	return ( dynamic_cast< conformation::symmetry::MirrorSymmetricConformation const * >( &conf ) );
}

bool
is_symmetric( conformation::symmetry::SymmetryInfo const & symminfo )
{
	return symminfo.get_use_symmetry();
}


/// @details Generate a symmetric conformation from a monomeric pose and symmetry information
/// stored in the SymmData object
conformation::symmetry::SymmetricConformationOP
setup_symmetric_conformation(
	conformation::Conformation & src_conformation,
	conformation::symmetry::SymmData & symmdata,
	std::map< core::Size, char > src_conf2pdb_chain
)
{
	if ( symmdata.get_num_components() > 1 ) {
		std::map<char,std::pair<Size,Size> > const chain2range = get_chain2range(src_conformation,src_conf2pdb_chain);
		core::kinematics::FoldTree newft = get_component_contiguous_foldtree(src_conformation.fold_tree(),chain2range);
		src_conformation.fold_tree( newft );
	}

	Size njump_orig = src_conformation.fold_tree().num_jump();

	// maybe a little inefficient: first build a standard conformation, then construct symmetric conformation at the end
	conformation::ConformationOP conf_op = src_conformation.clone();
	conformation::Conformation & conf = *conf_op;

	// recenter the pose to the origin?
	if ( symmdata.get_recenter() ) {
		core::conformation::symmetry::recenter( conf, symmdata );
	}

	// Setup temporary variables
	Size const nres_monomer( conf.size() );
	Size const njump_monomer( conf.fold_tree().num_jump() );

	// setup the pseudo residues
	std::map< std::string, conformation::symmetry::VirtualCoordinate > coords ( symmdata.get_virtual_coordinates() );
	std::map< std::string, std::pair< std::string, std::string > > virtual_connects( symmdata.get_virtual_connects() );
	std::map< Size, std::string > virtual_num_to_id (symmdata.get_virtual_num_to_id() );

	// Setup virtual residues
	bool has_mirror_operations(false);
	{
		// create the new residue
		chemical::ResidueTypeSetCOP rsd_set( src_conformation.residue_type_set_for_conf() );

		conformation::ResidueOP rsd( conformation::ResidueFactory::create_residue( *conformation::virtual_type_for_conf( src_conformation ) ) );
		conformation::ResidueOP irsd( conformation::ResidueFactory::create_residue( *conformation::inv_virtual_type_for_conf( src_conformation ) ) );

		for ( Size i =1; i <= virtual_num_to_id.size(); ++i ) {
			if ( virtual_num_to_id.find( i ) == virtual_num_to_id.end() ) {
				utility_exit_with_message( "[ERROR] Cannot find jump number..." );
			}
			std::string tag ( virtual_num_to_id.find( i )->second );
			if ( coords.find( tag ) == coords.end() ) {
				utility_exit_with_message( "[ERROR] Cannot find tag " + tag );
			}
			conformation::symmetry::VirtualCoordinate virt_coord( coords.find( tag )->second );
			numeric::xyzVector < core::Real > orig( virt_coord.get_origin() );
			if ( virt_coord.get_mirror_z() ) {
				has_mirror_operations = true; //There is at least one mirror operation in the set of mirror operations.
				irsd->set_xyz( "ORIG", orig );
				irsd->set_xyz( "X", virt_coord.get_x().normalized() + orig );
				irsd->set_xyz( "Y", virt_coord.get_y().normalized() + orig );
				conf.append_residue_by_jump( *irsd, conf.size() ); //append it to the end of the monomer
			} else {
				rsd->set_xyz( "ORIG", orig );
				rsd->set_xyz( "X", virt_coord.get_x().normalized() + virt_coord.get_origin() );
				rsd->set_xyz( "Y", virt_coord.get_y().normalized() + virt_coord.get_origin() );
				conf.append_residue_by_jump( *rsd, conf.size() ); //append it to the end of the monomer
			}
		}
	}

	// now insert the other conformations. This step will create a fold_tree that will be discarded at a later stage.
	// Optimally we want to set the fold tree directly but that turns out to be difficult. TODO change
	// kinematics::Jump const monomer_jump( conf.jump( njump_monomer+1 ) );
	for ( Size i=1; i< symmdata.get_subunits(); ++i ) {
		Size const insert_seqpos( nres_monomer * i + 1 ); // desired sequence position of residue 1 in inserted conformation
		Size const insert_jumppos( njump_monomer * i + 1 ); // desired jump number of 1st jump in inserted conformation
		Size const anchor_pseudo( nres_monomer * i + i + 1 ); // in current numbering
		Size const new_jump_number( njump_monomer*( i + 1 ) + i + 1 );
		conf.insert_conformation_by_jump( src_conformation, insert_seqpos, insert_jumppos, anchor_pseudo, new_jump_number);
	}

	core::kinematics::visualize_fold_tree(conf.fold_tree());

	// Now generate the fold_tree from the symmdata and set it into the conformation. In the future we want to be able to use a
	// atom tree for this.
	kinematics::FoldTree f ( core::conformation::symmetry::set_fold_tree_from_symm_data( src_conformation, symmdata, src_conf2pdb_chain ) );

	// set the root of the fold tree
	Size new_root ( conf.size() - coords.size() + symmdata.get_root() );
	f.reorder( new_root );
	TR.Debug << f << std::endl;
	conf.fold_tree( f );

	// now build the symmetry info object
	conformation::symmetry::SymmetryInfo symm_info_raw( symmdata, nres_monomer, njump_monomer+1-symmdata.get_num_components() );
	std::map<char,std::pair<Size,Size> > component_bounds;
	for ( Size ic = 1; ic <= src_conformation.num_chains(); ++ic ) {
		char chain = src_conf2pdb_chain[ic];
		component_bounds[chain].first = src_conformation.chain_begin(ic);
		component_bounds[chain].second = src_conformation.chain_end(ic);
	}
	symm_info_raw.set_multicomponent_info(
		symmdata.get_num_components(),
		symmdata.get_components(),
		component_bounds,
		symmdata.get_subunit_name_to_component(),
		symmdata.get_jump_name_to_components(),
		symmdata.get_jump_name_to_subunits()
	);
	conformation::symmetry::SymmetryInfo const & symm_info(symm_info_raw);

	if ( symmdata.get_num_components() > 1 ) {
		for ( std::map< Size, SymDof >::const_iterator j = symm_info.get_dofs().begin(); j != symm_info.get_dofs().end(); ++j ) {
			std::string const & dofname( symm_info.get_jump_name(j->first) );
			utility::vector1<char> compchild = symmdata.components_moved_by_jump(dofname);
			utility::vector1<Size> subchild = symmdata.subunits_moved_by_jump(dofname);
			TR << "MULTICOMPONENT " << "DOF " << dofname << std::endl;
			TR << "MULTICOMPONENT " << " moves these components:";
			for ( utility::vector1<char>::const_iterator i = compchild.begin(); i != compchild.end(); ++i ) {
				TR << " " << *i;
			}
			TR << std::endl;
			TR << "MULTICOMPONENT " << " moves these subunits:";
			for ( utility::vector1<Size>::const_iterator i = subchild.begin(); i != subchild.end(); ++i ) {
				TR << " " << *i;
			}
			TR << std::endl;
		}
		for ( utility::vector1<char>::const_iterator i = symm_info.get_components().begin(); i != symm_info.get_components().end(); ++i ) {
			TR << "MULTICOMPONENT " << *i << " " << symm_info.get_component_bounds().find(*i)->second.first << "-" << symm_info.get_component_bounds().find(*i)->second.second << std::endl;
			TR << "MULTICOMPONENT " << " moved by dofs:";
			for ( std::map< Size, SymDof >::const_iterator j = symm_info.get_dofs().begin(); j != symm_info.get_dofs().end(); ++j ) {
				std::string const & dofname( symm_info.get_jump_name(j->first) );
				utility::vector1<char> compchild = symmdata.components_moved_by_jump(dofname);
				if ( std::find(compchild.begin(),compchild.end(),*i) == compchild.end() ) continue;
				TR << " " << dofname;
			}
			TR << std::endl << "MULTICOMPONENT " << " contains virtuals:";
			for ( auto const & j : symm_info.get_subunit_name_to_component() ) {
				if ( j.second != *i ) continue;
				TR << " " << j.first;
			}
			TR << std::endl;
		}
	}

	// now create the symmetric conformation
	conformation::symmetry::SymmetricConformationOP symm_conf;

	if ( has_mirror_operations ) {
		symm_conf = conformation::symmetry::SymmetricConformationOP (
			new conformation::symmetry::MirrorSymmetricConformation( conf, symm_info ) );
	} else {
		symm_conf = conformation::symmetry::SymmetricConformationOP (
			new conformation::symmetry::SymmetricConformation( conf, symm_info ) );
	}

	// now we make the coordinates of the pose symmetric
	// -- apply independent jumps so the structure is created symmetrically
	for ( Size i=1; i<= f.num_jump(); ++i ) {
		if ( symm_info.jump_is_independent( i ) ) {
			symm_conf->set_jump( i, conf.jump( i ) );
		}
	}

	// Now, we need to update the residue identities of symmetric conformations (to ensure that ALA becomes DALA in mirrored subunits, for example):
	if ( has_mirror_operations ) {
		conformation::symmetry::MirrorSymmetricConformationOP mirror_conf( utility::pointer::dynamic_pointer_cast< conformation::symmetry::MirrorSymmetricConformation >( symm_conf ) );
		debug_assert( mirror_conf );
		mirror_conf->update_residue_identities();
	}

	// renumber the dof information to take the internal jumps into consideration
	core::conformation::symmetry::shift_jump_numbers_in_dofs( *symm_conf, (njump_orig+1-symmdata.get_num_components()) * symmdata.get_subunits() );

	TR << "=================== SYM FOLD TREE, jump notation: =symfixed= *indep* #symdof# jump[=follows] ========================\n"
		<< show_foldtree(*symm_conf,symmdata,get_chain2range(src_conformation,src_conf2pdb_chain)) << std::endl;

	return symm_conf;
}

// @details Make a fold tree from data stored in the SymmData object. This would have to make sense.
// this function would probably enjoy more checks...
kinematics::FoldTree
set_fold_tree_from_symm_data(
	conformation::Conformation & src_conformation,
	conformation::symmetry::SymmData & symmdata,
	std::map< core::Size, char > src_conf2pdb_chain
)
{
	std::map<char,std::pair<Size,Size> > const chain2range = get_chain2range(src_conformation,src_conf2pdb_chain);

	using namespace kinematics;
	FoldTree f, f_orig = src_conformation.fold_tree();
	f.clear();

	// Save number of monomer residues, subunits and virtuals
	Size num_res_subunit( src_conformation.size() );
	Size subunits( symmdata.get_subunits() );
	Size num_cuts_subunit( f_orig.num_cutpoint() );
	Size num_virtuals( symmdata.get_virtual_coordinates().size() );
	Size num_virtual_connects ( symmdata.get_virtual_connects().size() );
	Size num_res_real( num_res_subunit*subunits );
	Size anchor ( process_residue_request(src_conformation,symmdata.get_anchor_residue() ) );

	// Check that input data makes sense
	if ( anchor > num_res_subunit ) {
		utility_exit_with_message( "Anchor outside the subunit..." );
	}

	if ( symmdata.get_num_components() == 1 ) {


		Size num_jumps_subunit( f_orig.num_jump() );

		//store chemical edge information
		utility::vector1< core::kinematics::Edge > edges_subunit( f_orig.get_chemical_edges() );
		Size num_edges_subunit = edges_subunit.size();

		// Store number of jumps, cuts, and account for chemical edges.
		Size njumps( 0 );
		Size total_num_cuts( num_cuts_subunit*subunits + num_edges_subunit*subunits + subunits + num_virtuals - 1 );
		Size total_num_jumps( num_jumps_subunit*subunits + num_edges_subunit*subunits + num_virtual_connects );

		// store information from subunit foldtree
		utility::vector1< Size > cuts_subunit( f_orig.cutpoints() );
		ObjexxFCL::FArray1D< Size > cuts( total_num_cuts );
		ObjexxFCL::FArray2D< Size > jump_points_subunit( 2, num_jumps_subunit+num_edges_subunit ),
			jumps( 2, total_num_jumps );

		for ( Size i = 1; i<= f_orig.num_jump(); ++i ) {
			jump_points_subunit(1,i) = f_orig.upstream_jump_residue(i);
			jump_points_subunit(2,i) = f_orig.downstream_jump_residue(i);
		}

		// Now add jumps and cuts for the other subunits
		for ( Size i = 0; i < subunits; ++i ) {
			for ( Size j = 1; j <= num_jumps_subunit; ++j ) {
				++njumps;
				//int new_cut_pos( i*num_res_subunit + cuts_subunit.at( j ) );
				//cuts( njumps ) = new_cut_pos;
				int new_jump_pos1( i*num_res_subunit + jump_points_subunit(1,j) );
				int new_jump_pos2( i*num_res_subunit + jump_points_subunit(2,j) );
				//TR << " jump as jump " << new_jump_pos1 << " " << new_jump_pos2 << " " << new_cut_pos << std::endl;
				jumps(1, njumps ) = std::min(new_jump_pos1,new_jump_pos2);
				jumps(2, njumps ) = std::max(new_jump_pos1,new_jump_pos2);
			}
			//Add the edge points as jumps. This is a hack to get around tree_from_jumps_and_cuts ignoring edges.
			// They are replaced with true chemical edges later
			for ( Size j = 1; j <= num_edges_subunit; ++j ) {
				++njumps;
				//int new_cut_pos( i*num_res_subunit + edges_subunit.at( j ).stop()-1 );
				//cuts( njumps ) = new_cut_pos;
				int new_jump_pos1( i*num_res_subunit + edges_subunit.at(j).start() );
				int new_jump_pos2( i*num_res_subunit + edges_subunit.at(j).stop() );
				//TR << " edge as jump " << new_jump_pos1 << " " << new_jump_pos2 << " " << new_cut_pos << std::endl;
				jumps(1, njumps ) = std::min(new_jump_pos1,new_jump_pos2);
				jumps(2, njumps ) = std::max(new_jump_pos1,new_jump_pos2);
			}
			//Add the cuts
			for ( Size j=1; j<=cuts_subunit.size(); j++ ) {
				cuts( i*cuts_subunit.size()+j ) = i*num_res_subunit + cuts_subunit[j];
			}
		}

		// why copy the maps?
		std::map< std::string, std::pair< std::string, std::string > > const virtual_connect( symmdata.get_virtual_connects() );
		std::map< std::string, Size > virtual_id_to_num (symmdata.get_virtual_id_to_num() );
		std::map< Size, std::string > virtual_num_to_id (symmdata.get_virtual_num_to_id() );
		std::map< std::string, Size > virtual_id_to_subunit (symmdata.get_virt_id_to_subunit_num() );

		std::map< std::string, std::pair< std::string, std::string > >::const_iterator it_start = virtual_connect.begin();
		std::map< std::string, std::pair< std::string, std::string > >::const_iterator it_end = virtual_connect.end();
		for ( std::map< std::string, std::pair< std::string, std::string > >::const_iterator it = it_start; it != it_end; ++it ) {
			std::pair< std::string, std::string > connect( it->second );
			std::string pos_id1( connect.first );
			std::string pos_id2( connect.second );
			if ( pos_id2 == "SUBUNIT" ) {
				int subunit = virtual_id_to_subunit.find( pos_id1 )->second;
				++njumps;
				cuts( njumps ) = num_res_subunit*subunit;
				int jump_pos1( ( subunit-1 )*num_res_subunit + anchor );
				int jump_pos2( num_res_real + virtual_id_to_num.find( pos_id1 )->second );
				jumps(1, njumps ) = std::min(jump_pos1,jump_pos2);
				jumps(2, njumps ) = std::max(jump_pos1,jump_pos2);
			}
		}

		// Read jump structure from symmdata

		// Read the virtual connect data and add the new connections. They are stored in virtual_connects. We also need to
		// know the mapping between mapping between name of virtual residues and the jump number
		Size pos( 0 );
		for ( std::map< std::string, std::pair< std::string, std::string > >::const_iterator it = it_start; it != it_end; ++it ) {
			std::pair< std::string, std::string > connect( it->second );
			std::string pos_id1( connect.first );
			std::string pos_id2( connect.second );
			// We have already added the jumps from virtual residues to their corresponding subunits
			if ( pos_id2 == "SUBUNIT" ) continue;
			++njumps;
			++pos;
			debug_assert( virtual_id_to_num.find( pos_id1 ) != virtual_id_to_num.end() && virtual_id_to_num.find( pos_id2 ) != virtual_id_to_num.end() );
			Size pos1 ( virtual_id_to_num.find( pos_id1 )->second );
			Size pos2 ( virtual_id_to_num.find( pos_id2 )->second );
			cuts( njumps ) = num_res_real + pos;
			jumps(1, njumps ) = num_res_real + std::min(pos1,pos2);
			jumps(2, njumps ) = num_res_real + std::max(pos1,pos2);
		}

		// Now create foldtree
		f.tree_from_jumps_and_cuts( num_res_real + num_virtuals, njumps, jumps, cuts );

		//Replace the chemical edges into the fold tree before reordering
		for ( Size i=0; i < subunits; i++ ) {
			for ( Size j=1; j<=edges_subunit.size(); j++ ) {
				core::kinematics::Edge subunit_edge = edges_subunit[j];
				Size start = i*num_res_subunit + edges_subunit[j].start();
				Size stop = i*num_res_subunit + edges_subunit[j].stop();
				int jump_label = f.edge_label(start,stop);
				core::kinematics::Edge current_edge = f.jump_edge(jump_label);
				core::kinematics::Edge future_edge( start, stop, -2, subunit_edge.start_atom(), subunit_edge.stop_atom(), subunit_edge.keep_stub_in_residue());
				f.replace_edge(current_edge,future_edge);
				f.renumber_jumps_ordered();
			}
		}

		f.reorder( num_res_real + 1 );

		// now set jump atoms from f_orig
		{
			std::map<int,std::string> downstream_res_to_jump_atom;
			for ( Size i = 1; i <= f_orig.num_jump(); ++i ) {
				downstream_res_to_jump_atom[f_orig.downstream_jump_residue(i)] = f_orig.downstream_atom(i);
			}
			for ( Size i = 1; i <= f.num_jump(); ++i ) {
				int upres = f.upstream_jump_residue(i);
				int dnres = f.downstream_jump_residue(i);
				int upres1 = (upres-1) % num_res_subunit + 1;
				int dnres1 = (dnres-1) % num_res_subunit + 1;
				if ( dnres > (int)num_res_real ) continue;
				if ( downstream_res_to_jump_atom.find(dnres1) == downstream_res_to_jump_atom.end() ) continue;
				if ( downstream_res_to_jump_atom[dnres1] == "" ) continue;
				std::string upatom = upres > (int)num_res_real ? "X" : src_conformation.residue(upres1).atom_name(1);
				std::string dnatom = downstream_res_to_jump_atom[dnres1];
				TR << "set SYM jump atoms: " << upres << ":" << upatom << " " << dnres << ":" << dnatom << std::endl;
				f.set_jump_atoms(i, upatom, dnatom );
			}
		}

		return f;


	} else if ( symmdata.get_num_components() > 1 ) {
		using std::map;
		using std::endl;
		using std::string;
		using std::pair;
		using std::make_pair;
		using ObjexxFCL::format::I;

		Size num_jumps_subunit( f_orig.num_jump() + 1 - symmdata.get_num_components() );

		// Store number of jumps and cuts
		Size njumps( 0 ), ncuts(0);
		Size total_num_cuts( num_cuts_subunit*subunits + subunits + num_virtuals - 1 );
		Size total_num_jumps( num_jumps_subunit*subunits + num_virtual_connects );
		utility::vector1< Size > cuts_subunit( f_orig.cutpoints() );
		ObjexxFCL::FArray2D< Size > jumps( 2, total_num_jumps );

		// TR << "NUM COMPONENTS: " << symmdata.get_num_components() << std::endl;
		ObjexxFCL::FArray2D_int jump_points_subunit( 2, num_jumps_subunit );
		Size nsubjump = 0;
		for ( Size i = 1; i<= f_orig.num_jump(); ++i ) {
			Size const up = f_orig.upstream_jump_residue(i);
			Size const dn = f_orig.downstream_jump_residue(i);
			bool up_and_dn_within_component = is_jump_intracomponent(chain2range,up,dn);
			if ( !up_and_dn_within_component ) {
				TR << "this intrasub jump will be replaced with a SUBUNIT jump " << up << " " << dn << endl;
				continue;
			}
			nsubjump++;
			jump_points_subunit(1,nsubjump) = up;
			jump_points_subunit(2,nsubjump) = dn;
		}
		if ( nsubjump != num_jumps_subunit ) {
			TR << "f_orig.num_jump()             " << f_orig.num_jump() << std::endl;
			TR << "symmdata.get_num_components() " << symmdata.get_num_components() << std::endl;
			TR << "nsubjump                      " << nsubjump          << std::endl;
			TR << "num_jumps_subunit             " << num_jumps_subunit << std::endl;
			TR << f_orig << std::endl;
			utility_exit_with_message("intra-subunit jump mismatch!");
		}

		// Now add jumps and cuts for the other subunits
		for ( Size i = 0; i < subunits; ++i ) {
			for ( Size j = 1; j <= nsubjump; ++j ) {
				++njumps;
				int new_jump_pos1( i*num_res_subunit + jump_points_subunit(1,j) );
				int new_jump_pos2( i*num_res_subunit + jump_points_subunit(2,j) );
				jumps(1, njumps ) = std::min(new_jump_pos1,new_jump_pos2);
				jumps(2, njumps ) = std::max(new_jump_pos1,new_jump_pos2);
			}
		}


		std::map< std::string, std::pair< std::string, std::string > > const & virtual_connect( symmdata.get_virtual_connects() );
		std::map< std::string, Size > const & virtual_id_to_num (symmdata.get_virtual_id_to_num() );
		// std::map< Size, std::string > const & virtual_num_to_id (symmdata.get_virtual_num_to_id() );
		std::map< std::string, Size > const & virtual_id_to_subunit (symmdata.get_virt_id_to_subunit_num() );
		// for(map<string,Size>::const_iterator it = virtual_id_to_subunit.begin(); it != virtual_id_to_subunit.end(); ++it) {
		//  TR << "virtual_id_to_subunit " << it->first << " " << it->second << endl;
		// }
		map< string, char > const & virtual_id_to_subunit_chain(symmdata.get_virt_id_to_subunit_chain() );
		map< string, string > const & virtual_id_to_subunit_residue(symmdata.get_virt_id_to_subunit_residue() );
		// sanity check for SUBUNIT chain / residue
		// for(map<char,pair<Size,Size> >::const_iterator i = chain2range.begin(); i != chain2range.end(); ++i) {
		//  TR << "get_chain2range " << i->first << " " << i->second.first << "-" << i->second.second << endl;
		// }
		for ( auto const & elem : virtual_id_to_subunit_chain ) {
			string const & virt = elem.first;
			char chain = elem.second;
			if ( chain == (char)0 ) chain = src_conf2pdb_chain.count(1) ? src_conf2pdb_chain[1] : 'A';
			if ( chain2range.find(chain)==chain2range.end() ) {
				std::cerr << "missing chain in pdb: " << chain << std::endl;
				utility_exit_with_message("missing chain in symfile/pdb");
			}
			Size beg = chain2range.find(chain)->second.first;
			Size end = chain2range.find(chain)->second.second;
			Size resi = process_residue_request( src_conformation, virtual_id_to_subunit_residue.find(virt)->second, beg, end );
			if ( !( beg <= resi && resi <= end ) ) {
				TR << "SUBUNIT chain " << chain << " assertion fail: " << beg << " <= " << resi << " <= " << end << endl;
				utility_exit_with_message("requested residue out of bounds: "+virtual_id_to_subunit_residue.find(virt)->second);
			}
			// TR << "SUBUNIT chain: " << virt << " " << chain << " " << resi << endl;
		}

		for ( auto const & elem : virtual_connect ) {
			std::pair< std::string, std::string > connect( elem.second );
			std::string pos_id1( connect.first );
			std::string pos_id2( connect.second );
			if ( pos_id2 == "SUBUNIT" ) {
				char chain = virtual_id_to_subunit_chain.find(pos_id1)->second;
				if ( chain == (char)0 ) chain = src_conf2pdb_chain.count(1) ? src_conf2pdb_chain[1] : 'A';
				Size beg = chain2range.find(chain)->second.first;
				Size end = chain2range.find(chain)->second.second;
				Size resi = process_residue_request( src_conformation, virtual_id_to_subunit_residue.find(pos_id1)->second, beg, end );
				// TR << "set SUBUNIT jumps " << pos_id1 << " chain: " << chain << " anchor: " << resi << " range: " << beg << "-" << end << endl;
				int subunit = virtual_id_to_subunit.find( pos_id1 )->second;
				++njumps;
				Size this_anchor = resi ? resi : anchor;
				if ( this_anchor < beg || end < this_anchor ) utility_exit_with_message("bad anchor "+virtual_id_to_subunit_residue.find(pos_id1)->second);
				// cuts from subunit already
				// cuts(njumps) = (subunit-1)*num_res_subunit + end;
				int jump_pos1( (subunit-1)*num_res_subunit + this_anchor );
				int jump_pos2( num_res_real + virtual_id_to_num.find( pos_id1 )->second );
				jumps(1, njumps ) = std::min(jump_pos1,jump_pos2);
				jumps(2, njumps ) = std::max(jump_pos1,jump_pos2);
			}
		}

		// Read jump structure from symmdata

		// Read the virtual connect data and add the new connections. They are stored in virtual_connects. We also need to
		// know the mapping between mapping between name of virtual residues and the jump number
		Size pos( 0 );
		for ( auto const & elem : virtual_connect ) {
			std::pair< std::string, std::string > connect( elem.second );
			std::string pos_id1( connect.first );
			std::string pos_id2( connect.second );
			// We have already added the jumps from virtual residues to their corresponding subunits
			if ( pos_id2 == "SUBUNIT" ) continue;
			++njumps;
			++pos;
			debug_assert( virtual_id_to_num.find( pos_id1 ) != virtual_id_to_num.end() && virtual_id_to_num.find( pos_id2 ) != virtual_id_to_num.end() );
			Size pos1 ( virtual_id_to_num.find( pos_id1 )->second );
			Size pos2 ( virtual_id_to_num.find( pos_id2 )->second );
			jumps(1, njumps ) = num_res_real + std::min(pos1,pos2);
			jumps(2, njumps ) = num_res_real + std::max(pos1,pos2);
		}

		ObjexxFCL::FArray1D< Size > cuts( total_num_cuts );
		for ( Size i = 0; i < subunits; ++i ) {
			for ( Size j = 1; j <= f_orig.num_jump(); ++j ) {
				cuts( ++ncuts ) = i*num_res_subunit + cuts_subunit.at(j);
			}
			cuts(++ncuts) = num_res_subunit*(i+1);
		}
		for ( Size i = 1; i < num_virtuals; ++i ) {
			cuts(++ncuts) = i+num_res_real;
		}


		// // Now create foldtree
		// TR << "create new fold tree, oldnjump " << f.num_jump() << " newnjump " << njumps << std::endl;
		// for(int i = 1; i <= total_num_cuts; ++i){
		//  TR << i << " " << cuts(i) << std::endl;
		// }
		f.tree_from_jumps_and_cuts( num_res_real + num_virtuals, njumps, jumps, cuts, num_res_real+1 );
		f.reorder( num_res_real + 1 );

		// now set jump atoms from f_orig
		{
			std::map<int,std::string> downstream_res_to_jump_atom;
			for ( Size i = 1; i <= f_orig.num_jump(); ++i ) {
				downstream_res_to_jump_atom[f_orig.downstream_jump_residue(i)] = f_orig.downstream_atom(i);
			}
			for ( Size i = 1; i <= f.num_jump(); ++i ) {
				int upres = f.upstream_jump_residue(i);
				int dnres = f.downstream_jump_residue(i);
				int upres1 = (upres-1) % num_res_subunit + 1;
				int dnres1 = (dnres-1) % num_res_subunit + 1;
				if ( dnres > (int)num_res_real ) continue;
				if ( downstream_res_to_jump_atom.find(dnres1) == downstream_res_to_jump_atom.end() ) continue;
				if ( downstream_res_to_jump_atom[dnres1] == "" ) continue;
				std::string upatom = upres > (int)num_res_real ? "X" : src_conformation.residue(upres1).atom_name(1);
				std::string dnatom = downstream_res_to_jump_atom[dnres1];
				TR << "set SYM jump atoms: " << upres << ":" << upatom << " " << dnres << ":" << dnatom << std::endl;
				f.set_jump_atoms(i, upatom, dnatom );
			}
		}

		return f;

	} else {
		TR << "num_componints: " << symmdata.get_num_components() << std::endl;
		utility_exit_with_message("BAD num_componints");
	}

	return FoldTree(0); // this will never happen
}

// this function is STILL under construction...
kinematics::FoldTree
replaced_symmetric_foldtree_with_new_monomer(
	kinematics::FoldTree symm_f,
	conformation::symmetry::SymmetryInfo symmetry_info,
	kinematics::FoldTree monomer_f
)
{
	using namespace kinematics;
	FoldTree f=symm_f, f_new=monomer_f;

	// Save number of monomer residues, subunits and virtuals
	Size num_res_subunit( symmetry_info.num_independent_residues() );
	Size subunits( symmetry_info.subunits() );
	//Size anchor ( symmdata.get_anchor_residue() );


	// Now add jumps and cuts for the other subunits
	for ( Size i = 0; i < subunits; ++i ) {
		int begin ( i*num_res_subunit+1 );
		int end ( (i+1)*num_res_subunit );
		f.delete_segment( begin, end );
		//f.insert_fold_tree_by_jump(f, insert_seqpos, insert_jumppos, anchor_pos, anchor_jump_number, anchor_atom, root_atom);

	}

	return f;
}


// @details center the pose at the origin. Use the CA of the anchor residue
// as the anchor point
void
recenter(
	conformation::Conformation & src_conformation,
	conformation::symmetry::SymmData & symmdata
)
{
	conformation::Residue anchor ( src_conformation.residue( process_residue_request(src_conformation,symmdata.get_anchor_residue()) ) );
	Vector trans( anchor.xyz( "CA" ) );
	for ( Size i = 1; i <= src_conformation.size(); ++i ) {
		for ( Size j = 1; j <= src_conformation.residue_type(i).natoms(); ++j ) {
			id::AtomID id( j, i );
			src_conformation.set_xyz( id, src_conformation.xyz(id) - trans );
		}
	}
}

// @details shift jump numbers in dof
void
shift_jump_numbers_in_dofs(
	conformation::Conformation & conformation,
	Size shift
)
{
	using namespace core::conformation::symmetry;

	runtime_assert( is_symmetric( conformation ) );
	SymmetricConformation & symm_conf (
		dynamic_cast<SymmetricConformation & > ( conformation ) );
	SymmetryInfoOP symm_info( symm_conf.Symmetry_Info() );

	std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );
	std::map< Size, SymDof > dofs_shifted;
	std::map< Size, SymDof >::iterator it;
	std::map< Size, SymDof >::iterator it_begin = dofs.begin();
	std::map< Size, SymDof >::iterator it_end = dofs.end();

	for ( it = it_begin; it != it_end; ++it ) {
		int jump_nbr ( (*it).first + shift );
		dofs_shifted.insert( std::make_pair( jump_nbr, (*it).second ) );
	}
	symm_info->set_dofs( dofs_shifted );
}


kinematics::FoldTree
get_asymm_unit_fold_tree( core::conformation::Conformation const &conf ) {
	if ( !is_symmetric( conf ) ) {
		return conf.fold_tree();
	}

	// basic idea: only take edges with at least one end in the 1st subunit
	kinematics::FoldTree const &f = conf.fold_tree();
	kinematics::FoldTree f_new;

	conformation::symmetry::SymmetricConformation const & symm_conf ( dynamic_cast<conformation::symmetry::SymmetricConformation const & > ( conf ) );
	conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
	Size nres_subunit ( symm_info->num_independent_residues() );

	for ( auto const & it : f ) {
		if ( it.start() <= nres_subunit && it.stop() <= nres_subunit ) {
			f_new.add_edge( it );
		} else if ( it.stop() <= nres_subunit ) { // this is the jump to the subunit
			core::kinematics::Edge e_new = it;
			e_new.start() = nres_subunit + 1;
			f_new.add_edge( e_new );
		}
	}
	f_new.renumber_jumps_ordered();

	//TR.Error << "get_asymm_unit_fold_tree() called with " << f << std::endl;
	//TR.Error << "get_asymm_unit_fold_tree() returning " << f_new << std::endl;
	return f_new;
}

/// @brief    Converts an asymmetric foldtree (f) with virtual root into a
///           symmetric foldtree compatible with symmetric conformation (conf)
/// @param    conf - A symmetric conformation
/// @param    f - An asymmetric foldtree. This foldtree MUST have a virtual root
/// @details  This function does not require the symm data
void
symmetrize_fold_tree(
	core::conformation::Conformation const &conf,
	kinematics::FoldTree & f
) {

	if ( !is_symmetric( conf ) ) {
		return;
	}

	// basic idea: grabs jumps and cuts from monomer, symmetrize them,
	//    retain VRT jumps and cuts, generate a new fold tree
	kinematics::FoldTree f_orig = f;
	//TR.Error << "symmetrize_fold_tree() called with " << f << std::endl;
	//TR.Error << "reference fold tree =  " << p.fold_tree() << std::endl;

	conformation::symmetry::SymmetricConformation const & symm_conf (
		dynamic_cast<conformation::symmetry::SymmetricConformation const & > ( conf ) );
	conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
	Size nres_subunit ( symm_info->num_independent_residues() );
	Size nsubunits ( symm_info->subunits() );
	Size num_nonvrt = symm_info->num_total_residues_without_pseudo();

	// 1 - Get monomer jumps, cuts, and anchor
	Size new_anchor = 0;
	utility::vector1< Size > new_cuts;
	utility::vector1< std::pair<Size,Size> > new_jumps;

	for ( Size i = 1; i<= f_orig.num_jump(); ++i ) {
		Size down ( f_orig.downstream_jump_residue(i) );
		Size up ( f_orig.upstream_jump_residue(i) );
		if ( up > nres_subunit && new_anchor == 0 ) {
			new_anchor = down;
			TR << "symmetrize_fold_tree(): setting anchor to " << new_anchor << std::endl;
		} else {
			if ( up > nres_subunit && new_anchor != 0 ) {
				TR << "Warning!  Fold tree has multiple jumps from anchor, which is not allowed in symmetry.  Modifying." << std::endl;
				up = new_anchor;
			}

			// add this jump in each subunit
			for ( Size j=0;  j<nsubunits; ++j ) {
				// dont worry about order, they get sorted later
				new_jumps.push_back( std::pair<int,int>( j*nres_subunit+down, j*nres_subunit+up ) );
			}
		}
	}
	utility::vector1< int > cuts_vector( f_orig.cutpoints() );

	// 1A - symmetrize
	for ( Size i = 0;  i<nsubunits; ++i ) {
		for ( Size j = 1; j<=cuts_vector.size(); ++j ) {
			if ( i*nres_subunit + cuts_vector[j] != num_nonvrt ) {
				new_cuts.push_back( i*nres_subunit + cuts_vector[j] );
			}
		}
	}

	// 2 - Get symmetic jumps cuts and anchor
	// inter-VRT jumps
	kinematics::FoldTree const & f_pose = conf.fold_tree();
	for ( Size i = 1; i<= f_pose.num_jump(); ++i ) {
		Size down ( f_pose.downstream_jump_residue(i) );
		Size up ( f_pose.upstream_jump_residue(i) );
		// connections between VRTs are unchanged
		if ( up > num_nonvrt && down > num_nonvrt ) {
			new_jumps.push_back( std::pair<Size,Size>(up,down) );
		}
		// jumps to new anchor
		if ( up > num_nonvrt && down <= num_nonvrt ) {
			int subunit_i = symm_info->subunit_index(down);
			new_jumps.push_back( std::pair<Size,Size>( up, (subunit_i-1)*nres_subunit + new_anchor ) );
		}
		// jumps from new anchor
		if ( up <= num_nonvrt && down > num_nonvrt ) {
			int subunit_i = symm_info->subunit_index(up);
			new_jumps.push_back( std::pair<Size,Size>( (subunit_i-1)*nres_subunit + new_anchor , down ) );
		}
	}

	// cuts
	cuts_vector = f_pose.cutpoints();
	for ( Size i = 1; i<=cuts_vector.size(); ++i ) {
		if ( cuts_vector[i] >= (int) num_nonvrt ) {
			new_cuts.push_back( cuts_vector[i] );
		}
	}

	// 3 - combine
	ObjexxFCL::FArray1D< Size > cuts( new_cuts.size() );
	ObjexxFCL::FArray2D< Size > jumps( 2, new_jumps.size() );
	// Initialize jumps
	for ( Size i = 1; i<= new_jumps.size(); ++i ) {
		jumps(1,i) = std::min( new_jumps[i].first, new_jumps[i].second);
		jumps(2,i) = std::max( new_jumps[i].first, new_jumps[i].second);
		// DEBUG -- PRINT JUMPS AND CUTS
		//TR.Error << " jump " << i << " : " << jumps(1,i) << " , " << jumps(2,i) << std::endl;
	}
	for ( Size i = 1; i<= new_cuts.size(); ++i ) {
		cuts(i) = new_cuts[i];
		//TR.Error << " cut " << i << " : " << cuts(i) << std::endl;
	}

	// 4 make the foldtree
	f.clear();
	f.tree_from_jumps_and_cuts( conf.size(), new_jumps.size(), jumps, cuts, f_pose.root(), false ); // true );
	//TR.Error << "symmetrize_fold_tree() before reorder " << f << std::endl;
	//f.reorder( f_pose.root() );
	//TR.Error << "symmetrize_fold_tree() returning with " << f << std::endl;
}


void
set_asymm_unit_fold_tree( core::conformation::Conformation &conf, kinematics::FoldTree const &f) {
	if ( !is_symmetric( conf ) ) {
		conf.fold_tree( f );
	}

	kinematics::FoldTree f_new = f;
	symmetrize_fold_tree( conf, f_new );

	conf.fold_tree( f_new );
}

// this function is directly stolen from the docking code in protocols. The code duplication is introduced
// to avoid dependencies of core functionality from protocols
int
residue_center_of_mass(
	conformation::Conformation const & conformation,
	int const start,
	int const stop
)
{
	Vector center( 0.0 );
	for ( int i=start; i<=stop; ++i ) {
		//if( !pose.residue( i ).is_protein() ) continue;
		if ( !conformation.residue( i ).is_protein() ) {
			Vector ca_pos( conformation.residue( i ).nbr_atom_xyz() );
			center += ca_pos;
		} else {
			Vector ca_pos( conformation.residue( i ).atom( "CA" ).xyz() );
			center += ca_pos;
		}
	}
	center /= (stop-start+1);

	return core::conformation::symmetry::return_nearest_residue( conformation, start, stop, center );
}

// this function is directly stolen from the docking code in protocols. The code duplication is introduced
// to avoid dependencies of protocols on core functionality
int
return_nearest_residue(
	conformation::Conformation const & conformation,
	int const begin,
	int const end,
	Vector center
)
{
	Real min_dist = 9999.9;
	int res = 0;
	for ( int i=begin; i<=end; ++i ) {
		Vector ca_pos;
		if ( !conformation.residue( i ).is_protein() ) {
			ca_pos = conformation.residue( i ).nbr_atom_xyz();
		} else {
			ca_pos = conformation.residue( i ).atom( "CA" ).xyz() ;
		}

		ca_pos -= center;
		Real tmp_dist( ca_pos.length_squared() );
		if ( tmp_dist < min_dist ) {
			res = i;
			min_dist = tmp_dist;
		}
	}
	return res;
}

std::string
show_foldtree(
	core::conformation::symmetry::SymmetricConformation const & symm_conf,
	SymmData const & symmdata,
	std::map<char,std::pair<Size,Size> > const & chain2range
){
	Size nsub = symm_conf.Symmetry_Info()->subunits();
	Size nres_monomer = symm_conf.Symmetry_Info()->num_independent_residues();
	Size nreal = nsub * nres_monomer;

	// label real res by chain
	std::map<Size,std::string> labels;
	for ( auto const & i : chain2range ) {
		for ( Size is = 1; is <= nsub; ++is ) {
			for ( Size ir = i.second.first; ir <= i.second.second; ++ir ) {
				labels[(is-1)*nres_monomer+ir] = std::string("Sub") + ObjexxFCL::string_of(is) + std::string("") + i.first;
			}
		}
	}
	// label virtuals by name in symdef
	Size Nreal = symm_conf.Symmetry_Info()->num_total_residues_without_pseudo();
	for ( auto const & i : symmdata.get_virtual_num_to_id() ) {
		labels[i.first+Nreal] = i.second;
	}
	std::map<Size,char> mark_jump_to_res;
	for ( Size i = 1; i <= symm_conf.fold_tree().num_jump(); ++i ) {
		if ( symm_conf.Symmetry_Info()->jump_is_independent(i) ) {
			mark_jump_to_res[symm_conf.fold_tree().downstream_jump_residue(i)] = '*';
		}
	}
	for ( Size i = 1; i <= symm_conf.fold_tree().num_jump(); ++i ) {
		Size up = symm_conf.fold_tree().upstream_jump_residue(i);
		Size dn = symm_conf.fold_tree().downstream_jump_residue(i);
		if ( symm_conf.Symmetry_Info()->jump_is_independent(i) && up > nreal && dn > nreal ) {
			mark_jump_to_res[symm_conf.fold_tree().downstream_jump_residue(i)] = '=';
		}
	}
	std::map<Size,SymDof> dofs( symm_conf.Symmetry_Info()->get_dofs() );
	for ( std::map<Size,SymDof>::const_iterator i = dofs.begin(); i != dofs.end(); ++i ) {
		// TR << "SymDOF " << i->first << " " << i->second << std::endl;
		mark_jump_to_res[symm_conf.fold_tree().downstream_jump_residue(i->first)] = '#';
	}
	std::map<Size,Size> jump_follows;
	for ( Size i = 1; i <= symm_conf.fold_tree().num_jump(); ++i ) {
		if ( symm_conf.Symmetry_Info()->jump_follows(i)!=(Size)0 ) {
			jump_follows[i] = symm_conf.Symmetry_Info()->jump_follows(i);
		}
	}

	// for(Size i = 1; i <= symm_conf.fold_tree().num_jump(); ++i){
	//  using ObjexxFCL::format::I;
	//  Size up = symm_conf.fold_tree().upstream_jump_residue(i);
	//  Size dn = symm_conf.fold_tree().downstream_jump_residue(i);
	//  TR << "JUMP " << I(3,i) << " " << I(4,up) << " -> " << I(4,dn)
	//     << " " << (symm_conf.Symmetry_Info()->jump_is_independent(i) ? "I" : " ")
	//     << " " << I(3,symm_conf.Symmetry_Info()->jump_follows(i))
	//     << " " << labels[up] << " " << labels[dn]
	//     << " " << (dofs.find(i)!=dofs.end() ? "DOF" : "" )
	//     << std::endl;
	// }

	return core::kinematics::visualize_fold_tree(symm_conf.fold_tree(),labels,mark_jump_to_res,jump_follows);

}


} // symmetry
} // pose
} // core
