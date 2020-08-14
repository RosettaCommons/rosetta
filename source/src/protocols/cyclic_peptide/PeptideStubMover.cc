// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/PeptideStubMover.cc
/// @brief The PeptideStubMover prepends, appends, or inserts residues into an existing pose, or builds a new polymeric chain
/// @author Yifan Song
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitute.org): Added support for stripping
/// N-acetylation and C-methylamidation when appending residues, preserving phi and the previous
/// omega in the first case and psi and the following omega in the second.
/// @modified Jack Maguire, jackmaguire1444@gmail.com: breaking apply() into smaller methods

#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/cyclic_peptide/PeptideStubMoverCreator.hh>

#include <core/id/AtomID.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/variant_util.hh>
#include <core/id/AtomID.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/iterators.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>

#include <numeric/conversions.hh>

#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <algorithm>

static basic::Tracer TR( "protocols.cyclic_peptide.PeptideStubMover" );

using namespace core::select;

namespace protocols {
namespace cyclic_peptide {

PeptideStubMover::PeptideStubMover(){
	init();
}

PeptideStubMover::PeptideStubMover( PeptideStubMover const & ) = default;

PeptideStubMover::~PeptideStubMover()= default;

void PeptideStubMover::init() {
	reset_ = false;
	update_pdb_numbering_ = true;
}

///@brief if TR.Debug is visible, dump residue information to the console
///@author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
///@modified Jack Maguire, jackmaguire1444@gmail.com : moved to its own function
void dump_debug_output( core::pose::Pose const & pose ){
	//pose.dump_pdb("test.pdb");
	if ( TR.Debug.visible() ) {
		TR << pose.annotated_sequence() << std::endl;
		pose.fold_tree().show(TR.Debug);

		for ( core::Size ires : resids( pose ) ) {
			if ( pose.fold_tree().is_jump_point(ires) ) {
				TR.Debug << "jump: Residue " << ires << std::endl;
			}
			if ( pose.fold_tree().is_cutpoint(ires) ) {
				TR.Debug << "cut: Residue " << ires << std::endl;
			}

			for ( core::Size icon : utility::enumerate1(pose.residue_type(ires).n_possible_residue_connections()) ) {
				TR.Debug << "connection: Residue " << ires << " Atom " << pose.residue_type(ires).residue_connection(icon).atomno() << " to residue " << pose.residue(ires).connected_residue_at_resconn(icon) << ", connect id:" << pose.residue(ires).connect_map(icon).connid() << std::endl;
			}
		}
	}
}


void PeptideStubMover::apply( core::pose::Pose & pose )
{
	core::chemical::ResidueTypeSetCOP standard_residues = pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t );
	if ( reset_ ) {
		pose.clear();
	}

	using namespace core::chemical;

	utility::vector1< core::Size > resids_that_we_added;

	for ( core::Size istub : utility::indices1(stub_rsd_names_) ) {
		FirstResidAdded const first_resid =
			perform_single_iteration( pose, * standard_residues, istub );
		core::Size const n_residues_added =
			std::max( core::Size(1), core::Size(stub_rsd_repeat_[ istub ]) );

		for ( core::Size & r : resids_that_we_added ) {
			//Move downstream residues even further downstream.
			if ( r >= first_resid ) r += n_residues_added;
		}
		for ( core::Size ii = first_resid; ii < first_resid + n_residues_added; ++ii ) {
			resids_that_we_added.push_back( ii );
		}

		if ( TR.Debug.visible() ) TR.Debug << "Done appending " << stub_rsd_names_[istub] << std::endl;
	}// for istub

	dump_debug_output( pose );

	//If residues have been added, they need PDB chain IDs and numbers to be updated.
	if ( update_pdb_numbering_ ) update_pdb_numbering( pose, resids_that_we_added );
}


/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
PeptideStubMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & data
) {
	using namespace core;
	using namespace core::select::residue_selector;

	if ( tag->hasOption( "reset" ) ) set_reset_mode( tag->getOption< bool >( "reset" ) );
	if ( tag->hasOption( "update_pdb_numbering" ) ) set_update_pdb_numbering_mode(tag->getOption< bool >( "update_pdb_numbering" ) );

	if ( tag->hasOption( "label" ) ) {
		set_residue_label( tag->getOption< std::string >( "label" ) );
	}

	reset_mover_data();

	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	for ( tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it ) {
		utility::tag::Tag const & tag = **tag_it;

		core::Size anchor_rsd = 0;
		ResidueSelectorCOP anchor_rsd_selector = nullptr;

		if ( tag.hasOption( "anchor_rsd" ) ) {
			anchor_rsd = tag.getOption< core::Size >( "anchor_rsd" );
		} else if ( tag.hasOption( "anchor_rsd_selector" ) ) {
			anchor_rsd_selector = parse_residue_selector( *tag_it, data, "anchor_rsd_selector" );
		}

		add_residue(
			(*tag_it)->getName(),
			(*tag_it)->getOption<std::string>( "resname", "" ),
			(*tag_it)->getOption<core::Size>( "position", 0 ),
			(*tag_it)->getOption<bool>( "jump", false ),
			(*tag_it)->getOption<std::string>( "connecting_atom", "" ),
			(*tag_it)->getOption<core::Size>( "repeat", 1 ),
			anchor_rsd,
			anchor_rsd_selector,
			(*tag_it)->getOption<std::string>( "anchor_atom", "" )
		);
		/*if ( (*tag_it)->getName() == "Append" ) {
		stub_mode_.push_back(append);
		}
		else if ( (*tag_it)->getName() == "Prepend" ) {
		stub_mode_.push_back(prepend);
		}
		else if ( (*tag_it)->getName() == "Insert" ) {
		stub_mode_.push_back(insert);
		}
		stub_rsd_names_.push_back( (*tag_it)->getOption<std::string>( "resname", "" ) );
		stub_insert_pos_.push_back( (*tag_it)->getOption<core::Size>( "position", 0 ) );
		stub_rsd_jumping_.push_back( (*tag_it)->getOption<bool>( "jump", false ) );
		stub_rsd_connecting_atom_.push_back( (*tag_it)->getOption<std::string>( "connecting_atom", "" ) );
		stub_rsd_repeat_.push_back( (*tag_it)->getOption<core::Size>( "repeat", 1 ) );
		stub_anchor_rsd_.push_back( (*tag_it)->getOption<core::Size>( "anchor_rsd", 0 ) );
		stub_anchor_rsd_connecting_atom_.push_back( (*tag_it)->getOption<std::string>( "anchor_atom", "" ) ); */
	}
}


/// @brief Adds a residue to the list of residues to be appended, prepended, or inserted.
/// @details Calls add_residue() override that uses PSM_StubMode.
void PeptideStubMover::add_residue(
	std::string const &stubmode,
	std::string const &resname,
	core::Size const position,
	bool const jumpmode,
	std::string const &connecting_atom,
	core::Size const repeat,
	core::Size const anchor_rsd,
	core::select::residue_selector::ResidueSelectorCOP anchor_rsd_selector,
	std::string const &anchor_atom
) {

	runtime_assert_string_msg( stubmode=="Append" || stubmode=="Prepend" || stubmode=="Insert", "In PeptideStubMover::add_residue(): unrecognized stub creation mode.  Mode must be one of \"Append\", \"Prepend\", or \"Insert\"." );

	PSM_StubMode stubmode2;

	if ( stubmode == "Append" ) {
		stubmode2 = PSM_append;
	} else if ( stubmode == "Prepend" ) {
		stubmode2 = PSM_prepend;
	} else {
		debug_assert( stubmode == "Insert" ); //Should be true, given assertion above.
		stubmode2 = PSM_insert;
	}

	add_residue( stubmode2, resname, position, jumpmode, connecting_atom, repeat, anchor_rsd, anchor_rsd_selector, anchor_atom );
}

/// @brief Adds a residue to the list of residues to be appended, prepended, or inserted.
/// @details This version uses enums.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
PeptideStubMover::add_residue(
	PSM_StubMode const stubmode,
	std::string const &resname,
	core::Size const position,
	bool const jumpmode,
	std::string const &connecting_atom,
	core::Size const repeat,
	core::Size const anchor_rsd,
	core::select::residue_selector::ResidueSelectorCOP anchor_rsd_selector,
	std::string const &anchor_atom
) {
	runtime_assert_string_msg( stubmode == PSM_append || stubmode == PSM_prepend || stubmode == PSM_insert, "Error in PeptideStubMover::add_residue(): Invalid stub mode provided." );

	stub_mode_.push_back(stubmode);
	stub_rsd_names_.push_back( resname );
	stub_insert_pos_.push_back( position );
	stub_rsd_jumping_.push_back( jumpmode );
	stub_rsd_connecting_atom_.push_back( connecting_atom );
	stub_rsd_repeat_.push_back( repeat );
	stub_anchor_rsd_.push_back( anchor_rsd );
	anchor_rsd_selectors_.push_back( anchor_rsd_selector );
	stub_anchor_rsd_connecting_atom_.push_back( anchor_atom );
}


moves::MoverOP PeptideStubMover::clone() const { return utility::pointer::make_shared< PeptideStubMover >( *this ); }
moves::MoverOP PeptideStubMover::fresh_instance() const { return utility::pointer::make_shared< PeptideStubMover >(); }





//Private functions:

/// @brief Rebuilds all atoms that are dependent on bonds between residue_index and any other residues (including atoms on the other residues).
void PeptideStubMover::rebuild_atoms(
	core::pose::Pose &pose,
	core::Size const residue_index
) const {
	debug_assert(residue_index <= pose.size() && residue_index > 0);

	core::Size const nresconn = pose.residue(residue_index).n_possible_residue_connections();
	if ( nresconn == 0 ) return;

	for ( core::Size const ic : utility::enumerate1( nresconn ) ) {
		if ( !pose.residue(residue_index).connection_incomplete(ic) ) {
			core::Size const conn_at_index = pose.residue(residue_index).residue_connect_atom_index(ic); //The index of the connection atom
			for ( core::Size const ia : utility::enumerate1(pose.residue(residue_index).natoms() ) ) {
				if ( pose.residue(residue_index).icoor(ia).stub_atom(1).atomno()==conn_at_index ) { //If this is a child of the connection atom
					//TR << "Rebuilding rsd " << residue_index << " atom " << ia << " (" << pose.residue(residue_index).atom_name(ia) << ")" << std::endl; TR.flush();
					pose.conformation().set_xyz( core::id::AtomID( ia, residue_index ), pose.residue(residue_index).icoor(ia).build(pose.residue(residue_index), pose.conformation()) );
				}
			}
			core::Size const other_rsd_index = pose.residue(residue_index).connect_map(ic).resid(); //The residue index of the other residue that this one is connected to at connection ID "ic"
			core::Size const other_rsd_connect_index = pose.residue(residue_index).connect_map(ic).connid(); //The connection index on the other residue that links to THIS residue.
			core::Size const other_at_index = pose.residue(other_rsd_index).residue_connect_atom_index( other_rsd_connect_index ); //The index of the atom on the other residue that links to THIS residue.
			for ( core::Size const ia : utility::enumerate1(pose.residue(other_rsd_index).natoms() ) ) {
				if ( pose.residue(other_rsd_index).icoor(ia).stub_atom(1).atomno()==other_at_index ) { //If this is an immediate child of the other residue's atom that connects to this atom, rebuild it.
					pose.conformation().set_xyz( core::id::AtomID(ia, other_rsd_index), pose.residue(other_rsd_index).icoor(ia).build( pose.residue(other_rsd_index), pose.conformation() ) );
				}
			}//ia
		}//if ! connection_incomplete(ic)
	}//ic
}

/// @brief Updates the PDB numbering (PDB number/chain ID) as residues are added.
void
PeptideStubMover::update_pdb_numbering (
	core::pose::Pose &pose,
	utility::vector1< core::Size > & resids_that_we_added
) const {

	//Create pdbinfo if it does not exist
	if ( pose.pdb_info() == nullptr ) {
		pose.pdb_info( utility::pointer::make_shared< core::pose::PDBInfo >( pose, true ) );
	}
	core::pose::PDBInfoOP pdbinfo = pose.pdb_info();
	debug_assert( pdbinfo != nullptr );

	//Add labels
	bool const add_labels = !residue_label_.empty();
	if ( add_labels ) {
		for ( core::Size const resid : resids_that_we_added ) {
			runtime_assert( resid > 0 );
			runtime_assert( resid <= pose.size() );
			pdbinfo->add_reslabel( resid, residue_label_ );
		}
	}

	assign_chain_ids( pose );
}

std::string PeptideStubMover::get_name() const {
	return mover_name();
}

std::string PeptideStubMover::mover_name() {
	return "PeptideStubMover";
}

void PeptideStubMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "reset", xsct_rosetta_bool, "Delete the input pose" )
		+ XMLSchemaAttribute( "update_pdb_numbering", xsct_rosetta_bool, "Update the PDB chain IDs and numbers for this pose" )
		+ XMLSchemaAttribute::attribute_w_default( "label", xs_string, "Residue label for residues added by this mover. This requires update_pdb_numbering=\"true\"", "PEPTIDE_STUB_EXTENSION" );
	//go through subelements
	AttributeList subelement_attributes; //all have same attributes
	subelement_attributes
		+ XMLSchemaAttribute( "resname", xs_string, "Name of stub residue" )
		+ XMLSchemaAttribute( "position", xsct_non_negative_integer, "Position to insert stub residue" )
		+ XMLSchemaAttribute::attribute_w_default( "jump", xsct_rosetta_bool, "Append residue by jump?", "false" )
		+ XMLSchemaAttribute( "connecting_atom", xs_string, "Name of atom where residues should be connected" )
		+ XMLSchemaAttribute::attribute_w_default( "repeat", xsct_non_negative_integer, "Number of times to add this residue", "1" )
		+ XMLSchemaAttribute( "anchor_rsd", xsct_non_negative_integer, "Residue to which the stub residue should bond" )
		+ XMLSchemaAttribute( "anchor_rsd_selector", xs_string, "Selector alternative to \"anchor_rsd\"" )
		+ XMLSchemaAttribute( "anchor_atom", xs_string, "Atom to which the stub residue should bond" );
	XMLSchemaSimpleSubelementList subelements;
	//Append, Prepend, or Insert
	subelements
		.add_simple_subelement( "Append", subelement_attributes, "Specify a residue to append at a position" )
		.add_simple_subelement( "Prepend", subelement_attributes, "Specify a residue to prepend at a position" )
		.add_simple_subelement( "Insert", subelement_attributes, "Specify a residue to insert at a position" );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "Used to append/insert/prepend residues to a pose", attlist, subelements );
}

void
PeptideStubMover::assign_chain_ids(
	core::pose::Pose & pose
){
	if ( pose.pdb_info() == nullptr ) {
		pose.pdb_info( utility::pointer::make_shared< core::pose::PDBInfo >( pose, true ) );
	}
	runtime_assert( pose.pdb_info() != nullptr );
	core::pose::PDBInfo & pdbinfo = * pose.pdb_info();

	std::string const chainids = core::chemical::one_indexed_chr_chains();

	for ( core::Size const ir : resids( pose ) ) {
		core::Size const nlabels_before = pdbinfo.get_reslabels( ir ).size();

		char const chainid = chainids[ std::min< core::Size >( pose.chain(ir), 62 ) ];

		pdbinfo.set_resinfo(
			ir,
			chainid,
			static_cast<int>(ir),
			' ',//default
			"    ",//default
			false
		);

		core::Size const nlabels_after = pdbinfo.get_reslabels( ir ).size();
		runtime_assert( nlabels_before == nlabels_after );
	}
}

///@brief returns resid
FirstResidAdded
PeptideStubMover::perform_single_iteration(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & standard_residues,
	core::Size const istub
){

	core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( standard_residues.name_map(stub_rsd_names_[istub]) ) );

	// first stub always starts with a jump
	if ( istub == 1 && pose.size() == 0 ) {
		return add_residue_to_empty_pose( pose, *new_rsd, istub );
	}

	core::Size const anchor_rsd = get_anchor_rsd( pose, istub );

	if ( stub_rsd_jumping_[istub] ) { //Appending by jump.
		return append_by_jump( pose, anchor_rsd, *new_rsd, istub );
	}

	return append_by_bond( pose, anchor_rsd, *new_rsd, istub );
}

///@brief this function performs the early logic of append_by_bond()
///@details call handle_upper_terminus and handle_lower_terminus
void
PeptideStubMover::handle_termini_and_store_terminal_dihedrals(
	core::pose::Pose & pose,
	core::Size const anchor_rsd,
	core::Size const istub,
	core::Size const connecting_id,
	core::Size & anchor_connecting_id,
	core::Real & old_omega_minus1,
	core::Real & old_phi,
	core::Real & old_psi,
	core::Real & old_omega,
	bool & replace_upper_terminal_type,
	bool & replace_lower_terminal_type
){
	if ( stub_anchor_rsd_connecting_atom_[istub].empty() ) {
		if ( stub_mode_[istub] == PSM_append ) {
			handle_upper_terminus( pose, anchor_rsd, old_psi, old_omega, replace_upper_terminal_type );
			anchor_connecting_id = pose.residue_type(anchor_rsd).upper_connect_id();
		} else if ( stub_mode_[istub] == PSM_prepend ) {
			handle_lower_terminus( pose, anchor_rsd, old_phi, old_omega_minus1, replace_lower_terminal_type );
			anchor_connecting_id = pose.residue_type(anchor_rsd).lower_connect_id();
		}
	} else {
		core::Size const atomid(pose.residue_type(anchor_rsd).atom_index(stub_anchor_rsd_connecting_atom_[istub]));
		anchor_connecting_id = pose.residue_type(anchor_rsd).residue_connection_id_for_atom(atomid);

		if ( connecting_id == 0 ) {
			TR << "Error! Residue " << stub_rsd_names_[anchor_rsd] << " Atom " << stub_anchor_rsd_connecting_atom_[istub] << " cannot be connected" << std::endl;
		}
	}
	if ( (!replace_upper_terminal_type) && anchor_connecting_id == pose.residue_type(anchor_rsd).upper_connect_id() ) {
		handle_upper_terminus( pose, anchor_rsd, old_psi, old_omega, replace_upper_terminal_type );
	}
	if ( (!replace_lower_terminal_type) && anchor_connecting_id == pose.residue_type(anchor_rsd).lower_connect_id() ) {
		handle_lower_terminus( pose, anchor_rsd, old_phi, old_omega_minus1, replace_lower_terminal_type );
	}
}

///@brief determine how residues should connect when being appended by bond
core::Size
PeptideStubMover::get_connecting_id_for_append_by_bond(
	core::conformation::Residue const & new_rsd,
	core::Size const istub
){
	if ( stub_rsd_connecting_atom_[istub] == "" ) {
		core::Size const connecting_id = new_rsd.type().lower_connect_id();
		runtime_assert( connecting_id != 0 );
		return connecting_id;
	} else {
		core::Size atomid(new_rsd.atom_index(stub_rsd_connecting_atom_[istub]));
		core::Size const connecting_id = new_rsd.type().residue_connection_id_for_atom(atomid);
		if ( connecting_id == 0 ) {
			TR << "Error! Residue " << stub_rsd_names_[istub] << " Atom " << stub_rsd_connecting_atom_[istub] << " cannot be connected" << std::endl;
		}
		return connecting_id;
	}
}

///@brief Handle the case where the new residues start a new chain
///@eturns returns resid of the first new residue added
FirstResidAdded
PeptideStubMover::append_by_jump(
	core::pose::Pose & pose,
	core::Size const anchor_rsd,
	core::conformation::Residue & new_rsd,
	core::Size const istub
){
	runtime_assert_string_msg(stub_mode_[istub] == PSM_append, "Can only use append for jumps");
	pose.append_residue_by_jump(new_rsd, anchor_rsd);
	core::Size const first_resid_added = pose.size();

	for ( core::Size i_repeat = 2; i_repeat <= stub_rsd_repeat_[istub]; ++i_repeat ) {
		pose.append_residue_by_bond(new_rsd, true);
	}

	return first_resid_added;
}

///@brief Handle the case where the input pose has zero residues
///@eturns returns resid of the first new residue
FirstResidAdded
PeptideStubMover::add_residue_to_empty_pose(
	core::pose::Pose & pose,
	core::conformation::Residue & new_rsd,
	core::Size const istub
){
	runtime_assert_string_msg(stub_mode_[istub] == PSM_append, "Can only use append for the first residue");
	pose.append_residue_by_jump(new_rsd, 1);
	for ( core::Size i_repeat = 2; i_repeat <= stub_rsd_repeat_[istub]; ++i_repeat ) {
		pose.append_residue_by_bond(new_rsd, true);
	}
	return 1;
}

///@brief queries stub_anchor_rsd_ and anchor_rsd_selectors_ to find the residue to use as an achor
///@returns resid of anchor residue residue
core::Size
PeptideStubMover::get_anchor_rsd(
	core::pose::Pose const & pose,
	core::Size const istub
){
	if ( stub_anchor_rsd_[istub] != 0 ) return stub_anchor_rsd_[istub];

	//Residue Selector
	if ( anchor_rsd_selectors_[ istub ] != nullptr ) {
		core::Size selected_resid = 0;
		for ( core::Size resid : resids( pose, {anchor_rsd_selectors_[ istub ]} ) ) {
			if ( selected_resid != 0 ) utility_exit_with_message( "Multiple residues were selected to be the anchor residue" );

			selected_resid = resid;
			//don't break, keep looking to make sure there aren't multiple selected
		}
		return selected_resid;
	}

	if ( stub_mode_[istub] == PSM_prepend ) {
		return 1;
	} else {
		return pose.size();
	}
}

///@brief returns resid of new residue
FirstResidAdded
PeptideStubMover::append_by_bond(
	core::pose::Pose & pose,
	core::Size const anchor_rsd,
	core::conformation::Residue & new_rsd,
	core::Size const istub
){
	core::Size const connecting_id = get_connecting_id_for_append_by_bond( new_rsd, istub );

	core::Size anchor_connecting_id(0);
	core::Size new_index( anchor_rsd ); //The anchor residue's new index after appending or prepending or whatnot.
	core::Real old_omega_minus1(0), old_phi(0), old_psi(0), old_omega(0); //Only used to preserve phi and omega if we're replacing a C-methylamidation or an N-acetylation.
	bool replace_upper_terminal_type(false), replace_lower_terminal_type(false);
	handle_termini_and_store_terminal_dihedrals(
		pose,
		anchor_rsd,
		istub,
		connecting_id,
		anchor_connecting_id,
		old_omega_minus1,
		old_phi,
		old_psi,
		old_omega,
		replace_upper_terminal_type,
		replace_lower_terminal_type
	);

	FirstResidAdded resid_of_residue_being_added = 0;//This will be updated

	if ( stub_mode_[istub] == PSM_append ) {
		pose.append_residue_by_bond(new_rsd, true, connecting_id, anchor_rsd, anchor_connecting_id, false);
		new_index = anchor_rsd;
		resid_of_residue_being_added = new_index + 1;
		//rebuild the polymer bond dependent atoms:
		rebuild_atoms(pose, new_index);
	} else if ( stub_mode_[istub] == PSM_prepend ) {
		pose.prepend_polymer_residue_before_seqpos(new_rsd, anchor_rsd, true);
		new_index = anchor_rsd + 1; //The anchor residue's index has shifted by one.
		resid_of_residue_being_added = new_index - 1;
		rebuild_atoms(pose, new_index);
	} else if ( stub_mode_[istub] == PSM_insert ) {
		if ( stub_insert_pos_[istub] != 0 ) {
			pose.insert_residue_by_bond(new_rsd, stub_insert_pos_[istub], anchor_rsd, true, stub_anchor_rsd_connecting_atom_[istub], stub_rsd_connecting_atom_[istub], false, true);
			if ( stub_insert_pos_[istub] <= anchor_rsd ) {
				new_index = anchor_rsd + 1;
				resid_of_residue_being_added = new_index - 1;
			} else {
				new_index = anchor_rsd;
				resid_of_residue_being_added = new_index + 1;
			}
			//pose.dump_pdb("test.pdb");
		} else {
			pose.append_polymer_residue_after_seqpos(new_rsd, anchor_rsd, true);
			new_index = anchor_rsd;
			resid_of_residue_being_added = new_index + 1;
		}
		rebuild_atoms(pose, new_index);
	}

	//Update the omega-1 and phi (if we've replaced an N-acetylation) or the psi and omega (if we've replaced a C-methylamidation) to preserve these dihedral values:
	if ( !reset_ ) {
		preserve_old_mainchain_torsions( pose, new_index, old_omega_minus1, old_phi, old_psi, old_omega, replace_upper_terminal_type, replace_lower_terminal_type );
	}

	handle_repeats_in_append_by_bond( pose, anchor_rsd, new_rsd, istub );

	runtime_assert( resid_of_residue_being_added != 0 );
	return resid_of_residue_being_added;

}

///@brief add additional residues if needed
void
PeptideStubMover::handle_repeats_in_append_by_bond(
	core::pose::Pose & pose,
	core::Size const anchor_rsd,
	core::conformation::Residue & new_rsd,
	core::Size const istub
){
	for ( core::Size i_repeat = 2; i_repeat <= stub_rsd_repeat_[istub]; ++i_repeat ) {
		if ( stub_rsd_connecting_atom_[istub].empty() && stub_anchor_rsd_connecting_atom_[istub].empty() ) {
			if ( stub_mode_[istub] == PSM_append ) {
				pose.append_residue_by_bond(new_rsd, true);
			} else if ( stub_mode_[istub] == PSM_insert ) {
				if ( stub_insert_pos_[istub] == 0 ) {
					pose.append_polymer_residue_after_seqpos(new_rsd, anchor_rsd + i_repeat - 1, true);
				}
			} else if ( stub_mode_[istub] == PSM_prepend ) {
				pose.prepend_polymer_residue_before_seqpos(new_rsd, anchor_rsd, true);
			}
		} else {
			TR << "Cannot use repeat for non-canonical insertion" << std::endl;
			break;
		}
	}//for i_repeat
}

/// @brief Remove terminal types from the upper terminus.  Store the old psi and omega values.
/// @returns  Returns void, but if a terminal type was removed, old_psi will be set to the previous
/// psi value and old_omega will be set to the previous omega value.  The replace_upper_terminal_type var
/// is set to true if a terminal type was replaced and false otherwise.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
PeptideStubMover::handle_upper_terminus(
	core::pose::Pose & pose,
	core::Size const anchor_rsd,
	core::Real & old_psi,
	core::Real & old_omega,
	bool & replace_upper_terminal_type
) const {
	replace_upper_terminal_type = false;
	bool const is_valid_type( pose.residue_type(anchor_rsd).is_alpha_aa() || pose.residue_type(anchor_rsd).is_beta_aa() || pose.residue_type(anchor_rsd).is_peptoid() );
	if ( pose.residue_type( anchor_rsd ).is_upper_terminus() ||
			pose.residue_type( anchor_rsd ).has_variant_type( core::chemical::C_METHYLAMIDATION ) ||
			pose.residue_type( anchor_rsd ).has_variant_type(core::chemical::CUTPOINT_LOWER )
			) {
		if ( is_valid_type ) {
			replace_upper_terminal_type = true;
			old_psi=pose.psi(anchor_rsd);
			//old omega stored below:
		}
		if ( pose.residue_type( anchor_rsd ).has_variant_type( core::chemical::C_METHYLAMIDATION ) ) {
			if ( is_valid_type ) old_omega=pose.omega(anchor_rsd);
			core::pose::remove_variant_type_from_pose_residue(pose, core::chemical::C_METHYLAMIDATION, anchor_rsd);
		} else if ( pose.residue_type(anchor_rsd).has_variant_type(core::chemical::CUTPOINT_LOWER) ) {
			if ( is_valid_type ) old_omega=0.0;
			core::pose::remove_variant_type_from_pose_residue(pose, core::chemical::CUTPOINT_LOWER, anchor_rsd);
		} else {
			if ( is_valid_type ) old_omega=0.0;
			core::pose::remove_upper_terminus_type_from_pose_residue(pose, anchor_rsd);
		}
	}
}

/// @brief Remove terminal types from the lower terminus.  Store the old phi and omega_nminus1 values.
/// @returns  Returns void, but if a terminal type was removed, old_phi will be set to the previous
/// phi value and old_omega_nminus1 will be set to the previous upstream omega value.  The
/// replace_lower_terminal_type var is set to true if a terminal type was replaced and false otherwise.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
PeptideStubMover::handle_lower_terminus(
	core::pose::Pose & pose,
	core::Size const anchor_rsd,
	core::Real & old_phi,
	core::Real & old_omega_nminus1,
	bool & replace_lower_terminal_type
) const {
	replace_lower_terminal_type = false;
	bool const is_valid_type( pose.residue_type(anchor_rsd).is_alpha_aa() || pose.residue_type(anchor_rsd).is_beta_aa() || pose.residue_type(anchor_rsd).is_peptoid() );
	if ( pose.residue_type( anchor_rsd ).is_upper_terminus() ||
			pose.residue_type( anchor_rsd ).has_variant_type( core::chemical::N_ACETYLATION ) ||
			pose.residue_type( anchor_rsd ).has_variant_type(core::chemical::CUTPOINT_UPPER )
			) {
		if ( is_valid_type ) {
			replace_lower_terminal_type = true;
			old_phi=pose.phi(anchor_rsd);
			//old omega stored below:
		}
		if ( pose.residue_type( anchor_rsd ).has_variant_type( core::chemical::N_ACETYLATION ) ) {
			if ( is_valid_type ) {
				core::chemical::ResidueType const & restype( pose.residue_type(anchor_rsd) );
				old_omega_nminus1 = numeric::conversions::degrees( pose.conformation().torsion_angle(
					core::id::AtomID( restype.atom_index("CQ"), anchor_rsd ),
					core::id::AtomID( restype.atom_index("CP"), anchor_rsd ),
					core::id::AtomID( restype.atom_index("N"), anchor_rsd ),
					core::id::AtomID( restype.atom_index("CA"), anchor_rsd )
					) );
			}
			core::pose::remove_variant_type_from_pose_residue(pose, core::chemical::N_ACETYLATION, anchor_rsd);
		} else if ( pose.residue_type(anchor_rsd).has_variant_type(core::chemical::CUTPOINT_UPPER) ) {
			if ( is_valid_type ) old_omega_nminus1=0.0;
			core::pose::remove_variant_type_from_pose_residue(pose, core::chemical::CUTPOINT_UPPER, anchor_rsd);
		} else {
			if ( is_valid_type ) old_omega_nminus1=0.0;
			core::pose::remove_upper_terminus_type_from_pose_residue(pose, anchor_rsd);
		}
	}
}

/// @brief Update the omega-1 and phi (if we've replaced an N-acetylation) or the psi and omega
/// (if we've replaced a C-methylamidation) to preserve these dihedral values.
/// @details Builds a temporary foldtree rooted on the alpha carbon of the anchor residue.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
PeptideStubMover::preserve_old_mainchain_torsions(
	core::pose::Pose & pose,
	core::Size const anchor_res,
	core::Real const old_omega_minus1,
	core::Real const old_phi,
	core::Real const old_psi,
	core::Real const old_omega,
	bool const replace_upper_terminal_type,
	bool const replace_lower_terminal_type
) const {
	if ( !(replace_upper_terminal_type || replace_lower_terminal_type ) ) return; //Do nothing if there was no terminal type replaced.

	core::chemical::ResidueType const & restype( pose.residue_type(anchor_res) );
	if ( !(restype.is_alpha_aa() || restype.is_beta_aa() || restype.is_peptoid() ) ) return; //Do nothing if this isn't a suitable type.


	// Store the old foldtree:
	core::kinematics::FoldTree const saved_ft( pose.fold_tree() );

	// Temporarily append a virtual residue, which will be the root of the new foldtree.
	core::conformation::ResidueOP new_virt( core::conformation::ResidueFactory::create_residue( pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t )->name_map("VRT") ) );
	pose.append_residue_by_jump( *new_virt, anchor_res, "", "", true );
	core::Size const virtres( pose.total_residue() ); //Index of the new virtual.

	// Build a temporary foldtree:
	core::kinematics::FoldTree ft_copy;
	ft_copy.clear();
	core::Size jumpcount(1);
	core::Size const chaincount( pose.conformation().num_chains() );

	for ( core::Size chain(1); chain<chaincount; ++chain ) { //Loop through all chains
		core::Size const chnbgn( pose.conformation().chain_begin(chain) );
		core::Size const chnend( pose.conformation().chain_end(chain) );
		if ( anchor_res <= chnend && anchor_res >= chnbgn ) { //If the new root is in the current chain
			ft_copy.add_edge( virtres, anchor_res, jumpcount++, pose.residue_type(virtres).atom_name( pose.residue_type(virtres).root_atom()), "CA");
			if ( anchor_res < chnend ) ft_copy.add_edge( anchor_res, chnend, -1, "CA", pose.residue_type(chnend).atom_name(pose.residue_type(chnend).root_atom()) );
			if ( anchor_res > chnbgn ) ft_copy.add_edge( anchor_res, chnbgn, -1, "CA", pose.residue_type(chnbgn).atom_name(pose.residue_type(chnbgn).root_atom()) );
		} else { //Otherwise, just add the chain.
			ft_copy.add_edge( virtres, chnbgn, jumpcount++);
			ft_copy.add_edge( chnbgn, chnend, -1);
		}
	}
	ft_copy.reorder( virtres );
	if ( TR.visible() ) {
		TR << "Temporary fold tree: " << ft_copy << std::endl;
	}
	pose.fold_tree(ft_copy);

	if ( replace_upper_terminal_type ) {
		pose.set_psi( anchor_res, old_psi );
		pose.set_omega( anchor_res, old_omega );
	}
	if ( replace_lower_terminal_type ) {
		pose.set_phi( anchor_res, old_phi );
		//Ugh.  Need to march upstream to get atoms in prepended residue:
		core::conformation::Residue const & thisres( pose.residue(anchor_res) );
		core::Size const other_res_index( thisres.connected_residue_at_lower() ); //Index of residue connected at this residue's lower connection.
		runtime_assert( other_res_index != 0 ); //Should be true.
		core::conformation::Residue const & otherres( pose.residue(other_res_index) );
		core::Size const other_res_conn_atom( otherres.connect_atom(thisres) ); //The atom on the other residue that connects to this residue.
		runtime_assert( other_res_conn_atom != 0 ); //Should  be true.
		core::chemical::ResidueType const & otherrestype( otherres.type() );
		core::Size const other_res_conn_atom_parent( other_res_conn_atom == otherrestype.root_atom() ? 2 : otherrestype.icoor( other_res_conn_atom ).stub_atom1().atomno() );
		core::Size const this_res_downstream_atom( restype.lower_connect_atom() == restype.root_atom() /*I think this is always true, but better safe than sorry*/ ? 2 : restype.icoor(restype.lower_connect_atom()).stub_atom1().atomno() );

		pose.conformation().set_torsion_angle(
			core::id::AtomID( other_res_conn_atom_parent, other_res_index ),
			core::id::AtomID( other_res_conn_atom, other_res_index ),
			core::id::AtomID( restype.lower_connect_atom(), anchor_res ),
			core::id::AtomID( this_res_downstream_atom, anchor_res ),
			numeric::conversions::radians( old_omega_minus1 )
		);
	}

	// Set a foldtree without a JEDGE for removing the virtual (crashes otherwise):
	ft_copy.clear();
	jumpcount = 1;
	for ( core::Size chain : utility::enumerate1(pose.conformation().num_chains()) ) {
		if ( chain > 1 ) {
			ft_copy.add_edge( 1, pose.conformation().chain_begin(chain), jumpcount++ );
		}
		ft_copy.add_edge( pose.conformation().chain_begin(chain), pose.conformation().chain_end(chain), -1 );
	}
	pose.fold_tree(ft_copy);

	// Remove the virtual
	pose.delete_residue_slow( virtres );

	// Swap the old foldtree back in:
	pose.fold_tree( saved_ft );

	// Update residue neighbours:
	pose.update_residue_neighbors();

}

std::string PeptideStubMoverCreator::keyname() const {
	return PeptideStubMover::mover_name();
}

protocols::moves::MoverOP
PeptideStubMoverCreator::create_mover() const {
	return utility::pointer::make_shared< PeptideStubMover >();
}

void PeptideStubMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PeptideStubMover::provide_xml_schema( xsd );
}


} // moves
} // protocols
