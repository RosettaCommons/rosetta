// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/abscript/RigidChunkCMCM.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/RigidChunkCM.hh>
#include <protocols/abinitio/abscript/RigidChunkCMCreator.hh>

// Package headers
#include <core/environment/DofPassport.hh>

#include <protocols/environment/DofUnlock.hh>
#include <protocols/environment/EnvExcn.hh>
#include <protocols/environment/ProtectedConformation.hh>

#include <protocols/environment/claims/CutBiasClaim.hh>
#include <protocols/environment/claims/JumpClaim.hh>
#include <protocols/environment/claims/XYZClaim.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/AtomTree.hh>

// Project headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.tmpl.hh>
#include <protocols/loops/LoopsFileIO.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

#include <core/chemical/ResidueProperties.hh>

#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

//Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

#include <numeric/random/random.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/xyz.functions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers
#include <algorithm>

// ObjexxFCL Headers

//Req'd on WN32
#include <basic/datacache/WriteableCacheableMap.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer tr( "protocols.abinitio.abscript.RigidChunkCM", basic::t_info );

namespace protocols {
namespace abinitio {
namespace abscript {

using namespace core::environment;
using namespace protocols::environment;

// creator
// XRW TEMP std::string
// XRW TEMP RigidChunkCMCreator::keyname() const {
// XRW TEMP  return RigidChunkCM::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP RigidChunkCMCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new RigidChunkCM );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP RigidChunkCM::mover_name() {
// XRW TEMP  return "RigidChunkCM";
// XRW TEMP }

// This could be adapted to underpin the loops::Loops object, probably.
class ResidueChunkSelection {
	typedef utility::vector1< std::pair< core::Size, core::Size > > ChunkList;

public:
	typedef ChunkList::iterator iterator;
	typedef ChunkList::const_iterator const_iterator;

public:
	ResidueChunkSelection( utility::vector1< bool > sele ) {

		core::Size sele_start = 0;
		for ( core::Size i = 1; i <= sele.size(); ++i ) {
			if ( sele_start ) {
				if ( !sele[i] ) {
					add_chunk( sele_start, i-1 );
					sele_start = 0;
				} else {
					continue;
				}
			} else {
				if ( sele[i] ) {
					sele_start = i;
				} else {
					continue;
				}
			}
		}

		if ( sele_start ) {
			add_chunk( sele_start, sele.size() );
		}
	}

	void add_chunk( core::Size start, core::Size stop ){
		list_.push_back( std::make_pair( start, stop ) );
	}

	iterator begin() { return list_.begin(); }
	const_iterator begin() const { return list_.begin(); }

	iterator end() { return list_.end(); }
	const_iterator end() const { return list_.end(); }

private:
	ChunkList list_;
};


RigidChunkCM::RigidChunkCM():
	Parent(),
	sim_selector_( new core::select::residue_selector::TrueResidueSelector() ),
	templ_selector_( new core::select::residue_selector::TrueResidueSelector() ),
	xml_name_("")
{}

RigidChunkCM::RigidChunkCM(
	core::select::residue_selector::ResidueSelectorCOP selector,
	core::pose::Pose const& template_pose
):
	Parent(),
	template_( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose(template_pose) ) ) ),
	sim_selector_( selector ),
	templ_selector_( selector ),
	xml_name_("")
{}

loops::Loops read_rigid_core( std::string const& file){

	loops::PoseNumberedLoopFileReader reader;
	reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );

	std::ifstream infile( file.c_str() );

	if ( !infile.good() ) {
		tr.Error << "Error opening RBSeg file '" << file << "'" << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( "[ERROR] Error opening RBSeg file '" + file + "'" );
	}

	return loops::Loops( reader.read_pose_numbered_loops_file( infile, file, false ) );
}

void RigidChunkCM::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const&,
	protocols::moves::Movers_map const& movermap,
	core::pose::Pose const& ){

	using namespace basic::options;
	using namespace core::select::residue_selector;

	xml_name_ = std::string( tag->getOption< std::string >( "name" ) );

	std::string const SELECTOR = "selector";
	if ( tag->hasOption( SELECTOR ) ) {
		sim_selector( datamap.get_ptr< ResidueSelector const >( "ResidueSelector", tag->getOption<std::string>( SELECTOR ) ) );
	}

	std::string const APPLY_TO_TEMPLATE = "apply_to_template";
	utility::vector1< moves::MoverOP > apply_movers;
	if ( tag->hasOption( APPLY_TO_TEMPLATE ) ) {
		std::string const apply_to_templates = tag->getOption< std::string >( APPLY_TO_TEMPLATE );
		utility::vector1< std::string > movernames = utility::string_split( apply_to_templates, ',' );
		for ( utility::vector1< std::string >::const_iterator movername = movernames.begin();
				movername != movernames.end(); ++movername ) {
			if ( movermap.find( *movername ) != movermap.end() ) {
				apply_movers.push_back( movermap.find( *movername )->second );
			} else {
				std::ostringstream ss;
				ss << "In mover '" << xml_name_ << "', the 'apply_to_template' tag contained the mover '"
					<< *movername << "', which could not be found." << std::endl;
				throw utility::excn::EXCN_RosettaScriptsOption( ss.str() );
			}
		}
	}

	if ( tag->hasOption("template") ) {
		std::string file = tag->getOption< std::string >( "template" );
		if ( file == "INPUT" ) {
			if ( tag->hasOption( APPLY_TO_TEMPLATE ) ) {
				std::ostringstream ss;
				ss << "In " << this->get_name() << " the option '" << APPLY_TO_TEMPLATE
					<< "' is not combinable with input templates." << std::endl;
				throw utility::excn::EXCN_RosettaScriptsOption( ss.str() );
			}
			template_ = NULL;
		} else {
			core::pose::PoseOP p( new core::pose::Pose() );
			core::import_pose::pose_from_file( *p,
				*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ),
				file );

			for ( Size i = 1; i <= apply_movers.size(); ++i ) {
				if ( apply_movers[i] ) {
					apply_movers[i]->apply( *p );
					tr.Debug << "RigidChunkCM named " << this->get_name() << " applied " << apply_movers[i]->get_name() << std::endl;
				} else {
					std::ostringstream ss;
					ss << "RigidChunkCM named '" << this->get_name() << " couldn't apply one of its input movers because it doesn't exist.";
					throw utility::excn::EXCN_RosettaScriptsOption( ss.str() );
				}
			}
			template_ = p;
		}
	} else {
		throw utility::excn::EXCN_BadInput( "RigidChunkCM requires a template pdb with coordinates for rigid regions.");
	}

	// Determine which region of the template file we want to take.
	if ( tag->hasOption("region_file") ) {
		loops::Loops region_in = read_rigid_core( tag->getOption< std::string >( "region_file" ) );

		// Parse the loops object out into a ResidueIndexSelector. A better way to do this would be
		// to build a constructor for ResidueIndexSelector that takes a loops object, but the Loops
		// object is in protocols.3 and ResidueIndexSelector is in core.4. Lamesauce.
		std::stringstream ss;
		for ( loops::Loops::const_iterator loop = region_in.begin(); loop != region_in.end(); ++loop ) {
			ss << loop->start() << "-" << loop->stop() << ",";
		}
		ResidueSelectorCOP sele( new ResidueIndexSelector( ss.str() ) );
		templ_selector( sele );
	} else if ( tag->hasOption( "region_selector" ) ) {
		ResidueSelectorCOP sele( datamap.get_ptr< ResidueSelector const >( "ResidueSelector", tag->getOption<std::string>( "region_selector" ) ) );
		templ_selector( sele );
	} else if ( tag->hasOption( "region" ) ) {
		ResidueSelectorCOP sele( new ResidueIndexSelector( tag->getOption< std::string >( "region" ) ) );
		templ_selector( sele );
	} else {
		// I don't like this automatic defaulting to anything, but unfortunately some legacy applications
		// require this behavior. Lamesauce.

		tr.Warning << "RigidChunkCM named " << tag->getOption< std::string >( "name", "UNK" )
			<< " defaulted to ALL residues." << std::endl;
	}
}

core::Size find_disulfide_partner( core::pose::Pose const& pose,
	core::Size const resid ) {
	if ( pose.residue( resid ).type().has_variant_type( core::chemical::DISULFIDE ) ) {
		return 0;
	}

	using core::Size;

	if ( !pose.residue( resid ).has( "SG" ) ) {
		// no SG atom => no disulfide (e.g. because in centroid)
		return 0;
	}

	Size const cys_bound_atom( pose.residue( resid ).atom_index( "SG" ) );
	for ( Size jj = pose.residue( resid ).type().n_possible_residue_connections(); jj >= 1; --jj ) {
		if ( (Size) pose.residue( resid ).type().residue_connection( jj ).atomno() == cys_bound_atom ) {
			return pose.residue( resid ).connect_map( jj ).resid();
		}
	}

	return 0;
}

void RigidChunkCM::configure(
	core::pose::Pose const& in_p,
	utility::vector1< bool > const sim_selection
) {

	if ( !template_ ) {
		tr.Debug << "Building template from broker-time pose." << std::endl;
		template_ = core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( in_p ) ) );
	}

	utility::vector1< bool > templ_selection = templ_selector()->apply( templ() );

	{ // Debug output and error checking
		if ( tr.Debug.visible() ) {
			tr.Debug << "template selection for " << this->get_name() << ": ";
			for ( Size i = 1; i <= templ().size(); ++i ) {
				tr.Debug << ( templ_selection[i] ? "T" : "F" );
			}
			tr.Debug << std::endl;
			tr.Debug << "simulation selection for " << this->get_name() << ": ";
			for ( Size i = 1; i <= in_p.size(); ++i ) {
				tr.Debug << ( sim_selection[i] ? "T" : "F" );
			}
			tr.Debug << std::endl;

		}
		if ( templ_selection.index( true ) == 0 ) {
			std::ostringstream ss;
			ss << this->get_name() << " reports that its input loops file (aka rigid core file) contained no residues."
				<< " Check your input." << std::endl;
			throw utility::excn::EXCN_BadInput( ss.str() );
		} else if ( sim_selection.index( true ) == 0 ) {
			std::ostringstream ss;
			ss << this->get_name() << " reports that its selector (used on the simulation to determine where to insert residues)"
				<< " returned an empty selection. Check your input." << std::endl;
			throw utility::excn::EXCN_BadInput( ss.str() );
		}
	}

	core::Size templ_pos = 1;
	core::Size sim_pos = 1;
	std::map< core::Size, core::Size > templ_target;
	std::map< core::Size, core::Size > sim_origin;
	while ( templ_pos <= templ().size() &&
			sim_pos <= in_p.size() ) {
		if ( sim_selection[ sim_pos ] && templ_selection[ templ_pos ] ) {
			// this position is selected in both, it's a correspondence.
			templ_target[ templ_pos ] = sim_pos;
			sim_origin[ sim_pos ] = templ_pos;
			tr.Trace << "Added sim " << sim_pos << "/templ" << templ_pos << " to mapping." << std::endl;
			templ_pos += 1;
			sim_pos += 1;
		} else if ( sim_selection[ sim_pos ] ) {
			// this position is only selected in the simulation. Advance the template one.
			templ_pos += 1;
		} else if ( templ_selection[ templ_pos ] ) {
			// this position is only selected in the template. Advance the simulation one.
			sim_pos += 1;
		} else { // neither
			templ_pos += 1;
			sim_pos += 1;
		}
	}

	this->templ_target( templ_target );
	this->sim_origin( sim_origin );

	{ // Error checking and warnings.
		for ( core::Size i = templ_pos; i <= templ_selection.size(); ++i ) {
			if ( templ_selection[i] ) {
				tr.Warning << this->get_name() << " reports that "
					<< std::count( templ_selection.begin() + i - 1, templ_selection.end(), true )
					<< " residues beginning at " << templ().residue( i ).name3() << i
					<< "in the template do not fit in the selection given. "
					<< "This may or may not be problematic." << std::endl;
				break;
			}
		}
		for ( core::Size i = sim_pos; i <= sim_selection.size(); ++i ) {
			if ( sim_selection[i] ) {
				tr.Warning << this->get_name() << " reports that "
					<< std::count( sim_selection.begin() + i - 1, sim_selection.end(), true )
					<< " residues beginning at " << in_p.residue( i ).name3() << i
					<< " in the template do not fit in the selection given. "
					<< "This may or may not be problematic." << std::endl;
				break;
			}
		}

		for ( core::Size i = 1; i <= templ().size(); ++i ) {
			if ( templ_selection[i] &&
					templ().residue( i ).aa() == core::chemical::aa_cys ) {
				Size cys_partner = find_disulfide_partner( templ(), i );
				if ( cys_partner &&
						templ_target.find( cys_partner ) == templ_target.end() ) {
					std::ostringstream ss;
					ss << "RigidChunkClaimer " << this->get_name() << " reports failure on its selection/template "
						<< " combination because template residue " << templ().residue(cys_partner).name3()
						<< cys_partner << " is disulfide bonded with " << templ().residue( i ).name3() << i
						<< ". RigidChunkClaimer must include both or neither." << std::endl;
					throw utility::excn::EXCN_BadInput( ss.str() );
				}
			}
		}

		for ( core::Size i = 1; i <= in_p.size(); ++i ) {
			if ( sim_selection[i] &&
					in_p.residue( i ).aa() == core::chemical::aa_cys ) {
				Size cys_partner = find_disulfide_partner( in_p, i );
				if ( cys_partner &&
						sim_origin.find( cys_partner ) == sim_origin.end() ) {
					std::ostringstream ss;
					ss << "RigidChunkClaimer " << this->get_name() << " reports failure on its selection target "
						<< " selection because input pose residue " << in_p.residue(cys_partner).name3()
						<< cys_partner << " is disulfide bonded with " << in_p.residue( i ).name3() << i
						<< ". RigidChunkClaimer's rigid chunk must include both or neither." << std::endl;
					throw utility::excn::EXCN_BadInput( ss.str() );
				}
			}
		}
	}
}

claims::EnvClaims RigidChunkCM::yield_claims( core::pose::Pose const& in_p,
	basic::datacache::WriteableCacheableMapOP ){
	using namespace claims;
	EnvClaims claims;

	utility::vector1< bool > selection( in_p.size(), false );
	try {
		selection = sim_selector()->apply( in_p );
	} catch( utility::excn::EXCN_Msg_Exception& e ){
		std::ostringstream ss;
		ss << this->get_name() << " failed to apply its ResidueSelector in " << __FUNCTION__ << ".  ";
		ss << "Pose was length " << in_p.size() << ": " << in_p.sequence();
		e.add_msg(ss.str());
		throw;
	}
	configure( in_p, selection );

	ClientMoverOP this_ptr = utility::pointer::static_pointer_cast< ClientMover > ( get_self_ptr() );

	ResidueChunkSelection simulation_regions( selection );

	std::pair< core::Size, core::Size > prev_region = std::make_pair( 0, 0 );
	for ( ResidueChunkSelection::const_iterator loop_it = simulation_regions.begin();
			loop_it != simulation_regions.end(); ++loop_it ) {
		std::pair< core::Size, core::Size > const excl_region = *loop_it;

		// TODO: remove use of this XYZClaim constructor.
		XYZClaimOP xyz_claim( new XYZClaim( this_ptr, "BASE", excl_region ) );

		xyz_claim->strength( EXCLUSIVE, EXCLUSIVE );
		xyz_claim->set_relative( true ); //we don't care where it's located in space; just that the relative locations are ok.
		claims.push_back( xyz_claim );
		tr.Debug << this->get_name() << ": built EXCLUSIVE XYZClaim for " << excl_region.first << "-" << excl_region.second
			<< " in " << "BASE" << std::endl;

		// For initialization, we need to some of FDP's fixing, which requires some control outside the region
		std::pair< core::Size, core::Size > const supp_region = std::make_pair( std::max( Size( 1 ), excl_region.first-1 ),
			std::min( in_p.size(), excl_region.second+1 ) );
		XYZClaimOP support_claim( new XYZClaim( this_ptr, "BASE", supp_region ) );
		// These would be MUST_CONTROL, but sometimes the edges of the rigid chunk are cuts, and in this
		// case, we want to allow a directly abutting EXCLUSIVE claim. If it's a problem, we'll fail later.
		support_claim->strength( CAN_CONTROL, CAN_CONTROL );
		claims.push_back( support_claim );
		tr.Debug << this->get_name() << ": built support XYZClaim for " << supp_region.first << "-" << supp_region.second
			<< " in " << "BASE" << std::endl;

		claims.push_back( protocols::environment::claims::EnvClaimOP( new CutBiasClaim( this_ptr, "BASE", excl_region, 0.0 ) ) );

		if ( prev_region.first != 0 && prev_region.second != 0 ) {
			core::Size const jump_start = prev_region.first + ( ( prev_region.second - prev_region.first ) / 2 );
			core::Size const jump_end   = excl_region.first + ( ( excl_region.second - excl_region.first ) / 2 );

			JumpClaimOP j_claim( new JumpClaim( this_ptr,
				this->get_name()+"Jump"+utility::to_string( claims.size()/4 ),
				LocalPosition( "BASE", jump_start ),
				LocalPosition( "BASE", jump_end ) ) );
			j_claim->strength( EXCLUSIVE, EXCLUSIVE );
			j_claim->physical( false );

			tr.Debug << this->get_name() << ": built jump claim " << jump_start << "->" << jump_end << std::endl;

			claims.push_back( j_claim );
		}
		prev_region = excl_region;
	}

	return claims;
}
//
//bool missing_density( core::pose::Pose const& pose, loops::Loops const& core, core::Size region_offset ){
//  bool missing_density = false;
//
//  //sanity check: no missing density in backbon in any of the rigid_core residues?
//  for ( loops::Loops::const_iterator it = core.begin(); it!=core.end(); ++it ) {
//    for ( core::Size pos = it->start(); pos <=it->stop(); ++pos ) {
//      // Do we really have Sidechains ?
//      // check this my making sure that no SC atom is more than 20A (?) away from CA
//      numeric::xyzVector< core::Real> ca_pos = pose.residue( pos ).atom("CA").xyz();
//      numeric::xyzVector< core::Real> n_pos = pose.residue( pos ).atom("N").xyz();
//      numeric::xyzVector< core::Real> o_pos = pose.residue( pos ).atom("O").xyz();
//      if ( ( n_pos - ca_pos).length() > 20 || ( ( n_pos - o_pos ).length() > 20 ) ) {
//        tr.Error << "missing backbone in rigid-chunk at " << pos << std::endl;
//        missing_density = true;
//      }
//    }
//  }
//
//  return missing_density;
//}

/////////////// MAGIC FPD CODE FROM RIGID CHUNK CLAIMER /////////////////////////////
void fix_internal_coords_of_siblings( core::pose::Pose& pose,
	core::pose::Pose const& ref_pose,
	core::id::AtomID atom,
	core::id::AtomID ref_atom ) {
	using namespace core::id;

	runtime_assert( atom.rsd() >= 1 && atom.rsd() <= pose.size() );
	runtime_assert( pose.conformation().atom_tree().has( atom ) );
	runtime_assert( ref_pose.conformation().atom_tree().has( ref_atom ) );

	bool has_par1( pose.conformation().atom_tree().atom( atom ).parent() );
	bool ref_has_par1( ref_pose.conformation().atom_tree().atom( ref_atom ).parent() );

	//folding direction matters for the angle we have to set...hence find the parent atoms and get the angle
	AtomID par1O;
	AtomID ref_par1O;
	core::Size const offset = atom.rsd() - ref_atom.rsd();
	if ( has_par1 && ref_has_par1 ) {
		par1O=pose.conformation().atom_tree().atom( atom ).parent()->id();

		std::string const & aname( pose.residue( par1O.rsd() ).atom_name( par1O.atomno() ));
		core::Size const rsd_num = par1O.rsd() - offset;

		ref_par1O=core::id::AtomID( ref_pose.residue( rsd_num ).atom_index( aname ), rsd_num );
	} else {
		tr.Warning << "cannot fix internal coords of " << atom << " in RigidChunk because 1st parent is missing " << std::endl;
		return;
	}
	bool has_par2( pose.conformation().atom_tree().atom( par1O ).parent() );
	bool ref_has_par2( ref_pose.conformation().atom_tree().atom( ref_par1O ).parent() );
	core::id::AtomID par2O;
	core::id::AtomID ref_par2O;
	if ( has_par2 && ref_has_par2 ) {
		par2O=pose.conformation().atom_tree().atom( par1O ).parent()->id();

		std::string const & aname( pose.residue(par2O.rsd()).atom_name( par2O.atomno() ) );
		core::Size const rsd_num = par2O.rsd() - offset;

		ref_par2O=core::id::AtomID( ref_pose.residue( rsd_num ).atom_index( aname ), rsd_num );
	} else {
		tr.Warning << "cannot fix internal coords of " << atom << " in RigidChunk because 2nd parent is missing " << std::endl;
		return;
	}
	runtime_assert( ref_pose.conformation().atom_tree().has( ref_par1O ) );
	runtime_assert( ref_pose.conformation().atom_tree().has( ref_par2O ) );
	runtime_assert( pose.conformation().atom_tree().has( par1O ) );
	runtime_assert( pose.conformation().atom_tree().has( par2O ) );

	core::Real angle( numeric::angle_radians( ref_pose.xyz( ref_atom ), ref_pose.xyz( ref_par1O ), ref_pose.xyz( ref_par2O ) ) );
	tr.Trace << "ref angle direct: " << angle << std::endl;
	pose.conformation().set_bond_angle(  par2O, par1O, atom, angle );

	DOF_ID torsion_offset_dof( atom, PHI );
	DOF_ID ref_torsion_offset_dof( ref_atom, PHI );
	core::Real value( ref_pose.conformation().atom_tree().dof( ref_torsion_offset_dof ) );
	pose.conformation().set_dof( torsion_offset_dof, value );
}

class AtomPack {
public:
	AtomPack( core::Size upper_residue,
		core::conformation::Residue const& prev_rsd,
		core::conformation::Residue const& rsd ) :
		bbM1 ( prev_rsd.mainchain_atom( prev_rsd.n_mainchain_atoms()-2 ),  upper_residue-1 ),
		bb0  ( prev_rsd.mainchain_atom( prev_rsd.n_mainchain_atoms()-1 ),  upper_residue-1 ),
		bb1  ( prev_rsd.mainchain_atom( prev_rsd.n_mainchain_atoms()   ),  upper_residue-1 ),
		bb2  (      rsd.mainchain_atom( 1 ),  upper_residue   ),
		bb3  (      rsd.mainchain_atom( 2 ),  upper_residue   ),
		bb4  (      rsd.mainchain_atom( 3 ),  upper_residue   )
	{}
	core::id::AtomID bbM1;
	core::id::AtomID bb0;
	core::id::AtomID bb1;
	core::id::AtomID bb2;
	core::id::AtomID bb3;
	core::id::AtomID bb4;
};

void fix_mainchain_connect( core::pose::Pose& pose,
	core::Size global_upper,
	core::pose::Pose const& ref_pose,
	core::Size local_upper ) {

	//These guys just hold the various AtomIDs required to correclty address the angles that need to be repaired.
	AtomPack const local = AtomPack( local_upper,
		ref_pose.residue( local_upper-1),
		ref_pose.residue( local_upper ) );
	AtomPack const global = AtomPack( global_upper,
		pose.residue( global_upper-1),
		pose.residue( global_upper ) );

	core::conformation::Residue const & ref_resi = ref_pose.residue( local_upper );
	tr.Trace << "mainchain torsion: ref: " << ref_resi.mainchain_torsion( 1 ) << " atom-tree: "
		<< ref_pose.conformation().torsion_angle( local.bb1, local.bb2, local.bb3, local.bb4 ) << std::endl;

	core::conformation::Residue const & resi = pose.residue( global_upper );
	tr.Trace << "mainchain torsion (before): conf: " << resi.mainchain_torsion( 1 ) << " atom-tree: "
		<< pose.conformation().torsion_angle( global.bb1, global.bb2, global.bb3, global.bb4 ) << std::endl;

	pose.conformation().set_bond_length( global.bb1, global.bb2,
		ref_pose.conformation().bond_length( local.bb1, local.bb2 ) );
	pose.conformation().set_bond_angle ( global.bb0, global.bb1, global.bb2,
		ref_pose.conformation().bond_angle( local.bb0, local.bb1, local.bb2 ) );
	pose.conformation().set_bond_angle ( global.bb1, global.bb2, global.bb3,
		ref_pose.conformation().bond_angle( local.bb1, local.bb2, local.bb3 ) );
	pose.conformation().set_torsion_angle( global.bbM1, global.bb0, global.bb1, global.bb2,
		ref_pose.conformation().torsion_angle( local.bbM1, local.bb0, local.bb1, local.bb2 ) );
	pose.conformation().set_torsion_angle( global.bb0, global.bb1, global.bb2, global.bb3,
		ref_pose.conformation().torsion_angle( local.bb0, local.bb1, local.bb2, local.bb3 ) );
	pose.conformation().set_torsion_angle( global.bb1, global.bb2, global.bb3, global.bb4,
		ref_pose.conformation().torsion_angle( local.bb1, local.bb2, local.bb3, local.bb4 ) );

	core::conformation::Residue const & new_resi = pose.residue( global_upper ); //this should trigger update of coords and torsions
	tr.Trace << "mainchain torsion (after): conf: " << new_resi.mainchain_torsion( 1 ) << " atom-tree: "
		<< pose.conformation().torsion_angle( global.bb1, global.bb2, global.bb3, global.bb4 ) << std::endl;

	core::conformation::Residue const & prev_rsd( ref_pose.residue( local_upper-1 ) );
	core::conformation::Residue const &      rsd( ref_pose.residue( local_upper ) );

	if ( prev_rsd.has( "O" ) ) {
		core::id::AtomID ref_atomO( prev_rsd.atom_index( "O" ), local_upper-1 );
		core::id::AtomID atomO( pose.residue_type( global_upper-1 ).atom_index( "O" ), global_upper-1 );

		fix_internal_coords_of_siblings( pose, ref_pose, atomO, ref_atomO );
	}
	if ( rsd.has( "H" ) ) {
		core::id::AtomID ref_atomH( rsd.atom_index( "H" ), local_upper );
		core::id::AtomID atomH( new_resi.atom_index( "H" ), global_upper );
		runtime_assert( new_resi.has( "H" ) );

		fix_internal_coords_of_siblings( pose, ref_pose, atomH, ref_atomH );
	}

	if ( tr.Trace.visible() ) {
		bool ideal1( core::pose::is_ideal_position( local_upper, ref_pose ) );
		if ( ideal1 && !core::pose::is_ideal_position( global_upper, pose ) ) {
			tr.Warning << " pose in RigidChunkClaimer is not ideal at position " << global_upper << " although template pose was ideal there " << std::endl;
		}

		bool ideal2( core::pose::is_ideal_position( local_upper-1, ref_pose ) );
		if ( ideal2 && !core::pose::is_ideal_position( global_upper-1, pose ) ) {
			tr.Warning << " pose in RigidChunkClaimer is not ideal at position " << global_upper-1 << " although template pose was ideal there " << std::endl;
		}
	}
}
//////////////////////////////// END MAGIC FPD CODE ////////////////////////////////////////////

void RigidChunkCM::initialize( Pose& pose ){

	DofUnlock activation( pose.conformation(), passport() );

	//  if ( missing_density( pose, rigid_core(), region_offset_ ) ) {
	//    throw utility::excn::EXCN_BadInput( " missing density in backbone of rigid-chunk. Check your LOOP definitions.");
	//  }

	core::pose::Pose reference( pose );

	for ( Size sim_pos = 1; sim_pos <= pose.size(); ++sim_pos ) {

		if ( sim_origin().find( sim_pos ) != sim_origin().end() ) {
			Size const templ_pos = sim_origin().at( sim_pos );
			core::conformation::Residue const templ_res = templ().residue( templ_pos );

			try {
				tr.Trace << "Replacing simulation " << pose.residue( sim_pos ).name3() << sim_pos
					<< " with template " << templ_res.name3() << templ_pos << std::endl;

				ProtectedConformation const& conf = static_cast< ProtectedConformation const& >( pose.conformation() );
				pose.replace_residue( sim_pos, *conf.match_variants( sim_pos, templ_res ) , false );

				debug_assert( reference.residue( sim_pos ).is_lower_terminus() == pose.residue( sim_pos ).is_lower_terminus() );
				debug_assert( reference.residue( sim_pos ).is_upper_terminus() == pose.residue( sim_pos ).is_upper_terminus() );

			} catch ( EXCN_Env_Security_Exception& e ) {
				std::ostringstream ss;

				ss << this->get_name() << " failed trying to replace the simulation "
					<< pose.residue( sim_pos ).name3() << sim_pos << "(" << pose.residue( sim_pos ).natoms() << " atoms) with the template "
					<< templ().residue( templ_pos ).name3() << templ_pos <<  " ( "
					<< templ().residue( templ_pos ).natoms() << " atoms ).  ";
				if ( pose.residue( sim_pos ).name3() != templ().residue( templ_pos ).name3() ) {
					ss << "The problem is probably that the sequence shift is off. Surrounding sequences: " << std::endl
						<< "template:   ";

					for ( Size k = std::max( Size( 1 ), sim_origin().at( sim_pos ) - 5 ); k <= std::min( templ().size(), sim_origin().at( sim_pos ) + 5 ); ++k ) { ss << templ().residue( k ).name1(); }
					ss << std::endl
						<< "simulation: ";
					for ( Size k = std::max( Size( 1 ), sim_pos - 5 ); k <= std::min( pose.size(), sim_pos + 5 ); ++k ) { ss << pose.residue( k ).name1(); }
					ss << std::endl;
				} else if ( ( pose.is_centroid() && !templ().is_centroid() ) ||
						( pose.is_fullatom() && !templ().is_fullatom() ) ) {
					ss << "The problem is that the input is fullatom and it's being a applied to a centroid pose (or vice versa)--"
						<< "note the differing atom counts. If so, the 'apply_to_template' option to apply a "
						<< "SwitchResidueTypeSetMover to the template. There could also be a problem with differing variant." << std::endl;
				} else {
					ss << "The problem is that there are a differing number of atoms between the two residues, and replace_residue "
						<< "doesn't know how to handle the difference (differing degrees of freedom and all). This is usually from "
						<< "differing variants." << std::endl
						<< "Simulation: " << utility::to_string( pose.residue( sim_pos ).type().properties().get_list_of_variants() ) << std::endl
						<< "Template: " << utility::to_string( templ().residue( templ_pos ).type().properties().get_list_of_variants() ) << std::endl;
				}

				throw utility::excn::EXCN_BadInput( ss.str() );
			}
		}
	}

	if ( tr.Debug.visible() ) {
		tr.Debug << pose.fold_tree() << std::endl;
		tr.Debug << templ().fold_tree() << std::endl;
	}

	// correct mainchain connections

	for ( Size sim_pos = 1; sim_pos <= pose.size(); ++sim_pos ) {

		// If the residue before this is not in the region
		if ( sim_origin().find( sim_pos - 1 ) == sim_origin().end() &&
				sim_origin().find( sim_pos ) != sim_origin().end() ) {
			core::Size const templ_pos = sim_origin().at( sim_pos );

			bool lower_connect = ( sim_pos > 1  &&
				!pose.residue( sim_pos ).is_lower_terminus() &&
				!pose.fold_tree().is_cutpoint( sim_pos - 1 ) );
			try {
				if ( lower_connect &&
						( templ_pos - 1 < 1 ||
						templ().residue( templ_pos ).is_lower_terminus() ||
						templ().fold_tree().is_cutpoint( templ_pos -1 ) ) ) {
					// Here, the template doesn't have valid coordinates to connect to, so we have to do it
					// using reference values.
					tr.Debug << "fixing lower connection for " << sim_pos << " using non-template values." << std::endl;
					fix_mainchain_connect( pose, sim_pos, reference, sim_pos );
				} else if ( lower_connect ) {
					tr.Debug << "fixing lower connection for " << sim_pos << std::endl;
					fix_mainchain_connect( pose, sim_pos, templ(), templ_pos );
				} else {
					tr.Debug << "NOT fixing lower connection for " << sim_pos
						<< " ( lower_terminus : " << ( pose.residue( sim_pos ).is_lower_terminus() ? "T" : "F" )
						<< ", cutpoint " << ( pose.fold_tree().is_cutpoint( sim_pos - 1 ) ? "T" : "F" )
						<< std::endl;
				}
			} catch( core::environment::EXCN_Env_Exception& e ) {
				std::ostringstream ss;
				ss << this->get_name() << " couldn't repair the chunk's N-terminal connection to the mainchain at " << lower_connect
					<< " because it does not have DoF access at that position (which is one outside the selected region). "
					<< "Some other claimer must be claiming EXCLUSIVE at this position.";
				e.add_msg( ss.str() );
				throw;
			}
		}
		// If the residue after this is not in the region
		if ( sim_origin().find( sim_pos + 1 ) == sim_origin().end() &&
				sim_origin().find( sim_pos ) != sim_origin().end() ) {
			core::Size const templ_pos = sim_origin().at( sim_pos );

			bool upper_connect = ( sim_pos < pose.size() &&
				!pose.residue( sim_pos ).is_upper_terminus() &&
				!pose.fold_tree().is_cutpoint( sim_pos ) );

			try {
				if ( upper_connect &&
						( templ_pos + 1 > templ().size() ||
						templ().residue( templ_pos+1 ).is_upper_terminus() ||
						templ().fold_tree().is_cutpoint( templ_pos ) ) ) {
					// Here, the template doesn't have valid coordinates to connect to, so we have to do it
					// using reference values.
					tr.Debug << "fixing upper connection for " << sim_pos << " using non-template values." << std::endl;
					fix_mainchain_connect( pose, sim_pos+1, reference, sim_pos+1 );
				} else if ( upper_connect ) {
					tr.Debug << "fixing upper connection for " << sim_pos << std::endl;
					fix_mainchain_connect( pose, sim_pos+1, templ(), templ_pos+1 );
				} else {
					tr.Debug << "NOT fixing upper connection for " << sim_pos << std::endl;
				}
			} catch( core::environment::EXCN_Env_Exception& e ) {
				std::ostringstream ss;
				ss << this->get_name() << " couldn't repair the chunk's C-terminal connection to the mainchain at " << upper_connect
					<< " because it does not have DoF access at that position (which is one outside the selected region). "
					<< "Some other claimer must be claiming EXCLUSIVE at this position.";
				e.add_msg( ss.str() );
				throw;
			}
		}
	}


}


void RigidChunkCM::apply( core::pose::Pose& ){
}

loops::Loops RigidChunkCM::select_parts( loops::Loops const& rigid_core, core::Size random_grow_loops_by ) {
	loops::Loops current_rigid_core;

	if ( rigid_core.size() < 1 ) {
		tr.Error << "Given rigid core had size < 1. Check your loop definitions." << std::endl;
		throw utility::excn::EXCN_BadInput( "Given rigid core had size < 1. Check your loop definitions." );
	}

	for ( Size attempts = 1; attempts <= 50 && current_rigid_core.size() != 0; ++attempts ) {
		for ( loops::Loops::const_iterator it = rigid_core.begin(); it != rigid_core.end(); ++it ) {
			if ( numeric::random::rg().uniform() >= it->skip_rate() )  {
				current_rigid_core.push_back( *it );
			}
		}
	}

	if ( current_rigid_core.size() == 0 ) {
		current_rigid_core = rigid_core;
	}

	if ( random_grow_loops_by > 0 ) {
		core::Size nres( current_rigid_core[ current_rigid_core.size() ].stop() + 200 ); //it doesn't matter for this where exactly nres is.
		loops::Loops loops( current_rigid_core.invert( nres ) );
		loops.grow_all_loops( nres, random_grow_loops_by );
		tr.Info << "Enlarged loops: " << std::endl;
		tr.Info << loops << std::endl;
		current_rigid_core = loops.invert( nres );
	}

	return current_rigid_core;
}

void RigidChunkCM::sim_selector(
	core::select::residue_selector::ResidueSelectorCOP selector
) {
	if ( Parent::state_check( __FUNCTION__, ( selector.get() == sim_selector_.get() ) ) ) {
		sim_selector_ = selector;
	}
}

void RigidChunkCM::templ_selector(
	core::select::residue_selector::ResidueSelectorCOP selector
) {
	if ( Parent::state_check( __FUNCTION__, ( selector.get() == templ_selector_.get() ) ) ) {
		templ_selector_ = selector;
	}
}

core::select::residue_selector::ResidueSelectorCOP RigidChunkCM::sim_selector() const {
	debug_assert( sim_selector_ );
	return sim_selector_;
}

core::select::residue_selector::ResidueSelectorCOP RigidChunkCM::templ_selector() const {
	debug_assert( templ_selector_ );
	return templ_selector_;
}

// XRW TEMP std::string RigidChunkCM::get_name() const {
// XRW TEMP  std::ostringstream ss;
// XRW TEMP  ss << "RigidChunkCM(";
// XRW TEMP
// XRW TEMP  if ( xml_name_ == "" ) {
// XRW TEMP   ss << sim_selector()->get_name() << "+" << templ_selector()->get_name();
// XRW TEMP   return ss.str();
// XRW TEMP  } else {
// XRW TEMP   ss << xml_name_;
// XRW TEMP  }
// XRW TEMP
// XRW TEMP  ss << ")";
// XRW TEMP  return ss.str();
// XRW TEMP }

void RigidChunkCM::passport_updated() {
	if ( !has_passport() ) {
		templ_target_.clear();
		sim_origin_.clear();
	}
}

moves::MoverOP RigidChunkCM::clone() const {
	return moves::MoverOP( new RigidChunkCM( *this ) );
}

std::string RigidChunkCM::get_name() const {
	return mover_name();
}

std::string RigidChunkCM::mover_name() {
	return "RigidChunkCM";
}

void RigidChunkCM::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute(
		"name", xs_string,
		"Unique name for this client mover");
	attlist + XMLSchemaAttribute(
		"selector", xs_string,
		"Residue selector specifying the region over which this client mover will be applied");
	attlist + XMLSchemaAttribute(
		"apply_to_template", xs_string,
		"Comma-separated list of movers to apply to the template region");
	attlist + XMLSchemaAttribute(
		"template", xs_string,
		"PDB file to be used as a template for the specified region. If the argument 'INPUT' is supplied instead, it will fix the coordinates to their current positions");
	attlist + XMLSchemaAttribute(
		"region_file", xs_string,
		"Loops file specifying the region of the template pdb to use as a template for the pose");
	attlist + XMLSchemaAttribute(
		"region_selector", xs_string,
		"Residue selector specifying the region of the template pdb to use as a template for the pose");
	attlist + XMLSchemaAttribute(
		"region", xs_string,
		"String specifying indices of residues to be used as a template");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Holds a specified region of the pose constant, fixed to the coordinates in the template pdb file.",
		attlist );
}

std::string RigidChunkCMCreator::keyname() const {
	return RigidChunkCM::mover_name();
}

protocols::moves::MoverOP
RigidChunkCMCreator::create_mover() const {
	return protocols::moves::MoverOP( new RigidChunkCM );
}

void RigidChunkCMCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RigidChunkCM::provide_xml_schema( xsd );
}


} // abscript
} // abinitio
} // protocols
