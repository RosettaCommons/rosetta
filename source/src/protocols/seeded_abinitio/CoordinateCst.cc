// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file
/// @brief allows to add coordinate constraints to the pose based on parsed spans or residues
///     residues specified are parsed at runtime not at parse time, so it can be used somewhat used with length changes
///     currently it only adds CA constraints, but hopefully at some time point others too
///
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

// Unit headers
#include <protocols/seeded_abinitio/CoordinateCst.hh>
#include <protocols/seeded_abinitio/CoordinateCstCreator.hh>

#include <core/scoring/func/Func.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <numeric/xyzVector.hh>

#include <core/scoring/func/HarmonicFunc.hh>
#include <core/types.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <utility/string_util.hh>
#include <boost/foreach.hpp>
#include <utility/tag/Tag.hh>
// C++ headers
#include <iostream>
#include <string>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>

#include <core/pose/Pose.hh>

//pose
#include <core/pose/Pose.fwd.hh>

#include <utility/vector0.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>


namespace protocols {
namespace seeded_abinitio {

using namespace core;
using namespace scoring::constraints;
using namespace protocols::moves;

static THREAD_LOCAL basic::Tracer TR( "protocols.seeded_abinitio.CoordinateCst" );
static THREAD_LOCAL basic::Tracer TR_debug( "CoordinateCst.Debug" );

std::string
CoordinateCstCreator::keyname() const
{
	return CoordinateCstCreator::mover_name();
}

protocols::moves::MoverOP
CoordinateCstCreator::create_mover() const {
	return protocols::moves::MoverOP( new CoordinateCst );
}

std::string
CoordinateCstCreator::mover_name()
{
	return "CoordinateCst";
}

CoordinateCst::~CoordinateCst() = default;

CoordinateCst::CoordinateCst() :
	protocols::moves::Mover( CoordinateCstCreator::mover_name() )
{
	stddev_ =  3.0;
	use_jumps_ = true;
	anchor_res_.clear();
	span_vector_.clear();
	unparsed_residue_.clear();
	jump_ = 0;
	anchor_atom_id_ = ""; // by default, use the same atom name as atom_id_ -- denote their sameness by setting this to the empty string
	atom_id_ = "CA";
}

protocols::moves::MoverOP
CoordinateCst::clone() const {
	return( protocols::moves::MoverOP( new CoordinateCst( *this ) ) );
}

protocols::moves::MoverOP
CoordinateCst::fresh_instance() const {
	return protocols::moves::MoverOP( new CoordinateCst );
}

///parse residues at run time, in case there was a lenght change
void parse_spans(
	pose::Pose & pose,
	utility::vector1 < std::pair < std::string, std::string > > const span_vector,
	std::set < core::Size > resi_collection
)
{

	for ( Size iter = 1 ; iter <= span_vector.size() ; ++ iter ) {
		TR.Debug <<"sanity check, span_vector[iter].first " <<span_vector[iter].first <<std::endl;
		core::Size const begin = core::pose::parse_resnum( span_vector[iter].first, pose ) ;
		core::Size const end   = core::pose::parse_resnum( span_vector[iter].second, pose );
		//runtime_assert( end > begin );
		//runtime_assert( begin>=1);
		//runtime_assert( end<=pose.total_residue() );

		for ( Size resi = begin; resi <= end ; resi ++ ) {
			resi_collection.insert( resi );
			TR  << "runtime parsed of span, residue: " << resi << std::endl;
		}
	}
}

utility::vector1< core::Size >
adjust_residues(
	pose::Pose & pose,
	std::string design_residues
){

	utility::vector1< std::string > const design_keys( utility::string_split( design_residues, ',' ) );
	utility::vector1< core::Size > design_res;

	BOOST_FOREACH ( std::string const key, design_keys ) {
		core::Size const resnum( core::pose::parse_resnum( key, pose ));
		//TR.Debug<<"design within seed, residue: "<< key <<", parsed: "<< resnum <<std::endl;
		design_res.push_back( resnum);
		TR.Debug<<"parsed "<<key<<std::endl;
	}
	//TR.Debug<<"runtime designable: " << design_res <<std::endl;
	return design_res;
}//end parsing design residues


///residue parsing at runtime
void adjust_single_residues(
	pose::Pose & pose,
	std::string single_residues,
	std::set < core::Size > resi_collection
){

	utility::vector1< std::string > const single_keys( utility::string_split( single_residues, ',' ) );

	BOOST_FOREACH ( std::string const key, single_keys ) {
		core::Size const resnum( core::pose::parse_resnum( key, pose ));
		resi_collection.insert( resnum );
	}
}

void add_coordinate_constraints(
	pose::Pose & pose,
	std::set < core::Size > const & constrain_residues,
	core::Size const anchor_resnum,
	core::Real const coord_sdev,
	std::string anchor_atom_name, // if this is the empty string, use the same string as "atom_name"
	std::string atom_name,
	core::scoring::func::HarmonicFuncOP & coord_cst_func
)
{

	using namespace core::scoring::constraints;
	using namespace core::conformation;
	//should this rather be a ConstraintSet?
	ConstraintCOPs cst;

	if ( anchor_atom_name == "" ) {
		anchor_atom_name = atom_name;
	}

	if ( anchor_resnum != 0 ) {
		TR<<"Anchor residue for the coordinate constraint is "<< anchor_resnum <<std::endl;

		if ( anchor_atom_name == "CB" ) {
			if ( pose.residue( anchor_resnum ).aa() == core::chemical::aa_gly ) {
				anchor_atom_name = "CA";
			}
		}
	}
	if ( ! pose.residue( anchor_resnum ).has( anchor_atom_name ) ) {
		utility_exit_with_message( "Could not add coordinate constraint to anchor residue atom " + anchor_atom_name + " on residue "
			+ utility::to_string( anchor_resnum ) + " (" + pose.residue( anchor_resnum ).name()
			+ ") since that residue does not have an atom with that name" );
	}
	core::id::AtomID const anchor_atom( core::id::AtomID( pose.residue( anchor_resnum ).atom_index( anchor_atom_name ), anchor_resnum ) );

	if ( !coord_cst_func ) {
		coord_cst_func = core::scoring::func::HarmonicFuncOP( new core::scoring::func::HarmonicFunc( 0.0, 0.0 ) );
	}
	coord_cst_func->sd( coord_sdev );

	BOOST_FOREACH ( core::Size res, constrain_residues ) {
		// safety check.
		if ( ! pose.residue( res ).has( atom_name ) ) {
			utility_exit_with_message( "Could not add coordinate constraint to " + atom_name + " on residue " + utility::to_string( res ) + " (" +
				pose.residue( res ).name() + ") since that residue does not have an atom with that name" );
		}
		TR<<"Coordinate-constraining residue " << pose.residue( res ).name() << " " << res <<std::endl;
		///core::chemical::ResidueType rsd_type( pose.residue( res ).type() );
		Size atomindex =  pose.residue( res ).atom_index( atom_name );
		cst.push_back( core::scoring::constraints::ConstraintOP( new CoordinateConstraint(
			core::id::AtomID( atomindex, res ),
			anchor_atom,
			pose.residue( res ).xyz( atomindex ),
			coord_cst_func) ));
		cst = pose.add_constraints( cst );
	}
}

//option to use coordinate constraints on CA or on the functional group!(todo)
//for right now using CA

void
CoordinateCst::apply( pose::Pose & pose )
{

	using namespace core;

	//collect all residues that should have coordinate constraints
	std::set < core::Size > constrain_residues_set;

	if ( unparsed_residue_ != "" ) {
		adjust_single_residues( pose, unparsed_residue_, constrain_residues_set );
	}

	parse_spans( pose, span_vector_ , constrain_residues_set );

	TR.Debug << "constrain residue set size : " << constrain_residues_set.size() <<" -- if this is 0, you did NOT read in any residues "<< std::endl;

	Size anchor_res = 0;
	TR.Debug<< "anchor " << anchor_res_ <<std::endl;

	if ( !use_jumps_ && anchor_res_ != "" ) {
		//anchor needs to be parsed too, could do that explicitly or as a set, in case anyone wanted at some point more than one anchor...
		//std::set < core::Size > anchor_set;
		//adjust_single_residues( pose, anchor_res_ , anchor_set);
		utility::vector1 < core::Size > anchor_res_vec( adjust_residues( pose, anchor_res_ ) );
		//for now just the first elemen
		anchor_res = anchor_res_vec[1];
		TR << "-- anchor residue: " <<  anchor_res << " --" << std::endl;
	}

	if ( use_jumps_ ) {
		if ( pose.fold_tree().jump_edge( jump_ ).start()  == 0 ) {
			TR<<"no jump detected, defaulting anchor to 1" << std::endl;
			anchor_res = 1;
		}
		if ( constrain_residues_set.empty() ) { //size() == 0 ) {
			constrain_residues_set.insert( pose.fold_tree().jump_edge( jump_ ).stop() );
			TR<< "adding constraints to the downstream jump atom" << std::endl;
		}

		anchor_res = (pose.fold_tree().jump_edge( jump_ ).start());
		TR<< "setting anchor residue to upstream jump residue: " << anchor_res_ << std::endl;
	}


	core::scoring::func::HarmonicFuncOP coord_cst_func( new core::scoring::func::HarmonicFunc( 0.0, 0.0 ) );
	add_coordinate_constraints( pose, constrain_residues_set, anchor_res, stddev_, anchor_atom_id_, atom_id_, coord_cst_func );

}//end apply


std::string
CoordinateCst::get_name() const {
	return CoordinateCstCreator::mover_name();
}


void
CoordinateCst::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	Movers_map const &,
	Pose const & /*pose*/)
{

	TR<<"CoordinateCst mover has been invoked"<<std::endl;
	stddev_ = tag->getOption<core::Real>( "stddev", 0.5);

	TR<<"setting constraint standard deviation to "<< stddev_<< std::endl;
	use_jumps_ = true;

	atom_id_ = tag->getOption< std::string >( "atom" , "CA" );

	// if not provided, set anchor_atom_id_ to the empty string which will signal to the
	// add_coordinate_constraints function above that it should have the same value as the atom_id_ string.
	anchor_atom_id_ = tag->getOption< std::string > ("anchor_atom", "" );

	if ( tag->hasOption( "anchor_res" ) ) {
		anchor_res_ = tag->getOption< std::string >( "anchor_res" );
		TR<< "user specified following residue as anchor atom for the coordinate constraints: " << anchor_res_ << std::endl;

		if ( anchor_res_ != ""  ) {
			TR<< "WARNING, if jump is specified, it will not be used, since the specification of an anchor residue will override it" <<std::endl;
			use_jumps_ = false;
		}//end excluding jump
	}//end anchor

	if ( tag->hasOption( "jump" ) ) {
		jump_ = tag->getOption< core::Size >( "jump", 0 );
	}

	if ( use_jumps_ ) {
		TR<< "using jump atoms for anchor, that is the first one, for jump #: " << jump_ << std::endl;
	}

	//parsing option for CA or CB to use or whether just functional groups!

	//parsing branch tags
	utility::vector0< TagCOP > const & branch_tags( tag->getTags() );

	BOOST_FOREACH ( TagCOP const btag, branch_tags ) {

		if ( btag->getName() == "Span" || btag->getName() == "span" ||  btag->getName() == "Seeds" ) {
			std::string const beginS( btag->getOption<std::string>( "begin" ) );
			std::string const endS( btag->getOption<std::string>( "end" ) );
			std::pair <std::string,std::string> seedpair;
			seedpair.first  = beginS;
			TR.Debug<<"parsing loop verteces: " << beginS << " " <<endS <<std::endl;
			seedpair.second = endS;
			span_vector_.push_back( seedpair );

		}//end seeds or chain

		if ( btag->getName() == "residue" ) {
			unparsed_residue_ = tag->getOption< std::string >("residue");
		}//end residue tags

	}//end b-tags
}//end parse my tag

}//seeded_abinitio/
}//protocol
