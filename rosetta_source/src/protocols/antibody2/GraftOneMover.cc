// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody2/GraftOneMover.cc
/// @brief grafts a cdr onto the template of an antibody framework
/// @detailed
/// @author Jianqing Xu (xubest@gmail.com)

#include <protocols/antibody2/GraftOneMover.hh>
#include <core/conformation/Conformation.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/antibody2/AntibodyInfo.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.tmpl.hh>


#include <core/id/AtomID_Map.hh>
#include <core/pose/util.tmpl.hh>

#include <basic/Tracer.hh>


static basic::Tracer TRG("protocols.antibody2.GraftOneMover");

namespace protocols {
namespace antibody2 {
using namespace core;



GraftOneMover::GraftOneMover(){}

GraftOneMover::GraftOneMover( Size query_start, Size query_end, std::string template_name, scoring::ScoreFunctionOP scorefxn ) : Mover( "GraftOneMover" ), scorefxn_( scorefxn )
{
	query_start_ = query_start;
	query_end_ = query_end;
	set_default( template_name );
} // GraftOneMover default constructor



// GraftOneMover default destructor
GraftOneMover::~GraftOneMover() {}



void GraftOneMover::set_default( std::string template_name )
{
	std::string const path = basic::options::option[ basic::options::OptionKeys::in::path::path ]()[1];
	TRG << "Reading in template: " << path << template_name << ".pdb " << std::endl;
	core::import_pose::pose_from_pdb( template_pose_, path + template_name + ".pdb" );
	antibody2::AntibodyInfo ab_info( template_pose_, template_name );
	template_start_ = ab_info.current_start;
	template_end_ = ab_info.current_end;
	template_name_ = template_name;
} // GraftOneMover::set_default

std::string
GraftOneMover::get_name() const { return "GraftOneMover"; }

// copy ctor
GraftOneMover::GraftOneMover( GraftOneMover const & rhs ) {
    initForEqualOperatorAndCopyConstructor(*this, rhs);
}

///@brief assignment operator
GraftOneMover & GraftOneMover::operator=( GraftOneMover const & rhs ){
    //abort self-assignment
    if (this == &rhs) return *this;
    Mover::operator=(rhs);
    initForEqualOperatorAndCopyConstructor(*this, rhs);
    return *this;
}

void GraftOneMover::initForEqualOperatorAndCopyConstructor(GraftOneMover & lhs, GraftOneMover const & rhs) {

}



void GraftOneMover::apply( pose::Pose & pose_in )
{
	TRG<<"step1        I am here 7.4.1"<<std::endl;
	Size const nres( pose_in.total_residue() ); // Total residues
	Size query_size = ( query_end_ - query_start_ )+1;
	Size flank_size ( 5 );

	//set up packer
	pack::task::PackerTaskOP task;
	task = pack::task::TaskFactory::create_packer_task( pose_in );
	task->restrict_to_repacking();
	task->or_include_current( false );
	utility::vector1< bool > allow_repack( nres, false );
	for ( Size i = query_start_; i <= query_end_; ++i ) allow_repack[ i ] = true;
	task->restrict_to_residues( allow_repack );
	protocols::simple_moves::PackRotamersMoverOP packer = new protocols::simple_moves::PackRotamersMover( scorefxn_, task );

	// create a sub pose with  5 flanking residues on either side of CDR loop
	pose::Pose truncated_pose( pose_in, query_start_-flank_size, query_end_+flank_size );
	truncated_pose.conformation().delete_residue_range_slow( flank_size+1, ( query_size ) + flank_size );

	// create atom map for superimposing 2 flanking resiudes
	id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, template_pose_, id::BOGUS_ATOM_ID );

	for( Size start_stem = 1; start_stem < flank_size; ++start_stem ) {
		Size const ref_stem ( start_stem+1 ); // remember there are 5 flanking residues
		for( Size j=1; j <= 4; j++ ) {
			id::AtomID const id1( j, start_stem );
			id::AtomID const id2( j, ref_stem );
			atom_map[ id1 ] = id2;
		}
	}

	// start at the end of the actual loop
	for( Size end_stem = query_size+flank_size; end_stem < query_size+flank_size+4; ++end_stem ) {
		Size const ref_stem ( end_stem-query_size+1 );
		for( Size j=1; j <= 4; j++ ) {
			id::AtomID const id1( j, end_stem );
			id::AtomID const id2( j, ref_stem );
			atom_map[ id1 ] = id2;
		}
	}
	scoring::superimpose_pose( template_pose_, truncated_pose, atom_map );
	template_pose_.dump_pdb(template_name_);

	for ( Size i=query_start_; i <= query_end_; ++i ) {
		Size template_num ( i - (query_start_-5) );
		core::conformation::Residue const & source_rsd( pose_in.residue( i ) );
		core::conformation::Residue const & target_rsd( template_pose_.residue( template_num) );

		Size const natoms( source_rsd.natoms() );

		bool any_missing( false );
		id::AtomID_Mask missing( false );
		// dimension the missing-atom mask
		core::pose::initialize_atomid_map( missing, pose_in );

		if( source_rsd.name() != target_rsd.name() ) pose_in.set_secstruct( i, 'X' );
		for ( Size j=1; j<= natoms; ++j ) {
			std::string const & atom_name( source_rsd.atom_name(j) );
			if ( target_rsd.has( atom_name ) ) {
				pose_in.set_xyz( id::AtomID( source_rsd.atom_index( atom_name),i ), target_rsd.xyz( atom_name ) );
			} else {
				any_missing = true;
				missing[ id::AtomID( pose_in.residue_type(i).atom_index(source_rsd.atom_name(j)), i ) ] = true;
			}
		}

		if ( any_missing ) {
			pose_in.conformation().fill_missing_atoms( missing );
		}
	}
	packer->apply( pose_in );
	pose_in.dump_pdb(template_name_+"_graft");
} // GraftOneMover::apply



}  // namespace antibody2
}  // namespace protocols
