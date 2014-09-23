// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /devel/LoopExtend/LoopExtendMover.cc
/// @brief Extends a loop based on loop definition file with loop::extendby ALA after cutpoint using kinematic closure.
/// @author Daniel J. Mandell

// Unit Headers
#include <devel/loop_extend/LoopExtendMover.hh>

// Package Headers
#include <protocols/loops/Loops.hh>
//#include <protocols/loops/loops_main.hh> //loops_set_move_map

// Project Headers

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
// AUTO-REMOVED #include <core/conformation/util.hh>
// AUTO-REMOVED
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
#include <core/chemical/ChemicalManager.hh> //CENTROID, FA_STANDARD
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>

#include <protocols/loops/loops_main.hh>

#include <utility/vector1.hh>


//movers
// will need LoopMover_KIC

// option key includes
//#include <basic/options/keys/loops.OptionKeys.gen.hh>
//#include <basic/options/keys/run.OptionKeys.gen.hh>

using basic::T;

static thread_local basic::Tracer tr( "devel.LoopExtend.LoopExtendMover" );


namespace devel{
namespace loop_extend{



///@details Extends a loop
void LoopExtendMover::apply( core::pose::Pose & pose )
{
	using namespace core::chemical;
	using namespace core::conformation;
	using core::Real;
	using core::Size;
	//tr << "rocking loop:" << loop_ << std::endl;
	// setup the fold tree for extension
	core::kinematics::FoldTree f_new, f_orig=pose.fold_tree();
	protocols::loops::Loops loops; // fold_tree_from_loops needs loops even though we only use one loop
	loops.add_loop(loop_);
	protocols::loops::fold_tree_from_loops( pose, loops, f_new, true /* include terminal cutpoints */);
	pose.fold_tree( f_new );
	Size insert_point = loop_.cut()+1;
	// residue type set for making alanines from ResidueFactory
	ResidueTypeSetCOP residue_set
	( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

	for (Size i=loop_.start(); i<=loop_.stop(); i++) {
		for ( Size i = loop_.start(); i <= loop_.stop(); ++i ) {
			if ( i != loop_.start() )	pose.set_phi( i, init_phi_ );
			if ( i != loop_.stop() ) pose.set_psi( i, init_psi_ );
			if ( ( i != loop_.start() ) && ( i != loop_.stop() ) ) pose.set_omega( i, init_omega_ );
		}
	}

	for (Size i=1; i<=extend_len_; ++i) {
		tr << "Adding residue: " << i << " at position " << insert_point << std::endl;
		ResidueOP ala_rsd = ResidueFactory::create_residue( residue_set->name_map( "ALA" ) );
		//Real xpos = static_cast<Real> (ala_rsd->xyz(2).x());
		//tr << "extend res CA xpos before prepend: " << xpos << std::endl;
		// prepend the new alanine before the residue after the cutpoint
		pose.conformation().safely_prepend_polymer_residue_before_seqpos( *ala_rsd, insert_point, true );
		//xpos = static_cast<Real> (ala_rsd->xyz(2).x());
		//tr << "extend res CA xpos after prepend: " << xpos << std::endl;
		// set ideal omega angle at junction
		//extn.set_torsion( id::TorsionID( 1, id::BB, 3 ),  180.0);
	}

	return;
}//LoopDesignMover::apply

std::string
LoopExtendMover::get_name() const {
	return "LoopExtendMover";
}



}//LoopExtend
}//devel
