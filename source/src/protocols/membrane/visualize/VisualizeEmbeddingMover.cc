// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file     protocols/membrane/VisualizeEmbeddingMoverCreator.hh
/// @brief      Visualize Embedding normal and center with Virtual Residues
/// @details    Add a set of virtual residues as an additional chain to the
///    membrane pose. This tool is strictly for visualization of
///    the implicit membrane and should not be present in modeling.
///    Last Modified: 11/20/14
/// @author  JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_visualize_VisualizeEmbeddingMover_cc
#define INCLUDED_protocols_membrane_visualize_VisualizeEmbeddingMover_cc

// Unit Headers
#include <protocols/membrane/visualize/VisualizeEmbeddingMover.hh>
#include <protocols/membrane/visualize/VisualizeEmbeddingMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Project Headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.membrane.visualize.VisualizeEmbeddingMover" );

namespace protocols {
namespace membrane {
namespace visualize {

////////////////////
/// Constructors ///
////////////////////

/// @brief   Default Constructor
/// @details  Gets the embeddings from the pose and topology
VisualizeEmbeddingMover::VisualizeEmbeddingMover() : embeddings_( new protocols::membrane::geometry::Embedding() ) {
	register_options();
}

/// @brief   Constructor from embedding
/// @details  Visualizes defined embedding
VisualizeEmbeddingMover::VisualizeEmbeddingMover( protocols::membrane::geometry::EmbeddingOP embedding ) :
	embeddings_( embedding ) {
	register_options();
}

/// @brief Copy Constructor
/// @details Creates a deep copy of the visualize membrane mover class
VisualizeEmbeddingMover::VisualizeEmbeddingMover( VisualizeEmbeddingMover const & src ) :
	protocols::moves::Mover( src ),
	embeddings_( src.embeddings_ )
{}

/// @brief Assignment Operator
/// @details Overloads "=" assignemnt for deep copying
VisualizeEmbeddingMover &
VisualizeEmbeddingMover::operator=( VisualizeEmbeddingMover const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new VisualizeEmbeddingMover( *this ) );

}

/// @brief Destructor
VisualizeEmbeddingMover::~VisualizeEmbeddingMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
VisualizeEmbeddingMover::clone() const {
	return ( protocols::moves::MoverOP( new VisualizeEmbeddingMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
VisualizeEmbeddingMover::fresh_instance() const {
	return protocols::moves::MoverOP( new VisualizeEmbeddingMover() );
}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
VisualizeEmbeddingMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new VisualizeEmbeddingMover() );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
VisualizeEmbeddingMoverCreator::keyname() const {
	return VisualizeEmbeddingMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
VisualizeEmbeddingMoverCreator::mover_name() {
	return "VisualizeEmbeddingMover";
}

/////////////////////
/// Mover Methods ///
/////////////////////

void
VisualizeEmbeddingMover::apply( core::pose::Pose & pose ) {

	using namespace core;
	using namespace core::conformation;
	using namespace protocols::membrane::geometry;
	using namespace core::conformation::membrane;
	using namespace core::kinematics;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ OptionKeys::mp::visualize::embedding ].user() ) {

		TR << "Adding virtual residues to the pose to visualize embeddings of TMspans" << std::endl;

		// Check that I am a membrane pose
		if ( ! pose.conformation().is_membrane() ) {
			utility_exit_with_message("Cannot visualize a non-membrane pose!");
		}

		// Determine if the pose is fullatom (needed for residue typesets)
		bool fullatom = pose.is_fullatom();

		// get topology
		SpanningTopology topo = *pose.conformation().membrane_info()->spanning_topology();

		// get embedding object from pose and topology if not set previously
		if ( embeddings_->nspans() == 0 ) {
			EmbeddingOP emb( new Embedding( topo, pose ) );
			embeddings_ = emb;
		}
		embeddings_->show( TR );

		// create empty object to push embedding residues back into
		utility::vector1< ResidueOP > embedding_residues;

		// get thickness from MembraneInfo
		Size thickness = pose.conformation().membrane_info()->membrane_thickness();

		// grab residues out of embedding object
		for ( Size i = 1; i <= embeddings_->nspans(); ++i ) {

			// get center and normal POINTS
			Vector center = embeddings_->embedding( i )->center();
			Vector normal = center + embeddings_->embedding( i )->normal().normalize( thickness );

			// Create a new residue, append to list
			ResidueOP emb = create_embedding_virtual( center, normal, fullatom );

			// Append residue to list
			embedding_residues.push_back( emb );
		}

		// add total embedding
		Vector center_tot = embeddings_->total_embed()->center();
		Vector normal = embeddings_->total_embed()->normal().normalize( thickness );
		Vector normal_tot = center_tot + normal;

		// Create a new residue, append to list
		ResidueOP emb_tot = create_embedding_virtual( center_tot, normal_tot, fullatom );
		embedding_residues.push_back( emb_tot );

		// if first residue, create new chain
		bool is_first( true );

		// Append Residues to the pose
		for ( Size i = 1; i <= embedding_residues.size(); ++i ) {
			if ( is_first ) {
				pose.append_residue_by_jump( *embedding_residues[i], pose.total_residue(), "", "", true );
				is_first = false;
			} else {
				pose.append_residue_by_jump( *embedding_residues[i], pose.total_residue(), "", "", false );
			}
		}
	}
}

std::string
VisualizeEmbeddingMover::get_name() const {
	return "VisualizeEmbeddingMover";
}

//////////////////////
/// Helper Methods ///
//////////////////////

/// @brief Register Options with JD2
void
VisualizeEmbeddingMover::register_options() {

	using namespace basic::options;
	option.add_relevant( OptionKeys::mp::visualize::embedding );

}

/// @brief create virtual residue for visualization in PDB file
core::conformation::ResidueOP
VisualizeEmbeddingMover::create_embedding_virtual( core::Vector center, core::Vector normal, bool fullatom ) {

	using namespace core;
	using namespace core::conformation;
	using namespace core::chemical;

	// Grab the current residue typeset and create a new residue
	ResidueTypeSetCOP const & residue_set(
		core::chemical::ChemicalManager::get_instance()->residue_type_set( fullatom ? core::chemical::FA_STANDARD : core::chemical::CENTROID )
	);

	// Create a new Residue from rsd typeset of type MEM
	ResidueTypeCOP rsd_type( residue_set->get_representative_type_name3("EMB") );
	ResidueType const & membrane( *rsd_type );
	ResidueOP rsd( ResidueFactory::create_residue( membrane ) );

	// set something close to center for thickness, so it doesn't confuse the viewer
	Vector bogus( center.x()+0.001, center.y()+0.001, center.z()+0.001 );

	// Fill the Residue with normal/center info: AtomID, position
	rsd->set_xyz(1, center);
	rsd->set_xyz(2, normal);
	rsd->set_xyz(3, bogus);

	return rsd;
}


} // visualize
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_visualize_VisualizeEmbeddingMover_cc
