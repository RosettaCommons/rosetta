// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//////////////////////////////////////////////////////////////////////
///
/// @brief
/// Compare the sequences between a native and designed protein
///
/// @details
/// This is an implementation taken from Ron Jacak, Douglas Renfrew, Matt O Mera.
/// The main function that is called is the get_sequence_recovery() function. You can
/// pass this function a list of native pdbs and designed pdbs, or just 1 native and 1
/// designed pdb. The sequence recovery will be output in a file called sequencerecovery.txt
/// along with a substitution matrix in a file called submatrix.txt
///
///
/// @references
/// "Native sequences are close to optimal" paper
///
///
/// @author
/// Ron Jacak,
/// Douglas Renfrew (renfrew@nyu.edu) ( added rotamer recovery, cleanup )
/// Steven Combs (moved it into a general use class)
///
/////////////////////////////////////////////////////////////////////////

// Unit headers
//#include <devel/init.hh>

//project Headers
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <protocols/toolbox/pose_metric_calculators/SequenceComparison.hh>


// Utility Headers
#include <basic/Tracer.hh>
#include <basic/prof.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>


// Numeric Headers

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>

// C++ headers
#include <sstream>
#include <map>

//Auto Headers
#include <core/conformation/PointGraphData.hh>
#include <core/graph/UpperEdgeGraph.hh>
#include <utility/vector0.hh>


namespace protocols{
namespace toolbox{
namespace pose_metric_calculators{

static thread_local basic::Tracer TR( "seqrecovery" );

using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace ObjexxFCL::format;


/// @brief load custom TaskOperations according to an xml-like utility::tag file
core::pack::task::TaskFactoryOP SequenceComparison::setup_tf( core::pack::task::TaskFactoryOP task_factory_ ) {
	using namespace core::pack::task::operation;

		task_factory_->push_back( TaskOperationCOP( new pack::task::operation::InitializeFromCommandline ) );


	return task_factory_;

}

/// @brief return the set of residues that are designable based given pose
std::set< core::Size > SequenceComparison::fill_designable_set( core::pose::Pose & pose, core::pack::task::TaskFactoryOP & tf ) {

	//we need to score the pose for many of the task operations passed from cmd line
	std::set< Size > designable_set;
	core::pack::task::PackerTaskOP design_task( tf->create_task_and_apply_taskoperations( pose ) );

#ifndef NDEBUG
	TR<< "Task is: \n" << *(design_task)  << std::endl;
#endif

	// iterate over all residues
	for ( Size ii = 1; ii<= design_task->total_residue(); ++ii ) {
		if( design_task->being_designed( ii ) )
			designable_set.insert( ii );
	}

	return designable_set;

}


/// @brief helper method which uses the tenA nb graph in the pose object to fill a vector with nb counts
void SequenceComparison::fill_num_neighbors( pose::Pose & pose, utility::vector1< core::Size > & num_nbs ) {

	using core::conformation::PointGraph;
	using core::conformation::PointGraphOP;

	PointGraphOP pg( new PointGraph ); // create graph
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg ); // create vertices
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, 10.0 /* ten angstrom distance */ ); // create edges

	num_nbs.resize( pose.n_residue(), 0 );
	for ( core::Size ii=1; ii <= pose.total_residue(); ++ii ) {

		// a PointGraph is a typedef of UpperEdgeGraph< PointGraphVertexData, PointGraphEdgeData >
		// so any of the method in UpperEdgeGraph should be avail. here. The UpperEdgeGraph provides access to nodes
		// via a get_vertex() method, and each vertex can report back how many nbs it has.
		// So something that before was really complicated (nb count calculation) is done in <10 lines of code.
		// the assumption we're making here is that a pose residue position ii is the same index as the point graph vertex
		// that is indeed the case if you look at what the function residue_point_graph_from_pose().
		num_nbs[ ii ] = pg->get_vertex(ii).num_neighbors_counting_self();
	}

	return;
}


/// @brief iterates over all designed positions and determines identity to native. outputs recoveries to file.
void SequenceComparison::measure_sequence_recovery( utility::vector1<core::pose::Pose> & native_poses, utility::vector1<core::pose::Pose> & redesign_poses ) {

	// setup main arrays used for calculation
	utility::vector1< core::Size > n_correct( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_native( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_designed( chemical::num_canonical_aas, 0 );

	utility::vector1< core::Size > n_correct_core( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_native_core( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_designed_core( chemical::num_canonical_aas, 0 );

	utility::vector1< core::Size > n_correct_surface( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_native_surface( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_designed_surface( chemical::num_canonical_aas, 0 );

	ObjexxFCL::FArray2D_int sub_matrix( chemical::num_canonical_aas, chemical::num_canonical_aas, 0 );

	Size n_correct_total(0); Size n_total(0);
	Size n_correct_total_core(0); Size n_total_core(0);
	Size n_correct_total_surface(0); Size n_total_surface(0);

	Size surface_exposed_cutoff(surface_exposure_);
	Size core_cutoff(core_cutoff_);

	// iterate through all the structures
	utility::vector1< core::pose::Pose >::iterator native_itr( native_poses.begin() ), native_last( native_poses.end() );
	utility::vector1< core::pose::Pose >::iterator redesign_itr( redesign_poses.begin() ), redesign_last( redesign_poses.end() );

	while( ( native_itr != native_last ) && (redesign_itr != redesign_last ) ) {

		// get local copies of the poses
		core::pose::Pose native_pose( *native_itr );
		core::pose::Pose redesign_pose( *redesign_itr );

		// figure out the task & neighbor info
		core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory );
		std::set< Size > design_set;
		utility::vector1< core::Size > num_neighbors;

		// setup what residues we are going to look at...
		setup_tf( task_factory );
		design_set = fill_designable_set( native_pose, task_factory );
		fill_num_neighbors( native_pose, num_neighbors );

		// record native sequence
		// native_sequence vector is sized for the WHOLE pose not just those being designed
		// it doesn't matter because we only iterate over the number of designed positions
		Size const nres( native_pose.total_residue() );
		utility::vector1< chemical::AA > native_sequence( nres );

		// iterate over designable positions
		for ( std::set< core::Size >::const_iterator it = design_set.begin(), end = design_set.end(); it != end; ++it ) {

			if ( ! native_pose.residue(*it).is_protein() ) {
				native_sequence[ *it ] = chemical::aa_unk;
				continue;
			}

			native_sequence[ *it ] = native_pose.residue( *it ).aa();
			n_native[ native_pose.residue(*it).aa() ]++;

			//determine core/surface
			if ( num_neighbors[*it] >= core_cutoff ) {
				n_native_core[ native_pose.residue(*it).aa() ]++;
				n_total_core++;
			}

			if ( num_neighbors[*it] < surface_exposed_cutoff ) {
				n_native_surface[ native_pose.residue(*it).aa() ]++;
				n_total_surface++;
			}

		} // end finding native seq

		/// measure seq recov
		for ( std::set< core::Size >::const_iterator it = design_set.begin(), end = design_set.end(); it != end; ++it ) {

			// don't worry about recovery of non-protein residues
			if ( redesign_pose.residue( *it ).is_protein() ) {
				n_total++;

				// increment the designed count
				n_designed[ redesign_pose.residue(*it).aa() ]++;

				if ( num_neighbors[*it] >= core_cutoff ) { n_designed_core[ redesign_pose.residue(*it).aa() ]++; }
				if ( num_neighbors[*it] < surface_exposed_cutoff ) { n_designed_surface[ redesign_pose.residue(*it).aa() ]++; }

				// then check if it's the same
				if ( native_sequence[ *it ] == redesign_pose.residue(*it).aa() ) {
					n_correct[ redesign_pose.residue(*it).aa() ]++;

					if ( num_neighbors[*it] >= core_cutoff ) {
						n_correct_core[ redesign_pose.residue(*it).aa() ]++;
						n_correct_total_core++;
					}
					if ( num_neighbors[*it] < surface_exposed_cutoff ) {
						n_correct_surface[ redesign_pose.residue(*it).aa() ]++;
						n_correct_total_surface++;
					}
					n_correct_total++;
				}

				// set the substitution matrix for this go round
				sub_matrix( native_pose.residue(*it).aa(), redesign_pose.residue(*it).aa() )++;
			}

		} // end measure seq reovery

		// increment iterators
		native_itr++; redesign_itr++;
	}

	// open sequence recovery file stream
	utility::io::ozstream outputFile( "sequencerecovery.txt" ) ;

	// write header
	outputFile << "Residue\tNo.correct core\tNo.native core\tNo.designed core\tNo.correct/ No.native core\tNo.correct/ No.designed core\t"
				 << "No.correct\tNo.native\tNo.designed\tNo.correct/ No.native\tNo.correct/ No.designed\t"
				 << "Residue\tNo.correct surface\tNo.native surface\tNo.designed surface\tNo.correct/ No.native\tNo.correct/ No.designed" << std::endl;

	// write AA data
	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {

		outputFile << chemical::name_from_aa( chemical::AA(ii) ) << "\t"
					<< n_correct_core[ ii ] << "\t" << n_native_core[ ii ] << "\t" << n_designed_core[ ii ] << "\t";

		if ( n_native_core[ii] != 0 ) outputFile << F(4,2, (float)n_correct_core[ii]/n_native_core[ii] ) << "\t";
		else outputFile << "---\t";
		if ( n_designed_core[ii] != 0 ) outputFile << F(4,2, (float)n_correct_core[ii]/n_designed_core[ii] ) << "\t";
		else outputFile << "---\t";

		// debug
		//if ( n_native_core[ii] != 0 ) std::cout << F(4,2, (float)n_correct_core[ii]/n_native_core[ii] ) << "\t";
		//if ( n_designed_core[ii] != 0 ) std::cout << F(4,2, (float)n_correct_core[ii]/n_designed_core[ii] ) << "\t";

		outputFile << n_correct[ ii ] << "\t" << n_native[ ii ] << "\t" << n_designed[ ii ] << "\t";
		if ( n_native[ii] != 0 ) outputFile << F(4,2, (float)n_correct[ii]/n_native[ii] ) << "\t";
		else outputFile << "---\t";
		if ( n_designed[ii] != 0 ) outputFile << F(4,2, (float)n_correct[ii]/n_designed[ii] ) << "\t";
		else outputFile << "---\t";

		// debug
		//if ( n_native[ii] != 0 ) std::cout << F(4,2, (float)n_correct[ii]/n_native[ii] ) << "\t";
		//if ( n_designed[ii] != 0 ) std::cout << F(4,2, (float)n_correct[ii]/n_designed[ii] ) << "\t";

		outputFile << chemical::name_from_aa( chemical::AA(ii) ) << "\t"
							 << n_correct_surface[ ii ] << "\t" << n_native_surface[ ii ] << "\t" << n_designed_surface[ ii ] << "\t";

		if ( n_native_surface[ii] != 0 ) outputFile << F(4,2, (float)n_correct_surface[ii]/n_native_surface[ii] ) << "\t";
		else outputFile << "---\t";
		if ( n_designed_surface[ii] != 0 ) outputFile << F(4,2, (float)n_correct_surface[ii]/n_designed_surface[ii] ) << "\t";
		else outputFile << "---\t";

		// debug
		//if ( n_native_surface[ii] != 0 ) std::cout << F(4,2, (float)n_correct_surface[ii]/n_native_surface[ii] ) << "\t";
		//if ( n_designed_surface[ii] != 0 ) std::cout << F(4,2, (float)n_correct_surface[ii]/n_designed_surface[ii] ) << "\t";

		outputFile << std::endl;
	}

	// write totals
	outputFile << "Total\t"
				<< n_correct_total_core << "\t" << n_total_core << "\t\t" << F(5,3, (float)n_correct_total_core/n_total_core ) << "\t\t"
				<< n_correct_total << "\t" << n_total << "\t\t" << F(5,3, (float)n_correct_total/n_total ) << "\t\tTotal\t"
				<< n_correct_total_surface << "\t" << n_total_surface << "\t\t" << F(5,3, (float)n_correct_total_surface/n_total_surface )
				<< std::endl;


	// output the sequence substitution file
	utility::io::ozstream matrixFile( "submatrix.txt" ) ; //defaults to submatrix.txt

	// write the header
	matrixFile << "AA_TYPE" << "\t" ;
	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
		matrixFile << "nat_"<<chemical::name_from_aa( chemical::AA(ii) ) << "\t";
	}
	matrixFile<<std::endl;

	// now write the numbers
	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) { //redesigns
		matrixFile << "sub_" << chemical::name_from_aa( chemical::AA(ii) );
		for ( Size jj = 1; jj <= chemical::num_canonical_aas; ++jj ) { //natives
			//std::cout<<"Native: "<< jj << " Sub: " << ii << "  Value: "<<sub_matrix( jj, ii ) << std::endl;
			matrixFile<< "\t" << sub_matrix( jj, ii );
		}
		matrixFile << std::endl;
	}


}


void SequenceComparison::get_sequence_recovery(core::pose::Pose & native, core::pose::Pose & designed){
	utility::vector1<core::pose::Pose> native_poses;
	utility::vector1<core::pose::Pose> redesign_poses;

	native_poses.push_back(native);
	redesign_poses.push_back(designed);

	get_sequence_recovery(native_poses, redesign_poses);
}

//@brief main method for the sequence recovery protocol
void SequenceComparison::get_sequence_recovery(utility::vector1<core::pose::Pose> & native_poses, utility::vector1<core::pose::Pose> & redesign_poses){

	if ( native_poses.size() != redesign_poses.size() ) {
		utility_exit_with_message( "Size of native pdb list file   does not equal size of redesign pdb list! \n" );
	}

	TR << "Measuring sequence recovery" << std::endl;
	measure_sequence_recovery( native_poses, redesign_poses );

}


}
}
}
