// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/carbohydrates/glycan_relax.cc
/// @brief Application for Glycan Relax.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) and Jason W. Labonte (JWLabonte@gmail.com)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>
#include <protocols/moves/Mover.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/string_util.hh>

static THREAD_LOCAL basic::Tracer TR("glycan_info");


void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	option.add_relevant( in::file::s );
	option.add_relevant( in::file::l );

}

// Class Definitions //////////////////////////////////////////////////////////

/// @brief  Class to print out info specific to glycan poses.  Will be expanded as needed.
class GlycanInfoMover : public protocols::moves::Mover {
public:  // Standard methods
	/// @brief  Default constructor.
	GlycanInfoMover() : protocols::moves::Mover()
	{
		init();
	}

	/// @brief  Copy constructor.
	GlycanInfoMover( GlycanInfoMover const & object_to_copy ) : Mover( object_to_copy )
	{
		copy_data( *this, object_to_copy );
	}

	// Assignment operator
	GlycanInfoMover &
	operator=( GlycanInfoMover const & object_to_copy )
	{
		// Abort self-assignment.
		if ( this == &object_to_copy ) {
			return *this;
		}

		protocols::moves::Mover::operator=( object_to_copy );
		copy_data( *this, object_to_copy );
		return *this;
	}

	// Destructor
	virtual ~GlycanInfoMover() {}


public:  // Standard Rosetta methods



	/// @brief  Generate string representation of DockGlycansProtocol for debugging purposes.
	virtual
	void
	show( std::ostream & output=std::cout ) const
	{
		protocols::moves::Mover::show( output );  // name, type, tag
	}


	// Mover methods
	/// @brief  Return the name of the Mover.
	virtual
	std::string
	get_name() const
	{
		return type();
	}

	virtual
	protocols::moves::MoverOP
	clone() const
	{
		return protocols::moves::MoverOP( new GlycanInfoMover( *this ) );
	}

	virtual
	protocols::moves::MoverOP
	fresh_instance() const
	{
		return protocols::moves::MoverOP( new GlycanInfoMover );
	}

	std::string
	get_attachment_point_string( core::pose::Pose const & pose, core::Size resnum){
		using utility::to_string;
		using namespace core::chemical::carbohydrates;

		CarbohydrateInfoCOP info = pose.residue(resnum).carbohydrate_info();
		std::string outstring = "";
		std::string attach = "_->";

		if ( info->mainchain_glycosidic_bond_acceptor() ) {
			outstring = attach + to_string(info->mainchain_glycosidic_bond_acceptor());
		}

		for ( uint i = 1; i <= info->n_branches(); ++i ) {
			outstring = outstring + "," +attach + to_string( info->branch_point( i ));
		}
		return outstring;
	}
	/// @brief  Apply the corresponding protocol to <pose>.
	virtual
	void
	apply( core::pose::Pose & pose )
	{
		using namespace core::pose::carbohydrates;

		core::Size protein_branches = 0;
		core::Size carbohydrate_residues = 0;
		for ( core::Size resnum = 1; resnum <= pose.size(); ++resnum ) {
			if ( pose.residue( resnum ).is_carbohydrate() ) {
				std::string attachment_points = get_attachment_point_string( pose, resnum);
				core::Size parent_res = find_seqpos_of_saccharides_parent_residue( pose.residue( resnum ));
				bool bp = pose.residue( resnum ).is_branch_point();


				std::cout << "Carbohydrate: "<< resnum  << " Parent: " << parent_res << " BP: "<<bp << " CON: "
					<< utility::pad_right( attachment_points, 10) << pose.residue( resnum ).carbohydrate_info()->short_name() << std::endl;

				carbohydrate_residues += 1;

			} else if ( pose.residue( resnum ).is_branch_point() ) {
				std::cout << "Branch Point: " << pose.residue( resnum ).name3()<<" "<< resnum << std::endl;
				protein_branches += 1;

			}
		}
		std::cout << "Glycan Residues: " << carbohydrate_residues <<std::endl;
		std::cout << "Protein BPs: " << protein_branches << std::endl;

	}

	virtual bool
	reinitialize_for_each_job() const {
		return true;
	}


private:  // Private methods
	// Set command-line options.  (Called by init())
	//void
	//set_commandline_options()
	//{
	// using namespace basic::options;

	//}


	// Initialize data members from arguments.
	void
	init()
	{

		type( "GlycanInfoMover" );



	}

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void
	copy_data( GlycanInfoMover & /*object_to_copy_to*/, GlycanInfoMover const & /*object_to_copy_from*/ )
	{

	}



private:  // Private data



};

typedef utility::pointer::shared_ptr< GlycanInfoMover > GlycanInfoMoverOP;

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::jd2;

		devel::init( argc, argv );
		//register_options();

		if ( ( ! option [ in::file::l ].user() ) && ( ! option [ in::file::s ].user() ) ) {
			utility_exit_with_message("Please specify either -s or -l to specify the input PDB.");
		}



		// Make sure the default JobOutputter is SilentJobOutputter to ensure that when this
		// is called with default arguments is prints a proper scorefile and not the hacky thing that
		// the  JobOutputter scorefile() function produces (which for example skips Evaluators!!)

		// Set up a job outputter that writes a scorefile and no PDBs and no Silent Files.
		// Copied from score_jd2.
		SilentFileJobOutputterOP jobout( new SilentFileJobOutputter );
		jobout->set_write_no_structures();
		jobout->set_write_separate_scorefile(true);


		protocols::jd2::JobDistributor::get_instance()->set_job_outputter( JobDistributorFactory::create_job_outputter( jobout ));


		GlycanInfoMoverOP mover_protocol( new GlycanInfoMover() );

		protocols::jd2::JobDistributor::get_instance()->go( mover_protocol );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
