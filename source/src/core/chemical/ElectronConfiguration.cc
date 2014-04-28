// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
// This file is based off of code from the Biochemistry Library (BCL).
// The BCL is copyright Vanderbilt University (Meiler Lab), a RosettaCommons member

/// @file   core/chemical/gasteiger/ElectronConfiguration.hh
/// @brief  The data for the BCL electron configuration
/// @author To Rosetta transitioning: Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/ElectronConfiguration.hh>

#include <core/chemical/gasteiger/util.hh>

#include <numeric/util.hh>
#include <utility/string_util.hh>
#include <utility/tools/make_vector.hh>

#include <string>

namespace core {
namespace chemical {

using utility::tools::make_vector;
using core::Size;


std::vector<std::string> const & ElectronConfiguration::PrincipalQuantumNumber_strings()
{
	static std::vector<std::string> PrincipalQuantumNumber_strings_ =
		make_vector<std::string> (
			"1",
			"2",
			"3",
			"4",
			"5",
			"6",
			"7",
			"PrincipalQuantumNumber" //GetStaticClassName< PrincipalQuantumNumber>()
											   );
	return PrincipalQuantumNumber_strings_;
}

std::vector<std::string> const & ElectronConfiguration::AngularMomentumQuantumNumber_strings()
{
	static std::vector<std::string>	AngularMomentumQuantumNumber_strings_ =
		make_vector<std::string> (
			"s",
			"p",
			"d",
			"f",
			"AngularMomentumQuantumNumber" //GetStaticClassName< AngularMomentumQuantumNumber>()
												  );
	return AngularMomentumQuantumNumber_strings_;
}

    //! @brief PrincipalQuantumNumber as string
    //! @param NUM the PrincipalQuantumNumber desired
    //! @return the PrincipalQuantumNumber as string
    const std::string &ElectronConfiguration::get_descriptor( const PrincipalQuantumNumber &NUM)
    {
		return PrincipalQuantumNumber_strings()[ NUM ];
    }

    //! @brief PrincipalQuantumNumber from string
    //! @return NUM the PrincipalQuantumNumber desired
    //! @param the PrincipalQuantumNumber as string
    ElectronConfiguration::PrincipalQuantumNumber ElectronConfiguration::get_principal_quantum_number( std::string const &STR ) {
    	for( int ii = 0; ii < MaxPrincipleQuantumNumber; ++ii) {
    		if( PrincipalQuantumNumber_strings()[ ii ] == STR ) return PrincipalQuantumNumber(ii);
    	}
    	return MaxPrincipleQuantumNumber;
    }

    //! @brief AngularMomentumQuantumNumber as string
    //! @param NUM the AngularMomentumQuantumNumber desired
    //! @return the AngularMomentumQuantumNumber as string
    const std::string &ElectronConfiguration::get_descriptor( const AngularMomentumQuantumNumber &NUM)
    {
		return AngularMomentumQuantumNumber_strings()[ NUM];
    }

    //! @brief AngularMomentumQuantumNumber as string
    //! @return NUM the AngularMomentumQuantumNumber desired
    //! @param the AngularMomentumQuantumNumber as string
    ElectronConfiguration::AngularMomentumQuantumNumber ElectronConfiguration::get_angular_momentum_quantum_number( std::string const &STR ) {
    	for( int ii = 0; ii < MaxAngularMomentumQuantumNumber; ++ii) {
    		if( AngularMomentumQuantumNumber_strings()[ ii ] == STR ) return AngularMomentumQuantumNumber(ii);
    	}
    	return MaxAngularMomentumQuantumNumber;
    }

  //////////
  // data //
  //////////

//#    //! single instance of that class
//#    const util::SiPtr< const util::ObjectInterface> ElectronConfiguration::s_Instance
//#    (
//#      GetObjectInstances().AddInstance( new ElectronConfiguration())
//#    );

    // core::Size [ 7][ 4];
	std::vector< std::vector<Size> > const & ElectronConfiguration::s_MaxElectronsInOrbital()
    {
		static std::vector< std::vector<Size> > s_MaxElectronsInOrbital_ =
			make_vector< std::vector<Size> > (
				make_vector< Size > ( 2, 0, 0, 0),
				make_vector< Size > ( 2, 6, 0, 0),
				make_vector< Size > ( 2, 6, 0, 0),
				make_vector< Size > ( 2, 6, 10, 0),
				make_vector< Size > ( 2, 6, 10, 0),
				make_vector< Size > ( 2, 6, 10, 14),
				make_vector< Size > ( 2, 6, 10, 14)
				);

		return s_MaxElectronsInOrbital_;
	}

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    ElectronConfiguration::ElectronConfiguration() :
      valence_electrons_sp_( numeric::get_undefined_size() ),
      valence_electrons_spd_( numeric::get_undefined_size() ),
      valence_quantum_number_( numeric::get_undefined_size() )
    {
      for
      (
        core::Size shell_type( 0), num_shells( MaxPrincipleQuantumNumber);
        shell_type < num_shells;
        ++shell_type
      )
      {
        for
        (
          core::Size orbital_type( 0), num_orbitals( MaxAngularMomentumQuantumNumber);
          orbital_type < num_orbitals;
          ++orbital_type
        )
        {
          electrons_[ shell_type][ orbital_type] = numeric::get_undefined_size();
        }
      }
    }

    //! construct from actual number of electrons
    ElectronConfiguration::ElectronConfiguration
    (
      const core::Size VALENCE_ELECTRONS_SP,
      const core::Size VALENCE_ELECTRONS_SPD,
      const core::Size NUMBER_ELECTRONS_1S,
      const core::Size NUMBER_ELECTRONS_1P,
      const core::Size NUMBER_ELECTRONS_1D,
      const core::Size NUMBER_ELECTRONS_1F,
      const core::Size NUMBER_ELECTRONS_2S,
      const core::Size NUMBER_ELECTRONS_2P,
      const core::Size NUMBER_ELECTRONS_2D,
      const core::Size NUMBER_ELECTRONS_2F,
      const core::Size NUMBER_ELECTRONS_3S,
      const core::Size NUMBER_ELECTRONS_3P,
      const core::Size NUMBER_ELECTRONS_3D,
      const core::Size NUMBER_ELECTRONS_3F,
      const core::Size NUMBER_ELECTRONS_4S,
      const core::Size NUMBER_ELECTRONS_4P,
      const core::Size NUMBER_ELECTRONS_4D,
      const core::Size NUMBER_ELECTRONS_4F,
      const core::Size NUMBER_ELECTRONS_5S,
      const core::Size NUMBER_ELECTRONS_5P,
      const core::Size NUMBER_ELECTRONS_5D,
      const core::Size NUMBER_ELECTRONS_5F,
      const core::Size NUMBER_ELECTRONS_6S,
      const core::Size NUMBER_ELECTRONS_6P,
      const core::Size NUMBER_ELECTRONS_6D,
      const core::Size NUMBER_ELECTRONS_6F,
      const core::Size NUMBER_ELECTRONS_7S,
      const core::Size NUMBER_ELECTRONS_7P,
      const core::Size NUMBER_ELECTRONS_7D,
      const core::Size NUMBER_ELECTRONS_7F
    ) :
      valence_electrons_sp_( VALENCE_ELECTRONS_SP),
      valence_electrons_spd_( VALENCE_ELECTRONS_SPD),
      valence_quantum_number_( 0)
    {
      electrons_[ e_1][ e_S] = NUMBER_ELECTRONS_1S;
      electrons_[ e_1][ e_P] = NUMBER_ELECTRONS_1P;
      electrons_[ e_1][ e_D] = NUMBER_ELECTRONS_1D;
      electrons_[ e_1][ e_F] = NUMBER_ELECTRONS_1F;
      electrons_[ e_2][ e_S] = NUMBER_ELECTRONS_2S;
      electrons_[ e_2][ e_P] = NUMBER_ELECTRONS_2P;
      electrons_[ e_2][ e_D] = NUMBER_ELECTRONS_2D;
      electrons_[ e_2][ e_F] = NUMBER_ELECTRONS_2F;
      electrons_[ e_3][ e_S] = NUMBER_ELECTRONS_3S;
      electrons_[ e_3][ e_P] = NUMBER_ELECTRONS_3P;
      electrons_[ e_3][ e_D] = NUMBER_ELECTRONS_3D;
      electrons_[ e_3][ e_F] = NUMBER_ELECTRONS_3F;
      electrons_[ e_4][ e_S] = NUMBER_ELECTRONS_4S;
      electrons_[ e_4][ e_P] = NUMBER_ELECTRONS_4P;
      electrons_[ e_4][ e_D] = NUMBER_ELECTRONS_4D;
      electrons_[ e_4][ e_F] = NUMBER_ELECTRONS_4F;
      electrons_[ e_5][ e_S] = NUMBER_ELECTRONS_5S;
      electrons_[ e_5][ e_P] = NUMBER_ELECTRONS_5P;
      electrons_[ e_5][ e_D] = NUMBER_ELECTRONS_5D;
      electrons_[ e_5][ e_F] = NUMBER_ELECTRONS_5F;
      electrons_[ e_6][ e_S] = NUMBER_ELECTRONS_6S;
      electrons_[ e_6][ e_P] = NUMBER_ELECTRONS_6P;
      electrons_[ e_6][ e_D] = NUMBER_ELECTRONS_6D;
      electrons_[ e_6][ e_F] = NUMBER_ELECTRONS_6F;
      electrons_[ e_7][ e_S] = NUMBER_ELECTRONS_7S;
      electrons_[ e_7][ e_P] = NUMBER_ELECTRONS_7P;
      electrons_[ e_7][ e_D] = NUMBER_ELECTRONS_7D;
      electrons_[ e_7][ e_F] = NUMBER_ELECTRONS_7F;

      for( core::Size valence_number = 1, max_shells = core::Size( e_7) + 1; valence_number < max_shells; valence_number++)
      {
        if( electrons_[ valence_number][ e_S] > 0)
        {
          valence_quantum_number_ = ElectronConfiguration::PrincipalQuantumNumber( valence_number);
        }
      }
    }

//#    //! @brief virtual copy constructor
//#    ElectronConfiguration *ElectronConfiguration::Clone() const
//#    {
//#      return new ElectronConfiguration( *this);
//#    }

  /////////////////
  // data access //
  /////////////////

//#    //! @brief returns class name
//#    //! @return the class name as const ref std::string
//#    const std::string &ElectronConfiguration::GetClassIdentifier() const
//#    {
//#      return GetStaticClassName( *this);
//#    }

    //! @return number valence_electrons_sp
    core::Size ElectronConfiguration::valence_electrons_sp() const
    {
      return valence_electrons_sp_;
    }

    //! @return number valence_electrons_sp
    core::Size ElectronConfiguration::valence_electrons_spd() const
    {
      return valence_electrons_spd_;
    }

//! @return the maximum number of electrons in SP orbitals for the noble gas in this period
core::Size ElectronConfiguration::max_valence_electrons_sp() const
{
	return s_MaxElectronsInOrbital()[ valence_quantum_number_][ e_S] + s_MaxElectronsInOrbital()[ valence_quantum_number_][ e_P];
}

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ElectronConfiguration::read( std::istream &ISTREAM)
    {
      std::string tag;
      ISTREAM >> tag;
      if ( tag != "ElectronConfiguration:" ) {
    	  std::string line;
    	  getline( ISTREAM, line );
    	  std::cout << "$$$" << line << std::endl;
    	  getline( ISTREAM, line );
    	  std::cout << "$$$" << line << std::endl;
    	  utility_exit_with_message( "Malformatted Bcl elements file - 'ElectronConfiguration:' tag expected '" + tag + "' found." );
      }
      // read member
      gasteiger::safe_read( ISTREAM, valence_electrons_sp_); //io::Serialize::read( m_valence_electrons_sp, ISTREAM);
      gasteiger::safe_read( ISTREAM, valence_electrons_spd_); //io::Serialize::read( valence_electrons_spd_, ISTREAM);
      gasteiger::safe_read( ISTREAM, valence_quantum_number_); //io::Serialize::read( valence_quantum_number_, ISTREAM);

      // reset each state
      for
      (
        core::Size shell_type( 0), num_shells( MaxPrincipleQuantumNumber);
        shell_type < num_shells;
        ++shell_type
      )
      {
        for
        (
          core::Size orbital_type( 0), num_orbitals( MaxAngularMomentumQuantumNumber);
          orbital_type < num_orbitals;
          ++orbital_type
        )
        {
          core::Size nelec;
          gasteiger::safe_read( ISTREAM, nelec);
          electrons_[ shell_type][ orbital_type] = nelec;
        }
      }
//#
//#      // get how many occupied states there were
//#      core::Size num_occupied_states( 0);
//#      ISTREAM >> num_occupied_states; //io::Serialize::read( num_occupied_states, ISTREAM);
//#      for( core::Size curr_state( 0); curr_state < num_occupied_states; ++curr_state)
//#      {
//#        std::string state;
//#        ISTREAM >> state;
//#
//#        PrincipalQuantumNumber principle_n( get_principal_quantum_number(state.substr( 0, 1)) );
//#        AngularMomentumQuantumNumber angular_n( get_angular_momentum_quantum_number( state.substr( 1, 1)) );
//#
//#        core::Size read_electrons( utility::string2int( state.substr( 2) ) );
//#
//#        electrons_[ principle_n][ angular_n] = read_electrons;
//#      }

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &ElectronConfiguration::write( std::ostream &OSTREAM) const
    {
      OSTREAM << "    ElectronConfiguration: ";
      // write member
      gasteiger::safe_write( OSTREAM, valence_electrons_sp_); //io::Serialize::Write( m_valence_electrons_sp, OSTREAM, INDENT);
      gasteiger::safe_write( OSTREAM, valence_electrons_spd_); //io::Serialize::Write( valence_electrons_spd_, OSTREAM, INDENT);
      gasteiger::safe_write( OSTREAM, valence_quantum_number_); //io::Serialize::Write( valence_quantum_number_, OSTREAM, INDENT);

//#      // count the number of occupied states
//#      core::Size num_occupied_states( 0);
//#      for
//#      (
//#        core::Size shell_type( 0), num_shells( MaxPrincipleQuantumNumber);
//#        shell_type < num_shells;
//#        ++shell_type
//#      )
//#      {
//#        for
//#        (
//#          core::Size orbital_type( 0), num_orbitals( MaxAngularMomentumQuantumNumber);
//#          orbital_type < num_orbitals;
//#          ++orbital_type
//#        )
//#        {
//#          if( electrons_[ shell_type][ orbital_type] > 0)
//#          {
//#            ++num_occupied_states;
//#          }
//#        }
//#      }
//#
//#      // write how many occupied states there are, followed by each occupied state
//#      OSTREAM << num_occupied_states << ' '; //io::Serialize::Write( num_occupied_states, OSTREAM, INDENT) << ' ';
      for
      (
        core::Size shell_type( 0), num_shells( MaxPrincipleQuantumNumber);
        shell_type < num_shells;
        ++shell_type
      )
      {
        for
        (
          core::Size orbital_type( 0), num_orbitals( MaxAngularMomentumQuantumNumber);
          orbital_type < num_orbitals;
          ++orbital_type
        )
        {
//#          if( electrons_[ shell_type][ orbital_type] > 0)
//#          {
//#            OSTREAM << get_descriptor( PrincipalQuantumNumber( shell_type ) ) // PrincipalQuantumNumberEnum( PrincipalQuantumNumber( shell_type))
//#                    << get_descriptor( AngularMomentumQuantumNumber( orbital_type) ) // AngularMomentumQuantumNumberEnum( AngularMomentumQuantumNumber( orbital_type))
//#          }
        	gasteiger::safe_write( OSTREAM, electrons_[ shell_type][ orbital_type] );
        }
      }

      OSTREAM << std::endl;

      // return
      return OSTREAM;
    }


} // namespace core
} // namespace chemical
