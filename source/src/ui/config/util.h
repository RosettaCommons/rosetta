#pragma once

#include <QString>

namespace ui {
namespace config {


/// structure to hold user login+password
struct UserCredentials
{
	QString user, password;
};

/// retrive user credentials stored in user level config file
UserCredentials get_user_credentials();


/// store user credentials into user level config file
void set_user_credentials(UserCredentials const &c);



/// get/set preferered external PDB file viewer
QString get_pdb_viewer_path();
void set_pdb_viewer_path(QString const &);



/// retrive user credentials for particular host from user level config file
//UserCredentials get_user_credentials(QString const &host);

} // namespace config
} // namespace ui
