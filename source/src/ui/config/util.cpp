#include <ui/config/util.h>

#include <QSettings>

namespace ui {
namespace config {

/*UserCredentials get_user_credentials(QString const &host)
{
	QString prefix = "credentials/" + (host.isEmpty() ? "" : host + '_');

	QSettings settings;

	UserCredentials result;

	result.user     = settings.value(prefix + "user", "").toString();
	result.password = settings.value(prefix + "password", "").toString();

	return result;
}*/

static QString const _credentials_prefix_ = "credentials/";


UserCredentials get_user_credentials()
{
	QSettings settings;

	UserCredentials result;

	result.user     = settings.value(_credentials_prefix_ + "user", "").toString();
	result.password = settings.value(_credentials_prefix_ + "password", "").toString();

	return result;
}


void set_user_credentials(UserCredentials const &c)
{
	QSettings settings;

	settings.setValue(_credentials_prefix_ + "user", c.user);
	settings.setValue(_credentials_prefix_ + "password", c.password);
}


} // namespace config
} // namespace ui
