#ifndef CONFIG_DIALOG_H
#define CONFIG_DIALOG_H

#include <QDialog>

namespace Ui {
class ConfigDialog;
}

namespace ui {
namespace config {

class ConfigDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ConfigDialog(QWidget *parent = 0);
    ~ConfigDialog();


private Q_SLOTS:
	void update_ui_from_user_settings();
	void update_user_settings_from_ui();

private:
    Ui::ConfigDialog *ui;
};


} // namespace config
} // namespace ui

#endif // CONFIG_DIALOG_H
