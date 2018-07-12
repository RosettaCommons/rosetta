#ifndef FUNCTION_SETUP_DIALOG_H
#define FUNCTION_SETUP_DIALOG_H

#include <ui/network/argument.fwd.h>

#include <QDialog>

namespace Ui {
class FunctionSetupDialog;
}


namespace ui {
namespace network {

class FunctionSetupDialog : public QDialog
{
    Q_OBJECT

public:
    explicit FunctionSetupDialog(QString const & function_name, ArgumentsSP const &, QWidget *parent = 0);
    ~FunctionSetupDialog();

private Q_SLOTS:
	void ok_clicked();


private:
    Ui::FunctionSetupDialog *ui;
};

} // namespace network
} // namespace ui


#endif // FUNCTION_SETUP_DIALOG_H
