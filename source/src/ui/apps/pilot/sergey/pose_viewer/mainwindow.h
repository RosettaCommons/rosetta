#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <ui/network/bowman.h>


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private Q_SLOTS:
	void client_connected(quint64);
	void client_disconnected(quint64);

private:
	ui::network::BowmanThread bowman;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
