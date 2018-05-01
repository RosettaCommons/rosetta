#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <ui/network/bowman.h>

#include <QDebug>

using namespace ui::network;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

	bowman.start();
	connect(&bowman, &BowmanThread::client_connected,    this, &MainWindow::client_connected);
	connect(&bowman, &BowmanThread::client_disconnected, this, &MainWindow::client_disconnected);

	//setAttribute(Qt::WA_DeleteOnClose);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::client_connected(quint64 id)
{
	qDebug() << "MainWindow::client_connected(): " << id;
}

void MainWindow::client_disconnected(quint64 id)
{
	qDebug() << "MainWindow::client_disconnected(): " << id;
}
