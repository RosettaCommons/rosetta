#include "mainwindow.h"
#include <QApplication>

#include <core/init/init.hh>
//#include <protocols/init/init.hh>
#include <libxml/parser.h> // dummy include to make sure external excludes path is set

#include <QDebug>
#if defined(MAC) || defined(__APPLE__) || defined (__OSX__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

int main(int argc, char *argv[])
{
    qDebug() << "simple_pose_draw::main: Calling Rosetta init()...";

    core::init::init(argc, argv);

	glutInit(&argc, argv);

    QApplication a(argc, argv);
    MainWindow w;
	w.resize(600,800);
	w.setTitle(QString("Pose drawing demo"));
    w.show();

    return a.exec();
}
