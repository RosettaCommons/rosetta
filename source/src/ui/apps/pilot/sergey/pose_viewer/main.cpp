#include <ui/widgets/pose_editor.h>
#include <QApplication>

#include <protocols/init/init.hh>
#include <libxml/parser.h> // dummy include to make sure external excludes path is set
#include <ui/ui_lib_test.h>
#if defined(MAC) || defined(__APPLE__) || defined (__OSX__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#include <QDebug>

int main(int argc, char *argv[])
{
	qRegisterMetaType<JSON_SP> ("JSON_SP");
	qRegisterMetaType<JSON_CSP>("JSON_CSP");

	qRegisterMetaType<std::string>("std::string");

	//qRegisterMetaType<ui::network::FunctionID>("ui::network::FunctionID");


    qDebug() << "test::main: Calling Rosetta init()...";

    glutInit(&argc, argv);
    protocols::init::init(argc, argv);
    ui::ui_lib_test::test();

    QApplication a(argc, argv);
	ui::widgets::PoseEditor w;
    w.show();

    return a.exec();
}
