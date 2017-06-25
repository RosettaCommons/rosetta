#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QOpenGLWindow>

#include <core/pose/Pose.fwd.hh>

namespace Ui {
class MainWindow;
}

class MainWindow : public QOpenGLWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

public Q_SLOTS:

	void updateAnimation();

protected:

	void initializeGL() override;
	void resizeGL(int w, int h) override;
	void paintGL() override;
	void paintEvent(QPaintEvent* event) override;
	void resizeEvent(QResizeEvent* event) override;

private:

	void perspectiveGL( GLdouble fovY, GLdouble aspect, GLdouble zNear, GLdouble zFar ) const;

	void create_pose() const;

private:
    Ui::MainWindow *ui;

	QOpenGLContext* context_;

	QOpenGLFunctions* open_gl_functions_;

	core::pose::PoseOP pose_;

	float rotation_;
};

#endif // MAINWINDOW_H
