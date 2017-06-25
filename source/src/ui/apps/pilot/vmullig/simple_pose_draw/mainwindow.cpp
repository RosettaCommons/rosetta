#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QTimer>
#include <QSurfaceFormat>
#include <QOpenGLFunctions>
#include <QtOpenGL>

#if defined(MAC) || defined(__APPLE__) || defined (__OSX__)
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <core/pose/Pose.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>
#include <numeric/xyzVector.hh>
#include <core/chemical/AtomType.hh>

MainWindow::MainWindow(QWidget */*parent*/) :
	ui(new Ui::MainWindow),
	context_(new QOpenGLContext),
	open_gl_functions_(nullptr), /*Initialized below*/
	pose_(new core::pose::Pose ),
	rotation_(0)
{
	setSurfaceType( QWindow::OpenGLSurface );
	QSurfaceFormat format;
	format.setProfile(QSurfaceFormat::CompatibilityProfile);
	format.setVersion(2, 1);
	setFormat(format);

	context_->setFormat(format);
	context_->create();
	context_->makeCurrent(this);

	create_pose();

	open_gl_functions_ = context_->functions();

	QTimer *timer(new QTimer(this));
	connect(timer, SIGNAL(timeout()), this, SLOT(updateAnimation()));
	timer->start(30);
}

MainWindow::~MainWindow()
{
    delete ui;
	delete context_;
}

void
MainWindow::updateAnimation() {
	rotation_ += 1;
	if(rotation_ > 360.0) rotation_ -= 360.0;
	pose_->set_chi(1, 10, rotation_*2.0 );
	update();
}


void
MainWindow::initializeGL() {
	glEnable(GL_DEPTH_TEST);
	resizeGL(this->width(), this->height());

	glShadeModel(GL_SMOOTH);

	glEnable(GL_LIGHT1);
	GLfloat lightAmbient1[] = { 0.1f, 0.15f, 0.2f, 1.0f };
	GLfloat lightDiffuse1[] = { 1.0f, 0.95f, 0.92f, 1.0f };
	GLfloat lightPosition1[] = { -200.0f, 200.0f, -55.0f, 1.0f};
	glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmbient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDiffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, lightDiffuse1);
	glLightfv(GL_LIGHT1, GL_POSITION, lightPosition1);

	glEnable(GL_LIGHT2);
	GLfloat lightAmbient2[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	GLfloat lightDiffuse2[] = { 0.3f, 0.45f, 0.5f, 1.0f };
	GLfloat lightPosition2[] = { 200.0f, -200.0f, -15.0f, 1.0f};
	glLightfv(GL_LIGHT2, GL_AMBIENT, lightAmbient2);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, lightDiffuse2);
	glLightfv(GL_LIGHT2, GL_SPECULAR, lightDiffuse2);
	glLightfv(GL_LIGHT2, GL_POSITION, lightPosition2);
}

void
MainWindow::resizeGL(
	int w,
	int h
) {
	glViewport(0,0,w,h);
	qreal aspectRatio(qreal(w)/qreal(h));
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1*aspectRatio, 1*aspectRatio, -1, 1, 1, -1);
	perspectiveGL(75, 1.0, 0.1, 40000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

}

void
MainWindow::paintGL() {
	resizeGL( width(), height() );
	glClearColor(0.2f,0.3f,0.4f,1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glEnable(GL_LIGHTING);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, -2.0f, -25.0f);
	glRotatef(rotation_, 0.0f, 1.0f, 0.0f);
	glRotatef(-25.0f, 1.0f, 0.0f, 0.0f);
	glRotatef(25.0f, 0.0f, 0.0f, 1.0f);

	GLfloat colorvect[] = { 0.25f, 0.25f, 0.3f, 1.0f};
	GLfloat colorvectN[] = { 0.15f, 0.15f, 0.75f, 1.0f};
	GLfloat colorvectO[] = { 0.75f, 0.15f, 0.15f, 1.0f};
	GLfloat colorvectS[] = { 0.85f, 0.82f, 0.07f, 1.0f};
	GLfloat colorvectH[] = { 0.75f, 0.75f, 0.75f, 1.0f};
	GLfloat colorvectspec[] = { 0.2f, 0.2f, 0.2f, 1.0f};
	glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvect );
	glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, colorvectspec );
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 15.0f );

	for(core::Size ir(1), irmax(pose_->total_residue()); ir<=irmax; ++ir) {
		for(core::Size ia(1), iamax(pose_->residue(ir).natoms()); ia<=iamax; ++ia) {
			numeric::xyzVector< core::Real > xyz( pose_->xyz( core::id::AtomID(ia, ir) ) );

			glTranslatef( xyz.x(), xyz.y(), xyz.z() );
			std::string const elem( pose_->residue(ir).atom_type(ia).element() );
			if( !elem.compare("N") ) {
				glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectN );
			} else if( !elem.compare("S") ) {
				glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectS );
			} else if( !elem.compare("O") ) {
				glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectO );
			} else if( pose_->residue(ir).atom_type(ia).is_polar_hydrogen() || !pose_->residue(ir).atom_name(ia).compare(" HG ") ) {
				glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectH );
			} else {
				glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvect );
			}
			glutSolidSphere( 0.75*pose_->residue(ir).atom_type(ia).lj_radius(), 32, 16 );
			glTranslatef( -xyz.x(), -xyz.y(), -xyz.z() );
		}
	}

	glDisable(GL_LIGHTING);

	glFlush();

}

void
MainWindow::paintEvent(
	QPaintEvent* //event
) {
	paintGL();
}

void
MainWindow::resizeEvent(
	QResizeEvent* event
) {
	event->accept();
	resizeGL(this->width(), this->height());
}

void
MainWindow::perspectiveGL( GLdouble fovY, GLdouble aspect, GLdouble zNear, GLdouble zFar ) const
{
	const GLdouble pi = 3.1415926535897932384626433832795;
	GLdouble fW, fH;

	//fH = tan( (fovY / 2) / 180 * pi ) * zNear;
	fH = tan( fovY / 360 * pi ) * zNear;
	fW = fH * aspect;

	glFrustum( -fW, fW, -fH, fH, zNear, zFar );
}

void
MainWindow::create_pose() const {
	core::chemical::ResidueTypeSetCOP standard_residues( pose_->residue_type_set_for_pose( core::chemical::FULL_ATOM_t ) );
	for(core::Size i(1); i<=20; ++i) {
		core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( standard_residues->name_map( i==10 ? "CYS" : "ALA") ) );
		if(i==1) {
			pose_->append_residue_by_jump( *new_rsd, 1 );
		} else {
			pose_->append_residue_by_bond( *new_rsd, true );
		}
	}
	core::pose::add_variant_type_to_pose_residue( *pose_, core::chemical::LOWER_TERMINUS_VARIANT, 1 );
	core::pose::add_variant_type_to_pose_residue( *pose_, core::chemical::UPPER_TERMINUS_VARIANT, 20 );
	for(core::Size i(1), imax(pose_->total_residue()); i<=imax; ++i) {
		pose_->set_phi(i, -57.8);
		pose_->set_psi(i, -47.0);
		pose_->set_omega(i, 180.0);
	}
	pose_->set_chi(1, 10, 0.0);
	pose_->set_chi(2, 10, 180.0);

	//Compute center of mass:
	numeric::xyzVector< core::Real > xyz;
	xyz.x(0); xyz.y(0); xyz.z(0);
	for(core::Size ir(1), irmax(pose_->total_residue()); ir<=irmax; ++ir) {
		xyz.x( xyz.x() + pose_->residue(ir).xyz("CA").x() );
		xyz.y( xyz.y() + pose_->residue(ir).xyz("CA").y() );
		xyz.z( xyz.z() + pose_->residue(ir).xyz("CA").z() );
	}
	xyz *= 1.0/static_cast<core::Real>(pose_->total_residue());

	//Centering the pose:
	for(core::Size ir(1), irmax(pose_->total_residue()); ir<=irmax; ++ir) {
		for(core::Size ia(1), iamax(pose_->residue(ir).natoms()); ia<=iamax; ++ia) {
			core::id::AtomID curat(ia,ir);
			pose_->set_xyz(curat, pose_->xyz(curat) - xyz);
		}
	}

	pose_->update_residue_neighbors();
}
