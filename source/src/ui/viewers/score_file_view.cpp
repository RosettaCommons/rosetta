#include <ui/viewers/score_file_view.h>
#include "ui_score_file_view.h"

#include <ui/task/task_view.h>
#include <ui/task/file.h>
#include <ui/task/task.h>

#include <QScatterSeries>
#include <QFileInfo>
#include <QDir>
#include <QDebug>

namespace ui {
namespace viewers {

QT_CHARTS_USE_NAMESPACE

std::pair<std::vector<QString>,  std::vector<ScorePoint> > parse_score_file(QByteArray const &file)
{
	/* QList<QByteArray> */ auto lines = file.split('\n');

	std::vector<QString> keys;
	std::vector< std::map<QString, QString> > result;

	while( not lines.isEmpty() ) {
		QString line = lines.takeFirst();
		if( line.startsWith("SCORE: ") ) {
			QStringList fields = line.split(' ', QString::SkipEmptyParts).mid(1);
			//qDebug() << fields;

			if( keys.empty() ) {
				while( not fields.isEmpty() ) keys.push_back( fields.takeFirst() );
			} else {
				std::map<QString, QString> r;

				for(uint i = 0; not fields.isEmpty()  and  i < keys.size(); ++i) {
					r[ keys[i] ] = fields.takeFirst();
				}
				result.push_back( std::move(r) );
			}
		}
	}

	return std::make_pair(move(keys), move(result));
}

auto static QPointF_Less = [](QPointF const &lhs, QPointF const &rhs) -> bool {
	if( lhs.x() == rhs.x() ) return lhs.y() < rhs.y();
		else return lhs.x() < rhs.x();
};


ScoreFileView::ScoreFileView(std::pair<QString const, task::FileSP> const & score_file, task::TaskSP const &task,  QWidget *parent) :
    QWidget(parent),
	task_(task),
	score_file_(score_file),
	decoys_(QPointF_Less),
    ui(new Ui::ScoreFileView)
{
    ui->setupUi(this);

	QFileInfo fi(score_file_.first);
	decoys_path_ = fi.dir().path();

	update_ui_from_file_data();
}


void ScoreFileView::update_ui_from_file_data()
{
	auto chart = ui->chart_view->chart();
	chart->removeAllSeries();
	decoys_.clear();

	/* std::pair<std::vector<QString>,  std::vector<ScorePoint> > */ auto keys_sf = parse_score_file( score_file_.second->data() );

	score_terms_ = keys_sf.first;
	auto & sf = keys_sf.second;

    QScatterSeries *plot = new QScatterSeries();

    ui->chart_view->chart()->setTitle( QString("Score Plot for %1").arg(score_file_.first) );

    plot->setName(x_axis_ + " / " + y_axis_);
	plot->setMarkerSize(8.0);
	plot->setColor(Qt::blue);

	for( auto & point : sf ) {
		bool ok;
		double x = point[x_axis_].toDouble(&ok); if( not ok ) continue;
		double y = point[y_axis_].toDouble(&ok); if( not ok ) continue;

		QPointF p(x, y);
		*plot << p;

		decoys_.insert( std::make_pair(p, point) );
		//

		//qDebug() << "point:" << point;
		//qDebug() << "adding:" << p;
	}

	ui->chart_view->setRenderHint(QPainter::Antialiasing);
    ui->chart_view->chart()->addSeries(plot);
    //ui->chart_view->chart()->addSeries(m_scatter);
    ui->chart_view->chart()->createDefaultAxes();
    //ui->chart_view->chart()->axisX()->setRange(0, 4.5);
    //ui->chart_view->chart()->axisY()->setRange(0, 4.5);

	ui->chart_view->setRubberBand( QChartView::RectangleRubberBand );
	//ui->chart_view->chart()->zoom(0.5);

	connect(plot, &QScatterSeries::clicked, this, &ScoreFileView::point_clicked);
	//connect(plot, &QScatterSeries::doubleClicked, this, &ScoreFileView::point_double_clicked);

	connect(plot, &QScatterSeries::hovered, this, &ScoreFileView::hovered);


	for(auto const & st : score_terms_) {
		ui->x_axis->addItem(st);
		ui->y_axis->addItem(st);
	}
	ui->x_axis->setCurrentText(x_axis_);
	ui->y_axis->setCurrentText(y_axis_);
}


ScoreFileView::~ScoreFileView()
{
    delete ui;
}

void ScoreFileView::point_clicked(QPointF const & point)
{
	auto it = decoys_.find(point);
	if( it != decoys_.end() ) {
		QString description = it->second["description"];
		if( not description.isEmpty() ) {
			QString file_name = decoys_path_ + "/" + description + ".pdb";

			//file_name = "input/input.pdb";
			qDebug() << "ChartView::point_clicked:" << point << "->" << file_name;

			auto & files = task_->files();

			auto it = files.find(file_name);
			if( it == files.end() ) {
				qDebug() << "ChartView::point_clicked: could not find file:" << file_name;
			} else {
				open_file_viewer(*it, task_, this);
			}
		}
	}
}

void ScoreFileView::hovered(const QPointF &point, bool state)
{
	//qDebug() << "ScoreFileView::hovered:" << point << " " << state;
	if(state) {
		auto it = decoys_.find(point);
		if( it != decoys_.end() ) {
			QString text = QString("name: %1").arg(it->second["description"]);
			for(auto const & field : it->second)
				if( field.first != "description" ) text += QString("\n%1: %2").arg(field.first, field.second);
			setToolTip(text);
			//setToolTip(QString("%1 , %2").arg(point.x()).arg(point.y()));
		}
	} else setToolTip("");
}

void ScoreFileView::on_x_axis_activated(QString const &text)
{
	x_axis_ = text;
	update_ui_from_file_data();
}

void ScoreFileView::on_y_axis_activated(QString const &text)
{
	y_axis_ = text;
	update_ui_from_file_data();
}


// void ScoreFileView::point_double_clicked(QPointF const & point)
// {
// 	qDebug() << "ChartView::point_double_clicked:" << point;
// }


} // namespace viewers
} // namespace ui
