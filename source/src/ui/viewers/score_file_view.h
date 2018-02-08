#ifndef SCORE_FILE_VIEW_H
#define SCORE_FILE_VIEW_H

#include <ui/task/file.fwd.h>
#include <ui/task/task.fwd.h>

#include <QWidget>

#include <map>
#include <functional>

namespace Ui {
class ScoreFileView;
}

namespace ui {
namespace viewers {

using ScorePoint = std::map<QString, QString>;

class ScoreFileView : public QWidget
{
    Q_OBJECT

	using Decoys = std::map<QPointF, ScorePoint, std::function< bool(QPointF const &lhs, QPointF const &rhs) > >;

public:
    explicit ScoreFileView(std::pair<QString const, task::FileSP> const & score_file, task::TaskSP const &task, QWidget *parent = 0);
    ~ScoreFileView();

private Q_SLOTS:
	void update_ui_from_file_data();


	void point_clicked(QPointF const & point);
	//void point_double_clicked(QPointF const & point);
	void hovered(const QPointF &point, bool state);

	void on_x_axis_activated(QString const &text);
	void on_y_axis_activated(QString const &text);


private:
	task::TaskSP task_;

	std::pair<QString const, task::FileSP> score_file_;

	Decoys decoys_;
	QString decoys_path_;

	QString x_axis_="rms", y_axis_="total_score";

	std::vector<QString> score_terms_;

    Ui::ScoreFileView *ui;
};


} // namespace viewers
} // namespace ui

#endif // SCORE_FILE_VIEW_H
