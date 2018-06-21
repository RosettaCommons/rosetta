#ifndef BOWMAN_VIEW_H
#define BOWMAN_VIEW_H

#include <QWidget>

namespace Ui {
class BowmanView;
}

class BowmanView : public QWidget
{
    Q_OBJECT

public:
    explicit BowmanView(QWidget *parent = 0);
    ~BowmanView();

private:
    Ui::BowmanView *ui;
};

#endif // BOWMAN_VIEW_H
