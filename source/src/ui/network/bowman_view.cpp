#include "bowman_view.h"
#include "ui_bowman_view.h"

BowmanView::BowmanView(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::BowmanView)
{
    ui->setupUi(this);
}

BowmanView::~BowmanView()
{
    delete ui;
}
