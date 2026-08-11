#include "MainWindow.h"

#include "GeometryCanvas.h"

#include <QAction>
#include <QKeySequence>
#include <QStatusBar>
#include <QToolBar>

namespace gui {

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent),
      m_canvas(new GeometryCanvas(this)) {
    setWindowTitle(tr("Geometry Visualizer"));
    resize(960, 700);
    setCentralWidget(m_canvas);

    QToolBar* toolbar = addToolBar(tr("Shapes"));
    toolbar->setMovable(false);

    QAction* addRectangleAction = createAction(tr("Add Rectangle"), QKeySequence(Qt::CTRL | Qt::Key_R));
    connect(addRectangleAction, &QAction::triggered, m_canvas, &GeometryCanvas::addRectangle);
    toolbar->addAction(addRectangleAction);

    QAction* addEllipseAction = createAction(tr("Add Circle"), QKeySequence(Qt::CTRL | Qt::Key_E));
    connect(addEllipseAction, &QAction::triggered, m_canvas, &GeometryCanvas::addEllipse);
    toolbar->addAction(addEllipseAction);

    toolbar->addSeparator();

    QAction* deleteAction = createAction(tr("Delete Selected"), QKeySequence::Delete);
    connect(deleteAction, &QAction::triggered, m_canvas, &GeometryCanvas::removeSelected);
    toolbar->addAction(deleteAction);

    QAction* clearAction = createAction(tr("Clear Canvas"), QKeySequence(Qt::CTRL | Qt::Key_L));
    connect(clearAction, &QAction::triggered, m_canvas, &GeometryCanvas::clearShapes);
    toolbar->addAction(clearAction);

    statusBar()->showMessage(tr("Add a shape, then select and drag it to move."));
}

QAction* MainWindow::createAction(const QString& text, const QKeySequence& shortcut) {
    QAction* action = new QAction(text, this);
    action->setShortcut(shortcut);
    addAction(action);
    return action;
}

} // namespace gui
