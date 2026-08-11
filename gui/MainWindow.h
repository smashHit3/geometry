#ifndef GUI_MAIN_WINDOW_H
#define GUI_MAIN_WINDOW_H

#include <QKeySequence>
#include <QMainWindow>
#include <QString>

class QAction;

namespace gui {

class GeometryCanvas;

class MainWindow : public QMainWindow {
public:
    explicit MainWindow(QWidget* parent = nullptr);

private:
    GeometryCanvas* m_canvas;

    QAction* createAction(const QString& text, const QKeySequence& shortcut);
};

} // namespace gui

#endif // GUI_MAIN_WINDOW_H
