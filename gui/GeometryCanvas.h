#ifndef GUI_GEOMETRY_CANVAS_H
#define GUI_GEOMETRY_CANVAS_H

#include <QGraphicsScene>
#include <QGraphicsView>
#include <QRectF>

class QGraphicsItem;

namespace gui {

class GeometryCanvas : public QGraphicsView {
public:
    explicit GeometryCanvas(QWidget* parent = nullptr);

    void addRectangle();
    void addEllipse();
    void removeSelected();
    void clearShapes();

private:
    QGraphicsScene m_scene;
    int m_nextShapeOffset = 0;

    QRectF nextShapeBounds();
    void addShape(QGraphicsItem* item);
};

} // namespace gui

#endif // GUI_GEOMETRY_CANVAS_H
