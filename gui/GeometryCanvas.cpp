#include "GeometryCanvas.h"

#include <QGraphicsEllipseItem>
#include <QGraphicsRectItem>
#include <QPainter>

namespace {

constexpr qreal shapeSize = 100.0;
constexpr qreal shapeOffset = 24.0;

QPen shapePen() {
    return {QColor("#1d4ed8"), 2.0};
}

QBrush shapeBrush() {
    return {QColor(96, 165, 250, 90)};
}

} // namespace

namespace gui {

GeometryCanvas::GeometryCanvas(QWidget* parent)
    : QGraphicsView(parent) {
    setScene(&m_scene);
    setRenderHint(QPainter::Antialiasing);
    setDragMode(QGraphicsView::RubberBandDrag);
    setSceneRect(-400.0, -300.0, 800.0, 600.0);
    setBackgroundBrush(QColor("#f8fafc"));
}

void GeometryCanvas::addRectangle() {
    addShape(m_scene.addRect(nextShapeBounds(), shapePen(), shapeBrush()));
}

void GeometryCanvas::addEllipse() {
    addShape(m_scene.addEllipse(nextShapeBounds(), shapePen(), shapeBrush()));
}

void GeometryCanvas::removeSelected() {
    const auto selectedItems = m_scene.selectedItems();
    for (QGraphicsItem* item : selectedItems) {
        m_scene.removeItem(item);
        delete item;
    }
}

void GeometryCanvas::clearShapes() {
    m_scene.clear();
    m_nextShapeOffset = 0;
}

QRectF GeometryCanvas::nextShapeBounds() {
    const qreal offset = static_cast<qreal>(m_nextShapeOffset) * shapeOffset;
    m_nextShapeOffset = (m_nextShapeOffset + 1) % 8;
    return {-shapeSize / 2.0 + offset, -shapeSize / 2.0 + offset, shapeSize, shapeSize};
}

void GeometryCanvas::addShape(QGraphicsItem* item) {
    item->setFlag(QGraphicsItem::ItemIsMovable);
    item->setFlag(QGraphicsItem::ItemIsSelectable);
}

} // namespace gui
