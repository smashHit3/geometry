#include "MainWindow.h"

#include <QApplication>

int main(int argc, char* argv[]) {
    QApplication application(argc, argv);
    gui::MainWindow window;
    window.show();
    return application.exec();
}