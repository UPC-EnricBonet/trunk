#include <QCoreApplication>
#include <QTranslator>
#include <QStringList>
#include <QString>
#include <QDir>
#include <iostream>
using namespace std;

//#include "proostglobals.h"

extern int cheprooMain(QStringList& cheprooParams);


int main ( int argc, char *argv[] )
{
    QCoreApplication app(argc, argv);

    QTranslator translator;
    // For Automatic language selection, use: QLocale::system().name()
    QString translationFile = app.arguments()[0].section(QDir::separator() , 0, -2) +  QDir::separator() + "cheprooLib_es";
    bool translationLoaded = translator.load(translationFile);
    app.installTranslator(&translator);

    QStringList cheprooArgs = app.arguments();
    int i = cheprooMain(cheprooArgs);

    return(i);
}
