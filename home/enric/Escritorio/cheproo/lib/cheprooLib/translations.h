#ifndef TRANSLATIONS_H
#define TRANSLATIONS_H

#include "base.h"

#define PrTr(category, message) qPrintable(QCoreApplication::translate(category, message))

#endif