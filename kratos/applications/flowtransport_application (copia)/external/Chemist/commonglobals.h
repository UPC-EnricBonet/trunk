#ifndef COMMONGLOBALS_H
#define COMMONGLOBALS_H

#include <QtGlobal>
#include <QString>
#include <QDateTime>
#include <valarray>
#include <math.h>
#include <map>


#include <vector>
using namespace std;

#include "generalexception.h"

#define CHEPROOST_VERSION "0.01"

# define PROOST_EXEC "proost"
# define UPDATER_EXEC "ProostXMLUpdater"

static const char* INTERNAL_ERROR = QT_TRANSLATE_NOOP("prGlobals", "Proost internal error");
static const char* ERROR_SUB_NOT_IMPLEMENTED_YET = QT_TRANSLATE_NOOP("prGlobals", "ToDo Using a sub which has not been implemented yet");
static const char* ERROR_BADMETHODPAR = QT_TRANSLATE_NOOP("prGlobals", "Error: the input parameters of a method had unacceptable values");
static const char* ERROR_IN_READ =  QT_TRANSLATE_NOOP("prGlobals", "Error during reading of input file");
static const char* ERROR_BADMETHOD =  QT_TRANSLATE_NOOP("prGlobals", "Error: the method that was called is not implemented by the object that received the method call");
static const char* ERROR_FIND_ATTR =  QT_TRANSLATE_NOOP("prGlobals", "An obligatory xml attribute was missing");
static const char* ERROR_FIND_ELEM =  QT_TRANSLATE_NOOP("prGlobals", "An obligatory xml element was missing");
static const char* ERROR_NEEDS_UPDATE =  QT_TRANSLATE_NOOP("prGlobals", "The method where this error was generated needs to be updated and cannot be used currently");
static const char* ERROR_NOT_FINISHED =  QT_TRANSLATE_NOOP("prGlobals", "The method where this error was generated is not finished");
static const char* ERROR_TYPECAST =  QT_TRANSLATE_NOOP("prGlobals", "An object of the wrong type was provided");
static const char* ERROR_NOT_DOUBLE =  QT_TRANSLATE_NOOP("prGlobals", "Error: expected a double, but did not find it.");
static const char* ERROR_NOT_INTEGER =  QT_TRANSLATE_NOOP("prGlobals", "Error: expected an integer, but did not find it.");
static const QString ID_NOT_NEEDED = "ID_";




template<class T>
T FindInMap(std::map<QString, T>* mapToSearch, QString ObjectToFind, QString errorMessage)
{
    typename std::map<QString, T>::iterator it; 
    it = mapToSearch->find(ObjectToFind);
	if ( it == mapToSearch->end() )
    {
        throw GeneralException(errorMessage);
    }
    return(it->second);
}

template<class T>
T FindInMap(std::map<QString, T>* mapToSearch, QString ObjectToFind)
{
    typename std::map<QString, T>::iterator it; 
    it = mapToSearch->find(ObjectToFind);
	if ( it == mapToSearch->end() )
    {
        return 0;
    }
    else
    {
        return(it->second);
    }
}

template<class T>
T FindInMap(std::multimap<QString, T>* mapToSearch, QString ObjectToFind, QString errorMessage)
{
    typename std::multimap<QString, T>::iterator it; 
    it = mapToSearch->find(ObjectToFind);
	if ( it == mapToSearch->end() )
    {
        throw GeneralException(errorMessage);
    }
    return(it->second);
}

template<class T>
T FindInMap(std::multimap<QString, T>* mapToSearch, QString ObjectToFind)
{
    typename std::multimap<QString, T>::iterator it; 
    it = mapToSearch->find(ObjectToFind);
	if ( it == mapToSearch->end() )
    {
        return 0;
    }
    else
    {
        return(it->second);
    }
}

template<class T>
vector<T> FindAllInMultiMap(std::multimap<QString, T>* mapToSearch, QString ObjectToFind)
{
    vector<T> returnValues;
    typename std::multimap<QString, T>::iterator it;

    for(it = mapToSearch->equal_range(ObjectToFind).first; it != mapToSearch->equal_range(ObjectToFind).second; it++)
    {
        returnValues.push_back(it->second);
    }

    return returnValues;
};


#endif
