#ifndef GENERALEXCEPTION_H
#define GENERALEXCEPTION_H

#include <exception>
#include <QString>
#include <string>

class GeneralException :
	public std::exception
{
private:
	QString mMessage;
public:
	QString Message(){ return this->mMessage; }
    std::string MessageString(){ return this->mMessage.toStdString(); }
    void PrefixMessage(QString prefixString) {mMessage = prefixString + mMessage;}
public:
	bool operator== ( const char * other ) const { return (other == mMessage); }
    bool operator== ( QString& other ) const { return (other == mMessage); }
    bool operator== ( GeneralException& other ) const { return (other.Message() == mMessage); }
public:
    GeneralException(){}

    GeneralException(QString aMessage){ this->mMessage = aMessage; }

    GeneralException(QString aMessage, QString aMessage2){ this->mMessage = aMessage + "." + aMessage2; }

    GeneralException(char* aMessage){ this->mMessage = aMessage; }

    GeneralException(char* aMessage, char* aMessage2){ this->mMessage = QString(aMessage) + "." + QString(aMessage2); }

    virtual ~GeneralException() throw(){}
};

class NumericalException :
	public std::exception
{
private:
	QString mMessage;
public:
	QString Message(){ return this->mMessage; }
    std::string MessageString(){ return this->mMessage.toStdString(); }
    void PrefixMessage(QString prefixString) {mMessage = prefixString + mMessage;}
public:
	bool operator== ( const char * other ) const { return (other == mMessage); }
    bool operator== ( QString& other ) const { return (other == mMessage); }
    bool operator== ( GeneralException& other ) const { return (other.Message() == mMessage); }
public:
    NumericalException(){}

    NumericalException(QString aMessage){ this->mMessage = aMessage; }

    NumericalException(QString aMessage, QString aMessage2){ this->mMessage = aMessage + "." + aMessage2; }

    NumericalException(char* aMessage){ this->mMessage = aMessage; }

    NumericalException(char* aMessage, char* aMessage2){ this->mMessage = QString(aMessage) + "." + QString(aMessage2); }

    virtual ~NumericalException() throw(){}
};


#endif
