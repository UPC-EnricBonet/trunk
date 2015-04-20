#ifndef XMLQT_H
#define XMLQT_H

#define XML_NOATTRIBUTE "__XML_NOATTRIBUTE__"
#define XML_ATTR_TYPE "type"

#include <QtXml/QtXml>
#include <QString>
#include <QStringList>
#include <QVariant>
#include <QList>
#include <iostream>

#include <vector>
using namespace std;

static const char* ERROR_OPENING_INPUT_FILE = QT_TRANSLATE_NOOP("xmlqt", "Error opening inputfile");
static const char* ERROR_LOADING_XML_FILE = QT_TRANSLATE_NOOP("xmlqt", "Error loading xml file");


class XmlQt
{
    public:
        XmlQt();

		static void xmlParseFile(QString xmlFileName, QDomDocument* completeDoc);

		static QString getXMLTag(void* curElement);

		static QDomElement convertToQDomElement(void* sourceData);

		static QDomElement getFirstElemByTag(QDomElement parentElem, 
											 QString elementTag);

        static QDomElement getFirstChildElemByTag(const QDomElement parentElem, QString elementTag);

		static vector<QDomElement> getElementsByTagName(const QDomElement parentElem,
														 QString elementTag);

		static void getXMLAttributes(QDomNode node, 
									 QStringList attrs, 
									 QStringList &values);

		static void getXMLAttributes(QDomNode node, 
									 QStringList attrs, 
									 vector<int> &values);

		static void getXMLAttributes(QDomNode node, 
									 QStringList attrs, 
									 vector<double> &values);

        static void getXMLChildAttribute(QDomNode parentNode,
                                        QString elementTag,
                                        QString attrDescription,
                                        QVariant &attrValue);

        static void getXMLChildAttribute(QDomNode parentNode,
                                        QString elementTag,
                                        QString attrDescription,
                                        double &attrValue);

        static void getXMLChildAttribute(QDomNode parentNode,
                                        QString elementTag,
                                        QString attrDescription,
                                        int &attrValue);

        static void getXMLChildAttribute(QDomNode parentNode,
                                        QString elementTag,
                                        QString attrDescription,
                                        QString &attrValue);

        static void getXMLChildAttributes(QDomNode parentNode,
                                      QString elementTag,
                                      QStringList attrDescriptions,
                                      QStringList &attrValues);

        static void getXMLChildAttributes(QDomNode parentNode,
                                      QString elementTag,
                                      QStringList attrDescriptions,
                                      QList<QVariant> &attrValues);

        static void getXMLChildAttributes(QDomNode parentNode,
                                      QString elementTag,
                                      QStringList attrDescriptions,
                                      QList<QList<QVariant> > &attrValuesList);

        static void getXMLAttribute(QDomNode node, 
                   QString attr, 
                   int& value);

        static void getXMLAttribute(QDomNode node, 
                   QString attr, 
                   double& value);

        static void getXMLAttribute(QDomNode node, 
                   QString attr, 
                   bool& value);

        static void getXMLAttribute(QDomNode node, 
                   QString attr, 
                   QString& value);

		static void getXMLAttribute(QDomNode node, 
                   QString attr, 
				   std::string &value);

        static void getXMLAttribute(QDomNode node, 
                   QString attr,
                   QVariant& value);


        static bool hasXMLAttribute(QDomElement node, 
                   QString attr);

		static void includeXMLFile(QDomElement &xmlElement,
								   QString inclRef);

		static void includeAllXmlFiles(QDomElement &rootElement,
									   QString inclTag,
									   QString inclRef);

		/// \brief returns content b from a child element of the form <a> b <\a>
		/// if element <a> is not present, aVaule is not changed. 
		/// if b is not convertible to the output type, an exception is thrown. 
		static void getXMLChildContent(QDomElement aNode,
			QString aTagname,
			double& aValue,
			bool& isfound);



		/// \brief returns content b from a child element of the form <a> b <\a>
		static void getXMLChildContent(QDomElement aNode,
			QString aTagname,
			int & aValue,
			bool & isfound);

		/// \brief returns content b from a child element of the form <a> b <\a>
		/// if element <a> is not present, aVaule is not changed. 
		/// if b is not convertible to the output type, an exception is thrown. 
		static void getXMLChildContent(QDomElement aNode,
            QString aTagname,
            bool& aValue, 
            bool& isfound);

		/// \brief Reads a vector of double from xml content
		static void ReadXMLContentVector(QDomElement aNode, vector<double>&  aResult);

		/// \brief Reads a vector of int from xml content
		static void ReadXMLContentVector(QDomElement aNode, vector<int>&  aResult);

		static void ReadXMLContentVector(QDomElement aNode, int nrCol,vector<vector<double> >& aResult);

		/// \brief Get elements by name and value (name = NAME, attribute = VALUE
        static vector<QDomElement> GetElemsByNameAndAttrib(const QDomElement aNode,
            const QString aTagname, 
            const QString aAttribname,
            const QString aAttribValue);


        /// \brief Parses a file that optionally contains file references.
        static void xmlDeepParseFile(QString xmlFileName, QDomElement& completeDoc);

        ///\brief substitues a file reference for the tree in the file
        static void xmlSubstituteFile(QDomElement xTree, const QString currentPath);

        ///\brief substitues a file reference for the tree in the file
        static void xmlElemToFile(QDomElement xTree, QString filename);


        ///\brief substitues a file reference for the tree in the file
        static bool IsSameFile(QString filename1,QString filename2 );
};

#endif
