#include "xmlqt.h"
#include <assert.h>

#include "generalexception.h"
#include "xmlvector.h"
#include "translations.h"
static const char* ERROR_STR2NUM = QT_TRANSLATE_NOOP("xmlqt", "Error: Conversion of  string to numeric type failed");
#define XML_TAG_XMLFILE "xmlFile"
#define XML_ATTR_PATH "path"
#define XML_TAG_FRAGMENT "proost_fragment"
    
XmlQt::XmlQt()
{
}

void 
XmlQt::xmlParseFile(QString xmlFileName, QDomDocument* completeDoc)
{
	
	QString filePath = QDir::cleanPath(xmlFileName);
    
    QFile xmlFile;
    xmlFile.setFileName(filePath);

	if (!xmlFile.open(QIODevice::ReadOnly) )
	{
        QByteArray msg = qPrintable(ERROR_OPENING_INPUT_FILE + QString(": ") + filePath);
		throw GeneralException(PrTr("xmlqt", msg.constData()));
	}
	QString errorMsg = "";
    int errorLine;
    int errorColumn;
	if (!completeDoc->setContent(&xmlFile, &errorMsg, &errorLine, &errorColumn) )
    {
		 xmlFile.close();
         throw GeneralException(PrTr("xmlqt", ERROR_LOADING_XML_FILE) + QString(":%1 - Line: %2 Column: %3").arg(errorMsg).arg(errorLine).arg(errorColumn));
	}

	xmlFile.close();

	return;
}

QString 
XmlQt::getXMLTag(void* curElement)
{
	QDomElement* curElem = static_cast < QDomElement* > (curElement);
	return curElem->tagName();
}

QDomElement 
XmlQt::convertToQDomElement(void* sourceData)
{
	QDomElement* xmlFirstElemPtr = static_cast<QDomElement*> (sourceData);
	return (*xmlFirstElemPtr);
}

QDomElement 
XmlQt::getFirstElemByTag(const QDomElement parentElem, QString elementTag)
{
	QDomNodeList xmlNodeList = parentElem.elementsByTagName(elementTag);
    const int nodeListSize = xmlNodeList.size();
	if( nodeListSize > 0)
	{
		return (xmlNodeList.at(0).toElement());
	}
	else
	{
		return(QDomElement());
	}
}
QDomElement 
XmlQt::getFirstChildElemByTag(const QDomElement parentElem, QString elementTag)
{
    return  parentElem.firstChildElement(elementTag);

}

vector<QDomElement> 
XmlQt::getElementsByTagName(const QDomElement parentElem,
												 QString elementTag)
{
	vector<QDomElement> resElemList;
    assert(!parentElem.isNull());
	QDomNodeList allChildren = parentElem.elementsByTagName(elementTag);
    int numOfChildren = allChildren.size();
	for(int i = 0; i<numOfChildren;i++)
	{
        if(allChildren.at(i).parentNode() == parentElem)
        {
			resElemList.push_back(allChildren.at(i).toElement());
        }
	}
	return resElemList;
}

void 
XmlQt::getXMLChildAttribute(QDomNode parentNode,
                                         QString elementTag,
                                         QString attrDescription,
                                         QVariant &attrValue)
{
    QDomNode curNode = parentNode.firstChild();
    while( !curNode.isNull() )
    {
        if (curNode.isElement() )
        {
            QDomElement curElement = curNode.toElement();
            QString curTag = curElement.tagName();
            if (curTag == elementTag)
            {
                getXMLAttribute(curElement,
                                            attrDescription, 
                                            attrValue);
            }
        }
        curNode = curNode.nextSibling();
    }
}

void 
XmlQt::getXMLChildAttribute(QDomNode parentNode,
                                      QString elementTag,
                                      QString attrDescription,
                                      double &attrValue)
{
    QVariant tmpValue = attrValue;
    getXMLChildAttribute(parentNode, 
                         elementTag, 
                         attrDescription,
                         tmpValue);
    attrValue = tmpValue.toDouble();
}

void 
XmlQt::getXMLChildAttribute(QDomNode parentNode,
                                      QString elementTag,
                                      QString attrDescription,
                                      int &attrValue)
{
    QVariant tmpValue = attrValue;
    getXMLChildAttribute(parentNode, 
                         elementTag, 
                         attrDescription,
                         tmpValue);
    attrValue = tmpValue.toInt();
}

void 
XmlQt::getXMLChildAttribute(QDomNode parentNode,
                                      QString elementTag,
                                      QString attrDescription,
                                      QString &attrValue)
{
    QVariant tmpValue = attrValue;
    getXMLChildAttribute(parentNode, 
                         elementTag, 
                         attrDescription,
                         tmpValue);
    attrValue = tmpValue.toString();
}

void 
XmlQt::getXMLChildAttributes(QDomNode parentNode,
                                       QString elementTag,
                                       QStringList attrDescriptions,
                                       QList<QList<QVariant> > &attrValuesList)
{
    QList<QList<QVariant> > tmpValuesList;
    QDomNode curNode = parentNode.firstChild();
    while( !curNode.isNull() )
    {
        if (curNode.isElement() )
        {
            QDomElement curElement = curNode.toElement();
            QString curTag = curElement.tagName();
            if (curTag == elementTag)
            {
                QList<QVariant> tmpValues;
                tmpValues = attrValuesList.at(0);
                for(int i=0;i<attrDescriptions.size();++i)
                {
                    XmlQt::getXMLAttribute(curElement,
                                                attrDescriptions.at(i), 
                                                tmpValues[i]);
                }
                tmpValuesList.push_back(tmpValues);
            }
        }
        curNode = curNode.nextSibling();
    }
    attrValuesList = tmpValuesList;
}

void 
XmlQt::getXMLChildAttributes(QDomNode parentNode,
                                       QString elementTag,
                                       QStringList attrDescriptions,
                                       QStringList &attrValues)
{
	QDomElement curElement = parentNode.firstChildElement();
    while( !curElement.isNull() )
    {
        QString curTag = curElement.tagName();
        if (curTag == elementTag)
        {
            for(int i=0;i<attrDescriptions.size();i++)
            {
				attrValues << QString();
                XmlQt::getXMLAttribute(curElement,
                                       attrDescriptions.at(i), 
                                       attrValues[i]);
            }
        }
        curElement = curElement.nextSiblingElement();
    }
}

void 
XmlQt::getXMLChildAttributes(QDomNode parentNode,
                                       QString elementTag,
                                       QStringList attrDescriptions,
                                       QList<QVariant> &attrValues)
{
	QDomElement curElement = parentNode.firstChildElement();
    while( !curElement.isNull() )
    {
        QString curTag = curElement.tagName();
        if (curTag == elementTag)
        {
            for(int i=0;i<attrDescriptions.size();++i)
            {
				attrValues << QVariant();
                XmlQt::getXMLAttribute(curElement,
                                       attrDescriptions.at(i), 
                                       attrValues[i]);
            }
        }
        curElement = curElement.nextSiblingElement();
    }
}

void 
XmlQt::getXMLAttributes(QDomNode node, 
							 QStringList attrs, 
							 QStringList &values)
{
	values.clear();
	for(int i=0;i<attrs.size();++i)
	{
		values.push_back(QString());
		XmlQt::getXMLAttribute(node,
							   attrs[i], 
							   values[i]);
	}
}

void 
XmlQt::getXMLAttributes(QDomNode node, 
							 QStringList attrs, 
							 vector<int> &values)
{
	values.clear();
	for(int i=0;i<attrs.size();++i)
	{
		values.push_back(int());
		XmlQt::getXMLAttribute(node,
							   attrs[i], 
							   values[i]);
	}
}

void 
XmlQt::getXMLAttributes(QDomNode node, 
							 QStringList attrs, 
							 vector<double> &values)
{
	values.clear();
	for(int i=0;i<attrs.size();++i)
	{
		values.push_back(double());
		XmlQt::getXMLAttribute(node,
							   attrs[i], 
							   values[i]);
	}
}

void 
XmlQt::getXMLAttribute(QDomNode node, 
                    QString attr, 
                    int& value)
{
    QVariant curValue(value);
    getXMLAttribute(node, 
                    attr, 
                    curValue);
    value = curValue.toInt();
}

void 
XmlQt::getXMLAttribute(QDomNode node,
                    QString attr, 
                    double& value)
{
    QVariant curValue(value);
    getXMLAttribute(node, 
                    attr, 
                    curValue);
    value = curValue.toDouble();
}

void 
XmlQt::getXMLAttribute(QDomNode node, 
                     QString attr, 
                     bool& value)
{
    QString curValue;

    getXMLAttribute(node, 
                    attr,
                    curValue);
    curValue = curValue.trimmed().toLower();
    if (curValue == "true")
    {
        value = true;
    }
    else if  (curValue == "false")
    {
        value = false;
    }
    else if  (curValue == "")
    {
        return;
    }
    else
    {
        QString msg = "Error: trying to convert to boolean the value " + curValue + " of attribute " + attr + " of element " + node.nodeName() +
            " and it failed. Value should be ""true "" or ""false""";
        throw GeneralException("XmlQt", msg);
    }
}

void 
XmlQt::getXMLAttribute(QDomNode node, 
                    QString attr, 
                    QString& value)
{
    QVariant curValue(value);
    getXMLAttribute(node, 
                    attr, 
                    curValue);
    value = curValue.toString();
}

void 
XmlQt::getXMLAttribute(QDomNode node, 
                    QString attr, 
					std::string &value)
{
	QString tmpValue;
    getXMLAttribute(node, 
                    attr, 
                    tmpValue);
    value = tmpValue.toStdString();
}

void 
XmlQt::getXMLAttribute(QDomNode node, 
                    QString attr, 
                    QVariant& value)
{
	QVariant::Type curType = value.type();
    QString tagString = node.toElement().attribute(attr, XML_NOATTRIBUTE);
    if (tagString == XML_NOATTRIBUTE)
	{
        if(curType == QVariant::String 
           || curType == QVariant::Bool)
		{
			value = QString();
		}
	}
    else
    {
        if(curType == QVariant::String
           || curType == QVariant::Bool)
            value = tagString;
        else if(curType == QVariant::Double)
            value = tagString.toDouble();
        else if(curType == QVariant::Int)
            value = tagString.toInt();
    }
}

void 
XmlQt::includeXMLFile(QDomElement &xmlElement,
						   QString inclRef)
{
	QString xmlFile;
	XmlQt::getXMLAttribute(xmlElement, inclRef, xmlFile);
	QDomDocument newDoc;
	XmlQt::xmlParseFile(xmlFile, &newDoc);
	QDomElement firstElement = newDoc.documentElement();
	
	QDomDocument curDoc = xmlElement.ownerDocument();
	QDomNode newNode = curDoc.importNode(firstElement, true);
	QString testoutput = curDoc.toString();
	xmlElement.parentNode().appendChild(newNode);
}

void 
XmlQt::includeAllXmlFiles(QDomElement &rootElement,
							   QString inclTag,
							   QString inclRef)
{
	QDomElement curChild = rootElement.firstChildElement();
	while(!curChild.isNull())
	{
		QString curTag = curChild.tagName();
		if(curTag == inclTag)
		{
			XmlQt::includeXMLFile(curChild, inclRef);

		}
		includeAllXmlFiles(curChild, inclTag, inclRef);
		curChild = curChild.nextSiblingElement();
	}
}



void XmlQt::getXMLChildContent(QDomElement aNode,  QString aTagname, double& aValue, bool& isfound)
{

	vector<QDomElement> elem = XmlQt::getElementsByTagName(aNode, aTagname );
                  
	if (elem.size() == 0) 
		{isfound = false;
		return;}
	else if (elem.size() >1 )
		{throw;}
	else 
		{isfound =true;}

	QString content= elem[0].text();

	bool ok;

	aValue = content.toDouble(&ok);

	if (!ok) {
		throw;
	}

}


		/// \brief returns content b from a child element of the form <a> b <\a>
void XmlQt::getXMLChildContent(QDomElement aNode,  QString aTagname, int& aValue, bool& isfound)
{

	vector<QDomElement> elem = XmlQt::getElementsByTagName(aNode, aTagname );
                  
	if (elem.size() ==0) 
		{isfound = false;
		return;}
	else if (elem.size() >1 )
		{throw;}
	else 
		{isfound =true;}

	QString content= elem[0].text();

	bool ok;

	aValue = content.toDouble(&ok);

	if (!ok) {
		throw;
	}
};
void XmlQt::getXMLChildContent(QDomElement aNode,  QString aTagname, bool& aValue, bool& isfound)
{

	vector<QDomElement> elem = XmlQt::getElementsByTagName(aNode, aTagname );
                  
	if (elem.size() ==0) 
		{isfound = false;
		return;}
	else if (elem.size() >1 )
		{throw;}
	else 
		{isfound =true;}

	QString content= elem[0].text();

	content = content.toLower().trimmed();
    if (QString::compare(content,"true")==0)
    {
        aValue = true;
    }
    else if (QString::compare(content,"false")==0)
    {
        aValue = false;
    }
    else
    {throw;}

};
void
XmlQt::ReadXMLContentVector(QDomElement aNode, vector<double>& aResult )
{

	 //if we have a "vector" inside 
	vector <QDomElement> XVector = XmlQt::getElementsByTagName(aNode,XML_TAG_VECTOR);
	if (XVector.size() == 0)
	{
		QString dummy  = aNode.text().simplified();
		
		foreach(QString elem, dummy.split(" ", QString::SkipEmptyParts)){
            bool ok;
			aResult.push_back(elem.toDouble(&ok));
            if (!ok)
            {
                throw(ERROR_STR2NUM);
            }

		}	
	}
	else
	{
		for (int i=0 ; i < XVector.size() ; i++)
		{
			vector <double> x = XmlVector<double>::GetVector(XVector[i]);
			aResult.insert(aResult.end(),x.begin(),x.end());
		}
	}

};
void
XmlQt::ReadXMLContentVector(QDomElement aNode, vector<int>& aResult)
{

	// if we have a "vector" inside 
	QDomElement XVector = XmlQt::getFirstElemByTag(aNode,XML_TAG_VECTOR);
	if (XVector.isNull())
	{
        QString dummy1  = aNode.text();
        QString dummy2 = dummy1.simplified();
        QStringList temp = dummy2.split(" ", QString::SkipEmptyParts);
		
        foreach(QString elem,  temp){
			bool ok;
            aResult.push_back(elem.toInt(&ok));
            if (!ok)
            {
                throw(ERROR_STR2NUM);
            }
		}	
	}
	else
	{
		aResult = XmlVector<int>::GetVector(XVector);
	}
	

};
void 
XmlQt::ReadXMLContentVector(QDomElement aNode, int aNrCol,vector<vector<double> >& aResult)
{
    aResult.clear();
	vector<double> vector;
    XmlQt::ReadXMLContentVector(aNode,vector);
	
	if (aNrCol == 0) return;


	//nr elem in vector divided by nr o columns should be integer
	if (vector.size() % aNrCol != 0) return;

	int nrRow = vector.size() / aNrCol;
	aResult.resize(nrRow);

	int total = 0;
	for (int i = 0; i < nrRow;i++){
		aResult[i].resize(aNrCol);
		for (int j = 0;j < aNrCol; j++,total++ ){
			aResult[i][j] = vector.at(total);
		}
	}

};
vector<QDomElement>
XmlQt::GetElemsByNameAndAttrib(const QDomElement aNode, const QString aTagname,const QString aAttribname, const QString aAttribValue)
{
    vector<QDomElement> byName;
    vector<QDomElement> out;

    byName = XmlQt::getElementsByTagName(aNode, aTagname);

    for (unsigned i = 0; i < byName.size(); i++)
    {
        if ( byName[i].hasAttribute(aAttribname))
        {
            QString itsAttribval = byName[i].attribute(aAttribname);
            if (QString::compare(aAttribValue, itsAttribval) == 0)
            {
                out.push_back(byName[i]);
            }
        }
    }
    return out;

}
 bool 
XmlQt::hasXMLAttribute(QDomElement node, QString attr)
 {
     return node.hasAttribute(attr);
 }
void
XmlQt::xmlDeepParseFile(QString xmlFileName, QDomElement& completeDoc)
{
    
    // read file
    QDomDocument mydoc("doc");
    XmlQt::xmlParseFile(xmlFileName,&mydoc);
    completeDoc = mydoc.documentElement();

    // look for file tag and substitute
    XmlQt::xmlSubstituteFile(completeDoc, xmlFileName );
}

void 
XmlQt::xmlSubstituteFile(QDomElement xTree,  const QString currentPath)
{   

    QDomNodeList xChildren = xTree.childNodes();
    int nrChildren = (int) xChildren.size();
    
    //loop over elements contained in this one
    for  (int i = 0; i < nrChildren; i++)
    {

        QDomElement xChild = xChildren.at(i).toElement();
        QString dtagnam = xChild.tagName();

        // see if name of element is "XmlFile"
        if (xChild.tagName().toStdString() == XML_TAG_XMLFILE)
        {
            // get  the referenced file name
            QString path = xChild.text().trimmed();



            if ( IsSameFile(path, currentPath ) )
            {
                GeneralException e("Error:filename inside xmlFile element cannot be the name of the file it is in itself."); 
                throw(e);
            }
       
            cout << "now parsing..." << path.toStdString()<< endl; 
            QDomElement xFileTree;
            bool itSubstituted = false;
            XmlQt::xmlDeepParseFile(path, xFileTree);

            vector<QDomElement> nodes_to_add;

            // if the xFileTree contains an XmlFile element as its main tag
            if ( xFileTree.tagName() == XML_TAG_FRAGMENT)
            {
                QDomNodeList xdummy  = xFileTree.childNodes();
                for (int i = 0; i < xdummy.size(); i++)
                {
                    nodes_to_add.push_back(xdummy.at(i).toElement() );
                }
            }
            else
            {
                nodes_to_add.push_back(xFileTree);
            }


            // insert its contend and remove the "XmlFIle" element. 
            for(auto it = nodes_to_add.begin(); it != nodes_to_add.end(); it++)
            {
                xTree.insertAfter(*it,xChild );
            }

            // remove the original "XmlFile"  element 
            xTree.removeChild(xChild);



            cout << "done with..." << path.toStdString()<< endl; 
        }
        else
        {

            XmlQt::xmlSubstituteFile(xChild, currentPath);

        }

    }

}

void 
XmlQt::xmlElemToFile(QDomElement xTree, QString filename)
{
    QFile outfile;

        
    //now try to open the outfile in write mode.
    outfile.setFileName(filename); 
    if (! outfile.open(QIODevice::WriteOnly))
    {
        cout << "the output file you specified cannot be written to."<< endl;;
        exit(1);
    }     

    // now write output to file
    QTextStream outstream(&outfile);
    xTree.save(outstream,4);
    outstream.flush();
    outfile.close();
}


bool 
XmlQt::IsSameFile(QString filename1,QString filename2 )
{
    QDir d1(filename1);
    QDir d2(filename2);



    filename1 = QDir::cleanPath(d1.absolutePath());
    filename2 = QDir::cleanPath(d2.absolutePath());

#ifdef _WIN32
    filename1 = filename1.toLower();
    filename2 = filename2.toLower();    
#else
    // do nothing
#endif

    


    return (filename1 == filename2);
}