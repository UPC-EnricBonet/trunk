#ifndef CDESERIALIZATION_H
#define CDESERIALIZATION_H

#include <QString>
#include <vector>
#include <map>
#include "xmlqt.h"

using namespace std;

class CDeserialization
{
private:
	QDomElement mNode;
	QString mDesClassName;
public:
	CDeserialization(QString aClassName,QString* aIn)
	{
		mDesClassName = aClassName;
		QDomDocument* doc = new QDomDocument("mydocument");
		doc->setContent(*aIn);
		mNode = doc->documentElement();
		if (mNode.tagName() != this->mDesClassName )
		{
			mNode = XmlQt::getFirstElemByTag(mNode, this->mDesClassName);
			if ( mNode.isNull() ) throw;
		}
	}
	CDeserialization(QString aClassName,QDomElement* aNode)
	{
		mDesClassName = aClassName;
		mNode = *aNode;
		if (mNode.tagName() != this->mDesClassName )
		{
			mNode = XmlQt::getFirstElemByTag(mNode, this->mDesClassName);
			if ( mNode.isNull() ) throw;
		}
	}

public:

    static void StringToElement(QString& xmlString, QDomElement& node)
    {
        QDomDocument* doc = new QDomDocument("mydocument");
		doc->setContent(xmlString);
		node = doc->documentElement();
    }

    static void VerifyAttrName(QDomElement &elem, QString attrName)
    {
        QString tagName = elem.tagName();
        if (tagName != attrName )
		{
			elem = XmlQt::getFirstElemByTag(elem, attrName);
			if ( elem.isNull() ) throw;
		}
    }

	void Begin()
	{
		//do nothing
	}
	void End()
	{
		//do nothing
	}

    static void DeserializeAttr(QDomElement &elem, QString attrName, int& iAttribute)
	{
		QDomElement node = XmlQt::getFirstElemByTag(elem,attrName);
		if ( node.isNull() ) 
        { 
            throw;
        }
		iAttribute = node.text().toInt();
	}

    static void DeserializeAttr(QDomElement &elem, QString attrName, long& iAttribute)
	{
		QDomElement node = XmlQt::getFirstElemByTag(elem,attrName);
		if ( node.isNull() ) 
        { 
            throw;
        }
		iAttribute = node.text().toLong();
	}

    static void DeserializeAttr(QDomElement &elem, QString attrName, QString& iAttribute, bool getByTag = true)
	{
        QDomElement node = elem;
        if(getByTag)
        {
		    node = XmlQt::getFirstElemByTag(elem,attrName);
        }
		if ( node.isNull() ) 
        {
            throw;
        }
        iAttribute = node.text();
	}

    static void DeserializeAttr(QDomElement &elem, QString attrName, bool& iAttribute)
	{
		QDomElement node = XmlQt::getFirstElemByTag(elem,attrName);
		if ( node.isNull() ) 
        {
            throw;
        }
		iAttribute = (bool) node.text().toInt();
	}

    template<typename T>
    static void DeserializeAttr(QDomElement &elem, QString attrName, T* &iAttribute)
	{
        QDomElement node = XmlQt::getFirstElemByTag(elem,attrName);
        if ( node.isNull() ) 
        {
            iAttribute = 0;
        }
        else
        {
            iAttribute = new T();
            QString bla = node.text();
            iAttribute->Deserialize(node, attrName);
        }

    }

    static void DeserializeAttr(QDomElement &elem, QString attrName, vector<int>& iAttribute, bool getByTag = true)
	{
        iAttribute.clear();
        QDomElement node = elem;
        if(getByTag)
        {
            node = XmlQt::getFirstElemByTag(elem,attrName);
        }
        if ( node.isNull() ) 
        {
            throw;
        }
        QString nodeString = node.text();
        QStringList strList = nodeString.split(" ");
        // If there is one item and it is empty, do not execute the following statements
        if( ! (strList.size() == 1 && strList.at(0).isEmpty()))
        {
            iAttribute.resize(strList.size());
            for(int i=0; i< strList.size(); i++)
            {
                iAttribute[i] = strList[i].toInt();
            }
        }
	}

    static void DeserializeAttr(QDomElement &elem, QString attrName, vector<double>& iAttribute, bool getByTag = true)
	{
        iAttribute.clear();
        QDomElement node = elem;
        if(getByTag)
        {
            node = XmlQt::getFirstElemByTag(elem,attrName);
        }
        if ( node.isNull() ) 
        {
            throw;
        }
        QString nodeString = node.text();
        QStringList strList = nodeString.split(" ");
        // If there is one item and it is empty, do not execute the following statements
        if( ! (strList.size() == 1 && strList.at(0).isEmpty()))
        {
            iAttribute.resize(strList.size());
            for(int i=0; i< strList.size(); i++)
            {
                iAttribute[i] = strList[i].toDouble();
            }
        }
	}

    template<typename T>
    static void DeserializeAttr(QDomElement &elem, QString attrName, vector<vector<T> >& iAttribute, bool getByTag = true)
	{
        iAttribute.clear();
        QDomElement node = elem;
        if(getByTag)
        {
            node = XmlQt::getFirstElemByTag(elem,attrName);
        }
		if ( node.isNull() ) 
        {
            throw;
        }
		QString allText = node.text();
        QStringList strList = allText.split(";");
        iAttribute.resize(strList.size());
        for(int i=0; i< strList.size(); i++)
        {
            QDomDocument doc( "tmpdoc" );
            QDomElement tmpNode = doc.createElement("tmp");
            tmpNode.appendChild(doc.createTextNode(strList[i]));   
            DeserializeAttr(tmpNode, "", iAttribute[i], false);
        }
	}

    template<typename T>
    static void DeserializeAttr(QDomElement &elem, QString attrName, vector<vector<vector<T> > >& iAttribute, bool getByTag = true)
	{
        iAttribute.clear();
        QDomElement node = elem;
        if(getByTag)
        {
            node = XmlQt::getFirstElemByTag(elem,attrName);
        }
		if ( node.isNull() ) 
        {
            throw;
        }
		QString allText = node.text();
        QStringList strList = allText.split(":");
        iAttribute.resize(strList.size());
        for(int i=0; i< strList.size(); i++)
        {
            QDomDocument doc( "tmpdoc" );
            QDomElement tmpNode = doc.createElement("tmp");
            tmpNode.appendChild(doc.createTextNode(strList[i]));   
            DeserializeAttr(tmpNode, "", iAttribute[i], false);
        }
	}

    template<typename T>
    static void DeserializeAttr(QDomElement &elem, QString attrName, vector<vector<vector<vector<T> > > >& iAttribute, bool getByTag = true)
	{
        iAttribute.clear();
        QDomElement node = elem;
        if(getByTag)
        {
            node = XmlQt::getFirstElemByTag(elem,attrName);
        }
		if ( node.isNull() ) 
        {
            throw;
        }
		QString allText = node.text();
        QStringList strList = allText.split("_");
        iAttribute.resize(strList.size());
        for(int i=0; i< strList.size(); i++)
        {
            QDomDocument doc( "tmpdoc" );
            QDomElement tmpNode = doc.createElement("tmp");
            tmpNode.appendChild(doc.createTextNode(strList[i]));   
            DeserializeAttr(tmpNode, "", iAttribute[i], false);
        }
	}

    template<typename T>
    static void DeserializeAttr(QDomElement &elem, QString attrName, QString itemType, vector<T*> &iAttribute)
	{
        iAttribute.clear();
        QDomElement node = XmlQt::getFirstElemByTag(elem,attrName);
		if (! node.isNull() )
        {
            QDomElement nextItem = node.firstChildElement();
		    while(!nextItem.isNull())
            {
                if( nextItem.tagName() == itemType + "Ref")
                {
                    T* newItem = new T();
                    newItem->SetName(nextItem.text());
                    iAttribute.push_back(newItem);
                }
                else if(nextItem.tagName() == itemType)
                {
                    T* newItem = 0;
                    if(nextItem.hasChildNodes())
                    {
                        newItem = new T();
                        newItem->Deserialize(nextItem, itemType);
                    }
                    iAttribute.push_back( newItem );
                }
                nextItem = nextItem.nextSiblingElement();
            }

        }
    }

	template<typename T>
	static void DeserializePtrMap(QDomElement &elem, QString attrName, map<QString,T*> &aMap)
	{
		QDomElement node = XmlQt::getFirstElemByTag(elem,attrName);
		if ( node.isNull() ) 
        {
            throw;
        }
		vector<QDomElement> itemTags = XmlQt::getElementsByTagName(node,"mapItem");
		for ( unsigned i = 0 ; i < itemTags.size() ; i ++ )
		{
			QDomElement first = XmlQt::getFirstElemByTag(itemTags[i],"first");
			QString key = first.text();
			QDomElement second = XmlQt::getFirstElemByTag(itemTags[i],"second");
			T* value = new T;
			value->Deserialize(second,"second");
			aMap.insert( make_pair(key,value) );
		}
	}

};

#endif
