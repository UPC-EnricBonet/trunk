#include "base.h"

#include <QFile>
#include <QTextStream>

//#include "xmlconstants.h"
#include "serialization.h"
#include "deserialization.h"
#include "objectmanager.h"

Base::Base()
{
    this->mClassName = "Base";
	this->mName = "";
    this->mIsAllRead = false;
}

Base::Base(Base* aProost)
{
    this->mClassName = "Base";
    this->mName = aProost->mName;
    this->mIsAllRead = false;
}

Base::Base(QDomElement aNode)
{
    this->mClassName = "Base";
	this->mNodeFound = aNode;
	this->mIsAllRead = false;
	this->mName = "";
	if ( ! aNode.isNull() )
	{
        QString atrName;
        XmlQt::getXMLAttribute(aNode, "name", atrName);
		if ( atrName.isEmpty() )
        {
            this->mName = this->GenerateIdentifier(this->mClassName);
		}
        else
        {
		    this->mName = atrName.trimmed();
        }	
	}
}

void 
Base::PrThrow(const char* message)
{
    QString curObject = mClassName + "(" + mName + ").";
    throw GeneralException(QCoreApplication::translate(qPrintable(mClassName), message));
}

void 
Base::PrThrow(GeneralException &aException)
{
    QString curObject = mClassName + "(" + mName + ").";
    aException.PrefixMessage(curObject);  
    throw aException;
}

void 
Base::Read()
{
	if ( this->mNodeFound.isNull() ) 
	{
		this->PrThrow(ERROR_UNINITIALIZED_XML_NODE);
	}
    try
	{
        if(!this->mIsAllRead) this->Read(this->mNodeFound);
    }
    catch(GeneralException& ex)
	{
		this->PrThrow(ex);
	}
}

void 
Base::Read(const QDomElement aNode)
{
	try
	{
        pRead(aNode);
        this->mIsAllRead = true;
     }
    catch(GeneralException& ex)
	{
		this->PrThrow(ex);
	}
}

void 
Base::pRead(const QDomElement aNode)
{
    if ( aNode.isNull() ) this->PrThrow(ERROR_XML_NODE_CANNOT_BE_NULL);

    if ( this->mNodeFound != aNode || !mIsAllRead  )
	{
		this->mNodeFound = aNode;

        // check if the name has not been set already. 
        if ( this->mName.isEmpty())
        {
            QString atrName;
            XmlQt::getXMLAttribute(aNode, XML_ATTR_NAME, atrName);


		    this->mName = atrName.trimmed();
		    if ( this->mName == "" ) this->PrThrow(ERROR_NAME_ATTRIBUTE_CANNOT_BE_VOID);
        }
	}


	
}

void 
Base::ReadReferences(const QDomElement aNode,
                        QString refType, 
                        vector<QDomElement> &refTags,
                        vector<QString>* aListOfIDs)
{
    vector<QDomElement> allRefTags = XmlQt::getElementsByTagName(aNode, XML_TAG_REF);
    for (unsigned i = 0 ; i < allRefTags.size(); i++ )
	{
        QString attrRefType;
        XmlQt::getXMLAttribute(allRefTags[i], XML_ATTR_TYPE, attrRefType);
        attrRefType = attrRefType.trimmed();
        if(attrRefType.isEmpty())
		{
            this->PrThrow(ERROR_REF_HAS_NO_TYPE);
		}
        if(attrRefType == refType)
        {
            refTags.push_back(allRefTags[i]);
            QString atrID;
		    XmlQt::getXMLAttribute(allRefTags[i], XML_ATTR_ID, atrID);
            if(aListOfIDs != 0)
            {
                aListOfIDs->push_back(atrID);
            }
        }
    }
}

void 
Base::ReadFirstReferences(const QDomElement aNode,
                             QString refType, 
                             QDomElement &refTag,
                             QString* atrID)
{
    vector<QDomElement> allRefTags = XmlQt::getElementsByTagName(aNode, XML_TAG_REF);
    for (unsigned i = 0 ; i < allRefTags.size(); i++ )
	{
        QString attrRefType;
        XmlQt::getXMLAttribute(allRefTags[i], XML_ATTR_TYPE, attrRefType);
        attrRefType = attrRefType.trimmed();
        if(attrRefType.isEmpty())
		{
            this->PrThrow(ERROR_REF_HAS_NO_TYPE);
		}

        if(attrRefType == refType)
        {
            refTag = allRefTags[i];
            if(atrID != 0)
            {
		        XmlQt::getXMLAttribute(allRefTags[i], XML_ATTR_ID, *atrID);
            }
            break;
        }
    }
}

void 
Base::Copy(Base &aTarget, Base &aSource)
{
    aTarget.mClassName = aSource.mClassName;
    aTarget.mIsAllRead = aSource.mIsAllRead;
    aTarget.mName      = aSource.mName;
    aTarget.mNodeFound = aSource.mNodeFound;
}

QString
Base::GenerateIdentifier(const QString& root)
{
    gTimeStamp++;
    QString newname = root +QString::number(gTimeStamp);
    return newname;
}



bool 
Base::addRef(Base* proostObj)
{
	this->PrThrow(ERROR_THIS_ISNOT_IMPLEMENT_TO_THIS_CLASS_GUI);

	return false;
}

GuiItem*
Base::getGuiItem()
{
	return this->mGuiItem;
}

void
Base::setGuiItem(GuiItem* aGuiItem)
{
	this->mGuiItem = aGuiItem;
}


