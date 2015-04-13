#include "objectmanager.h"



ObjectManager* ObjectManager::mSingleton = 0;
long ObjectManager::gTimeStemp = 0;
double ObjectManager::mTimeTolerance = 1e-5;
double ObjectManager::mSpaceTolerance = 1e-5;


ObjectManager::ObjectManager()
{ 
    mAllLists = new map<QString,ObjectList<Base>*>; 
}

ObjectManager::~ObjectManager()
{
    map<QString,ObjectList<Base>*>::iterator ilist;
    for (ilist =this->mAllLists->begin(); ilist != this->mAllLists->end();ilist++)
    {
        if (ilist->second != 0)
        {
            delete ilist->second;
        }
    }
    delete mAllLists;
}

ObjectManager* 
ObjectManager::instance()
{
        // prevent from having more than one singleton
        if ( mSingleton == 0 ) 
        {
                mSingleton = new ObjectManager();
                mSingleton->CreateNewFactory();
        }

        return mSingleton;
}


void 
ObjectManager::instDestroy()
{
    if (mSingleton != 0)
    {
        if (mSingleton->mFactory != 0) delete mSingleton->mFactory;
        delete mSingleton;
        mSingleton = 0;
    }
}




void ObjectManager::UnregisterProostClass(QString className)
{
    ObjectManager::mFactory->Unregister(className);
}


void ObjectManager::UnregisterAllProostClasses()
{
        mFactory->UnregisterAll();
}


Base* ObjectManager::newInstanceOfClassName(QString aClassName, QString aName, bool mustExist )
{
    // El nodo no puede ser nulo
    if ( aName == "" ) throw GeneralException(ERROR_NAME_ATTRIBUTE_CANNOT_BE_VOID);

    // Obtiene la nueva instancia con creación sin nodo
    QDomElement aNullNode;
    Base *ret = this->mFactory->Create(aClassName, aNullNode);

    //Coloca el aName en el name de la nueva instancia.
    if ( ret != 0)
    {
            if ( aName == "" ) throw GeneralException(ERROR_NAME_ATTRIBUTE_CANNOT_BE_VOID);
            ret->SetName(aName);    
    }
    else
    {
        if ( mustExist )
        {
                throw GeneralException(ERROR_CLASS_NOT_REGISTERED + aClassName);
        }
    }
    return ret;
}





void ObjectManager::EraseProostListBaseClass(void)
{
	map<QString,ObjectList<Base>* >::iterator it;
	for ( it = mAllLists->begin() ; it != mAllLists->end() ; it ++ )
	{
		ObjectList<Base>* plist = (*it).second;
		delete plist;
	}
	mAllLists->erase( mAllLists->begin(), mAllLists->end() );
}


    Base* ObjectManager::newInstanceFromXML(QDomElement aNode, bool mustExist)
{
	// Obtiene el atributo type que indica el nombre de la clase derivada a instancias.
    QString atrType = aNode.attribute(XML_ATTR_TYPE, XML_NOATTRIBUTE).trimmed();
    atrType = atrType.trimmed();
    if ( atrType == XML_NOATTRIBUTE )
    {
        QString atrName = aNode.attribute("name", XML_NOATTRIBUTE).trimmed();
        atrName = atrName.trimmed();
        QString message;
        message.append(ERROR_TYPE_ATTRIBUTE_NOT_FOUND).append("name = ").append(atrName);
        throw GeneralException(message);
    }

    Base* ret = this->CreateNewObjectInstance(atrType, aNode);

	if ( ret == 0 && mustExist )
	{
		throw GeneralException(ERROR_CLASS_NOT_REGISTERED + atrType);
	}
	return ret;
}

Base* ObjectManager::newInstanceFromXML(QDomElement aNode, QString aClassName, bool mustExist )
{
    // Obtiene el atributo type que en este caso es la clase parent, porque no tiene especializaciones.
	QString atrType = aClassName;

    Base* ret = this->CreateNewObjectInstance(atrType, aNode);

	if ( ret == 0 && mustExist )
	{
		throw GeneralException(ERROR_CLASS_NOT_REGISTERED + atrType);
	}
	return ret;
}

Base* ObjectManager::FindInstance (QString aClassName, QString aName)
{
    ObjectList<Base>* curList;
    QString errorMsg = "Error: a proostList with objects of type " + aClassName + "could not be found";
    curList = FindInMap(mAllLists, aClassName, errorMsg);

	Base* inst = curList->Find(aName);
	if ( inst != 0 )
	{
		if ( ! inst->isAllRead() )
		{
			inst->Read();
		}
	}
	else
	{
        throw GeneralException(QString("Name not found") + " in list: \"" + aClassName + "\"" + " named: " + aName);
	}
	return inst;
}

bool ObjectManager::IsInList(QString aClassName, QString aName)
{
    ObjectList<Base>* curList;
    curList = FindInMap(mAllLists, aClassName, ERROR_CLASS_NAME_NOT_FOUND);

	Base* inst = curList->Find(aName);
	if ( inst != 0 )
	{
        return true;
	}
	return false;
}

Base* ObjectManager::AddToListOfClassName(QString aClassName, QString aName)
{
    if ( aClassName == "" ) throw GeneralException("AddToListOfClassName: class name cannot be empty");
	if ( aName == "" ) throw GeneralException("AddToListOfClassName: object name cannot be empty");

    // Find the correct list
    ObjectList<Base>* aProostList;
    aProostList = FindInMap(mAllLists, aClassName, ERROR_CLASS_NAME_NOT_FOUND);

	// Add a new instance
	Base* retVal = this->newInstanceOfClassName(aClassName, aName, true);
	if ( retVal )
	{
        aProostList->Add(retVal->name(), retVal);
	}
	else
	{
		throw GeneralException("Name not found");
	}

	return retVal;
}

Base* ObjectManager::AddToListOfClassName(QString aListName, QString aClassName, QString aName)
{
	if ( aClassName == "" ) throw GeneralException("AddToListOfClassName: class name cannot be empty");
	if ( aName == "" ) throw GeneralException("AddToListOfClassName: object name cannot be empty");

    // Find the correct list
    ObjectList<Base>* aProostList;
    aProostList = FindInMap(mAllLists, aListName, ERROR_CLASS_NAME_NOT_FOUND);

	// Add a new instance
	Base* retVal = this->newInstanceOfClassName(aClassName, aName, true);
	if ( retVal )
	{
		aProostList->Add(retVal->name(), retVal);
	}
	else
	{
		throw GeneralException("Name not found");
	}

	return retVal;
}

void ObjectManager::AddToListOfClassName(QString aClassName, Base* aObjectToAdd)
{
	if ( aClassName == "" ) throw GeneralException("AddToListOfClassName: class name cannot be empty");
        
    // Find the correct list
    ObjectList<Base>* aProostList;
    aProostList = FindInMap(mAllLists, aClassName, ERROR_CLASS_NAME_NOT_FOUND);

	// Add a new instance
	aProostList->Add(aObjectToAdd->name(),aObjectToAdd);
}

Base* ObjectManager::AddToListFromXML(QString aClassName, QDomElement aInitialNode)
{
	if ( aClassName == "" ) throw GeneralException("AddToListOfClassName: class name cannot be empty");

	// Find correct list
    ObjectList<Base>* aProostList;
    QString msg = "Error: looking for ProostList with name "+ aClassName + " but did not find it";
    aProostList = FindInMap(mAllLists, aClassName, msg);

	// Add new instance
	Base* retVal = this->newInstanceFromXML( aInitialNode );
	if ( retVal )
	{
		aProostList->Add(retVal->name(), retVal);
		retVal->Read(aInitialNode);
	}
	else
	{
		throw GeneralException("Name not found");
	}

    return retVal;
}


void ObjectManager::EraseAllLists()
{
    try
    {
	    for ( map<QString,ObjectList<Base>*>::iterator it = mAllLists->begin(); it != mAllLists->end() ; it++ )
	    {

		    (*it).second->EraseList();
                 
            delete (*it).second;
	    }
	    mAllLists->erase( mAllLists->begin() , mAllLists->end() );
    }
    catch(...)
    {
        cout<< "Proost appears to have finished comnputations correctly, but an error ocurred while freeing memory. This should be debugged."; 
    }
}


void ObjectManager::FindFirstList(const QDomElement aNode,
                                QString listName, 
                                QDomElement &listTag)
{
    vector<QDomElement> allListTags = XmlQt::getElementsByTagName(aNode, XML_TAG_PROOSTLIST);
    for(unsigned i = 0 ; i < allListTags.size(); i++ )
	{
        QString attrType;
        XmlQt::getXMLAttribute(allListTags[i], XML_ATTR_TYPE, attrType);
        attrType = attrType.trimmed();
        if(attrType.isEmpty())
		{
            throw GeneralException("the type attribute of a Proostlist is absent");
		}
        if(attrType == listName)
        {
            listTag = allListTags[i];
            break;
        }
    }
}


 void  ObjectManager::FillAllListsFromXML(QDomElement aInitialNode, bool typeAttributeMustExist, bool mustDelete)
{
    // if we must delete  all, do it straight away
    if (mustDelete)
    {
	    for ( map<QString,ObjectList<Base>*>::iterator it = mAllLists->begin(); it != mAllLists->end() ; it++ )
	    {
            it->second->EraseList();
	    }
    }

    // then read all, but do not delete 

	for ( map<QString,ObjectList<Base>*>::iterator it = mAllLists->begin(); it != mAllLists->end() ; it++ )
	{
		FillListFromXML((*it).second,aInitialNode, typeAttributeMustExist, false);
	}

	// Now that it has all instances allocated let's read all that must be read

	for (map<QString,ObjectList<Base>*>::iterator it = mAllLists->begin(); it != mAllLists->end() ; it++ )
	{
		ObjectList<Base>* pl = (*it).second;
		if ( pl->mMustRead )
		{
			pl->ReadList();
		}
	}
}


void  ObjectManager::FillListFromXML(ObjectList<Base>* aProostList, QDomElement aInitialNode, bool typeAttrMustExist, bool mustDelete)
{
	if ( aProostList == 0 ) throw GeneralException("Internal error");

	if ( mustDelete == true )
	{
		aProostList->EraseList();
	}

	// Si m_tagListName no es nulo ("") entonces busca el primer nodo descendiente de aInitialNode con 
	// el tag m_tagListName y sino se asume que el nodo es aInitialNode.

    QDomElement currentList;
    if( !aProostList->mTagListName.isEmpty() )
	{   
        QString initNodeTag = aInitialNode.tagName();
        this->FindFirstList(aInitialNode, aProostList->mTagListName, currentList);
        QString curListTag = currentList.tagName();
    }
    if(currentList.isNull())
    {
        return;
    }

    vector<QDomElement> listElements;
    listElements = XmlQt::getElementsByTagName(currentList, aProostList->mTagName);
	for(unsigned i = 0 ; i < listElements.size() ;i++)
	{
        if ( typeAttrMustExist )
        {
		    Base* retVal = this->newInstanceFromXML(listElements.at(i));
		    if(retVal != 0)
		    {
			    aProostList->Add(retVal->name(), retVal);
		    }
		    else
		    {
			    throw GeneralException("Name not found");
		    }
        }
        else
        {
            Base* retVal = this->newInstanceFromXML(listElements.at(i), aProostList->mTagName);
		    if(retVal != 0)
		    {
			    aProostList->Add(retVal->name(), retVal);
		    }
		    else
		    {
			    throw GeneralException("Name not found");
		    }
        }
 

	}
}
