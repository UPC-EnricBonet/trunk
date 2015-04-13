#ifndef OBJECTMANAGER_H
#define OBJECTMANAGER_H



#include <QString>

#include "base.h"
#include "translations.h"
#include "generalexception.h"
#include "objectfactory.h"
#include "objectlist.h"
#include "commonglobals.h"
 
static const char* ERROR_CLASS_NOT_REGISTERED = QT_TRANSLATE_NOOP("ObjectManager", "Class not registered:");
static const char* ERROR_CLASS_NAME_NOT_FOUND = QT_TRANSLATE_NOOP("ObjectManager", "Class name not found");


    /// \brief 
/// This class can be used to manage the creation, storage and destruction of any number of elements
/// that derive from "Base". The objects are stored on objectlists. 



class ObjectManager
{
    
private:
    // singleton variable
        static ObjectManager* mSingleton;
   protected:
    // Proost object factory
    ObjectFactory<Base*(QDomElement),QString>* mFactory;
public:

    /// \brief timestamp 
    static long gTimeStemp;    
    static double mTimeTolerance;
    static double mSpaceTolerance;

    ObjectManager();

    virtual ~ObjectManager();
    
    /// \brief get an instance of the object manager
    static ObjectManager* instance();

    /// \brief destroy the instance of the object manager
    static void instDestroy();

    /// \brief Register a class 
    template<class classType>
    void RegisterProostClass(QString className)
    {
        mFactory->Register<classType>(className);
    }


    /// \brief Unregister a class 
    void UnregisterProostClass(QString className);

    /// \brief Unregister all classes 
    void UnregisterAllProostClasses();

    /// \brief create a new instance of a class given the name of it
    Base* newInstanceOfClassName(QString aClassName, QString aName, bool mustExist = true);
 

    /// \brief Creates a proostlist for objects of a certain type, if it does not exist already.
    template<class classType>
	void RegisterProostListBaseClass(QString className, QString tagListName, QString tagName, bool aMustRead = false)
    {
	    // If it doesn't exist create proostList and add to the map of maps
        bool listDoesNotExists = (mAllLists->find(className) == mAllLists->end());
	    if(listDoesNotExists) 
	    {
		    ObjectList<Base>* plist = (ObjectList<Base>*) new ObjectList<classType>(tagListName,tagName, aMustRead);
		    mAllLists->insert( make_pair( className , plist ) );
	    }
    }


    /// \brief remove the proostlist of a certain object type
    void EraseProostListBaseClass(void);

    /// \brief Create new instance (and put it on list) given XML node
    Base* newInstanceFromXML(QDomElement aNode, bool mustExist = true);

    /// \brief Create new instance (and put it on list) given XML node and class name
    Base* newInstanceFromXML(QDomElement aNode, QString aClassName, bool mustExist = true);

    /// \brief Finds instance in lists given class name and object name
    Base* FindInstance (QString aClassName, QString aName);

    /// \brief Indicates if item with given class name and object name is on the proostlist. 
   bool IsInList(QString aClassName, QString aName);

    /// \brief creates an object on a proostlist with a given name and class type. 
    Base* AddToListOfClassName(QString aClassName, QString aName);

    /// \brief creates an object on a proostlist with a given name and class type.
    Base* AddToListOfClassName(QString aListName, QString aClassName, QString aName);

    /// \brief puts an already created object on the proper proostlist. 
    void AddToListOfClassName(QString aClassName, Base* aObjectToAdd);


    /// \brief puts an item on a proostlist given its XML node
    Base* AddToListFromXML(QString aClassName, QDomElement aInitialNode);


    /// \brief Erases all lists
    virtual void EraseAllLists();



    /// \brief find proostlist of objects of a given  type
    void FindFirstList(const QDomElement aNode,
                                  QString listName, 
                                  QDomElement &listTag);


    /// \brief reads  input from xml files (it stores name and a pointer to the xml node.)
    virtual void FillAllListsFromXML(QDomElement aInitialNode, bool typeAttributeMustExist, bool mustDelete);

 
	//brief return proostList (use in GUI)
	map<QString, ObjectList<Base>* >* getList() {return mAllLists;};
	
	// GUI method
	void ReplaceInstance (QString aClassName, QString aOldName, QString aNewName)
    {
        ObjectList<Base>* curList;
        curList = FindInMap(mAllLists, aClassName, ERROR_CLASS_NAME_NOT_FOUND);

	    Base* inst = curList->Find(aOldName);
	    if ( inst == 0 )
            throw GeneralException(QString("Name not found") + " in list: \"" + aClassName + "\"" + " named: " + aOldName);
		inst->SetName(aNewName);
		curList->Replace(aOldName,aNewName);
    }
	
	

protected:
    /// \brief Add all objects of a single list, from an xml input file
    void FillListFromXML(ObjectList<Base>* aProostList, QDomElement aInitialNode, bool typeAttrMustExist, bool mustDelete);




protected:
    
    /// \brief Map of maps: Dictionary (indexed by class name) of instances (indexed by name)
	map<QString, ObjectList<Base>* >* mAllLists;

    /// \brief  creates an object factory
    void CreateNewFactory(void) {mFactory = new ObjectFactory<Base*(QDomElement),QString>;}


    /// \brief  creates a proost object of a given type. 
    Base* CreateNewObjectInstance(QString atrType, QDomElement aNode) {return(this->mFactory->Create(atrType,aNode));}

};


#endif
