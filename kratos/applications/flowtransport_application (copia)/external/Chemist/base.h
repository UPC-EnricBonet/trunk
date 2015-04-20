#ifndef BASE_H
#define BASE_H

#define gTimeStamp ObjectManager::gTimeStemp
#define gTimeTolerance ObjectManager::mTimeTolerance
#define gSpaceTolerance ObjectManager::mSpaceTolerance

#include <QString>
#include <vector>

#include "commonxmlconstants.h"
#include "generalexception.h"
#include "xmlqt.h"
#include "guiitem.h"

static const char* ERROR_NAME_ATTRIBUTE_CANNOT_BE_VOID = QT_TRANSLATE_NOOP("Base", "The name cannot be empty");
static const char* ERROR_NAME_ATTRIBUTE_NOT_FOUND = QT_TRANSLATE_NOOP("Base", "The XML element does not have the attribute 'name'");
static const char* ERROR_TYPE_ATTRIBUTE_NOT_FOUND = QT_TRANSLATE_NOOP("Base", "The XML element does not have the attribute 'type'");
static const char* ERROR_REF_HAS_NO_TYPE = QT_TRANSLATE_NOOP("Base", "The XML ref element does not have a 'type'");
static const char* ERROR_UNINITIALIZED_XML_NODE = QT_TRANSLATE_NOOP("Base", "The XML node has not been initialized");
static const char* ERROR_XML_NODE_CANNOT_BE_NULL = QT_TRANSLATE_NOOP("Base", "The XML node cannot be null");
static const char* ERROR_DERIVATIVES_CONSTANT_FIELD = QT_TRANSLATE_NOOP("Base", "A constant field must not calculate any kind of derivatives");
static const char* ERROR_THIS_ISNOT_IMPLEMENT_TO_THIS_CLASS_GUI = QT_TRANSLATE_NOOP("Base", "This functions is not implement to this class (GUI)");
static const char* ERROR_THIS_PROOSTOBJECTNAME_ALREADY_EXISTS = QT_TRANSLATE_NOOP("Base", "Name already exists (GUI)");
/*!
 \brief Virtual class from which most of Proost classes derive.

Usually an "Object Factory" is called to create any instance of a derived class (registered) on it. 

For each instance it retains its identifier (name), optionally the XML node in which it is located (m_nodeFound),
and a value indicating whether it has been read completely (mIsAllRead). 

*/
class Base
{

protected:

    QString mName;							///< \brief Instance name. Usually a unique identifier used when the derived class belongs to a Proost list (see proostList \<T\>)

	QDomElement mNodeFound;                 ///< \brief Node in the input XML file where the instance. If zero indicates that it has been created dynamically within the code

    bool mIsAllRead;                        ///< \brief Indicates whether it has read the entire instance attributes (or just the name and node).

	GuiItem* mGuiItem;						///< \brief Storage all information from GUI GUI

	QDomElement mNodeFoundParent;			///< \brief Parent Node in the input XML file where the instance  GUI



public:



    QString mClassName;

	/// \brief Instance name. Usually a unique identifier used when the derived class belongs to a Proost list (see proostList \<T\>)
	QString name() { return this->mName; } 

    void SetName(QString newName) { this->mName = newName; } 

    /// \brief Indicates whether it has read the entire instance attributes (or just the name and node).
	bool isAllRead() { return this->mIsAllRead; }

    /// \brief Indicates whether it has read the entire instance attributes (or just the name and node).
	void  HasRead(bool isRead) {  this->mIsAllRead = isRead; }
    

    /// \brief fgenerates an identifier (eg, for name or ID) that is unique,it starts with root and then it has the timestamp
    static QString GenerateIdentifier(const QString& root);
    
    ///	\brief Constructor
	Base();

    ///	\brief Destructor
    virtual ~Base(){}
    

    ///	\brief Copy Constructor
    Base(Base* aProost);


    ///	\brief Constructor from an XML node
	Base(QDomElement aNode               ///< XML node with CProost attributes, i.e. "name" and "type".
            );

    GuiItem* getGuiItem();			 ///< \brief get GuiItem GUI

	void setGuiItem(GuiItem* aGuiItem);	 ///< \brief set GuiItem GUI

	QString className() { return this->mClassName; }


    /// \brief Copies a CProost to another
    void Copy(Base &aTarget, Base &aSource);


    ///	\brief Reads attributes. Must be used after having set XML node attribute by its constructor.
	void Read();

    ///	\brief Reads attributes from an XML node.
	void Read(const QDomElement aNode       ///< XML node with CProost attributes, i.e. "name" and "type".
              );

    ///	\brief Reads attributes from an XML node.
	virtual void pRead(
					const QDomElement aNode ///< XML node with CProost attributes, i.e. "name" and "type".
                    );

    ///	\brief filters the child nodes of aNode. Returns those that ar refs to objects of a certain type.
    virtual void ReadReferences(const QDomElement aNode,
                                QString refType, 
                                vector<QDomElement> &refTags,
                                vector<QString>* aListOfIDs = 0);

    ///	\brief filters the child nodes of aNode. Returns the first ref to objects of a certain type.
    virtual void ReadFirstReferences(const QDomElement aNode,
                                QString refType, 
                                QDomElement &refTag,
                                QString* atrID = 0);

    ///	\brief filters the child nodes of aNode. Returns the object pointed to by 
    ///  the first reference of a certain type.
    template<typename OMType, typename ObjectToRetrieveType>
    void ReadFirstRefObject(const QDomElement aNode,
                          QString refType,
                          QString className,
                          const char* castingError,
                          ObjectToRetrieveType** objectToRetrieve,
                          const char* notFoundError = 0)
    {
        OMType* om = OMType::instance();
        QDomElement objectNode;
        this->ReadFirstReferences(aNode,
                                  refType, 
                                  objectNode);
        if(objectNode.isNull())
        {
            if(notFoundError == 0)
            {
                // If no error message is provided, just set the pointer 
                // to the object to retrieve  to 0.
                *objectToRetrieve = 0;
            }
            else
            {
                this->PrThrow(notFoundError);
            }
        }
        else
        {
            QString objectName = objectNode.text().trimmed();
            *objectToRetrieve = dynamic_cast<ObjectToRetrieveType*>(om->FindInstance(className, objectName));
            if(*objectToRetrieve == 0) 
            {
                this->PrThrow(castingError);
            }
        }
    }

    ///	\brief filters the child nodes of aNode. Returns the object pointed to by 
    ///  the first reference of a certain type WITH a given ID
    template<typename OMType, typename ObjectToRetrieveType>
    void ReadFirstRefObject(const QDomElement aNode,
                          QString refType,
                          QString className,
                          QString IDvalue,
                          const char* castingError,
                          ObjectToRetrieveType** objectToRetrieve,
                          const char* notFoundError = 0)
    {
        OMType* om = OMType::instance();
        QDomElement objectNode;
        
        // get references
        vector<QDomElement> refs = XmlQt::getElementsByTagName(aNode, XML_TAG_REF);
        
        // select those with desired class and ID
        vector<QDomElement> selectedRefs;
        for (int i = 0; i < refs.size(); i++)
        {
            // get class name
            QString refclass ;
            XmlQt::getXMLAttribute(refs[i], XML_ATTR_TYPE,refclass);
    
            // get ID name
            QString refID ;
            XmlQt::getXMLAttribute(refs[i], XML_ATTR_ID,refID);
    
            if (refclass == className && refID == IDvalue) selectedRefs.push_back(refs[i]);
    
         }
    
        // see if we have exactly one hit
        if (selectedRefs.size() != 1);
        {
            if(notFoundError == 0)
            {
                // If no error message is provided, just set the pointer 
                // to the object to retrieve  to 0.
                *objectToRetrieve = 0;
            }
            else
            {
                this->PrThrow(notFoundError);
            }
        }


        objectNode = selectedRefs[0];
        QString objectName = objectNode.text().trimmed();
        *objectToRetrieve = dynamic_cast<ObjectToRetrieveType*>(om->FindInstance(className, objectName));
        if(*objectToRetrieve == 0) 
        {
            this->PrThrow(castingError);
        }
    }

    virtual void PrThrow(const char* message);
    
    virtual void PrThrow(GeneralException &aException);

	// \brief add ref to proost object
	virtual bool addRef(Base*Obj);

	/// \brief set XMLNode object GUI
	void setNodeFound(QDomElement aNodeFound) {this->mNodeFound=aNodeFound; }
	/// \brief return parentXMLNode object GUI
	QDomNode parentNodeFound() {return this->mNodeFoundParent; }
	/// \brief set parentXMLNode object GUI
	void setParentNodeFound(QDomElement aParentNodeFound) {this->mNodeFoundParent=  aParentNodeFound; }

	/// \brief return XMLNode object GUI
	QDomNode nodeFound() {return this->mNodeFound; }



};
#endif
