#ifndef GUIITEM_H
#define GUIITEM_H

#include <vector>
#include <QString>
//#include <QToolTip>

#include "generalexception.h"
#include "xmlqt.h"

/*!
 \brief Falta QToolTip.

COMMENT

*/
class GuiItem	
{
private:

	QString mClassName;

protected:

	QString mName;
	QString mValue;
	QString mText;
	QString mId;
	int mMin,mMax,mGui;
	QString mWidget;
	QString mGroupBox;
	//QToolTip mInfo;
	vector <GuiItem*>* mAttributes;
	vector <GuiItem*>* mRefs;
	vector <QString> mAdmitedValues;

	 /// iTemplate stuff
	 //QDomNode* iTemplate;
    
    ///	\brief Constructor
	

public:	
	
	GuiItem();

	void addAttributes(vector <GuiItem*>* aAttributes);
	void addAttribute (GuiItem* aAttribute);

	void addRefs(vector <GuiItem*>* aRefs);
	void addRef (GuiItem* aRef);

	GuiItem* createAttribute (QString name, QString value);
	//brief <tagName name="aName" type="aType"/>
	GuiItem* createNameNode (QString tagName, QString name, QString type);
	static GuiItem* createRefNode (QString type, QString text);
	//brief <tagName>text</tagName>
	GuiItem* createTextNode (QString tagName, QString text);

	//Read
	static GuiItem* XMLToGuiItem (QString aClassName,const QDomElement aParentNode, const QDomElement aNode, const QDomElement aSkeletonRoot);
	static void passAttributes (const QDomElement aNode,const QDomElement aNodeParent,const QDomElement aSkeletonRoot,GuiItem* aGuiItem);
	void XMLSkeletonToGuiItem (QString aSkeletonClassTagName, QString aSkeletonItemTagName, const QDomElement aSkeletonRoot);
	static QDomElement findXMLNode(QString aXMLClassTagName, const QDomElement aXMLRoot);

	//Write
	void GuiItemToXML (QDomDocument& doc, QDomElement& parent);
	//void fillNode (QDomNode* childNode, QDomElement& root);

	///////////////////////////////////////////////
	bool hasAttributes();
	bool hasChilds();
	bool hasText();
    static void guiNumberItems(GuiItem* aItem, vector <GuiItem*>* aTreeItems, int aGuiNumber);
	//////////////////////////////////////////////

	QString getValueFromAttribute(QString aName);
	GuiItem* getAttributeByName (QString aName);
	void getSetGuiItemByValue (const QString& aOldValue, const QString& aNewValue,const QString& aLabelName,const QString& aWidget);
	GuiItem* getGuiItemByNameAndValue (const QString& aName, const QString& aValue);
	//static void getGuiItemByValue (const QString& aValue,GuiItem* guiItemToParse, GuiItem* searchedGuiItem);
	void replaceCreateIDAttribute (const QString& aName, const QString& aOldValue, const QString& aNewValue);

	void replaceNameRefGuiItem (const QString& aName, const QString& aOldValue, const QString& aNewValue);
	void replaceNameGuiItem (const QString& aNewValue);

	QString getName()
	{
		return mName;
	}

	void setName(QString aName)
	{
		mName=aName;
	}

	QString getValue()
	{
		return mValue;
	}

	void setValue(QString aValue)
	{
			mValue=aValue;
	}

	QString getText()
	{
		return mText;
	}

	void setText(QString aText)
	{
		mText=aText;
	}

	QString getId()
	{
		return mId;
	}

	void setId(QString aId)
	{
		 mId=aId;
	}

	int  getMin()
	{
		return mMin;
	}

	void setMin(int aMin)
	{
		 mMin=aMin;
	}

	int  getMax()
	{
		return mMax;
	}

	void setMax(int aMax)
	{
		 mMax=aMax;
	}

	int  getGui()
	{
		return mGui;
	}

	void setGui(int aGui)
	{
		 mGui=aGui;
	}
	QString getWidget()
	{
		return mWidget;
	}

	void setWidget(QString aWidget)
	{
		mWidget=aWidget;
	}
	QString getGroupBox()
	{
		return mGroupBox;
	}

	void setGroupBox(QString aGroupBox)
	{
		mGroupBox=aGroupBox;
	}

	vector <GuiItem*>* getAttributes ()
	{
		return mAttributes;
	}

	void setAttributes (vector <GuiItem*>* aAttributes)
	{
		mAttributes=aAttributes;
	}

	vector <GuiItem*>* getRefs ()
	{
		return mRefs;
	}

	void setRefs (vector <GuiItem*>* aRefs)
	{
		mRefs=aRefs;
	}

	vector <QString> getAdmitedValues ()
	{
		return mAdmitedValues;
	}	

	QString getClassName ()
	{
		return this->mClassName;
	}

	void setClassName (QString aClassName)
	{
		mClassName=aClassName;
	}

	/*QDomNode* getiTemplate()
	{
		return this->iTemplate;
	}

	void setiTemplate (QDomNode* aDomNode)
	{
		iTemplate=aDomNode;
	}*/

	/*QDomNode getDomNode()
	{
		return this->mDomNode;
	}

	void setDomNode (QDomNode aDomNode)
	{
		mDomNode=aDomNode;
	}*/

};

#endif