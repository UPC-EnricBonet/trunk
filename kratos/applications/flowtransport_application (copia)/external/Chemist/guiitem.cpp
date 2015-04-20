#include "guiitem.h" 

#include <QFile>
#include <QTextStream>

#include "commonxmlconstants.h"
#include "serialization.h"
#include "deserialization.h"

GuiItem::GuiItem()
{
	this->mName="";
	this->mValue="";
	this->mText="";
	this->mGui =0;
	this->mGroupBox="";
	this->mAttributes = new vector<GuiItem*>;
	this->mRefs = new vector<GuiItem*>;
	//this->mDomNode=new QDomNode;
	//this->iTemplate=new QDomElement;
}

void 
GuiItem::addAttributes(vector <GuiItem*>* aAttributes)
{
	for(int i=aAttributes->size()-1;i>=0;i--)
	 {
		 this->mAttributes->push_back(aAttributes->at(i));
	 }
}

void 
GuiItem::addAttribute (GuiItem* aAttribute)
{
	if( this->mAttributes == 0 )
    {
		 this->mAttributes = new vector<GuiItem*>;
	}	 
	this->mAttributes->push_back(aAttribute);
}

void 
GuiItem::addRefs (vector <GuiItem*>* aRefs)
{
	 for(int i=aRefs->size()-1;i>=0;i--)
	 {
		 this->mRefs->push_back(aRefs->at(i));
	 }
}

void 
GuiItem::addRef (GuiItem* aRef)
{
	/*if( this->mRefs == 0 )
    {
		 this->mRefs = new vector<GuiItem*>;
	}	 */
	this->mRefs->push_back(aRef);

}

GuiItem* 
GuiItem::createAttribute (QString name, QString value)
{
	GuiItem* attribute=new GuiItem();
	attribute->setName(name);
	attribute->setValue(value);

	return attribute;
}

GuiItem* 
GuiItem::createRefNode (QString type, QString text)
{
		 GuiItem* refNode=new GuiItem();
		 refNode->setName("ref");
		 refNode->setText(text);

		 GuiItem* attributeType=new GuiItem();
		 attributeType->setName("type");
		 attributeType->setValue(type);

		 refNode->addAttribute(attributeType);

		return refNode;

}

GuiItem* 
GuiItem::createNameNode (QString tagName, QString name, QString type)
{
		 GuiItem* nameNode=new GuiItem();
		 nameNode->setName(tagName);

		 GuiItem* attributeName=new GuiItem();
		 attributeName->setName("name");
		 attributeName->setValue(name);

		 GuiItem* attributeType=new GuiItem();
		 attributeType->setName("type");
		 attributeType->setValue(type);

		 nameNode->addAttribute(attributeType);

		return nameNode;
}

GuiItem* 
GuiItem::createTextNode (QString tagName, QString text)
{
		 GuiItem* textNode=new GuiItem();
		 textNode->setName(tagName);
		 textNode->setValue(text);

		 return textNode;
}

GuiItem*
GuiItem::XMLToGuiItem (QString aClassName,const QDomElement aParentNode, const QDomElement aNode, const QDomElement aSkeletonRoot)
{
	GuiItem* guiItem=new GuiItem();
	guiItem->setClassName(aClassName.trimmed());

	QDomElement parentNode=aParentNode;
	QDomElement node=aNode;

	//Set Name & value
	guiItem->setName(node.nodeName().trimmed());
	if(!node.toElement().nodeValue().isEmpty())
		guiItem->setValue(node.nodeValue().trimmed());

	// Set Text, if it necesary
	if(node.firstChild().isText())
		guiItem->setText(node.toElement().text().trimmed());

	//Set Attributes
	if(node.hasAttributes())
		GuiItem::passAttributes(node,parentNode,aSkeletonRoot,guiItem);

	//Set skeleton
	QString skeletonClassTagName=aClassName;

	QString skeletonItemTagName;
	if(guiItem->getName()!="ref")
		skeletonItemTagName=parentNode.nodeName()+aNode.nodeName();
	else
		skeletonItemTagName= parentNode.nodeName()+aNode.nodeName()+guiItem->getAttributeByName("type")->getValue();

	guiItem->XMLSkeletonToGuiItem(skeletonClassTagName, skeletonItemTagName, aSkeletonRoot);

	//Set Childs
	if(aNode.hasChildNodes() && !aNode.childNodes().item(0).isText())//Avoid Text like node
	{
		for (unsigned i=0;i< aNode.childNodes().length();i++)
		{
			GuiItem* guiItemChild=new GuiItem();

			guiItemChild = XMLToGuiItem (aClassName,node,node.childNodes().item(i).toElement(),aSkeletonRoot);

			guiItem->addRef(guiItemChild);
		}
	}


	return guiItem;
}

QDomElement 
GuiItem::findXMLNode(QString aXMLClassTagName, const QDomElement aXMLRoot)
{
	QDomElement classNode;

	for (unsigned i =0; i< aXMLRoot.childNodes().length(); i++)
		if(aXMLRoot.childNodes().item(i).nodeName()==aXMLClassTagName)
		{
			classNode=aXMLRoot.childNodes().item(i).toElement();
			break;
		}

	return classNode;
}

void
GuiItem::passAttributes  (const QDomElement aNode,const QDomElement aNodeParent,const QDomElement aSkeletonRoot,GuiItem* aGuiItem)
{

		QDomNamedNodeMap attributesMap;
		attributesMap=aNode.attributes();
		int jj=attributesMap.length();

		for(unsigned i=0;i < attributesMap.length();i++)
		{
			GuiItem* attribute=new GuiItem();
			attribute->setClassName(aGuiItem->getClassName().trimmed());
			attribute->setName(attributesMap.item(i).toAttr().name().trimmed());
			attribute->setValue(attributesMap.item(i).toAttr().value().trimmed());
			attribute->XMLSkeletonToGuiItem(aGuiItem->getClassName(), aNode.nodeName()+attributesMap.item(i).toAttr().name(),aSkeletonRoot);
			aGuiItem->addAttribute(attribute);
		}	
}

void 
GuiItem::XMLSkeletonToGuiItem (QString aSkeletonClassTagName, QString aSkeletonItemTagName, const QDomElement aSkeletonRoot)
{

	QDomElement classNode=GuiItem::findXMLNode(aSkeletonClassTagName,aSkeletonRoot);

	QDomNode skeletonNode = classNode.toElement().elementsByTagName(aSkeletonItemTagName).at(0);

	this->setMin (skeletonNode.attributes().namedItem("min").nodeValue().trimmed().toInt());
	this->setMax (skeletonNode.attributes().namedItem("max").nodeValue().trimmed().toInt());
	this->setGui (skeletonNode.attributes().namedItem("gui").nodeValue().trimmed().toInt());
	this->setWidget (skeletonNode.attributes().namedItem("widget").nodeValue().trimmed());
	this->setGroupBox (skeletonNode.attributes().namedItem("groupBox").nodeValue().trimmed());
}


void
GuiItem::guiNumberItems(GuiItem* aItem, vector <GuiItem*>* aTreeItems, int aGuiNumber)
{	
	if( aTreeItems == 0 )
    {
		 aTreeItems = new vector<GuiItem*>;
	}

	if(aItem->getGui()!=0 && aItem->getGui() == aGuiNumber)
		aTreeItems->push_back(aItem);

	//Rastrear atributos también y meterlos en los gui items IMPORTANTE
	if(aItem->hasAttributes())
	{
	vector <GuiItem*>* attributes = new vector<GuiItem*>;
	attributes=aItem->getAttributes();
		for(unsigned i=0;i<attributes->size();i++)
			if(attributes->at(i)->getGui()!=0 && attributes->at(i)->getGui() == aGuiNumber)
				aTreeItems->push_back(attributes->at(i));	
	}

	if(aItem->hasChilds())
	{
		vector <GuiItem*>* childs = aItem->getRefs();
		for(unsigned i=0;i<childs->size();i++)
			GuiItem::guiNumberItems(childs->at(i), aTreeItems,aGuiNumber);	
	}
}

void
GuiItem::GuiItemToXML (QDomDocument& doc, QDomElement& parent)
{
	//Set Atributes rootNode
	if(hasAttributes())
		for(unsigned i=0;i<mAttributes->size();i++)
		{
			GuiItem* guiItemAttr = new GuiItem;
			guiItemAttr = getAttributes()->at(i);

			parent.setAttribute(guiItemAttr->getName(),guiItemAttr->getValue());
		}
 
	//Set TextNode rootNode
	if(hasText())
	{

		QDomText childDomNode =doc.createTextNode(mText);
		parent.appendChild(childDomNode);
	}

	//Set childsNode rootNode
	if(hasChilds())
		for(unsigned i=0;i<mRefs->size();i++)
		{
			GuiItem* guiItemChild = new GuiItem;
			guiItemChild = getRefs()->at(i);

			QDomElement childDomNode =doc.createElement(guiItemChild->getName());
			parent.appendChild(childDomNode);
			guiItemChild->GuiItemToXML(doc,childDomNode);
		}
}

bool 
GuiItem::hasChilds()
{
	if(this->getRefs()->empty())// childs->empty con "#text" vector size ????
		return false;
	else
		return true;
}

bool 
GuiItem::hasAttributes()
{

	if(this->getAttributes()->empty())
		return false;
	else
		return true;
}

bool 
GuiItem::hasText()
{
	if(this->getText().isEmpty())
		return false;
	else
		return true;
}

GuiItem* 
GuiItem::getAttributeByName (QString aName)
	{
		 for (int i = 0; i < mAttributes->size(); ++i) {
			 if (this->mAttributes->at(i)->getName() == aName)
				 return mAttributes->at(i);
		 }
	}

QString
GuiItem::getValueFromAttribute(QString aName)
{
	GuiItem* item = this->getAttributeByName(aName);
	return item->getValue();
}

void
GuiItem::getSetGuiItemByValue (const QString& aOldValue, const QString& aNewValue,const QString& aLabelName,const QString& aWidget)
{
	if(this->getName()!="ref" && this->getName() == aLabelName && this->getWidget() == aWidget && this->getValue() == aOldValue)
	{
			this->setValue(aNewValue);
			return;
	}
	if(this->getName()!="ref" && this->getName() == aLabelName && this->getWidget() == aWidget && this->getText() == aOldValue)
	{
			this->setText(aNewValue);
			return;
	}
	else if(this->getName()=="ref" && this->getAttributeByName("type")->getValue() == aLabelName && this->getWidget() == aWidget && this->getText() == aOldValue)
	{
		this->setText(aNewValue);
		return;
	}
		
	for (unsigned i=0;i<getAttributes()->size();i++)
		if(getAttributes()->at(i)->getName() == aLabelName && getAttributes()->at(i)->getWidget() == aWidget && getAttributes()->at(i)->getValue()==aOldValue)
			{
				this->getAttributes()->at(i)->setValue(aNewValue);
				return;
			}

	for(unsigned i=0;i<getRefs()->size();i++)
		this->getRefs()->at(i)->getSetGuiItemByValue(aOldValue,aNewValue,aLabelName,aWidget);

}

void
GuiItem::replaceCreateIDAttribute (const QString& aName, const QString& aOldValue, const QString& aNewValue)
{
	if(!aOldValue.isEmpty())
	{
		std::vector <GuiItem*>* refVec= this->getRefs();
		for(unsigned i=0;i<refVec->size();i++)
			if(refVec->at(i)->getText()==aName)
				refVec->at(i)->getAttributeByName("ID")->setValue(aNewValue);
	}
	else
	{ 
		GuiItem *IDattribute=this->createAttribute("ID",aNewValue);
			std::vector <GuiItem*>* refVec= this->getRefs();
			for(unsigned i=0;i<refVec->size();i++)
				if(refVec->at(i)->getText()==aName)
					refVec->at(i)->addAttribute(IDattribute);
	}
}

void
GuiItem::replaceNameRefGuiItem (const QString& aName, const QString& aOldValue, const QString& aNewValue)
{
		std::vector <GuiItem*>* refVec= this->getRefs();
		for(unsigned i=0;i<refVec->size();i++)
			if(refVec->at(i)->getText()==aName)
				refVec->at(i)->setText(aNewValue);
}

void
GuiItem::replaceNameGuiItem (const QString& aNewValue)
{
		this->getAttributeByName("name")->setValue(aNewValue);
}


GuiItem*
GuiItem::getGuiItemByNameAndValue (const QString& aName, const QString& aValue)
{
	if(this->mName==aName && this->mValue==aValue || this->mText==aValue)
		return this;
	
	for(unsigned i=0;i<this->mAttributes->size();i++)
		if(this->mAttributes->at(i)->mName==aName && this->mAttributes->at(i)->mValue==aValue)
			return  this->mAttributes->at(i);

	for(unsigned i=0;i<this->mRefs->size();i++)
		if(this->mRefs->at(i)->mName==aName && this->mRefs->at(i)->mValue==aValue || this->mRefs->at(i)->mText==aValue)
			 this->mRefs->at(i)->getGuiItemByNameAndValue(aName,aValue);

}


// Otra opción, as static function
/*void 
GuiItem::getGuiItemByValue (const QString& aValue,GuiItem* guiItemToParse, GuiItem* searchedGuiItem)
{
	if(guiItemToParse->getValue() == aValue)
		{
			searchedGuiItem=guiItemToParse;
			return;
		}	

	for (unsigned i=0;i<guiItemToParse->getAttributes()->size();i++)
		if(guiItemToParse->getAttributes()->at(i)->getValue()==aValue)
			{
				searchedGuiItem= guiItemToParse->getAttributes()->at(i);
				return;
			}
		else
		continue;

	for(unsigned i=0;i<guiItemToParse->getRefs()->size();i++)
		if(guiItemToParse->getRefs()->at(i)->getValue()==aValue)
				{
					searchedGuiItem= guiItemToParse->getRefs()->at(i);
					return;
				}

		else
			GuiItem::getGuiItemByValue(aValue,guiItemToParse->getRefs()->at(i),searchedGuiItem);

}*/


/*void
GuiItem::fillNode(QDomNode* childNode, QDomElement& root)
{
	
	//Set TagName rootNode
	//childNode->toElement().setTagName(getName());
	 childNode->setNodeValue(getName());
	 QString dddDdd=childNode->toElement().nodeName();
	 QString dddDddd=childNode->toElement().nodeValue();
	 QDomDocument dddDdssdd=childNode->ownerDocument();QString hhhhh=dddDdssdd.firstChild().toElement().tagName();

	//Set Atributes rootNode
	if(hasAttributes())
		for(unsigned i=0;i<getAttributes()->size();i++)
		{
			GuiItem*curGuiItem = getAttributes()->at(i);

			QDomAttr attribute=childNode->toElement().attributeNode(curGuiItem->getName());
			attribute.setValue(curGuiItem->getValue());
		}

	//Set TextNode rootNode
	if(hasText())
	{
		QDomText* text=new QDomText;
		text->setNodeValue(getText());QString hhhh=this->getText();QString jj=childNode->toElement().tagName();
		childNode->appendChild(text->toElement());
	}

	if(hasChilds())
		for(unsigned i=0;i<getRefs()->size();i++)
		{
			QDomNode* curNode=new QDomNode;
			getRefs()->at(i)->fillNode(curNode,root);
			childNode->appendChild(curNode->toElement());
		}


}*/
