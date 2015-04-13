#ifndef XMLVECTOR
#define XMLVECTOR

#include "xmlqt.h"
#include <QString>

static const QString XML_TAG_CONSTANT = "constant";
static const QString XML_TAG_LINEAR = "linear";
static const QString XML_ATTR_NRELEM = "nrelem";
static const QString XML_ATTR_VALUE = "value";
static const QString XML_ATTR_STARTVALUE = "startvalue";
static const QString XML_ATTR_INCREMENT = "increment";
static const QString XML_TAG_VECTOR = "XmlVector";

/* \brief 
This static class allows providing numerical data in a more convenient format in the xml file.
Currently one type of input vector is supported: the constant vector. 
instead of
<SomeTag> 1 1 1 1 1 1 1 1 1 1 1 </SomeTag>

we can give 
<SomeTag><XmlVector nrElem = "11" > <Constant value = "4" /></XmlVector></SomeTag>


*/
template <typename T> class XmlVector
{
public:

	static  vector<T>  GetVector(QDomElement aNode);
};
template <typename T>  vector<T> XmlVector<T>::GetVector(QDomElement aNode)
{
	vector<T> result;

	// get nr elements 
	int nrelem = 0;
	XmlQt::getXMLAttribute(aNode, XML_ATTR_NRELEM, nrelem);

	
	// get fill
	QDomElement XConstant = XmlQt::getFirstElemByTag(aNode,XML_TAG_CONSTANT);
	if (!XConstant.isNull())
	{
		T value = 0;
		XmlQt::getXMLAttribute(XConstant, XML_ATTR_VALUE, value);

		if (nrelem > 0)
		{
			result.resize(nrelem,value);
		}
        return result;
	}

	// get fill
	QDomElement XLinear = XmlQt::getFirstElemByTag(aNode,XML_TAG_LINEAR);
	if (!XLinear.isNull())
    {
		T startvalue = 0;
		XmlQt::getXMLAttribute(XLinear, XML_ATTR_STARTVALUE, startvalue);

		T increment = 0;
		XmlQt::getXMLAttribute(XLinear, XML_ATTR_INCREMENT, increment);

        result.resize(nrelem);
		for ( int ielem = 0; ielem < nrelem; ielem++)
		{
			result[ielem] = startvalue + ielem *increment;
		}
        return result;
        
    }

    

	//see if we have a content string
	QString vectext = aNode.text();
	if (!vectext.isEmpty())
	{
		//NOU !!!!?????????? MIRAR
		// description of real numbers
		QString pattern =  "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?";
		QRegExp rx = QRegExp(pattern);


		 int count = 0;
		 int pos = 0;
		 int endpos = 0;
		 QString num;

		 while ((pos = rx.indexIn(vectext, pos)) != -1)
		 {
     
			num = vectext.mid(pos, + rx.matchedLength());
			double interm = num.toDouble();
			result.push_back( (T) interm);
			pos +=  rx.matchedLength();

		}
		 return result;
	}

	throw("Error in XMlVector ");
    return result;	
}
#endif
