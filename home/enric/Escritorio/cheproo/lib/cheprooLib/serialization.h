#ifndef SERIALIZATION_H
#define SERIALIZATION_H

#include <QString>
#include <QTextStream>
#include <vector>
#include <map>
#include <valarray>
#include <set>

using namespace std;

class Serialization
{
private:
    QString mSerClassName;
    QTextStream* mOut;
    QString *mS;
public:
    Serialization(QString aClassName, QTextStream* aOut)
    {
        mSerClassName = aClassName; mOut = aOut;
        mS = mOut->string();
    }

    void Begin()
    {
        *mOut << "<" << mSerClassName << ">";
    }
    void End()
    {
        *mOut << "</" << mSerClassName << ">";
    }

    template<typename T>
    void SerializeAttribute(QString attrName, T iAttribute, bool addTag = true)
    {
        if(addTag)
        {
            *mOut << "<" << attrName << ">";
        }
        *mOut << iAttribute;
        if(addTag)
        {
            *mOut << "</" << attrName << ">";
        }
    }

    template<typename T>
    void SerializeAttribute(QString attrName, T* iAttribute, bool goDeep)
    {
        if(iAttribute == 0)
        {
            this->SerializeAttribute(attrName, "");
        }
        else
        {
            if(goDeep)
            {
                iAttribute->Serialize(*mS, attrName, goDeep);
            }
            else
            {
                iAttribute->SerializeName(*mS, attrName);
            }
        }
    }

    void SerializeAttribute(QString attrName, bool iAttribute, bool addTag = true)
    {
        if(addTag)
        {
            *mOut << "<" << attrName << ">";
        }
        
        if(iAttribute)
        {
            iAttribute = 1;
        }
        *mOut << iAttribute;

        if(addTag)
        {
            *mOut << "</" << attrName << ">";
        }
    }

    template<typename T>
    void SerializeAttribute(QString attrName, vector<T> iAttribute, bool addTag = true)
    {
        if(addTag)
        {
            *mOut << "<" << attrName << ">";
        }
        int size =  iAttribute.size();
        if(size>0)
        {
            int i;
            for(i=0; i < size-1 ; i++)
            {
                this->SerializeAttribute("", iAttribute.at(i), false);
                *mOut << " ";
            }
            *mOut << iAttribute.at(i);
        }
        if(addTag)
        {
            *mOut << "</" << attrName << ">";
        }
    }

    template<typename T>
    void SerializeAttribute(QString attrName, vector<vector<T> > iAttribute, bool addTag = true)
    {
        if(addTag)
        {
            *mOut << "<" << attrName << ">";
        }
        int size =  iAttribute.size();
        if(size > 0)
        {
            int i;
            for(i=0; i < size-1 ; i++)
            {
                this->SerializeAttribute("", iAttribute.at(i), false);
                *mOut << ";";
            }
            this->SerializeAttribute("", iAttribute.at(i), false);
        }
        if(addTag)
        {
            *mOut << "</" << attrName << ">";
        }
    }

    template<typename T>
    void SerializeAttribute(QString attrName, vector<vector<vector<T> > > iAttribute, bool addTag = true)
    {
        if(addTag)
        {
            *mOut << "<" << attrName << ">";
        }
        int size =  iAttribute.size();
        if(size > 0)
        {
            int i;
            for(i=0; i < size-1 ; i++)
            {
                this->SerializeAttribute("", iAttribute.at(i), false);
                *mOut << ":";
            }
            this->SerializeAttribute("", iAttribute.at(i), false);
        }
        if(addTag)
        {
            *mOut << "</" << attrName << ">";
        }
    }

    template<typename T>
    void SerializeAttribute(QString attrName, vector<vector<vector<vector<T> > > > iAttribute, bool addTag = true)
    {
        if(addTag)
        {
            *mOut << "<" << attrName << ">";
        }
        int size =  iAttribute.size();
        if(size > 0)
        {
            int i;
            for(i=0; i < size-1 ; i++)
            {
                this->SerializeAttribute("", iAttribute.at(i), false);
                *mOut << "_";
            }
            this->SerializeAttribute("", iAttribute.at(i), false);
        }
        if(addTag)
        {
            *mOut << "</" << attrName << ">";
        }
    }

    template<typename T>
    void SerializeAttribute(QString attrName, QString itemType, vector<T*> iAttribute, bool goDeep)
    {
        *mOut << "<" << attrName << ">";
        foreach(T* curItem, iAttribute)
        {
            this->SerializeAttribute(itemType, curItem, goDeep);
        }
        *mOut << "</" << attrName << ">";
    }

    template<typename T>
    void SerializeAttribute(QString attrName, map<QString,T* > aMap, bool goDeep)
    {
        *mOut << "<" << attrName << ">";
        
        typename map<QString, T* >::iterator it;
        for ( it = aMap.begin(); it != aMap.end() ; it++ )
        {
            *mOut << "<mapItem>";
            SerializeAttribute("first", it->first);

            T* curObject = it->second;
            if(curObject)
            {
                if ( goDeep ) 
                {
                    
                    curObject->Serialize(*mS, "second", true);
                }
                else
                {
                    curObject->SerializeName(*mS, "second");
                }
            }

            *mOut << "</mapItem>";
        }
        *mOut << "</" << attrName << ">";
    }

    template<typename T>
    void SerializeAttribute(QString attrName, multimap<QString,T* > aMap, bool goDeep)
    {
        *mOut << "<" << attrName << ">";
        
        typename multimap<QString, T* >::iterator it;
        for ( it = aMap.begin(); it != aMap.end() ; it++ )
        {
            *mOut << "<mapItem>";
            SerializeAttribute("first", it->first);

            T* curObject = it->second;
            if(curObject)
            {
                if ( goDeep ) 
                {
                    
                    curObject->Serialize(*mS, "second", true);
                }
                else
                {
                    curObject->SerializeName(*mS, "second");
                }
            }

            *mOut << "</mapItem>";
        }
        *mOut << "</" << attrName << ">";
    }

    template<typename T>
    void SerializeMap(QString attrName, map<QString,T> aMap)
    {
        *mOut << "<" << attrName << ">";
        typename map<QString,T>::iterator it;
        for ( it = aMap.begin(); it != aMap.end() ; it++ )
        {
            *mOut << "<mapItem>";
            SerializeAttribute("first", (*it).first);
            SerializeAttribute<T>("second", (*it).second);

            *mOut << "</mapItem>";
        }
        *mOut << "</" << attrName << ">";
    }

    template<typename T>
    void SerializeQt(QString attrName, T aQtVar)
    {
        *mOut << "<" << attrName << ">";
        QString myString;
        QDebug myDebugStream(&myString);
        myDebugStream << aQtVar;
        *mOut << myString;
        *mOut << "</" << attrName << ">";
    }

    template<typename T>
    void SerializeAttribute(QString attrName, valarray<T> iAttribute, bool addTag = true)
    {
        if(addTag)
        {
            *mOut << "<" << attrName << ">";
        }
        int size =  iAttribute.size();
        if(size>0)
        {
            int i;
            for(i=0; i < size-1 ; i++)
            {
                this->SerializeAttribute("", iAttribute[i], false);
                *mOut << " ";
            }
            *mOut << iAttribute[i];
        }
        if(addTag)
        {
            *mOut << "</" << attrName << ">";
        }
    }
    
    template<typename T>
    void SerializeAttribute(QString attrName, set<T> iAttribute,  bool goDeep)
    {
       *mOut << "<" << attrName << ">";

       int size =  iAttribute.size();

       if (goDeep)
       {
            for(typename set<T>::iterator i= iAttribute.begin(); i != iAttribute.end() ; i++)
            {
                //this->SerializeAttribute("", i, false);
            }
       }


        *mOut << "</" << attrName << ">";
        
    }

    template<typename T>
    void SerializeAttribute(QString attrName, valarray<valarray<T> > iAttribute, bool addTag = true)
    {
        if(addTag)
        {
            *mOut << "<" << attrName << ">";
        }
        int size =  iAttribute.size();
        if(size > 0)
        {
            int i;
            for(i=0; i < size-1 ; i++)
            {
                this->SerializeAttribute("", iAttribute[i], false);
                *mOut << ";";
            }
            this->SerializeAttribute("", iAttribute[i], false);
        }
        if(addTag)
        {
            *mOut << "</" << attrName << ">";
        }
    }

    template<typename T>
    void SerializeAttribute(QString attrName, valarray<valarray<valarray<T> > > iAttribute, bool addTag = true)
    {
        if(addTag)
        {
            *mOut << "<" << attrName << ">";
        }
        int size =  iAttribute.size();
        if(size > 0)
        {
            int i;
            for(i=0; i < size-1 ; i++)
            {
                this->SerializeAttribute("", iAttribute[i], false);
                *mOut << ":";
            }
            this->SerializeAttribute("", iAttribute[i], false);
        }
        if(addTag)
        {
            *mOut << "</" << attrName << ">";
        }
    }
};

#endif
