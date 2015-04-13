#ifndef OBJECTLIST_H
#define OBJECTLIST_H

#include <map>
#include <QString>
#include <vector>
using namespace std;

#include "generalexception.h"
#include "xmlqt.h"
#include "serialization.h"
#include "deserialization.h"


#define XML_TAG_PROOSTLIST "proostList"


/*
	usefull comments. 
*/
#define vectype vector<pair<keyType,valueType> >

template <typename keyType, typename valueType> class Vec2map_wrapper
{
private:
     vectype mInstances;

public:
     typename vectype::iterator begin();

     typename vectype::iterator  end();

     void erase(typename vectype::iterator first, typename vectype::iterator last );

     void insert(keyType t, valueType s);

	 void replace(keyType oldKey, valueType oldValue,keyType newKey, valueType newValue);

     typename vectype::iterator find(const keyType& tofind );





};

 template <typename keyType, typename valueType>  typename vectype::iterator Vec2map_wrapper<keyType, valueType>::begin()
{
    return mInstances.begin();
}
 template <typename keyType, typename valueType>  typename vectype::iterator Vec2map_wrapper<keyType, valueType>::end()
{
    return mInstances.end();
}

template <typename keyType, typename valueType>   void  Vec2map_wrapper<keyType, valueType>::erase(typename vectype::iterator first,typename vectype::iterator last )
{

    this->mInstances.erase(first, last);

}
template <typename keyType, typename valueType>   void   Vec2map_wrapper<keyType, valueType>::insert(keyType t, valueType s)
{
    pair<keyType, valueType> dummy;
    dummy =  std::make_pair(t, s);
    this->mInstances.push_back(dummy);
}

template <typename keyType, typename valueType>  typename vectype::iterator    Vec2map_wrapper<keyType, valueType>::find(const keyType& tofind )
{
    for( typename vectype::iterator it = this->mInstances.begin(); it != this->mInstances.end(); it++)
    {
        if((*it).first == tofind)
        {
            return it;
        }
    }
    return this->mInstances.end();
}

template<typename T> class ObjectList
{
public:

    /// \brief name of list
    QString mTagListName;
    
    /// \brief name of objects on list    
    QString mTagName;

    bool mMustRead;

    Vec2map_wrapper<QString, T*> mInstances;

public:
    ObjectList(){}

    ObjectList(const QString aTagListName, const QString aTagName, const bool aMustRead)
    {
        mTagListName = aTagListName; mTagName = aTagName; mMustRead = aMustRead;
    }

    void Add(QString aName, T* instance)
    {
        if ( instance != 0  )
        {
            if ( aName == "" ) 
            {
                instance->SetName(Base::GenerateIdentifier(instance->mClassName));

                aName = instance->name();
            }


            if ( mInstances.find(aName) == mInstances.end())
            {
                mInstances.insert( aName, instance  );
            }
            else
            {
                 QString msg = "element: " + aName + " already in list: " + this->mTagName;
                 cout << msg.toStdString();
                 throw(GeneralException(msg));
             }
        }

    }

    T* Find(QString aName)
    {
        typename  vector<pair<QString, T*> >::iterator iter = mInstances.find(aName.trimmed() );

        if ( iter == mInstances.end() ) return 0;

        return (*iter).second;
    }

	    void Replace(QString aOldName,QString aNewName)
    {
        typename  vector<pair<QString, T*> >::iterator iter = mInstances.find(aOldName.trimmed() );

        if ( iter == mInstances.end() ) return;

		(*iter).first=aNewName;
    }

    void Destroy(QString aName)
    {
        T* instance = Find(aName);
        if ( instance )
        {
            delete instance;
        }
        else
        {
            throw GeneralException("Cannot destroy, instance not found");
        }
    }



    T* Read(QDomElement aNode)
    {

        T* retVal = Find(aNode.text());
        if ( retVal != 0 ) throw GeneralException("Cannot read, instance not found");

        if ( retVal->mNodeFound != 0  &&  retVal->mIsAllRead == false ) 
        {
            retVal->Read();
        }

        return retVal;
    }

    void EraseList( )
    {
        typename  vector<pair<QString,T*> >::iterator it;
        for ( it = mInstances.begin(); it != mInstances.end() ; it++ )
        {
            T* elem = (*it).second;
            delete elem;
        }

        mInstances.erase( mInstances.begin() , mInstances.end() );
    }


    void ReadList()
    {
        typename  vector<pair<QString,T*> >::iterator it;
        for ( it = mInstances.begin(); it != mInstances.end() ; it++ )
        {
            T* elem = (*it).second;
            elem->Read();
        }
    }
	    vector<T*> iterator()
    {
		vector<T*> retVec;
        typename  vector<pair<QString,T*> >::iterator it;
        for ( it = mInstances.begin(); it != mInstances.end() ; it++ )
        {
			T* elem = (*it).second;
			retVec.push_back(elem);
        }
		return retVec;
    }
};
#endif
