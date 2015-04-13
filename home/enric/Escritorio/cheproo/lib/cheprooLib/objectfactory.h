#ifndef OBJECTFACTORY_H
#define OBJECTFACTORY_H

#include <map>

template <typename CtorSignature, typename UniqueIdType> class ObjectFactory;

/// \brief creates an instance of a specialization class and returns it as an instance of a baseclass
template<typename BaseClassType, typename ClassType>
BaseClassType CreateObject()
{
   return new ClassType();
}


/// \brief the object factory contains a map relating class type (specified by some constant such as the class name)
/// and a pointer to the constructor of that class. It can therefore be used to create instances by giving it the desired type of object. 
// before putting it to use we must register the classes, giving each a name (key) and a pointer to its constructor as value. 

template<typename BaseClassType, typename UniqueIdType>
class ObjectFactory<BaseClassType (), UniqueIdType>
{
protected:
    //  the string "BaseClassType" is now called (*CreateObjectFunc)();
    typedef BaseClassType (*CreateObjectFunc)();

public:
   typedef typename std::map<UniqueIdType, CreateObjectFunc>::const_iterator ConstIterator;
   typedef typename std::map<UniqueIdType, CreateObjectFunc>::iterator Iterator;


/// \brief stores a function pointer to constructor in its  map (class name) to (constructor function)
   template<typename ClassType>
   bool Register(UniqueIdType unique_id)
   {
      if (m_object_creator.find(unique_id) != m_object_creator.end())
               return false;

      m_object_creator[unique_id] = &CreateObject<BaseClassType, ClassType>();

      return true;
   }

/// \brief removes a function pointer to constructor from its internal map of (class name) to (constructor function)
   bool Unregister(UniqueIdType unique_id)
   {
      return (m_object_creator.erase(unique_id) == 1);
   }

/// \brief unregisters all
   void UnregisterAll(void)
   {
      m_object_creator.erase( m_object_creator.begin(), m_object_creator.end() );
   }


   /// \brief Creates 
   BaseClassType Create(UniqueIdType unique_id)
   {
      Iterator iter = m_object_creator.find(unique_id);

      if (iter == m_object_creator.end())
         return NULL;

      return ((*iter).second)();
   }

   /// \brief returns const iterator to begin of  m_object_creator
   ConstIterator GetBegin() const
   {
      return m_object_creator.begin();
   }

   /// \brief returns  iterator to begin of  m_object_creator
   Iterator GetBegin()
   {
      return m_object_creator.begin();
   }


   /// \brief returns const  iterator to end of  m_object_creator
   ConstIterator GetEnd() const
   {
      return m_object_creator.end();
   }



   /// \brief returns iterator to end of  m_object_creator
   Iterator GetEnd()
   {
      return m_object_creator.end();
   }

protected:

   /// \brief map of (class name) to (pointer to construcotr of that class)
   std::map<UniqueIdType, CreateObjectFunc> m_object_creator;
};

template<typename BaseClassType, typename Param1Type, typename ClassType>
BaseClassType CreateObject(Param1Type param1)
{
   return new ClassType(param1);
}



/// \brief this object factory class does the same as the previous one, but it 
/// passes on a parameter to the constructor of any generated object. 
template<typename BaseClassType, typename Param1Type, typename UniqueIdType>
class ObjectFactory<BaseClassType (Param1Type),  UniqueIdType>
{
protected:
   typedef BaseClassType (*CreateObjectFunc)(Param1Type);

public:
   typedef typename std::map<UniqueIdType, CreateObjectFunc>::const_iterator ConstIterator;
   typedef typename std::map<UniqueIdType, CreateObjectFunc>::iterator Iterator;


   template<typename ClassType>
   bool Register(UniqueIdType unique_id)
   {
      if (m_object_creator.find(unique_id) != m_object_creator.end())
               return false;

      m_object_creator[unique_id] = &CreateObject<BaseClassType, Param1Type, ClassType>;

      return true;
   }

   bool Unregister(UniqueIdType unique_id)
   {
      return (m_object_creator.erase(unique_id) == 1);
   }

   void UnregisterAll(void)
   {
       if(m_object_creator.size() > 0)
       {
            m_object_creator.erase( m_object_creator.begin(), m_object_creator.end() );
       }
   }

   BaseClassType Create(UniqueIdType unique_id, Param1Type param1)
   {
      Iterator iter = m_object_creator.find(unique_id);

      if (iter == m_object_creator.end())
         return NULL;

      return ((*iter).second)(param1);
   }

   ConstIterator GetBegin() const
   {
      return m_object_creator.begin();
   }

   Iterator GetBegin()
   {
      return m_object_creator.begin();
   }

   ConstIterator GetEnd() const
   {
      return m_object_creator.end();
   }

   Iterator GetEnd()
   {
      return m_object_creator.end();
   }

protected:
   std::map<UniqueIdType, CreateObjectFunc> m_object_creator;
};


#endif
