#ifndef CCHEPROOBASE_H
#define CCHEPROOBASE_H

#include "base.h"

class CCheprooBase : public Base
{

public:

    /// \brief Constructor.
	CCheprooBase();

    /// \brief Constructor from a XML file.
	CCheprooBase(
		   QDomElement aNode ///< XML node with attributes
		   );
    
    virtual ~CCheprooBase();

};

#endif
