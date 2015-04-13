#include "ccheproobase.h"

CCheprooBase::CCheprooBase()
{
    this->mClassName = "CCheprooBase";
}

CCheprooBase::CCheprooBase(QDomElement aNode) : Base(aNode)
{
    this->mClassName = "CCheprooBase";

   /* if ( this->name() == "")
    {
        this->PrThrow(ERROR_NAME_ATTRIBUTE_CANNOT_BE_VOID);
    } */
}

CCheprooBase::~CCheprooBase()
{
}

