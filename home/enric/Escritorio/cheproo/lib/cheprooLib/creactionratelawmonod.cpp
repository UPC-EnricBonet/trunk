#include "creactionratelawmonod.h"


CReactionRateLawMonod::CReactionRateLawMonod()
{
    this->mClassName = "CReactionRateLawMonod";
}

CReactionRateLawMonod::CReactionRateLawMonod(QDomElement aNode)
: CReactionRateLaw(aNode)
{
    this->mClassName = "CReactionRateLawMonod";
}

CReactionRateLawMonod::~CReactionRateLawMonod()
{
}

void
CReactionRateLawMonod::Read(QDomElement aNode)
{
    CReactionRateLaw::Read(aNode);
}
