#ifndef GEND_H
#define GEND_H

#include "objectmanager.h"

void gEnd()
{
	ObjectManager* om = ObjectManager::instance();

	om->EraseProostListBaseClass();
	om->UnregisterAllProostClasses();

}
#endif