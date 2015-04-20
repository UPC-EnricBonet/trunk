#ifndef GINIT_CHEPROO_H
#define GINIT_CHEPROO_H

#include "classnameconstants.h"
#include "xmlconstants.h"
#include "objectmanager.h"


#include "caqueous.h"
#include "cmineral.h"
#include "cmineralpure.h"
#include "caqueousdebyehuckel.h"
#include "caqueousideal.h"
#include "caqueouspitzer.h"
#include "ccheprooplusplus.h"
#include "cgas.h"
#include "cgasbinarymixture.h"
#include "cgasideal.h"
#include "cglobalchemicalsystem.h"
#include "clocalchemicalsystem.h"
#include "clocalchemicalsystemconstprimary.h"
#include "clocalchemicalsystemsaaltink.h"
#include "cchemicalcomposition.h"
#include "cphase.h"
#include "creaction.h"
#include "creactionratelaw.h"
#include "creactionratelawlasaga.h"
#include "creactionratelawmonod.h"
#include "cspecies.h"
#include "csurface.h"
#include "csurfacecexch.h"

void gInit_cheproo()
{
    ObjectManager* om = ObjectManager::instance();

    om->RegisterProostClass<CAqueous>(CLASS_NAME_AQUEOUS);
    om->RegisterProostClass<CAqueousDebyeHuckel>(CLASS_NAME_AQUEOUSDEBYEHUCKEL);
    om->RegisterProostClass<CAqueousIdeal>(CLASS_NAME_AQUEOUSIDEAL);
    om->RegisterProostClass<CAqueousIdeal>(CLASS_NAME_AQUEOUSPITZER);
    om->RegisterProostClass<CCheprooPlusPlus>(CLASS_NAME_CHEPROOPLUSPLUS);
    om->RegisterProostClass<CGas>(CLASS_NAME_GAS);
    om->RegisterProostClass<CGasBinaryMixture>(CLASS_NAME_GASBINARYMIXTURE);
    om->RegisterProostClass<CGasIdeal>(CLASS_NAME_GASIDEAL);
    om->RegisterProostClass<CGlobalChemicalSystem>(CLASS_NAME_GLOBALCHEMICALSYSTEM);
    om->RegisterProostClass<CLocalChemicalSystem>(CLASS_NAME_LOCALCHEMICALSYSTEM);
    om->RegisterProostClass<CLocalChemicalSystemConstPrimary>(CLASS_NAME_LOCALCHEMICALSYSTEMCONSTPRIMARY);
    om->RegisterProostClass<CLocalChemicalSystemSaaltink>(CLASS_NAME_LOCALCHEMICALSYSTEMSAALTINK);
    om->RegisterProostClass<CMineral>(CLASS_NAME_MINERAL);
    om->RegisterProostClass<CMineralPure>(CLASS_NAME_MINERALPURE);
    om->RegisterProostClass<CChemicalComposition>(CLASS_NAME_CHEMICALCOMPOSITION);
    om->RegisterProostClass<CReaction>(CLASS_NAME_REACTION);
    om->RegisterProostClass<CReactionRateLaw>(CLASS_NAME_REACTIONRATELAW);
    om->RegisterProostClass<CReactionRateLawLasaga>(CLASS_NAME_REACTIONRATELAWLASAGA);
    om->RegisterProostClass<CReactionRateLawMonod>(CLASS_NAME_REACTIONRATELAWMONOD);
    om->RegisterProostClass<CSpecies>(CLASS_NAME_SPECIES);
    om->RegisterProostClass<CSpecies>(CLASS_NAME_SURFACE);
    om->RegisterProostClass<CSpecies>(CLASS_NAME_SURFACECEXCHANGE);
    om->RegisterProostClass<CPhase>(CLASS_NAME_PHASE);

    // Proost lists (only Base Clases that have lists (e.g. "phases")).

    om->RegisterProostListBaseClass<CChemicalComposition>(CLASS_NAME_CHEMICALCOMPOSITION, XML_TAG_CHEMCOMP, XML_TAG_CHEMCOMP);

}

#endif