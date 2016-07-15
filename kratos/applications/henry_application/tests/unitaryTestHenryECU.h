/* *********************************************************
 *
 *   Last Modified by:    $Author: VictorBez $
 *   Date:                $Date: 2015        $
 *   Revision:            $Revision: 1.0     $
 *
 * ***********************************************************/


#if !defined(KRATOS_UNITARY_TEST_UTILITIES)
#define  KRATOS_UNITARY_TEST_UTILITIES


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"
#include <string>
#include <stdio.h>
//Inheritance
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_elements/flowpressuretrans_2d.h"



//Application
#include "henry_application.h"

using namespace std;

namespace Kratos
{
    
class UnitaryTestHenryECU
{
public:
    
   
    //Doc-> constructor
    UnitaryTestHenryECU(ModelPart& model_part,unsigned int aDomainSize = 2):mr_model_part(model_part)
    {
    }

    ~UnitaryTestHenryECU(){}
 


    void UnitaryTest()
    {
         //get the "model element" named "BaseTemperatureElement"
         //Element const& r_reference_element = KratosComponents<Element>::Get("FlowPressureTrans2D"); 
         
         //loop to make the copy
         /*for(ElementsContainerType::iterator it = (this->model_part).ElementsBegin(); it!=(this->model_part).ElementsEnd(); it++)
         {
               //Geometry< Node<3> >::Pointer pgeom = it->pGetGeometry(); //get the connectivity of the origin element
               //unsigned int new_id = it->Id(); //we want to create a new element with the same Id as the old one
               //PropertiesType::Pointer pprop = it->GetProperties(); //get the properties (of the old model part)
               //create a copy using the "reference element" as a model
              // Element::Pointer p_element = rReferenceElement.Create(new_id, pgeom ,properties);
               //fill the new model part 
               //new_model_part.Elements().push_back(p_element);
               it->CalculateDensityElement();
         }*/

         /*ModelPart::ElementsContainerType& rElements = mr_model_part.Elements();
         for (unsigned int rElements; rElements < mr_model_part.Elements().size(); rElements++)
         {
              rElements.CalculateDensityElement();
         }*/
         //ModelPart::ElementsContainerType& rElements = mr_model_part.Elements();
         for (int i=0; i < mr_model_part.Elements().size(); i++)
        {
            ModelPart::ElementsContainerType::iterator iel = mr_model_part.ElementsBegin()+i;
            //ModelPart::ElementsContainerType::iterator itElem = GetModelPart().ElementsBegin()+i;
            FlowPressureTrans2D* i_floweleme = dynamic_cast<FlowPressureTrans2D*>(&*iel);
            i_floweleme->CalculateDensityElement();
            //iel->CalculateUnitaryFunctions();//obtenemos valores calculados en cada una de las fucniones que hacemos testeo unitario
            //KRATOS_WATCH(elem);
        }
         //KRATOS_WATCH("hola_UnitaryTest()")
    }


    
protected:
    


private:    

    unsigned int mDomainSize;
    ModelPart& mr_model_part;

}; /* Class UnitaryTestUtilities */

} /* namespace Kratos.*/

#endif /* KRATOS_UNITARY_TEST_UTILITIES  defined */

