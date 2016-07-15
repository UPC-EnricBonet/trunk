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
    UnitaryTestHenryECU(ModelPart& model_part,
    unsigned int aDomainSize = 2)
    {
         this->mDomainSize = aDomainSize;
        ProcessInfo& CurrentProcessInfo = model_part.GetProcessInfo();
    }

    ~UnitaryTestHenryECU(){}
 


    void UnitaryTest()
    {
         KRATOS_WATCH("hola")
    }


    
protected:
    


private:    

    unsigned int mDomainSize;

}; /* Class UnitaryTestUtilities */

} /* namespace Kratos.*/

#endif /* KRATOS_UNITARY_TEST_UTILITIES  defined */

