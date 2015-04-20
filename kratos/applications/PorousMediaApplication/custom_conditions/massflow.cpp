// Project includes 
#include "includes/define.h"
#include "custom_conditions/massflow.h"
#include "porous_media_application.h"
#include "utilities/math_utils.h"


namespace Kratos
{
    
        typedef GeometryData::KratosGeometryType KratosGeometryType;
        
	//************************************************************************************
	//************************************************************************************
	MassFlow::MassFlow(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	MassFlow::MassFlow(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
            
	}
	Condition::Pointer MassFlow::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new MassFlow(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	MassFlow::~MassFlow()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void MassFlow::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
            KRATOS_TRY
            switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
            {
                case 1:
               {
                   break;
                }
                case 2:
                {
                   this->CalculateRightHandSideConcentration(rRightHandSideVector, rCurrentProcessInfo);
                   break;
                }
                default:
                {
                    KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
                }
            }
            KRATOS_CATCH("")
           
        }

        void MassFlow::CalculateRightHandSideConcentration(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
            //calculation flags
            bool CalculateStiffnessMatrixFlag = false;
            bool CalculateResidualVectorFlag = true;
            MatrixType temp = Matrix();
            
            this->CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
        }
        
        void MassFlow::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
            KRATOS_TRY
                    
            switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
            {
                case 1:
                {
                    break;
                }
                case 2:
                {
                    this->CalculateLocalSystemConcentration(rLeftHandSideMatrix,rRightHandSideVector, rCurrentProcessInfo);
                    break;
                }
                default:
                {
                    KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
                }
            }
            
            KRATOS_CATCH("")
                    

        }
        
        void MassFlow::CalculateLocalSystemConcentration(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
        {
            //calculation flags
            bool CalculateStiffnessMatrixFlag = true;
            bool CalculateResidualVectorFlag = true;
            
            this->CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
        }
        
        void MassFlow::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
            switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
            {
                case 1:
                {
                    break;
                }
                case 2:
                {
                    this->EquationIdVectorConcentration(rResult, rCurrentProcessInfo);
                    break;
                }
                default:
                {
                    KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
                }
            }
	}
        
        void MassFlow::EquationIdVectorConcentration(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
                int number_of_nodes = GetGeometry().PointsNumber();
                
		unsigned int index;
		unsigned int dim = 1;
		rResult.resize(number_of_nodes*dim);
		for (int i=0;i<number_of_nodes;i++)
		{
			index = i*dim;
			rResult[index] = (GetGeometry()[i].GetDof(CONCENTRATION).EquationId());			
		}
        }

	//************************************************************************************
	//************************************************************************************
	void MassFlow::GetDofList(DofsVectorType& rConditionalDofList,ProcessInfo& rCurrentProcessInfo)
	{
   
            switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
            {
                case 1:
                {
                    break;
                }
                case 2:
                {
                    this->GetDofListConcentration(rConditionalDofList, rCurrentProcessInfo);
                    break;
                }
                default:
                {
                            KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
                }
            }
         
	}
        void MassFlow::GetDofListConcentration(DofsVectorType& rConditionalDofList,ProcessInfo& rCurrentProcessInfo)
	{
   
            unsigned int dim = 1;
            rConditionalDofList.resize(GetGeometry().size()*dim);
            unsigned int index;
            for (unsigned int i=0;i<GetGeometry().size();i++)
            {

                index = i*dim;
                rConditionalDofList[index] = (GetGeometry()[i].pGetDof(CONCENTRATION));
            }
         
	}
        
	//************************************************************************************
	//************************************************************************************
	void MassFlow::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag)
	{
		KRATOS_TRY

                //Set Matrix & vector size
                unsigned int number_of_nodes = GetGeometry().size();
                unsigned int MatSize=number_of_nodes;
            
                const double prescribedValue = GetProperties()[PRESCRIBED_VALUE];
                const double prescribedDensity = 1000.0 + 700.0*(prescribedValue);
                
                double inputFlow_i = 0.0;
                if(GetGeometry()[0].FastGetSolutionStepValue(SINKSOURCE_BALANCE) > 0.0)
                    inputFlow_i= GetGeometry()[0].FastGetSolutionStepValue(SINKSOURCE_BALANCE);
                
                if (CalculateStiffnessMatrixFlag==true)
                {
                    if(rLeftHandSideMatrix.size1() != MatSize)
                        rLeftHandSideMatrix.resize(MatSize,MatSize,false);
                    noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize);
                    
                    rLeftHandSideMatrix(0,0) = prescribedDensity*inputFlow_i;
                }
                
                if (CalculateResidualVectorFlag==true)
                {
                    if(rRightHandSideVector.size() != MatSize)
                        rRightHandSideVector.resize(MatSize,false);
                    
                    const double currentSVValue = GetGeometry()[0].FastGetSolutionStepValue(CONCENTRATION);
                    rRightHandSideVector[0] = prescribedDensity*inputFlow_i*(prescribedValue - currentSVValue);
                } 

                KRATOS_CATCH("");
	}
          
	//************************************************************************************
	//************************************************************************************
} // Namespace Kratos
