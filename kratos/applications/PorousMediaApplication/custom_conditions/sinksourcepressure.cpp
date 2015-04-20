// Project includes 
#include "includes/define.h"
#include "custom_conditions/sinksourcepressure.h"
#include "porous_media_application.h"
#include "utilities/math_utils.h"


namespace Kratos
{
    
        typedef GeometryData::KratosGeometryType KratosGeometryType;
        
	//************************************************************************************
	//************************************************************************************
	SinkSourcePressure::SinkSourcePressure(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	SinkSourcePressure::SinkSourcePressure(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{

	}
	Condition::Pointer SinkSourcePressure::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new SinkSourcePressure(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	SinkSourcePressure::~SinkSourcePressure()
	{
	}

	//************************************************************************************
	//************************************************************************************
        void SinkSourcePressure::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
            KRATOS_TRY
            switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
            {
                case 1:
               {
                   this->CalculateRightHandSidePressure(rRightHandSideVector, rCurrentProcessInfo);
                   break;
                }
                case 2:
                {
                   break;
                }
                default:
                {
                    KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
                }
            }
            KRATOS_CATCH("")
        }
	void SinkSourcePressure::CalculateRightHandSidePressure(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
            KRATOS_TRY
            
            bool CalculateStiffnessMatrixFlag =false;
            MatrixType temp = Matrix();
            
            this->CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag);
            
            KRATOS_CATCH("")
            
        }
        
        void SinkSourcePressure::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
            KRATOS_TRY
            switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
            {
                case 1:
               {
                   this->CalculateLocalSystemPressure(rLeftHandSideMatrix,rRightHandSideVector, rCurrentProcessInfo);
                    break;
                }
                case 2:
                {
                    break;
                }
                default:
                {
                    KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
                }
            }
            KRATOS_CATCH("")
		
        }
        void SinkSourcePressure::CalculateLocalSystemPressure(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{

            KRATOS_TRY

            bool CalculateStiffnessMatrixFlag =true;
            
            this->CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag);
	
             KRATOS_CATCH("")
            
        }
        void SinkSourcePressure::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
            switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
            {
                case 1:
               {
                   this->EquationIdVectorPressure(rResult, rCurrentProcessInfo);
                    break;
                }
                case 2:
                {
                    break;
                }
                default:
                {
                    KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
                }
            }

	}
        void SinkSourcePressure::EquationIdVectorPressure(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
            	int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int index;
		unsigned int dim = 1;
		rResult.resize(number_of_nodes*dim);
		for (int i=0;i<number_of_nodes;i++)
		{
			index = i*dim;
			rResult[index] = (GetGeometry()[i].GetDof(PRESSURE).EquationId());			
		}
	}

	//************************************************************************************
	//************************************************************************************
	void SinkSourcePressure::GetDofList(DofsVectorType& rConditionalDofList,ProcessInfo& rCurrentProcessInfo)
	{
            switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
            {
                case 1:
               {
                   this->GetDofListPressure(rConditionalDofList, rCurrentProcessInfo);
                    break;
                }
                case 2:
                {
                    break;
                }
                default:
                {
                    KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
                }
            }

	}
        void SinkSourcePressure::GetDofListPressure(DofsVectorType& rConditionalDofList,ProcessInfo& rCurrentProcessInfo)
	{
            unsigned int dim = 1;
            rConditionalDofList.resize(GetGeometry().size()*dim);
            unsigned int index;
            for (unsigned int i=0;i<GetGeometry().size();i++)
            {

                index = i*dim;
                rConditionalDofList[index] = (GetGeometry()[i].pGetDof(PRESSURE));
            }
        }
	//************************************************************************************
	//************************************************************************************
	void SinkSourcePressure::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool aCalculateStiffnessMatrixFlag)
	{
            
		KRATOS_TRY

                //Set Matrix & vector size
                unsigned int number_of_nodes = GetGeometry().size();
                unsigned int MatSize=number_of_nodes;
            
                double sinkSource = GetProperties()[SINK_SOURCE];
                
                if (aCalculateStiffnessMatrixFlag==true)
                {
                    if(rLeftHandSideMatrix.size1() != MatSize)
                        rLeftHandSideMatrix.resize(MatSize,MatSize,false);
                    noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize);
                }
                
                if(rRightHandSideVector.size() != MatSize)
                    rRightHandSideVector.resize(MatSize,false);

                
                KratosGeometryType typeOfElement= GetGeometry().GetGeometryType();
                
                switch ( typeOfElement ) {
                    case  GeometryData::Kratos_Point2D:
                  this->CalculateLocalSystemPoint(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, sinkSource);
                  break;
                    case GeometryData::Kratos_Line2D2:
                  this->CalculateLocalSystemLine(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, sinkSource);
                  break;
                    case GeometryData::Kratos_Triangle2D3:
                  this->CalculateLocalSystemTriangle(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, sinkSource);
                  break;
                default:
                  KRATOS_ERROR(std::logic_error, "This element is not yet implement in order to solve source term!!!!! ERROR", "");
                  break;
                }    
                
                KRATOS_CATCH("");
              
            
	}
        
        
        void SinkSourcePressure::CalculateLocalSystemPoint(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource)
        {
            
            rRightHandSideVector[0] = aSinkSource;
        }
        
        void SinkSourcePressure::CalculateLocalSystemLine(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource)
        {
                //calculate lenght
                double xlenght = GetGeometry()[1].X() - GetGeometry()[0].X();
                double ylenght = GetGeometry()[1].Y() - GetGeometry()[0].Y();
                double lenght = xlenght*xlenght + ylenght*ylenght;
                lenght = sqrt(lenght);
		
		rRightHandSideVector[0] = (aSinkSource*lenght)/2;
                rRightHandSideVector[1] = (aSinkSource*lenght)/2;
		
        }
        
        void SinkSourcePressure::CalculateLocalSystemTriangle(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double aSinkSource)
        {
                //calculate area
                double Area = GetGeometry().Area();
		
		rRightHandSideVector[0] = (aSinkSource*Area)/3;
                rRightHandSideVector[1] = (aSinkSource*Area)/3;
                rRightHandSideVector[2] = (aSinkSource*Area)/3;
                
                
        }


        

         
        
        
	//************************************************************************************
	//************************************************************************************
} // Namespace Kratos
