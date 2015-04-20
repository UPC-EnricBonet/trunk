// Project includes 
#include "includes/define.h"
#include "custom_conditions/leakage.h"
#include "flowtransport_application.h"
#include "utilities/math_utils.h"


namespace Kratos
{
    
        typedef GeometryData::KratosGeometryType KratosGeometryType;
        
	//************************************************************************************
	//************************************************************************************
	Leakage::Leakage(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Leakage::Leakage(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
            //this->mSV= GetProperties()[SV];
            
            this->mNumber_of_nodes = GetGeometry().PointsNumber();
            
            this->mLeakage_coeff = GetProperties()[LEAKAGE_COEFFICIENT];
            
            this->mLevel = GetProperties()[LEVEL];
	}
	Condition::Pointer Leakage::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new Leakage(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	Leakage::~Leakage()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void Leakage::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
                //calculation flags
                bool CalculateStiffnessMatrixFlag = false;
                bool CalculateResidualVectorFlag = true;
                MatrixType temp = Matrix();

                CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
                
		
	}

	//************************************************************************************
	//************************************************************************************
	void Leakage::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		
            //calculation flags
            bool CalculateStiffnessMatrixFlag = true;
            bool CalculateResidualVectorFlag = true;

            CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
                
               
	}
        
        void Leakage::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag)
        {
            KRATOS_TRY
                
                //Set Matrix & vector size
                unsigned int number_of_nodes=this->mNumber_of_nodes;
            
                if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
                {
                    if(rLeftHandSideMatrix.size1() != number_of_nodes)
                            rLeftHandSideMatrix.resize(number_of_nodes,number_of_nodes,false);
                    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_nodes,number_of_nodes);
                }
                
                if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
                {
                    if(rRightHandSideVector.size() != number_of_nodes)
			rRightHandSideVector.resize(number_of_nodes,false);
                }
                
                
                //Catch type of element
                KratosGeometryType typeOfElement= GetGeometry().GetGeometryType();
                 
                switch ( typeOfElement ) {
                    case  GeometryData::Kratos_Point2D:
                  this->CalculateLocalSystemPoint(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
                  break;
                    case GeometryData::Kratos_Line2D2:
                  this->CalculateLocalSystemLine(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
                  break;
                    case GeometryData::Kratos_Triangle2D3:
                  this->CalculateLocalSystemTriangle(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
                  break;
                default:
                  KRATOS_ERROR(std::logic_error, "This element is not yet implement in order to solve source term!!!!! ERROR", "");
                  break;
                }
                
                KRATOS_CATCH("")
        }
        
        
        void Leakage::CalculateLocalSystemPoint(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag)
        {
		
                double currentHeadLevel = GetGeometry()[0].FastGetSolutionStepValue(HEAD_LEVEL);
                ////double stepHeadLevel = GetGeometry()[0].GetSolutionStepValue(HEAD_LEVEL,0);
		
                if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
                    rLeftHandSideMatrix(0,0) = this->mLeakage_coeff*currentHeadLevel;
                
                if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
                    rRightHandSideVector[0] = this->mLeakage_coeff*this->mLevel;
                	
        }
        
        void Leakage::CalculateLocalSystemLine(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag)
        {
            double xlenght = GetGeometry()[1].X() - GetGeometry()[0].X();
            double ylenght = GetGeometry()[1].Y() - GetGeometry()[0].Y();
            double lenght = xlenght*xlenght + ylenght*ylenght;
            lenght = sqrt(lenght);   
            
                
                if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
                {
                    double currentHeadLevel0 = GetGeometry()[0].FastGetSolutionStepValue(HEAD_LEVEL);
                    double currentHeadLevel1 = GetGeometry()[1].FastGetSolutionStepValue(HEAD_LEVEL);
                    rLeftHandSideMatrix(0,0) = lenght*this->mLeakage_coeff*currentHeadLevel0;
                    rLeftHandSideMatrix(1,1) = lenght*this->mLeakage_coeff*currentHeadLevel1;
                }
		
                if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
                {
                    rRightHandSideVector[0] = lenght*this->mLeakage_coeff*this->mLevel;
                    rRightHandSideVector[1] = lenght*this->mLeakage_coeff*this->mLevel;
                }
		
        }
        
        void Leakage::CalculateLocalSystemTriangle(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag)
        {

            double Area = GetGeometry().Area();
            
                if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
                {
                    double currentHeadLevel0 = GetGeometry()[0].FastGetSolutionStepValue(HEAD_LEVEL);
                    double currentHeadLevel1 = GetGeometry()[1].FastGetSolutionStepValue(HEAD_LEVEL);
                    double currentHeadLevel2 = GetGeometry()[2].FastGetSolutionStepValue(HEAD_LEVEL);
                
                    rLeftHandSideMatrix(0,0) = Area*this->mLeakage_coeff*currentHeadLevel0;
                    rLeftHandSideMatrix(1,1) = Area*this->mLeakage_coeff*currentHeadLevel1;
                    rLeftHandSideMatrix(2,2) = Area*this->mLeakage_coeff*currentHeadLevel2;
                }
                
		
                if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
                {
                    rRightHandSideVector[0] = Area*this->mLeakage_coeff*this->mLevel;
                    rRightHandSideVector[1] = Area*this->mLeakage_coeff*this->mLevel;
                    rRightHandSideVector[2] = Area*this->mLeakage_coeff*this->mLevel;
                }
                
        }

	//************************************************************************************
	//************************************************************************************
	void Leakage::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
                rResult.resize(mNumber_of_nodes);
                
                unsigned int index = 0;
                typedef Node<3>::DofsContainerType Dofs;
                typedef Dofs::iterator Dof;

                for (unsigned int node=0;node<mNumber_of_nodes;node++)
		{
                    Node<3>::DofsContainerType& node_all_dofs = GetGeometry()[node].GetDofs();
                    typename Dofs::iterator Dof_begin=node_all_dofs.begin();
                    typename Dofs::iterator Dof_end=node_all_dofs.end();
                    
                    //Loop over all Dof in the node
                    for(Dof currentDof = Dof_begin; currentDof != Dof_end; currentDof++)
                    {
                        std::string currentDofname = currentDof->GetVariable().Name();
                        //if(currentDofname == this->mSV)
                      //  {
                            rResult[index++] = (GetGeometry()[node].GetDof((currentDof)->GetVariable()).EquationId());
                       // }
                    }
                }
            
	}

	//************************************************************************************
	//************************************************************************************
	  void Leakage::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
	{
              	unsigned int index = 0;
                typedef Node<3>::DofsContainerType Dofs;
                typedef Dofs::iterator Dof;
                
                for (unsigned int node=0;node<mNumber_of_nodes;node++)
		{
                    Node<3>::DofsContainerType& node_all_dofs = GetGeometry()[node].GetDofs();
                    typename Dofs::iterator Dof_begin=node_all_dofs.begin();
                    typename Dofs::iterator Dof_end=node_all_dofs.end();
                    
                    //Loop over all Dof in the node
                    for(Dof currentDof = Dof_begin; currentDof != Dof_end; currentDof++)
                    {
                        std::string currentDofname = currentDof->GetVariable().Name();
                        //if(currentDofname == this->mSV)
                       // {
                          //  ConditionalDofList[index++] = (GetGeometry()[node].pGetDof((currentDof)->GetVariable()));
                       // }
                        
                    }
                }
	}
} // Namespace Kratos
