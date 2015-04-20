// Project includes 
#include "includes/define.h"
#include "custom_conditions/sinksource.h"
#include "flowtransport_application.h"
#include "utilities/math_utils.h"


namespace Kratos
{
    
        typedef GeometryData::KratosGeometryType KratosGeometryType;
        
	//************************************************************************************
	//************************************************************************************
	SinkSource::SinkSource(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	SinkSource::SinkSource(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
            
           /* std::vector <std::string> SVTypes(4);
            SVTypes[0]= GetProperties()[SV];
            SVTypes[1]= GetProperties()[SV_X];
            SVTypes[2]= GetProperties()[SV_Y];
            SVTypes[3]= GetProperties()[SV_Z];
            
            mSinkSourceIsVectorial = false;
            if(SVTypes[0].empty())
                mSinkSourceIsVectorial = true;

            for ( std::vector<std::string>::iterator it = SVTypes.begin(); it != SVTypes.end(); ++it )
            {
                if(!(*it).empty())
                    mSVs.push_back(*it);
            }
 
            mNumber_of_nodes = GetGeometry().PointsNumber();
            
            */mSinkSource = GetProperties()[SINK_SOURCE];
            
	}
	Condition::Pointer SinkSource::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new SinkSource(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	SinkSource::~SinkSource()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void SinkSource::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
            bool CalculateStiffnessMatrixFlag =false;
            MatrixType temp = Matrix();
            
            this->CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag);
        }

        void SinkSource::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
            bool CalculateStiffnessMatrixFlag =true;
            
            this->CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag);
        }
        
        void SinkSource::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const unsigned int conditionDofs = mSVs.size(); 
                rResult.resize(mNumber_of_nodes*conditionDofs);
                 
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
                        // Loop over Dof at the condition
                        for (std::vector<std::string>::iterator SV_i = mSVs.begin() ; SV_i != mSVs.end(); ++SV_i)
                        {
                            if(currentDofname == *SV_i)
                            {
                                rResult[index++] = (GetGeometry()[node].GetDof((currentDof)->GetVariable()).EquationId());
                            }
                        }
                    }
                }
            
               /*
                const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		rResult.resize(number_of_nodes*this->mNumberOfUnknows);
		unsigned int index;
                
                for (unsigned int i=0;i<number_of_nodes;i++)
		{
                        unsigned int nodesNumberOfDof =GetGeometry()[i].GetDofs().size();
			index = i*this->mNumberOfUnknows;
			rResult[index] = (GetGeometry()[i].GetDof(this->mSV).EquationId());	//DofsContainerType		
		}
                */
                
            //Get idEquation of all Dofs
            /*
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
                typedef Node<3>::DofsContainerType Dofs;
                typedef Dofs::iterator Dof;
                
                for (unsigned int node=0;node<number_of_nodes;node++)
		{
                    Node<3>::DofsContainerType& node_dofs = GetGeometry()[node].GetDofs();
                    typename Dofs::iterator Dof_begin=node_dofs.begin();
                    typename Dofs::iterator Dof_end=node_dofs.end();
                    
                    for(Dof it = Dof_begin; it != Dof_end; it++)
                        rResult.push_back(GetGeometry()[node].GetDof((it)->GetVariable()).EquationId());
                }
                
             */
	}

	//************************************************************************************
	//************************************************************************************
	void SinkSource::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
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
                        // Loop over Dof at the condition
                        for (std::vector<std::string>::iterator SV_i = mSVs.begin() ; SV_i != mSVs.end(); ++SV_i)
                        {
                            if(currentDofname == *SV_i)
                            {
                                ConditionalDofList[index++] = (GetGeometry()[node].pGetDof((currentDof)->GetVariable()));
                            }
                        }
                    }
                }
                
                
                //No Abstration code
                /*
		ConditionalDofList.resize(GetGeometry().PointsNumber()*this->mNumberOfUnknows);
		unsigned int index;
		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			index = i*this->mNumberOfUnknows;
			ConditionalDofList[index] = (GetGeometry()[i].pGetDof(this->mSV));
		}
                */
         
	}
        
	//************************************************************************************
	//************************************************************************************
	void SinkSource::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool aCalculateStiffnessMatrixFlag)
	{
		KRATOS_TRY

                //Set Matrix & vector size
                const unsigned int conditionDofs = mSVs.size();
                const unsigned int MatSize=mNumber_of_nodes*conditionDofs;
            
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
                    this->CalculateLocalSystemPoint(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
                    break;
                  case GeometryData::Kratos_Line2D2:
                    this->CalculateLocalSystemLine(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
                    break;
                  case GeometryData::Kratos_Triangle2D3:
                    this->CalculateLocalSystemTriangle(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
                    break;
                  default:
                    KRATOS_ERROR(std::logic_error, "This element is not yet implement in order to solve source term!!!!! ERROR", "");
                    break;
                }    
                
                KRATOS_CATCH("");
	}
        
        
        void SinkSource::CalculateLocalSystemPoint(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
        {
            
		unsigned int index = 0;
                for (unsigned int node=0;node<mNumber_of_nodes;node++)
		{
                    for (std::vector<std::string>::iterator SV_i = mSVs.begin() ; SV_i != mSVs.end(); ++SV_i)
                    {
                        rRightHandSideVector[index++] = this->mSinkSource/double(mNumber_of_nodes);
                    }
                }
        }
        
        void SinkSource::CalculateLocalSystemLine(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
        {

            double xlenght = GetGeometry()[1].X() - GetGeometry()[0].X();
            double ylenght = GetGeometry()[1].Y() - GetGeometry()[0].Y();
            double lenght = xlenght*xlenght + ylenght*ylenght;
            lenght = sqrt(lenght);
                
            if(this->mSinkSourceIsVectorial)
            {
                array_1d<double,3> Normal;
                this->CalculateNormalFace2D(Normal);
                Normal /= (-1.0*lenght);//Inputs?
                
                
                unsigned int index = 0;
                for (unsigned int node=0;node<mNumber_of_nodes;node++)
		{
                    for (std::vector<std::string>::iterator SV_i = mSVs.begin() ; SV_i != mSVs.end(); ++SV_i)
                    {
                        int dimIndex =SV_i - mSVs.begin();
                        rRightHandSideVector[index++] = (this->mSinkSource*lenght*Normal[dimIndex])/double(mNumber_of_nodes);
                    }
                } 
            }
            else
            {
		
                unsigned int index = 0;
                for (unsigned int node=0;node<mNumber_of_nodes;node++)
		{
                    for (std::vector<std::string>::iterator SV_i = mSVs.begin() ; SV_i != mSVs.end(); ++SV_i)
                    {
                        rRightHandSideVector[index++] = (this->mSinkSource*lenght)/double(mNumber_of_nodes);
                    }
                }
               
            }
		
        }
        
        void SinkSource::CalculateLocalSystemTriangle(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
        {
            double Area = GetGeometry().Area();
            
            if(this->mSinkSourceIsVectorial)
            {
		 array_1d<double,3> Normal;
                this->CalculateNormalFace3D(Normal);
                double lenght = std::sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);
                Normal /= lenght;

                unsigned int index = 0;
                for (unsigned int node=0;node<mNumber_of_nodes;node++)
		{
                    for (std::vector<std::string>::iterator SV_i = mSVs.begin() ; SV_i != mSVs.end(); ++SV_i)
                    {
                        int dimIndex =SV_i - mSVs.begin();
                        rRightHandSideVector[index++] = (this->mSinkSource*Area*Normal[dimIndex])/double(mNumber_of_nodes);
                    }
                }
            }
            else
            {
                
                unsigned int index = 0;
                for (unsigned int node=0;node<mNumber_of_nodes;node++)
		{
                    for (std::vector<std::string>::iterator SV_i = mSVs.begin() ; SV_i != mSVs.end(); ++SV_i)
                    {
                        rRightHandSideVector[index++] = (this->mSinkSource*Area)/double(mNumber_of_nodes);
                    }
                }
                
            }
                
        }

        void SinkSource::CalculateNormalFace2D(array_1d<double,3>& aN)
        {
            Geometry<Node<3> >& pGeometry = this->GetGeometry();

            aN[0] =   pGeometry[1].Y() - pGeometry[0].Y();
            aN[1] = - (pGeometry[1].X() - pGeometry[0].X());
            aN[2] =    0.00;
                
        }
        
        
         
        void SinkSource::CalculateNormalFace3D(array_1d<double,3>& aN)
        {
            Geometry<Node<3> >& pGeometry = this->GetGeometry();

            array_1d<double,3> v1,v2;
            v1[0] = pGeometry[1].X() - pGeometry[0].X();
            v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
            v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

            v2[0] = pGeometry[2].X() - pGeometry[0].X();
            v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
            v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

            MathUtils<double>::CrossProduct(aN,v1,v2);
            aN *= 0.5;
        }
         
        
        
	//************************************************************************************
	//************************************************************************************
} // Namespace Kratos
