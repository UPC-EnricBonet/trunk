//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/flowtrans_2d.h"
#include "porous_media_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

/*#include "utilities/enrichment_utilities.h"

#include "processes/process.h"
#include "includes/node.h"
//#include "includes/element.h"
#include "includes/model_part.h"*/

namespace Kratos
{
	


	//************************************************************************************
	//************************************************************************************
	FlowTrans2D::FlowTrans2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	FlowTrans2D::FlowTrans2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
            this->mNumber_of_nodes = GetGeometry().PointsNumber();
            this->SetPropertiesElement();
	}

	Element::Pointer FlowTrans2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new FlowTrans2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	FlowTrans2D::~FlowTrans2D()
	{
	}

        void FlowTrans2D::SetPropertiesElement()
        { 
           this->mDensity = 1;//GetProperties()[DENSITY];
           this->mPorosity = GetProperties()[POROSITY];
           this->mPermeability = GetProperties()[PERMEABILITY_WATER];
           this->mSpecificStorage = GetProperties()[SPECIFIC_STORAGE];
           this->mDiffusionCoefficient = GetProperties()[DIFFUSION_COEFFICIENT];

        }
        
	//************************************************************************************
	//************************************************************************************
        void FlowTrans2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
		case 1:
		{
			this->CalculateFlow(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
			break;
		}
		case 2:
		{
                        this->CalculateDarcyFlow(rCurrentProcessInfo);
			this->CalculateTransport(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
			break;
		}
		default:
		{
			KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
		}

		KRATOS_CATCH("");
	}
        //************************************************************************************
	//************************************************************************************
	void FlowTrans2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}	 

	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void FlowTrans2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}
        /*void FlowTrans2D::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;
                if (rCurrentProcessInfo[FRACTIONAL_STEP]==2)
                    this->CalculateDarcyFlow(rCurrentProcessInfo);
		else
                    KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);


		KRATOS_CATCH("");
	}*/
        void FlowTrans2D::CalculateDarcyFlow(ProcessInfo& rCurrentProcessInfo)
        {
            //Area
            double Area;
            const int Ndim = 2;
            //GradN
            boost::numeric::ublas::bounded_matrix<double, 3, 2 > msDN_DX;
            //N
            array_1d<double, 3 > msN;
            //Geometry staff
	    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);       
            
            //Compute gradient of head
            array_1d<double, Ndim > headGrad;
            

            for (unsigned int node = 0; node < this->mNumber_of_nodes; node++)
            {
               double currentHeadLevel =  GetGeometry()[node].FastGetSolutionStepValue(HEAD_LEVEL);
               for (unsigned int dim = 0; dim < Ndim; dim++)
                {
                   headGrad[dim]+=currentHeadLevel*msDN_DX.at_element(node,dim);
                }
               
            }

            /*
            //if we have a  buoyancy term, substract it from the gradient
            if ( isBuoyancy)
            {
                array_1d<double, Ndim > grav;
                for (unsigned int dim = 0; dim < Ndim; dim++)
                {
                    if(dim= Ndim-1)
                         grav[dim]=-9.8;
                     else
                         grav[dim]=0.0;

                    headGrad[dim] -= grav[dim] ;
                }
            }
            */

            array_1d<double, 2 > darcyFlowVec;
	    double darcyFlowX, darcyFlowY;
            
            //K*gradh
            boost::numeric::ublas::bounded_matrix<double,2,2> mK = ZeroMatrix(2,2);
            mK(0,0)=this->mPermeability;
            mK(1,1)=this->mPermeability;
                 
            darcyFlowVec=-1*prod(mK,trans(headGrad));
            darcyFlowX=darcyFlowVec[0];
            darcyFlowY=darcyFlowVec[1];
            
            rCurrentProcessInfo.SetValue(DARCY_FLOW_X,darcyFlowX);
            rCurrentProcessInfo.SetValue(DARCY_FLOW_Y,darcyFlowY);
            
            
        }
	//************************************************************************************
	//************************************************************************************
        void FlowTrans2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
		case 1:
		{
			this->FlowEquationEquationIdVector(rResult,rCurrentProcessInfo);
			break;
		}
		case 2:
		{
			this->TransportEquationIdVector(rResult,rCurrentProcessInfo);
			break;
		}
		default:
		{
			KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
		}

		KRATOS_CATCH("");
	}
        void FlowTrans2D::FlowEquationEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		
		if(rResult.size() != this->mNumber_of_nodes)
			rResult.resize(this->mNumber_of_nodes,false);	

		for (unsigned int i=0;i<this->mNumber_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(HEAD_LEVEL).EquationId();
	}
        void FlowTrans2D::TransportEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		
		if(rResult.size() != this->mNumber_of_nodes)
			rResult.resize(this->mNumber_of_nodes,false);	

		for (unsigned int i=0;i<this->mNumber_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(CONCENTRATION).EquationId();
	}
	//************************************************************************************
	//************************************************************************************
	void FlowTrans2D::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
		case 1:
		{
			this->GetFlowEquationDofList(rElementalDofList,rCurrentProcessInfo);
			break;
		}
		case 2:
		{
			this->GetTransportEquationDofList(rElementalDofList,rCurrentProcessInfo);
			break;
		}
		default:
		{
			KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
		}

		KRATOS_CATCH("");
	}
        void FlowTrans2D::GetFlowEquationDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		
		if(ElementalDofList.size() != this->mNumber_of_nodes)
			ElementalDofList.resize(this->mNumber_of_nodes);	
                
		for (unsigned int i=0;i<this->mNumber_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(HEAD_LEVEL);
	}
        void FlowTrans2D::GetTransportEquationDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		
		if(ElementalDofList.size() != this->mNumber_of_nodes)
			ElementalDofList.resize(this->mNumber_of_nodes);	
                
		for (unsigned int i=0;i<this->mNumber_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(CONCENTRATION);
	}
        //************************************************************************************
	//************************************************************************************
	void FlowTrans2D::CalculateFlow(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
                 this->mDeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
                 this->mStationary = rCurrentProcessInfo[STATIONARY];
                 
                //dimension of our local matrix (right & left)
                if (rLeftHandSideMatrix.size1() != this->mNumber_of_nodes)
                    rLeftHandSideMatrix.resize(this->mNumber_of_nodes, this->mNumber_of_nodes, false);
                if (rRightHandSideVector.size() != this->mNumber_of_nodes)
                    rRightHandSideVector.resize(this->mNumber_of_nodes, false);
                
                // Storage triangle factor Matrix
                boost::numeric::ublas::bounded_matrix<double, 3, 3 > msMassFactors = 1.0 / 3.0 * IdentityMatrix(3, 3);
                
                //GradN
                boost::numeric::ublas::bounded_matrix<double, 3, 2 > msDN_DX;
                //N?
                array_1d<double, 3 > msN;
                
                //Storage head level from n, n+1 time to interpolate at the storage term
                array_1d<double, 3 > ms_head_np_vec;
                
                //unknow at the previously time (ms_temp in purediffusion))
                array_1d<double,3> ms_head; //dimension = number of nodes //ms_temp=ms_head
                
                //Diffusion matrix constant, two components for each node?? (msD en purediffusion... conductivity)
                boost::numeric::ublas::bounded_matrix<double,2,2> mK = ZeroMatrix(2,2); //initializing the matrix as zero msD = mK 
                
		//getting data for the given geometry (N,gradN & Area)
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);       
                
		//reading properties and conditions
		mK(0,0)= this->mPermeability;
		mK(1,1)= this->mPermeability;
 
                std::ofstream myfile;
                myfile.open ("MatrixSystem.txt", std::ios::app);
                
                ////Contribution from storage to LHS
                double dt_inv = 1.0/ this->mDeltaTime;
                if(this->mStationary!=1)
                    noalias(rLeftHandSideMatrix) = dt_inv * this->mSpecificStorage * msMassFactors;
                
                //Contribution from darcyFlow to LHS
		noalias(rLeftHandSideMatrix) += prod(msDN_DX,Matrix(prod(mK,trans(msDN_DX))));  // Bt D B
                //noalias(rLeftHandSideMatrix) = mPermeability * prod(msDN_DX, trans(msDN_DX);
                
                rLeftHandSideMatrix *= Area;
                
                
                KRATOS_WATCH(rLeftHandSideMatrix);
                myfile << "[Ae]  "<< rLeftHandSideMatrix ;    
                
                ////Contribution from storage to RHS
                if(this->mStationary!=1)
                {
                    array_1d<double, 3 > step_unknown;
                    for (unsigned int iii = 0; iii < this->mNumber_of_nodes; iii++)
                        step_unknown[iii] =  GetGeometry()[iii].FastGetSolutionStepValue(HEAD_LEVEL, 1);
                    noalias(rRightHandSideVector) = dt_inv * this->mSpecificStorage * Area * prod(msMassFactors, step_unknown);
                    
                    myfile << "    [hn]   "<< step_unknown ;
                }
          
		//subtracting the dirichlet term at RHS
		// RHS -= LHS*DUMMY_UNKNOWNs
		for(unsigned int iii = 0; iii<this->mNumber_of_nodes; iii++)
			ms_head[iii] = GetGeometry()[iii].FastGetSolutionStepValue(HEAD_LEVEL);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_head);
                
                KRATOS_WATCH(rRightHandSideVector);
                myfile << "    [be]   "<< rRightHandSideVector;
                
                myfile << "    [unkonws]   "<< ms_head << std::endl << std::endl;
               
                myfile.close();
                
		KRATOS_WATCH(rRightHandSideVector);
                
	}

        void FlowTrans2D::CalculateTransport(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
                 this->mDeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
                 this->mStationary = rCurrentProcessInfo[STATIONARY];
                 
                //dimension of our local matrix (right & left)
                if (rLeftHandSideMatrix.size1() != this->mNumber_of_nodes)
                    rLeftHandSideMatrix.resize(this->mNumber_of_nodes, this->mNumber_of_nodes, false);
                if (rRightHandSideVector.size() != this->mNumber_of_nodes)
                    rRightHandSideVector.resize(this->mNumber_of_nodes, false);
                
                // Storage triangle factor Matrix
                boost::numeric::ublas::bounded_matrix<double, 3, 3 > msMassFactors = 1.0 / 3.0 * IdentityMatrix(3, 3);
                
                //GradN
                boost::numeric::ublas::bounded_matrix<double, 3, 2 > msDN_DX;
                //N?
                array_1d<double, 3 > msN;
                array_1d<double, 3 > ms_u_DN;
                
                array_1d<double, 2 > ms_darcyFlow_vec;
                
                //Storage concentration from n, n+1 time to interpolate at the storage term
                array_1d<double, 3 > ms_concentration_np_vec;
                
                array_1d<double,3> ms_concentration; //dimension = number of nodes //ms_temp=ms_head
                
                //Diffusion coeffient matrix
                boost::numeric::ublas::bounded_matrix<double,2,2> mD = ZeroMatrix(2,2); //initializing the matrix as zero msD = mK
                
                //Velocity vector
                //array_1d<double, 2 > ms_vel_gauss;
                
                //??
                //array_1d<double, 3 > ms_u_DN;
                //array_1d<double, 2 > grad_g;
       

		//getting data for the given geometry (N,gradN & Area)
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);       
                
		mD(0,0)=this->mDiffusionCoefficient;
		mD(1,1)=this->mDiffusionCoefficient;
		
                double darcyFlowX = rCurrentProcessInfo.GetValue(DARCY_FLOW_X);
                double darcyFlowY = rCurrentProcessInfo.GetValue(DARCY_FLOW_Y);
                
                //Ojo!!
                ms_darcyFlow_vec[0]=darcyFlowX;
                ms_darcyFlow_vec[1]=darcyFlowY;
                
                
                std::ofstream myfile;
                myfile.open ("MatrixSystem.txt", std::ios::app);
                
                //Advective contribution to the stiffness matrix
                noalias(ms_u_DN) = prod(msDN_DX, ms_darcyFlow_vec);
                noalias(rLeftHandSideMatrix) = this->mDensity * outer_prod(msN, ms_u_DN);
                
                ////Contribution from storage to LHS
                double dt_inv = 1.0/ this->mDeltaTime;
                if(this->mStationary!=1)
                    noalias(rLeftHandSideMatrix) += dt_inv * this->mDensity * this->mPorosity * msMassFactors;
                
                //Contribution from diffusion to LHS
		noalias(rLeftHandSideMatrix) += this->mDensity* prod(msDN_DX,Matrix(prod(mD,trans(msDN_DX))));  // Bt D B
                //noalias(rLeftHandSideMatrix) = mPermeability * prod(msDN_DX, trans(msDN_DX);
                
                rLeftHandSideMatrix *= Area;
                
                KRATOS_WATCH(rLeftHandSideMatrix);
                myfile << "[Ae]  "<< rLeftHandSideMatrix ;    
                
                ////Contribution from storage to RHS
                if(this->mStationary!=1)
                {
                    array_1d<double, 3 > ms_concentration_vec_np;
                    for (unsigned int iii = 0; iii < this->mNumber_of_nodes; iii++)
                        ms_concentration_vec_np[iii] =  GetGeometry()[iii].FastGetSolutionStepValue(CONCENTRATION, 1);
                    noalias(rRightHandSideVector) = dt_inv * this->mDensity * this->mPorosity * Area * prod(msMassFactors, ms_concentration_vec_np);
                    
                    myfile << "    [hnC]   "<< ms_concentration_vec_np ;
                }
               
                
		//subtracting the dirichlet term at RHS
		// RHS -= LHS*DUMMY_UNKNOWNs
		for(unsigned int iii = 0; iii< this->mNumber_of_nodes; iii++)
			ms_concentration[iii] = GetGeometry()[iii].FastGetSolutionStepValue(CONCENTRATION);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_concentration);
                
                KRATOS_WATCH(rRightHandSideVector);
                myfile << "    [beC]   "<< rRightHandSideVector;
                
                myfile << "    [unkonwsC]   "<< ms_concentration << std::endl << std::endl;
               
                myfile.close();
                
		KRATOS_WATCH(rRightHandSideVector);  
        }
                /*void FlowTrans2D::GetScalarDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);

	}
        void FlowTrans2D::GetVectorDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 3;
		const SizeType LocalSize = 6;
		GeometryType& rGeom = this->GetGeometry();

		if (ElementalDofList.size() != LocalSize)
			ElementalDofList.resize(LocalSize);

		SizeType LocalIndex = 0;

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
		}
	}
         */
         /*void FlowTrans2D::VectorUnkownEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 3;
		const SizeType LocalSize = 6;
		GeometryType& rGeom = this->GetGeometry();

		SizeType LocalIndex = 0;

		if (rResult.size() != LocalSize)
			rResult.resize(LocalSize, false);

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_X).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_Y).EquationId();
		}
	}*/

} // Namespace Kratos
