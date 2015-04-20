//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/flow_2d.h"
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
	Flow2D::Flow2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Flow2D::Flow2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer Flow2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new Flow2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	Flow2D::~Flow2D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void Flow2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

                const double delta_t = rCurrentProcessInfo.GetValue(DELTA_TIME);
                const int Stationary = rCurrentProcessInfo[STATIONARY];
                //const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

                //get number of points of element... Generic, but its always a triangle
                //if(GetGeometry().Dimension()!=2)
                    //KRATOS_ERROR(std::logic_error,  "Element dimension is wrong" , "");
                //unsigned int Ndim=GetGeometry().Dimension();
                const unsigned int number_of_points = GetGeometry().size(); 
                
                //dimension of our local matrix (right & left)
                if (rLeftHandSideMatrix.size1() != number_of_points)
                    rLeftHandSideMatrix.resize(number_of_points, number_of_points, false);
                if (rRightHandSideVector.size() != number_of_points)
                    rRightHandSideVector.resize(number_of_points, false);
                
                const unsigned int number_of_points1=3;
                const unsigned int NDim=2;
                
                // Storage triangle factor Matrix
                boost::numeric::ublas::bounded_matrix<double, number_of_points1, number_of_points1 > msMassFactors = 1.0 / 3.0 * IdentityMatrix(number_of_points1, number_of_points1);
                
                //GradN
                boost::numeric::ublas::bounded_matrix<double, number_of_points1, NDim > msDN_DX;
                //N?
                array_1d<double, number_of_points1 > msN;
                
                //Storage head level from n, n+1 time to interpolate at the storage term
                array_1d<double, number_of_points1 > ms_head_np_vec;
                
                //unknow at the previously time (ms_temp in purediffusion))
                //array_1d<double, number_of_points1 > ms_temp_vec_np;
                array_1d<double,number_of_points1> ms_head; //dimension = number of nodes //ms_temp=ms_head
                
                //Diffusion matrix constant, two components for each node?? (msD en purediffusion... conductivity)
                //boost::numeric::ublas::bounded_matrix<double, 2, 2 > Identity = IdentityMatrix(2, 2);
                //boost::numeric::ublas::bounded_matrix<double, 2, 2 > First;
                //boost::numeric::ublas::bounded_matrix<double, 2, 2 > Second;
                //boost::numeric::ublas::bounded_matrix<double, 2, 3 > Third;
                boost::numeric::ublas::bounded_matrix<double,NDim,NDim> mK = ZeroMatrix(NDim,NDim); //initializing the matrix as zero msD = mK
                
                //Velocity vector
                //array_1d<double, 2 > ms_vel_gauss;
                
                //??
                //array_1d<double, 3 > ms_u_DN;
                //array_1d<double, 2 > grad_g;
       

		//getting data for the given geometry (N,gradN & Area)
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);       
                
		//reading properties and conditions
		double permeability = GetProperties()[PERMEABILITY_WATER];
		mK(0,0)=permeability;
		mK(1,1)=permeability;
                
		//reading properties and conditions
		double specificStorage = GetProperties()[SPECIFIC_STORAGE];
		
		// main loop	
		//const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

                
                /*if(Stationary!=1){
                     noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * specificStorage * msMassFactors;

                     for (unsigned int iii = 0; iii < number_of_points; iii++)
                            ms_head_np_vec[iii] = BDFcoeffs[1] * GetGeometry()[iii].FastGetSolutionStepValue(HEAD_LEVEL, 1);
                                for (unsigned int step = 2; step < BDFcoeffs.size(); step++)
                        {
                        for (unsigned int iii = 0; iii < number_of_points; iii++)
                            ms_head_np_vec[iii] += BDFcoeffs[step] * GetGeometry()[iii].FastGetSolutionStepValue(HEAD_LEVEL, step);
                         }
                        noalias(rRightHandSideVector) -= prod(msMassFactors, ms_head_np_vec * specificStorage);
                    }*/
                
                std::ofstream myfile;
                myfile.open ("MatrixSystem.txt", std::ios::app);
                
                ////Contribution from storage to LHS
                double dt_inv = 1.0/ delta_t;
                if(Stationary!=1)
                    noalias(rLeftHandSideMatrix) = dt_inv * specificStorage * msMassFactors;
                
                //Contribution from darcyFlow to LHS
		noalias(rLeftHandSideMatrix) += prod(msDN_DX,Matrix(prod(mK,trans(msDN_DX))));  // Bt D B
                //noalias(rLeftHandSideMatrix) = permeability * prod(msDN_DX, trans(msDN_DX);
                
                rLeftHandSideMatrix *= Area;
                
                KRATOS_WATCH(rLeftHandSideMatrix);
                myfile << "[Ae]  "<< rLeftHandSideMatrix ;    
                
                ////Contribution from storage to RHS
                if(Stationary!=1)
                {
                    array_1d<double, 3 > step_unknown;
                    for (unsigned int iii = 0; iii < number_of_points; iii++)
                        step_unknown[iii] =  GetGeometry()[iii].FastGetSolutionStepValue(HEAD_LEVEL, 1);
                    noalias(rRightHandSideVector) = dt_inv * specificStorage * Area * prod(msMassFactors, step_unknown);
                    
                    myfile << "    [hn]   "<< step_unknown ;
                }
               
                
		//subtracting the dirichlet term at RHS
		// RHS -= LHS*DUMMY_UNKNOWNs
		for(unsigned int iii = 0; iii<number_of_points; iii++)
			ms_head[iii] = GetGeometry()[iii].FastGetSolutionStepValue(HEAD_LEVEL);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_head);
                
                KRATOS_WATCH(rRightHandSideVector);
                myfile << "    [be]   "<< rRightHandSideVector;
                
                myfile << "    [unkonws]   "<< ms_head << std::endl << std::endl;
               
                myfile.close();
                
		KRATOS_WATCH(rRightHandSideVector);
                
                CalculateDarcyFlow(rCurrentProcessInfo);
                
		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
        void Flow2D::CalculateDarcyFlow(ProcessInfo& rCurrentProcessInfo)
        {
            //Area
            double Area;
           
            //if(GetGeometry().Dimension()!=2)
                //KRATOS_ERROR(std::logic_error,  "Element dimension is wrong" , "");
            //unsigned int Ndim=GetGeometry().Dimension();
            //const unsigned int number_of_nodes = aNPoints; 
            const unsigned int number_of_nodes = 3;//GetGeometry().PointsNumber();
            const unsigned int Ndim=2;//aDim;
             
            //GradN
            boost::numeric::ublas::bounded_matrix<double, number_of_nodes, Ndim > msDN_DX;
            //N
            array_1d<double, number_of_nodes > msN;
            //Geometry staff
	    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);       
            
            //Compute gradient of head
            array_1d<double, Ndim > headGrad;

            for (unsigned int node = 0; node < number_of_nodes; node++)
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
            boost::numeric::ublas::bounded_matrix<double,Ndim,Ndim> mK = ZeroMatrix(Ndim,Ndim);
            double permeability = GetProperties()[PERMEABILITY_WATER];
            mK(0,0)=permeability;
            mK(1,1)=permeability;
                 
            darcyFlowVec=prod(mK,trans(headGrad));
            darcyFlowX=darcyFlowVec[0];
            darcyFlowY=darcyFlowVec[1];
            //array_1d<double, 2 > darcyFlowVec;
	    //const double darcyFlowX, darcyFlowX;
            
            //noalias(rLeftHandSideMatrix) += prod(msDN_DX,Matrix(prod(mK,trans(msDN_DX)))); 
            
            rCurrentProcessInfo.SetValue(DARCY_FLOW_X,darcyFlowX);
            rCurrentProcessInfo.SetValue(DARCY_FLOW_Y,darcyFlowY);
            double dX = rCurrentProcessInfo.GetValue(DARCY_FLOW_X);
            double dY = rCurrentProcessInfo.GetValue(DARCY_FLOW_Y);
            
            
            
        }
	void Flow2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}	 

	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void Flow2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void Flow2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(HEAD_LEVEL).EquationId();
	}

	//************************************************************************************
	//************************************************************************************
	  void Flow2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	
                
		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(HEAD_LEVEL);
                

	}



} // Namespace Kratos
