#ifndef CONTAINERTOOLS_H
#define CONTAINERTOOLS_H

#include <vector>

  template< class elemtype> class ContainerTools
{
  public:
    /// \brief sums a 3d container to a similarly shaped 3d container
  static void Sum
    (
       vector< vector<vector< elemtype > > >& target, ///< container we sum into
        vector< vector<vector< elemtype > > >& src,    ///< container we sum
        elemtype factor);                                          ///< factor multiplying src if any


   /// \brief computes A = w1*B + w2*C where A,B and C are equally shaped 4d containers, and w1 and w2 are scalars
  static void LinearCombination
   (
       vector< vector<vector<vector<elemtype > > > >& result, ///< we equal this to w1*B + w2*C
       vector< vector<vector<vector<elemtype > > > >& B,
       elemtype w1,
       vector< vector<vector<vector<elemtype > > > >& C,
       elemtype w2);


  /// \brief computes A = w1*B  where A and B  are equally shaped 4d containers, and w1 a scalar
  static void LinearCombination
  (
      vector< vector<vector<vector<elemtype > > > >& result, ///< we equal this to w1*B + w2*C
      vector< vector<vector<vector<elemtype > > > >& B,
      elemtype w1);

    /// \brief Sets a 2d array to zero
  static void SetToZero
    (vector< vector<elemtype> >& A);  
	
	    /// \brief Sets a 2d array to zero
  static void SetToZero
    (vector< vector< vector< vector< elemtype> > > >& A);  


  static elemtype Distance(const vector<elemtype>& p1,const vector<elemtype>& p2);

  static vector<elemtype> UnitVec(const vector<elemtype>& p1);


  static double norm(const vector<elemtype>& p1);

 /// \Set all the nodes of a subset
 //  This method obtains a matrix of nodes called curCoord from two vectors of nodes, where the matrix curCoordJ and curCoordS contains several vectors, 
 //  the method obtains for each node of a vector the possible combinations with all nodes of the other vector. 

  static void setMatrixFromVectors (vector<vector<elemtype>> &curCoord, int i, int nodeIndex, vector<vector<elemtype>> &curCoordS, vector<vector<elemtype>> &curCoordJ);
};


template< class elemtype> void
    ContainerTools<elemtype>::Sum(
      vector< vector<vector< elemtype > > >& target,
      vector< vector<vector< elemtype > > >& src,
      elemtype factor)
{
    unsigned idim1 = target.size();
    for (unsigned itgt = 0; itgt < idim1; itgt++)
    {
        vector<vector< elemtype> >* tgtdim1 = &target[itgt];
        vector<vector< elemtype> >* srcdim1 = &src[itgt];

        unsigned idim2 =  tgtdim1->size();
        for (unsigned j = 0; j < idim2; j++)
        {
             vector< elemtype>* tgtdim2 = &((*tgtdim1)[j]);
             vector< elemtype>* srcdim2 = &((*srcdim1)[j]);

             unsigned idim3 = tgtdim2->size();
             for (unsigned k = 0; k < idim3; k++)
             {
                 (*tgtdim2)[k]+= srcdim2->at(k)* factor;
             }

        }
    }
}
template< class elemtype>
void ContainerTools<elemtype>::LinearCombination
 (
     vector< vector<vector<vector<elemtype > > > >& result, ///< we equal this to w1*B + w2*C
     vector< vector<vector<vector<elemtype > > > >& B,
     elemtype w1,
     vector< vector<vector<vector<elemtype > > > >& C,
     elemtype w2)
{
    unsigned idim1 = result.size();
    for (unsigned itgt = 0; itgt < idim1; itgt++)
    {
        vector<vector<vector< elemtype> > >* resultdim1 = &result[itgt];
        vector<vector<vector< elemtype> > >* Bdim1 = &B[itgt];
        vector<vector<vector< elemtype> > >* Cdim1 = &C[itgt];

        unsigned idim2 =  resultdim1->size();
        for (unsigned j = 0; j < idim2; j++)
        {
            vector<vector< elemtype> >* resultdim2 = &(*resultdim1)[j];
            vector<vector< elemtype> >* Bdim2 = &(*Bdim1)[j];
            vector<vector< elemtype> >* Cdim2 = &(*Cdim1)[j];

             unsigned idim3 = resultdim2->size();
             for (unsigned k = 0; k < idim3; k++)
             {

                 vector< elemtype>* resultdim3 = &(*resultdim2)[k];
                 vector< elemtype>* Bdim3 = &(*Bdim2)[k];
                 vector< elemtype>* Cdim3 = &(*Cdim2)[k];

                 unsigned idim4 = resultdim3->size();
                 for (unsigned l = 0; l < idim4; l++)
                 {
                     (*resultdim3)[l] = w1* Bdim3->at(l) + w2* Cdim3->at(l);
                 }

             }

        }
    }

}

template<class elemtype>
void  ContainerTools<elemtype>::LinearCombination
 (
     vector< vector<vector<vector<elemtype > > > >& result, ///< we equal this to w1*B + w2*C
     vector< vector<vector<vector<elemtype > > > >& B,
     elemtype w1)
{
    unsigned idim1 = result.size();
    for (unsigned itgt = 0; itgt < idim1; itgt++)
    {
        vector<vector<vector< elemtype> > >* resultdim1 = &result[itgt];
        vector<vector<vector< elemtype> > >* Bdim1 = &B[itgt];

        unsigned idim2 =  resultdim1->size();
        for (unsigned j = 0; j < idim2; j++)
        {
            vector<vector< elemtype> >* resultdim2 = &(*resultdim1)[j];
            vector<vector< elemtype> >* Bdim2 = &(*Bdim1)[j];

             unsigned idim3 = resultdim2->size();
             for (unsigned k = 0; k < idim3; k++)
             {

                 vector< elemtype>* resultdim3 = &(*resultdim2)[k];
                 vector< elemtype>* Bdim3 = &(*Bdim2)[k];

                 unsigned idim4 = resultdim3->size();
                 for (unsigned l = 0; l < idim4; l++)
                 {
                     (*resultdim3)[l] = w1* Bdim3->at(l);
                 }

             }

        }
    }

}
    /// \brief Sets a 2d array to zero
template<class elemtype>
void  ContainerTools<elemtype>::SetToZero
    (
       vector< vector< elemtype > > & A)
{
    unsigned idim1 = A.size();
    for (unsigned itgt = 0; itgt < idim1; itgt++)
    {
        vector<elemtype>* B = &A[itgt];
        unsigned idim2 =  B->size();
        for (unsigned j = 0; j < idim2; j++)
        {
            (*B)[j] = 0;
        }
    }
}
    /// \brief Sets a 2d array to zero
template<class elemtype>
void  ContainerTools<elemtype>::SetToZero
    (  vector< vector< vector< vector<  elemtype > > > >& A)
{
   typename vector< vector< vector< vector<  elemtype > > > >::iterator A_1;
   typename vector< vector< vector<  elemtype > > >::iterator A_2;
   typename vector< vector<  elemtype > >::iterator A_3;
   
   for(A_1 = A.begin(); A_1 != A.end(); A_1++)
   {
	   for(A_2 = A_1->begin(); A_2 != A_1->end(); A_2++)
	   {
			for (A_3 = A_2->begin(); A_3 != A_2->end(); A_3++)
			{
				fill(A_3->begin(), A_3->end(), 0);
			}
	   }	
   }
}

	///! \brief Inline function to compute distance between two points
template<class elemtype>
elemtype ContainerTools<elemtype>::Distance(const vector<elemtype>& p1, const vector<elemtype>& p2){
		elemtype out= 0;
		for (unsigned i = 0;i < p1.size();i++){
			out += pow(p1[i]-p2[i],2);
		}
		out = sqrt(out);
		return out;
	};


template<class elemtype>
vector<elemtype> ContainerTools<elemtype>::UnitVec(const vector<elemtype>& p1)
{
    // make unit vector in direction
    const int dim = p1.size();
    vector<double> udir(dim ) ;
    
    double norm= ContainerTools<elemtype>::norm(p1);

    if (norm == 0)
    {
        udir.resize(0);
        return udir;
    }
    
    for (int i = 0; i < dim; i++)
    {
        udir[i] =  p1[i]/norm;
    }

    return udir; 
}
template<class elemtype>
double ContainerTools<elemtype>::norm(const vector<elemtype>& p1)
{
    const int dim = p1.size();

    
    double dnorm= 0;
    for (int i = 0; i < dim; i++)
    {
         dnorm += p1[i]*p1[i];
    }
     dnorm = sqrt( dnorm);   
   
    return  dnorm; 
}

template<class elemtype>
void ContainerTools<elemtype>::setMatrixFromVectors (vector<vector<elemtype>> &curCoord, int i, int nodeIndex, vector<vector<elemtype>> &curCoordS, vector<vector<elemtype>> &curCoordJ)
{
	vector<elemtype>curCoordE1 = (curCoordJ)[i];
	vector<elemtype>curCoordE2 = (curCoordS)[i+1];
	int nIndex = nodeIndex;
	nodeIndex=0;
	if(nodeIndex==0)nIndex=curCoordE1.size();
	(curCoord).resize(nIndex*curCoordE2.size());
	for (unsigned n =0 ; n < curCoordE2.size() ;n++)
	{
		for(int m = 0; m<nIndex; m++)
	   {
		   elemtype v2 = curCoordE2[n];
		   if(i==0)
		   {
			   vector<elemtype> v;
			   elemtype v1 = curCoordE1[m];
			   v.push_back(v1);
			   v.push_back(v2);
			   (curCoord)[nodeIndex]=v;
		   }
		   else
		   {
		   vector<elemtype> v=(curCoordJ)[m];
		   v.push_back(v2);
		   (curCoord)[nodeIndex]=v;
		   }
		   nodeIndex++;
	   }
	}
    (curCoordJ)=(curCoord);
}
#endif

