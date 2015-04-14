#ifndef VALARRAYTOOLS_H
#define VALARRAYTOOLS_H

#include <QString>
#include <QStringList>

/// \brief This class contains methods that provide some functionalituy for valarrays, 
/// such ac copying them, allocating them etcetera. 
class ValarrayTools
{
public:
     

    /// \brief this method sets a valarray<valarray<double>> that is already allocated to zero

    static void SetToZero(valarray<double> &in)
    {
        int idim1 = in.size();
        for (int i = 0; i < idim1;i++)
        {
            in[i] = 0;           
        }
    }

/// \brief this method sets a valarray<valarray<double>> that is already allocated to zero

    static void SetToZero(valarray<valarray<double> > &in)
    {
        int idim1 = in.size();
        for (int i = 0; i < idim1;i++)
        {
            int idim2 = in[i].size();
            for (int j = 0; j < idim2; j++)
            {
                in[i][j] = 0;
            }
        }
    }


/// \brief sets valarray to zero 
    static void SetToZero(valarray<valarray<valarray<double> > > &in)
    {
        unsigned int idim1 = in.size();
        for (unsigned int i = 0; i < idim1;i++)
        {
            unsigned int idim2 = in[i].size();
            for (unsigned int j = 0; j < idim2; j++)
            {
                    in[i][j] = 0;
            }

        }
    }

/// \brief sets valarray to zero 
    static void SetToZero(valarray<valarray<valarray<valarray<double> > > > &in)
    {
        unsigned int idim1 = in.size();
        for (unsigned int i = 0; i < idim1;i++)
        {
            unsigned int idim2 = in[i].size();
            for (unsigned int j = 0; j < idim2; j++)
            {

                unsigned int idim3 = in[i][j].size();
                for (unsigned int k = 0; k < idim3; k++)
                {

                    in[i][j][k] = 0;
                }

            }

        }
    }

/// \brief this method sets a valarray<valarray<double>> that is already allocated to a value

    static void SetToValue(valarray<valarray<double> > &in, double aValue)
    {
        int idim1 = in.size();
        for (int i = 0; i < idim1;i++)
        {
            int idim2 = in[i].size();
            for (int j = 0; j < idim2; j++)
            {
                in[i][j] =aValue;
            }
        }
    }

// this method resizes a valarray<valarray<double>> to match the shape of another one
//and fills it with the value "val"
    static void
    Allocate(valarray<valarray<double> > &out,const valarray<valarray<double> > &in,  const double val)
    {
        int idim1 = in.size();
        out.resize(idim1);

        for (int i = 0; i< idim1;i++)
        {
            int idim2 = in[i].size();
            out[i].resize(idim2, val);      
        }
    }

	static void
    Allocate(valarray<valarray<valarray<valarray<double> > > > &out, const valarray<valarray<valarray<valarray<double> > > > &in,  const double val)
    {
        int idim1 = in.size();
        out.resize(idim1);

        for (int i = 0; i< idim1;i++)
        {
            int idim2 = in[i].size();
			out[i].resize(idim2);
			for(int j = 0; j< idim2;j++)
			{
				int idim3 = in[i][j].size();
				out[i][j].resize(idim3);
				for(int k = 0; k<idim3;k++)
				{
					int idim4 = in[i][j][k].size();
					out[i][j][k].resize(idim4, val);      
				}
			}
        }
    }

/// \brief copies a valarray into another one. Does not allocate. 
    static void
    Copy(valarray<valarray<double> > &out ,const  valarray<valarray<double> > &in)
    {
        int idim1 = in.size();

        for (int i = 0; i< idim1;i++)
        {
            out[i]= in[i];
        }
    }

/// \brief computes out = weight1* in1 + weight2 * in2
    static void
    LinearCombination(valarray<valarray<double> > &out
                ,const double w1, 
                const  valarray<valarray<double> > &in1
                ,const double w2, 
                const  valarray<valarray<double> > &in2)
    {

        int idim1 = in1.size();
        for (int i = 0; i < idim1;i++)
        {
            int idim2 = in1[i].size();
            for (int j = 0; j < idim2; j++)
            {
                out[i][j] = w1 *  in1[i][j] + w2 *  in2[i][j] ;
            }
        }
    }

/// \brief computes out = out + weight1* in1
    static void
    Add(valarray<valarray<double> > &out
                ,const double w1, 
                const  valarray<valarray<double> > &in1
               )
    {

        int idim1 = in1.size();
        for (int i = 0; i < idim1;i++)
        {
            int idim2 = in1[i].size();
            for (int j = 0; j < idim2; j++)
            {
                out[i][j] += w1 *  in1[i][j];
            }
        }
    }    

/// \brief computes out = weight1* in1 + weight2 * in2
    static void
    LinearCombination(valarray<valarray<valarray<valarray<double> > > > &out
                ,const double w1, 
                const  valarray<valarray<valarray<valarray<double> > > > &in1
                ,const double w2, 
                const  valarray<valarray<valarray<valarray<double> > > > &in2)
    {

        int idim1 = in1.size();
        for (int i = 0; i < idim1;i++)
        {
            int idim2 = in1[i].size();
            for (int j = 0; j < idim2; j++)
            {
				int idim3 = in1[i][j].size();
				for (int k= 0; k < idim3;k++)
				{
					int idim4 = in1[i][j][k].size();
					for (int l = 0; l < idim4;l++)
					{
						out[i][j][k][l] = w1 *  in1[i][j][k][l] + w2 *  in2[i][j][k][l] ;
					}
				}
            }
        }
    }
/// \brief computes out = weight1* in1 
    static void
    LinearCombination(valarray<valarray<valarray<valarray<double> > > > &out
                ,const double w1, 
                const  valarray<valarray<valarray<valarray<double> > > > &in1)
    {

        int idim1 = in1.size();
        for (int i = 0; i < idim1;i++)
        {
            int idim2 = in1[i].size();
            for (int j = 0; j < idim2; j++)
            {
				int idim3 = in1[i][j].size();
				for (int k= 0; k < idim3;k++)
				{
					int idim4 = in1[i][j][k].size();
					for (int l = 0; l < idim4;l++)
					{
						out[i][j][k][l] = w1 *  in1[i][j][k][l] ;
					}
				}
            }
        }
    }

// \brief computes out = weight1* in1
    static void
    LinearCombination(valarray<valarray<double> > &out
                ,const double w1, 
                const  valarray<valarray<double> > &in1)
    {

        int idim1 = in1.size();
        for (int i = 0; i < idim1;i++)
        {
            int idim2 = in1[i].size();
            for (int j = 0; j < idim2; j++)
            {
                out[i][j] = w1 *  in1[i][j] ;
            }
        }
    }

// \brief checks if 2 are equal
    static bool
    IsEqual(const valarray<valarray<double> > &in1
           ,const valarray<valarray<double> > &in2
           ,const double absmaxdiff)
    {

        int i1dim1 = in1.size();
        int i2dim1 = in2.size(); 
        if ( i1dim1 != i2dim1) return false;


        for (int i = 0; i < i1dim1;i++)
        {
            int idim2 = in1[i].size();
            if ( in2[i].size() != idim2 ) return false;
            for (int j = 0; j < idim2; j++)
            {
                if (in1[i][j] != in2[i][j])
                {
                    if (abs(in1[i][j] - in2[i][j]) > absmaxdiff) return false;
                }
            }
        }

        return true;
    }

// \brief checks if 2 are equal
    static bool
    IsEqual(const valarray<valarray<valarray<valarray<double> > > > &in1
           ,const valarray<valarray<valarray<valarray<double> > > > &in2
           ,const double absmaxdiff)
    {

        int i1dim1 = in1.size();
        int i2dim1 = in2.size(); 
        if ( i1dim1 != i2dim1) return false;

        for (int i = 0; i < i1dim1;i++)
        {
            int idim2 = in1[i].size();
            if ( in2[i].size() != idim2 ) return false;
            for (int j = 0; j < idim2; j++)
            {
				int idim3 = in1[i][j].size();
				if ( in2[i][j].size() != idim3 ) return false;
				for (int k = 0; k < idim3; k++)
				{
					int idim4 = in1[i][j][k].size();
					if ( in2[i][j][k].size() != idim4 ) return false;
					for (int l = 0; l < idim4; l++)
					{	
						if (in1[i][j][k][l] != in2[i][j][k][l])
						{
							if (abs(in1[i][j][k][l] - in2[i][j][k][l]) > absmaxdiff) return false;
						}
					}
				}
            }
        }

        return true;
    }



// \brief reads one from file
    static void 
    fromString(const QString str,  valarray<valarray<double> > &out)
    {
        QStringList strList = str.split(";");
        out.resize(strList.size());
        for (unsigned i = 0; i < out.size(); i++)
        {
            out[i].resize(1);
            out[i][0] = strList[i].toDouble();
        }
    }


	// \brief reads one from file
    static void 
    fromString(const QString str,  valarray<valarray<valarray<valarray<double> > > > &out)
    {
        QStringList strList1 = str.split(";");
        out.resize(strList1.size());

        for (unsigned i = 0; i < out.size(); i++)
        {
			QStringList strList2 = strList1[i].split("*");
			out[i].resize(strList2.size());
			
			for (unsigned j=0; j < out[i].size();j++)
			{
				QStringList strList3 = strList2[j].split(",");
				out[i][j].resize(strList3.size());
			
				for (unsigned k=0; k < out[i][j].size();k++)
				{

					QStringList strList4 = strList3[k].split("|");
					out[i][j][k].resize(strList4.size());				

					for (unsigned m=0;m <  out[i][j][k].size();m++)
					{
						out[i][j][k][m]= strList4[m].toDouble();
					}
				}
			}
        }
    }

};
#endif