/* 
 * File:   LibChemist.h
 * Author: enric
 *
 * Created on 20 de noviembre de 2014, 12:37
 */

#ifndef LIBCHEMIST_H
#define	LIBCHEMIST_H
#include <math.h>
#include <vector>

using namespace std;

class LibChemist
{
public:
    class T{};
    int *hi;

    void AssignChemistValues(vector<double> vector, const double valor);
    double getTuputaMadre();
};


#endif	/* LIBCHEMIST_H */

