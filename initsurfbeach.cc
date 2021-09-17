#include "initsurfbeach.h"
#include "globals.h"
#include "math.h"
////////////////////////////////////////
// class CInitSurfBeach
//

CInitSurfBeach::CInitSurfBeach(const dunepar& par, string prefix)
{
	m_h= par.getdefault<double>( prefix+"beach.h", duneglobals::dx());
    double m_slope= duneglobals::slope();//par.getdefault<double>( prefix+"beach.angle", 45);
    
    if (m_slope > 0) {
        m_l = m_h/(duneglobals::dx()*m_slope);
    }
}

void CInitSurfBeach::init_2d_scal(TFktScal& array)
{
    int x, y;

    for( x= 0; x< duneglobals::nx(); ++x )
        for( y= 0; y< duneglobals::ny(); ++y ){
            if (x < m_l) {
//                double a = 0*(x-0.5*m_l)*duneglobals::dx()/3.;
                array(x, y) =  m_h / m_l * x; // + 0.2*exp(-a*a); // m_h * (1 - (1-x*1.0/m_l)*(1-x*1.0/m_l));
            } else {
                array(x, y) = m_h;
            }
            //           double a = (x-0.5*m_l)*duneglobals::dx()/3.;
            //           array(x, y) = m_h * (1 - exp(- x /m_l)) + 0.2*exp(-a*a);
        }
}

/* test for oblique wind on rotated profile
void CInitSurfBeach::init_2d_scal(TFktScal& array)
{
    int x, y, theta, x_0, x_new, y_dim;

    theta = 45/360; // in radians
    x_0 = 20; // L_veg in m
    x_new = 32; // till when to incline wind
    y_dim = 32; // y limit size

    // to be commented before run
    for(x=0; x<static_cast<int>(duneglobals::nx()/20); ++x)
        for(y=0; y<static_cast<int>(x*tan(theta)); ++y){
            if (x < m_l) {
                array(x, y) = m_h / m_l * x;
            } else {
                array(x, y) = m_h;
            }

        }
    // to be commented before run
    for(x=0; x<static_cast<int>(duneglobals::nx()/20); ++x)
        for(y=static_cast<int>(x*tan(theta)); y<duneglobals::ny(); ++y){

                array(x, y) = 0;

        }

    for(x= 0; x< duneglobals::nx(); ++x )
        for( y= 0; y< duneglobals::ny(); ++y ){
            if (x < (x_0 - y*(x_0-x_new)/y_dim)) {
                array(x, y) = m_h/m_l*x/3;
            } else {
                array(x, y) = m_h;
            }
        }
}
*/

/* Part 2 - test initial profile
void CInitSurfBeach::init_2d_scal(TFktScal& array)
{
	int x, y;
	
	for( x= 0; x< duneglobals::nx(); ++x )
		for( y= 0; y< duneglobals::ny(); ++y ){
           if (x < m_l) {
                double a = 0*(x-0.5*m_l)*duneglobals::dx()/3.;
                array(x, y) = m_h * (1 - (1-x*1.0/m_l)*(1-x*1.0/m_l)); // m_h / m_l * x; // + 0.2*exp(-a*a);
            } else {
                array(x, y) = m_h;
            }
           double a = (x-0.5*m_l)*duneglobals::dx()/3.;
            array(x, y) = m_h * (1 - exp(- x /m_l)) + 0.2*exp(-a*a);
        }
}
*/

/* Original do not delete
void CInitSurfBeach::init_2d_scal(TFktScal& array)
{
    int x, y;

    for( x= 0; x< duneglobals::nx(); ++x )
        for( y= 0; y< duneglobals::ny(); ++y ){
            if (x < m_l) {
//                double a = 0*(x-0.5*m_l)*duneglobals::dx()/3.;
                array(x, y) =  m_h / m_l * x; // + 0.2*exp(-a*a); // m_h * (1 - (1-x*1.0/m_l)*(1-x*1.0/m_l));
            } else {
                array(x, y) = m_h;
            }
            //           double a = (x-0.5*m_l)*duneglobals::dx()/3.;
            //           array(x, y) = m_h * (1 - exp(- x /m_l)) + 0.2*exp(-a*a);
        }
}
*/