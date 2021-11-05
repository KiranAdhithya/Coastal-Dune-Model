/******************************************************************************
 $Id: storm.cc,v 1.9 2004/09/23 12:41:36 schatz Exp $
 ******************************************************************************/

#include <math.h>
#include <random>
#include <iostream>
#include <algorithm>

#include "globals.h"
#include "func.h"
#include "storm.h"
#include "avalanche.h"


//*****************************************************************************
//  class storm

storm::storm(const dunepar& p) : dunedata(p)
{
    m_avalanche = avalanche::create(p);
    
    m_sflux.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), 0.0);
    m_hst.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx(), 0.0);
    
    m_storm_iter = p.getdefault("storm.totaliter", 50);
    m_Q = p.getdefault("storm.step", 0.1);
    m_scalefactor = p.getdefault("storm.scalefactor", 1.0);

    m_overwash = p.getdefault("storm.overwash", 1);
    m_seed = p.getdefault("storm.seed", 1);
    m_frequency = p.getdefault("storm.frequency", 1.0);
    m_intensity = p.getdefault("storm.intensity", 0.35);

    m_Sdt = duneglobals::secyear()*duneglobals::timefrac()/m_frequency;

    cout << "!! STORM = " << m_Sdt << endl;

    m_shore_HMWL = duneglobals::HMWL();
    m_watertable = duneglobals::MSL();
    m_slope = duneglobals::slope(); //p.getdefault("beach.slope0", 0.5);

    m_Smax = p.getdefault("storm.Smax", 10000.0);
    
    /*
    // Frequency
    double m_freq = p.getdefault("storm.freq", 18.0);
    m_Sdt = duneglobals::secyear()*duneglobals::timefrac()/m_freq;
    
    cout << "!! STORM = " << m_Sdt << endl;
    */

    // Uncomment for marked Poisson HWE
    // Poisson distribution

//    poisson_param = m_frequency/m_Td; // frequency in events/year

    poisson_param = m_frequency; // frequency in events/year
    inter_event_time = 0.0;
    event_time = 0.0; //42290000.0/duneglobals::secyear()/duneglobals::timefrac();  // change to make it storm time to start HWE after dune growth (in years)
    rand_n = 0.0;
    //event_tstep = 0.0;
    //end_time = 0.0;

    unsigned seed = m_seed;

    std::default_random_engine generator (seed);
    std::uniform_real_distribution<double> distribution(0.0926,1.0); //0.0926 to get minimum of 2 days

    for(int i = 0; i< 1000;i++){

        rand_n = distribution(generator);
        inter_event_time = -log(1.0 - rand_n) / poisson_param;
        event_time = event_time + inter_event_time; //  adding two days as minimum inter-arrival time
        event_times[i] = event_time*365*24*60*60; //time in seconds
        //event_tstep = event_time*365*24*60*60/1000;
        //event_tstep = std::round(event_tstep);
    }

    // STORM INITIALIZATION
    /// ERLANG (GAMMA) DISTRIBUTION!!!

    double intensity = 1.0/m_intensity; // Exponential parameter in m^(-1)

    unsigned seed2 = m_seed;

    std::default_random_engine generator2 (seed2);
    std::exponential_distribution<double> distribution2(intensity);

    //int shape = 3;
    //double scalefactor0 = 0.4;
    //std::gamma_distribution<double> distribution2(shape,scalefactor0);

    for(int i = 0; i < 1000; i++) {
        surge[i] = distribution2(generator2);
    }

/*    // test with constant intensity
    double H_0 = 3.9504; //6.5402; //4.6849; //3.9504;  // in m
    double intensity = H_0; // in m
    //double intens[10] = {0.7};

    for(int i = 0; i < 10; i++) {

        //        cout << surge[i] << "# StormT" << endl;
        surge[i] = (10-i)*intensity/10;
    }
*/
    stormindex = 0;
    
}

double storm::impact( double shoreline, TFktScal &h, TFktScal &h_nonerod, TFktScal &overwash, double* add_great2)
{
    // read the shoreline position
    m_shoreline = shoreline;
    
    overwash.SetAll(0.0);

    cout << "!! STORMINDEX = " << stormindex << " storm intensity " << surge[stormindex] << " \n" << endl;

    int iterstorm0 = 1; //10 // iteration is done for wet avalanche stability
    double isurge = m_scalefactor * surge[stormindex];//iterstorm0 - 100 + rand() % 300;

    *add_great2 = surge[stormindex];

    // avoid inundation
    m_Tsurge = m_shore_HMWL + (isurge < m_Smax ? isurge : m_Smax);

    //save_2d_scalarray( "h_prestorm", h);

    //for (int i=0; i<iterstorm0; i++) {
            calc(h, overwash);
            //m_avalanche->calc(h, h_nonerod); // comment if no avalanche while iterating erosion
    //}

    // save_2d_scalarray( "h_aval", h);

    m_avalanche->calc(h, h_nonerod); // comment if no avalanche after erosion

    //save_2d_scalarray( "h_poststorm", h);

    double HMax = h.GetMax();

    if (HMax < m_shore_HMWL + 0.3)
    {
        m_Tsurge *= -1.0;
    }

    return m_Tsurge;
}

void storm::stop( double time, double timestep, bool &calc_storm, double* add_great, double* add_small)
{
    
    if( calc_storm > 0 ){
        calc_storm = 0;
        *add_great = surge[stormindex];
        *add_small = time;

    } else {
        double tstep = time / timestep;

        // for Poisson storms
        //bool exists = std::find(std::begin(event_tsteps), std::end(event_tsteps), tstep) != std::end(event_tsteps);

        calc_storm = (time >= event_times[stormindex+1] ? 1:0);

/*        // for cyclic storms
        double Sstep = m_Sdt / timestep;
        int tnextstorm = (int) tstep % (int) Sstep; //750; //1000; //500;
        calc_storm = (tnextstorm == 0 && tstep >= 0 ? 1 : 0);
*/


/*       // for single fixed storm
        double Sstep = 120000;
        calc_storm = (((int) tstep % (int) Sstep == 0) && tstep >= 0 ? 1 : 0);
*/
        if (calc_storm > 0) {
            stormindex++;
            // cout << "!! STORM = " << stormindex << ' '  << tstep << " storm time step \n" << endl;

            //        cout << "!! STORM = " << calc_storm << ' ' << tstep << ' ' << stormindex << ' ' << tnextstorm << endl; // for cyclic storms

        }
        
    }
}

void storm::calc( TFktScal &h, TFktScal &overwash )
{
    //for (int iter=0; iter < 1; iter++) { // iter<factor*m_storm_iter * m_Tsurge * m_Tsurge; factor 1,100/200 to increase storm erosion time
        Step(h, overwash);
    //}
    //cout << "!! SLOPE = " << m_slope << endl;
}

void storm::Step(TFktScal &h, TFktScal &overwash)
{
    //double hnext, hi, hprev, Sfactor, hx, hxx, divq,q,dx = duneglobals::dx();
    double hi;
    int max_x,x;


    for (int y=0; y < duneglobals::ny(); y++) {
        bool cont = true;

        double hmax=0.0;

        for (x =m_shoreline+1; x<duneglobals::nx();x++){
            if(h(x,y)>hmax){
                hmax = h(x,y);
                max_x = x;
            }
        }

        if(hmax<m_Tsurge){
            for (x =m_shoreline+1; x<duneglobals::nx();x++){

                h(x,y) = m_shore_HMWL;
                overwash(x, y) = 1;
            }
        } else if(m_overwash == 0){

            for (x =m_shoreline+1; x<max_x && cont;x++) { // start from next point to shoreline // x = m_shoreline +1;

                hi = h(x, y);
                if (hi < m_Tsurge) {
                    // surge < h -> h(x+1) = surge and STOP
                    h(x, y) = m_shore_HMWL;
                    overwash(x, y) = 1;
                    }
                else {
                    cont = false;
                }
            }
        }
    }

    //    cout << "!! SURGE = " << m_Tsurge << ' ' << x << endl;

}

/*!  Saves the arrays m_u and m_rho.  */

void storm::save_arrays()
{
    //	save_2d_vecarray( "grad_h", m_grad_h_up );
 //   save_2d_scalarray( "h_st", m_hst);
//    save_2d_scalarray( "", m_div_q );
//	save_2d_scalarray( "flux_beach", m_sflux );
}
