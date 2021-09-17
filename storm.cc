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
    
    m_shore_HMWL = duneglobals::HMWL();
    m_watertable = duneglobals::MSL();
    m_slope = duneglobals::slope(); //p.getdefault("beach.slope0", 0.5);

    m_Smax = p.getdefault("storm.Smax", 10000.0);
    
    // Frequency
    double m_freq = p.getdefault("storm.freq", 18.0);
    m_Sdt = duneglobals::secyear()*duneglobals::timefrac()/m_freq;
    
    cout << "!! STORM = " << m_Sdt << endl;

/*
    // Uncomment for marked Poisson HWE
    // Poisson distribution
    poisson_param = 7.8; // frequency in events/year
    inter_event_time = 0.0;
    event_time = 0.0; //42290000.0/365/24/60/60;  // change to make it storm time to start HWE after dune growth (in years)
    event_tstep = 0.0;
    end_time = 0.0;
    rand_n = 0.0;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    for(int i = 0; i< 5000000;i++){
        rand_n = distribution(generator);
        inter_event_time = -log(1.0-rand_n)/poisson_param;
        event_time = event_time + inter_event_time;
        event_tstep = event_time*365*24*60*60/1000;
        event_tstep = std::round(event_tstep);
        event_tsteps[i] = event_tstep;
    }


    // STORM INITIALIZATION
    /// ERLANG (GAMMA) DISTRIBUTION!!!

    double intensity = 1.0/0.353; // intensity in m

    std::default_random_engine generator2;
    std::exponential_distribution<double> distribution2(intensity);

    //int shape = 3;
    //double scalefactor0 = 0.4;
    //std::gamma_distribution<double> distribution2(shape,scalefactor0);

    for(int i = 0; i < 5000000; i++) {
        surge[i] = distribution2(generator2);
        // cout << surge[i] << "# StormT" << endl;
    }
*/


    // test with constant intensity
    double H_0 = 3.53; // in m
    double intensity = 0.5*H_0; // in m

    for(int i = 0; i < 10; i++) {
        surge[i] = intensity;
        //        cout << surge[i] << "# StormT" << endl;
    }

    stormindex = 0;
    
}

double storm::impact( double shoreline, TFktScal &h, TFktScal &h_nonerod, TFktScal &overwash, double* add_great2)
{
    // read the shoreline position
    m_shoreline = shoreline;
    
    overwash.SetAll(0.0);

    cout << "!! STORMINDEX = " << stormindex << " storm intensity " << surge[stormindex] << " \n" << endl;

    int iterstorm0 = 10;
    double isurge = m_scalefactor * surge[stormindex];//iterstorm0 - 100 + rand() % 300;

    *add_great2 = surge[stormindex];

    // avoid inundation
    m_Tsurge = m_shore_HMWL + (isurge < m_Smax ? isurge : m_Smax);

    double HMax = h.GetMax();

    // update to have erosion only when dune is overtopped (include if condition)
    //if(m_Tsurge>HMax){
        for (int i=0; i<iterstorm0; i++) {
            calc(h, overwash);
            // m_avalanche->calc(h, h_nonerod); // comment if no avalanche while iterating erosion
        }
    //}

    m_avalanche->calc(h, h_nonerod); // comment if no avalanche after erosion

    HMax = h.GetMax();

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

/*      // for Poisson storms
        bool exists = std::find(std::begin(event_tsteps), std::end(event_tsteps), tstep) != std::end(event_tsteps);
        calc_storm = exists;
*/

/*
        // for cyclic storms
        double Sstep = m_Sdt / timestep;
        int tnextstorm = (int) tstep % (int) Sstep; //750; //1000; //500;
        calc_storm = (tnextstorm == 0 && tstep >= 0 ? 1 : 0);
*/

        // for single fixed storm
        double Sstep = 10000;
        calc_storm = (((int) tstep % (int) Sstep == 0) && tstep >= 0 ? 1 : 0);

        if (calc_storm > 0) {
            stormindex++;
            cout << "!! STORM = " << stormindex << ' '  << tstep << " storm time step \n" << endl;

            //        cout << "!! STORM = " << calc_storm << ' ' << tstep << ' ' << stormindex << ' ' << tnextstorm << endl; // for cyclic storms

        }
        
    }
}

void storm::calc( TFktScal &h, TFktScal &overwash )
{
    for (int iter=0; iter < m_storm_iter * m_Tsurge * m_Tsurge; iter++) {
        Step(h, overwash);
    }
}

void storm::Step(TFktScal &h, TFktScal &overwash)
{
    double hnext, hi, hprev, Sfactor, hx, hxx, divq;
    double dx = duneglobals::dx();
       
    for (int y=0; y < duneglobals::ny(); y++) {
        bool cont = true;

        for (int x=m_shoreline; x < duneglobals::nx() && cont; x++) {

            hi = h(x, y);
            // definition of B.C
            // left: h = MHWL
            hprev = (x == m_shoreline ? m_shore_HMWL : h(x - 1, y));
            // right: depend on storm surge
            // surge > h -> h(x+1) = h(x) (hx = 0)
            hnext = (x == duneglobals::nx() - 1 ? m_shore_HMWL : h(x + 1, y));
            if (hi < m_Tsurge && h(x + 1, y) > m_Tsurge) {
                // surge < h -> h(x+1) = surge and STOP
                hnext = m_Tsurge;
                cont = false;
            }

            // Auxiliar
            hx = 0.5 * (hnext - hprev) / dx;
            hxx = (hnext - 2 * hi + hprev) / dx / dx;
            Sfactor = (m_Tsurge - hi);

            // Flux
            // in
            //   m_sflux(x,y) = (m_slope - hx) * Sfactor * Sfactor;
            // div Q
            divq = Sfactor * (hxx * Sfactor + 2 * (m_slope - hx) * hx);

            // Evol
            h(x, y) += m_Q * divq / m_Tsurge / m_Tsurge;

            // Submerge index
            overwash(x, y) = 1;


        }
    }

//    cout << "!! SURGE = " << m_Tsurge << ' ' << x << endl;

}

/*!  Saves the arrays m_u and m_rho.  */

void storm::save_arrays()
{
    //	save_2d_vecarray( "grad_h", m_grad_h_up );
//    save_2d_scalarray( "h_st", m_hst);
//    save_2d_scalarray( "", m_div_q );
//	save_2d_scalarray( "flux_beach", m_sflux );
}
