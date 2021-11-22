/******************************************************************************
 $Id: flowbeach.h,v 1.5 2004/09/23 12:41:36 schatz Exp $
 ******************************************************************************/

#ifndef STORM_H
#define STORM_H

#include "globals.h"
#include "func.h"
#include "PTG_FileNames.h"
#include "vec.h"

class dunepar;
class avalanche;

/*! Continuum aproach of sand relaxation by coasts (including coasts over solid substrates) */

class storm : public dunedata
{

public:
    // construction
    storm(const dunepar& p);
    virtual ~storm() {}
        
    void CalcGradUp(TFktScal &h);
    void Step(TFktScal &h, TFktScal &overwash);

    virtual double impact(double shoreline, TFktScal &h, TFktScal &h_nonerod, TFktScal &overwash, double* add_great2);
    virtual void stop( double time, double timestep, bool &calc_storm, double* add_great, double* add_small);

    virtual void calc( TFktScal &h, TFktScal &overwash );
  
    virtual void save_arrays();
        
private:
 
    /*!  Avalanche relaxation object.  */
    avalanche *m_avalanche;

    TFktScal m_sflux;
    TFktScal m_hst;
    
    int m_storm_iter;

    double m_Smax;
    double m_Sdt;
    double m_Q, m_scalefactor, m_shore_HMWL, m_watertable, m_surge, m_Tsurge;
    double m_slope;
    double m_storm_start;
    
    int m_shoreline;
    
    double surge[5000];
    int stormindex;
    int m_overwash;
    int m_seed;

    double m_frequency;
    double m_intensity;
    double poisson_param;
    double inter_event_time;
    double event_time;
    double event_tstep;
    double end_time;
    double rand_n;
    double event_times[5000];
};


#endif

