/* this is a generic soil cohort model to track organics */

#include "stdio.h"
#include "math.h"
#include "share_sedlib.h"
#include "fctns_sedlib.h"
#include "stdlib.h"

void sedlib_core_soil_cohort()

{
     int  i, j, k;

     double veg_source, veg_sink, veg_biomass_old;
     double grain_mass[50], grain_mass_old[50];
     double layer_mass[50], layer_mass_old[50];
     double new_bl_thick, total_bed_thick;
     double bl_prop, bl_prop_old, org_prop, net_org_prop;

     bl_prop = 0.0;
     bl_prop_old = 0.0;
     org_prop = 0.0;
     net_org_prop = 0.0;

/* update vegetation density */

     veg_source = (max_veg_source/rho) * MAX((limiting_veg_depth - MAX(depth_new,0.0)),0.0) / limiting_veg_depth;
     veg_sink = max_veg_source / max_eq_mass;

     veg_biomass_old = 0.0;
     for (k = 0; k < n_sed; k++)
      {
        if (k == roots_constituent) 
          {
          for (i = 0; i < n_blay-1; i ++)
           {
             veg_biomass_old += sed_sg[k] * (1. - bl_por[i]) * bl_thick[i] * bl_dist[i][k];
           }
           veg_biomass_old /= root_to_shoot_ratio;
          }
      }
     veg_biomass = (veg_biomass_old + veg_source*sed_time_step)/(1. + veg_sink*sed_time_step);
/* now update each layer.  */
/* dont include active layer: load veg mass in next bl */
  total_bed_thick = 0.0;
  for ( i = n_blay-2; i > -1; i -- )
    {
     layer_mass[i] = 0.0;
     layer_mass_old[i] = 0.0;
     for (k = 0; k < n_sed; k++)
      {
       grain_mass[k] = sed_sg[k] * (1. - bl_por[i]) * bl_thick[i] * bl_dist[i][k];
       grain_mass_old[k] = grain_mass[k];
       layer_mass_old[i] += grain_mass_old[k];
       if (k == roots_constituent)
         {
          grain_mass[k] = (1./.99) * veg_biomass * root_to_shoot_ratio * 
                                        (pow(0.01,MIN((total_bed_thick/root_thickness_limit),1.0)) -
                                         pow(0.01,MIN(((total_bed_thick + bl_thick[i])/root_thickness_limit),1.0)));
         }
       if (k == refractory_constituent)
         {
          grain_mass[k] += (1./.99) * veg_biomass * (1. + root_to_shoot_ratio) * veg_sink 
                                    * sed_time_step * (1. - labile_fraction) *  
                                        (pow(0.01,MIN((total_bed_thick/root_thickness_limit),1.0)) -
                                         pow(0.01,MIN(((total_bed_thick + bl_thick[i])/root_thickness_limit),1.0)));
         }
         layer_mass[i] += grain_mass[k];
       }
       total_bed_thick += bl_thick[i]; 
/* now update porosity and other erosion parameters */
    for (j =0; j < 4; j ++)
      {
      if(sed_type[k] == 1 || j == 0){
        if (j == 0) bl_prop_old = bl_por[i];
        if (j == 1) bl_prop_old = bl_ces[i];
        if (j == 2) bl_prop_old = bl_erc[i];
        if (j == 3) bl_prop_old = bl_ere[i];
        net_org_prop = 0.0;
        for (k = 0; k < n_sed; k++)
         {
           if (j == 0) org_prop = sed_por[k];
           if (j == 1) org_prop = sed_ces[k];
           if (j == 2) org_prop = sed_erc[k];
           if (j == 3) org_prop = sed_ere[k];
           {
             net_org_prop += org_prop * (grain_mass[k] - grain_mass_old[k]);
           }
         }
         if (layer_mass[i] > 1.E-16)
           bl_prop = (layer_mass_old[i]*bl_prop_old + net_org_prop) / layer_mass[i];
         if (j == 0) bl_por[i] = bl_prop;
         if (j == 1) bl_ces[i] = bl_prop;
         if (j == 2) bl_erc[i] = bl_prop;
         if (j == 3) bl_ere[i] = bl_prop;
        }
      }
/* now update bed layer thickness and distribution */
      new_bl_thick = 0.0;
      for (k = 0; k < n_sed; k++)
       {
         new_bl_thick += grain_mass[k]/(sed_sg[k]*(1. - bl_por[i]));
       }
      for (k = 0; k < n_sed; k++)
       {
         if (new_bl_thick > 1.E-16)
           bl_dist[i][k] = grain_mass[k]/(sed_sg[k]*(1. - bl_por[i])*new_bl_thick);
       }

       disp += (new_bl_thick - bl_thick[i]);
       as_ceiling += (new_bl_thick - bl_thick[i]);
       bl_thick[i] = new_bl_thick;
     }
/* now accumulate root and refractory biomass for output (alreadt calculated veg biomass above */
     root_biomass = 0.0;
     refractory_biomass = 0.0;
     for (k = 0; k < n_sed; k++)
      {
        for (i = 0; i < n_blay-1; i ++)
         {
           if (k == roots_constituent) root_biomass += rho * sed_sg[k] * (1. - bl_por[i]) * bl_thick[i] * bl_dist[i][k];
           if (k == refractory_constituent) refractory_biomass += rho * sed_sg[k] * (1. - bl_por[i]) * bl_thick[i] * bl_dist[i][k];
         }
      }
      veg_biomass *= rho;
}
