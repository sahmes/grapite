#ifndef GRAPE6_H
#define GRAPE6_H

#if defined(__cplusplus)
extern "C"
{
#endif
    /* C interface */

    /*
     * standard functions.
     * the number of the cards is hidden to the user.
     */
    void g6_open_all(void);
    void g6_close_all(void);
    int g6_set_j_particle_all(int address, int index, double tj, double dtj, double mass,
                              double a2by18[3], double a1by6[3], double aby2[3], double v[3], double x[3]);
    int g6_set_j_particle_mxonly_all(int address, int index, double mass, double x[3]);
    void g6_set_ti_all(double ti);
    void g6calc_firsthalf_all(int nj, int ni, int index[], double xi[][3], double vi[][3], 
                              double fold[][3], double jold[][3], double phiold[], double eps2, double h2[]);
    void g6calc_firsthalf0_all(int nj, int ni, int index[], double xi[][3], double vi[][3], 
                               double fold[][3], double jold[][3], double phiold[], double *eps2, double h2[], int mode);
    int g6calc_lasthalf_all(int nj, int ni, int index[], double xi[][3], double vi[][3], 
                            double eps2, double h2[], double acc[][3], double jerk[][3], double pot[]);
    int g6calc_lasthalf0_all(int nj, int ni, int index[], double xi[][3], double vi[][3], 
                             double *eps2, double h2[], double acc[][3], double jerk[][3], double pot[], int mode);
    int g6calc_lasthalf2_all(int nj, int ni, int index[], double xi[][3], double vi[][3],
                             double eps2, double h2[], double acc[][3], double jerk[][3], double pot[], int nnbindex[]);
    int g6_read_neighbour_list_all(void);
    int g6_get_neighbour_list_all(int ipipe, int maxlength, int *nblen, int nbl[]);
    void g6_set_nip_all(int nip);
    void g6_set_njp_all(int njp);
    void g6_set_i_particle_scales_from_real_value_all(int address, double acc[3], double jerk[3], double phi,
                                                      double jfactor, double ffactor);
    void g6_set_i_particle_all(int address, int index, double x[3], double v[3], double eps2, double h2);
    int g6_get_force_all(double acc[][3], double jerk[][3], double phi[], int flag[]);
    int g6_get_force_etc_all(double acc[][3], double jerk[][3], double phi[], int nnbindex[], int flag[]);
    void g6_get_predicted_j_particles_all(int nj, double (*x)[3], double (*v)[3]);
    int g6_getnjmax_all(void);

    /*
     * primitive functions to control multiple cards individually.
     * the user needs to specify card's device id explicitly.
     */
    void g6_open(int clusterid);
    void g6_close(int clusterid);
    void g6_set_tunit(int newtunit);
    void g6_set_xunit(int newxunit);
    int g6_set_j_particle(int clusterid, int address, int index, double tj, double dtj, double mass,
                          double a2by18[3], double a1by6[3], double aby2[3], double v[3], double x[3]);
    int g6_set_j_particle_mxonly(int  clusterid, int address, int index, double mass, double x[3]);
    void g6_set_ti(int clusterid, double ti);
    void g6calc_firsthalf(int clusterid, int nj, int ni, int index[], double xi[][3], double vi[][3], 
                          double fold[][3], double jold[][3], double phiold[], double eps2, double h2[]);
    void g6calc_firsthalf0(int clusterid, int nj, int ni, int index[], double xi[][3], double vi[][3], 
                           double fold[][3], double jold[][3], double phiold[], double *eps2, double h2[], int mode);
    int g6calc_lasthalf(int clusterid, int nj, int ni, int index[], double xi[][3], double vi[][3], 
                        double eps2, double h2[], double acc[][3], double jerk[][3], double pot[]);
    int g6calc_lasthalf0(int clusterid, int nj, int ni, int index[], double xi[][3], double vi[][3], 
                         double *eps2, double h2[], double acc[][3], double jerk[][3], double pot[], int mode);
    int g6calc_lasthalf2(int clusterid, int nj, int ni, int index[], double xi[][3], double vi[][3],
                         double eps2, double h2[], double acc[][3], double jerk[][3], double pot[], int nnbindex[]);
    int g6_read_neighbour_list(int clusterid);
    int g6_get_neighbour_list(int clusterid, int ipipe, int maxlength, int *nblen, int nbl[]);
    void g6_set_neighbour_list_sort_mode(int mode);
    int g6_get_neighbour_list_sort_mode(void);
    int g6_npipes(void);
    void g6_set_nip(int clusterid, int nip);
    void g6_set_njp(int clusterid, int njp);
    void g6_set_i_particle_scales_from_real_value(int clusterid, int address, double acc[3], double jerk[3], double phi,
                                                  double jfactor, double ffactor);
    void g6_set_i_particle(int clusterid, int address, int index, double x[3], double v[3], double eps2, double h2);
    int g6_get_force(int clusterid, double acc[][3], double jerk[][3], double phi[], int flag[]);
    int g6_get_force_etc(int clusterid, double acc[][3], double jerk[][3], double phi[], int nnbindex[], int flag[]);
    void g6_get_predicted_j_particles(int clusterid, int nj, double x[][3], double v[][3]);
    int g6_getnjmax(int clusterid);

    int grapite_tag_particles(int N, double x[][3], double v[][3], int skip, double fraction);

#if defined(__cplusplus)
}
#endif

#endif /* GRAPE6_H */
