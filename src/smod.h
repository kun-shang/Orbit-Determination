
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

#ifndef _SMOD_H_
    #define _SMOD_H_

    #ifndef REAL128
       #define REAL128
       typedef long double real128;
    #endif

    #ifndef __STDIO_H__
        #include <stdio.h>
    #endif

    #ifdef MSDOS
        #ifndef __CONIO_H__
            #include <conio.h>
        #endif
    #endif

    #ifdef LINUX
        void getch(void) {}
    #endif


    #ifndef __MATH_H__
        #include <math.h>
    #endif

    #ifndef __STRING_H__
        #include <string.h>
    #endif

    #ifndef __STDLIB_H__
        #include <stdlib.h>
    #endif

    #ifndef __CTYPE_H__
        #include <ctype.h>
    #endif

    #ifndef __TIME_H__
        #include <time.h>
    #endif

    #ifndef _NOVAS_H_
        #include "novas.h"
    #endif

    #ifndef __EPHMAN_H_
        #include "eph_manager.h"
    #endif


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


    typedef struct InfStruct
    {
        int leaps;
        double jd0;
        double jdt;
        double mjd;
        double gmst;
        double deltat;
        double utc;
        double ut1;
        double tt;
        double xp;
        double yp;
        double ut1_utc;
        double dx;
        double dy;
        double c_ie[9];
        double c_ei[9];
    }InfStruct;

    typedef struct OTStruct
    {
        double ds;
        int argn[6];
        char name;
        int n;
        int m;
        double cp;
        double sp;
        double cm;
        double sm;
    }OTStruct;


    #define DYNPAR 4
    
    short int CENTER, GRAVDEGREE, GRAVORDER,
        ACCURACY, LEAPSECS, PERB[11];
    int DIM_TE, DIM_OR;
    double JD0, STEP_OR, STEP_TE, *TE_EPH, *OR_EPH, GMDE[11], 
        SRPJ, ASRP, CONS, DCON, BIAS, DBIA, SATMASS, SATAREA;
    char FILE_GRV[80], FILE_EOP[80];
    double *cn0, *cnm, *snm, gma[2];
    double j2, j3, j4, jc21, js21, jc22, js22;

    OTStruct *otfes;

    int PERMT, STIDE, OTIDE, NBODY, RELTIV, NFES, NMAX, MMAX, NSMAX, NOMAX, NEOP,
        CT;
    double AMR, *COEFG, *COEFS, *COEFO, *EOPMT, GMCT, RCT;

    int IBASB, IBAST, ISRPB, ISRPT, IK2;
    
    double SRPB, SRPT, BASB, BAST, K2;
    
    int MEST, SLOVEFOR, MSTA, MOBS, MSRP, MGCS, MTK2, MDYN, MSOL; 
    
    #define I_TITAN 20

    typedef struct CSStruct
    {
        int n;
        int m;
        int cs;
        double initv;
        double dadcsn[3];
        double dadcse[3];
        double dadcs[3];
    }CSStruct;


    CSStruct CSinfo[1000];
#ifdef _cplusplus
extern "C"
{
#endif

extern void DOODSN (double *jd_tdb, double *gmst, int *argn, int *ncon,
        double *doodarg, double *ang);
#define DOODSN doodsn_
#ifdef _cplusplus
}
#endif


double obs_vel (double jd, double utc, double *obs, int part, double *bmat);
void initsolvefor (double *xsm, double *x);
void updsolvefor (double *x);
void getsolvefor ();

double stidecs_k2(InfStruct *info, double k2, double *stcs, int *body, int nbody);
double cs2ada (double *llr, double *cs, int nmax, double *acc,
                int part, double *dadr, int flag);
int brank(double *a, int m, int n);

void choldc(double *a, int n, double p[]);

void cholsl(double *a, int n, double p[], double b[], double x[]);

void solvels_chol(double *a, int n, double *y, double *x, int nocov);

    void solvegaus(double *a, int n, double *y, double *x);

    void pt_orb (double ts_orb, double te_orb, double step_orb, int dim_eph);

    double obs_alt (double jd, double utc, double *obs, int part, double *bmat);
    double obs_dsn (double jd, double utc, double *obs, int part, double *bmat);

    int get_ephemeris (double tjd[2], int to, int from, double *x);

    double accel_pm_part (double *tjd, double *x, double *fnt, 
        int part, double *dadr);



    double accel_nb_part (double *tjd, double *xic, double *fnt, 
        int part, double *dadr);

    double f_bcrs (double *jd, double *xi, int exclude, 
                   double *acc, int part, double *dadr);

    double accel_sr_part (double *tjd, double *xic, double *acc,
        int part, double *dadr, double *dadpb, double *dadpt);

    double accel_gt_part (double *tjd, double *xic, double *ag, 
        int part, double *dadr, double *dadk2);



    double mgrn1(double u,double g,double *r);

    double stfrqdep(double jdt, double gmst, double *c20f, double *c21f, double *s21f, double *c22f, double *s22f);

    int openotcs (char *infile);

    double otidecs(double jdt, double gmst, int nmax, double *coef);
    
    double getinfo(double *tjd, InfStruct *info);

//    double stidecs_Anelastic(double *tjd, double *c_ie, int id_perm, double *stcs);
    double stidecs_Anelastic(InfStruct *info, int id_perm, double *stcs);

    double lgdr(double t, int nmax, int m, double *pbar);

    double accel_gravtide (double *tjd, double *xic, 
            double *ag, double *at, double *ao);

    double stidecs(double *tjd, double *c_ie, int id_perm, double *stcs);

    double accel_pmiers (double *tjd, double *x, double *fnt, double *fgr);

    void aei2xyz (double ele[6], double pos[3], double vel[3]);

    double kepler(double M,double e);

    void xyz2llh (double *vt, double *llh);

    double earth_fullaccel (double jd, double tt, double *xic, double *fxic);

    double accel_point (double *tjd, double *x, double *fnt, double *fgr);

    double accel_slrad (double *tjd, double *xic, double *acc);

    double accel_nbody (double *tjd, double *xic, double *fnt, double *fgr);

    double force_bcrs (double *jd, double *xi, short int exclude, 
                   double *fnt, double *fgr);

    void mt (double *a, int m, int n, double *b);

    void brmul (double *a, double *b, int m,int n, int k,double *c);

    double rkf78 (double jd, double t, double h, double *x, int dim, 
              void (*fun)(int, double, double,double *,double *));
//              double (*fun)( double, double,double *,double *));
    
    void fun_accel (int dim, double jd, double tt, double *xic, double *fxic);

    double rkf78_auto (double h, double t, double *x, int dim, double err, 
              double (*fun)(double,double *,double *), int autoadjust);

    double opengrv (char file_grv[2][200], double *coef, int nmax, int mmax);

    double earth_pointmass (double jd, double tdbs, double *x, double *f);

    double modvect (double *v);

    double dotvect (double *v1, double *v2);

    void crsvect (double *v1, double *v2, double *v);

    double chosephase (double sinvalue, double cosvalue);

    double xyz2aei(double ele[6], double pos[3], double vel[3]);

    double accel_gravt (double *tjd, double *xic, double *a4);

    double cs2acc (double *llr, double *cs, double gm, double a, int nmax, 
               double *acc);

    double lgdr2(double t, int nmax, int m, 
             double *pbar, double *pbar1, double *pbar2);

    void openeop (char file_eop[200], int mjds, int num, double *eopmat);

    void geteop (double mjd, double *xp, double *yp, 
               double *ut1_utc, double *dx, double *dy);

    void itrf2gcrf(double jd, double utc, double *vt, double *vc);

    void gcrf2itrf(double jd, double utc, double *vc, double *vt);

    int getlps (double jdutc);



/* Reference epoch (J2000), Julian Date */
    #define DJ00 (2451545.0)
    
/* Days per Julian millennium */
    #define DJM (365250.0)
    
/* Degrees to radians */
    #define DD2R (1.745329251994329576923691e-2)
    
/* 2Pi */
    #define D2PI (6.283185307179586476925287)


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


/*
#ifdef _cplusplus
extern "C" 
{
#endif 
	void _stdcall DPLEPH (double *tjdtdb, int *targe, int *center,
		double *posvel);   
#ifdef _cplusplus
}
#endif
*/


void xyz2rtn(double *x, double *v, double *xyz, double *rtn);

double stidecs_old(double *tjd, double gma1, double k2, 
               double *c20, double *c21, double *s21, double *c22, double *s22);


double al_part (double *xc2, double *dodx, double *dodp);


double simula_altim (double tdb, double *calculable, short int part, double *bmat);

void azelev (double jd_ut1, double delta_t, short int accuracy,
             double x, double y, double *llh, double ra,
             double dec, double *zd, double *az);

void xyz2llh (double *vt, double *llh);

double lagrange (double *y, int dim_y, int dim_x, double t, double *z);

double delta_iid (double *jd, double *xi, double *ii, double *id);

double delta_tdb (double *txice, double *txics, double *deltat);

double fun_fullaccel (double tdbs, double *xic, double *fxic);

double fun_fullstate (double tdbs, double *state, double *fstate);

double accel_bcrs (double *jd, double *xi, short int part, short int exclude, 
                   double *acc, double *dadr, double *dadp);

double accel_ntrel (double *tjd, double *xic, short int part, 
                    double *acc, double *dadr, double *dadp);

double accel_nonsp (double *tjd, double *xic, short int part, 
                    double *acc, double *dadr, double *dadp);

double accel_radpr (double *tjd, double *xic, short int part, 
                    double *acc, double *dadr, double *dadp);

void in2pa (double *jd, double *te);


double itrf2icrf (double jd, double utc, double *vt, double *vc);

double fun_pointmass (double tdbs, double *x, double *f);


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double simula_phase (double utc3, double utc0, double *station3, 
                     short int uplink, double *station1, short int genrel,
                     double *calculable, double *azimuth, double *elevation, 
                     short int part, double *bmat);

double simula_dople (double utc3, double tc, double *station3, 
                     short int uplink, double *station1, short int genrel,
                     double *calculable, double *azimuth, double *elevation, 
                     short int part, double *bmat);


double simula_range (double utc3, double *station3, short int uplink, 
                     double *station1, short int genrel, 
                     double *calculable, double *azimuth, double *elevation, 
                     short int part, double *bmat);

double ltsolution (double utc_3, double *station3, short int uplink,  
                   double *station1, short int genrel, real128 *lt, 
                   double *azimuth, double *elevation, short int part,
                   double *bmat, double *txic);

real128 lt_form (double tdb3, double tdb2, double *re3, double *rp2, 
                   int genrel, real128 *rs3, real128 *rs2);

double lt_part (real128 *rs3, real128 *rs2, real128 *rs1, int uplink, 
                  double *dodx, double *dodp);

double opengravfile (double *cn0, double *cnm, double *snm, double *gma);


void iau_pns (double *jd, double *te, int cent);

void iau_s (double *jd, double *tes, int cent);

void iau_pn (double *jd, double *tes, int cent);

short int mbf2cel (double *jd_tdb, double *te);

void in2me (double *jd, double *te, short int derivation);

void me2pa (double *te);

void rotmatx (double rad, double *matx, short int deri);

void rotmaty (double rad, double *maty, short int deri);

void rotmatz (double rad, double *matz, short int deri);

double enlgr (double *x, double *y, int n, double t);



int bssgj (double *a,int n);

int brinv (double *a,int n);

double iauDtdb(double date1, double date2,
               double ut, double elong, double u, double v);





#endif


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

