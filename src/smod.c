
#ifndef _SMOD_H_
    #include "smod.h"
#endif



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* brinv - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int brinv (double *a,int n)
{ 
    int *is, *js, i, j, k, l, u, v, rk;
    double d, p;

    is = (int *)malloc (n * sizeof(int));
    js = (int *)malloc (n * sizeof(int));
    for (k = 0; k <= n - 1; k++)
    { 
        d=0.0;

        for (i = k; i <= n - 1; i++)
        {
            for (j = k; j <= n - 1; j++)
            { 
                l = i * n + j; 
                p = fabs (a[l]);
                if (p > d) 
                { 
                    d = p; 
                    is[k] = i; 
                    js[k] = j;
                }
            }
        }

        if (d + 1.0 == 1.0)
        { 
            free (is); 
            free (js); 
            rk = brank(a, n, n);
            printf ("error: Matrix is ill-conditioned:\n");
            printf ("      not inv\n");
            printf ("       rank = %d\n", rk);
            exit(0);
        }

        if (is[k] != k)
        {
            for (j = 0; j <= n - 1; j++)
            { 
                u = k * n + j; 
                v = is[k] * n + j;
                p = a[u];
                a[u] = a[v];
                a[v] = p;
            }
        }

        if (js[k] != k)
        {
            for (i = 0; i <= n - 1; i++)
            { 
                u = i * n + k; 
                v = i * n + js[k];
                p = a[u];
                a[u] = a[v];
                a[v] = p;
            }
        }

        l = k * n + k;
        a[l] = 1.0 / a[l];
        
        for (j = 0; j <= n - 1; j++)
        {
            if (j != k)
            { 
                u = k * n + j; 
                a[u] = a[u] * a[l];
            }
        }
        
        for (i = 0; i <= n - 1; i++)
        {
            if (i != k)
                for (j = 0; j <= n - 1; j++)
                    if (j != k)
                    { 
                        u = i * n + j;
                        a[u] = a[u] - a[i * n + k] * a[k * n + j];
                    }
        }

        for (i = 0; i <= n - 1; i++)
        {
            if (i != k)
            { 
                u = i * n + k; 
                a[u] = - a[u] * a[l];
            }
        }
    }
    
    for (k = n - 1; k >= 0; k--)
    { 
        if (js[k] != k)
        {
            for (j = 0; j <= n - 1; j++)
            { 
                u = k * n + j; 
                v = js[k] * n + j;
                p = a[u]; 
                a[u] = a[v]; 
                a[v] = p;
            }
        }
        
        if (is[k] != k)
        {
            for (i = 0; i <= n - 1; i++)
            { 
                u = i * n + k; 
                v = i * n + is[k];
                p = a[u]; 
                a[u] = a[v]; 
                a[v] = p;
            }
        }
    }
    
    free(is); 
    free(js);
    return(1);
}







int brank(double *a, int m, int n)
  { int i,j,k,nn,is,js,l,ll,u,v;
    double q,d;
    nn=m;
    if (m>=n) nn=n;
    k=0;
    for (l=0; l<=nn-1; l++)
      { q=0.0;
        for (i=l; i<=m-1; i++)
        for (j=l; j<=n-1; j++)
          { ll=i*n+j; d=fabs(a[ll]);
        if (d>q) { q=d; is=i; js=j;}
          }
        if (q+1.0==1.0) return(k);
        k=k+1;
        if (is!=l)
          { for (j=l; j<=n-1; j++)
              { u=l*n+j; v=is*n+j;
                d=a[u]; a[u]=a[v]; a[v]=d;
              }
          }
        if (js!=l)
          { for (i=l; i<=m-1; i++)
              { u=i*n+js; v=i*n+l;
                d=a[u]; a[u]=a[v]; a[v]=d;
              }
          }
        ll=l*n+l;
        for (i=l+1; i<=n-1; i++)
          { d=a[i*n+l]/a[ll];
            for (j=l+1; j<=n-1; j++)
              { u=i*n+j;
                a[u]=a[u]-d*a[l*n+j];
              }
          }
      }
    return(k);
  }



void choldc(double *a, int n, double p[])
/*
 * Given a positive-definite symmetric matrix a[1..n][1..n], 
 * this routine constructs its Cholesky decomposition, A = L · LT . 
 * On input, only the upper triangle of a need be given; it is not modified. 
 * The Cholesky factor L is returned in the lower triangle of a, 
 * except for its diagonal elements which are returned in p[1..n].
 *
*/
{
    int i,j,k, rk;
    double sum;

    for (i=0;i<n;i++) 
    {
        for (j=i;j<n;j++) 
        {
            sum=a[i * n + j];
            for (k=i-1;k>=0;k--) 
                sum -= a[i * n + k]*a[j * n + k];
            if (i == j) 
            {
//                printf("i = %d\tj = %d\tsum = %e\n", i, j, sqrt(sum));
                if (sum <= 0.0)
                {
                    rk = brank(a, n, n);
                    printf("error: Matrix is not positive definite:\n");
                    printf("       i = %d\tj = %d\tsum = %e\n", i, j, sum);
                    printf("       rank = %d\n", rk);
                    exit(0);
                }
                p[i]=sqrt(sum);
            } 
            else a[j * n + i]=sum/p[i];
        }
    }


}


void cholsl(double *a, int n, double p[], double b[], double x[])
/*
 * Solves the set of n linear equations A · x = b, 
 * where a is a positive-definite symmetric matrix. 
 * a[1..n][1..n] and p[1..n] are input as the output of the routine choldc. 
 * Only the lower subdiagonal portion of a is accessed. 
 * b[1..n] is input as the right-hand side vector. 
 * The solution vector is returned in x[1..n]. 
 * a, n, and p are not modified and can be left in place for 
 * successive calls with different right-hand sides b. 
 * b is not modified unless you identify b and x in the calling sequence,
 * which is allowed.
 *
*/
{
    int i,k;
    double sum;

    for (i=0;i<n;i++) 
    {
        for (sum=b[i],k=i-1;k>=0;k--) 
            sum -= a[i * n + k]*x[k];
        x[i]=sum/p[i];
    }
    for (i=n-1;i>=0;i--) 
    {
        for (sum=x[i],k=i+1;k<n;k++) 
            sum -= a[k * n + i]*x[k];
        x[i]=sum/p[i];
    }


}






///////////********************////////////////
void solvels_chol(double *a, int n, double *y, double *x, int nocov)
{
    double *p, sum;
    int i, k, j;

    p = (double *)calloc (n, sizeof(double));

    choldc(a, n, p);
    cholsl(a, n, p, y, x);

    if (nocov == 1) 
    {
        free (p);
        return;
    }

    for (i=0;i<n;i++) 
    { 
        a[i * n + i]=1.0/p[i];
        for (j=i+1;j<n;j++) 
        {
            sum=0.0;
            for (k=i;k<j;k++) 
                sum -= a[j * n + k]*a[k * n + i]; 
            a[j * n + i]=sum/p[j];
        } 
    }
    
    for (i = 0; i <= n - 1; i++)
    {
        for (j = i; j <= n - 1; j++)
        {
            sum = 0.0;
            for (k = j; k <= n - 1; k++)
                sum = sum + a[k * n + i] * a[k * n + j];
            a[i * n + j] = sum;
        }
    }

    for (i = 0; i <= n - 1; i++)
    {
        for (j = 0; j <= i - 1; j++)
        {
            a[i * n + j] = 0;
        }
    }

    free(p);
    return;
}    















void solvegaus(double *a, int n, double *y, double *x)
{
    brinv(a, n);
    brmul(a, y, n, n, 1, x);
}












void pt_orb (double ts_orb, double te_orb, double step_orb, int dim_eph)
{
    FILE *fp_fxyz, *fp_faei, *fp_frtn, *fp_fllh;
    int i, n;
    double tt, lps, utc, xtm[6], *eph, dist, velt, tp, rtn_p[3], rtn_v[3], 
        ele[6], llh[3];

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    /*--print orbit--*/
    if((fp_fxyz=fopen("forb.xyz","w"))==NULL)
    {
        printf("Cannot write fort.xyz!\n");
      	getch();
        exit(0);
    }

    if((fp_faei=fopen("forb.aei","w"))==NULL)
    {
        printf("Cannot write fort.aei!\n");
      	getch();
        exit(0);
    }

    if((fp_frtn=fopen("forb.rtn","w"))==NULL)
    {
        printf("Cannot write fort.rtn!\n");
      	getch();
        exit(0);
    }
    if((fp_fllh=fopen("forb.llh","w"))==NULL)
    {
        printf("Cannot write fort.llh!\n");
      	getch();
        exit(0);
    }

        
    eph = (double *) calloc (dim_eph - 1, sizeof(double));

    i = 0;
    for (utc = ts_orb; utc <= te_orb; utc = ts_orb + step_orb * i)
    {
        lps = getlps (JD0 + ts_orb/86400.0);
        tt = utc + (lps + 32.184);

        lagrange (OR_EPH, DIM_OR, dim_eph, tt, eph);
        for (n = 0; n < 6; n++)
            xtm[n] = eph[n];


        dist = sqrt (xtm[0] * xtm[0] + xtm[1] * xtm[1] + xtm[2] * xtm[2]);
        velt = sqrt (xtm[3] * xtm[3] + xtm[4] * xtm[4] + xtm[5] * xtm[5]);

        tp = xyz2aei(ele, &xtm[0], &xtm[3]);

        xyz2rtn(&xtm[0], &xtm[3], &xtm[0], rtn_p);
        xyz2rtn(&xtm[0], &xtm[3], &xtm[3], rtn_v);
    
        xyz2llh(xtm, llh);

        fprintf (fp_fxyz, "%14.4f  %14.6f  %26.14f  %26.14f  %26.14f  ",
            JD0, utc, xtm[0], xtm[1], xtm[2]);
        fprintf (fp_fxyz, "%24.16f  %24.16f  %24.16f  %16.4f  %14.6f \n", 
            xtm[3], xtm[4], xtm[5], dist, velt);

        fprintf (fp_faei, "%14.4f  %14.6f  %26.14f  %10.6f  %12.6f  ",
            JD0, utc, ele[0], ele[1], ele[2]);
        fprintf (fp_faei, "%12.6f  %12.6f  %12.4f  %12.4f \n", 
            ele[3], ele[4], ele[5], tp);

        fprintf (fp_frtn, "%14.4f  %14.6f  %16.4f  %16.4f  %16.4f  ",
            JD0, utc, rtn_p[0], rtn_p[1], rtn_p[2]);
        fprintf (fp_frtn, "%14.6f  %14.6f  %14.6f  %16.4f  %14.6f \n", 
            rtn_v[0], rtn_v[1], rtn_v[2], dist, velt);

        fprintf (fp_fllh, "%14.4f  %14.6f  %26.14f  %26.14f  %26.14f\n",
            JD0, utc, llh[0], llh[1], llh[2] - RCT);
        i++;

    }
    fclose(fp_fxyz);
    fclose(fp_faei);
    fclose(fp_frtn);
    fclose(fp_fllh);
    free (eph);

    return;
}





  double mgrn1(double u,double g,double *r)
  { int i,m;
    double s,w,v,t;
    s=65536.0; w=2053.0; v=13849.0;
    t=0.0;
    for (i=1; i<=12; i++)
      { *r=(*r)*w+v; m=(int)(*r/s);
        *r=*r-m*s; t=t+(*r)/s;
      }
    t=u+g*(t-6.0);
    return(t);
  }


    

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
 * geteop - interpolation of eop
 * mjd: double, input MJD
 * xp, yp, ut1_utc, dx, dy: output EOP
 
 * http://hpiers.obspm.fr/iers/eop/eopc04_05/eopc04_IAU2000.62-now
 * version: 20 Aug 2010
 */
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void geteop (double mjd, double *xp, double *yp, 
               double *ut1_utc, double *dx, double *dy)
{
    int i, mjdi;
    double x1, y1, dt1, dx1, dy1, x2, y2, dt2, dx2, dy2; 

            
    for (i = 0; i < NEOP; i ++)
    {
        mjdi = (int)EOPMT[i * 6 + 0];

        if (mjdi == (int)mjd)
        {
            x1  = EOPMT[i * 6 + 1];
            y1  = EOPMT[i * 6 + 2];
            dt1 = EOPMT[i * 6 + 3];
            dx1 = EOPMT[i * 6 + 4];
            dy1 = EOPMT[i * 6 + 5];
            i++;
            x2  = EOPMT[i * 6 + 1];
            y2  = EOPMT[i * 6 + 2];
            dt2 = EOPMT[i * 6 + 3];
            dx2 = EOPMT[i * 6 + 4];
            dy2 = EOPMT[i * 6 + 5];

            break;
        }
    }    

    *xp      = x1 + (x2-x1) * (mjd - mjdi);
    *yp      = y1 + (y2-y1) * (mjd - mjdi);
    *ut1_utc = dt1 + (dt2-dt1) * (mjd - mjdi);
    *dx      = dx1 + (dx2-dx1) * (mjd - mjdi);
    *dy      = dy1 + (dy2-dy1) * (mjd - mjdi);

}





double getinfo(double *tjd, InfStruct *info)
{
    int n;
    double gmsth, ux[3] = {1,0,0}, uy[3] = {0,1,0}, uz[3] = {0,0,1}, 
        tx[3], ty[3], tz[3], xp, yp, ut1_utc, dx, dy;


    info->jd0 = tjd[0];
    info->tt = tjd[1] * 86400.0;

    info->jdt = info->jd0 + info->tt / 86400.0;
        
    info->leaps = getlps (info->jdt);

    info->utc = info->tt - info->leaps - 32.184;

    if (CT != 2)
    {
        iau_pns(tjd, info->c_ei, CT);
        mt(info->c_ei, 3, 3, info->c_ie);     
        return 1;
    }

    info->mjd = info->jd0 - 2400000.5 + info->utc/86400.0;
    geteop (info->mjd, &xp, &yp, &ut1_utc, &dx, &dy);	

    info->xp      = xp;
    info->yp      = yp;
    info->ut1_utc = ut1_utc;
    info->dx      = dx;
    info->dy      = dy;


    info->deltat = 32.184 + info->leaps - info->ut1_utc;
    info->ut1 = info->utc + info->ut1_utc;

    sidereal_time (info->jd0, info->ut1/86400.0, info->deltat,0,1,1, &gmsth);

    info->gmst = gmsth / 24 * 360.0 * DEG2RAD;


    cel_pole (info->jdt, 2, info->dx * 1e3, info->dy * 1e3);

    cel2ter (info->jd0, info->ut1 / 86400.0, info->deltat, 1, 1, 0,
        info->xp, info->yp, ux, tx);
    cel2ter (info->jd0, info->ut1 / 86400.0, info->deltat, 1, 1, 0,
        info->xp, info->yp, uy, ty);
    cel2ter (info->jd0, info->ut1 / 86400.0, info->deltat, 1, 1, 0,
        info->xp, info->yp, uz, tz);

    for (n = 0; n < 3; n++)
    {
        info->c_ie[n*3] = tx[n];
        info->c_ie[n*3+1] = ty[n];
        info->c_ie[n*3+2] = tz[n];
    }
    
    mt(info->c_ie, 3, 3, info->c_ei);

//        printf ("%d\tjd0 = %.10f\t ut1 = %.10f\t utc = %.10f\n", i,info[i].jd0, info[i].ut1, info[i].utc );
//        for (n = 0; n < 9; n++)
//            printf ("%e\n", info[i].c_ie[n]);


    return 0;
}












/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* itrf2gcrf(icrf) - from earth fixed to earth inertial
* jd: double, integral part of JD day, unit: day
* utc: double, fractional part of JD day, unit: seconds
* @param2: description of param2 
        
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void itrf2gcrf(double jd, double utc, double *vt, double *vc)
{
    double lps, mjd, xp, yp, ut1_utc, dx, dy, delta_t, ut1, tt;

    mjd = jd - 2400000.5;
    geteop (mjd + utc/86400.0, &xp, &yp, &ut1_utc, &dx, &dy);	
    lps = getlps (jd + utc/86400.0);
    delta_t = 32.184 + lps - ut1_utc;
    ut1 = utc + ut1_utc;
    tt = utc + (lps + 32.184); 	

    cel_pole (jd + tt / 86400.0, 2, dx * 1e3, dy * 1e3);

    ter2cel (jd, ut1 / 86400.0, delta_t, 1, 1, 0,
        xp, yp, vt, vc); /*--vc unit: m--*/    

}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
 * * gcrf(icrf)2itrf - from earth fixed to earth inertial
 * * jd: double, integral part of JD day, unit: day
 * * utc: double, fractional part of JD day, unit: seconds
 * * @param2: description of param2 
 *         
 *         * version: 20 Aug 2010
 *         */
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void gcrf2itrf(double jd, double utc, double *vc, double *vt)
{
    double lps, mjd, xp, yp, ut1_utc, dx, dy, delta_t, ut1, tt;

    mjd = jd - 2400000.5;
    geteop (mjd + utc/86400.0, &xp, &yp, &ut1_utc, &dx, &dy);
    lps = getlps (jd + utc/86400.0);
    delta_t = 32.184 + lps - ut1_utc;
    ut1 = utc + ut1_utc;
    tt = utc + (lps + 32.184);

    cel_pole (jd + tt / 86400.0, 2, dx * 1e3, dy * 1e3);

    cel2ter (jd, ut1 / 86400.0, delta_t, 1, 1, 0,
        xp, yp, vc, vt); /*--vc unit: m--*/

}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
 *
 * getlps - get the leap seconds value for input JD
 * 
 * jdutc: double, Julian Day of UTC 
 * return: short int, leap seconds
 *
 * version: Mar 2013
 * 
 */
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int getlps (double jdutc)
{
/*
 *
1972 JAN  1 =JD 2441317.5  TAI-UTC=  10.0       S + (MJD - 41317.) X 0.0      S
1972 JUL  1 =JD 2441499.5  TAI-UTC=  11.0       S + (MJD - 41317.) X 0.0      S
1973 JAN  1 =JD 2441683.5  TAI-UTC=  12.0       S + (MJD - 41317.) X 0.0      S
1974 JAN  1 =JD 2442048.5  TAI-UTC=  13.0       S + (MJD - 41317.) X 0.0      S
1975 JAN  1 =JD 2442413.5  TAI-UTC=  14.0       S + (MJD - 41317.) X 0.0      S
1976 JAN  1 =JD 2442778.5  TAI-UTC=  15.0       S + (MJD - 41317.) X 0.0      S
1977 JAN  1 =JD 2443144.5  TAI-UTC=  16.0       S + (MJD - 41317.) X 0.0      S
1978 JAN  1 =JD 2443509.5  TAI-UTC=  17.0       S + (MJD - 41317.) X 0.0      S
1979 JAN  1 =JD 2443874.5  TAI-UTC=  18.0       S + (MJD - 41317.) X 0.0      S
1980 JAN  1 =JD 2444239.5  TAI-UTC=  19.0       S + (MJD - 41317.) X 0.0      S
1981 JUL  1 =JD 2444786.5  TAI-UTC=  20.0       S + (MJD - 41317.) X 0.0      S
1982 JUL  1 =JD 2445151.5  TAI-UTC=  21.0       S + (MJD - 41317.) X 0.0      S
1983 JUL  1 =JD 2445516.5  TAI-UTC=  22.0       S + (MJD - 41317.) X 0.0      S
1985 JUL  1 =JD 2446247.5  TAI-UTC=  23.0       S + (MJD - 41317.) X 0.0      S
1988 JAN  1 =JD 2447161.5  TAI-UTC=  24.0       S + (MJD - 41317.) X 0.0      S
1990 JAN  1 =JD 2447892.5  TAI-UTC=  25.0       S + (MJD - 41317.) X 0.0      S
1991 JAN  1 =JD 2448257.5  TAI-UTC=  26.0       S + (MJD - 41317.) X 0.0      S
1992 JUL  1 =JD 2448804.5  TAI-UTC=  27.0       S + (MJD - 41317.) X 0.0      S
1993 JUL  1 =JD 2449169.5  TAI-UTC=  28.0       S + (MJD - 41317.) X 0.0      S
1994 JUL  1 =JD 2449534.5  TAI-UTC=  29.0       S + (MJD - 41317.) X 0.0      S
1996 JAN  1 =JD 2450083.5  TAI-UTC=  30.0       S + (MJD - 41317.) X 0.0      S
1997 JUL  1 =JD 2450630.5  TAI-UTC=  31.0       S + (MJD - 41317.) X 0.0      S
1999 JAN  1 =JD 2451179.5  TAI-UTC=  32.0       S + (MJD - 41317.) X 0.0      S
2006 JAN  1 =JD 2453736.5  TAI-UTC=  33.0       S + (MJD - 41317.) X 0.0      S
2009 JAN  1 =JD 2454832.5  TAI-UTC=  34.0       S + (MJD - 41317.) X 0.0      S
2012 JUL  1 =JD 2456109.5  TAI-UTC=  35.0       S + (MJD - 41317.) X 0.0      S
 *
*/

    short int lps;
    double jd = jdutc;

    if (jd > 2456109.5)
        lps = 35;
    else if (jd > 2454832.5)
        lps = 34;
    else if (jd > 2453736.5)
        lps = 33;
    else if (jd > 2451179.5)
        lps = 32;
    else if (jd > 2450630.5)
        lps = 31;
    else if (jd > 2450083.5)
        lps = 30;
    else if (jd > 2449534.5)
        lps = 29;
    else if (jd > 2449169.5)
        lps = 28;
    else if (jd > 2448804.5)
        lps = 27;
    else if (jd > 2448257.5)
        lps = 26;
    else if (jd > 2447892.5)
        lps = 25;
    else if (jd > 2447161.5)
        lps = 24;
    else if (jd > 2446247.5)
        lps = 23;
    else if (jd > 2445516.5)
        lps = 22;
    else if (jd > 2445151.5)
        lps = 21;
    else if (jd > 2444786.5)
        lps = 20;
    else 
    {
        printf ("No leapsecond configured before 1981 JUL  1 =JD 2444786.5\n");
        exit (0);
    }
    return lps;
  

}

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


void openeop (char file_eop[200], int mjds, int num, double *eopmat)
{
    FILE *fp_eop;
    int i, mjdi;
    char string[160];
   
    if ((fp_eop = fopen (file_eop,"r")) == NULL)
    {
        printf ("Cannot open eop file?\n");
        exit (0);
    }
        
//    for (i = 0; i < 13;i++)
//        fgets (string, 160, fp_eop);
    while (feof(fp_eop) == 0)
    {
        fgets (string, 160, fp_eop);
        sscanf (string, "%*d%*d%*d%d", &mjdi);
        if (mjdi == mjds - 1)
        {
            for (i = 0; i < num; i ++)
            {
                fgets (string, 160, fp_eop);
                sscanf (string, "%*d%*d%*d%lf%lf%lf%lf%*lf%lf%lf",
                    &eopmat[i * 6 + 0], &eopmat[i * 6 + 1], &eopmat[i * 6 + 2], 
                    &eopmat[i * 6 + 3], &eopmat[i * 6 + 4], &eopmat[i * 6 + 5]);
//                printf("mjd = %f\n", eopmat[i * 6 + 0]);
            }
            break;
        }
    }    

    fclose (fp_eop);
}





void aei2xyz (double ele[6], double pos[3], double vel[3])
{
    double a, e, i, omega, w, M, E,r, P[3], Q[3], n, GM, radius;
    int x;

    GM = GMCT;
    radius = RCT;
    a = ele[0]; 
    e = ele[1];
    i = ele[2] * DEG2RAD;
    omega = ele[3] * DEG2RAD;
    w = ele[4] * DEG2RAD;
    M = ele[5] * DEG2RAD;

    n=sqrt(GM/(a*a*a));

    E=kepler(M,e);
    P[0]=cos(omega)*cos(w)-sin(omega)*sin(w)*cos(i);
    P[1]=sin(omega)*cos(w)+cos(omega)*sin(w)*cos(i);
    P[2]=sin(w)*sin(i);
    Q[0]=-cos(omega)*sin(w)-sin(omega)*cos(w)*cos(i);
    Q[1]=-sin(omega)*sin(w)+cos(omega)*cos(w)*cos(i);
    Q[2]=cos(w)*sin(i);
    for(x=0;x<3;x++)
    {
        pos[x]=a*(cos(E)-e)*P[x]+a*sqrt(1-e*e)*sin(E)*Q[x];
    }
    r = modvect (pos);

    if (r <= radius)
    {
        printf("error: r <= radius ! in aei2xyz \n");
    }

    for(x=0;x<3;x++)
    {
        vel[x]=-a*a*n/r*sin(E)*P[x]+a*a*n/r*sqrt(1-e*e)*cos(E)*Q[x];
    }
}




double kepler(double M,double e)
{
    double E0,E1=M;
    do
    {
        E0=E1;
        E1=M+e*sin(E0);
    }
    while(fabs(E0-E1)>=1e-10);
    return(E1);
}


double xyz2aei(double ele[6], double pos[3], double vel[3])
{
    double a, e, omega, i, w, E, M, r, v, h, HV[3], n, 
        GM, radius, Pz, Qz;
    
    GM = GMCT;
    radius = RCT;
    r = modvect (pos);
    v = modvect (vel);
    
    if (r <= radius)
    {
        printf("error: r <= radius ! in xyz2aei \n");
    }

    a = 1.0 / (2.0 / r - v * v / GM);

    crsvect (pos, vel, HV);
    h = modvect(HV);
    e = sqrt (1.0 - h * h / GM / a);
    i = acos (HV[2] / h);     //unit: rad
    omega = chosephase (HV[0] / h / sin(i), - HV[1] / h / sin(i));   //unit: rad

    if(a <= 0)
    {    
        ele[0]=a,ele[1]=e,ele[2]=i * RAD2DEG;
        ele[3]=omega * RAD2DEG,ele[4]=0,ele[5]=0;
//        printf("error: a <= 0 !\n");
        return 0;
    }
    
    if(a <= radius)
    {    
        printf("warning: a <= radius !\n");
    }
    


    n = sqrt ( GM / (a*a*a) );

    E = chosephase ( dotvect(pos, vel) / (a * a * n * e), (1.0 - r / a) / e);  //unit: rad

    M = E - e * sin(E);        //unit: rad

    Pz = (cos(E) / r * pos[2] - sin(E) / n / a * vel[2]);
    Qz = (sin(E) / r / sqrt(1.0-e*e) * pos[2] + (cos(E) - e) / n / a / sqrt(1.0-e*e) * vel[2]);

    w = chosephase ( Pz / sin(i), Qz /sin(i));  //unit: rad

    ele[0] = a;
    ele[1] = e;
    ele[2] = i * RAD2DEG;
    ele[3] = omega * RAD2DEG;
    ele[4] = w * RAD2DEG;
    ele[5] = M * RAD2DEG;
                          
    return TWOPI / n;                              
}
                




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* opengravfile ¨C 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double opengrv (char file_grv[2][200], double *coef, int nmax, int mmax)
{
    FILE *fp_grv;
    double c,s;
    int n,m, l, ind;
    char string[200], name[20];

    if ((fp_grv = fopen (file_grv[0],"r")) == NULL)
    {
        printf ("Cannot open gravity file?\n");
        exit (0);
    }


//    coef[0] = 1;    // include zero degree term
    coef[0] = 0;    // exclude zero degree term

    while (1)
    {
        if (fgets (string, 200, fp_grv) == NULL) break;
//        sscanf (string, "%d%d%lf%lf", &n, &m, &c, &s);    

        if (strlen(file_grv[1])==0)
        {
            sscanf (string, "%d%d%lf%lf", &n, &m, &c, &s);  
        }
        else 
        {
            sscanf (string, "%s", name);    
            if (strcmp (name,file_grv[1]) ==0)  
            {
                sscanf (string, "%*s%d%d%lf%lf", &n, &m, &c, &s);
//                printf ("n = %d m = %d c = %e s = %e\n", n, m, c, s);
            }
        }


//        if (n > nmax || n < 0)
        if (n > nmax || n < 2 || m > mmax)   // permanently exclude degree 1 @7/24/2012
            continue;
        else if (m == 0)
        {
            coef[n] = c;
        }
        else 
        {
            l = nmax - m + 1;
            ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
            coef[ind + n - m] = c;
            coef[ind + n - m + l] = s;
//            coef[ind + n - m] = 0;
//            coef[ind + n - m + l] = 0;
        }
    }

    printf ("coef[2] = %e\n", coef[2]);
    if (PERMT == 1)
    {
        coef[2]  = coef[2] - 4.201e-9; //tn32
//        coef[2]  = coef[2] - 4.1736e-9;  //tn36
    }

    printf ("coef[2] = %e\n", coef[2]);
    fclose(fp_grv);
    return 0;
}












/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* xyz2llh - xyz to latitude, longitude, height
* @param1: description of param1
* @param2: description of param2

* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void xyz2llh (double *vt, double *llh)
{
    double r, cosphi, phi, costhe, sinthe;
    r = sqrt (vt[0] * vt[0] + vt[1] * vt[1] + vt[2] * vt[2]);

    cosphi = vt[2] / r;
    phi = acos(cosphi) ;
    costhe = vt[0] / r / sin(phi);
    sinthe = vt[1] / r / sin(phi);
    llh[2] = r;
    llh[1] = chosephase(sinthe, costhe) * RAD2DEG;
    llh[0] = 90.0 - phi * RAD2DEG;
}





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* fun_pointmass - abandoned
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double earth_pointmass (double jd, double tdbs, double *x, double *f)
{
    double GM, radius, gmde, fnt[3], fgr[3], r, s2, rrd, a, b;
    int n, gamma;

    GM = GMCT;
    radius = RCT;

    gmde = GM * 86400.0 * 86400.0 / AU / AU / AU;
    gamma = 1;

    f[0] = x[3]; 
    f[1] = x[4]; 
    f[2] = x[5];

    r = sqrt (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    s2 = x[3] * x[3] + x[4] * x[4] + x[5] * x[5];
    rrd = x[0] * x[3] + x[1] * x[4] + x[2] * x[5];
	
    a = 2 * (1 + gamma) * gmde / r - gamma * s2;
    b = 2 * (1 + gamma) * rrd;
    for (n = 0; n < 3; n++)
        fgr[n] =  gmde / C_AUDAY / C_AUDAY / r / r / r 
        * ( a * x[n] + b * x[n+3] );

    fnt[0] = - gmde / (r*r*r) * x[0];
    fnt[1] = - gmde / (r*r*r) * x[1];
    fnt[2] = - gmde / (r*r*r) * x[2];

    for (n = 0; n < 3; n++)
    {
        f[3 + n] = fnt[n]; 
    }

	return 0;
}



double accel_point (double *tjd, double *x, double *fnt, double *fgr)
{
    double GM, radius, gmde, r, s2, rrd, a, b;
    int n, gamma;

    GM = GMCT;
    radius = RCT;

    gmde = GM * 86400.0 * 86400.0 / AU / AU / AU;
    gamma = 1;

    r = sqrt (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    s2 = x[3] * x[3] + x[4] * x[4] + x[5] * x[5];
    rrd = x[0] * x[3] + x[1] * x[4] + x[2] * x[5];
	
    a = 2 * (1 + gamma) * gmde / r - gamma * s2;
    b = 2 * (1 + gamma) * rrd;
    for (n = 0; n < 3; n++)
        fgr[n] =  gmde / C_AUDAY / C_AUDAY / r / r / r 
        * ( a * x[n] + b * x[n+3] );

    fnt[0] = - gmde / (r*r*r) * x[0];
    fnt[1] = - gmde / (r*r*r) * x[1];
    fnt[2] = - gmde / (r*r*r) * x[2];

	return 0;
}












double accel_pmiers (double *tjd, double *x, double *fnt, double *fgr)
{
    double GME, GMS, J, Jv[3], beta, gamma, r, v2, pv, pJ, a, b, p[3], v[3],
        pxv[3], vxJ[3], ps[3], vs[3], rs, vsxps[3], vsxpsxv[3],
        term1[3], term2[3], term3[3];
    int n;
    short int sun = 10;

    GME = GMCT; //m^3/s^2
    J = 9.8e8; //m^2/s
    gamma = 1;
    beta = 1;
    GMS = 1.32712442076e20; //m^3/s^2

    GME = GME * 86400.0 * 86400.0 / AU / AU / AU;
    GMS = GMS * 86400.0 * 86400.0 / AU / AU / AU;
    J = J * 86400.0 / AU / AU; 
    
    Jv[0] = 0; Jv[1] = 0; Jv[2] = J;
    p[0] = x[0]; p[1] = x[1]; p[2] = x[2];
    v[0] = x[3]; v[1] = x[4]; v[2] = x[5];

    r = modvect(p);
    v2 = dotvect(v, v);
    pv = dotvect(p, v);
    pJ = dotvect(p, Jv);
    crsvect(p, v, pxv);
    crsvect(v, Jv, vxJ);    

    planet_ephemeris (tjd, CT, sun, ps, vs);
    rs = modvect(ps);
    crsvect(vs, ps, vsxps);
    crsvect(vsxps, v, vsxpsxv);    

    a = 2 * (beta + gamma) * GME / r - gamma * v2;
    b = 2 * (1 + gamma) * pv;
    for (n = 0; n < 3; n++)
    {
        term1[n] = GME / C_AUDAY / C_AUDAY / r / r / r *
                ( a * p[n] + b * v[n] );
        term2[n] = GME / C_AUDAY / C_AUDAY / r / r / r * (1 + gamma) * 
                ( 3/r/r * pxv[n] * pJ + vxJ[n] );
        term3[n] = - GMS / C_AUDAY / C_AUDAY / rs / rs / rs * (1 + 2 * gamma) *
                vsxpsxv[n];

        fgr[n] = term1[n]
                + term2[n] + term3[n];
//        printf ("%15.12f\t%15.12f\t%15.12f\n", term1[n],term2[n],term2[n]);
    }


    fnt[0] = - GME / (r*r*r) * p[0];
    fnt[1] = - GME / (r*r*r) * p[1];
    fnt[2] = - GME / (r*r*r) * p[2];

    return 0;
}





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* fun_fullstate -transition matrix(36), orbit(6), sensitivity matrix(6*DYNPAR)
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void fun_accel (int dim, double jd, double tt, double *state, double *fstate)
{
    int n, i,k, part;
    double tjd[2], xic[6], dfdx[36], dxdx0[36],
        acc1[3], dadr1[9],
        acc2[3], dadr2[9],
        acc3[3], dadr3[9], dadsrpb[3], dadsrpt[3],
        acc4[4], dadr4[9], dadk2[3], 
        acc[3], dadr[9], 
        fxic[6], fdxdx0[36];

    double ap[3], an[3], ar[3], ag[3], apgr[3], angr[3], at[3], ao[3];
    double *dadp, *dxdp, *dfdpp, *dfdp, *fdxdp;

    tjd[0] = jd;    tjd[1] = tt / 86400.0;
//    tjd[0] = jd;    tjd[1] = tt;

    if (dim < 6)
    {
        printf ("error: fun_accel input dim < 6!\n");
        exit (0);
    }
    else if (dim == 6)
        part = 0;
    else if (dim > 6)
    {
        part = 1;    
    }
    for (n = 0; n < 6; n++)
    {
        xic[n] = state[n];
    }   

/* acc, partial to xyz: dadr, partial to parameters dadp*/
//    accel_ntrel (tjd, xic, part, acc1, dadr1, dadp1);
//    accel_nonsp (tjd, xic, part, acc2, dadr2, dadp2);
//    accel_radpr (tjd, xic, part, acc3, dadr3, dadp3);
/*todo: air drag acc & partial to vxvyvz dadv*/

    accel_pm_part (tjd, xic, ap, part, dadr1);
    accel_nb_part (tjd, xic, an, part, dadr2);
    accel_sr_part (tjd, xic, ar, part, dadr3, dadsrpb, dadsrpt);
    accel_gt_part (tjd, xic, ag, part, dadr4, dadk2);



    for (n = 0; n <= 2; n++)
    {
        acc[n] = ap[n] + an[n] + ar[n] + ag[n];
    }

    fxic[0] = xic[3];
    fxic[1] = xic[4];
    fxic[2] = xic[5];
    fxic[3] = acc[0];
    fxic[4] = acc[1];
    fxic[5] = acc[2];

    for (n = 0; n < 6; n++)
    {
        fstate[n] = fxic[n];
    }

    if (part == 0)
    {
        return;
    }

    
    for (n = 0; n < 36; n++)
    {
        dxdx0[n] = state[n + 6];
    }

    for (n = 0; n <= 8; n++)
    {
        dadr[n] = dadr1[n] + dadr2[n] + dadr3[n] + dadr4[n];
//        dadr[n] = dadr1[n];
    }


    for (n = 0; n < 36; n++)
    {
        dfdx[n] = 0;
    }
    dfdx[3]  = 1; 
    dfdx[10] = 1; 
    dfdx[17] = 1;
    for (n = 0; n < 3; n++)
    {
        dfdx[n + 18] = dadr[n];
        dfdx[n + 24] = dadr[n + 3];
        dfdx[n + 30] = dadr[n + 6];
    }
    brmul(dfdx, dxdx0, 6, 6, 6, fdxdx0);


    for (n = 0; n < 36; n++)
    {
        fstate[n + 6] = fdxdx0[n];
    }


    if (MDYN == 0)
        return;


    dadp = (double *) calloc ( 3 * MDYN, sizeof(double));
    dxdp = (double *) calloc ( 6 * MDYN, sizeof(double));
    dfdpp = (double *) calloc ( 6 * MDYN, sizeof(double));
    dfdp = (double *) calloc ( 6 * MDYN, sizeof(double));
    fdxdp = (double *) calloc ( 6 * MDYN, sizeof(double));


    for (n = 0; n < 6 * MDYN; n++)
    {
        dxdp[n] = state[n + 42];
    }



    i = 0;
    if (MSRP > 0)
    {     
        for (n = 0; n < 3; n++)
        {
            dadp[n * MDYN + i] = dadsrpb[n];
        }
        i++;
    }
    if (MSRP > 1)
    {     
        for (n = 0; n < 3; n++)
        {
            dadp[n * MDYN + i] = dadsrpt[n];
        }
        i++;
    }
    if (MTK2 > 0)
    {     
        for (n = 0; n < 3; n++)
        {
            dadp[n * MDYN + i] = dadk2[n];
        }
        i++;
    }
    if (MGCS > 0)
    {     
        for (k = 0; k < MGCS; k ++)
        {
            for (n = 0; n < 3; n++)
            {
                dadp[n * MDYN + i] = CSinfo[k].dadcs[n];
            }
            i++;
        }
    }





    brmul(dfdx, dxdp, 6, 6, MDYN, dfdpp);
    for (n = 0; n < 3 * MDYN; n++)
    {
        dfdp[n] = 0;
        dfdp[n + 3 * MDYN] = dadp[n];
    }
    for (n = 0; n < 6 * MDYN; n++)
    {
        fdxdp[n] = dfdpp[n] + dfdp[n];
    }
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/


    for (n = 0; n < 6 * MDYN; n++)
    {
        fstate[n + 42]= fdxdp[n];
    }

    free (dadp);
    free (dxdp);
    free (dfdpp);
    free (dfdp);
    free (fdxdp);

    return;
}










double accel_pm_part (double *tjd, double *x, double *acc, int part, double *dadr)
{
    double GME, GMS, J, Jv[3], beta, gamma, r, v2, pv, pJ, a, b, p[3], v[3],
        pxv[3], vxJ[3], ps[3], vs[3], xsc[6], rs, vsxps[3], vsxpsxv[3],
        unit[9], ppt[9], r5, r3, fgr[3], fnt[3], term1[3], term2[3], term3[3];
    int n;
    short int sun = 10;

    GME = GMCT; //m^3/s^2

//    GME = GME * 86400.0 * 86400.0 / AU / AU / AU;
    p[0] = x[0]; p[1] = x[1]; p[2] = x[2];
    v[0] = x[3]; v[1] = x[4]; v[2] = x[5];

    r = modvect(p);

    acc[0] = - GME / (r*r*r) * p[0];
    acc[1] = - GME / (r*r*r) * p[1];
    acc[2] = - GME / (r*r*r) * p[2];

    
    if (part == 1)
    {
        unit[0] = 1; unit[1] = 0; unit[2] = 0;
        unit[3] = 0; unit[4] = 1; unit[5] = 0;
        unit[6] = 0; unit[7] = 0; unit[8] = 1;
        r5 = pow (r, 5);
        r3 = pow (r, 3);
        brmul (p, p, 3,1,3, ppt);
        for (n = 0; n <= 8; n++)
        {
            dadr[n] = 3 * GME * ppt[n] / r5
            - GME * unit[n] / r3;
        }
    }



    if (RELTIV == 0)
        return 0;

    J = 9.8e8; //m^2/s
    gamma = 1;
    beta = 1;
    GMS = 1.32712442076e20; //m^3/s^2
//    GMS = GMS * 86400.0 * 86400.0 / AU / AU / AU;
//    J = J * 86400.0 / AU / AU; 
    
    Jv[0] = 0; Jv[1] = 0; Jv[2] = J;

    v2 = dotvect(v, v);
    pv = dotvect(p, v);
    pJ = dotvect(p, Jv);
    crsvect(p, v, pxv);
    crsvect(v, Jv, vxJ);    

//    planet_ephemeris (tjd, CT, sun, ps, vs);
    get_ephemeris (tjd, CT, sun, xsc);
    for (n = 0; n < 3; n++)
    {
        ps[n] = xsc[n];
        vs[n] = xsc[n + 3];
    }

    rs = modvect(ps);
    crsvect(vs, ps, vsxps);
    crsvect(vsxps, v, vsxpsxv);    

    a = 2 * (beta + gamma) * GME / r - gamma * v2;
    b = 2 * (1 + gamma) * pv;
    for (n = 0; n < 3; n++)
    {
        term1[n] = GME / C / C / r / r / r *
                ( a * p[n] + b * v[n] );
        term2[n] = GME / C / C / r / r / r * (1 + gamma) * 
                ( 3/r/r * pxv[n] * pJ + vxJ[n] );
        term3[n] = - GMS / C / C / rs / rs / rs * (1 + 2 * gamma) *
                vsxpsxv[n];

        fgr[n] = term1[n]
                + term2[n] + term3[n];
//        printf ("%15.12f\t%15.12f\t%15.12f\n", term1[n],term2[n],term2[n]);
    }


    acc[0] = acc[0] + fgr[0];
    acc[1] = acc[1] + fgr[1];
    acc[2] = acc[2] + fgr[2];

    return 0;




}






int get_ephemeris (double tjd[2], int to, int from, double *x)
{
    double jd0 = 2451545.00000000, lt, tdbj2000, fromTtoS[6], pos[3], vel[3];
    int n;
    short int center, target;


    if (from <= 12 && to <= 12)
    {
        center = (short int)from;
        target = (short int)to;
        planet_ephemeris (tjd, target, center, pos, vel);
        x[0] = pos[0] * AU;
        x[1] = pos[1] * AU;
        x[2] = pos[2] * AU;
        x[3] = vel[0] * AU / 86400.0;
        x[4] = vel[1] * AU / 86400.0;
        x[5] = vel[2] * AU / 86400.0;
    }
    else if (from == I_TITAN || to == I_TITAN) //titan
    {

        tdbj2000 = ((tjd[0] - jd0) + tjd[1]) * 86400.0;

        spkezr_c ("SATURN", tdbj2000, "J2000", "NONE", "TITAN", fromTtoS, &lt);
/*
Procedure

   void spkezr_c ( ConstSpiceChar     *targ,
                   SpiceDouble         et,
                   ConstSpiceChar     *ref,
                   ConstSpiceChar     *abcorr,
                   ConstSpiceChar     *obs,
                   SpiceDouble         starg[6],
                   SpiceDouble        *lt        )
 
 
   Return the state (position and velocity) of a target body 
   relative to an observing body, optionally corrected for light 
   time (planetary aberration) and stellar aberration. 
*/

        if (from == I_TITAN) center = (short int) to;
        else center = (short int) from;
        planet_ephemeris (tjd, center, 5, pos, vel);

        for (n = 0; n < 3; n++)
        {    
            x[n] = pos[n] * AU + fromTtoS[n] * 1000.0;
            x[n + 3] = vel[n] * AU / 86400.0 + fromTtoS[n + 3] * 1000.0;
        }
        if (to == I_TITAN)
        {
            for (n = 0; n < 6; n++)
            {
                x[n] = - x[n];
            }
        }    
    }
    else 
    {
        printf ("error in get_ephemeris: from = %d\t to = %d\n", from, to);
    }
    return 0;
}









double accel_nb_part (double *tjd, double *xic, double *acc, int part, double *dadr)
{
    int n;
    short int ssbary = 11;
    double xcb[6], xib[6], acb[3], aib[3], dadr1[9], dadrc[9], dadri[9];

    if (NBODY == 0)
    {
        if (part == 1)
        {
            for (n = 0; n <= 8; n++)
            {
                dadr[n] = 0;
            }
        }

        for (n = 0; n <= 2; n++)
        {
            acc[n] = 0;
        }

        return 0;
    }


//    planet_ephemeris (tjd, CT, ssbary, &xcb[0], &xcb[3]);
    get_ephemeris (tjd, CT, ssbary, xcb);
    f_bcrs (tjd, xcb, CT, acb, part, dadrc);
    for (n = 0; n <= 5; n++)
    {
        xib[n] = xic[n] + xcb[n];
    }
    f_bcrs (tjd, xib, CT, aib, part, dadri);

    for (n = 0; n <= 2; n++)
    {
        acc[n] = aib[n] - acb[n];
    }

    if (part == 1)
    {
        for (n = 0; n <= 8; n++)
        {
            dadr[n] = dadri[n];
        }
    }

    return 0;
}











/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double f_bcrs (double *tjd, double *xi, int exclude, 
                   double *acc, int part, double *dadr)
{
    double fnt[3], fgr[3], xj[11][6], xij[11][6], rij[11], xjk[6], rjk, 
        xddj[3], sumil, sumjk, sdi2, sdj2, rdirdj, rrrdr2, rjirdd, 
        rij5, rij3, xijt[9], gra, grb, beta, gamma, unit[9];
    short int ssbary, l, k, j, n, flag_gr;
    
    ssbary = 11;
    gamma = 1.0;
    beta = 1.0;
    unit[0] = 1; unit[1] = 0; unit[2] = 0;
    unit[3] = 0; unit[4] = 1; unit[5] = 0;
    unit[6] = 0; unit[7] = 0; unit[8] = 1;

    for (j = 0; j <= 10; j++)
    {
//        planet_ephemeris (jd, j, ssbary, &xj[j][0], &xj[j][3]);
        get_ephemeris (tjd, j, ssbary, xj[j]);
        for (n = 0; n < 6; n++)
        {
            xij[j][n] = xi[n] - xj[j][n];
        }
        rij[j] = sqrt (xij[j][0] * xij[j][0] 
            + xij[j][1] * xij[j][1] + xij[j][2] * xij[j][2]);
    }
    
    flag_gr = 0;
    for (n = 0; n < 3; n ++)
        fnt[n] = 0;
    for (j = 0; j <= 10; j++)
    {
        if (PERB[j] == 2)
            flag_gr = 1;
        if (PERB[j] == 0)
            continue;
        if (j == exclude)
            continue;
        for (n = 0; n < 3; n++)
            fnt[n] = fnt[n] 
            - GMDE[j] / (rij[j] * rij[j] * rij[j]) * xij[j][n];
    }

    if (part == 1)
    {
        for (n = 0; n <= 8; n++)
        {
            dadr[n] = 0;
        }
        for (j = 0; j <= 10; j++)
        {
            if (j == exclude)
                continue;
            rij5 = pow (rij[j], 5);
            rij3 = pow (rij[j], 3);
            brmul (xij[j], xij[j], 3,1,3, xijt);
            for (n = 0; n <= 8; n++)
            {
                dadr[n] = dadr[n] + 3 * GMDE[j] * xijt[n] / rij5
                - GMDE[j] * unit[n] / rij3;
            }
        }
    }

    if (flag_gr == 0)
    {
        for (n = 0; n < 3; n++)
            acc[n] =  fnt[n];
        return 0;
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    sdi2 = xi[3] * xi[3] + xi[4] * xi[4] + xi[5] * xi[5];
    sumil = 0;
    for (l = 0; l < 11; l ++)
    {
        if ( l == exclude)
            continue;
        if (PERB[l] != 2)
            continue;
        sumil = sumil + GMDE[l] / rij[l];
    }

    for (n = 0; n < 3; n ++)
        fgr[n] = 0;
    for (j = 0; j < 11; j ++)
    {
        if (PERB[j] != 2)
            continue;
        if (j == exclude)
            continue;
        sumjk = 0;
        for (n = 0; n < 3; n ++)
            xddj[n] = 0;
        for (k = 0; k < 11; k ++)
        {
            if (k == j)	
                continue;	//k!=j
            if (PERB[k] != 2)
                continue;
            for (n = 0; n < 3; n++)
                xjk[n] = xj[j][n] - xj[k][n];
            rjk = sqrt (xjk[0] * xjk[0] + xjk[1] * xjk[1] + xjk[2] * xjk[2]);
            sumjk = sumjk + GMDE[k] / rjk;
            for (n = 0; n < 3; n ++)
                xddj[n] = xddj[n] - GMDE[k] / (rjk * rjk * rjk) * xjk[n];
        }
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
        sdj2 = xj[j][3] * xj[j][3] + xj[j][4] * xj[j][4] 
            + xj[j][5] * xj[j][5];
        rdirdj = xi[3] * xj[j][3] + xi[4] * xj[j][4] + xi[5] * xj[j][5];
        rrrdr2 = pow( ( xij[j][0] * xj[j][3] + xij[j][1] * xj[j][4] 
            + xij[j][2] * xj[j][5]) / rij[j], 2);
        rjirdd = - ( xij[j][0] * xddj[0] + xij[j][1] * xddj[1] 
            + xij[j][2] * xddj[2]);
        
        gra = - 2 * (beta + gamma) * sumil - (2 * beta -1) * sumjk 
            + gamma * sdi2 + (1 + gamma) * sdj2
            - 2 * (1 + gamma) * rdirdj - 1.5 * rrrdr2 + 0.5 * rjirdd;

        grb = xij[j][0] * ((2+2*gamma) * xi[3] - (1+2*gamma) * xj[j][3])
            + xij[j][1] * ((2+2*gamma) * xi[4] - (1+2*gamma) * xj[j][4])
            + xij[j][2] * ((2+2*gamma) * xi[5] - (1+2*gamma) * xj[j][5]);

        for (n = 0; n < 3; n ++)
        {
            fgr[n] = fgr[n] 
                + GMDE[j] / (rij[j] * rij[j] * rij[j]) 
                * ( - xij[j][n]) * gra / C / C
                + GMDE[j] / (rij[j] * rij[j] * rij[j]) 
                * xij[j][n + 3] * grb / C / C
                + GMDE[j] / rij[j] * (3 + 4 * gamma) * 0.5   
                * xddj[n] / C / C;
        }
    }

    for (n = 0; n < 3; n++)
        acc[n] =  fgr[n] + fnt[n];
    return 1;
}



















double accel_sr_part (double *tjd, double *xic, double *acc, int part,
        double *dadr, double *dadsrpb, double *dadsrpt)
{
    double j, c1, ap, m, rsp, usp[3], xis[6], xsc[6], f, 
        xist[9], unit[9], rsp3;
    short int n, sun;
    short int ssbary = 11;


    if (AMR == 0)
    {
        if (part == 1)
        {
            for (n = 0; n <= 8; n++)
            {
                dadr[n] = 0;
            }
            for (n = 0; n <= 3; n++)
            {
                dadsrpb[n] = 0;
                dadsrpt[n] = 0;
            }
        }

        for (n = 0; n <= 2; n++)
        {
            acc[n] = 0;
        }

        return 0;
    }


    sun = 10;
    unit[0] = 1; unit[1] = 0; unit[2] = 0;
    unit[3] = 0; unit[4] = 1; unit[5] = 0;
    unit[6] = 0; unit[7] = 0; unit[8] = 1;

//    planet_ephemeris (tjd, sun, CT, &xsc[0], &xsc[3]);    
    get_ephemeris (tjd, sun, CT, xsc);
    for (n = 0; n <= 5; n++)
    {
        xis[n] = xic[n] - xsc[n];
    }
    rsp = sqrt (xis[0] * xis[0] + xis[1] * xis[1] + xis[2] * xis[2]);
    usp[0] = xis[0] / rsp;
    usp[1] = xis[1] / rsp;
    usp[2] = xis[2] / rsp;

    j  = 1352.5;   //kg/s3
//    j  = 1359.4;   //kg/s3
//    m  = SATMASS;     //kg
//    ap = SATAREA;        //m2
    c1 = j / C * AU * AU;   //kg/s2/m*au*au
    f  = c1 * AMR  / rsp / rsp;    
//    f  = c1 * ap / m  / rsp / rsp;    
//kg/s2/m*au*au * m2 / kg / au  / au = m/s2
//    f = f / AU * 86400.0 * 86400.0;


//    printf ("SRPB = %f\t SRPT = %f\n", SRPB, SRPT);
    for (n = 0; n < 3; n++)
    {
       acc[n] = f * usp[n] * (1 + SRPB + SRPT * tjd[1]);
    }

    if (part == 0)
        return 1;

    rsp3 = rsp * rsp * rsp;
    brmul (xis, xis, 3,1,3, xist);
    for (n = 0; n <= 8; n++)
    {
        dadr[n] = - f * (3 * xist[n] / rsp3 - unit[n] / rsp) * 
            (1 + SRPB + SRPT * tjd[1]);
    }

    for (n = 0; n <= 2; n++)
    {
        dadsrpb[n] = f * usp[n];      
        dadsrpt[n] = f * usp[n] * tjd[1];   
    }

    return 0;
}




















double accel_gt_part (double *tjd, double *xic, double *acc, int part,
        double *dadr, double *dadk2)
{
    int n,k, lps, ntide, blst[12], nb;
    double GM, radius, *tmp, pi[3], pe[3], llr[3], c_ie[9], c_ei[9], ae[3], ai[3];
    double jd0, tt, utc, te[9], tx[3], ty[3], tz[3], ao[3], as[3], ag[3], dadk2e[3],
        vx[3] = {1,0,0}, vy[3] = {0,1,0}, vz[3] = {0,0,1}, dadre[9], dadres[9], dadrei[9];
    InfStruct info;

    GM = GMCT;
    radius = RCT;

    for (n = 0; n <= 2; n++)
    {
        acc[n] = 0;
    }

    if (NMAX < 2)
    {
        if (part == 1)
        {
            for (n = 0; n <= 8; n++)
            {
                dadr[n] = 0;
            }
            for (n = 0; n <= 3; n++)
            {
                dadk2[n] = 0;
            }
        }

        return 0;
    }


    for (n = 0; n < 3; n++)
    {
        pi[n] = xic[n];
    }

    getinfo(tjd, &info);

    brmul(info.c_ie, pi, 3, 3, 1, pe);  
    
    xyz2llh(pe, llr);
        
//    cs2acc (llr, COEFG, GM, radius, NMAX, ae);
    cs2ada (llr, COEFG, NMAX, ae, part, dadre, 1);

    brmul(info.c_ei, ae, 3, 3, 1, ag);

    for (n = 0; n < 3; n++)
    {
        acc[n] = ag[n];
    }

    if (STIDE != 0)
    {

        if (STIDE == 3 && CT == 2)
            stidecs_Anelastic(&info, 1, COEFS);
        else if (STIDE == 2 && CT == 2)
            stidecs(tjd, info.c_ie, 1, COEFS);
        else if (CT == 2)
        {
            blst[0] = 10; blst[1] = 9; nb = 2;
            stidecs_k2 (&info, K2, COEFS, blst, nb);
        }
        else if (CT == 9)
        {
            blst[0] = 10; blst[1] = 2; nb = 2;
            stidecs_k2 (&info, K2, COEFS, blst, nb);
        }
        else if (CT == 20)
        {
            blst[0] = 10; blst[1] = 5; nb = 2;
            stidecs_k2 (&info, K2, COEFS, blst, nb);
        }
        else 
        {
            blst[0] = 10; nb = 1;
            stidecs_k2 (&info, K2, COEFS, blst, nb);
        }
//        cs2acc (llr, COEFS, GM, radius, NSMAX, ae);

        cs2ada (llr, COEFS, NSMAX, ae, part, dadres, 0);

        brmul(info.c_ei, ae, 3, 3, 1, as);
    
        for (n = 0; n < 3; n++)
        {
            acc[n] = acc[n] + as[n];
        }   
        for (n = 0; n <= 8; n++)
        {
            dadre[n] = dadre[n] + dadres[n];
        }
    
    }


    if (OTIDE != 0) // N.A.
    {
        otidecs(info.jdt, info.gmst, NOMAX, COEFO);

        cs2acc (llr, COEFO, GM, radius, NOMAX, ae);

        brmul(info.c_ei, ae, 3, 3, 1, ao);

        for (n = 0; n < 3; n++)
        {
            acc[n] = acc[n] + ao[n];
        }   
    }

    if (part == 0)
        return 1;


    brmul(dadre, info.c_ie, 3, 3, 3, dadrei);			
    brmul(info.c_ei, dadrei, 3, 3, 3, dadr);	

    if (MTK2 == 1)
    {
        stidecs_k2 (&info, 1, COEFS, blst, nb);
        cs2ada (llr, COEFS, NSMAX, dadk2e, 0, dadres, 0);
        brmul(info.c_ei, dadk2e, 3, 3, 1, dadk2);
    }

    if (MGCS > 0)
    {
        for (k = 0; k < MGCS; k ++)
        {
            brmul(info.c_ei, CSinfo[k].dadcse, 3, 3, 1, CSinfo[k].dadcs);
        }    
    }
    return 0;

}





//    nmax = 4;      
//    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
double stidecs_k2(InfStruct *info, double k2, double *stcs, int *body, int nbody)
{
    int sun, i;

    double xs[6], gms2e, tjd[2], 
        pse[3], llrs[3], pbar[4], t,
        p20s, p30s, p21s, p31s, p22s, p32s, p33s,
        rers, c20, c21, s21, c22, s22;

    tjd[0] = info->jd0;
    tjd[1] = info->tt/86400.0;


    c20 = 0; c21 = 0; s21 = 0; c22 = 0; s22 = 0;
    for (i = 0; i < nbody; i++)
    {
        sun = body[i]; 
        if (sun > 12)
            continue;
    
        get_ephemeris (tjd, sun, CT, xs);
//        planet_ephemeris (tjd, sun, earth, ps, vs);

        brmul (info->c_ie, xs, 3, 3, 1, pse);  
    
        xyz2llh(pse, llrs);

        t = sin(llrs[0] * DEG2RAD);
        lgdr(t, 3, 0, pbar); p20s = pbar[2]; p30s = pbar[3];
        lgdr(t, 3, 1, pbar); p21s = pbar[1]; p31s = pbar[2];
        lgdr(t, 3, 2, pbar); p22s = pbar[0]; p32s = pbar[1];
        lgdr(t, 3, 3, pbar); p33s = pbar[0];


        gms2e  =  GMDE[sun]/GMCT;
        rers = RCT / llrs[2];

        c20 += k2/5.0 * ( gms2e * pow(rers, 3) * p20s );
        c21 += k2/5.0 * ( gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) );
        s21 += k2/5.0 * ( gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );
        c22 += k2/5.0 * ( gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) );
        s22 += k2/5.0 * ( gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );
    }
    
    stcs[0]  = 0; //c00; 
    stcs[1]  = 0; //c10; 
    stcs[2]  = c20; 
    stcs[3]  = 0; 
    stcs[4]  = c21; 

    stcs[5]  = 0; //s11; 
    stcs[6]  = s21; 
    stcs[7]  = c22; 
    stcs[8]  = s22; 
    
    return 0;

}








/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* fun_fullaccel - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double earth_fullaccel (double jd, double tt, double *xic, double *fxic)
{
    int n;
    short int ssbary = 11, part = 0;
    double tjd[2], ap[3], an[3], ar[3], ag[3], apgr[3], angr[3], at[3], ao[3],
        acc[3];

    tjd[0] = jd;
    tjd[1] = tt;
//    printf("%f\t%f\n", jd, tt);

//    accel_point (tjd, xic, a1, a1gr);
    accel_pmiers (tjd, xic, ap, apgr);


    if (NBODY == 1)
        accel_nbody (tjd, xic, an, angr);
    else 
    {
        for (n = 0; n < 3; n++)
        {
            an[n] = 0;
            angr[n] = 0;
        }
    }



    if (AMR != 0 )
        accel_slrad (tjd, xic, ar);
    else 
    {
        for (n = 0; n < 3; n++)
        {
            ar[n] = 0;
        }
    }



    if (NMAX >= 2 )
    {
//        if (STIDE == 1)
//        {
            accel_gravtide (tjd, xic, ag, at, ao);
//        }
//        else
//        {
//            accel_gravt (tjd, xic, ag);
//            at[0] = 0; at[1] = 0; at[2] = 0;
//        }
    }
    else 
    {
        for (n = 0; n < 3; n++)
        {
            ag[n] = 0;
            at[n] = 0;
            ao[n] = 0;
        }
    }



    if (RELTIV == 0)
    {
        for (n = 0; n < 3; n++)
        {
            apgr[n] = 0;
            angr[n] = 0;
        }

    }
    else if (RELTIV == 1)
    {
        for (n = 0; n < 3; n++)
        {
            angr[n] = 0;
        }
    }



    for (n = 0; n < 3; n++)
    {
//        a1[n] = 0; a1gr[n] = 0;
//        a2gr[n] = 0;
//        ao[n] = 0;
        acc[n] = ap[n] + an[n] + ar[n] + ag[n] + at[n] + ao[n] 
            + apgr[n] + angr[n];
//        printf ("ag=%e at=%e ao=%e \n", ag[n], at[n], ao[n]);
    }

//    exit(0);

    fxic[0] = xic[3];
    fxic[1] = xic[4];
    fxic[2] = xic[5];
    fxic[3] = acc[0];
    fxic[4] = acc[1];
    fxic[5] = acc[2];

    return 0;
}



double accel_gravtide (double *tjd, double *xic, double *ag, double *as, double *ao)
{
    int n, lps, ntide;
    double GM, radius, *stcs, pi[3], pe[3], llr[3], c_ie[9], c_ei[9], ae[3], ai[3];
    double jd0, tt, utc, te[9], tx[3], ty[3], tz[3],
        vx[3] = {1,0,0}, vy[3] = {0,1,0}, vz[3] = {0,0,1};
    InfStruct info;

    GM = GMCT;
    radius = RCT;

 
    for (n = 0; n < 3; n++)
    {
        pi[n] = xic[n] * AU;
    }


    getinfo(tjd, &info);
/*
        printf ("jd0 = %.10f\t lps = %d\t gmt = %.10f\n", info.jd0, info.leaps, info.gmst );
        for (n = 0; n < 9; n++)
            printf ("%e\n", info.c_ie[n]);

    


    jd0 = tjd[0];
    tt = tjd[1] * 86400.0;

    lps = getlps (jd0 + tt/86400.0);
    utc = tt - (lps + 32.184);

//    printf("%f\t%f\n", jd0, utc);
    itrf2gcrf(jd0, utc, vx, tx);
    itrf2gcrf(jd0, utc, vy, ty);
    itrf2gcrf(jd0, utc, vz, tz);


    for (n = 0; n < 3; n++)
    {
        c_ei[n*3] = tx[n];
        c_ei[n*3+1] = ty[n];
        c_ei[n*3+2] = tz[n];  //c_ei mean ITRF2ICRF
    }

    mt(c_ei, 3, 3, c_ie);

        for (n = 0; n < 9; n++)
            printf ("%e\n", c_ie[n]);

    exit(0);
*/

    brmul(info.c_ie, pi, 3, 3, 1, pe);  
    
    xyz2llh(pe, llr);
//    xyz2llh(pi, llr);

        
    cs2acc (llr, COEFG, GM, radius, NMAX, ae);
//    cs2acc (llr, COEFG, GM, radius, NMAX, ai);
//    printf("ae = %f\t%f\t%f\n", ae[0], ae[1], ae[2]);
//    exit(0);

    brmul(info.c_ei, ae, 3, 3, 1, ai);
    
    for (n = 0; n < 3; n++)
    {
        ag[n] = ai[n] / AU * 86400 * 86400;
    }


    if (STIDE == 0 && OTIDE == 0)
    {
        for (n = 0; n < 3; n++)
        {
            as[n] = 0;
            ao[n] = 0;
        }
        return 0;
    }



//    ntide = 4;

//    COEFS = (double *) calloc ( (ntide + 1) * (ntide + 1), sizeof(double));

//    id_perm = PERM;
//    stidecs(tjd, c_ie, 1, stcs);
    if (STIDE == 1)
        stidecs_Anelastic(&info, 1, COEFS);
    if (STIDE == 2)
        stidecs(tjd, info.c_ie, 1, COEFS);

    cs2acc (llr, COEFS, GM, radius, NSMAX, ae);

    brmul(info.c_ei, ae, 3, 3, 1, ai);

    for (n = 0; n < 3; n++)
    {
        as[n] = ai[n] / AU * 86400 * 86400;
    }


    if (OTIDE == 0)
    {
        for (n = 0; n < 3; n++)
        {
            ao[n] = 0;
        }
        return 0;
    }

//    NOMAX = 4;
//    COEFO = (double *) calloc ( (NOMAX + 1) * (NOMAX + 1), sizeof(double));

    otidecs(info.jdt, info.gmst, NOMAX, COEFO);
   
//    printf ("jdt = %e\t gmst = %f\t NOMAX = %d\t COEFO = %e\n", info.jdt, info.gmst, NOMAX, COEFO[5]);
//    for (n=0;n<25;n++)
//        printf ("COEFO = %e\t COEFS = %e\n", COEFO[n], COEFS[n]);

    cs2acc (llr, COEFO, GM, radius, NOMAX, ae);

    brmul(info.c_ei, ae, 3, 3, 1, ai);

    for (n = 0; n < 3; n++)
    {
        ao[n] = ai[n] / AU * 86400 * 86400;
    }

    if (STIDE == 0)
    {
        for (n = 0; n < 3; n++)
        {
            as[n] = 0;
        }
        return 0;
    }

    return 0;

}




int openotcs (char *infile)
{
    FILE *fp_ot;
    int i;
    char string[100];

    if ((fp_ot = fopen (infile,"r")) == NULL)
    {
        printf ("Cannot open otide file?\n");
        exit (0);
    }

//    fgets (string, 100, fp_ot);
//    fgets (string, 100, fp_ot);
//    fgets (string, 100, fp_ot);
//    fgets (string, 100, fp_ot);

    i = 0;
    while (1)
    {
        if (fgets (string, 100, fp_ot) == NULL) break;
        sscanf (string, "%lf%s%d%d%lf%lf%lf%lf", 
            &otfes[i].ds, &otfes[i].name, &otfes[i].n, &otfes[i].m, &otfes[i].cp, &otfes[i].sp, &otfes[i].cm, &otfes[i].sm);    

        otfes[i].argn[0] = (int)(otfes[i].ds/100)%10;
        otfes[i].argn[1] = (int)(otfes[i].ds/10)%10 - 5;
        otfes[i].argn[2] = (int)(otfes[i].ds/1)%10 - 5;
        otfes[i].argn[3] = (int)(otfes[i].ds*10)%10 - 5;
        otfes[i].argn[4] = (int)(otfes[i].ds*100)%10 - 5;
        otfes[i].argn[5] = (int)(otfes[i].ds*1000)%10 - 5;
        i++;
    }

//    (*n) = i;
    fclose(fp_ot);

    return 0;

}


double otidecs(double jdt, double gmst, int nmax, double *coef)
{
    double doodarg[6], ang, cp, sp, cm, sm;
    int i, ncon = 1, n,m, l, ind;

    for (i = 0; i < NFES; i++)
    {
        if (otfes[i].n > nmax)
        {
            continue;
        }

        n = otfes[i].n;
        m = otfes[i].m;
        cp = otfes[i].cp;
        sp = otfes[i].sp;
        cm = otfes[i].cm;
        sm = otfes[i].sm;

//        DOODSN(&jdt, &gmst, otfes[i].argn, &ncon, doodarg, &ang);
//  ang=0;
//        printf ("cp = %e\t sp = %e\t ang = %e\t doodarg = %e\n", cp, sp, ang, doodarg[3]);

        if (m == 0)
        {
            coef[n] = coef[n] + 1e-11 * ((cp+cm) * cos(ang) + (sp+sm)*sin(ang));
        }
        else 
        {
            l = nmax - m + 1;
            ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
            coef[ind + n - m] = coef[ind + n - m] + 1e-11 * ((cp+cm) * cos(ang) + (sp+sm)*sin(ang));
            coef[ind + n - m + l] = coef[ind + n - m + l] + 1e-11 * ((sp-sm) * cos(ang) - (cp-cm)*sin(ang));
        }
    }



    return 0;

}







//    nmax = 4;      
//    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
double stidecs_Anelastic(InfStruct *info, int id_perm, double *stcs)
{

//    double gms2e  =  332946.048166;                   
//    double gmm2e  =  1/81.3005690699;
//    double gms2e  =  332946.0487185;                   
//    double gmm2e  =  1/81.3005538970823;

    double GMsun = 1.32712442076e20;
    double gms2e, gmm2e  = 0.0123000383;

//    double c20pt = -4.1736e-9;
    double c20pt = -4.201e-9;

    double k20   =  0.29525;
    double k21   =  0.29470;
    double k22   =  0.29801;

    double REk20 =  0.30190;
    double REk21 =  0.29830;
    double REk22 =  0.30102;
    double IMk21 = -0.00144;
    double IMk22 = -0.00130;
    double k20pa = -0.00089;
    double k21pa = -0.00080;
    double k22pa = -0.00057;

    double k20p  = -0.00087;
    double k21p  = -0.00079;
    double k22p  = -0.00057;

    double k30   =  0.093;
    double k31   =  0.093;
    double k32   =  0.093;
    double k33   =  0.094;

    short int moon = 9, earth = 2, sun = 10, n;

    double ps[3], vs[3], pm[3], vm[3], tjd[2],
        pse[3], pme[3], llrs[3], llrm[3], pbar[4], t,
        p20m, p30m, p21m, p31m, p22m, p32m, p33m, 
        p20s, p30s, p21s, p31s, p22s, p32s, p33s,
        rerm, rers, c20, c30, c40, c21, s21, c22, s22, c31, s31, 
        c32, s32, c33, s33, c41, s41, c42, s42, 
        c20f, c21f, s21f, c22f, s22f;

    double GM, radius;

    GM = 398600.44180E+09;
    radius = 6378136.6;

    gms2e = GMsun/GM;

    tjd[0] = info->jd0;
    tjd[1] = info->tt/86400.0;

// Luni-solar ephemeris

    planet_ephemeris (tjd, sun, earth, ps, vs);
    planet_ephemeris (tjd, moon, earth, pm, vm);
    for (n = 0; n < 3; n++)
    {
        ps[n] = ps[n] * AU;
        pm[n] = pm[n] * AU;
    }
 
    
//    icrf2itrf(num, ps, pse);
//    icrf2itrf(num, pm, pme);

    brmul (info->c_ie, ps, 3, 3, 1, pse);   //inertial to fixed matrix gmat = rmat*tbt
    brmul (info->c_ie, pm, 3, 3, 1, pme);   //inertial to fixed matrix gmat = rmat*tbt


    xyz2llh(pse, llrs);
    xyz2llh(pme, llrm);

    t = sin(llrm[0] * DEG2RAD);
    lgdr(t, 3, 0, pbar); p20m = pbar[2]; p30m = pbar[3];
    lgdr(t, 3, 1, pbar); p21m = pbar[1]; p31m = pbar[2];
    lgdr(t, 3, 2, pbar); p22m = pbar[0]; p32m = pbar[1];
    lgdr(t, 3, 3, pbar); p33m = pbar[0];
    t = sin(llrs[0] * DEG2RAD);
    lgdr(t, 3, 0, pbar); p20s = pbar[2]; p30s = pbar[3];
    lgdr(t, 3, 1, pbar); p21s = pbar[1]; p31s = pbar[2];
    lgdr(t, 3, 2, pbar); p22s = pbar[0]; p32s = pbar[1];
    lgdr(t, 3, 3, pbar); p33s = pbar[0];


    rerm = radius / llrm[2];
    rers = radius / llrs[2];

// Frequency Independent Terms

// C20
    c20 = REk20/5.0 * ( gmm2e * pow(rerm, 3) * p20m 
                    + gms2e * pow(rers, 3) * p20s );
// C21/S21
    c21 = + REk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) )
          + IMk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD) 
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );


    s21 = - IMk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) )
          + REk21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD) 
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );


// C22/S22
    c22 = + REk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) )
          + IMk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );

    s22 = - IMk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) )
          + REk22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );


// C30
    c30 = k30/7.0 * ( gmm2e * pow(rerm, 4) * p30m 
                    + gms2e * pow(rers, 4) * p30s );
// C31/S31
    c31 = k31/7.0 * ( gmm2e * pow(rerm, 4) * p31m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 4) * p31s * cos(llrs[1] * DEG2RAD) );
    s31 = k31/7.0 * ( gmm2e * pow(rerm, 4) * p31m * sin(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 4) * p31s * sin(llrs[1] * DEG2RAD) );
// C32/S32
    c32 = k32/7.0 * ( gmm2e * pow(rerm, 4) * p32m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 4) * p32s * cos(llrs[1] * DEG2RAD * 2.0) );
    s32 = k32/7.0 * ( gmm2e * pow(rerm, 4) * p32m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 4) * p32s * sin(llrs[1] * DEG2RAD * 2.0) );
// C33/S33
    c33 = k33/7.0 * ( gmm2e * pow(rerm, 4) * p33m * cos(llrm[1] * DEG2RAD * 3.0)
                    + gms2e * pow(rers, 4) * p33s * cos(llrs[1] * DEG2RAD * 3.0) );
    s33 = k33/7.0 * ( gmm2e * pow(rerm, 4) * p33m * sin(llrm[1] * DEG2RAD * 3.0)
                    + gms2e * pow(rers, 4) * p33s * sin(llrs[1] * DEG2RAD * 3.0) );

// C40
    c40 = k20pa/5.0* ( gmm2e * pow(rerm, 3) * p20m 
                    + gms2e * pow(rers, 3) * p20s );
// C41/S41
    c41 = k21pa/5.0* ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) );
    s41 = k21pa/5.0* ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );
// C42/S42
    c42 = k22pa/5.0* ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) );
    s42 = k22pa/5.0* ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );

    
    
    
    stcs[0]  = 0; //c00; 
    stcs[1]  = 0; //c10; 
    stcs[2]  = c20; 
    stcs[3]  = c30; 
    stcs[4]  = c40; 

    stcs[5]  = 0; //c11; 
    stcs[6]  = c21; 
    stcs[7]  = c31; 
    stcs[8]  = c41; 
    stcs[9]  = 0; //s11; 
    stcs[10] = s21; 
    stcs[11] = s31; 
    stcs[12] = s41; 

    stcs[13] = c22; 
    stcs[14] = c32; 
    stcs[15] = c42; 
    stcs[16] = s22; 
    stcs[17] = s32; 
    stcs[18] = s42; 

    stcs[19] = c33; 
    stcs[20] = 0; //c43; 
    stcs[21] = s33; 
    stcs[22] = 0; //s43; 

    stcs[23] = 0; //c44; 
    stcs[24] = 0; //s44; 
  
    
// Frequency Dependent Terms

    c20f = 0; c21f = 0; s21f = 0; c22f = 0; s22f = 0;
//    stfrqdep(info->jdt, info->gmst, &c20f, &c21f, &s21f, &c22f, &s22f);


    stcs[2]  = c20 + c20f; 
    stcs[6]  = c21 + c21f; 
    stcs[10] = s21 + s21f; 
    stcs[13] = c22 + c22f; 
    stcs[16] = s22 + s22f; 

    if(id_perm==1) 
    {
        stcs[2]  = c20 + c20f - c20pt; 
    }
    return 0;

}





//    nmax = 4;      
//    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
double stidecs(double *tjd, double *c_ie, int id_perm, double *stcs)
{

//    double gms2e  =  332946.048166;                   
//    double gmm2e  =  1/81.3005690699;
//    double gms2e  =  332946.0487185;                   
//    double gmm2e  =  1/81.3005538970823;

    double GMsun = 1.32712442076e20;
    double gms2e, gmm2e  = 0.0123000383;

//    double c20pt = -4.1736e-9;
    double c20pt = -4.201e-9;

    double k20   =  0.29525;
    double k21   =  0.29470;
    double k22   =  0.29801;
/*
    double REk20 =  0.30190;
    double REk21 =  0.29830;
    double REk22 =  0.30102;
    double IMk21 = −0.00144;
    double IMk22 = −0.00130;
    double k20pa = −0.00089;
    double k21pa = −0.00080;
    double k22pa = −0.00057;
*/
    double k20p  = -0.00087;
    double k21p  = -0.00079;
    double k22p  = -0.00057;

    double k30   =  0.093;
    double k31   =  0.093;
    double k32   =  0.093;
    double k33   =  0.094;

    short int moon = 9, earth = 2, sun = 10, n;

    double ps[3], vs[3], pm[3], vm[3],
        pse[3], pme[3], llrs[3], llrm[3], pbar[4], t,
        p20m, p30m, p21m, p31m, p22m, p32m, p33m, 
        p20s, p30s, p21s, p31s, p22s, p32s, p33s,
        rerm, rers, c20, c30, c40, c21, s21, c22, s22, c31, s31, 
        c32, s32, c33, s33, c41, s41, c42, s42, 
        c20f, c21f, s21f, c22f, s22f;

    double GM, radius;

    GM = 398600.44180E+09;
    radius = 6378136.6;

    gms2e = GMsun/GM;

//    tjd[0] = info[num].jd0;
//    tjd[1] = info[num].tt/86400.0;

// Luni-solar ephemeris

    planet_ephemeris (tjd, sun, earth, ps, vs);
    planet_ephemeris (tjd, moon, earth, pm, vm);
    for (n = 0; n < 3; n++)
    {
        ps[n] = ps[n] * AU;
        pm[n] = pm[n] * AU;
    }
 
    
//    icrf2itrf(num, ps, pse);
//    icrf2itrf(num, pm, pme);

    brmul (c_ie, ps, 3, 3, 1, pse);   //inertial to fixed matrix gmat = rmat*tbt
    brmul (c_ie, pm, 3, 3, 1, pme);   //inertial to fixed matrix gmat = rmat*tbt


    xyz2llh(pse, llrs);
    xyz2llh(pme, llrm);

    t = sin(llrm[0] * DEG2RAD);
    lgdr(t, 3, 0, pbar); p20m = pbar[2]; p30m = pbar[3];
    lgdr(t, 3, 1, pbar); p21m = pbar[1]; p31m = pbar[2];
    lgdr(t, 3, 2, pbar); p22m = pbar[0]; p32m = pbar[1];
    lgdr(t, 3, 3, pbar); p33m = pbar[0];
    t = sin(llrs[0] * DEG2RAD);
    lgdr(t, 3, 0, pbar); p20s = pbar[2]; p30s = pbar[3];
    lgdr(t, 3, 1, pbar); p21s = pbar[1]; p31s = pbar[2];
    lgdr(t, 3, 2, pbar); p22s = pbar[0]; p32s = pbar[1];
    lgdr(t, 3, 3, pbar); p33s = pbar[0];


    rerm = radius / llrm[2];
    rers = radius / llrs[2];

// Frequency Independent Terms

// C20
    c20 = k20/5.0 * ( gmm2e * pow(rerm, 3) * p20m 
                    + gms2e * pow(rers, 3) * p20s );
// C21/S21
    c21 = k21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) );
    s21 = k21/5.0 * ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD) 
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );
// C22/S22
    c22 = k22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) );
    s22 = k22/5.0 * ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );


// C30
    c30 = k30/7.0 * ( gmm2e * pow(rerm, 4) * p30m 
                    + gms2e * pow(rers, 4) * p30s );
// C31/S31
    c31 = k31/7.0 * ( gmm2e * pow(rerm, 4) * p31m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 4) * p31s * cos(llrs[1] * DEG2RAD) );
    s31 = k31/7.0 * ( gmm2e * pow(rerm, 4) * p31m * sin(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 4) * p31s * sin(llrs[1] * DEG2RAD) );
// C32/S32
    c32 = k32/7.0 * ( gmm2e * pow(rerm, 4) * p32m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 4) * p32s * cos(llrs[1] * DEG2RAD * 2.0) );
    s32 = k32/7.0 * ( gmm2e * pow(rerm, 4) * p32m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 4) * p32s * sin(llrs[1] * DEG2RAD * 2.0) );
// C33/S33
    c33 = k33/7.0 * ( gmm2e * pow(rerm, 4) * p33m * cos(llrm[1] * DEG2RAD * 3.0)
                    + gms2e * pow(rers, 4) * p33s * cos(llrs[1] * DEG2RAD * 3.0) );
    s33 = k33/7.0 * ( gmm2e * pow(rerm, 4) * p33m * sin(llrm[1] * DEG2RAD * 3.0)
                    + gms2e * pow(rers, 4) * p33s * sin(llrs[1] * DEG2RAD * 3.0) );

// C40
    c40 = k20p/5.0* ( gmm2e * pow(rerm, 3) * p20m 
                    + gms2e * pow(rers, 3) * p20s );
// C41/S41
    c41 = k21p/5.0* ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) );
    s41 = k21p/5.0* ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );
// C42/S42
    c42 = k22p/5.0* ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) );
    s42 = k22p/5.0* ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );

    
    
    
    stcs[0]  = 0; //c00; 
    stcs[1]  = 0; //c10; 
    stcs[2]  = c20; 
    stcs[3]  = c30; 
    stcs[4]  = c40; 

    stcs[5]  = 0; //c11; 
    stcs[6]  = c21; 
    stcs[7]  = c31; 
    stcs[8]  = c41; 
    stcs[9]  = 0; //s11; 
    stcs[10] = s21; 
    stcs[11] = s31; 
    stcs[12] = s41; 

    stcs[13] = c22; 
    stcs[14] = c32; 
    stcs[15] = c42; 
    stcs[16] = s22; 
    stcs[17] = s32; 
    stcs[18] = s42; 

    stcs[19] = c33; 
    stcs[20] = 0; //c43; 
    stcs[21] = s33; 
    stcs[22] = 0; //s43; 

    stcs[23] = 0; //c44; 
    stcs[24] = 0; //s44; 
  
    
// Frequency Dependent Terms

    c20f = 0; c21f = 0; s21f = 0; c22f = 0; s22f = 0;
//    stfrqdep(info.jdt, info.gmst, &c20f, &c21f, &s21f, &c22f, &s22f);


    stcs[2]  = c20 + c20f; 
    stcs[6]  = c21 + c21f; 
    stcs[10] = s21 + s21f; 
    stcs[13] = c22 + c22f; 
    stcs[16] = s22 + s22f; 

    if(id_perm==1) 
    {
        stcs[2]  = c20 + c20f - c20pt; 
    }
    return 0;

}



double stfrqdep(double jdt, double gmst, double *c20f, double *c21f, double *s21f, double *c22f, double *s22f)
{
    double sets[71][8] = {
     0,5,5,5,6,5,  16.6e-12,  -6.7e-12, 
     0,5,5,5,7,5,  -0.1e-12,   0.1e-12, 
     0,5,6,5,5,4,  -1.2e-12,   0.8e-12, 
     0,5,7,5,5,5,  -5.5e-12,   4.3e-12, 
     0,5,7,5,6,5,   0.1e-12,  -0.1e-12, 
     0,5,8,5,5,4,  -0.3e-12,   0.2e-12, 
     0,6,3,6,5,5,  -0.3e-12,   0.7e-12, 
     0,6,5,4,4,5,   0.1e-12,  -0.2e-12, 
     0,6,5,4,5,5,  -1.2e-12,   3.7e-12, 
     0,6,5,4,6,5,   0.1e-12,  -0.2e-12, 
     0,6,5,6,5,5,   0.1e-12,  -0.2e-12, 
     0,7,3,5,5,5,   0.0e-12,   0.6e-12, 
     0,7,5,3,5,5,   0.0e-12,   0.3e-12, 
     0,7,5,5,5,5,   0.6e-12,   6.3e-12, 
     0,7,5,5,6,5,   0.2e-12,   2.6e-12, 
     0,7,5,5,7,5,   0.0e-12,   0.2e-12, 
     0,8,3,6,5,5,   0.1e-12,   0.2e-12, 
     0,8,5,4,5,5,   0.4e-12,   1.1e-12, 
     0,8,5,4,6,5,   0.2e-12,   0.5e-12, 
     0,9,3,5,5,5,   0.1e-12,   0.2e-12, 
     0,9,5,3,5,5,   0.1e-12,   0.1e-12, 
     1,2,5,7,5,5,  -0.1e-12,   0.0e-12, 
     1,2,7,5,5,5,  -0.1e-12,   0.0e-12, 
     1,3,5,6,4,5,  -0.1e-12,   0.0e-12, 
     1,3,5,6,5,5,  -0.7e-12,   0.1e-12, 
     1,3,7,4,5,5,  -0.1e-12,   0.0e-12, 
     1,4,5,5,4,5,  -1.3e-12,   0.1e-12, 
     1,4,5,5,5,5,  -6.8e-12,   0.6e-12, 
     1,4,7,5,5,5,   0.1e-12,   0.0e-12, 
     1,5,3,6,5,5,   0.1e-12,   0.0e-12, 
     1,5,5,4,4,5,   0.1e-12,   0.0e-12, 
     1,5,5,4,5,5,   0.4e-12,   0.0e-12, 
     1,5,5,6,5,5,   1.3e-12,  -0.1e-12, 
     1,5,5,6,6,5,   0.3e-12,   0.0e-12, 
     1,5,7,4,5,5,   0.3e-12,   0.0e-12, 
     1,5,7,4,6,5,   0.1e-12,   0.0e-12, 
     1,6,2,5,5,6,  -1.9e-12,   0.1e-12, 
     1,6,3,5,4,5,   0.5e-12,   0.0e-12, 
     1,6,3,5,5,5, -43.4e-12,   2.9e-12, 
     1,6,4,5,5,4,   0.6e-12,   0.0e-12, 
     1,6,4,5,5,6,   1.6e-12,  -0.1e-12, 
     1,6,5,3,4,5,   0.1e-12,   0.0e-12, 
     1,6,5,5,3,5,   0.1e-12,   0.0e-12, 
     1,6,5,5,4,5,  -8.8e-12,   0.5e-12, 
     1,6,5,5,5,5, 470.9e-12, -30.2e-12, 
     1,6,5,5,6,5,  68.1e-12,  -4.6e-12, 
     1,6,5,5,7,5,  -1.6e-12,   0.1e-12, 
     1,6,6,4,5,5,   0.1e-12,   0.0e-12, 
     1,6,6,5,4,4,  -0.1e-12,   0.0e-12, 
     1,6,6,5,5,4, -20.6e-12,  -0.3e-12, 
     1,6,6,5,5,6,   0.3e-12,   0.0e-12, 
     1,6,6,5,6,4,  -0.3e-12,   0.0e-12, 
     1,6,7,3,5,5,  -0.2e-12,   0.0e-12, 
     1,6,7,3,6,5,  -0.1e-12,   0.0e-12, 
     1,6,7,5,5,5,  -5.0e-12,   0.3e-12, 
     1,6,7,5,6,5,   0.2e-12,   0.0e-12, 
     1,6,8,5,5,4,  -0.2e-12,   0.0e-12, 
     1,7,3,6,5,5,  -0.5e-12,   0.0e-12, 
     1,7,3,6,6,5,  -0.1e-12,   0.0e-12, 
     1,7,5,4,4,5,   0.1e-12,   0.0e-12, 
     1,7,5,4,5,5,  -2.1e-12,   0.1e-12, 
     1,7,5,4,6,5,  -0.4e-12,   0.0e-12, 
     1,8,3,5,5,5,  -0.2e-12,   0.0e-12, 
     1,8,5,3,5,5,  -0.1e-12,   0.0e-12, 
     1,8,5,5,5,5,  -0.6e-12,   0.0e-12, 
     1,8,5,5,6,5,  -0.4e-12,   0.0e-12, 
     1,8,5,5,7,5,  -0.1e-12,   0.0e-12, 
     1,9,5,4,5,5,  -0.1e-12,   0.0e-12, 
     1,9,5,4,6,5,  -0.1e-12,   0.0e-12, 
     2,4,5,6,5,5,  -0.3e-12,   0.0e-12, 
     2,5,5,5,5,5,  -1.2e-12,   0.0e-12  
    };

    double doodarg[6], ang, c20 = 0, c21 = 0, s21 = 0, c22 = 0, s22 = 0;
    int i, nsets = 71, argn[6], ncon = 1;

    for (i=0;i<nsets;i++)
    {
        argn[0] = (int)sets[i][0];
        argn[1] = (int)sets[i][1] - 5;
        argn[2] = (int)sets[i][2] - 5;
        argn[3] = (int)sets[i][3] - 5;
        argn[4] = (int)sets[i][4] - 5;
        argn[5] = (int)sets[i][5] - 5;
        
//        DOODSN(&info[num].jdt, &info[num].gmst, argn, &ncon, doodarg, &ang);
//        DOODSN(&jdt, &gmst, otfes[i].argn, &ncon, doodarg, &ang);

//  C20 correction: Long period tidal constituent

        if(argn[0]==0) 
        {
           c20 = c20 + sets[i][6]*cos(ang) - sets[i][7]*sin(ang);
        }

//  C21/S21 correction: Diurnal period tidal constituent

        if(argn[0]==1) 
        {
            c21 = c21 + sets[i][6]*sin(ang) + sets[i][7]*cos(ang);
            s21 = s21 + sets[i][6]*cos(ang) - sets[i][7]*sin(ang);
        }

//  C22/S22 correction: Semi-diurnal period tidal constituent

        if(argn[0]==2) 
        {
            c22 = c22 + sets[i][6]*cos(ang);
            s22 = s22 - sets[i][6]*sin(ang);
        }
    }

    *c20f = c20;
    *c21f = c21;
    *s21f = s21;
    *c22f = c22;
    *s22f = s22;



    return 0;

}



double accel_gravt (double *tjd, double *xic, double *a4)
{
    int n, lps;
    double GM, radius, pi[3], pe[3], llr[3], c_ie[9], c_ei[9], ae[3], ai[3];
    double jd0, tt, utc, te[9], tx[3], ty[3], tz[3],
        vx[3] = {1,0,0}, vy[3] = {0,1,0}, vz[3] = {0,0,1};

    GM = 398600.44150E+09;
    radius = 6378136.3;

 
    for (n = 0; n < 3; n++)
    {
        pi[n] = xic[n] * AU;
    }

    jd0 = tjd[0];
    tt = tjd[1] * 86400.0;

    lps = getlps (jd0 + tt/86400.0);
    utc = tt - (lps + 32.184);

//    printf("%f\t%f\n", jd0, utc);
    itrf2gcrf(jd0, utc, vx, tx);
    itrf2gcrf(jd0, utc, vy, ty);
    itrf2gcrf(jd0, utc, vz, tz);


    for (n = 0; n < 3; n++)
    {
        c_ie[n*3] = tx[n];
        c_ie[n*3+1] = ty[n];
        c_ie[n*3+2] = tz[n];
    }

    mt(c_ie, 3, 3, c_ei);

    brmul(c_ei, pi, 3, 3, 1, pe);  
    
    xyz2llh(pe, llr);

        
    cs2acc (llr, COEFG, GM, radius, NMAX, ae);
//    printf("ae = %f\t%f\t%f\n", ae[0], ae[1], ae[2]);
//    exit(0);

    brmul(c_ie, ae, 3, 3, 1, ai);
    
    for (n = 0; n < 3; n++)
    {
        a4[n] = ai[n] / AU * 86400 * 86400;
    }

    return 0;

}








double cs2ada (double *llr, double *cs, int nmax, double *ae,
                int part, double *dadre, int flagdadcs)
{
    int n, m, k, l, ind, ic, is, label;

    double slat, clat, slon, clon, sclt, cclt, *cosml, *sinml, 
        *aprn, *pbar, *pbar1, *pbar2, *pt, *pt1, *pt2, lat, lon, r, vi, t, 
        gm, a, an[3], c_en[9], c_ne[9], dadrn[9], dadrne[9], dadre1[9], 
        dadre2[9], dcdrn[27], dcdre[27], dcdxe[9], dcdye[9], dcdze[9], 
        dadrecx[3], dadrecy[3], dadrecz[3], dadrec[9], dadrea[9];

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    gm = GMCT;
    a = RCT; 

    lat = llr[0];
    lon = llr[1];
    r = llr[2];

    slat = sin(lat * DEG2RAD);
    clat = cos(lat * DEG2RAD);
    slon = sin(lon * DEG2RAD);
    clon = cos(lon * DEG2RAD);
    sclt = clat;
    cclt = slat;

    cosml = (double *) calloc ( nmax + 1, sizeof(double)); //cos(m*lamta)
    sinml = (double *) calloc ( nmax + 1, sizeof(double)); //sin(m*lamta)
    aprn = (double *) calloc ( nmax + 1, sizeof(double));  //sin(m*lamta)
    pbar = (double *) calloc ( nmax + 1, sizeof(double));
    pbar1 = (double *) calloc ( nmax + 1, sizeof(double));
    pbar2 = (double *) calloc ( nmax + 1, sizeof(double));
    pt  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    pt1  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    pt2  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

    cosml[0] = 1; sinml[0] = 0;
    for (m = 1; m <= nmax; m++)
    {
        cosml[m] = cos(m * lon * DEG2RAD);
        sinml[m] = sin(m * lon * DEG2RAD);
    }
    
    for (n = 0; n <= nmax; n++)
    {
        aprn[n] = pow (a / r, n) * gm / r;
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    an[0] = 0; an[1] = 0; an[2] = 0;
    t = cclt; vi = 0;
    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;
        lgdr2(t, nmax, m, pbar, pbar1, pbar2);
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
                n = k;
                ic = n;
                is = 0;
            }
            else
            {
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                n = k + m;
                ic = ind + n - m;
                is = ind + n - m + l;
            }
            pt[ic] = aprn[n] * pbar[k] * cosml[m];
            pt[is] = aprn[n] * pbar[k] * sinml[m];
            pt1[ic] = aprn[n] * pbar1[k] * cosml[m];
            pt1[is] = aprn[n] * pbar1[k] * sinml[m];
            pt2[ic] = aprn[n] * pbar2[k] * cosml[m];
            pt2[is] = aprn[n] * pbar2[k] * sinml[m];

            vi = vi + pt[ic] * cs[ic] +  pt[is] * cs[is];
            an[0] = an[0] + pt1[ic] * cs[ic] + pt1[is] * cs[is];
            an[1] = an[1] - m *  pt[is] * cs[ic] + m *  pt[ic] * cs[is];          
            an[2] = an[2] + (n+1) * pt[ic] * cs[ic] + (n+1) * pt[is] * cs[is];
//            an[2] = an[2] + (n+1) * ( pt[ic] * cs[ic] +  pt[is] * cs[is]);
        }
    }


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    free (pbar);
    free (pbar1);
    free (pbar2);

    free (cosml);
    free (sinml);
    free (aprn);


    an[0] = - an[0] / r;
    an[1] = + an[1] / r / sclt;
    an[2] = + an[2] / r;


    c_ne[0] = - slat * clon;      //nsys: North-East-Down
    c_ne[1] = - slon;
    c_ne[2] = - clat * clon;
    c_ne[3] = - slat * slon;
    c_ne[4] = clon;
    c_ne[5] = - clat * slon;
    c_ne[6] = clat;
    c_ne[7] = 0;
    c_ne[8] = - slat;



    brmul(c_ne, an, 3, 3, 1, ae);  //from n-sys to e-sys

    for (n = 0; n < 9; n++)
        dadre[n] = 0;

    if (part == 0)
    {
        free (pt);
        free (pt1);
        free (pt2);
        return vi;
    }



    for (n = 0; n < 9; n++)
        dadrn[n] = 0;
    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
                n = k;
                ic = n;
                is = 0;
            }
            else
            {
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                n = k + m;
                ic = ind + n - m;
                is = ind + n - m + l;
            }
            
            dadrn[0] += pt2[ic] * cs[ic] + pt2[is] * cs[is];
            dadrn[3] += m * ( - pt[is] * cs[ic] + pt[ic] * cs[is]) * cclt / sclt / sclt
                      - m * ( - pt1[is] * cs[ic] + pt1[ic] * cs[is]) / sclt;
            dadrn[6] -= (n+1) * (pt1[ic] * cs[ic] + pt1[is] * cs[is]);

            dadrn[1] -= m * ( - pt1[is] * cs[ic] + pt1[ic] * cs[is]) / sclt;
            dadrn[4] -= m * m * (pt[ic] * cs[ic] + pt[is] * cs[is]) / sclt / sclt;
            dadrn[7] += m * (n+1) * ( - pt[is] * cs[ic] + pt[ic] * cs[is]) / sclt;

            dadrn[2] -= (n+2) * (pt1[ic] * cs[ic] + pt1[is] * cs[is]);
            dadrn[5] += m * (n+2) * ( - pt[is] * cs[ic] + pt[ic] * cs[is]) / sclt;
            dadrn[8] += (n+2) * (n+1)* (pt[ic] * cs[ic] + pt[is] * cs[is]);
        }
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    for (n = 0; n < 9; n++)
        dadrn[n] = dadrn[n] / r / r;

    mt(c_ne, 3, 3, c_en);
    brmul(dadrn, c_en, 3, 3, 3, dadrne);			
    brmul(c_ne, dadrne, 3, 3, 3, dadrea);	




    for (n = 0; n < 27; n++)
        dcdrn[n] = 0;

    dcdrn[0]  = - sclt * clon;  dcdrn[1]  = cclt * slon / sclt;
    dcdrn[3]  = 0;              dcdrn[4]  = - clon / sclt;
    dcdrn[6]  = cclt * clon;/**/dcdrn[7]  = slon;

    dcdrn[9]  = - sclt * slon;  dcdrn[10] = - cclt * clon / sclt;
    dcdrn[12] = 0;              dcdrn[13] = - slon / sclt;
    dcdrn[15] = cclt * slon;    dcdrn[16] = - clon;

    dcdrn[18] = - cclt;         dcdrn[19] = 0;
    dcdrn[21] = 0;              dcdrn[22] = 0;
    dcdrn[24] = - sclt;         dcdrn[25] = 0;

    for (n = 0; n < 27; n++)
        dcdrn[n] = dcdrn[n] / r;

    brmul(dcdrn, c_en, 9, 3, 3, dcdre);		

    for (n = 0; n < 9; n++)	
    {
        dcdxe[n] = dcdre[n*3];
        dcdye[n] = dcdre[n*3+1];
        dcdze[n] = dcdre[n*3+2];
//        printf ("%e\n", dcdze[n]);
    }


    brmul(dcdxe, an, 3, 3, 1, dadrecx);	
    brmul(dcdye, an, 3, 3, 1, dadrecy);  
    brmul(dcdze, an, 3, 3, 1, dadrecz);  

    for (n = 0; n < 3; n++)	
    {
        dadrec[n*3] = dadrecx[n]; 
        dadrec[n*3+1] = dadrecy[n]; 
        dadrec[n*3+2] = dadrecz[n];
    }




    for (n = 0; n <= 8; n++)
    {
        dadre[n] = dadrec[n] + dadrea[n]; 
    }


    if (flagdadcs == 0)
    {
        free (pt);
        free (pt1);
        free (pt2);
        return vi;
    }




// for dadcs_nm
//




    for (k = 0; k < MGCS; k ++)
    {   
        n = CSinfo[k].n; m = CSinfo[k].m; label = CSinfo[k].cs;

        if (m == 0)
        {
            ic = n;
            CSinfo[k].dadcsn[0] = - pt1[ic] / r;
            CSinfo[k].dadcsn[1] = 0;          
            CSinfo[k].dadcsn[2] = (n+1) * pt[ic] / r;
        }
        else 
        {
            l = nmax - m + 1;
            ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
            ic = ind + n - m;
            is = ind + n - m + l;
            if (label == 1)
            {
                CSinfo[k].dadcsn[0] = - pt1[ic] / r; 
                CSinfo[k].dadcsn[1] = - m *  pt[is] / r / sclt;          
                CSinfo[k].dadcsn[2] = (n+1) * pt[ic] / r;
            }

            if (label == -1)
            {
                CSinfo[k].dadcsn[0] = - pt1[is] / r;
                CSinfo[k].dadcsn[1] =   m *  pt[ic] / r / sclt;          
                CSinfo[k].dadcsn[2] = (n+1) * pt[is] / r;
            }
        }

        brmul(c_ne, CSinfo[k].dadcsn, 3, 3, 1, CSinfo[k].dadcse);	
    }






    free (pt);
    free (pt1);
    free (pt2);
    return vi;

    

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
































/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double cs2acc (double *llr, double *cs, double gm, double a, int nmax, 
               double *acc)
{
    int n, m, k, l, ind;

    double sinf, cosf, sinlon, coslon, sincolat, coscolat, *cosml, *sinml, 
        *aprn, *pbar, *pbar1, *pbar2, accn[3], c_ei[9], c_en[9], c_in[9],
        *pt, *ptt, lat, lon, r, vi, dvdr, dvdcolat, dvdlon, t;


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    lat = llr[0];
    lon = llr[1];
    r = llr[2];

    sinf = sin(lat * DEG2RAD);
    cosf = cos(lat * DEG2RAD);
    sinlon = sin(lon * DEG2RAD);
    coslon = cos(lon * DEG2RAD);
    sincolat = cosf;
    coscolat = sinf;

//    #pragma omp parallel private(cosml, sinml, aprn, pbar, pbar1, pbar2, n, m, k, l, ind, sinf, cosf)
    cosml = (double *) calloc ( nmax + 1, sizeof(double)); //cos(m*lamta)
    sinml = (double *) calloc ( nmax + 1, sizeof(double)); //sin(m*lamta)
    aprn = (double *) calloc ( nmax + 1, sizeof(double));  //sin(m*lamta)
    pbar = (double *) calloc ( nmax + 1, sizeof(double));
    pbar1 = (double *) calloc ( nmax + 1, sizeof(double));
    pbar2 = (double *) calloc ( nmax + 1, sizeof(double));
    pt  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
    ptt  = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));

    for (m = 0; m <= nmax; m++)
    {
        cosml[m] = cos(m * lon * DEG2RAD);
        sinml[m] = sin(m * lon * DEG2RAD);
    }
    
    for (n = 0; n <= nmax; n++)
    {
        aprn[n] = pow (a / r, n) * gm / r;
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    t = coscolat; vi = 0; dvdlon = 0; dvdcolat = 0; dvdr = 0;
    for (m = 0; m <= nmax; m ++)
    {
        l = nmax - m + 1;
        lgdr2(t, nmax, m, pbar, pbar1, pbar2);
//        lgdr(t, nmax, m, pbar);
        for (k = 0; k < l; k++)
        {
            if (m==0)
            {
//                ind = 0;
                n = k + m;
                pt[k] = aprn[n] * pbar[k];
                ptt[k] = aprn[n] * pbar1[k];
                vi = vi + pt[k] * cs[k];
            
//                if (n>=2)
                {
                    dvdr = dvdr + (n+1) * pt[k] * cs[k];
                    dvdcolat = dvdcolat + ptt[k] * cs[k];
                }
            }
            else
            {
                ind = nmax + 1 + (2 * nmax - m + 2) * (m - 1);
                n = k + m;
                pt[ind + n - m] = aprn[n] * pbar[k] * cosml[m];
                pt[ind + n - m + l] = aprn[n] * pbar[k] * sinml[m];
                ptt[ind + n - m] = aprn[n] * pbar1[k] * cosml[m];
                ptt[ind + n - m + l] = aprn[n] * pbar1[k] * sinml[m];
                vi = vi + pt[ind + n - m] * cs[ind + n - m];
                vi = vi + pt[ind + n - m + l] * cs[ind + n - m + l];

                dvdcolat = dvdcolat + ptt[ind + n - m] * cs[ind + n - m];
                dvdcolat = dvdcolat + ptt[ind + n - m + l] * cs[ind + n - m + l];

                dvdlon = dvdlon - m * pt[ind + n - m + l] * cs[ind + n - m];
                dvdlon = dvdlon + m * pt[ind + n - m] * cs[ind + n - m + l];          

                dvdr = dvdr + (n+1) * pt[ind + n - m] * cs[ind + n - m];
                dvdr = dvdr + (n+1) * pt[ind + n - m + l] * cs[ind + n - m + l];
            }
        }
    }

//    dvdcolat = - dvdcolat * sincolat; //tmd!!
    dvdcolat = dvdcolat;
    dvdlon = + dvdlon;
    dvdr = - dvdr / r;


    accn[0] = - dvdcolat / r;
    accn[1] = + dvdlon / r / sincolat;
    accn[2] = - dvdr;


    c_en[0] = - sinf * coslon;      //from fixed to up-east-north system: rmat
    c_en[1] = - sinlon;
    c_en[2] = - cosf * coslon;
    c_en[3] = - sinf * sinlon;
    c_en[4] = coslon;
    c_en[5] = - cosf * sinlon;
    c_en[6] = cosf;
    c_en[7] = 0;
    c_en[8] = - sinf;


//    mt(info[num].c_ie, 3, 3, info[num].c_ei);
//    brmul (info[num].c_ei, c_en, 3, 3, 3, c_in);  //inertial to fixed matrix gmat = rmat*tbt

//    brmul(c_in, accn, 3, 3, 1, acc);  //from fixed acc to inertial acc

    brmul(c_en, accn, 3, 3, 1, acc);  //from fixed acc to inertial acc

//    *v = vi;
//    *dvdt = - ANGVEL * dvdlon;

    free (pbar);
    free (pbar1);
    free (pbar2);
    free (pt);
    free (ptt);

    free (cosml);
    free (sinml);
    free (aprn);
    return 1;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/




double lgdr(double t, int nmax, int m, double *pbar)

/*
! THIS CALCULATES THE FULLY NORMALIZED LEGENDRE FUNCTION WITH GIVEN ORDER(M),
! MAXIMUM DEGREE (NMAX), AND GIVEN EVALUATION POINT, T (COSINES OF COLATITUDE).
! THIS RETURNS ALL Pn,m, P'n,m, AND P''n,m (m=<n<=Nmax).
! THE RECURSION FORMULAR FOR THE FUNCTION ITSELF IS GIVEN IN JEKELI(1996).
! THE RECURSION FORMULAR FOR THE 1ST DERIVATIVE IS GIVEN IN TSCHERNING, ET AL(1983).
! THE FORMULAR FOR THE 2ND DERIVATIVE IS FROM THE ASSOCIATE LEGENDRE EQUATION.
! NOTE : EQUATIONS GIVEN IN TSCHERNING, ET AL(1983) HAVE ERRATA.
!
! S.C. Han, 1/24/01 (MODIFIED FOR CRAY T94 2/13/01)
!
*/
{
    int i;
//REAL*8 :: PBAR(NMAX-M+1),PBAR1(NMAX-M+1),PBAR2(NMAX-M+1),T,P00,P11,C,D
    double p00, p11, c, d;
//! THE FULLY NORMALIZED ASSOCIATED LEGENDRE FUNCTION
//! Pm,m : JEKEIL (A.3c) & (A.3d) , P'm,m : TSCHERNING (7)

    p00 = 1.0; 
    p11 = sqrt (3.0*(1.0-t*t));
    if (m>=1)
    {
        pbar[0] = p11; 

        for (i = 2; i <= m; i++)
        {
            pbar[0] = sqrt((2.0*i+1.0)/(2.0*i)*(1.0-t*t))*pbar[0];
        }
    }
    else 
    {
        pbar[0]=p00; 
    }

    if (nmax - m + 1 >= 2)
    {
        pbar[1] = sqrt(2.0*m +3.0) * t * pbar[0];
    }

    for(i = 3; i <= nmax-m+1; i++)
    {
        c=((2.0*m+2.0*i-3.0) * (2.0*m + 2.0*i-1.0)) / ((i-1.0)*(2.0*m+i-1.0));
        d=((2.0*m+2.0*i-1.0)*(2.0*m+i-2.0)*(i-2.0)) 
            / ((2.0*m+2.0*i-5.0)*(i-1.0)*(2.0*m+i-1.0));
        pbar[i-1] = sqrt(c)*t*pbar[i-2] - sqrt(d) * pbar[i-3];
    }

    return 0;
}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double lgdr2(double t, int nmax, int m, 
             double *pbar, double *pbar1, double *pbar2)

/*
! THIS CALCULATES THE FULLY NORMALIZED LEGENDRE FUNCTION WITH GIVEN ORDER(M),
! MAXIMUM DEGREE (NMAX), AND GIVEN EVALUATION POINT, T (COSINES OF COLATITUDE).
! THIS RETURNS ALL Pn,m, P'n,m, AND P''n,m (m=<n<=Nmax).
! THE RECURSION FORMULAR FOR THE FUNCTION ITSELF IS GIVEN IN JEKELI(1996).
! THE RECURSION FORMULAR FOR THE 1ST DERIVATIVE IS GIVEN IN TSCHERNING, ET AL(1983).
! THE FORMULAR FOR THE 2ND DERIVATIVE IS FROM THE ASSOCIATE LEGENDRE EQUATION.
! NOTE : EQUATIONS GIVEN IN TSCHERNING, ET AL(1983) HAVE ERRATA.
!
! S.C. Han, 1/24/01 (MODIFIED FOR CRAY T94 2/13/01)
!
*/
{
    int i;
//REAL*8 :: PBAR(NMAX-M+1),PBAR1(NMAX-M+1),PBAR2(NMAX-M+1),T,P00,P11,C,D
    double p00, p11, c, d;
//! THE FULLY NORMALIZED ASSOCIATED LEGENDRE FUNCTION
//! Pm,m : JEKEIL (A.3c) & (A.3d) , P'm,m : TSCHERNING (7)

    p00 = 1.0; 
    p11 = sqrt (3.0*(1.0-t*t));
    if (m>=1)
    {
        pbar[0] = p11; 
        pbar1[0] = sqrt(3.0) * t;

        for (i = 2; i <= m; i++)
        {
            pbar1[0] = sqrt((2.0*i+1.0)/(2.0*i))*(sqrt(1.0-t*t)*pbar1[0]+t*pbar[0]);
//            pbar1[0] = sqrt((2.0*i+1.0)/(2.0*i))*(sqrt(1.0-t*t)*pbar1[0]+t*pbar[0]/(-sqrt(1.0-t*t)));
            pbar[0] = sqrt((2.0*i+1.0)/(2.0*i)*(1.0-t*t))*pbar[0];
        }
    }
    else 
    {
        pbar[0]=p00; 
        pbar1[0]=0.0;
    }

//    ! Pm+1,m : JEKEIL (A.3b)

    if (nmax - m + 1 >= 2)
    {
        pbar[1] = sqrt(2.0*m +3.0) * t * pbar[0];
    }

//  ! Pn,m (n>=m+2) : JEKEIL (A.3a)

    for(i = 3; i <= nmax-m+1; i++)
    {
        c=((2.0*m+2.0*i-3.0) * (2.0*m + 2.0*i-1.0)) / ((i-1.0)*(2.0*m+i-1.0));
        d=((2.0*m+2.0*i-1.0)*(2.0*m+i-2.0)*(i-2.0))/((2.0*m+2.0*i-5.0)*(i-1.0)*(2.0*m+i-1.0));
        pbar[i-1] = sqrt(c)*t*pbar[i-2] - sqrt(d) * pbar[i-3];
    }

//  ! THE FULLY NORMALIZED ASSOCIATED LEGENDRE FUNCTION - 1ST DERIVATIVE
//  ! P'n,m (n>=m+1) : TSCHERNING (8)
    for (i=2; i<=nmax-m+1; i++)
    {
        c = 1.0/sqrt(1.0-t*t)*t*(m+i-1);
        d = 1.0/sqrt(1.0-t*t)*sqrt((((m+i-1)*(m+i-1)-m*m)*(2.0*(m+i-1)+1.0))/(2.0*(m+i-1)-1.0));
//!! found it different from TSCHERNING (8),dcl-2010-2-14
//!! Jianbin confirms code is correct, dcl-2010-2-15
//!!      D=1D0/SQRT(1D0-T**2)/SQRT((((M+I-1)**2-M**2)*(2D0*(M+I-1)+1D0))/(2D0*(M+I-1)-1D0))
        pbar1[i-1] = c * pbar[i-1] - d * pbar[i-2];
    }

//! THE FULLY NORMALIZED ASSOCIATED LEGENDRE FUNCTION - 2ND DERIVATIVE
//! P''n,m (n>=m) : ASSOCIATE LEGENDRE EQUATION (2ND ORDER DIFFERENTIAL EQN.)

    for (i=1;i<=nmax-m+1;i++)
    {
        pbar2[i-1] = (-t/sqrt(1.0-t*t)) * pbar1[i-1] 
            - ((m+i-1)*(m+i)-m*m/(1.0-t*t)) * pbar[i-1];
    }
    return 0;
}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/





double accel_slrad (double *tjd, double *xic, double *acc)
{
    double j, c1, ap, m, rsp, usp[3], xis[6], xsc[6], f, 
        xist[9], unit[9], rsp3;
    short int n, sun;
    short int ssbary = 11;

    sun = 10;
    unit[0] = 1; unit[1] = 0; unit[2] = 0;
    unit[3] = 0; unit[4] = 1; unit[5] = 0;
    unit[6] = 0; unit[7] = 0; unit[8] = 1;

    planet_ephemeris (tjd, sun, CT, &xsc[0], &xsc[3]);
    for (n = 0; n <= 5; n++)
    {
        xis[n] = xic[n] - xsc[n];
    }
    rsp = sqrt (xis[0] * xis[0] + xis[1] * xis[1] + xis[2] * xis[2]);
    usp[0] = xis[0] / rsp;
    usp[1] = xis[1] / rsp;
    usp[2] = xis[2] / rsp;

    j  = 1352.5;   //kg/s3
//    j  = 1359.4;   //kg/s3
//    m  = SATMASS;     //kg
//    ap = SATAREA;        //m2
    c1 = j / C * 1 * 1;   //kg/s2/m*au*au
    f  = c1 * AMR  / rsp / rsp;    
//    f  = c1 * ap / m  / rsp / rsp;    
//kg/s2/m*au*au * m2 / kg / au  / au = m/s2
    f = f / AU * 86400.0 * 86400.0;

    acc[0] = f * usp[0];
    acc[1] = f * usp[1];
    acc[2] = f * usp[2]; 


    return 0;
}





double accel_nbody (double *tjd, double *xic, double *fnt, double *fgr)
{
    int n;
    short int ssbary = 11;

    double xcb[6], xib[6], fnti[3], fntb[3], fgri[3], fgrb[3];

    planet_ephemeris (tjd, CT, ssbary, &xcb[0], &xcb[3]);
    force_bcrs (tjd, xcb, CT, fntb, fgrb);
    for (n = 0; n <= 5; n++)
    {
        xib[n] = xic[n] + xcb[n];
    }
    force_bcrs (tjd, xib, CT, fnti, fgri);
//    force_bcrs (tjd, xib, 99, fnti, fgri);
    for (n = 0; n <= 2; n++)
    {
        fnt[n] = fnti[n] - fntb[n];
        fgr[n] = fgri[n] - fgrb[n];

    }

    return 0;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double force_bcrs (double *jd, double *xi, short int exclude, 
                   double *fnt, double *fgr)
{
    double xj[11][6], xij[11][6], rij[11], xjk[6], rjk, 
        xddj[3], sumil, sumjk, sdi2, sdj2, rdirdj, rrrdr2, rjirdd, gm[11], GMDE[11],
        rij5, rij3, xijt[9], gra, grb, beta, gamma, unit[9],gm2de;
    short int ssbary, l, k, j, n, flag_gr;
    
    gm[0] =   2.203208082807623e+13;
    gm[1] =      3.248586038641429e+14;
    gm[2] =     398600.44150E+09;
    gm[3] =     4.28283719012840e+13;
    gm[4] =      1.267127698227696e+17;
    gm[5] =     3.794062664949063e+16;
    gm[6] =      5.794549096929744e+15;
    gm[7] =     6.836534169987595e+15;
    gm[8] =    9.816009029289940e+11;
    gm[9] =      4.902801056E+12;
    gm[10] =      1.32712442076e20;

    gm2de = 86400.0 * 86400.0 / AU / AU / AU;
    for (n = 0; n <= 10; n++)
        GMDE[n] = gm[n] * gm2de;

    ssbary = 11;
//    ssbary = 10;
    gamma = 1.0;
    beta = 1.0;
    unit[0] = 1; unit[1] = 0; unit[2] = 0;
    unit[3] = 0; unit[4] = 1; unit[5] = 0;
    unit[6] = 0; unit[7] = 0; unit[8] = 1;

    for (j = 0; j <= 10; j++)
    {
        planet_ephemeris (jd, j, ssbary, &xj[j][0], &xj[j][3]);
        for (n = 0; n < 6; n++)
        {
            xij[j][n] = xi[n] - xj[j][n];
        }
        rij[j] = sqrt (xij[j][0] * xij[j][0] 
            + xij[j][1] * xij[j][1] + xij[j][2] * xij[j][2]);
    }
    
    for (n = 0; n < 3; n ++)
        fnt[n] = 0;
    for (j = 0; j <= 10; j++)
    {
        if (j == exclude)
            continue;
        for (n = 0; n < 3; n++)
            fnt[n] = fnt[n] 
            - GMDE[j] / (rij[j] * rij[j] * rij[j]) * xij[j][n];
    }


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    sdi2 = xi[3] * xi[3] + xi[4] * xi[4] + xi[5] * xi[5];
    sumil = 0;
    for (l = 0; l < 11; l ++)
    {
        if ( l == exclude)
            continue;
        sumil = sumil + GMDE[l] / rij[l];
    }

    for (n = 0; n < 3; n ++)
        fgr[n] = 0;
    for (j = 0; j < 11; j ++)
    {
        if (j == exclude)
            continue;
        sumjk = 0;
        for (n = 0; n < 3; n ++)
            xddj[n] = 0;
        for (k = 0; k < 11; k ++)
        {
            if (k == j)	
                continue;	//k!=j
            for (n = 0; n < 3; n++)
                xjk[n] = xj[j][n] - xj[k][n];
            rjk = sqrt (xjk[0] * xjk[0] + xjk[1] * xjk[1] + xjk[2] * xjk[2]);
            sumjk = sumjk + GMDE[k] / rjk;
            for (n = 0; n < 3; n ++)
                xddj[n] = xddj[n] - GMDE[k] / (rjk * rjk * rjk) * xjk[n];
        }
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
        sdj2 = xj[j][3] * xj[j][3] + xj[j][4] * xj[j][4] 
            + xj[j][5] * xj[j][5];
        rdirdj = xi[3] * xj[j][3] + xi[4] * xj[j][4] + xi[5] * xj[j][5];
        rrrdr2 = pow( ( xij[j][0] * xj[j][3] + xij[j][1] * xj[j][4] 
            + xij[j][2] * xj[j][5]) / rij[j], 2);
        rjirdd = - ( xij[j][0] * xddj[0] + xij[j][1] * xddj[1] 
            + xij[j][2] * xddj[2]);
        
        gra = - 2 * (beta + gamma) * sumil - (2 * beta -1) * sumjk 
            + gamma * sdi2 + (1 + gamma) * sdj2
            - 2 * (1 + gamma) * rdirdj - 1.5 * rrrdr2 + 0.5 * rjirdd;

        grb = xij[j][0] * ((2+2*gamma) * xi[3] - (1+2*gamma) * xj[j][3])
            + xij[j][1] * ((2+2*gamma) * xi[4] - (1+2*gamma) * xj[j][4])
            + xij[j][2] * ((2+2*gamma) * xi[5] - (1+2*gamma) * xj[j][5]);

        for (n = 0; n < 3; n ++)
        {
            fgr[n] = fgr[n] 
                + GMDE[j] / (rij[j] * rij[j] * rij[j]) 
                * ( - xij[j][n]) * gra / C_AUDAY / C_AUDAY
                + GMDE[j] / (rij[j] * rij[j] * rij[j]) 
                * xij[j][n + 3] * grb / C_AUDAY / C_AUDAY
                + GMDE[j] / rij[j] * (3 + 4 * gamma) * 0.5   
                * xddj[n] / C_AUDAY / C_AUDAY;
        }
    }
    return 1;
}










double modvect (double *v)
{
    return  sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}


double dotvect (double *v1, double *v2)
{
    return  v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}


void crsvect (double *v1, double *v2, double *v)
{
    v[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v[2] = v1[0] * v2[1] - v1[1] * v2[0]; 
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* chosephase 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double chosephase (double sinvalue, double cosvalue)
{
    double sv = sinvalue, cv = cosvalue;
    if (sv >= 0 && cv >= 0) 
        return (asin (sv));
    if (sv > 0 && cv < 0) 
        return (acos (cv));
    if (sv < 0 && cv < 0) 
        return ( - asin (sv) + TWOPI / 2.0);
    else 
        return (asin (sv) + TWOPI);
}






/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/****************************************************************************/
/*                                                                          */
/*		Functions for Runge-Kutta integrator                                */
/*                                                                          */
/*      Version:    2009-9-8                                                */
/*                                                                          */
/*      Copyright (c) 2009 shangkun@shao.ac.cn All Right Reserved           */
/*                                                                          */
/****************************************************************************/

/*
  Version: 2009-9-8 
  Version: 2009-9-13 integrate forwards & backwards
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double rkf78 (double jd, double t, double h, double *x, int dim, 
              void (*fun)(int, double, double,double *,double *))
//              double (*fun)(double, double,double *,double *))
/*
    purpose: auto-adjusted Runge-Kutta-Ful... integrator 
    input:  double h				integration step
            double t                integrate from t to t+h
            double *x               x(t)
            int dim	                dim(x)
            double err              tolerance of step control
            double (*fun)()			right(force) function
    output: double *x				x(t+h)
    return: h                       new step after adjustment
*/
{	
    int i, j, n, flag = 0;
    double *y, *k, *f, d = 0, tn;
    double a[13] = { 0, 2.0/27, 1.0/9, 1.0/6, 5.0/12, 1.0/2, 5.0/6, 1.0/6, 
        2.0/3, 1.0/3, 1.0, 0, 1.0 };
    double c[13] = { 0, 0, 0, 0, 0, 34.0/105, 9.0/35, 9.0/35, 9.0/280, 
        9.0/280, 0, 41.0/840, 41.0/840 };
    double b[13][12] = 
    {
        {0},
        {2.0/27},
        {1.0/36,1.0/12},
        {1.0/24,0,1.0/8},
        {5.0/12,0,-25.0/16,25.0/16},
        {1.0/20,0,0,1.0/4,1.0/5},
        {-25.0/108,0,0,125.0/108,-65.0/27,125.0/54},
        {31.0/300,0,0,0,61.0/225,-2.0/9,13.0/900},
        {2.0,0,0,-53.0/6,704.0/45,-107.0/9,67.0/90,3.0},
        {-91.0/108,0,0,23.0/108,-976.0/135,311.0/54,-19.0/60,17.0/6,-1.0/12},
        {2383.0/4100,0,0,-341.0/164,4496.0/1025,-301.0/82,2133.0/4100,
        45.0/82,45.0/164,18.0/41},
        {3.0/205,0,0,0,0,-6.0/41,-3.0/205,-3.0/41,3.0/41,6.0/41},
        {-1777.0/4100,0,0,-341.0/164,4496.0/1025,-289.0/82,2193.0/4100,
        51.0/82,33.0/164,12.0/41,0,1.0}
    };
    
    y = (double *) calloc (dim, sizeof(double));
    k = (double *) calloc (dim*13, sizeof(double));
    f = (double *) calloc (dim, sizeof(double));
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    
    do
    {	
        for (i = 0; i <= 12; i++)
        {			
            tn = t + a[i] * h;
            for (n = 0; n <= dim - 1; n++)
            {	
                y[n] = x[n];
                for (j = 0; j <= i-1; j++)
                    y[n] = y[n] + h * b[i][j] * k[n*13+j];
            }
            fun (dim, jd, tn, y, f);
//            fun (jd, tn, y, f);
            for (n = 0; n <= dim - 1; n++)
            {
                k[n*13+i] = f[n];
            }
        }
        d = 0;
        for (n = 0; n <= dim - 1; n++)
        {
            d = d + fabs (41.0 / 840 * (k[n*13+0] + k[n*13+10] 
                - k[n*13+11] - k[n*13+12]) * h);
        }
		
        flag = 0;
    }while (flag == 1);

    for (n = 0; n <= dim - 1; n++)
    {
        for (i = 0; i <= 12; i++)
            x[n] = x[n] + h * c[i] * k[n*13+i];	
    }
	
    free (y);
    free (f);
    free (k);
    return h;
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* mt - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void mt (double *a, int m, int n, double *b)
{
    int i, j;
    for (i = 0; i <= m - 1; i++)
    {
        for (j = 0; j <= n - 1; j++)
            b[j * m + i] = a[i * n + j];
    }
    return;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* brmul - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void brmul (double *a, double *b, int m,int n, int k,double *c)
{ 
    int i, j, l, u;
    for (i = 0; i <= m - 1; i++)
    {
        for (j = 0; j <= k - 1; j++)
        {
            u = i * k + j; 
            c[u] = 0.0;
            for (l = 0; l <= n - 1; l++)
                c[u] = c[u] + a[i * n + l] * b[l * k + j];
        }
    }
    return;
}







void xyz2rtn(double *x, double *v, double *xyz, double *rtn)
{ 
    double scal_x, scal_v, vr[3], vn[3], vt[3];

    scal_x=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    scal_v=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

// c...unit vector in R direction
    vr[0]=x[0]/scal_x;
    vr[1]=x[1]/scal_x;
    vr[2]=x[2]/scal_x;
 
// c...unit direction in N direction
    vn[0]=(vr[1]*v[2]-vr[2]*v[1])/scal_v;
    vn[1]=(vr[2]*v[0]-vr[0]*v[2])/scal_v;
    vn[2]=(vr[0]*v[1]-vr[1]*v[0])/scal_v;
        
// c...unit direction in T direction
    vt[0]=(vn[1]*vr[2]-vn[2]*vr[1]);
    vt[1]=(vn[2]*vr[0]-vn[0]*vr[2]);
    vt[2]=(vn[0]*vr[1]-vn[1]*vr[0]);


// drtn(i,4)=dsqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2])
    rtn[0]=xyz[0]*vr[0]+xyz[1]*vr[1]+xyz[2]*vr[2];
    rtn[1]=xyz[0]*vt[0]+xyz[1]*vt[1]+xyz[2]*vt[2];
    rtn[2]=xyz[0]*vn[0]+xyz[1]*vn[1]+xyz[2]*vn[2];

    return;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* lagrange interpolation order = 6, 2*order points
* @param1: description of param1
* @param2: description of param2
* todo
    order = input parameter    
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double lagrange (double *y, int dim_y, int dim_x, double t, double *z)
{ 
    int i, j, k, m, dim, order = 8;
    double s;
    
    i = 0;
    while ((y[i * dim_x] < t) && (i < dim_y)) 
        i = i + 1;
    k = i - order;
    if (k < 0) 
        k = 0;
    m = i + order - 1;
    if (m > dim_y - 1) 
        m = dim_y - 1;

    for (dim = 0; dim < dim_x - 1; dim++)
    {
        z[dim] = 0;
    }

    for (i = k; i <= m; i++)
    { 
        s = 1.0; 
        for (j = k; j <= m; j++)
        {
            if (j != i) 
            {
                s = s * (t - y[j * dim_x]) / (y[i * dim_x] - y[j * dim_x]);
            }
        }        
        for (dim = 0; dim < dim_x - 1; dim++)
        {
            z[dim] = z[dim] + s * y[i * dim_x + dim + 1];
        }
    }
    return 0;
}











double obs_alt (double jd, double utc, double *obs, int part, double *bmat)
{
    int n, lps, i;
    double r, h, ref, *eph, *dxdp, *dodpo, *dodpd, *dodpp, xc2[6], dxdx0[36], dodx[6], 
        dodx0[6], tt, tjd[2], xsc[6], dx[3];

    ref = RCT; //should be topography height

    lps = getlps (JD0 + utc/86400.0);
    tt = utc + (lps + 32.184);

//    tjd[0] = JD0;    tjd[1] = tt / 86400.0;
//    get_ephemeris (tjd, 2, CT, xsc);
 
    if (part == 0)
    {
        lagrange (OR_EPH, DIM_OR, 7, tt, xc2);
    }
    if (part == 1)
    {
        eph = (double *) calloc (42 + 6 * MDYN, sizeof(double));
        lagrange (OR_EPH, DIM_OR, 42 + 6 * MDYN + 1, tt, eph);
        for (n = 0; n < 6; n++)
            xc2[n] = eph[n];
        for (n = 0; n < 36; n++)
            dxdx0[n] = eph[n + 6];
        if (MDYN > 0)
        {
            dxdp = (double *) calloc (6 * MDYN, sizeof(double));
            dodpd = (double *) calloc (MDYN, sizeof(double));
            for (n = 0; n < 6 * MDYN; n++)
                dxdp[n] = eph[n + 42];
        }
        if (MOBS > 0)
            dodpo = (double *) calloc (MOBS, sizeof(double));
        free (eph);
    }

    h = modvect(xc2) - RCT;

    for (n = 0; n < 3; n++)
        dx[n] = xc2[n];
    r = modvect(dx);

    *obs = r - ref + BASB + BAST * utc;

    if (part == 0)
        return h;

    for (n = 0; n < 3; n++)
    {
        dodx[n] = (xc2[n])/r;
        dodx[n + 3] = 0;
    }

    brmul (dodx, dxdx0, 1, 6, 6, dodx0);
    for (n = 0; n < 6; n++)
        bmat[n] = dodx0[n];
    
    if (MEST == 0) 
        return h;

    i = 0;
    if (MOBS > 0)
    {     
        dodpo[i] = 1;
        i++;
    }
    if (MOBS > 1)
    {
        dodpo[i] = utc;
        i++;
    }

    if (MDYN > 0)
    {
        brmul (dodx, dxdp, 1, 6, MDYN, dodpd);
    }
    if (MOBS > 0)
    for (n = 0; n < MOBS; n++)
        bmat[6 + n] = dodpo[n];
    if (MDYN > 0)
    for (n = 0; n < MDYN; n++)
        bmat[6 + MOBS + n] = dodpd[n];
    
    if (MDYN > 0)
    {
        free (dxdp);
        free (dodpd);
    }
    if (MOBS > 0)
        free (dodpo);
    return h;


}






double obs_vel (double jd, double utc, double *obs, int part, double *bmat)
{
    int n, lps, i;
    double v, h, *eph, *dxdp, *dodpo, *dodpd, *dodpp, xc2[6], dxdx0[36], dodx[6], 
        dodx0[6], tt, tjd[2], xsc[6], dv[3];

    lps = getlps (JD0 + utc/86400.0);
    tt = utc + (lps + 32.184);

    tjd[0] = JD0;    tjd[1] = tt / 86400.0;
    get_ephemeris (tjd, 2, CT, xsc);
 
    if (part == 0)
    {
        lagrange (OR_EPH, DIM_OR, 7, tt, xc2);
    }
    if (part == 1)
    {
        eph = (double *) calloc (42 + 6 * MDYN, sizeof(double));
        lagrange (OR_EPH, DIM_OR, 42 + 6 * MDYN + 1, tt, eph);
        for (n = 0; n < 6; n++)
            xc2[n] = eph[n];
        for (n = 0; n < 36; n++)
            dxdx0[n] = eph[n + 6];
        if (MDYN > 0)
        {
            dxdp = (double *) calloc (6 * MDYN, sizeof(double));
            dodpd = (double *) calloc (MDYN, sizeof(double));
            for (n = 0; n < 6 * MDYN; n++)
                dxdp[n] = eph[n + 42];
        }
        if (MOBS > 0)
            dodpo = (double *) calloc (MOBS, sizeof(double));
        free (eph);
    }

    h = modvect(xc2) - RCT;
        
    for (n = 0; n < 3; n++)
        dv[n] = xc2[n + 3] - xsc[n + 3];
    v = modvect(dv);

    *obs = v + BASB + BAST * utc;

    if (part == 0)
        return h;

    for (n = 0; n < 3; n++)
    {
        dodx[n] = 0;
        dodx[n + 3] = (xc2[n + 3] - xsc[n + 3])/v;
    }

    brmul (dodx, dxdx0, 1, 6, 6, dodx0);
    for (n = 0; n < 6; n++)
        bmat[n] = dodx0[n];
    
    if (MEST == 0) 
        return h;

    i = 0;
    if (MOBS > 0)
    {     
        dodpo[i] = 1;
        i++;
    }
    if (MOBS > 1)
    {
        dodpo[i] = utc;
        i++;
    }

    if (MDYN > 0)
    {
        brmul (dodx, dxdp, 1, 6, MDYN, dodpd);
    }
    if (MOBS > 0)
    for (n = 0; n < MOBS; n++)
        bmat[6 + n] = dodpo[n];
    if (MDYN > 0)
    for (n = 0; n < MDYN; n++)
        bmat[6 + MOBS + n] = dodpd[n];
    
    if (MDYN > 0)
    {
        free (dxdp);
        free (dodpd);
    }
    if (MOBS > 0)
        free (dodpo);
    return h;


}









double obs_dsn (double jd, double utc, double *obs, int part, double *bmat)
{
    int n, lps, i;
    double r, h, *eph, *dxdp, *dodpo, *dodpd, *dodpp, xc2[6], dxdx0[36], dodx[6], 
        dodx0[6], tt, tjd[2], xsc[6], dx[3];

    lps = getlps (JD0 + utc/86400.0);
    tt = utc + (lps + 32.184);

    tjd[0] = JD0;    tjd[1] = tt / 86400.0;
    get_ephemeris (tjd, 2, CT, xsc);
 
    if (part == 0)
    {
        lagrange (OR_EPH, DIM_OR, 7, tt, xc2);
    }
    if (part == 1)
    {
        eph = (double *) calloc (42 + 6 * MDYN, sizeof(double));
        lagrange (OR_EPH, DIM_OR, 42 + 6 * MDYN + 1, tt, eph);
        for (n = 0; n < 6; n++)
            xc2[n] = eph[n];
        for (n = 0; n < 36; n++)
            dxdx0[n] = eph[n + 6];
        if (MDYN > 0)
        {
            dxdp = (double *) calloc (6 * MDYN, sizeof(double));
            dodpd = (double *) calloc (MDYN, sizeof(double));
            for (n = 0; n < 6 * MDYN; n++)
                dxdp[n] = eph[n + 42];
        }
        if (MOBS > 0)
            dodpo = (double *) calloc (MOBS, sizeof(double));
        free (eph);
    }

    h = modvect(xc2) - RCT;
        
    for (n = 0; n < 3; n++)
        dx[n] = xc2[n] - xsc[n];
    r = modvect(dx);

    *obs = r + BASB + BAST * utc;

    if (part == 0)
        return h;

    for (n = 0; n < 3; n++)
    {
        dodx[n] = (xc2[n] - xsc[n])/r;
        dodx[n + 3] = 0;
    }

    brmul (dodx, dxdx0, 1, 6, 6, dodx0);
    for (n = 0; n < 6; n++)
        bmat[n] = dodx0[n];
    
    if (MEST == 0) 
        return h;

    i = 0;
    if (MOBS > 0)
    {     
        dodpo[i] = 1;
        i++;
    }
    if (MOBS > 1)
    {
        dodpo[i] = utc;
        i++;
    }

    if (MDYN > 0)
    {
        brmul (dodx, dxdp, 1, 6, MDYN, dodpd);
    }
    if (MOBS > 0)
    for (n = 0; n < MOBS; n++)
        bmat[6 + n] = dodpo[n];
    if (MDYN > 0)
    for (n = 0; n < MDYN; n++)
        bmat[6 + MOBS + n] = dodpd[n];
    
    if (MDYN > 0)
    {
        free (dxdp);
        free (dodpd);
    }
    if (MOBS > 0)
        free (dodpo);
    return h;


}







void getsolvefor ()
{
       
//    MOBS = 0;   //2; 
//    MSRP = 0;   //2; 
//    MTK2 = 0;   //1; 
//    MGCS = 1;   //6;

    MDYN = MSRP + MTK2 + MGCS;                // dim of sensitivity matrix
    MSOL = MOBS + MDYN + 6;                       // dim of regress matrix
    MEST = MOBS + MDYN;  //= MOBS + MDYN - 6

    MSTA = 42 + MDYN * 6;

/*
    if (MGCS > 0)    
    {
        CSinfo = (CSStruct *) calloc ( MGCS, sizeof(CSStruct));

        CSinfo[0].n = 2; CSinfo[0].m = 0; CSinfo[0].cs = 0;
//        CSinfo[1].n = 3; CSinfo[1].m = 0; CSinfo[1].cs = 0; 
    }
*/
    return;

}








void initsolvefor (double *xsm, double *x)
{
    int i, k, n, m, ind, l, ic, is, label;

    for (k = 0; k < 6; k ++)
    {
        x[k] = xsm[k];
    }
    i = 6;
    if (MOBS > 0)
    {
        x[i] = BASB;
        i++;
    }
    if (MOBS > 1)
    {
        x[i] = BAST;
        i++;
    }
    if (MSRP > 0)
    {
        x[i] = SRPB;
        i++;
    }
    if (MSRP > 1)
    {
        x[i] = SRPT ;
        i++;
    }
    if (MTK2 > 0)
    {
        x[i] = K2;
        i++;
    }

    if (MGCS > 0)
    for (k = 0; k < MGCS; k ++)
    {
        n = CSinfo[k].n; m = CSinfo[k].m; label = CSinfo[k].cs;

        if (m == 0)
        {
            x[i] = COEFG[n] + CSinfo[k].initv;
            COEFG[n] = x[i];
        }
        else
        {
            l = NMAX - m + 1;
            ind = NMAX + 1 + (2 * NMAX - m + 2) * (m - 1);
            ic = ind + n - m;
            is = ind + n - m + l;
            if (label == 1)
            {
                x[i] = COEFG[ic] + CSinfo[k].initv;
                COEFG[ic] = x[i];
            }

            if (label == -1)
            {
                x[i] = COEFG[is] + CSinfo[k].initv;
                COEFG[is] = x[i];
            }
        }
        i++;
    }

    return;

}






void updsolvefor (double *x)
{
    int i, k, n, m, ind, l, ic, is, label;

    i = 6;
    if (MOBS > 0)
    {
        BASB = x[i];
        i++;
    }
    if (MOBS > 1)
    {
        BAST = x[i];
        i++;
    }
    if (MSRP > 0)
    {
        SRPB = x[i];
        i++;
    }
    if (MSRP > 1)
    {
        SRPT = x[i];
        i++;
    }
    if (MTK2 > 0)
    {
        K2 = x[i];
        i++;
    }

    if (MGCS > 0)
    for (k = 0; k < MGCS; k ++)
    {
        n = CSinfo[k].n; m = CSinfo[k].m; label = CSinfo[k].cs;

        if (m == 0)
        {
            COEFG[n] = x[i];
        }
        else
        {
            l = NMAX - m + 1;
            ind = NMAX + 1 + (2 * NMAX - m + 2) * (m - 1);
            ic = ind + n - m;
            is = ind + n - m + l;
            if (label == 1)
            {
                COEFG[ic] = x[i];
            }

            if (label == -1)
            {
                COEFG[is] = x[i];
            }
        }
        i++;
    }



}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* simula_phase - simulate total phase count observable
* @param1: description of param1
* @param2: description of param2
* todo: 
        1 one-way doppler deltat accumlated error
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double simula_phase (double utc3, double utc0, double *station3, 
                     short int uplink, double *station1, short int genrel,
                     double *calculable, double *azimuth, double *elevation, 
                     short int part, double *bmat)
{
    double txice[7], txics[7], *bmats, *bmate, deltat;
    int n;
    
    real128 lts2[3], lte2[3];

    bmats = (double *) calloc ( SLOVEFOR, sizeof(double));
    bmate = (double *) calloc ( SLOVEFOR, sizeof(double));


    ltsolution (utc0, station3, uplink, station1, genrel, lts2, 
        azimuth, elevation, part, bmats, txics);
    ltsolution (utc3, station3, uplink, station1, genrel, lte2, 
        azimuth, elevation, part, bmate, txice);
  
    *calculable = (lte2[2] - lts2[2]) * C;

    if (uplink == 0)    //one-way doppler deltat, time correction
    {
        delta_tdb (txice, txics, &deltat);
        *calculable = *calculable + deltat * (txice[0] - txics[0]) * C;
    }

    if (part == 1)
    {
        for (n = 0; n < 6 + DYNPAR; n++)
        {
            bmat[n] = bmate[n] - bmats[n];
        }
        bmat[8] = 1;
        bmat[9] = utc3;
    }
    *calculable = *calculable + BIAS + DBIA * utc3;
    
    free (bmats); 
    free (bmate);

    return 0;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* simula_dople - simulate doppler observable
* @param1: description of param1
* @param2: description of param2
* todo: 
        1 one-way doppler deltat accumlated error
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double simula_dople (double utc3, double tc, double *station3, 
                     short int uplink, double *station1, short int genrel,
                     double *calculable, double *azimuth, double *elevation, 
                     short int part, double *bmat)
{
    double txice[7], txics[7], *bmats, *bmate, deltat;
    double dop_old, dop_new;
    int n;
    real128 lts2[3], lte2[3], dlt;

    bmats = (double *) calloc ( SLOVEFOR, sizeof(double));
    bmate = (double *) calloc ( SLOVEFOR, sizeof(double));


    ltsolution (utc3 + tc / 2, station3, uplink, station1, genrel, lte2, 
        azimuth, elevation, part, bmate, txice);
    ltsolution (utc3 - tc / 2, station3, uplink, station1, genrel, lts2, 
        azimuth, elevation, part, bmats, txics);

    dop_old = (double) ((lte2[1] - lts2[1]) / (real128)tc);

    dlt = lte2[2] - lts2[2];
    dop_new = (double) (dlt / (real128)tc * (real128)C);
    *calculable = dop_new;

    if (uplink == 0)     //one-way doppler deltat, proper time correction
    {
        delta_tdb (txice, txics, &deltat);
        *calculable = *calculable + deltat * C;
    }

    if (part == 1)
    {
        for (n = 0; n < 6 + DYNPAR; n++)
        {
            bmat[n] = (bmate[n] - bmats[n]) / tc;
        }
        bmat[8] = 1;
        bmat[9] = utc3;
    }
    *calculable = *calculable + BIAS + DBIA * utc3;
    free (bmats); 
    free (bmate);
    return 0;
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* delta_tdb - one-way doppler deltat, proper time correction
* txice: [0]:  satellite TDB time(s), [1]~[7]satellite coordinates(AU, day)
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double delta_tdb (double *txice, double *txics, double *deltat)
{
    double ie, is, t, ied, isd, tjde[2], tjds[2], xcbe[6], xcbs[6],
        xibe[6], xibs[6];
    short int ssbary, n;
    
    ssbary = 11;
    t = txice[0] - txics[0];
    tjde[0] = JD0;
    tjde[1] = txice[0] / 86400.0;
    tjds[0] = JD0;
    tjds[1] = txics[0] / 86400.0;

    planet_ephemeris (tjde, CENTER, ssbary, &xcbe[0], &xcbe[3]);
    for (n = 0; n < 6; n++)
        xibe[n] = txice[n + 1] + xcbe[n];
    planet_ephemeris (tjds, CENTER, ssbary, &xcbs[0], &xcbs[3]);
    for (n = 0; n < 6; n++)
        xibs[n] = txics[n + 1] + xcbe[n];

    delta_iid (tjde, xibe, &ie, &ied);
    delta_iid (tjds, xibs, &is, &isd);

    *deltat = 1.0 / 2.0 * (ie + is) * t - 1.0 / 12.0 * (ied - isd) * t * t;
    return 0;
}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* delta_tdb -  one-way light-time proper time correction
* txice: [0]:  satellite TDB time(s), [1]~[7]satellite coordinates(AU, day)
* @param2: description of param2
* todo: 
        1 Uobl, time correction due to non-spherical potential 
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double delta_iid (double *jd, double *xi, double *ii, double *id)
{
    double uu, vv, L, ud, vd,
    rdd[3], xj[11][6], xij[11][6], rij[11], rijd[11]; 
    short int ssbary, j, n;
    
    ssbary = 11;
    L = 1.550520e-8;

    uu = 0;
    ud = 0;
    for (n = 0; n < 3; n ++)
        rdd[n] = 0;
    for (j = 0; j <= 10; j++)
    {
        if (PERB[j] == 0)
            continue;
        planet_ephemeris (jd, j, ssbary, &xj[j][0], &xj[j][3]);
        for (n = 0; n < 6; n++)
        {
            xij[j][n] = xi[n] - xj[j][n];
        }
        rij[j] = sqrt (xij[j][0] * xij[j][0] 
            + xij[j][1] * xij[j][1] + xij[j][2] * xij[j][2]);
        rijd[j] = (xij[j][0] * xij[j][3] + xij[j][1] * xij[j][4] 
            + xij[j][2] * xij[j][5]) / rij[j];
        uu = uu + GMDE[j] / rij[j];
        ud = ud - GMDE[j] / rij[j] / rij[j] * rijd[j];
        for (n = 0; n < 3; n++)
            rdd[n] = rdd[n] 
            - GMDE[j] / (rij[j] * rij[j] * rij[j]) * xij[j][n];
    }

    vv = xi[3] * xi[3] + xi[4] * xi[4] + xi[5] * xi[5];
    vd = 2.0 * (xi[3] * rdd[0] + xi[4] * rdd[1] + xi[5] * rdd[2]);

    *ii = (uu + vv / 2.0) / C_AUDAY / C_AUDAY - L;
    *id = (ud + vd / 2.0) / C_AUDAY / C_AUDAY / 86400.0;

    return 0;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* simula_range - simulate doppler observable
* @param1: description of param1
* @param2: description of param2
* todo: 

* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double simula_range (double utc3, double *station3, short int uplink, 
                     double *station1, short int genrel, 
                     double *calculable, double *azimuth, double *elevation, 
                     short int part, double *bmat)
{
    double txic[7];
    real128 lt[3];

    ltsolution (utc3, station3, uplink, station1, genrel, lt, 
        azimuth, elevation, part, bmat, txic);

    if (uplink == 0)    //one-way range: only work for near-earth satellite
    {
        *calculable = (double) lt[0];
    }
    if (uplink == 1)    //two/three-way range: work for deep space
    {
        *calculable = (double) (lt[2] * (real128)C);
    }

    if (part == 1)
    {
        bmat[8] = 1;    //DBIAS
        bmat[9] = utc3;
    }
    *calculable = *calculable + BIAS + DBIA * utc3;
    return 0;
}






/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*

* ltsolution - light time solution
* @param:
    utc3        : unit: day; 
    station3    : receive station;
    uplink      : no uplink == 0; yes uplink == 1 
    station1    : transmit station;		
    genrel      : no general relativity correction == 0; yes == 1
    calculable  : ;
    azimuth     : ;
    elevation   : ;
    partial     : no partial == 0 ; yes partial == 1 
    bmat        : (partial == 0: satellite coordinates(6), partial == 1: partial)
    txic        : satellite coordinates(t2)

* 
* version: 20 Aug 2010

*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double ltsolution (double utc_3, double *station3, short int uplink,  
                   double *station1, short int genrel, real128 *lt, 
                   double *azimuth, double *elevation, short int part,
                   double *bmat, double *txic)
{
    double re3fi[3], re3[3], re1fi[3], re1[3], vec[3], xe3[6], re3n[3], 
        xe1[6], re1n[3], xc2[6], 
        secdiff, 
        ra, dec, zd, az, secdiff3, secdiff1,
        utc_1, ut1_3, ut1_1, tt_3, tt_1, tdb_3, tdb_2, tdb_1, 
        ut1_utc, xp, yp, xp3, yp3, dx, dy, delta_t3, t, elong, u, v, 
        dxdx0[36], dxdp[6 * DYNPAR], dodx[6], dodp[DYNPAR], 
        dodx0[6], dodpp[DYNPAR], eph[42 + 6 * DYNPAR], 
        te[9], llh3[3], llh1[3];
    real128 tao231, tao232, tao121, tao122, taoerr, r23,  r12, xb3[6], xb2[6], xb1[6];
    int n, flag, dim_par;
	
    taoerr   = 1.0e-12L; //1nanosec;
//    taoerr   = 1.0e-8L; //1nanosec;
    dim_par = 42 + 6 * DYNPAR + 1;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*--light time iteration 2 -> 3 --*/

    re3fi[0] = station3[0];
    re3fi[1] = station3[1];
    re3fi[2] = station3[2];
    xyz2llh(re3fi, llh3);

    elong = llh3[1] * DEG2RAD;
    u = sqrt (re3fi[0] * re3fi[0] + re3fi[1] * re3fi[1]) / 1000.0;
    v = re3fi[2] / 1000.0;

/*--time scales transformation --*/
    geteop (utc_3, &xp3, &yp3, &ut1_utc, &dx, &dy);	
    delta_t3 = 32.184 + LEAPSECS - ut1_utc;
    ut1_3 = utc_3 + ut1_utc;
    tt_3 = utc_3 + (LEAPSECS + 32.184); 	
    secdiff3 = iauDtdb (JD0, tt_3 / 86400.0, ut1_3 / 86400.0, elong, u, v);
    tdb_3 = tt_3 + secdiff3;

/*--station coordinate interpolation--*/    
    lagrange (TE_EPH, DIM_TE, 10, utc_3, te);        
    brmul (te, re3fi, 3, 3, 1, re3);
    lagrange (TE_EPH, DIM_TE, 10, utc_3 + 1.0, te);        
    brmul (te, re3fi, 3, 3, 1, re3n);
    for (n = 0; n < 3; n++)
    {
        xe3[n] = re3[n] / AU;
        xe3[n + 3] = (re3n[n] - re3[n]) / AU * 86400;
    }

/*--satellite coordinate interpolation--*/    
    if (part == 0)
    {
        lagrange (OR_EPH, DIM_OR, 7, tdb_3, xc2);
    }
    if (part == 1)
    {
        lagrange (OR_EPH, DIM_OR, dim_par, tdb_3, eph);
        for (n = 0; n < 6; n++)
            xc2[n] = eph[n + 36];
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*--iteration--*/    
    r23 = lt_form (tdb_3, tdb_3, xe3, xc2, genrel, xb3, xb2);
    tao231 = (real128)r23 * (real128)AU_SEC;
    tdb_2 = tdb_3 - (double) tao231;

    flag = -1;
    do
    {
        flag++;
        tao232 = tao231;
        if (part == 0)
        {
            lagrange (OR_EPH, DIM_OR, 7, tdb_2, xc2);
        }
        if (part == 1)
        {
            lagrange (OR_EPH, DIM_OR, dim_par, tdb_2, eph);
            for (n = 0; n < 6; n++)
                xc2[n] = eph[n + 36];
        }
        
        r23 = lt_form (tdb_3, tdb_2, xe3, xc2, genrel, xb3, xb2);
        tao231 = (real128)r23 * (real128)AU_SEC;
        tdb_2 = tdb_3 - (double) tao231;
    }while (fabsl (tao232-tao231) > taoerr);



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*--light time iteration 1 -> 2--*/    


    if (uplink == 1)
    {
        re1fi[0] = station1[0];
        re1fi[1] = station1[1];
        re1fi[2] = station1[2];
        xyz2llh(re1fi, llh1);

        elong = llh1[1] * DEG2RAD;
        u = sqrt (re1fi[0] * re1fi[0] + re1fi[1] * re1fi[1]) / 1000.0;
        v = re1fi[2] / 1000.0;

    /*--time scales transformation --*/
        tdb_1 = tdb_2; //unit: s
        tdb2tt (JD0 + tdb_1 / 86400.0, &t, &secdiff);
        tt_1 = tdb_1 - secdiff;
        utc_1 = tt_1 - (LEAPSECS + 32.184); 	
        geteop (utc_1, &xp, &yp, &ut1_utc, &dx, &dy);	
        ut1_1 = utc_1 + ut1_utc;	        
        secdiff1 = iauDtdb (JD0, tdb_1 / 86400.0, ut1_1 / 86400.0, 
            elong, u, v);
        tt_1 = tdb_1 - secdiff1;
        utc_1 = tt_1 - (LEAPSECS + 32.184); 	

    /*--station coordinate interpolation--*/    
        lagrange (TE_EPH, DIM_TE, 10, utc_1, te);        
        brmul (te, re1fi, 3, 3, 1, re1);
        lagrange (TE_EPH, DIM_TE, 10, utc_1 + 1.0, te);        
        brmul (te, re1fi, 3, 3, 1, re1n);
        for (n = 0; n < 3; n++)
        {
            xe1[n] = re1[n] / AU;
            xe1[n + 3] = (re1n[n] - re1[n]) / AU * 86400;
        }

        r12 = lt_form (tdb_1, tdb_2, xe1, xc2, genrel, xb1, xb2);
        tao121 = (real128)r12 * (real128)AU_SEC;
        tdb_1 = tdb_2 - (double) tao121;
     
    /*--iteration--*/    
        flag = -1;
        do
        {
            flag++;
            tao122 = tao121;

        /*--time scales transformation --*/
            tdb2tt (JD0 + tdb_1 / 86400.0, &t,&secdiff);
            tt_1 = tdb_1 - secdiff;
            utc_1 = tt_1 - (LEAPSECS + 32.184); 
            geteop (utc_1, &xp, &yp, &ut1_utc, &dx, &dy);	

            ut1_1 = utc_1 + ut1_utc;
            secdiff1 = iauDtdb (JD0, tdb_1 / 86400.0, ut1_1 / 86400.0, 
                elong, u, v);
            tt_1 = tdb_1 - secdiff1;
            utc_1 = tt_1 - (LEAPSECS + 32.184); 	

        /*--station coordinate interpolation--*/    
            lagrange (TE_EPH, DIM_TE, 10, utc_1, te);        
            brmul (te, re1fi, 3, 3, 1, re1);
            lagrange (TE_EPH, DIM_TE, 10, utc_1 + 1.0, te);        
            brmul (te, re1fi, 3, 3, 1, re1n);
            for (n = 0; n < 3; n++)
            {
                xe1[n] = re1[n] / AU;
                xe1[n + 3] = (re1n[n] - re1[n]) / AU * 86400;
            }
            
            r12 = lt_form (tdb_1, tdb_2, xe1, xc2, genrel, xb1, xb2);
            tao121 = (real128)r12 * (real128)AU_SEC;
            tdb_1 = tdb_2 - (double) tao121;
        }while (fabsl (tao122-tao121) > taoerr);
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*-- partial time: tdb2 --*/    
    if (part == 0)
    {
        lagrange (OR_EPH, DIM_OR, 7, tdb_2, bmat);

    }
    if (part == 1)
    {
        lagrange (OR_EPH, DIM_OR, dim_par, tdb_2, eph);

        for (n = 0; n < 36; n++)
            dxdx0[n] = eph[n];
        for (n = 0; n < 6 * DYNPAR; n++)
            dxdp[n] = eph[n + 42];
        lt_part (xb3, xb2, xb1, uplink, dodx, dodp);
        brmul (dodx, dxdx0, 1, 6, 6, dodx0);
        brmul (dodx, dxdp, 1, 6, DYNPAR, dodpp);
        for (n = 0; n < DYNPAR; n++)
            dodp[n] = dodpp[n] + dodp[n];

        for (n = 0; n < 3; n++)
        {
            bmat[n] = dodx0[n];
            bmat[n + 3] = dodx0[n + 3] * 86400.0;
        }
        bmat[6] = dodp[0] * AU; // l/c: au, 
        bmat[7] = dodp[1] * AU * 86400.0; // l/(c*d-1): au*d, 
    }

    txic[0] = tdb_2;
    for (n = 1; n < 7; n++)	
    {
        txic[n] = xc2[n - 1];
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*--calculate light time observable --*/
    if (uplink == 0)
    {
        lt[0] = (real128)r23 * (real128)AU;			
        lt[1] = ((real128)utc_3 - (real128)tdb_2) * (real128)C;
        lt[2] = (real128)tao231 - (real128)secdiff3 - (real128)(LEAPSECS + 32.184);
    }

    if (uplink == 1)
    {
        lt[0] = ((real128)r12 + (real128)r23) * (real128)AU;			
        lt[1] = ((real128)utc_3 - (real128)utc_1) * (real128)C;			
        lt[2] = (real128)tao231 + (real128)tao121 + (real128)secdiff1 - (real128)secdiff3;			
//        lt[2] = tao231 + tao121 + secdiff1 - secdiff3;			
    }

    for (n = 0; n < 3; n++)
        vec[n] = xb2[n] - xb3[n];
    vector2radec (vec, &ra,&dec);

    azelev (ut1_3 / 86400.0 + JD0, delta_t3, ACCURACY,
               xp3, yp3, llh3, ra, dec, &zd, &az);

    *azimuth = az;
    *elevation = 90.0 - zd;

    return 0;
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* lt_part - partial of light time 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double lt_part (real128 *xb3, real128 *xb2, real128 *xb1, int uplink, 
                 double *dodx, double *dodp)
{
    double r23, r12, p23, p12, pt2, pt1, rp12, rpt1, dt2dx[3], dt1dx[3];
    int n;
    
    r23 = sqrt ((xb2[0] - xb3[0]) * (xb2[0] - xb3[0]) 
        + (xb2[1] - xb3[1]) * (xb2[1] - xb3[1])
        + (xb2[2] - xb3[2]) * (xb2[2] - xb3[2]));

    p23 = ((xb3[0] - xb2[0]) * xb2[3] + (xb3[1] - xb2[1]) * xb2[4]
        + (xb3[2] - xb2[2]) * xb2[5]) / r23;

    pt2 = (1 - p23 / C_AUDAY);

    dt2dx[0] = (xb3[0] - xb2[0]) / r23 / pt2;
    dt2dx[1] = (xb3[1] - xb2[1]) / r23 / pt2;
    dt2dx[2] = (xb3[2] - xb2[2]) / r23 / pt2;

    if (uplink == 0)
    {
        dodx[0] = - dt2dx[0];
        dodx[1] = - dt2dx[1];
        dodx[2] = - dt2dx[2];
    }
    if (uplink == 1)
    {
        r12 = sqrt ((xb2[0] - xb1[0]) * (xb2[0] - xb1[0]) 
            + (xb2[1] - xb1[1]) * (xb2[1] - xb1[1])
            + (xb2[2] - xb1[2]) * (xb2[2] - xb1[2]));

        p12 = ((xb2[0] - xb1[0]) * xb1[3] + (xb2[1] - xb1[1]) * xb1[4]
            + (xb2[2] - xb1[2]) * xb1[5]) / r12;
        
        pt1 = (1 - p12 / C_AUDAY);

        rp12 = ((xb2[0] - xb1[0]) * (xb2[3] - xb1[3]) 
              + (xb2[1] - xb1[1]) * (xb2[4] - xb1[4]) 
              + (xb2[2] - xb1[2]) * (xb2[5] - xb1[5])) / r12;

        rpt1 = (1 - (rp12 + p12) / C_AUDAY);


        for (n = 0; n < 3; n++)
        {
            dt1dx[n] = (dt2dx[n] * rpt1 - (xb2[n] - xb1[n]) / r12) / pt1;
        }

        for (n = 0; n < 3; n++)
        {
            dodx[n] = - dt1dx[n];
        }

    }

    dodx[3] = 0; 
    dodx[4] = 0; 
    dodx[5] = 0;
    
    dodp[0] = 0; //
    dodp[1] = 0; //
    
    return 0;
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* lt_form - calculate the light time equation
* @param:
    tdb3, re3[3]    : station time, coordinates (unit:AU)
    tdb2, rp2[3]    : satellite time, coordinates (AU)
    genrel          : 
    *rs3            : output: station coordinates to SSB (AU)
    *rs2            : output: satellite coordinates to SSB (AU)
    return          : light time solution (AU)
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
real128 lt_form (double tdb3, double tdb2, double *re3, double *rp2, 
                int genrel, real128 *rs3, real128 *rs2)
{
    double gamma, tjd2[2], tjd3[2], re[3], ve[3], 
        rp[3], vp[3], xe[6], xp[6];
    real128 rlight, rgen, ri, rj, r12, rse[3], rsp[3], r23[3], rlt[3];
    short int earth = 2, sun = 10, j;
    int n;
    gamma = 1;
    tjd2[0] = JD0;
    tjd3[0] = JD0;
    tjd2[1] = tdb2 / 86400.0;
    tjd3[1] = tdb3 / 86400.0;

    planet_ephemeris (tjd3, earth, sun, &xe[0], &xe[3]);
    for (n = 0; n < 6; n++)
        rs3[n] = (real128)re3[n] + (real128)xe[n];
    ri = sqrtl (rs3[0] * rs3[0] + rs3[1] * rs3[1] + rs3[2] * rs3[2]);

    planet_ephemeris (tjd2, CENTER, sun, &xp[0], &xp[3]);
    for (n = 0; n < 6; n++)
        rs2[n] = (real128)rp2[n] + (real128)xp[n];
    rj = sqrtl (rs2[0] * rs2[0] + rs2[1] * rs2[1] + rs2[2] * rs2[2]);


    for (n = 0; n < 3; n++)
        r23[n] = (real128)xe[n] - (real128)xp[n];
    for (n = 0; n < 3; n++)
        rlt[n] = (real128)r23[n] + (real128)re3[n] - (real128)rp2[n];


//    rlight = sqrtl (((real128)rs3[0] - (real128)rs2[0]) * ((real128)rs3[0] - (real128)rs2[0]) 
//        + ((real128)rs3[1] - (real128)rs2[1]) * ((real128)rs3[1] - (real128)rs2[1])
//        + ((real128)rs3[2] - (real128)rs2[2]) * ((real128)rs3[2] - (real128)rs2[2]));

     rlight = sqrtl( rlt[0] * rlt[0] +  rlt[1] * rlt[1] +  rlt[2] * rlt[2]);
     rlight = sqrtl ((rs3[0] - rs2[0]) * (rs3[0] - rs2[0]) 
        + (rs3[1] - rs2[1]) * (rs3[1] - rs2[1])
        + (rs3[2] - rs2[2]) * (rs3[2] - rs2[2]));


    if (genrel == 1)
    {
        rgen = (1L + (real128)gamma) * (real128)GMDE[10] / (real128)C_AUDAY / (real128)C_AUDAY 
            * logl ((ri + rj + rlight 
            + (1L + (real128)gamma) * (real128)GMDE[10] / (real128)C_AUDAY / (real128)C_AUDAY) 
            / (ri + rj - rlight 
            + (1L + (real128)gamma) * (real128)GMDE[10] / (real128)C_AUDAY / (real128)C_AUDAY));
        rlight = rlight + rgen;
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
        for (j = 0; j <= 9; j++)
        {
//			if (j ==9) continue;		
            planet_ephemeris (tjd3, earth, j, re, ve);
            for (n = 0; n < 3; n++)
                rse[n] = (real128)re3[n] + (real128)re[n];
            ri = sqrtl (rse[0] * rse[0] + rse[1] * rse[1] + rse[2] * rse[2]);

            planet_ephemeris (tjd2, CENTER, j, rp, vp);
            for (n = 0; n < 3; n++)
                rsp[n] = (real128)rp2[n] + (real128)rp[n];
            rj = sqrtl (rsp[0] * rsp[0] + rsp[1] * rsp[1] + rsp[2] * rsp[2]);
	
            r12 = sqrtl((rse[0] - rsp[0]) * (rse[0] - rsp[0]) 
                + (rse[1] - rsp[1]) * (rse[1] - rsp[1])
                + (rse[2] - rsp[2]) * (rse[2] - rsp[2]));
            rgen = (1L + (real128)gamma) * (real128)GMDE[j] / (real128)C_AUDAY / (real128)C_AUDAY 
                * logl ((ri + rj + r12 ) / (ri + rj - r12));
            rlight = rlight + rgen;
        }
    }
    return rlight;
}





void azelev (double jd_ut1, double delta_t, short int accuracy,
              double x, double y, double *llh, double ra,
              double dec, double *zd, double *az)

{

   double sinlat, coslat, sinlon, coslon, sindc, cosdc, sinra, cosra,
      uze[3], une[3], uwe[3], uz[3], un[3], uw[3], p[3], pz, pn, pw,
      proj;

/*
   Preliminaries.
*/

   sinlat = sin (llh[0] * DEG2RAD);
   coslat = cos (llh[0] * DEG2RAD);
   sinlon = sin (llh[1] * DEG2RAD);
   coslon = cos (llh[1] * DEG2RAD);
   sindc = sin (dec * DEG2RAD);
   cosdc = cos (dec * DEG2RAD);
   sinra = sin (ra * 15.0 * DEG2RAD);
   cosra = cos (ra * 15.0 * DEG2RAD);

/*
   Set up orthonormal basis vectors in local Earth-fixed system.

   Define vector toward local zenith in Earth-fixed system (z axis).
*/
   uze[0] = coslat * coslon;
   uze[1] = coslat * sinlon;
   uze[2] = sinlat;

/*
   Define vector toward local north in Earth-fixed system (x axis).
*/

   une[0] = -sinlat * coslon;
   une[1] = -sinlat * sinlon;
   une[2] = coslat;

/*
   Define vector toward local west in Earth-fixed system (y axis).
*/

   uwe[0] = sinlon;
   uwe[1] = -coslon;
   uwe[2] = 0.0;

/*
   Obtain vectors in celestial system.

   Rotate Earth-fixed orthonormal basis vectors to celestial system
   (wrt equator and equinox of date).
*/

   ter2cel (jd_ut1,0.0,delta_t,1,accuracy,1,x,y,uze, uz);
   ter2cel (jd_ut1,0.0,delta_t,1,accuracy,1,x,y,une, un);
   ter2cel (jd_ut1,0.0,delta_t,1,accuracy,1,x,y,uwe, uw);

/*
   Define unit vector 'p' toward object in celestial system
   (wrt equator and equinox of date).
*/

   p[0] = cosdc * cosra;
   p[1] = cosdc * sinra;
   p[2] = sindc;

/*
   Compute coordinates of object wrt orthonormal basis.

   Compute components of 'p' - projections of 'p' onto rotated
   Earth-fixed basis vectors.
*/

   pz = p[0] * uz[0] + p[1] * uz[1] + p[2] * uz[2];
   pn = p[0] * un[0] + p[1] * un[1] + p[2] * un[2];
   pw = p[0] * uw[0] + p[1] * uw[1] + p[2] * uw[2];

/*
   Compute azimuth and zenith distance.
*/

   proj = sqrt (pn * pn + pw * pw);

   if (proj > 0.0)
      *az = -atan2 (pw, pn) * RAD2DEG;

   if (*az < 0.0)
      *az += 360.0;

   if (*az >= 360.0)
      *az -= 360.0;

   *zd = atan2 (proj, pz) * RAD2DEG;
}





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* fun_pointmass - abandoned
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

double fun_pointmass (double tdbs, double *x, double *f)
{
    double fnt[3], fgr[3], r, s2, rrd, a, b;
    int n, gamma;

    gamma = 1;

    f[0] = x[3]; 
    f[1] = x[4]; 
    f[2] = x[5];

    r = sqrt (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    s2 = x[3] * x[3] + x[4] * x[4] + x[5] * x[5];
    rrd = x[0] * x[3] + x[1] * x[4] + x[2] * x[5];
	
    a = 2 * (1 + gamma) * GMDE[CENTER] / r - gamma * s2;
    b = 2 * (1 + gamma) * rrd;
    for (n = 0; n < 3; n++)
        fgr[n] =  GMDE[CENTER] / C_AUDAY / C_AUDAY / r / r / r 
        * ( a * x[n] + b * x[n+3] );

    fnt[0] = - GMDE[CENTER] / (r*r*r) * x[0];
    fnt[1] = - GMDE[CENTER] / (r*r*r) * x[1];
    fnt[2] = - GMDE[CENTER] / (r*r*r) * x[2];

    for (n = 0; n < 3; n++)
    {
        f[3 + n] = fnt[n] + fgr[n]; 
    }

	return 0;
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* fun_fullaccel - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double fun_fullaccel (double tdbs, double *xic, double *fxic)
{
    int n;
    short int ssbary = 11, part = 0;
    double tjd[2], acc1[3], acc2[3], acc3[3], acc[3], dum1[1], dum2[1];

    tjd[0] = JD0;
    tjd[1] = tdbs;

    accel_ntrel (tjd, xic, part, acc1, dum1, dum2);
    accel_nonsp (tjd, xic, part, acc2, dum1, dum2);
    accel_radpr (tjd, xic, part, acc3, dum1, dum2);

    for (n = 0; n <= 2; n++)
    {
        acc[n] = acc1[n] + acc2[n] + acc3[n];
    }

    fxic[0] = xic[3];
    fxic[1] = xic[4];
    fxic[2] = xic[5];
    fxic[3] = acc[0];
    fxic[4] = acc[1];
    fxic[5] = acc[2];

    return 0;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* fun_fullstate -transition matrix(36), orbit(6), sensitivity matrix(6*DYNPAR)
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double fun_fullstate (double tdbs, double *state, double *fstate)
{
    int n;
    short int ssbary = 11, part = 1;
    double tjd[2], xic[6], dfdx[36], dxdx0[36], dfdp[6 * DYNPAR], 
        dfdpp[6 * DYNPAR], dxdp[6 * DYNPAR],
        acc1[3], dadr1[9], dadp1[3 * DYNPAR],
        acc2[3], dadr2[9], dadp2[3 * DYNPAR],
        acc3[3], dadr3[9], dadp3[3 * DYNPAR],
        acc[3], dadr[9], dadp[3 * DYNPAR],
        fxic[6], fdxdx0[36], fdxdp[6 * DYNPAR];

    tjd[0] = JD0;
    tjd[1] = tdbs;

    for (n = 0; n < 36; n++)
    {
        dxdx0[n] = state[n];
    }
    for (n = 0; n < 6; n++)
    {
        xic[n] = state[n + 36];
    }
    for (n = 0; n < 6 * DYNPAR; n++)
    {
        dxdp[n] = state[n + 42];
    }

/* acc, partial to xyz: dadr, partial to parameters dadp*/
    accel_ntrel (tjd, xic, part, acc1, dadr1, dadp1);
    accel_nonsp (tjd, xic, part, acc2, dadr2, dadp2);
    accel_radpr (tjd, xic, part, acc3, dadr3, dadp3);
/*todo: air drag acc & partial to vxvyvz dadv*/

    for (n = 0; n <= 2; n++)
    {
        acc[n] = acc1[n] + acc2[n] + acc3[n];
    }
    for (n = 0; n <= 8; n++)
    {
        dadr[n] = dadr1[n] + dadr2[n] + dadr3[n];
    }
    for (n = 0; n < 3 * DYNPAR; n++)
    {
        dadp[n] = dadp1[n] + dadp2[n] + dadp3[n];
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    for (n = 0; n < 36; n++)
    {
        dfdx[n] = 0;
    }
    dfdx[3]  = 1; 
    dfdx[10] = 1; 
    dfdx[17] = 1;
    for (n = 0; n < 3; n++)
    {
        dfdx[n + 18] = dadr[n];
        dfdx[n + 24] = dadr[n + 3];
        dfdx[n + 30] = dadr[n + 6];
    }
    brmul(dfdx, dxdx0, 6, 6, 6, fdxdx0);

    fxic[0] = xic[3];
    fxic[1] = xic[4];
    fxic[2] = xic[5];
    fxic[3] = acc[0];
    fxic[4] = acc[1];
    fxic[5] = acc[2];

    brmul(dfdx, dxdp, 6, 6, DYNPAR, dfdpp);
    for (n = 0; n < 3 * DYNPAR; n++)
    {
        dfdp[n] = 0;
        dfdp[n + 3 * DYNPAR] = dadp[n];
    }
    for (n = 0; n < 6 * DYNPAR; n++)
    {
        fdxdp[n] = dfdpp[n] + dfdp[n];
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    for (n = 0; n < 36; n++)
    {
        fstate[n] = fdxdx0[n];
    }
    for (n = 0; n < 6; n++)
    {
        fstate[n + 36] = fxic[n];
    }
    for (n = 0; n < 6 * DYNPAR; n++)
    {
        fstate[n + 42]= fdxdp[n];
    }

    return 0;
}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* accel_ntrel - Newtonian + Relativistic acceleration
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double accel_ntrel (double *tjd, double *xic, short int part, 
                    double *acc, double *dadr, double *dadp)
{
    int n;
    short int ssbary = 11;
    double xcb[6], acb[3], xib[6], aib[3], dadr1[9], dum[9], 
        dadp0[3*DYNPAR], dadp1[3*DYNPAR];

    planet_ephemeris (tjd, CENTER, ssbary, &xcb[0], &xcb[3]);
    accel_bcrs (tjd, xcb, part, CENTER, acb, dum, dadp0);
    for (n = 0; n <= 5; n++)
    {
        xib[n] = xic[n] + xcb[n];
    }
    accel_bcrs (tjd, xib, part, 99, aib, dadr1, dadp1);
    for (n = 0; n <= 2; n++)
    {
        acc[n] = aib[n] - acb[n];
    }

    if (part == 1)
    {
        for (n = 0; n <= 8; n++)
        {
            dadr[n] = dadr1[n];
        }
        for (n = 0; n <= 3 * DYNPAR - 1; n++)
        {
            dadp[n] = dadp1[n] - dadp0[n];
        }
    }
    return 0;
}






/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* accel_bcrs - Newtonian + Relativistic acceleration
* @param1: description of param1
* @param2: description of param2
* todo: 
        1 partial to parameters
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double accel_bcrs (double *jd, double *xi, short int part, short int exclude, 
                   double *acc, double *dadr, double *dadp)
{
    double fnt[3], fgr[3], xj[11][6], xij[11][6], rij[11], xjk[6], rjk, 
        xddj[3], sumil, sumjk, sdi2, sdj2, rdirdj, rrrdr2, rjirdd, 
        rij5, rij3, xijt[9], gra, grb, beta, gamma, unit[9];
    short int ssbary, l, k, j, n, flag_gr;
    
    ssbary = 11;
    gamma = 1.0;
    beta = 1.0;
    unit[0] = 1; unit[1] = 0; unit[2] = 0;
    unit[3] = 0; unit[4] = 1; unit[5] = 0;
    unit[6] = 0; unit[7] = 0; unit[8] = 1;

    for (j = 0; j <= 10; j++)
    {
        planet_ephemeris (jd, j, ssbary, &xj[j][0], &xj[j][3]);
        for (n = 0; n < 6; n++)
        {
            xij[j][n] = xi[n] - xj[j][n];
        }
        rij[j] = sqrt (xij[j][0] * xij[j][0] 
            + xij[j][1] * xij[j][1] + xij[j][2] * xij[j][2]);
    }
    
    flag_gr = 0;
    for (n = 0; n < 3; n ++)
        fnt[n] = 0;
    for (j = 0; j <= 10; j++)
    {
        if (PERB[j] == 2)
            flag_gr = 1;
        if (PERB[j] == 0)
            continue;
        if (j == exclude)
            continue;
        for (n = 0; n < 3; n++)
            fnt[n] = fnt[n] 
            - GMDE[j] / (rij[j] * rij[j] * rij[j]) * xij[j][n];
    }

    if (part == 1)
    {
        for (n = 0; n <= 3 * DYNPAR - 1; n++)
        {
            dadp[n] = 0;
        }
        for (n = 0; n <= 8; n++)
        {
            dadr[n] = 0;
        }
        for (j = 0; j <= 10; j++)
        {
            if (j == exclude)
                continue;
            rij5 = pow (rij[j], 5);
            rij3 = pow (rij[j], 3);
            brmul (xij[j], xij[j], 3,1,3, xijt);
            for (n = 0; n <= 8; n++)
            {
                dadr[n] = dadr[n] + 3 * GMDE[j] * xijt[n] / rij5
                - GMDE[j] * unit[n] / rij3;
            }
        }
    }

    if (flag_gr == 0)
    {
        for (n = 0; n < 3; n++)
            acc[n] =  fnt[n];
        return 0;
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    sdi2 = xi[3] * xi[3] + xi[4] * xi[4] + xi[5] * xi[5];
    sumil = 0;
    for (l = 0; l < 11; l ++)
    {
        if ( l == exclude)
            continue;
        if (PERB[l] != 2)
            continue;
        sumil = sumil + GMDE[l] / rij[l];
    }

    for (n = 0; n < 3; n ++)
        fgr[n] = 0;
    for (j = 0; j < 11; j ++)
    {
        if (PERB[j] != 2)
            continue;
        if (j == exclude)
            continue;
        sumjk = 0;
        for (n = 0; n < 3; n ++)
            xddj[n] = 0;
        for (k = 0; k < 11; k ++)
        {
            if (k == j)	
                continue;	//k!=j
            if (PERB[k] != 2)
                continue;
            for (n = 0; n < 3; n++)
                xjk[n] = xj[j][n] - xj[k][n];
            rjk = sqrt (xjk[0] * xjk[0] + xjk[1] * xjk[1] + xjk[2] * xjk[2]);
            sumjk = sumjk + GMDE[k] / rjk;
            for (n = 0; n < 3; n ++)
                xddj[n] = xddj[n] - GMDE[k] / (rjk * rjk * rjk) * xjk[n];
        }
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
        sdj2 = xj[j][3] * xj[j][3] + xj[j][4] * xj[j][4] 
            + xj[j][5] * xj[j][5];
        rdirdj = xi[3] * xj[j][3] + xi[4] * xj[j][4] + xi[5] * xj[j][5];
        rrrdr2 = pow( ( xij[j][0] * xj[j][3] + xij[j][1] * xj[j][4] 
            + xij[j][2] * xj[j][5]) / rij[j], 2);
        rjirdd = - ( xij[j][0] * xddj[0] + xij[j][1] * xddj[1] 
            + xij[j][2] * xddj[2]);
        
        gra = - 2 * (beta + gamma) * sumil - (2 * beta -1) * sumjk 
            + gamma * sdi2 + (1 + gamma) * sdj2
            - 2 * (1 + gamma) * rdirdj - 1.5 * rrrdr2 + 0.5 * rjirdd;

        grb = xij[j][0] * ((2+2*gamma) * xi[3] - (1+2*gamma) * xj[j][3])
            + xij[j][1] * ((2+2*gamma) * xi[4] - (1+2*gamma) * xj[j][4])
            + xij[j][2] * ((2+2*gamma) * xi[5] - (1+2*gamma) * xj[j][5]);

        for (n = 0; n < 3; n ++)
        {
            fgr[n] = fgr[n] 
                + GMDE[j] / (rij[j] * rij[j] * rij[j]) 
                * ( - xij[j][n]) * gra / C_AUDAY / C_AUDAY
                + GMDE[j] / (rij[j] * rij[j] * rij[j]) 
                * xij[j][n + 3] * grb / C_AUDAY / C_AUDAY
                + GMDE[j] / rij[j] * (3 + 4 * gamma) * 0.5   
                * xddj[n] / C_AUDAY / C_AUDAY;
        }
    }

    for (n = 0; n < 3; n++)
        acc[n] =  fgr[n] + fnt[n];
    return 1;
}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* accel_radpr - solar radiation press & partial to srp coefficients 
* @param1: description of param1
* @param2: description of param2
* todo: 
        1 earth shadow 
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double accel_radpr (double *tjd, double *xic, short int part, 
                    double *acc, double *dadr, double *dadp)
{
    double j, c1, ap, m, rsp, usp[3], xis[6], xsc[6], f, 
        xist[9], unit[9], rsp3;
    short int n, sun;

    sun = 10;
    unit[0] = 1; unit[1] = 0; unit[2] = 0;
    unit[3] = 0; unit[4] = 1; unit[5] = 0;
    unit[6] = 0; unit[7] = 0; unit[8] = 1;

    planet_ephemeris (tjd, sun, CENTER, &xsc[0], &xsc[3]);
    for (n = 0; n <= 5; n++)
    {
        xis[n] = xic[n] - xsc[n];
    }
    rsp = sqrt (xis[0] * xis[0] + xis[1] * xis[1] + xis[2] * xis[2]);
    usp[0] = xis[0] / rsp;
    usp[1] = xis[1] / rsp;
    usp[2] = xis[2] / rsp;

    j  = 1352.5;   //kg/s3
    m  = SATMASS;     //kg
    ap = SATAREA;        //m2
    c1 = j / C * 1 * 1;   //kg/s2/m*au*au
    f  = c1 * ap / m  / rsp / rsp;    
//kg/s2/m*au*au * m2 / kg / au  / au = m/s2
    f = f / AU * 86400.0 * 86400.0;

//    acc[0] = f * usp[0] * (1 + CONS + DCON * tjd[1]);
//    acc[1] = f * usp[1] * (1 + CONS + DCON * tjd[1]);
//    acc[2] = f * usp[2] * (1 + CONS + DCON * tjd[1]);
    acc[0] = f * usp[0];
    acc[1] = f * usp[1];
    acc[2] = f * usp[2]; 

    if (part == 0)
        return 1;

    rsp3 = rsp * rsp * rsp;
    brmul (xis, xis, 3,1,3, xist);
    for (n = 0; n <= 8; n++)
        dadr[n] = - f * (3 * xist[n] / rsp3 - unit[n] / rsp) ;
//        * (1 + CONS + DCON * tjd[1]);

//    for (n = 0; n <= 2; n++)
//    {
//        dadp[n * DYNPAR] = f * usp[n];      
//        dadp[n * DYNPAR + 1] = dadp[n * DYNPAR] * tjd[1];   
//        dadp[n * DYNPAR] = 0;      
//        dadp[n * DYNPAR + 1] = 0;   

//    }

        for (n = 0; n <= 3 * DYNPAR - 1; n++)
        {
            dadp[n] = 0;
        }


    return 0;
}






//    nmax = 4;      
//    stcs = (double *) calloc ( (nmax + 1) * (nmax + 1), sizeof(double));
double stidecs_old(double *tjd, double gma1, double k2, 
               double *c20, double *c21, double *s21, double *c22, double *s22)
{

    double gms2e;                   
    double gmm2e;
    

//    short int moon = 9, earth = 2, sun = 10;
    short int moon = 2, earth = 9, sun = 10;

    double ps[3], vs[3], pm[3], vm[3],
        pse[3], pme[3], llrs[3], llrm[3], pbar[4], t,
        p20m, p30m, p21m, p31m, p22m, p32m, p33m, 
        p20s, p30s, p21s, p31s, p22s, p32s, p33s,
        rerm, rers, tb[9], tbt[9];


// Luni-solar ephemeris

    planet_ephemeris (tjd, sun, CENTER, ps, vs);
    planet_ephemeris (tjd, moon, CENTER, pm, vm);

    iau_pns (tjd, tb, CENTER);	
    mt (tb, 3, 3, tbt);			
    brmul (tbt,ps,3,3,1,pse);	
    brmul (tbt,pm,3,3,1,pme);	

    
    xyz2llh(pse, llrs);
    xyz2llh(pme, llrm);

    t = sin(llrm[0] * DEG2RAD);
    lgdr(t, 3, 0, pbar); p20m = pbar[2]; p30m = pbar[3];
    lgdr(t, 3, 1, pbar); p21m = pbar[1]; p31m = pbar[2];
    lgdr(t, 3, 2, pbar); p22m = pbar[0]; p32m = pbar[1];
    lgdr(t, 3, 3, pbar); p33m = pbar[0];
    t = sin(llrs[0] * DEG2RAD);
    lgdr(t, 3, 0, pbar); p20s = pbar[2]; p30s = pbar[3];
    lgdr(t, 3, 1, pbar); p21s = pbar[1]; p31s = pbar[2];
    lgdr(t, 3, 2, pbar); p22s = pbar[0]; p32s = pbar[1];
    lgdr(t, 3, 3, pbar); p33s = pbar[0];


    gms2e  =  GMDE[sun]/GMDE[CENTER];                   
//    gmm2e  =  GMDE[moon]/GMDE[CENTER]; 
    gmm2e = 0;

    rerm = gma1 / llrm[2];
    rers = gma1 / llrs[2];

// Frequency Independent Terms

// C20
    *c20 = k2/5.0 * ( gmm2e * pow(rerm, 3) * p20m 
                    + gms2e * pow(rers, 3) * p20s );
// C21/S21
    *c21 = k2/5.0 * ( gmm2e * pow(rerm, 3) * p21m * cos(llrm[1] * DEG2RAD)
                    + gms2e * pow(rers, 3) * p21s * cos(llrs[1] * DEG2RAD) );
    *s21 = k2/5.0 * ( gmm2e * pow(rerm, 3) * p21m * sin(llrm[1] * DEG2RAD) 
                    + gms2e * pow(rers, 3) * p21s * sin(llrs[1] * DEG2RAD) );
// C22/S22
    *c22 = k2/5.0 * ( gmm2e * pow(rerm, 3) * p22m * cos(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * cos(llrs[1] * DEG2RAD * 2.0) );
    *s22 = k2/5.0 * ( gmm2e * pow(rerm, 3) * p22m * sin(llrm[1] * DEG2RAD * 2.0)
                    + gms2e * pow(rers, 3) * p22s * sin(llrs[1] * DEG2RAD * 2.0) );
  
    
    return 0;

}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* accel_nonsp - non-spherical force 
* @param1: description of param1
* @param2: description of param2
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double accel_nonsp (double *tjd, double *xic, short int part, 
                    double *acc, double *dadr, double *dadp)
{
    int n, m, i;
    static int flag = 0;
//    static double *cn0, *cnm, *snm, gma[2];

    double xfc[3] = {0}, tb[9] = {0}, tbt[9] ={0}, rmat[9] = {0}, 
        gmat[9] = {0}, gmatt[9] = {0}, r, rxy, sinf, cosf, sinl, cosl, 
        lamta, *pn, *pnm, *pnp, *pnmp, *pnpp, *pnmpp,
        frj[3] = {0}, frcs[3] = {0}, fr[3] = {0}, *cosml, *sinml, *aprn,
        peprp[27], pepr[27], prtpx[9], prtpy[9], prtpz[9], pgtpx[9], 
        pgtpy[9], pgtpz[9],	pgx[3], pgy[3], pgz[3], part1[9], prjpx[3], 
        prjpy[3], prjpz[3], prcspx[3], prcspy[3], prcspz[3], prpr[9], 
        gtpr[9], part2[9], prtpxx[9], prtpyy[9], prtpzz[9];

    double cunit, dfd2r[3], dfd2[3], dfdkr[3], dfdk[3], k2, c20, c21, s21, c22, s22,
        unit, nup, ndown;

    if (part == 1)
    {
        for (n = 0; n <= 3 * DYNPAR - 1; n++)
        {
            dadp[n] = 0;
        }
        for (n = 0; n <= 8; n++)
        {
            dadr[n] = 0;
        }
    }

    if (GRAVDEGREE < 2)
    {
        for (n = 0; n <= 2; n++)
        {
            acc[n] = 0;
        }
        return 1;
    }

    if (flag != 9)			//
    {
        cn0  = (double *) calloc (GRAVDEGREE, sizeof(double));
        cnm  = (double *) calloc (GRAVDEGREE * GRAVDEGREE, sizeof(double));
        snm  = (double *) calloc (GRAVDEGREE * GRAVDEGREE, sizeof(double));
        opengravfile (cn0, cnm, snm, gma);
        flag = 9;
    }


    k2 = CONS;
    stidecs_old(tjd, gma[1], k2, &c20, &c21, &s21, &c22, &s22);

//    cn0[1] = (-8.745054708184200e-04 + c20 + DCON * 1.0e-16 * tjd[1]) * sqrt(5);
    n = 2;
    cn0[n-1] = (j2 + c20 + DCON * 1.0e-8) * sqrt(2*n+1);
    n = 3;
    cn0[n-1] = (j3 + BIAS * 1.0e-8) * sqrt(2*n+1);
    n = 4;
    cn0[n-1] = (j4 + DBIA * 1.0e-8) * sqrt(2*n+1);


    n = 2; m = 1;
    {
            nup = 1.0;
            ndown = 1.0;
            for (i = 1; i <= n - m; i++)
                nup = nup * i;
            for (i = 1; i <= n + m; i++)
                ndown = ndown * i;
            unit = sqrt (2 * (2*n+1.0)*nup/ndown);
            cnm[(n-1)*GRAVDEGREE + (m-1)] = (jc21 + c21)* unit;
            snm[(n-1)*GRAVDEGREE + (m-1)] = (js21 + s21)* unit;
    }

    n = 2; m = 2;
    {
            nup = 1.0;
            ndown = 1.0;
            for (i = 1; i <= n - m; i++)
                nup = nup * i;
            for (i = 1; i <= n + m; i++)
                ndown = ndown * i;
            unit = sqrt (2 * (2*n+1.0)*nup/ndown);
            cnm[(n-1)*GRAVDEGREE + (m-1)] = (jc22 + c22)* unit;
            snm[(n-1)*GRAVDEGREE + (m-1)] = (js22  + s22)* unit;
    }

    //    cn0[2] = -1.188691064601560e-05 * sqrt(7) * (1 + DCON);

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*****************************rotation matrix ********************************/
    iau_pns (tjd, tb, CENTER);	//tb
    mt (tb, 3, 3, tbt);			//tbt
    brmul (tbt,xic,3,3,1,xfc);	//rb
	
    r = sqrt (xfc[0] * xfc[0] + xfc[1] * xfc[1] + xfc[2] * xfc[2]);	
    //define up-east-north system 
    rxy  = sqrt (xfc[0] * xfc[0] + xfc[1] * xfc[1]);
    sinf = xfc[2] / r;
    cosf = rxy / r;
    sinl = xfc[1] / rxy;
    cosl = xfc[0] / rxy;

    rmat[0] = cosf * cosl;		//from fixed to up-east-north system: rmat
    rmat[1] = cosf * sinl;
    rmat[2] = sinf;
    rmat[3] = -sinl;
    rmat[4] = cosl;
    rmat[5] = 0;
    rmat[6] = -sinf * cosl;
    rmat[7] = -sinf * sinl;
    rmat[8] = cosf;

    brmul (rmat,tbt,3,3,3,gmat);	//inertial to fixed matrix gmat = rmat*tbt
    mt (gmat, 3, 3, gmatt);		//fixed to inertial matrix gmatt

    lamta = chosephase (sinl, cosl); //rad

    cosml = (double *) calloc ( GRAVDEGREE, sizeof(double)); //cos(m*lamta)
    sinml = (double *) calloc ( GRAVDEGREE, sizeof(double)); //sin(m*lamta)
    aprn = (double *) calloc ( GRAVDEGREE, sizeof(double));  //sin(m*lamta)

    for (m = 1; m <= GRAVDEGREE; m++)
    {
        cosml[m-1] = cos(m*lamta);
        sinml[m-1] = sin(m*lamta);
    }
	
    for (n = 1; n <= GRAVDEGREE; n++)
    {
        aprn[n-1] = pow (gma[1] / r, n);
    }

/******************************************************************/

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/****************Legendre Polynomial**********************************/
    pn    = (double *) calloc ( GRAVDEGREE, sizeof(double));				
    //Pn
    pnp   = (double *) calloc ( GRAVDEGREE, sizeof(double));				
    //Pn'
    pnm   = (double *) calloc ( GRAVDEGREE * GRAVDEGREE, sizeof(double));	
    //secf*Pmn
    pnmp  = (double *) calloc ( GRAVDEGREE * GRAVDEGREE, sizeof(double));	
    //cosf*Pmn'
    pnpp  = (double *) calloc ( GRAVDEGREE, sizeof(double));				
    //Pn''
    pnmpp = (double *) calloc ( GRAVDEGREE * GRAVDEGREE, sizeof(double));	
    //cos2fPmn''

    pn[0]   = sinf; 
    pnp[0]  = 1; 
    pnpp[0] = 0;
    pn[1]   = 3.0/2.0*sinf*sinf - 1.0/2.0;
    pnp[1]  = sinf + 2 * sinf;
    pnpp[1] = 3;
    
    for (n = 3; n <= GRAVDEGREE; n++)
    {
        pn[n-1] = (2 * n - 1.0) / n * sinf * pn[n-2] 
            - (n - 1.0) / n * pn[n-3]; //tmd!!!
        pnp[n-1] = sinf * pnp[n-2] + n * pn[n-2];
        pnpp[n-1] = sinf * pnpp[n-2] + (n+1) * pnp[n-2];
    }

    pnm[0] = 1; //secfP11 = 1
    for (n = 2; n <= GRAVDEGREE; n++)
    {
        pnm[(n-1) * GRAVDEGREE + n - 1] 
            = (2 * n - 1.0) * cosf * pnm[(n - 2) * GRAVDEGREE + (n - 2)];
    }
	
    pnm[GRAVDEGREE] = (2 * 2.0 - 1.0) / (2 - 1.0) * sinf * pnm[0];	
    //secfP21 = pnm[GRAVDEGREE]
    for (n = 3; n <= GRAVDEGREE; n++)
    {
        for (m = 1; m < n; m++)
        {
            pnm[(n-1) * GRAVDEGREE + (m-1)] = (2 * n - 1.0) / (n-m) 
                * sinf * pnm[(n-2) * GRAVDEGREE + (m-1)]
                - (n + m - 1.0) / (n - m) 
                * pnm[(n - 3) * GRAVDEGREE + (m - 1)];
//			printf ("%d\t%d\t%f\n", n, m, pnm[(n-1)*n2 + (m-1)]);
        }
    }

    pnmp[0] = -sinf * pnm[0];		//cosfP11'
    for (n = 2; n <= GRAVDEGREE; n++)
    {
        for (m = 1; m <= n; m++)
        {
            pnmp[(n - 1) * GRAVDEGREE + (m - 1)] = 
                - n * sinf * pnm[(n - 1) * GRAVDEGREE + (m - 1)] 
                + (n + m) * pnm[(n - 2) * GRAVDEGREE + (m - 1)];
        }
    }
    
    pnmpp[0] = sinf * pnmp[0] / cosf - pnm[0] * cosf;	
    //cos2fP11''
    pnmpp[GRAVDEGREE] = sinf * pnmp[GRAVDEGREE] / cosf 
        - pnm[GRAVDEGREE] * cosf - 3 * sinf * pnm[GRAVDEGREE + 1];	
    //cos2fP21'' = pnmpp[GRAVDEGREE]
    pnmpp[GRAVDEGREE + 1] = - 2 * pnm[GRAVDEGREE + 1] * cosf;	
    //cos2fP22'' = pnmpp[GRAVDEGREE+1]
    for (n = 3; n <= GRAVDEGREE; n++)
    {
        for (m = 1; m <= n; m++)
        {
            if (m == 1)
            {
                pnmpp[(n - 1) * GRAVDEGREE + (m - 1)] = 
                    + sinf * pnmp[(n - 1) * GRAVDEGREE + (m - 1)] / cosf 
                    - pnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosf
                    - 3 * sinf * pnm[(n - 1) * GRAVDEGREE + (m - 1) + 1] 
                    + pnm[(n - 1) * GRAVDEGREE + (m - 1) + 2] * cosf;
            }
            else
            {
                pnmpp[(n-1)*GRAVDEGREE + (m-1)] = - (n - 2) * sinf 
                    * pnmp[(n - 1) * GRAVDEGREE + (m-1)] / cosf
                    + (n + m) * pnmp[(n - 2) * GRAVDEGREE + (m-1)] / cosf 
                    - n * pnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosf;
            }
        }
    }
/*****************Legendre Polynomial**********************************/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

/**********************************************************/
    for (n = 1; n <= GRAVDEGREE; n++)
    {
        frj[0] = frj[0] + (-cn0[n - 1]) * aprn[n - 1] * (n + 1) * pn[n - 1];
        frj[1] = frj[1] + 0;
        frj[2] = frj[2] + (-cn0[n - 1]) * aprn[n - 1] * (-cosf) * pnp[n - 1];
    }
    for (n = 1; n <= GRAVDEGREE; n++)
    {
        for (m = 1; m <= n; m++)
        {
            if ( n == GRAVDEGREE && m > GRAVORDER)
            {
//				printf ("%d\t%d\n",n,m);
                break;
            }
            frcs[0] = frcs[0] + aprn[n - 1] 
                * ( - (n + 1)) * pnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosf
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);

            frcs[1] = frcs[1] + aprn[n - 1] 
                * m * pnm[(n - 1) * GRAVDEGREE + (m - 1)] 
                * (-cnm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1]);
            frcs[2] = frcs[2] + aprn[n-1] 
                * pnmp[(n - 1) * GRAVDEGREE + (m - 1)] 
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);
        }
    }

    for (n = 0; n < 3; n++)
    {
        fr[n] = (frj[n] + frcs[n]) * gma[0] / r / r;
    }

    brmul(gmatt,fr,3,3,1,acc);  //from fixed acc to inertial acc

    if (part == 0)
    {
        free (pn);
        free (pnp);
        free (pnm);
        free (pnmp);
        free (pnpp);
        free (pnmpp);
        free (cosml);
        free (sinml);
        free (aprn);
        return 1;
    }


/*************************************************************/

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/************************partial part1***********************************/

    peprp[0]  = 0;	
    peprp[1]  = - sinl / r;					
    peprp[2]  = - sinf * cosl / r;	//(573)
    peprp[3]  = 0;	
    peprp[4]  = cosl / r;					
    peprp[5]  = - sinf * sinl / r;	//(574)		//11.10change
    peprp[6]  = 0;	
    peprp[7]  = 0;							
    peprp[8]  = cosf / r;			//(575)
    peprp[9]  = 0;	
    peprp[10] = - cosl / r / cosf;			
    peprp[11] = 0;					//(576)
    peprp[12] = 0;	
    peprp[13] = - sinl / r / cosf;			
    peprp[14] = 0;					//(577)
    peprp[15] = 0;	
    peprp[16] = 0;							
    peprp[17] = 0;					//(578)
    peprp[18] = 0;	
    peprp[19] = sinf * sinl / r / cosf;		
    peprp[20] = - cosf * cosl / r;	//(579)
    peprp[21] = 0;	
    peprp[22] = - sinf * cosl / r / cosf;	
    peprp[23] = - cosf * sinl / r;	//(580)
    peprp[24] = 0;	
    peprp[25] = 0;							
    peprp[26] = - sinf / r;			//(581)

    brmul(peprp, gmat, 9, 3, 3, pepr);		//571

    for (n = 0; n < 9; n++)					//570
    {
        prtpx[n] = pepr[n*3];
        prtpy[n] = pepr[n*3+1];
        prtpz[n] = pepr[n*3+2];
    }

    mt (prtpx, 3, 3, prtpxx);			//OH MY GOD!!!11.11
    mt (prtpy, 3, 3, prtpyy);			//OH MY GOD!!!11.11
    mt (prtpz, 3, 3, prtpzz);			//OH MY GOD!!!11.11


    brmul(tb, prtpxx, 3, 3, 3, pgtpx);		//568
    brmul(tb, prtpyy, 3, 3, 3, pgtpy);  
    brmul(tb, prtpzz, 3, 3, 3, pgtpz);  

    brmul(pgtpx, fr, 3, 3, 1, pgx);			//558 first term	//11.10 change
    brmul(pgtpy, fr, 3, 3, 1, pgy);  
    brmul(pgtpz, fr, 3, 3, 1, pgz);  

    for (n = 0; n < 3; n++)	
    {
        part1[n*3] = pgx[n]; 
        part1[n*3+1] = pgy[n]; 
        part1[n*3+2] = pgz[n];
    }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/************************partial part2*************************************/

    for (n = 0; n < 3; n++)
    {
        prjpx[n]  = 0;	
        prjpy[n]  = 0;	
        prjpz[n]  = 0;
        prcspx[n] = 0;	
        prcspy[n] = 0;	
        prcspz[n] = 0;
    }

    for (n = 1; n <= GRAVDEGREE; n++)
    {
        prjpx[0] = prjpx[0] - (-cn0[n - 1]) * aprn[n - 1] 
            * (n + 1) * pn[n - 1] * (n + 2);        //561
        prjpx[2] = prjpx[2] - (-cn0[n - 1]) * aprn[n - 1] 
            * (-cosf) * pnp[n - 1] * (n + 2);       //561
        prjpz[0] = prjpz[0] + (-cn0[n - 1]) * aprn[n - 1] 
            * (n + 1) * cosf * pnp[n - 1];      //563
		prjpz[2] = prjpz[2]	+ (-cn0[n - 1]) * aprn[n - 1] 
            * ( sinf * pnp[n - 1] - cosf * cosf * pnpp[n - 1] );//563
    }
	
    for (n = 1; n <= GRAVDEGREE; n++)
    {
        for (m = 1; m <= n; m++)
        {
            if ( n == GRAVDEGREE && m > GRAVORDER)
            {
//				printf ("%d\t%d\n",n,m);
                break;
            }
//from 564 to 566
            prcspx[0] = prcspx[0] - aprn[n - 1] 
                * ( - (n + 1)) * pnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosf
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]) 
                * (n + 2);			
            prcspx[1] = prcspx[1] - aprn[n - 1] 
                * m * pnm[(n - 1) * GRAVDEGREE + (m - 1)] 
                * ( - cnm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1]) 
                * (n + 2);
            prcspx[2] = prcspx[2] - aprn[n - 1] 
                * pnmp[(n - 1) * GRAVDEGREE + (m - 1)] 
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1])
                * (n + 2);

            prcspy[0] = prcspy[0] + m * aprn[n - 1] 
                * (n + 1) * pnm[(n - 1) * GRAVDEGREE + (m - 1)]
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1] 
                - snm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1]);
            prcspy[1] = prcspy[1] + m * aprn[n - 1] 
                * m * pnm[(n - 1) * GRAVDEGREE + (m - 1)] / cosf
                * ( - cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                - snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);				
            prcspy[2] = prcspy[2] + m * aprn[n - 1] 
                * pnmp[(n - 1) * GRAVDEGREE + (m - 1)] / cosf
                * ( - cnm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1]);

            prcspz[0] = prcspz[0] + aprn[n - 1] 
                * ( - (n + 1)) * pnmp[(n - 1) * GRAVDEGREE + (m - 1)]
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);
            prcspz[1] = prcspz[1] + aprn[n - 1] 
                * m * (sinf * pnm[(n - 1) * GRAVDEGREE + (m - 1)] 
                + pnmp[(n - 1) * GRAVDEGREE + (m - 1)]) / cosf
                * ( - cnm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1]);			
            prcspz[2] = prcspz[2] + aprn[n - 1] 
                * ( pnmpp[(n - 1) * GRAVDEGREE + (m - 1)] 
                - pnmp[(n - 1) * GRAVDEGREE + (m - 1)] * sinf / cosf)
                * ( cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);
        }
    }

    for (n = 0; n < 3; n++)	
    {
        prpr[n*3]   = (prjpx[n] + prcspx[n]) * gma[0] / r / r / r;
        prpr[n*3+1] = (prjpy[n] + prcspy[n]) * gma[0] / r / r / r;
        prpr[n*3+2] = (prjpz[n] + prcspz[n]) * gma[0] / r / r / r;
    }

    brmul(prpr, gmat, 3, 3, 3, gtpr);			
    brmul(gmatt, gtpr, 3, 3, 3, part2);	

    for (n = 0; n <= 8; n++)
    {
        dadr[n] = part1[n] + part2[n]; 
    }

/*****************************************************************************/




    n = 2;
    dfd2r[0] =  (- 1.0e-8 * sqrt(2*n+1)) * aprn[n - 1] * (n + 1) * pn[n - 1];
    dfd2r[1] =  0;
    dfd2r[2] =  (- 1.0e-8 * sqrt(2*n+1)) * aprn[n - 1] * (-cosf) * pnp[n - 1];

    for (n = 0; n < 3; n++)
    {
        dfd2r[n] = dfd2r[n] * gma[0] / r / r;
    }
  
    brmul(gmatt,dfd2r,3,3,1,dfd2);  


    for (n = 0; n <= 2; n++)
    {
//        dadp[n * DYNPAR] = dfd2[n];      
//        dadp[n * DYNPAR + 1] = dfd3[n];   
        dadp[n * DYNPAR + 1] = dfd2[n];   
    }

    n = 3;
    dfd2r[0] =  (- 1.0e-8 * sqrt(2*n+1)) * aprn[n - 1] * (n + 1) * pn[n - 1];
    dfd2r[1] =  0;
    dfd2r[2] =  (- 1.0e-8 * sqrt(2*n+1)) * aprn[n - 1] * (-cosf) * pnp[n - 1];

    for (n = 0; n < 3; n++)
    {
        dfd2r[n] = dfd2r[n] * gma[0] / r / r;
    }
  
    brmul(gmatt,dfd2r,3,3,1,dfd2);  


    for (n = 0; n <= 2; n++)
    {
//        dadp[n * DYNPAR] = dfd2[n];      
//        dadp[n * DYNPAR + 1] = dfd3[n];   
        dadp[n * DYNPAR + 2] = dfd2[n];   
    }

    n = 4;
    dfd2r[0] =  (- 1.0e-8 * sqrt(2*n+1)) * aprn[n - 1] * (n + 1) * pn[n - 1];
    dfd2r[1] =  0;
    dfd2r[2] =  (- 1.0e-8 * sqrt(2*n+1)) * aprn[n - 1] * (-cosf) * pnp[n - 1];

    for (n = 0; n < 3; n++)
    {
        dfd2r[n] = dfd2r[n] * gma[0] / r / r;
    }
  
    brmul(gmatt,dfd2r,3,3,1,dfd2);  


    for (n = 0; n <= 2; n++)
    {
//        dadp[n * DYNPAR] = dfd2[n];      
//        dadp[n * DYNPAR + 1] = dfd3[n];   
        dadp[n * DYNPAR + 3] = dfd2[n];   
    }


    k2 = 1;
    stidecs_old(tjd, gma[1], k2, &c20, &c21, &s21, &c22, &s22);

    n = 2;
    cn0[n-1] = (c20) * sqrt(2*n+1);

    dfdkr[0] =  (-cn0[n - 1]) * aprn[n - 1] * (n + 1) * pn[n - 1];
    dfdkr[1] =  0;
    dfdkr[2] =  (-cn0[n - 1]) * aprn[n - 1] * (-cosf) * pnp[n - 1];


    n = 2; m = 1;
    {
            nup = 1.0;
            ndown = 1.0;
            for (i = 1; i <= n - m; i++)
                nup = nup * i;
            for (i = 1; i <= n + m; i++)
                ndown = ndown * i;
            unit = sqrt (2 * (2*n+1.0)*nup/ndown);
            cnm[(n-1)*GRAVDEGREE + (m-1)] = (  c21)* unit;
            snm[(n-1)*GRAVDEGREE + (m-1)] = (  s21)* unit;

            dfdkr[0] = dfdkr[0] + aprn[n - 1] 
                * ( - (n + 1)) * pnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosf
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);

            dfdkr[1] = dfdkr[1] + aprn[n - 1] 
                * m * pnm[(n - 1) * GRAVDEGREE + (m - 1)] 
                * (-cnm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1]);
            dfdkr[2] = dfdkr[2] + aprn[n-1] 
                * pnmp[(n - 1) * GRAVDEGREE + (m - 1)] 
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);

    }

    n = 2; m = 2;
    {
            nup = 1.0;
            ndown = 1.0;
            for (i = 1; i <= n - m; i++)
                nup = nup * i;
            for (i = 1; i <= n + m; i++)
                ndown = ndown * i;
            unit = sqrt (2 * (2*n+1.0)*nup/ndown);
            cnm[(n-1)*GRAVDEGREE + (m-1)] = (   c22)* unit;
            snm[(n-1)*GRAVDEGREE + (m-1)] = (   s22)* unit;

            dfdkr[0] = dfdkr[0] + aprn[n - 1] 
                * ( - (n + 1)) * pnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosf
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);

            dfdkr[1] = dfdkr[1] + aprn[n - 1] 
                * m * pnm[(n - 1) * GRAVDEGREE + (m - 1)] 
                * (-cnm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1]);
            dfdkr[2] = dfdkr[2] + aprn[n-1] 
                * pnmp[(n - 1) * GRAVDEGREE + (m - 1)] 
                * (cnm[(n - 1) * GRAVDEGREE + (m - 1)] * cosml[m - 1] 
                + snm[(n - 1) * GRAVDEGREE + (m - 1)] * sinml[m - 1]);

    }


    for (n = 0; n < 3; n++)
    {
        dfdkr[n] = dfdkr[n] * gma[0] / r / r;
    }
  
    brmul(gmatt,dfdkr,3,3,1,dfdk);  


    for (n = 0; n <= 2; n++)
    {
        dadp[n * DYNPAR] = dfdk[n];      
//        dadp[n * DYNPAR + 1] = dfd3[n];   
//        dadp[n * DYNPAR + 1] = dfd2[n];   
    }



/*****************************************************************************/

    free (pn);
    free (pnp);
    free (pnm);
    free (pnmp);
    free (pnpp);
    free (pnmpp);
    free (cosml);
    free (sinml);
    free (aprn);
    return 0;
}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* opengravfile - open gravity field
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double opengravfile (double *cn0, double *cnm, double *snm, double *gma)
{
    FILE *fp_gra;
    double value, c,s, nup, ndown, unit;
    int n,m, i;
    char string[100], name[20];

    if ((fp_gra = fopen (FILE_GRV,"r")) == NULL)
    {
        printf ("Cannot open gravity file?\n");
      	getch();
        exit (0);
    }

    while (feof (fp_gra) == 0)
    {
        fgets (string, 100, fp_gra);
        sscanf (string, "%s%lf", name, &value);	
        if (strcmp (name,"Gm") ==0)	
        {
            gma[0] = value / AU / AU / AU * 86400.0 * 86400.0;
        }
        if (strcmp (name,"RefDistance") ==0)
        {
            gma[1] = value / AU;
        }
        if (strcmp (name,"BEGIN") ==0)	
            break;
    }
   
    while (feof(fp_gra) == 0)
    {
        n = 999;
        m = 999;
        fgets (string, 100, fp_gra);
        sscanf (string, "%d%d%lf%lf", &n, &m, &c, &s);	
        if (n > GRAVDEGREE)
            continue;
        else if (m == 0)
        {
            unit = sqrt (2 * n + 1.0);
            cn0[n-1] = c * unit;
            if (n == 2) j2 = c;
            if (n == 3) j3 = c;
            if (n == 4) j4 = c;

        }
        else 
        {
            if (n == 2 && m == 1) {jc21 = c; js21 = s; }
            if (n == 2 && m == 2) {jc22 = c; js22 = s; }

            nup = 1.0;
            ndown = 1.0;
            for (i = 1; i <= n - m; i++)
                nup = nup * i;
            for (i = 1; i <= n + m; i++)
                ndown = ndown * i;
            unit = sqrt (2 * (2*n+1.0)*nup/ndown);
            cnm[(n-1)*GRAVDEGREE + (m-1)] = c * unit;
            snm[(n-1)*GRAVDEGREE + (m-1)] = s * unit;
        }
    }

    fclose(fp_gra);
    return 0;
}







/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* itrf2icrf - from earth fixed to earth inertial
* @param1: description of param1
* @param2: description of param2 
        
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double itrf2icrf(double jd, double utc, double *vt, double *vc)
{
    double xp = 0, yp = 0, ut1_utc = 0, dx = 0, dy = 0, delta_t, ut1, tt;

    geteop (utc, &xp, &yp, &ut1_utc, &dx, &dy);	

    delta_t = 32.184 + LEAPSECS - ut1_utc;
    ut1 = utc + ut1_utc;
    tt = utc + (LEAPSECS + 32.184); 	

    cel_pole (jd + tt / 86400.0, 2, dx * 1e3, dy * 1e3);

    ter2cel (jd, ut1 / 86400.0, delta_t, 1, ACCURACY, 0,
        xp, yp, vt, vc); /*--vc unit: m--*/    

    return 0;


}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* iau_pns - planet fixed to J2000 inertial (for gravity field)
        Report of the IAU/IAGWorking Group on cartographic
        coordinates and rotational elements: 2006
* @param1: description of param1
* @param2: description of param2
* todo: 
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void iau_pns (double *jd, double *te, int cent) 
{
    double tes[9] = {0}, tepn[9] ={0}, tb[9], utc;

    double vx[3] = {1,0,0}, vy[3] = {0,1,0}, vz[3] = {0,0,1}, te2[9];
//    double ty[3], tz[3], te1[9], tx[3];
    int i;

    if (cent == 2)
    {
        utc = jd[1] * 86400 - (LEAPSECS + 32.184);  //jd[1]: tt(tdt)
        lagrange (TE_EPH, DIM_TE, 10, utc, te2);
/*
        itrf2icrf(jd[0], utc, vx, tx);
        itrf2icrf(jd[0], utc, vy, ty);
        itrf2icrf(jd[0], utc, vz, tz);
        for (i = 0; i < 3; i++)
        {
            te1[i*3] = tx[i];
            te1[i*3+1] = ty[i];
            te1[i*3+2] = tz[i];
        }
*/
        for (i = 0; i < 9; i++)
            te[i] = te2[i];


    }
    else if (cent == 9)
    {
        mbf2cel (jd, te);
//        in2pa (jd, tb);
//        mt (tb, 3, 3, te);
    }
    else
    {
        cent = cent +1;
        iau_s (jd, tes, cent);			//IAU fixed to IAU inertial
        iau_pn (jd, tepn, cent);		//IAU inertial to J2000 inertial
        brmul (tepn,tes,3,3,3,te);
    }
    return;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* iau_s - from IAU fixed to IAU inertial, true-of-date equator and equinox
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void iau_s (double *jd, double *tes, int cent)
{
    double d, str, cosst, sinst;	
    
    d = jd[0] - 2451545.0;
    d = d + jd[1];

    switch (cent)		//sun0, mercury1, ..., pluto9 
    {
    case 0 : str = 84.176 + 14.1844000 * d; break;
    case 1 : str = 329.548 + 6.1385025 * d; break;
    case 2 : str =  160.20 - 1.4813688 * d; break;
    case 3 : str = 190.147 + 360.9856235 * d; break;
    case 4 : str = 176.63 + 350.89198226 * d; break;	
    case 5 : str = 284.95 + 870.5366420 * d; break;
    case 6 : str =  38.90 + 810.7939024 * d; break;
    case 7 : str = 203.81 - 501.1600928 * d; break;
    case 8 : str = 253.18 + 536.3128492 * d - 
                 0.48 * sin ((357.85 + 52.316 * d / 36525.0 ) * DEG2RAD);
        break;
    case 9 : str = 237.305 - 56.3625225 * d; break;
    case I_TITAN : str = 186.5855 + 22.5769768 * d; break;
    }

    cosst = cos (str * DEG2RAD);
    sinst = sin (str * DEG2RAD);
	
    tes[0] = cosst;
    tes[1] = -sinst;
    tes[2] = 0;
    tes[3] = sinst;
    tes[4] = cosst;
    tes[5] = 0;
    tes[6] = 0;
    tes[7] = 0;
    tes[8] = 1;
    return;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* iau_pn - from IAU inertial (for all planets) to J2000 inertial
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void iau_pn (double *jd, double *tes, int cent)
{	
    double ra0, dec0, jcent, cr, sr, cd, sd;
    
    jcent = jd[0] - 2451545.0;
    jcent = (jcent + jd[1]) / 36525.0;
    
    switch (cent)		//sun0, mercury1, ..., pluto9 
    {
    case 0 : 
        ra0  = 286.13;
        dec0 = 63.87; 
        break;		
    case 1 : 
        ra0  = 281.01 - 0.033 * jcent;
        dec0 = 61.45 - 0.005 * jcent; 
        break;			
    case 2 : 
        ra0  = 272.76;
        dec0 = 67.16; 
        break;	
    case 3 : 
        ra0  = 0.00 - 0.641 * jcent;
        dec0 = 90.0 - 0.557 * jcent; 
        break;	
    case 4 : 
        ra0  = 317.68143 - 0.1061 * jcent;
        dec0 = 52.88650 - 0.0609 * jcent; 
        break;	
    case 5 : 
        ra0  = 268.05 - 0.009 * jcent;
        dec0 = 64.49 + 0.003 * jcent; 
        break;	
    case 6 : 
        ra0  = 40.589 - 0.036 * jcent;
        dec0 = 83.537 - 0.004 * jcent; 
        break;
    case 7 : 
        ra0  = 257.311;
        dec0 = -15.175; 
        break;
    case 8 : 
        ra0  = 299.36 + 0.70 * sin ((357.85 + 52.316 * jcent) * DEG2RAD);
        dec0 = 43.46 - 0.51 * cos ((357.85 + 52.316 * jcent) * DEG2RAD); 
        break;
    case 9 : 
        ra0 = 313.02;
        dec0 = 9.09; 
        break;
    case I_TITAN : 
        ra0 = 39.4827;
        dec0 = 83.4279;        
        break;
    }

    cr = cos (ra0 * DEG2RAD);
    sr = sin (ra0 * DEG2RAD);
    cd = cos (dec0 * DEG2RAD);
    sd = sin (dec0 * DEG2RAD);

    tes[0] = -sr;
    tes[1] = -cr * sd;
    tes[2] = cr * cd;
    tes[3] = cr;
    tes[4] = -sr * sd;
    tes[5] = sr * cd;
    tes[6] = 0;
    tes[7] = cd;
    tes[8] = sd;
    return;
}





/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* in2pa - 	from inertial to moon fixed (PA)
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void in2pa(double *jd, double *te)
{
    double lib[6] = {0}, tb1[9], tb2[9], tb3[9], tb32[9];
	int target, center;

	target = 15;
	center = 0;

//	DPLEPH(jd, &target, &center, lib);

	rotmatz (lib[0], tb1, 0);
	rotmatx (lib[1], tb2, 0);
	rotmatz (lib[2], tb3, 0);

	brmul(tb3, tb2, 3, 3, 3, tb32);
	brmul(tb32, tb1, 3, 3, 3, te);
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* mbf2cel - simulate doppler observable
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
short int mbf2cel (double *jd_tdb, double *te)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function rotates a vector from the moon body-fixed system 
	  to the celestial system.  

   REFERENCES:
      P. Kenneth Seidelmann et. al. (2007). Report of the IAU/IAGWorking 
	  Group on cartographic coordinates and rotational elements: 2006

   INPUT
   ARGUMENTS:
      jd_tdb[2] (double)
		 TDB Julian date.
         High-order part (jd_tdb[0]) & Low-order part (jd_tdb[0]).
      method (short int)
         Selection for method
            = 0 ... IAU report formulae
            = 1 ... NASA/JPL DE/LE ephemeris 
      ref_sys (short int)
         Reference system in which moon body-fixed system is given
            = 0 ... Mean Earth/polar axis (ME) system
            = 1 ... Principal Axis (PA) system
      derivation (short int)
         Seclection derivation of parameters
            = 0 ... No derivation, vecc is normal
            = 1 ... fisrt parameter derivation, vecc is derivation
			= 2 ... second 
			= 3 ... third
      vecm[3] (double)
         Position vector referred to moon body-fixed system

   OUTPUT
   ARGUMENTS:
      vecc[3] (double)
         Position vector referred to ICRF axes (celestial system)

   RETURNED
   VALUE:
      =  0  ... everything is ok.
      =  1  ... invalid value of 'ref_sys'
      =  2  ... invalid value of 'method'

   GLOBALS
   USED:

   FUNCTIONS
   CALLED:

   VER./DATE/
   PROGRAMMER:
      V1.0/03-10/ (SHAO).

   NOTES:

------------------------------------------------------------------------
*/
{
    short int error = 0;
    double tb1[9], tb2[9], tbt1[9], tbt2[9];

/*
IAU report formulae
*/
    me2pa(tb1);
    mt(tb1, 3, 3, tbt1);

    in2me(jd_tdb, tb2, 0);	
    mt(tb2, 3, 3, tbt2);
    brmul(tbt2,tbt1,3,3,3,te);	//ME2ICRF

    return error;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* in2me - from inertial to moon fixed (ME)
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void in2me (double *jd, double *te, short int derivation)
{
    double ra, dec, w, lib[3], d, T, tb1[9], tb2[9], tb3[9], tb32[9], 
        E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11, E12, E13;
	
    d = jd[0] - 2451545.0 + jd[1];
    T = d / 36525.0;

    E1  = 125.045 - 0.0529921 * d;	
    E2  = 250.089 - 0.1059842 * d;	
    E3  = 260.008 + 13.0120009 * d;
    E4  = 176.625 + 13.3407154 * d;	
    E5  = 357.529 + 0.9856003 * d;	
    E6  = 311.589 + 26.4057084 * d;
    E7  = 134.963 + 13.0649930 * d;	
    E8  = 276.617 + 0.3287146 * d;	
    E9  = 34.226 + 1.7484877 * d;
    E10 = 15.134 - 0.1589763 * d; 	
    E11 = 119.743 + 0.0036096 * d;	
    E12 = 239.961 + 0.1643573 * d;
    E13 = 25.053 + 12.9590088 * d;

    ra  = 269.9949 + 0.0031 * T - 3.8787 * sin (E1 * DEG2RAD) 
        - 0.1204 * sin (E2 * DEG2RAD) + 0.0700 * sin (E3 * DEG2RAD) 
        - 0.0172 * sin (E4 * DEG2RAD) + 0.0072 * sin (E6 * DEG2RAD) 
        - 0.0052 * sin (E10 * DEG2RAD) + 0.0043 * sin (E13 * DEG2RAD);
    dec = 66.5392 + 0.0130 * T + 1.5419 * cos (E1 * DEG2RAD) 
        + 0.0239 * cos (E2 * DEG2RAD) - 0.0278 * cos (E3 * DEG2RAD) 
        + 0.0068 * cos (E4 * DEG2RAD) - 0.0029 * cos (E6 * DEG2RAD)
        + 0.0009 * cos (E7 * DEG2RAD) + 0.0008 * cos (E10 * DEG2RAD) 
        - 0.0009 * cos (E13 * DEG2RAD);
    w   = 38.3213 + 13.17635815 * d - 1.4e-12 * d * d 
        + 3.5610 * sin (E1 * DEG2RAD) + 0.1208 * sin (E2 * DEG2RAD) 
        - 0.0642 * sin (E3 * DEG2RAD) + 0.0158 * sin (E4 * DEG2RAD)
        + 0.0252 * sin (E5 * DEG2RAD) - 0.0066 * sin (E6 * DEG2RAD) 
        - 0.0047 * sin (E7 * DEG2RAD) - 0.0046 * sin (E8 * DEG2RAD) 
        + 0.0028 * sin (E9 * DEG2RAD) + 0.0052 * sin (E10 * DEG2RAD)
        + 0.0040 * sin (E11 * DEG2RAD) + 0.0019 * sin (E12 * DEG2RAD) 
        - 0.0044 * sin (E13 * DEG2RAD);

    lib[0] = (90.0 + ra) * DEG2RAD;
    lib[1] = (90.0 - dec) * DEG2RAD;
    lib[2] = w * DEG2RAD;
	
    rotmatz (lib[0], tb1, 0);
    rotmatx (lib[1], tb2, 0);
    rotmatz (lib[2], tb3, 0);

    if (derivation == 1)
        rotmatz (lib[0], tb1, 1);
    if (derivation == 2)
        rotmatx (lib[1], tb2, 1);
    if (derivation == 3)
        rotmatz (lib[2], tb3, 1);
    
    brmul(tb3, tb2, 3, 3, 3, tb32);
    brmul(tb32, tb1, 3, 3, 3, te);
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* me2pa - simulate doppler observable
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void me2pa (double *te)
{
    double tb1[9], tb2[9], tb3[9], tb32[9];

    rotmatx ( 0.1462 * ASEC2RAD, tb1, 0);
    rotmaty (79.0768 * ASEC2RAD, tb2, 0);
    rotmatz (63.8986 * ASEC2RAD, tb3, 0);
    brmul(tb3, tb2, 3, 3, 3, tb32);
    brmul(tb32, tb1, 3, 3, 3, te);
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* rotmatx - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void rotmatx (double rad, double *matx, short int deri)
{
    double cosst, sinst;

    cosst = cos(rad);
    sinst = sin(rad);

    matx[0] = 1;
    matx[1] = 0;
    matx[2] = 0;
    matx[3] = 0;
    matx[4] = cosst;
    matx[5] = sinst;
    matx[6] = 0;
    matx[7] = -sinst;
    matx[8] = cosst;

    if (deri == 1)
    {
        matx[0] = 0;
        matx[1] = 0;
        matx[2] = 0;
        matx[3] = 0;
        matx[4] = -sinst;
        matx[5] = cosst;
        matx[6] = 0;
        matx[7] = -cosst;
        matx[8] = -sinst;
    }
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* rotmaty - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void rotmaty (double rad, double *maty, short int deri)
{
    double cosst, sinst;

    cosst = cos(rad);
    sinst = sin(rad);

    maty[0] = cosst;
    maty[1] = 0;
    maty[2] = -sinst;
    maty[3] = 0;
    maty[4] = 1;
    maty[5] = 0;
    maty[6] = sinst;
    maty[7] = 0;
    maty[8] = cosst;

    if (deri == 1)
    {
        maty[0] = -sinst;
        maty[1] = 0;
        maty[2] = -cosst;
        maty[3] = 0;
        maty[4] = 0;
        maty[5] = 0;
        maty[6] = cosst;
        maty[7] = 0;
        maty[8] = -sinst;
    }

}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* rotmatz - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
void rotmatz (double rad, double *matz, short int deri)
{
    double cosst, sinst;

    cosst = cos(rad);
    sinst = sin(rad);

    matz[0] = cosst;
    matz[1] = sinst;
    matz[2] = 0;
    matz[3] = -sinst;
    matz[4] = cosst;
    matz[5] = 0;
    matz[6] = 0;
    matz[7] = 0;
    matz[8] = 1;

    if (deri == 1)
    {
        matz[0] = -sinst;
        matz[1] = cosst;
        matz[2] = 0;
        matz[3] = -cosst;
        matz[4] = -sinst;
        matz[5] = 0;
        matz[6] = 0;
        matz[7] = 0;
        matz[8] = 0;
    }
}




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/****************************************************************************/
/*                                                                          */
/*		Functions for Runge-Kutta integrator                                */
/*                                                                          */
/*      Version:    2009-9-8                                                */
/*                                                                          */
/*      Copyright (c) 2009 shangkun@shao.ac.cn All Right Reserved           */
/*                                                                          */
/****************************************************************************/

/*
  Version: 2009-9-8 
  Version: 2009-9-13 integrate forwards & backwards
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double rkf78_auto (double h, double t, double *x, int dim, double err, 
              double (*fun)(double,double *,double *), int autoadjust)
/*
    purpose: auto-adjusted Runge-Kutta-Ful... integrator 
    input:  double h				integration step
            double t                integrate from t to t+h
            double *x               x(t)
            int dim	                dim(x)
            double err              tolerance of step control
            double (*fun)()			right(force) function
    output: double *x				x(t+h)
    return: h                       new step after adjustment
*/
{	
    int i, j, n, flag = 0;
    double *y, *k, *f, d = 0, tn;
    double a[13] = { 0, 2.0/27, 1.0/9, 1.0/6, 5.0/12, 1.0/2, 5.0/6, 1.0/6, 
        2.0/3, 1.0/3, 1.0, 0, 1.0 };
    double c[13] = { 0, 0, 0, 0, 0, 34.0/105, 9.0/35, 9.0/35, 9.0/280, 
        9.0/280, 0, 41.0/840, 41.0/840 };
    double b[13][12] = 
    {
        {0},
        {2.0/27},
        {1.0/36,1.0/12},
        {1.0/24,0,1.0/8},
        {5.0/12,0,-25.0/16,25.0/16},
        {1.0/20,0,0,1.0/4,1.0/5},
        {-25.0/108,0,0,125.0/108,-65.0/27,125.0/54},
        {31.0/300,0,0,0,61.0/225,-2.0/9,13.0/900},
        {2.0,0,0,-53.0/6,704.0/45,-107.0/9,67.0/90,3.0},
        {-91.0/108,0,0,23.0/108,-976.0/135,311.0/54,-19.0/60,17.0/6,-1.0/12},
        {2383.0/4100,0,0,-341.0/164,4496.0/1025,-301.0/82,2133.0/4100,
        45.0/82,45.0/164,18.0/41},
        {3.0/205,0,0,0,0,-6.0/41,-3.0/205,-3.0/41,3.0/41,6.0/41},
        {-1777.0/4100,0,0,-341.0/164,4496.0/1025,-289.0/82,2193.0/4100,
        51.0/82,33.0/164,12.0/41,0,1.0}
    };
    
    y = (double *) calloc (dim, sizeof(double));
    k = (double *) calloc (dim*13, sizeof(double));
    f = (double *) calloc (dim, sizeof(double));
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
    
    do
    {	
        for (i = 0; i <= 12; i++)
        {			
            tn = t + a[i] * h;
            for (n = 0; n <= dim - 1; n++)
            {	
                y[n] = x[n];
                for (j = 0; j <= i-1; j++)
                    y[n] = y[n] + h * b[i][j] * k[n*13+j];
            }
            fun (tn,y,f);
            for (n = 0; n <= dim - 1; n++)
            {
                k[n*13+i] = f[n];
            }
        }
        d = 0;
        for (n = 0; n <= dim - 1; n++)
        {
            d = d + fabs (41.0 / 840 * (k[n*13+0] + k[n*13+10] 
                - k[n*13+11] - k[n*13+12]) * h);
        }
		
        flag = 0;
        if (autoadjust == 1)
        {
            if (d > err) //adapting step h
            {
                h = h/2.0; 
                flag = 1; 	
            }
            if ( (d < err * 1e-4) && (h < 5e-3))
            {
                h = h*2.0; 
                flag = 2;
            }
        }
    }while (flag == 1);

    for (n = 0; n <= dim - 1; n++)
    {
        for (i = 0; i <= 12; i++)
            x[n] = x[n] + h * c[i] * k[n*13+i];	
    }
	
    free (y);
    free (f);
    free (k);
    return h;
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* enlgr - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double enlgr (double *x, double *y, int n, double t)
{ 
    int i, j, k, m;
    double z, s;
    
    z = 0.0;
    if (n < 1) 
        return (z);
    if (n == 1) 
    { 
        z = y[0];
        return (z);
    }
    if (n == 2)
    { 
        z = (y[0] * (t - x[1]) - y[1] * (t - x[0])) / (x[0] - x[1]);
        return(z);
    }
    i = 0;
    while ((x[i] < t) && (i < n)) 
        i = i + 1;
    k = i - 4;
    if (k < 0) 
        k = 0;
    m = i + 3;
    if (m > n - 1) 
        m = n - 1;
    for (i = k; i <= m; i++)
    { 
        s = 1.0;
        for (j = k; j <= m; j++)
        {
            if (j != i) 
                s = s * (t - x[j]) / (x[i] - x[j]);
        }
        z = z + s * y[i];
    }
    return (z);
}



/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* bssgj - 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
int bssgj (double *a,int n)
{
    int i, j, k, m;
    double w, g, *b;

    b = (double *)malloc (n * sizeof(double));
    for (k = 0; k <= n - 1; k++)
    {
        w = a[0];
        if (fabs (w) + 1.0 == 1.0)
        { 
            free (b); 
            printf ("fail\n"); 
            return (-2);
        }
        m = n - k - 1;
        for (i = 1; i <= n - 1; i++)
        { 
            g = a[i * n]; 
            b[i] = g / w;
            if (i <= m) 
                b[i] = - b[i];
            for (j = 1; j <= i; j++)
                a[(i - 1) * n + j - 1] = a[i * n + j] + g * b[j];
        }
        a[n * n - 1] = 1.0 / w;
        for (i = 1; i <= n - 1; i++)
            a[(n - 1) * n + i - 1] = b[i];
    }
	
    for (i=0; i<=n-2; i++)
    {
        for (j=i+1; j<=n-1; j++)
            a[i*n+j]=a[j*n+i];
    }
     
    free(b);
    return(2);
}


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*
* iauDtdb - precise tdb-tt correction: better than +/- 3 nanoseconds 
* @param1: description of param1
* @param2: description of param2
*
* version: 20 Aug 2010
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
double iauDtdb(double date1, double date2,
               double ut, double elong, double u, double v)
/*
**  - - - - - - - -
**   i a u D t d b
**  - - - - - - - -
**
**  An approximation to TDB-TT, the difference between barycentric
**  dynamical time and terrestrial time, for an observer on the Earth.
**
**  The different time scales - proper, coordinate and realized - are
**  related to each other:
**
**            TAI             <-  physically realized
**             :
**          offset            <-  observed (nominally +32.184s)
**             :
**            TT              <-  terrestrial time
**             :
**    rate adjustment (L_G)   <-  definition of TT
**             :
**            TCG             <-  time scale for GCRS
**             :
**      "periodic" terms      <-  iauDtdb  is an implementation
**             :
**    rate adjustment (L_C)   <-  function of solar-system ephemeris
**             :
**            TCB             <-  time scale for BCRS
**             :
**    rate adjustment (-L_B)  <-  definition of TDB
**             :
**            TDB             <-  TCB scaled to track TT
**             :
**      "periodic" terms      <-  -iau_DTDB is an approximation
**             :
**            TT              <-  terrestrial time
**
**  Adopted values for the various constants can be found in the IERS
**  Conventions (McCarthy & Petit 2003).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2   double  date, TDB (Notes 1-3)
**     ut            double  universal time (UT1, fraction of one day)
**     elong         double  longitude (east positive, radians)
**     u             double  distance from Earth spin axis (km)
**     v             double  distance north of equatorial plane (km)
**
**  Returned (function value):
**                   double  TDB-TT (seconds)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**     Although the date is, formally, barycentric dynamical time (TDB),
**     the terrestrial dynamical time (TT) can be used with no practical
**     effect on the accuracy of the prediction.
**
**  2) TT can be regarded as a coordinate time that is realized as an
**     offset of 32.184s from International Atomic Time, TAI.  TT is a
**     specific linear transformation of geocentric coordinate time TCG,
**     which is the time scale for the Geocentric Celestial Reference
**     System, GCRS.
**
**  3) TDB is a coordinate time, and is a specific linear transformation
**     of barycentric coordinate time TCB, which is the time scale for
**     the Barycentric Celestial Reference System, BCRS.
**
**  4) The difference TCG-TCB depends on the masses and positions of the
**     bodies of the solar system and the velocity of the Earth.  It is
**     dominated by a rate difference, the residual being of a periodic
**     character.  The latter, which is modeled by the present function,
**     comprises a main (annual) sinusoidal term of amplitude
**     approximately 0.00166 seconds, plus planetary terms up to about
**     20 microseconds, and lunar and diurnal terms up to 2 microseconds.
**     These effects come from the changing transverse Doppler effect
**     and gravitational red-shift as the observer (on the Earth's
**     surface) experiences variations in speed (with respect to the
**     BCRS) and gravitational potential.
**
**  5) TDB can be regarded as the same as TCB but with a rate adjustment
**     to keep it close to TT, which is convenient for many applications.
**     The history of successive attempts to define TDB is set out in
**     Resolution 3 adopted by the IAU General Assembly in 2006, which
**     defines a fixed TDB(TCB) transformation that is consistent with
**     contemporary solar-system ephemerides.  Future ephemerides will
**     imply slightly changed transformations between TCG and TCB, which
**     could introduce a linear drift between TDB and TT;  however, any
**     such drift is unlikely to exceed 1 nanosecond per century.
**
**  6) The geocentric TDB-TT model used in the present function is that of
**     Fairhead & Bretagnon (1990), in its full form.  It was originally
**     supplied by Fairhead (private communications with P.T.Wallace,
**     1990) as a Fortran subroutine.  The present C function contains an
**     adaptation of the Fairhead code.  The numerical results are
**     essentially unaffected by the changes, the differences with
**     respect to the Fairhead & Bretagnon original being at the 1e-20 s
**     level.
**
**     The topocentric part of the model is from Moyer (1981) and
**     Murray (1983), with fundamental arguments adapted from
**     Simon et al. 1994.  It is an approximation to the expression
**     ( v / c ) . ( r / c ), where v is the barycentric velocity of
**     the Earth, r is the geocentric position of the observer and
**     c is the speed of light.
**
**     By supplying zeroes for u and v, the topocentric part of the
**     model can be nullified, and the function will return the Fairhead
**     & Bretagnon result alone.
**
**  7) During the interval 1950-2050, the absolute accuracy is better
**     than +/- 3 nanoseconds relative to time ephemerides obtained by
**     direct numerical integrations based on the JPL DE405 solar system
**     ephemeris.
**
**  8) It must be stressed that the present function is merely a model,
**     and that numerical integration of solar-system ephemerides is the
**     definitive method for predicting the relationship between TCG and
**     TCB and hence between TT and TDB.
**
**  References:
**
**     Fairhead, L., & Bretagnon, P., Astron.Astrophys., 229, 240-247
**     (1990).
**
**     IAU 2006 Resolution 3.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Moyer, T.D., Cel.Mech., 23, 33 (1981).
**
**     Murray, C.A., Vectorial Astrometry, Adam Hilger (1983).
**
**     Seidelmann, P.K. et al., Explanatory Supplement to the
**     Astronomical Almanac, Chapter 2, University Science Books (1992).
**
**     Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G. & Laskar, J., Astron.Astrophys., 282, 663-683 (1994).
**
**  This revision:  2008 May 24
**
**  Copyright (C) 2008 IAU SOFA Review Board.  See notes at end.
*/
{
   double t, tsol, w, elsun, emsun, d, elj, els, wt, w0, w1, w2, w3, w4,
          wf, wj;
   int j;
 
/*
** =====================
** Fairhead et al. model
** =====================
**
** 787 sets of three coefficients.
**
** Each set is
**    amplitude (microseconds)
**      frequency (radians per Julian millennium since J2000)
**      phase (radians)
**
** Sets   1-474 are the T**0 terms
**  "   475-679  "   "  T**1
**  "   680-764  "   "  T**2
**  "   765-784  "   "  T**3
**  "   785-787  "   "  T**4
*/
 
   static const double fairhd[787][3] = {
   /* 1, 10 */
      { 1656.674564e-6,     6283.075849991,  6.240054195 },
      {   22.417471e-6,     5753.384884897,  4.296977442 },
      {   13.839792e-6,    12566.151699983,  6.196904410 },
      {    4.770086e-6,      529.690965095,  0.444401603 },
      {    4.676740e-6,     6069.776754553,  4.021195093 },
      {    2.256707e-6,      213.299095438,  5.543113262 },
      {    1.694205e-6,      -3.523118349,   5.025132748 },
      {    1.554905e-6,    77713.771467920,  5.198467090 },
      {    1.276839e-6,     7860.419392439,  5.988822341 },
      {    1.193379e-6,     5223.693919802,  3.649823730 },
   /* 11, 20 */
      {    1.115322e-6,     3930.209696220,  1.422745069 },
      {    0.794185e-6,    11506.769769794,  2.322313077 },
      {    0.447061e-6,       26.298319800,  3.615796498 },
      {    0.435206e-6,     -398.149003408,  4.349338347 },
      {    0.600309e-6,     1577.343542448,  2.678271909 },
      {    0.496817e-6,     6208.294251424,  5.696701824 },
      {    0.486306e-6,     5884.926846583,  0.520007179 },
      {    0.432392e-6,       74.781598567,  2.435898309 },
      {    0.468597e-6,     6244.942814354,  5.866398759 },
      {    0.375510e-6,     5507.553238667,  4.103476804 },
   /* 21, 30 */
      {    0.243085e-6,     -775.522611324,  3.651837925 },
      {    0.173435e-6,    18849.227549974,  6.153743485 },
      {    0.230685e-6,     5856.477659115,  4.773852582 },
      {    0.203747e-6,    12036.460734888,  4.333987818 },
      {    0.143935e-6,     -796.298006816,  5.957517795 },
      {    0.159080e-6,    10977.078804699,  1.890075226 },
      {    0.119979e-6,       38.133035638,  4.551585768 },
      {    0.118971e-6,     5486.777843175,  1.914547226 },
      {    0.116120e-6,     1059.381930189,  0.873504123 },
      {    0.137927e-6,    11790.629088659,  1.135934669 },
   /* 31, 40 */
      {    0.098358e-6,     2544.314419883,  0.092793886 },
      {    0.101868e-6,    -5573.142801634,  5.984503847 },
      {    0.080164e-6,      206.185548437,  2.095377709 },
      {    0.079645e-6,     4694.002954708,  2.949233637 },
      {    0.062617e-6,       20.775395492,  2.654394814 },
      {    0.075019e-6,     2942.463423292,  4.980931759 },
      {    0.064397e-6,     5746.271337896,  1.280308748 },
      {    0.063814e-6,     5760.498431898,  4.167901731 },
      {    0.048042e-6,     2146.165416475,  1.495846011 },
      {    0.048373e-6,      155.420399434,  2.251573730 },
   /* 41, 50 */
      {    0.058844e-6,      426.598190876,  4.839650148 },
      {    0.046551e-6,       -0.980321068,  0.921573539 },
      {    0.054139e-6,    17260.154654690,  3.411091093 },
      {    0.042411e-6,     6275.962302991,  2.869567043 },
      {    0.040184e-6,       -7.113547001,  3.565975565 },
      {    0.036564e-6,     5088.628839767,  3.324679049 },
      {    0.040759e-6,    12352.852604545,  3.981496998 },
      {    0.036507e-6,      801.820931124,  6.248866009 },
      {    0.036955e-6,     3154.687084896,  5.071801441 },
      {    0.042732e-6,      632.783739313,  5.720622217 },
   /* 51, 60 */
      {    0.042560e-6,   161000.685737473,  1.270837679 },
      {    0.040480e-6,    15720.838784878,  2.546610123 },
      {    0.028244e-6,    -6286.598968340,  5.069663519 },
      {    0.033477e-6,     6062.663207553,  4.144987272 },
      {    0.034867e-6,      522.577418094,  5.210064075 },
      {    0.032438e-6,     6076.890301554,  0.749317412 },
      {    0.030215e-6,     7084.896781115,  3.389610345 },
      {    0.029247e-6,   -71430.695617928,  4.183178762 },
      {    0.033529e-6,     9437.762934887,  2.404714239 },
      {    0.032423e-6,     8827.390269875,  5.541473556 },
   /* 61, 70 */
      {    0.027567e-6,     6279.552731642,  5.040846034 },
      {    0.029862e-6,    12139.553509107,  1.770181024 },
      {    0.022509e-6,    10447.387839604,  1.460726241 },
      {    0.020937e-6,     8429.241266467,  0.652303414 },
      {    0.020322e-6,      419.484643875,  3.735430632 },
      {    0.024816e-6,    -1194.447010225,  1.087136918 },
      {    0.025196e-6,     1748.016413067,  2.901883301 },
      {    0.021691e-6,    14143.495242431,  5.952658009 },
      {    0.017673e-6,     6812.766815086,  3.186129845 },
      {    0.022567e-6,     6133.512652857,  3.307984806 },
   /* 71, 80 */
      {    0.016155e-6,    10213.285546211,  1.331103168 },
      {    0.014751e-6,     1349.867409659,  4.308933301 },
      {    0.015949e-6,     -220.412642439,  4.005298270 },
      {    0.015974e-6,    -2352.866153772,  6.145309371 },
      {    0.014223e-6,    17789.845619785,  2.104551349 },
      {    0.017806e-6,       73.297125859,  3.475975097 },
      {    0.013671e-6,     -536.804512095,  5.971672571 },
      {    0.011942e-6,     8031.092263058,  2.053414715 },
      {    0.014318e-6,    16730.463689596,  3.016058075 },
      {    0.012462e-6,      103.092774219,  1.737438797 },
   /* 81, 90 */
      {    0.010962e-6,        3.590428652,  2.196567739 },
      {    0.015078e-6,    19651.048481098,  3.969480770 },
      {    0.010396e-6,      951.718406251,  5.717799605 },
      {    0.011707e-6,    -4705.732307544,  2.654125618 },
      {    0.010453e-6,     5863.591206116,  1.913704550 },
      {    0.012420e-6,     4690.479836359,  4.734090399 },
      {    0.011847e-6,     5643.178563677,  5.489005403 },
      {    0.008610e-6,     3340.612426700,  3.661698944 },
      {    0.011622e-6,     5120.601145584,  4.863931876 },
      {    0.010825e-6,      553.569402842,  0.842715011 },
   /* 91, 100 */
      {    0.008666e-6,     -135.065080035,  3.293406547 },
      {    0.009963e-6,      149.563197135,  4.870690598 },
      {    0.009858e-6,     6309.374169791,  1.061816410 },
      {    0.007959e-6,      316.391869657,  2.465042647 },
      {    0.010099e-6,      283.859318865,  1.942176992 },
      {    0.007147e-6,     -242.728603974,  3.661486981 },
      {    0.007505e-6,     5230.807466803,  4.920937029 },
      {    0.008323e-6,    11769.853693166,  1.229392026 },
      {    0.007490e-6,    -6256.777530192,  3.658444681 },
      {    0.009370e-6,   149854.400134205,  0.673880395 },
   /* 101, 110 */
      {    0.007117e-6,       38.027672636,  5.294249518 },
      {    0.007857e-6,    12168.002696575,  0.525733528 },
      {    0.007019e-6,     6206.809778716,  0.837688810 },
      {    0.006056e-6,      955.599741609,  4.194535082 },
      {    0.008107e-6,    13367.972631107,  3.793235253 },
      {    0.006731e-6,     5650.292110678,  5.639906583 },
      {    0.007332e-6,       36.648562930,  0.114858677 },
      {    0.006366e-6,     4164.311989613,  2.262081818 },
      {    0.006858e-6,     5216.580372801,  0.642063318 },
      {    0.006919e-6,     6681.224853400,  6.018501522 },
   /* 111, 120 */
      {    0.006826e-6,     7632.943259650,  3.458654112 },
      {    0.005308e-6,    -1592.596013633,  2.500382359 },
      {    0.005096e-6,    11371.704689758,  2.547107806 },
      {    0.004841e-6,     5333.900241022,  0.437078094 },
      {    0.005582e-6,     5966.683980335,  2.246174308 },
      {    0.006304e-6,    11926.254413669,  2.512929171 },
      {    0.006603e-6,    23581.258177318,  5.393136889 },
      {    0.005123e-6,       -1.484472708,  2.999641028 },
      {    0.004648e-6,     1589.072895284,  1.275847090 },
      {    0.005119e-6,     6438.496249426,  1.486539246 },
   /* 121, 130 */
      {    0.004521e-6,     4292.330832950,  6.140635794 },
      {    0.005680e-6,    23013.539539587,  4.557814849 },
      {    0.005488e-6,       -3.455808046,  0.090675389 },
      {    0.004193e-6,     7234.794256242,  4.869091389 },
      {    0.003742e-6,     7238.675591600,  4.691976180 },
      {    0.004148e-6,     -110.206321219,  3.016173439 },
      {    0.004553e-6,    11499.656222793,  5.554998314 },
      {    0.004892e-6,     5436.993015240,  1.475415597 },
      {    0.004044e-6,     4732.030627343,  1.398784824 },
      {    0.004164e-6,    12491.370101415,  5.650931916 },
   /* 131, 140 */
      {    0.004349e-6,    11513.883316794,  2.181745369 },
      {    0.003919e-6,    12528.018664345,  5.823319737 },
      {    0.003129e-6,     6836.645252834,  0.003844094 },
      {    0.004080e-6,    -7058.598461315,  3.690360123 },
      {    0.003270e-6,       76.266071276,  1.517189902 },
      {    0.002954e-6,     6283.143160294,  4.447203799 },
      {    0.002872e-6,       28.449187468,  1.158692983 },
      {    0.002881e-6,      735.876513532,  0.349250250 },
      {    0.003279e-6,     5849.364112115,  4.893384368 },
      {    0.003625e-6,     6209.778724132,  1.473760578 },
   /* 141, 150 */
      {    0.003074e-6,      949.175608970,  5.185878737 },
      {    0.002775e-6,     9917.696874510,  1.030026325 },
      {    0.002646e-6,    10973.555686350,  3.918259169 },
      {    0.002575e-6,    25132.303399966,  6.109659023 },
      {    0.003500e-6,      263.083923373,  1.892100742 },
      {    0.002740e-6,    18319.536584880,  4.320519510 },
      {    0.002464e-6,      202.253395174,  4.698203059 },
      {    0.002409e-6,        2.542797281,  5.325009315 },
      {    0.003354e-6,   -90955.551694697,  1.942656623 },
      {    0.002296e-6,     6496.374945429,  5.061810696 },
   /* 151, 160 */
      {    0.003002e-6,     6172.869528772,  2.797822767 },
      {    0.003202e-6,    27511.467873537,  0.531673101 },
      {    0.002954e-6,    -6283.008539689,  4.533471191 },
      {    0.002353e-6,      639.897286314,  3.734548088 },
      {    0.002401e-6,    16200.772724501,  2.605547070 },
      {    0.003053e-6,   233141.314403759,  3.029030662 },
      {    0.003024e-6,    83286.914269554,  2.355556099 },
      {    0.002863e-6,    17298.182327326,  5.240963796 },
      {    0.002103e-6,    -7079.373856808,  5.756641637 },
      {    0.002303e-6,    83996.847317911,  2.013686814 },
   /* 161, 170 */
      {    0.002303e-6,    18073.704938650,  1.089100410 },
      {    0.002381e-6,       63.735898303,  0.759188178 },
      {    0.002493e-6,     6386.168624210,  0.645026535 },
      {    0.002366e-6,        3.932153263,  6.215885448 },
      {    0.002169e-6,    11015.106477335,  4.845297676 },
      {    0.002397e-6,     6243.458341645,  3.809290043 },
      {    0.002183e-6,     1162.474704408,  6.179611691 },
      {    0.002353e-6,     6246.427287062,  4.781719760 },
      {    0.002199e-6,     -245.831646229,  5.956152284 },
      {    0.001729e-6,     3894.181829542,  1.264976635 },
   /* 171, 180 */
      {    0.001896e-6,    -3128.388765096,  4.914231596 },
      {    0.002085e-6,       35.164090221,  1.405158503 },
      {    0.002024e-6,    14712.317116458,  2.752035928 },
      {    0.001737e-6,     6290.189396992,  5.280820144 },
      {    0.002229e-6,      491.557929457,  1.571007057 },
      {    0.001602e-6,    14314.168113050,  4.203664806 },
      {    0.002186e-6,      454.909366527,  1.402101526 },
      {    0.001897e-6,    22483.848574493,  4.167932508 },
      {    0.001825e-6,    -3738.761430108,  0.545828785 },
      {    0.001894e-6,     1052.268383188,  5.817167450 },
   /* 181, 190 */
      {    0.001421e-6,       20.355319399,  2.419886601 },
      {    0.001408e-6,    10984.192351700,  2.732084787 },
      {    0.001847e-6,    10873.986030480,  2.903477885 },
      {    0.001391e-6,    -8635.942003763,  0.593891500 },
      {    0.001388e-6,       -7.046236698,  1.166145902 },
      {    0.001810e-6,   -88860.057071188,  0.487355242 },
      {    0.001288e-6,    -1990.745017041,  3.913022880 },
      {    0.001297e-6,    23543.230504682,  3.063805171 },
      {    0.001335e-6,     -266.607041722,  3.995764039 },
      {    0.001376e-6,    10969.965257698,  5.152914309 },
   /* 191, 200 */
      {    0.001745e-6,   244287.600007027,  3.626395673 },
      {    0.001649e-6,    31441.677569757,  1.952049260 },
      {    0.001416e-6,     9225.539273283,  4.996408389 },
      {    0.001238e-6,     4804.209275927,  5.503379738 },
      {    0.001472e-6,     4590.910180489,  4.164913291 },
      {    0.001169e-6,     6040.347246017,  5.841719038 },
      {    0.001039e-6,     5540.085789459,  2.769753519 },
      {    0.001004e-6,     -170.672870619,  0.755008103 },
      {    0.001284e-6,    10575.406682942,  5.306538209 },
      {    0.001278e-6,       71.812653151,  4.713486491 },
   /* 201, 210 */
      {    0.001321e-6,    18209.330263660,  2.624866359 },
      {    0.001297e-6,    21228.392023546,  0.382603541 },
      {    0.000954e-6,     6282.095528923,  0.882213514 },
      {    0.001145e-6,     6058.731054289,  1.169483931 },
      {    0.000979e-6,     5547.199336460,  5.448375984 },
      {    0.000987e-6,    -6262.300454499,  2.656486959 },
      {    0.001070e-6,  -154717.609887482,  1.827624012 },
      {    0.000991e-6,     4701.116501708,  4.387001801 },
      {    0.001155e-6,      -14.227094002,  3.042700750 },
      {    0.001176e-6,      277.034993741,  3.335519004 },
   /* 211, 220 */
      {    0.000890e-6,    13916.019109642,  5.601498297 },
      {    0.000884e-6,    -1551.045222648,  1.088831705 },
      {    0.000876e-6,     5017.508371365,  3.969902609 },
      {    0.000806e-6,    15110.466119866,  5.142876744 },
      {    0.000773e-6,    -4136.910433516,  0.022067765 },
      {    0.001077e-6,      175.166059800,  1.844913056 },
      {    0.000954e-6,    -6284.056171060,  0.968480906 },
      {    0.000737e-6,     5326.786694021,  4.923831588 },
      {    0.000845e-6,     -433.711737877,  4.749245231 },
      {    0.000819e-6,     8662.240323563,  5.991247817 },
   /* 221, 230 */
      {    0.000852e-6,      199.072001436,  2.189604979 },
      {    0.000723e-6,    17256.631536341,  6.068719637 },
      {    0.000940e-6,     6037.244203762,  6.197428148 },
      {    0.000885e-6,    11712.955318231,  3.280414875 },
      {    0.000706e-6,    12559.038152982,  2.824848947 },
      {    0.000732e-6,     2379.164473572,  2.501813417 },
      {    0.000764e-6,    -6127.655450557,  2.236346329 },
      {    0.000908e-6,      131.541961686,  2.521257490 },
      {    0.000907e-6,    35371.887265976,  3.370195967 },
      {    0.000673e-6,     1066.495477190,  3.876512374 },
   /* 231, 240 */
      {    0.000814e-6,    17654.780539750,  4.627122566 },
      {    0.000630e-6,       36.027866677,  0.156368499 },
      {    0.000798e-6,      515.463871093,  5.151962502 },
      {    0.000798e-6,      148.078724426,  5.909225055 },
      {    0.000806e-6,      309.278322656,  6.054064447 },
      {    0.000607e-6,      -39.617508346,  2.839021623 },
      {    0.000601e-6,      412.371096874,  3.984225404 },
      {    0.000646e-6,    11403.676995575,  3.852959484 },
      {    0.000704e-6,    13521.751441591,  2.300991267 },
      {    0.000603e-6,   -65147.619767937,  4.140083146 },
   /* 241, 250 */
      {    0.000609e-6,    10177.257679534,  0.437122327 },
      {    0.000631e-6,     5767.611978898,  4.026532329 },
      {    0.000576e-6,    11087.285125918,  4.760293101 },
      {    0.000674e-6,    14945.316173554,  6.270510511 },
      {    0.000726e-6,     5429.879468239,  6.039606892 },
      {    0.000710e-6,    28766.924424484,  5.672617711 },
      {    0.000647e-6,    11856.218651625,  3.397132627 },
      {    0.000678e-6,    -5481.254918868,  6.249666675 },
      {    0.000618e-6,    22003.914634870,  2.466427018 },
      {    0.000738e-6,     6134.997125565,  2.242668890 },
   /* 251, 260 */
      {    0.000660e-6,      625.670192312,  5.864091907 },
      {    0.000694e-6,     3496.032826134,  2.668309141 },
      {    0.000531e-6,     6489.261398429,  1.681888780 },
      {    0.000611e-6,  -143571.324284214,  2.424978312 },
      {    0.000575e-6,    12043.574281889,  4.216492400 },
      {    0.000553e-6,    12416.588502848,  4.772158039 },
      {    0.000689e-6,     4686.889407707,  6.224271088 },
      {    0.000495e-6,     7342.457780181,  3.817285811 },
      {    0.000567e-6,     3634.621024518,  1.649264690 },
      {    0.000515e-6,    18635.928454536,  3.945345892 },
   /* 261, 270 */
      {    0.000486e-6,     -323.505416657,  4.061673868 },
      {    0.000662e-6,    25158.601719765,  1.794058369 },
      {    0.000509e-6,      846.082834751,  3.053874588 },
      {    0.000472e-6,   -12569.674818332,  5.112133338 },
      {    0.000461e-6,     6179.983075773,  0.513669325 },
      {    0.000641e-6,    83467.156352816,  3.210727723 },
      {    0.000520e-6,    10344.295065386,  2.445597761 },
      {    0.000493e-6,    18422.629359098,  1.676939306 },
      {    0.000478e-6,     1265.567478626,  5.487314569 },
      {    0.000472e-6,      -18.159247265,  1.999707589 },
   /* 271, 280 */
      {    0.000559e-6,    11190.377900137,  5.783236356 },
      {    0.000494e-6,     9623.688276691,  3.022645053 },
      {    0.000463e-6,     5739.157790895,  1.411223013 },
      {    0.000432e-6,    16858.482532933,  1.179256434 },
      {    0.000574e-6,    72140.628666286,  1.758191830 },
      {    0.000484e-6,    17267.268201691,  3.290589143 },
      {    0.000550e-6,     4907.302050146,  0.864024298 },
      {    0.000399e-6,       14.977853527,  2.094441910 },
      {    0.000491e-6,      224.344795702,  0.878372791 },
      {    0.000432e-6,    20426.571092422,  6.003829241 },
   /* 281, 290 */
      {    0.000481e-6,     5749.452731634,  4.309591964 },
      {    0.000480e-6,     5757.317038160,  1.142348571 },
      {    0.000485e-6,     6702.560493867,  0.210580917 },
      {    0.000426e-6,     6055.549660552,  4.274476529 },
      {    0.000480e-6,     5959.570433334,  5.031351030 },
      {    0.000466e-6,    12562.628581634,  4.959581597 },
      {    0.000520e-6,    39302.096962196,  4.788002889 },
      {    0.000458e-6,    12132.439962106,  1.880103788 },
      {    0.000470e-6,    12029.347187887,  1.405611197 },
      {    0.000416e-6,    -7477.522860216,  1.082356330 },
   /* 291, 300 */
      {    0.000449e-6,    11609.862544012,  4.179989585 },
      {    0.000465e-6,    17253.041107690,  0.353496295 },
      {    0.000362e-6,    -4535.059436924,  1.583849576 },
      {    0.000383e-6,    21954.157609398,  3.747376371 },
      {    0.000389e-6,       17.252277143,  1.395753179 },
      {    0.000331e-6,    18052.929543158,  0.566790582 },
      {    0.000430e-6,    13517.870106233,  0.685827538 },
      {    0.000368e-6,    -5756.908003246,  0.731374317 },
      {    0.000330e-6,    10557.594160824,  3.710043680 },
      {    0.000332e-6,    20199.094959633,  1.652901407 },
   /* 301, 310 */
      {    0.000384e-6,    11933.367960670,  5.827781531 },
      {    0.000387e-6,    10454.501386605,  2.541182564 },
      {    0.000325e-6,    15671.081759407,  2.178850542 },
      {    0.000318e-6,      138.517496871,  2.253253037 },
      {    0.000305e-6,     9388.005909415,  0.578340206 },
      {    0.000352e-6,     5749.861766548,  3.000297967 },
      {    0.000311e-6,     6915.859589305,  1.693574249 },
      {    0.000297e-6,    24072.921469776,  1.997249392 },
      {    0.000363e-6,     -640.877607382,  5.071820966 },
      {    0.000323e-6,    12592.450019783,  1.072262823 },
   /* 311, 320 */
      {    0.000341e-6,    12146.667056108,  4.700657997 },
      {    0.000290e-6,     9779.108676125,  1.812320441 },
      {    0.000342e-6,     6132.028180148,  4.322238614 },
      {    0.000329e-6,     6268.848755990,  3.033827743 },
      {    0.000374e-6,    17996.031168222,  3.388716544 },
      {    0.000285e-6,     -533.214083444,  4.687313233 },
      {    0.000338e-6,     6065.844601290,  0.877776108 },
      {    0.000276e-6,       24.298513841,  0.770299429 },
      {    0.000336e-6,    -2388.894020449,  5.353796034 },
      {    0.000290e-6,     3097.883822726,  4.075291557 },
   /* 321, 330 */
      {    0.000318e-6,      709.933048357,  5.941207518 },
      {    0.000271e-6,    13095.842665077,  3.208912203 },
      {    0.000331e-6,     6073.708907816,  4.007881169 },
      {    0.000292e-6,      742.990060533,  2.714333592 },
      {    0.000362e-6,    29088.811415985,  3.215977013 },
      {    0.000280e-6,    12359.966151546,  0.710872502 },
      {    0.000267e-6,    10440.274292604,  4.730108488 },
      {    0.000262e-6,      838.969287750,  1.327720272 },
      {    0.000250e-6,    16496.361396202,  0.898769761 },
      {    0.000325e-6,    20597.243963041,  0.180044365 },
   /* 331, 340 */
      {    0.000268e-6,     6148.010769956,  5.152666276 },
      {    0.000284e-6,     5636.065016677,  5.655385808 },
      {    0.000301e-6,     6080.822454817,  2.135396205 },
      {    0.000294e-6,     -377.373607916,  3.708784168 },
      {    0.000236e-6,     2118.763860378,  1.733578756 },
      {    0.000234e-6,     5867.523359379,  5.575209112 },
      {    0.000268e-6,  -226858.238553767,  0.069432392 },
      {    0.000265e-6,   167283.761587465,  4.369302826 },
      {    0.000280e-6,    28237.233459389,  5.304829118 },
      {    0.000292e-6,    12345.739057544,  4.096094132 },
   /* 341, 350 */
      {    0.000223e-6,    19800.945956225,  3.069327406 },
      {    0.000301e-6,    43232.306658416,  6.205311188 },
      {    0.000264e-6,    18875.525869774,  1.417263408 },
      {    0.000304e-6,    -1823.175188677,  3.409035232 },
      {    0.000301e-6,      109.945688789,  0.510922054 },
      {    0.000260e-6,      813.550283960,  2.389438934 },
      {    0.000299e-6,   316428.228673312,  5.384595078 },
      {    0.000211e-6,     5756.566278634,  3.789392838 },
      {    0.000209e-6,     5750.203491159,  1.661943545 },
      {    0.000240e-6,    12489.885628707,  5.684549045 },
   /* 351, 360 */
      {    0.000216e-6,     6303.851245484,  3.862942261 },
      {    0.000203e-6,     1581.959348283,  5.549853589 },
      {    0.000200e-6,     5642.198242609,  1.016115785 },
      {    0.000197e-6,      -70.849445304,  4.690702525 },
      {    0.000227e-6,     6287.008003254,  2.911891613 },
      {    0.000197e-6,      533.623118358,  1.048982898 },
      {    0.000205e-6,    -6279.485421340,  1.829362730 },
      {    0.000209e-6,   -10988.808157535,  2.636140084 },
      {    0.000208e-6,     -227.526189440,  4.127883842 },
      {    0.000191e-6,      415.552490612,  4.401165650 },
   /* 361, 370 */
      {    0.000190e-6,    29296.615389579,  4.175658539 },
      {    0.000264e-6,    66567.485864652,  4.601102551 },
      {    0.000256e-6,    -3646.350377354,  0.506364778 },
      {    0.000188e-6,    13119.721102825,  2.032195842 },
      {    0.000185e-6,     -209.366942175,  4.694756586 },
      {    0.000198e-6,    25934.124331089,  3.832703118 },
      {    0.000195e-6,     4061.219215394,  3.308463427 },
      {    0.000234e-6,     5113.487598583,  1.716090661 },
      {    0.000188e-6,     1478.866574064,  5.686865780 },
      {    0.000222e-6,    11823.161639450,  1.942386641 },
   /* 371, 380 */
      {    0.000181e-6,    10770.893256262,  1.999482059 },
      {    0.000171e-6,     6546.159773364,  1.182807992 },
      {    0.000206e-6,       70.328180442,  5.934076062 },
      {    0.000169e-6,    20995.392966449,  2.169080622 },
      {    0.000191e-6,    10660.686935042,  5.405515999 },
      {    0.000228e-6,    33019.021112205,  4.656985514 },
      {    0.000184e-6,    -4933.208440333,  3.327476868 },
      {    0.000220e-6,     -135.625325010,  1.765430262 },
      {    0.000166e-6,    23141.558382925,  3.454132746 },
      {    0.000191e-6,     6144.558353121,  5.020393445 },
   /* 381, 390 */
      {    0.000180e-6,     6084.003848555,  0.602182191 },
      {    0.000163e-6,    17782.732072784,  4.960593133 },
      {    0.000225e-6,    16460.333529525,  2.596451817 },
      {    0.000222e-6,     5905.702242076,  3.731990323 },
      {    0.000204e-6,      227.476132789,  5.636192701 },
      {    0.000159e-6,    16737.577236597,  3.600691544 },
      {    0.000200e-6,     6805.653268085,  0.868220961 },
      {    0.000187e-6,    11919.140866668,  2.629456641 },
      {    0.000161e-6,      127.471796607,  2.862574720 },
      {    0.000205e-6,     6286.666278643,  1.742882331 },
   /* 391, 400 */
      {    0.000189e-6,      153.778810485,  4.812372643 },
      {    0.000168e-6,    16723.350142595,  0.027860588 },
      {    0.000149e-6,    11720.068865232,  0.659721876 },
      {    0.000189e-6,     5237.921013804,  5.245313000 },
      {    0.000143e-6,     6709.674040867,  4.317625647 },
      {    0.000146e-6,     4487.817406270,  4.815297007 },
      {    0.000144e-6,     -664.756045130,  5.381366880 },
      {    0.000175e-6,     5127.714692584,  4.728443327 },
      {    0.000162e-6,     6254.626662524,  1.435132069 },
      {    0.000187e-6,    47162.516354635,  1.354371923 },
   /* 401, 410 */
      {    0.000146e-6,    11080.171578918,  3.369695406 },
      {    0.000180e-6,     -348.924420448,  2.490902145 },
      {    0.000148e-6,      151.047669843,  3.799109588 },
      {    0.000157e-6,     6197.248551160,  1.284375887 },
      {    0.000167e-6,      146.594251718,  0.759969109 },
      {    0.000133e-6,    -5331.357443741,  5.409701889 },
      {    0.000154e-6,       95.979227218,  3.366890614 },
      {    0.000148e-6,    -6418.140930027,  3.384104996 },
      {    0.000128e-6,    -6525.804453965,  3.803419985 },
      {    0.000130e-6,    11293.470674356,  0.939039445 },
   /* 411, 420 */
      {    0.000152e-6,    -5729.506447149,  0.734117523 },
      {    0.000138e-6,      210.117701700,  2.564216078 },
      {    0.000123e-6,     6066.595360816,  4.517099537 },
      {    0.000140e-6,    18451.078546566,  0.642049130 },
      {    0.000126e-6,    11300.584221356,  3.485280663 },
      {    0.000119e-6,    10027.903195729,  3.217431161 },
      {    0.000151e-6,     4274.518310832,  4.404359108 },
      {    0.000117e-6,     6072.958148291,  0.366324650 },
      {    0.000165e-6,    -7668.637425143,  4.298212528 },
      {    0.000117e-6,    -6245.048177356,  5.379518958 },
   /* 421, 430 */
      {    0.000130e-6,    -5888.449964932,  4.527681115 },
      {    0.000121e-6,     -543.918059096,  6.109429504 },
      {    0.000162e-6,     9683.594581116,  5.720092446 },
      {    0.000141e-6,     6219.339951688,  0.679068671 },
      {    0.000118e-6,    22743.409379516,  4.881123092 },
      {    0.000129e-6,     1692.165669502,  0.351407289 },
      {    0.000126e-6,     5657.405657679,  5.146592349 },
      {    0.000114e-6,      728.762966531,  0.520791814 },
      {    0.000120e-6,       52.596639600,  0.948516300 },
      {    0.000115e-6,       65.220371012,  3.504914846 },
   /* 431, 440 */
      {    0.000126e-6,     5881.403728234,  5.577502482 },
      {    0.000158e-6,   163096.180360983,  2.957128968 },
      {    0.000134e-6,    12341.806904281,  2.598576764 },
      {    0.000151e-6,    16627.370915377,  3.985702050 },
      {    0.000109e-6,     1368.660252845,  0.014730471 },
      {    0.000131e-6,     6211.263196841,  0.085077024 },
      {    0.000146e-6,     5792.741760812,  0.708426604 },
      {    0.000146e-6,      -77.750543984,  3.121576600 },
      {    0.000107e-6,     5341.013788022,  0.288231904 },
      {    0.000138e-6,     6281.591377283,  2.797450317 },
   /* 441, 450 */
      {    0.000113e-6,    -6277.552925684,  2.788904128 },
      {    0.000115e-6,     -525.758811831,  5.895222200 },
      {    0.000138e-6,     6016.468808270,  6.096188999 },
      {    0.000139e-6,    23539.707386333,  2.028195445 },
      {    0.000146e-6,    -4176.041342449,  4.660008502 },
      {    0.000107e-6,    16062.184526117,  4.066520001 },
      {    0.000142e-6,    83783.548222473,  2.936315115 },
      {    0.000128e-6,     9380.959672717,  3.223844306 },
      {    0.000135e-6,     6205.325306007,  1.638054048 },
      {    0.000101e-6,     2699.734819318,  5.481603249 },
   /* 451, 460 */
      {    0.000104e-6,     -568.821874027,  2.205734493 },
      {    0.000103e-6,     6321.103522627,  2.440421099 },
      {    0.000119e-6,     6321.208885629,  2.547496264 },
      {    0.000138e-6,     1975.492545856,  2.314608466 },
      {    0.000121e-6,      137.033024162,  4.539108237 },
      {    0.000123e-6,    19402.796952817,  4.538074405 },
      {    0.000119e-6,    22805.735565994,  2.869040566 },
      {    0.000133e-6,    64471.991241142,  6.056405489 },
      {    0.000129e-6,      -85.827298831,  2.540635083 },
      {    0.000131e-6,    13613.804277336,  4.005732868 },
   /* 461, 470 */
      {    0.000104e-6,     9814.604100291,  1.959967212 },
      {    0.000112e-6,    16097.679950283,  3.589026260 },
      {    0.000123e-6,     2107.034507542,  1.728627253 },
      {    0.000121e-6,    36949.230808424,  6.072332087 },
      {    0.000108e-6,   -12539.853380183,  3.716133846 },
      {    0.000113e-6,    -7875.671863624,  2.725771122 },
      {    0.000109e-6,     4171.425536614,  4.033338079 },
      {    0.000101e-6,     6247.911759770,  3.441347021 },
      {    0.000113e-6,     7330.728427345,  0.656372122 },
      {    0.000113e-6,    51092.726050855,  2.791483066 },
   /* 471, 480 */
      {    0.000106e-6,     5621.842923210,  1.815323326 },
      {    0.000101e-6,      111.430161497,  5.711033677 },
      {    0.000103e-6,      909.818733055,  2.812745443 },
      {    0.000101e-6,     1790.642637886,  1.965746028 },
 
   /* T */
      {  102.156724e-6,     6283.075849991,  4.249032005 },
      {    1.706807e-6,    12566.151699983,  4.205904248 },
      {    0.269668e-6,      213.299095438,  3.400290479 },
      {    0.265919e-6,      529.690965095,  5.836047367 },
      {    0.210568e-6,       -3.523118349,  6.262738348 },
      {    0.077996e-6,     5223.693919802,  4.670344204 },
   /* 481, 490 */
      {    0.054764e-6,     1577.343542448,  4.534800170 },
      {    0.059146e-6,       26.298319800,  1.083044735 },
      {    0.034420e-6,     -398.149003408,  5.980077351 },
      {    0.032088e-6,    18849.227549974,  4.162913471 },
      {    0.033595e-6,     5507.553238667,  5.980162321 },
      {    0.029198e-6,     5856.477659115,  0.623811863 },
      {    0.027764e-6,      155.420399434,  3.745318113 },
      {    0.025190e-6,     5746.271337896,  2.980330535 },
      {    0.022997e-6,     -796.298006816,  1.174411803 },
      {    0.024976e-6,     5760.498431898,  2.467913690 },
   /* 491, 500 */
      {    0.021774e-6,      206.185548437,  3.854787540 },
      {    0.017925e-6,     -775.522611324,  1.092065955 },
      {    0.013794e-6,      426.598190876,  2.699831988 },
      {    0.013276e-6,     6062.663207553,  5.845801920 },
      {    0.011774e-6,    12036.460734888,  2.292832062 },
      {    0.012869e-6,     6076.890301554,  5.333425680 },
      {    0.012152e-6,     1059.381930189,  6.222874454 },
      {    0.011081e-6,       -7.113547001,  5.154724984 },
      {    0.010143e-6,     4694.002954708,  4.044013795 },
      {    0.009357e-6,     5486.777843175,  3.416081409 },
   /* 501, 510 */
      {    0.010084e-6,      522.577418094,  0.749320262 },
      {    0.008587e-6,    10977.078804699,  2.777152598 },
      {    0.008628e-6,     6275.962302991,  4.562060226 },
      {    0.008158e-6,     -220.412642439,  5.806891533 },
      {    0.007746e-6,     2544.314419883,  1.603197066 },
      {    0.007670e-6,     2146.165416475,  3.000200440 },
      {    0.007098e-6,       74.781598567,  0.443725817 },
      {    0.006180e-6,     -536.804512095,  1.302642751 },
      {    0.005818e-6,     5088.628839767,  4.827723531 },
      {    0.004945e-6,    -6286.598968340,  0.268305170 },
   /* 511, 520 */
      {    0.004774e-6,     1349.867409659,  5.808636673 },
      {    0.004687e-6,     -242.728603974,  5.154890570 },
      {    0.006089e-6,     1748.016413067,  4.403765209 },
      {    0.005975e-6,    -1194.447010225,  2.583472591 },
      {    0.004229e-6,      951.718406251,  0.931172179 },
      {    0.005264e-6,      553.569402842,  2.336107252 },
      {    0.003049e-6,     5643.178563677,  1.362634430 },
      {    0.002974e-6,     6812.766815086,  1.583012668 },
      {    0.003403e-6,    -2352.866153772,  2.552189886 },
      {    0.003030e-6,      419.484643875,  5.286473844 },
   /* 521, 530 */
      {    0.003210e-6,       -7.046236698,  1.863796539 },
      {    0.003058e-6,     9437.762934887,  4.226420633 },
      {    0.002589e-6,    12352.852604545,  1.991935820 },
      {    0.002927e-6,     5216.580372801,  2.319951253 },
      {    0.002425e-6,     5230.807466803,  3.084752833 },
      {    0.002656e-6,     3154.687084896,  2.487447866 },
      {    0.002445e-6,    10447.387839604,  2.347139160 },
      {    0.002990e-6,     4690.479836359,  6.235872050 },
      {    0.002890e-6,     5863.591206116,  0.095197563 },
      {    0.002498e-6,     6438.496249426,  2.994779800 },
   /* 531, 540 */
      {    0.001889e-6,     8031.092263058,  3.569003717 },
      {    0.002567e-6,      801.820931124,  3.425611498 },
      {    0.001803e-6,   -71430.695617928,  2.192295512 },
      {    0.001782e-6,        3.932153263,  5.180433689 },
      {    0.001694e-6,    -4705.732307544,  4.641779174 },
      {    0.001704e-6,    -1592.596013633,  3.997097652 },
      {    0.001735e-6,     5849.364112115,  0.417558428 },
      {    0.001643e-6,     8429.241266467,  2.180619584 },
      {    0.001680e-6,       38.133035638,  4.164529426 },
      {    0.002045e-6,     7084.896781115,  0.526323854 },
   /* 541, 550 */
      {    0.001458e-6,     4292.330832950,  1.356098141 },
      {    0.001437e-6,       20.355319399,  3.895439360 },
      {    0.001738e-6,     6279.552731642,  0.087484036 },
      {    0.001367e-6,    14143.495242431,  3.987576591 },
      {    0.001344e-6,     7234.794256242,  0.090454338 },
      {    0.001438e-6,    11499.656222793,  0.974387904 },
      {    0.001257e-6,     6836.645252834,  1.509069366 },
      {    0.001358e-6,    11513.883316794,  0.495572260 },
      {    0.001628e-6,     7632.943259650,  4.968445721 },
      {    0.001169e-6,      103.092774219,  2.838496795 },
   /* 551, 560 */
      {    0.001162e-6,     4164.311989613,  3.408387778 },
      {    0.001092e-6,     6069.776754553,  3.617942651 },
      {    0.001008e-6,    17789.845619785,  0.286350174 },
      {    0.001008e-6,      639.897286314,  1.610762073 },
      {    0.000918e-6,    10213.285546211,  5.532798067 },
      {    0.001011e-6,    -6256.777530192,  0.661826484 },
      {    0.000753e-6,    16730.463689596,  3.905030235 },
      {    0.000737e-6,    11926.254413669,  4.641956361 },
      {    0.000694e-6,     3340.612426700,  2.111120332 },
      {    0.000701e-6,     3894.181829542,  2.760823491 },
   /* 561, 570 */
      {    0.000689e-6,     -135.065080035,  4.768800780 },
      {    0.000700e-6,    13367.972631107,  5.760439898 },
      {    0.000664e-6,     6040.347246017,  1.051215840 },
      {    0.000654e-6,     5650.292110678,  4.911332503 },
      {    0.000788e-6,     6681.224853400,  4.699648011 },
      {    0.000628e-6,     5333.900241022,  5.024608847 },
      {    0.000755e-6,     -110.206321219,  4.370971253 },
      {    0.000628e-6,     6290.189396992,  3.660478857 },
      {    0.000635e-6,    25132.303399966,  4.121051532 },
      {    0.000534e-6,     5966.683980335,  1.173284524 },
   /* 571, 580 */
      {    0.000543e-6,     -433.711737877,  0.345585464 },
      {    0.000517e-6,    -1990.745017041,  5.414571768 },
      {    0.000504e-6,     5767.611978898,  2.328281115 },
      {    0.000485e-6,     5753.384884897,  1.685874771 },
      {    0.000463e-6,     7860.419392439,  5.297703006 },
      {    0.000604e-6,      515.463871093,  0.591998446 },
      {    0.000443e-6,    12168.002696575,  4.830881244 },
      {    0.000570e-6,      199.072001436,  3.899190272 },
      {    0.000465e-6,    10969.965257698,  0.476681802 },
      {    0.000424e-6,    -7079.373856808,  1.112242763 },
   /* 581, 590 */
      {    0.000427e-6,      735.876513532,  1.994214480 },
      {    0.000478e-6,    -6127.655450557,  3.778025483 },
      {    0.000414e-6,    10973.555686350,  5.441088327 },
      {    0.000512e-6,     1589.072895284,  0.107123853 },
      {    0.000378e-6,    10984.192351700,  0.915087231 },
      {    0.000402e-6,    11371.704689758,  4.107281715 },
      {    0.000453e-6,     9917.696874510,  1.917490952 },
      {    0.000395e-6,      149.563197135,  2.763124165 },
      {    0.000371e-6,     5739.157790895,  3.112111866 },
      {    0.000350e-6,    11790.629088659,  0.440639857 },
   /* 591, 600 */
      {    0.000356e-6,     6133.512652857,  5.444568842 },
      {    0.000344e-6,      412.371096874,  5.676832684 },
      {    0.000383e-6,      955.599741609,  5.559734846 },
      {    0.000333e-6,     6496.374945429,  0.261537984 },
      {    0.000340e-6,     6055.549660552,  5.975534987 },
      {    0.000334e-6,     1066.495477190,  2.335063907 },
      {    0.000399e-6,    11506.769769794,  5.321230910 },
      {    0.000314e-6,    18319.536584880,  2.313312404 },
      {    0.000424e-6,     1052.268383188,  1.211961766 },
      {    0.000307e-6,       63.735898303,  3.169551388 },
   /* 601, 610 */
      {    0.000329e-6,       29.821438149,  6.106912080 },
      {    0.000357e-6,     6309.374169791,  4.223760346 },
      {    0.000312e-6,    -3738.761430108,  2.180556645 },
      {    0.000301e-6,      309.278322656,  1.499984572 },
      {    0.000268e-6,    12043.574281889,  2.447520648 },
      {    0.000257e-6,    12491.370101415,  3.662331761 },
      {    0.000290e-6,      625.670192312,  1.272834584 },
      {    0.000256e-6,     5429.879468239,  1.913426912 },
      {    0.000339e-6,     3496.032826134,  4.165930011 },
      {    0.000283e-6,     3930.209696220,  4.325565754 },
   /* 611, 620 */
      {    0.000241e-6,    12528.018664345,  3.832324536 },
      {    0.000304e-6,     4686.889407707,  1.612348468 },
      {    0.000259e-6,    16200.772724501,  3.470173146 },
      {    0.000238e-6,    12139.553509107,  1.147977842 },
      {    0.000236e-6,     6172.869528772,  3.776271728 },
      {    0.000296e-6,    -7058.598461315,  0.460368852 },
      {    0.000306e-6,    10575.406682942,  0.554749016 },
      {    0.000251e-6,    17298.182327326,  0.834332510 },
      {    0.000290e-6,     4732.030627343,  4.759564091 },
      {    0.000261e-6,     5884.926846583,  0.298259862 },
   /* 621, 630 */
      {    0.000249e-6,     5547.199336460,  3.749366406 },
      {    0.000213e-6,    11712.955318231,  5.415666119 },
      {    0.000223e-6,     4701.116501708,  2.703203558 },
      {    0.000268e-6,     -640.877607382,  0.283670793 },
      {    0.000209e-6,     5636.065016677,  1.238477199 },
      {    0.000193e-6,    10177.257679534,  1.943251340 },
      {    0.000182e-6,     6283.143160294,  2.456157599 },
      {    0.000184e-6,     -227.526189440,  5.888038582 },
      {    0.000182e-6,    -6283.008539689,  0.241332086 },
      {    0.000228e-6,    -6284.056171060,  2.657323816 },
   /* 631, 640 */
      {    0.000166e-6,     7238.675591600,  5.930629110 },
      {    0.000167e-6,     3097.883822726,  5.570955333 },
      {    0.000159e-6,     -323.505416657,  5.786670700 },
      {    0.000154e-6,    -4136.910433516,  1.517805532 },
      {    0.000176e-6,    12029.347187887,  3.139266834 },
      {    0.000167e-6,    12132.439962106,  3.556352289 },
      {    0.000153e-6,      202.253395174,  1.463313961 },
      {    0.000157e-6,    17267.268201691,  1.586837396 },
      {    0.000142e-6,    83996.847317911,  0.022670115 },
      {    0.000152e-6,    17260.154654690,  0.708528947 },
   /* 641, 650 */
      {    0.000144e-6,     6084.003848555,  5.187075177 },
      {    0.000135e-6,     5756.566278634,  1.993229262 },
      {    0.000134e-6,     5750.203491159,  3.457197134 },
      {    0.000144e-6,     5326.786694021,  6.066193291 },
      {    0.000160e-6,    11015.106477335,  1.710431974 },
      {    0.000133e-6,     3634.621024518,  2.836451652 },
      {    0.000134e-6,    18073.704938650,  5.453106665 },
      {    0.000134e-6,     1162.474704408,  5.326898811 },
      {    0.000128e-6,     5642.198242609,  2.511652591 },
      {    0.000160e-6,      632.783739313,  5.628785365 },
   /* 651, 660 */
      {    0.000132e-6,    13916.019109642,  0.819294053 },
      {    0.000122e-6,    14314.168113050,  5.677408071 },
      {    0.000125e-6,    12359.966151546,  5.251984735 },
      {    0.000121e-6,     5749.452731634,  2.210924603 },
      {    0.000136e-6,     -245.831646229,  1.646502367 },
      {    0.000120e-6,     5757.317038160,  3.240883049 },
      {    0.000134e-6,    12146.667056108,  3.059480037 },
      {    0.000137e-6,     6206.809778716,  1.867105418 },
      {    0.000141e-6,    17253.041107690,  2.069217456 },
      {    0.000129e-6,    -7477.522860216,  2.781469314 },
   /* 661, 670 */
      {    0.000116e-6,     5540.085789459,  4.281176991 },
      {    0.000116e-6,     9779.108676125,  3.320925381 },
      {    0.000129e-6,     5237.921013804,  3.497704076 },
      {    0.000113e-6,     5959.570433334,  0.983210840 },
      {    0.000122e-6,     6282.095528923,  2.674938860 },
      {    0.000140e-6,      -11.045700264,  4.957936982 },
      {    0.000108e-6,    23543.230504682,  1.390113589 },
      {    0.000106e-6,   -12569.674818332,  0.429631317 },
      {    0.000110e-6,     -266.607041722,  5.501340197 },
      {    0.000115e-6,    12559.038152982,  4.691456618 },
   /* 671, 680 */
      {    0.000134e-6,    -2388.894020449,  0.577313584 },
      {    0.000109e-6,    10440.274292604,  6.218148717 },
      {    0.000102e-6,     -543.918059096,  1.477842615 },
      {    0.000108e-6,    21228.392023546,  2.237753948 },
      {    0.000101e-6,    -4535.059436924,  3.100492232 },
      {    0.000103e-6,       76.266071276,  5.594294322 },
      {    0.000104e-6,      949.175608970,  5.674287810 },
      {    0.000101e-6,    13517.870106233,  2.196632348 },
      {    0.000100e-6,    11933.367960670,  4.056084160 },
 
   /* T^2 */
      {    4.322990e-6,     6283.075849991,  2.642893748 },
   /* 681, 690 */
      {    0.406495e-6,        0.000000000,  4.712388980 },
      {    0.122605e-6,    12566.151699983,  2.438140634 },
      {    0.019476e-6,      213.299095438,  1.642186981 },
      {    0.016916e-6,      529.690965095,  4.510959344 },
      {    0.013374e-6,       -3.523118349,  1.502210314 },
      {    0.008042e-6,       26.298319800,  0.478549024 },
      {    0.007824e-6,      155.420399434,  5.254710405 },
      {    0.004894e-6,     5746.271337896,  4.683210850 },
      {    0.004875e-6,     5760.498431898,  0.759507698 },
      {    0.004416e-6,     5223.693919802,  6.028853166 },
   /* 691, 700 */
      {    0.004088e-6,       -7.113547001,  0.060926389 },
      {    0.004433e-6,    77713.771467920,  3.627734103 },
      {    0.003277e-6,    18849.227549974,  2.327912542 },
      {    0.002703e-6,     6062.663207553,  1.271941729 },
      {    0.003435e-6,     -775.522611324,  0.747446224 },
      {    0.002618e-6,     6076.890301554,  3.633715689 },
      {    0.003146e-6,      206.185548437,  5.647874613 },
      {    0.002544e-6,     1577.343542448,  6.232904270 },
      {    0.002218e-6,     -220.412642439,  1.309509946 },
      {    0.002197e-6,     5856.477659115,  2.407212349 },
   /* 701, 710 */
      {    0.002897e-6,     5753.384884897,  5.863842246 },
      {    0.001766e-6,      426.598190876,  0.754113147 },
      {    0.001738e-6,     -796.298006816,  2.714942671 },
      {    0.001695e-6,      522.577418094,  2.629369842 },
      {    0.001584e-6,     5507.553238667,  1.341138229 },
      {    0.001503e-6,     -242.728603974,  0.377699736 },
      {    0.001552e-6,     -536.804512095,  2.904684667 },
      {    0.001370e-6,     -398.149003408,  1.265599125 },
      {    0.001889e-6,    -5573.142801634,  4.413514859 },
      {    0.001722e-6,     6069.776754553,  2.445966339 },
   /* 711, 720 */
      {    0.001124e-6,     1059.381930189,  5.041799657 },
      {    0.001258e-6,      553.569402842,  3.849557278 },
      {    0.000831e-6,      951.718406251,  2.471094709 },
      {    0.000767e-6,     4694.002954708,  5.363125422 },
      {    0.000756e-6,     1349.867409659,  1.046195744 },
      {    0.000775e-6,      -11.045700264,  0.245548001 },
      {    0.000597e-6,     2146.165416475,  4.543268798 },
      {    0.000568e-6,     5216.580372801,  4.178853144 },
      {    0.000711e-6,     1748.016413067,  5.934271972 },
      {    0.000499e-6,    12036.460734888,  0.624434410 },
   /* 721, 730 */
      {    0.000671e-6,    -1194.447010225,  4.136047594 },
      {    0.000488e-6,     5849.364112115,  2.209679987 },
      {    0.000621e-6,     6438.496249426,  4.518860804 },
      {    0.000495e-6,    -6286.598968340,  1.868201275 },
      {    0.000456e-6,     5230.807466803,  1.271231591 },
      {    0.000451e-6,     5088.628839767,  0.084060889 },
      {    0.000435e-6,     5643.178563677,  3.324456609 },
      {    0.000387e-6,    10977.078804699,  4.052488477 },
      {    0.000547e-6,   161000.685737473,  2.841633844 },
      {    0.000522e-6,     3154.687084896,  2.171979966 },
   /* 731, 740 */
      {    0.000375e-6,     5486.777843175,  4.983027306 },
      {    0.000421e-6,     5863.591206116,  4.546432249 },
      {    0.000439e-6,     7084.896781115,  0.522967921 },
      {    0.000309e-6,     2544.314419883,  3.172606705 },
      {    0.000347e-6,     4690.479836359,  1.479586566 },
      {    0.000317e-6,      801.820931124,  3.553088096 },
      {    0.000262e-6,      419.484643875,  0.606635550 },
      {    0.000248e-6,     6836.645252834,  3.014082064 },
      {    0.000245e-6,    -1592.596013633,  5.519526220 },
      {    0.000225e-6,     4292.330832950,  2.877956536 },
   /* 741, 750 */
      {    0.000214e-6,     7234.794256242,  1.605227587 },
      {    0.000205e-6,     5767.611978898,  0.625804796 },
      {    0.000180e-6,    10447.387839604,  3.499954526 },
      {    0.000229e-6,      199.072001436,  5.632304604 },
      {    0.000214e-6,      639.897286314,  5.960227667 },
      {    0.000175e-6,     -433.711737877,  2.162417992 },
      {    0.000209e-6,      515.463871093,  2.322150893 },
      {    0.000173e-6,     6040.347246017,  2.556183691 },
      {    0.000184e-6,     6309.374169791,  4.732296790 },
      {    0.000227e-6,   149854.400134205,  5.385812217 },
   /* 751, 760 */
      {    0.000154e-6,     8031.092263058,  5.120720920 },
      {    0.000151e-6,     5739.157790895,  4.815000443 },
      {    0.000197e-6,     7632.943259650,  0.222827271 },
      {    0.000197e-6,       74.781598567,  3.910456770 },
      {    0.000138e-6,     6055.549660552,  1.397484253 },
      {    0.000149e-6,    -6127.655450557,  5.333727496 },
      {    0.000137e-6,     3894.181829542,  4.281749907 },
      {    0.000135e-6,     9437.762934887,  5.979971885 },
      {    0.000139e-6,    -2352.866153772,  4.715630782 },
      {    0.000142e-6,     6812.766815086,  0.513330157 },
   /* 761, 770 */
      {    0.000120e-6,    -4705.732307544,  0.194160689 },
      {    0.000131e-6,   -71430.695617928,  0.000379226 },
      {    0.000124e-6,     6279.552731642,  2.122264908 },
      {    0.000108e-6,    -6256.777530192,  0.883445696 },
 
   /* T^3 */
      {    0.143388e-6,     6283.075849991,  1.131453581 },
      {    0.006671e-6,    12566.151699983,  0.775148887 },
      {    0.001480e-6,      155.420399434,  0.480016880 },
      {    0.000934e-6,      213.299095438,  6.144453084 },
      {    0.000795e-6,      529.690965095,  2.941595619 },
      {    0.000673e-6,     5746.271337896,  0.120415406 },
   /* 771, 780 */
      {    0.000672e-6,     5760.498431898,  5.317009738 },
      {    0.000389e-6,     -220.412642439,  3.090323467 },
      {    0.000373e-6,     6062.663207553,  3.003551964 },
      {    0.000360e-6,     6076.890301554,  1.918913041 },
      {    0.000316e-6,      -21.340641002,  5.545798121 },
      {    0.000315e-6,     -242.728603974,  1.884932563 },
      {    0.000278e-6,      206.185548437,  1.266254859 },
      {    0.000238e-6,     -536.804512095,  4.532664830 },
      {    0.000185e-6,      522.577418094,  4.578313856 },
      {    0.000245e-6,    18849.227549974,  0.587467082 },
   /* 781, 787 */
      {    0.000180e-6,      426.598190876,  5.151178553 },
      {    0.000200e-6,      553.569402842,  5.355983739 },
      {    0.000141e-6,     5223.693919802,  1.336556009 },
      {    0.000104e-6,     5856.477659115,  4.239842759 },
 
   /* T^4 */
      {    0.003826e-6,     6283.075849991,  5.705257275 },
      {    0.000303e-6,    12566.151699983,  5.407132842 },
      {    0.000209e-6,      155.420399434,  1.989815753 }
   };


/* Time since J2000.0 in Julian millennia. */
   t = ((date1 - DJ00) + date2) / DJM;
 
/* ================= */
/* Topocentric terms */
/* ================= */
 
/* Convert UT to local solar time in radians. */
   tsol = fmod(ut, 1.0) * D2PI + elong;
 
/* FUNDAMENTAL ARGUMENTS:  Simon et al. 1994. */
 
/* Combine time argument (millennia) with deg/arcsec factor. */
   w = t / 3600.0;
 
/* Sun Mean Longitude. */
   elsun = fmod(280.46645683 + 1296027711.03429 * w, 360.0) * DD2R;
 
/* Sun Mean Anomaly. */
   emsun = fmod(357.52910918 + 1295965810.481 * w, 360.0) * DD2R;
 
/* Mean Elongation of Moon from Sun. */
   d = fmod(297.85019547 + 16029616012.090 * w, 360.0) * DD2R;
 
/* Mean Longitude of Jupiter. */
   elj = fmod(34.35151874 + 109306899.89453 * w, 360.0) * DD2R;
 
/* Mean Longitude of Saturn. */
   els = fmod(50.07744430 + 44046398.47038 * w, 360.0) * DD2R;
 
/* TOPOCENTRIC TERMS:  Moyer 1981 and Murray 1983. */
   wt =   +  0.00029e-10 * u * sin(tsol + elsun - els)
          +  0.00100e-10 * u * sin(tsol - 2.0 * emsun)
          +  0.00133e-10 * u * sin(tsol - d)
          +  0.00133e-10 * u * sin(tsol + elsun - elj)
          -  0.00229e-10 * u * sin(tsol + 2.0 * elsun + emsun)
          -  0.02200e-10 * v * cos(elsun + emsun)
          +  0.05312e-10 * u * sin(tsol - emsun)
          -  0.13677e-10 * u * sin(tsol + 2.0 * elsun)
          -  1.31840e-10 * v * cos(elsun)
          +  3.17679e-10 * u * sin(tsol);
 
/* ===================== */
/* Fairhead et al. model */
/* ===================== */
 
/* T**0 */
   w0 = 0;
   for (j = 473; j >= 0; j--) {
      w0 += fairhd[j][0] * sin(fairhd[j][1] * t + fairhd[j][2]);
   }
 
/* T**1 */
   w1 = 0;
   for (j = 678; j >= 474; j--) {
      w1 += fairhd[j][0] * sin(fairhd[j][1] * t + fairhd[j][2]);
   }
 
/* T**2 */
   w2 = 0;
   for (j = 763; j >= 679; j--) {
      w2 += fairhd[j][0] * sin(fairhd[j][1] * t + fairhd[j][2]);
   }
 
/* T**3 */
   w3 = 0;
   for (j = 783; j >= 764; j--) {
      w3 += fairhd[j][0] * sin(fairhd[j][1] * t + fairhd[j][2]);
   }
 
/* T**4 */
   w4 = 0;
   for (j = 786; j >= 784; j--) {
      w4 += fairhd[j][0] * sin(fairhd[j][1] * t + fairhd[j][2]);
   }
 
/* Multiply by powers of T and combine. */
   wf = t * (t * (t * (t * w4 + w3) + w2) + w1) + w0;
 
/* Adjustments to use JPL planetary masses instead of IAU. */
   wj =   0.00065e-6 * sin(6069.776754 * t + 4.021194) +
          0.00033e-6 * sin( 213.299095 * t + 5.543132) +
        (-0.00196e-6 * sin(6208.294251 * t + 5.696701)) +
        (-0.00173e-6 * sin(  74.781599 * t + 2.435900)) +
          0.03638e-6 * t * t;
 
/* ============ */
/* Final result */
/* ============ */
 
/* TDB-TT in seconds. */
   w = wt + wf + wj;

   return w;

/*-----------------------------------------------------------------------
**
**  Copyright (C) 2008
**  Standards Of Fundamental Astronomy Review Board
**  of the International Astronomical Union.
**
**  =====================
**  SOFA Software License
**  =====================
**
**  NOTICE TO USER:
**
**  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
**  WHICH APPLY TO ITS USE.
**
**  1. The Software is owned by the IAU SOFA Review Board ("the Board").
**
**  2. Permission is granted to anyone to use the SOFA software for any
**     purpose, including commercial applications, free of charge and
**     without payment of royalties, subject to the conditions and 
**     restrictions listed below.
**
**  3. You (the user) may copy and adapt the SOFA software and its 
**     algorithms for your own purposes and you may copy and distribute
**     a resulting "derived work" to others on a world-wide, royalty-free 
**     basis, provided that the derived work complies with the following
**     requirements: 
**
**     a) Your work shall be marked or carry a statement that it (i) uses
**        routines and computations derived by you from software provided 
**        by SOFA under license to you; and (ii) does not contain
**        software provided by SOFA or software that has been distributed
**        by or endorsed by SOFA.
**
**     b) The source code of your derived work must contain descriptions
**        of how the derived work is based upon and/or differs from the
**        original SOFA software.
**
**     c) The name(s) of all routine(s) that you distribute shall differ
**        from the SOFA names, even when the SOFA content has not been
**        otherwise changed.
**
**     d) The routine-naming prefix "iau" shall not be used.
**
**     e) The origin of the SOFA components of your derived work must not
**        be misrepresented;  you must not claim that you wrote the
**        original software, nor file a patent application for SOFA
**        software or algorithms embedded in the SOFA software.
**
**     f) These requirements must be reproduced intact in any source
**        distribution and shall apply to anyone to whom you have granted 
**        a further right to modify the source code of your derived work.
**
**  4. In any published work or commercial products which includes
**     results achieved by using the SOFA software, you shall acknowledge
**     that the SOFA software was used in obtaining those results.
**
**  5. You shall not cause the SOFA software to be brought into
**     disrepute, either by misuse, or use for inappropriate tasks, or by
**     inappropriate modification.
**
**  6. The SOFA software is provided "as is" and the Board makes no 
**     warranty as to its use or performance.   The Board does not and 
**     cannot warrant the performance or results which the user may obtain 
**     by using the SOFA software.  The Board makes no warranties, express 
**     or implied, as to non-infringement of third party rights,
**     merchantability, or fitness for any particular purpose.  In no
**     event will the Board be liable to the user for any consequential,
**     incidental, or special damages, including any lost profits or lost
**     savings, even if a Board representative has been advised of such
**     damages, or for any claim by any third party.
**
**  7. The provision of any version of the SOFA software under the terms 
**     and conditions specified herein does not imply that future
**     versions will also be made available under the same terms and
**     conditions.
**
**  Correspondence concerning SOFA software should be addressed as
**  follows:
**
**     Internet email: sofa@rl.ac.uk
**     Postal address: IAU SOFA Center
**                     Rutherford Appleton Laboratory
**                     Chilton, Didcot, Oxon OX11 0QX
**                     United Kingdom
**
**---------------------------------------------------------------------*/
}


