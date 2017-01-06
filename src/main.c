/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*

    Solar-system Model and Orbit Determination
    SMOD
  
    Version:    1008mars (20 Aug 2010)
    Version:    smod_1201 (1 Jan 2012)

    Copyright (c) 2010 SHANG Kun (shangkun@shao.ac.cn) All Right Reserved 
    Copyright (c) 2012 Kun Shang (shang.34@osu.edu) All Right Reserved 

    NOT FULLY VALIDATED -- SUBJECT TO CORRECTION

*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

/*#define MSDOS */
#define LINUX



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eph_manager.h" 
#include "novas.h"
#include "smod.h"
#include "SpiceUsr.h"

//#define orbformat  

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
main (int argc, char *argv[])
{
    FILE *fp_stdin, *fp_simu, *fp_result, *fp_res;
    short int year, month, day, error, de_num;
    int uplink, iter, part, all, teout;
    int n, i, nsta, ndata, dim_par, sta1, sta3, consfg, dconfg, biasfg,slvls, genrel, maxiter,
        dbiafg,  type, sta_out, type_sel, sta3_sel, sta1_sel, lps, mjds, kpl, ntide, indicator;
    double pi[3], pe[3], llh[3], tt, tp, h, hcut, sigma2 = 0, sigma = 0;

    double tdb, tdbs, tdbe, utc, utcs, utce, ts_sim, te_sim, t_sim,
        start, finish, duration, xsm[6], xau[6], xtm[6], x2[6], dist, velt, 
        calculable, noise, azimuth, elevation, elcutoff, data_value, 
        rmsold, rmsnew, rms2n, step_day, jd_sel, ts_sel, te_sel,
        *state, *bi, *bit, *btbi, *btyi, *BTry, *BTrB, 
        *p, *px, *dx, *x0, *xnew, *xold, weight, dyi[1], tcdop, 
        jd_beg, jd_end, criteria, gm2de, correl, ts_orb, te_orb, step_orb, 
        ts_sta, te_sta, step_sta, step_sim, varcons, vardcon, varbias, 
        vardbia, tx[3], ty[3], tz[3], vx[3] = {1,0,0}, vy[3] = {0,1,0},
        vz[3] = {0,0,1}, te[9], gm[11]= {0}, var[6], vt[3], vc[3], 
        geo_loc[10][3], rtn_p[3], rtn_v[3], ele[6], nele, errstd = 0, seekrd;
    char line[200], card[20], f_eop[200], f_eph[200], f_ot[200], f_obs[200],
        f_grv[2][200]={"\0", "\0"}, stdname[100];  
    double data_chk64;
    real128 data_chk128;



    srand (time(NULL));
    seekrd = (rand() %100)/100.0;

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
 
    printf("Checking Quadruple-precision...\n\n");
        
    printf("shortI = %2d bytes\n", sizeof(short int));
    printf("int    = %2d bytes\n", sizeof(int));
    printf("longI  = %2d bytes\n", sizeof(long int));
    printf("float  = %2d bytes\n", sizeof(float));
    printf("double = %2d bytes\n", sizeof(double));
    printf("longD  = %2d bytes\n", sizeof(real128));

    data_chk64 = 0.12345678901234567890123456789012345678901234567890L;
    data_chk128 = 0.12345678901234567890123456789012345678901234567890L;
    printf("Input   : 0.12345678901234567890123456789012345678901234567890\n");
    printf("Output2 : %52.50f\n", data_chk64);
    printf("Output4 : %52.50Lf\n", data_chk128);
    printf("\n");

// 1AU = 1.5e11m = 0.15e12m = 0.15e15mm = 0.15e18micormeter

   
	if (argc == 1)
	{
		printf ("input file name: ");
		scanf ("%s", stdname);
	}
	else if (argc == 2)
	{
		sscanf (argv[1], "%s", stdname);
	}
	else
	{
        printf ("input argument error!\n");
       	getch ();
        exit (0);
	}
    if ( (fp_stdin = fopen (stdname,"r")) == NULL)
    {
        printf ("Cannot open stdin file!\n");
       	getch();
        exit (0);
    }
    printf ("\nreading input file \"%s\"...\n\n", stdname);


    MOBS = 0;  
    MSRP = 0;  
    MTK2 = 0;  
    MGCS = 0;  

    all = 1;
    teout = 0;
    while (feof(fp_stdin) == 0)
    {
        card[1] = '\0';
        fgets (line, 100, fp_stdin);
        sscanf (line, "%s", card);	

        if (strcmp (card, "planet") ==0)	
        {
            sscanf (line, "%*s%*s%d", &i);
//            i = i - 1;
            sscanf (line, "%*s%*s%*d%d%lf", &PERB[i], &gm[i]);
        }
		
        if (strcmp (card, "elcutoff") ==0)	
        {
            sscanf (line, "%*s%lf", &elcutoff);
        }
        
        if (strcmp (card, "hcut") ==0)	
        {
            sscanf (line, "%*s%lf", &hcut);
        }
        
        if (strcmp (card, "station") ==0)
        {
            sscanf (line, "%*s%*s%d", &i);
            i = i - 1;
            sscanf (line, "%*s%*s%*d%lf%lf%lf", &geo_loc[i][0], 
                &geo_loc[i][1], &geo_loc[i][2]);
        }

        if (strcmp (card, "indicator") ==0)
        {
            sscanf (line, "%*s%d", &indicator);
        }
        if (strcmp (card, "center") ==0)
        {
            sscanf (line, "%*s%d%lf%lf", &CT, &GMCT, &RCT);
//            CT = i - 1;
        }

        if (strcmp (card, "epoch") ==0)
        {
            sscanf (line, "%*s%hd%hd%hd", &year, &month, &day);
        }
        if (strcmp (card, "utcsec") ==0)
        {
            sscanf (line, "%*s%lf%lf", &utcs, &utce);
        }

        if (strcmp (card, "kepler") ==0)
        {
            sscanf (line, "%*s%d", &kpl);
        }

        if (strcmp (card, "element1") ==0)
        {
            sscanf (line, "%*s%lf%lf%lf", &ele[0], &ele[1], &ele[2]);
        }
        if (strcmp (card, "element2") ==0)
        {
            sscanf (line, "%*s%lf%lf%lf", &ele[3], &ele[4], &ele[5]);
        }

        if (strcmp (card, "intstep") ==0)	
        {
            sscanf (line, "%*s%lf", &STEP_OR);
        }
        if (strcmp (card, "orbfile") ==0)	
        {
            sscanf (line, "%*s%lf%lf%lf", &ts_orb, &te_orb, &step_orb);
        }
        if (strcmp (card, "f2istep") ==0)	
        {
            sscanf (line, "%*s%lf", &STEP_TE);
        }
        if (strcmp (card, "stafile") ==0)	
        {
            teout = 1;
            sscanf (line, "%*s%d%lf%lf%lf", 
                &sta_out, &ts_sta, &te_sta, &step_sta);
            sta_out = sta_out - 1;
        }

        if (strcmp (card, "simulate") ==0)	
        {
            sta3 = 0; sta1 = 0;
            sscanf (line, "%*s%d%lf%lf%lf%lf%d%d%lf", 
                &type, &noise, &ts_sim, &te_sim, 
                &step_sim, &sta3, &sta1, &t_sim);
            sta3 = sta3 - 1;
            sta1 = sta1 - 1;
        }
        if (strcmp (card, "relativ") ==0)	
        {
   	        sscanf (line, "%*s%d", &genrel);
        }
        if (strcmp (card, "rmsconv") ==0)	
        {
   	        sscanf (line, "%*s%lf", &criteria);
        }
        if (strcmp (card, "itermax") ==0)	
        {
        sscanf (line, "%*s%d", &maxiter);
        }
        if (strcmp (card, "varcov1") ==0)	
        {
   	        sscanf (line, "%*s%lf%lf%lf", &var[0], &var[1], &var[2]);
        }
        if (strcmp (card, "varcov2") ==0)	
        {
   	        sscanf (line, "%*s%lf%lf%lf", &var[3], &var[4], &var[5]);
        }


        if (strcmp (card, "solrad") ==0)	
        {
  	        sscanf (line, "%*s%lf%d%lf", &CONS, &consfg, &varcons);
        }
        if (strcmp (card, "solraddot") ==0)	
        {
  	        sscanf (line, "%*s%lf%d%lf", &DCON, &dconfg, &vardcon);
        }
        if (strcmp (card, "bias") ==0)	
        {
            sscanf (line, "%*s%lf%d%lf", &BIAS, &biasfg, &varbias);
        }
        if (strcmp (card, "biasdot") ==0)	
        {
            sscanf (line, "%*s%lf%d%lf", &DBIA, &dbiafg, &vardbia);
        }
        if (strcmp (card, "select") ==0)	
        {
            all = 0;
            sta3_sel = 0; sta1_sel = 0;
            sscanf (line, "%*s%d%lf%lf%d%d", 
                &type_sel, &ts_sel, &te_sel, &sta3_sel, &sta1_sel);
            sta3_sel = sta3_sel - 1;
            sta1_sel = sta1_sel - 1;
        }
        if (strcmp (card, "solvels") == 0)
        {
            sscanf (line, "%*s%d", &slvls);
        }


        if (strcmp (card, "gravity") ==0)   
        {
            sscanf (line, "%*s%d%d", &NMAX, &MMAX);
        }
        if (strcmp (card, "amratio") ==0)   
        {
            sscanf (line, "%*s%lf", &AMR);
        }
        if (strcmp (card, "nbody") ==0) 
        {
            sscanf (line, "%*s%d", &NBODY);
        }
        if (strcmp (card, "reltiv") ==0)    
        {
            sscanf (line, "%*s%d", &RELTIV);
        }
        if (strcmp (card, "permtide") ==0)
        {
            sscanf (line, "%*s%d", &PERMT);
        }
        if (strcmp (card, "bodytide") ==0)
        {
            sscanf (line, "%*s%d", &STIDE);
        }
        if (strcmp (card, "oceantide") ==0)
        {
            sscanf (line, "%*s%d", &OTIDE);
        }

        if (strcmp (card,"ephfile") ==0)
        {
            sscanf (line, "%*s%s", f_eph);
        }
        if (strcmp (card,"eopfile") ==0)
        {
            sscanf (line, "%*s%s", f_eop);
        }
        if (strcmp (card,"grvfile") ==0)
        {
            sscanf (line, "%*s%s", f_grv[0]);
        }
        if (strcmp (card,"oceanfile") ==0)
        {
            sscanf (line, "%*s%s", f_ot);
        }

        if (strcmp (card,"obsfile") ==0)
        {
            sscanf (line, "%*s%s", f_obs);
        }


        if (strcmp (card, "MOBS") ==0)
        {
            sscanf (line, "%*s%d%lf%lf", &MOBS, &BASB, &BAST);
        }

        if (strcmp (card, "MSRP") ==0)
        {
            sscanf (line, "%*s%d%lf%lf", &MSRP, &SRPB, &SRPT);
        }


        if (strcmp (card, "MTK2") ==0)
        {
            sscanf (line, "%*s%d%lf", &MTK2, &K2);
        }

        if (strcmp (card, "MGCS") ==0)
        {
            sscanf (line, "%*s%d%d%d%lf", &CSinfo[MGCS].n, &CSinfo[MGCS].m, 
                &CSinfo[MGCS].cs, &CSinfo[MGCS].initv);
            MGCS ++;
        }


    }
    
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

    if ((error = ephem_open (f_eph, &jd_beg,&jd_end,&de_num)) != 0)
    {
        if (error == 1)
            printf ("JPL ephemeris file not found.\n");
        else
            printf ("Error reading JPL ephemeris file header.\n");
        getch();
        exit (error);
    }
    else
    {
        printf ("JPL ephemeris DE%d open. Start JD = %10.2f  End JD = %10.2f\n",
            de_num, jd_beg, jd_end);
        printf ("\n");
    }

    if (CT == 20)
    {
        furnsh_c ( "inputs/sat052.bsp" );
        printf ("SPICE ephemeric file sat052 open\n\n");
    }


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

//    gm2de = 86400.0 * 86400.0 / AU / AU / AU;
    for (n = 0; n <= 10; n++)
//        GMDE[n] = gm[n] * gm2de;
        GMDE[n] = gm[n];
    JD0 = julian_date (year, month, day,0);



    mjds = (int)(JD0 - 2400000.5 + utcs/86400.0) - 2;
    NEOP = (int)(JD0 - 2400000.5 + utce/86400.0 - mjds) + 3;

    EOPMT  = (double *) calloc (NEOP * 6, sizeof(double));
    openeop (f_eop, mjds, NEOP, EOPMT);




    if (NMAX >= 2)
    {
        COEFG  = (double *) calloc ((NMAX + 1) * (NMAX + 1), sizeof(double));
        opengrv (f_grv, COEFG, NMAX, MMAX);

        NSMAX = 2;
        if (CT == 2 && STIDE > 1)
            NSMAX = 4;
        COEFS = (double *) calloc ( (NSMAX + 1) * (NSMAX + 1), sizeof(double));

        NOMAX = 4;
        COEFO = (double *) calloc ( (NOMAX + 1) * (NOMAX + 1), sizeof(double));

        NFES = 59462;
        otfes = (OTStruct *) calloc ( NFES, sizeof(OTStruct));
        openotcs (f_ot);

    }



    if (kpl == 0)
    {
        for (n = 0; n < 6; n++)
        {
            xsm[n] = ele[n];
        }
    }
    else if (kpl == 1)
    {
        aei2xyz (ele, &xsm[0], &xsm[3]);
    }
    else
    {
        printf ("kepler should be equal to 0 or 1 !\n");
        exit(0);
    }

    nsta = 1;

    ACCURACY = 1;

//    BASB = 2; BAST = 0.03; SRPB = 1; SRPT = 0.004; K2 = 0.0; 
//    COEFG[2] = COEFG[2] + 2e-7; 
//    COEFG[3] = -4e-6;
//    BASB = 0; BAST = 0.00; SRPB = 0; SRPT = 0.000; K2 = 0.3;
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/*--orbit integration & observable simulation--*/

    if (indicator == 1 || indicator == 2)
    {
        start = clock();  
        DIM_OR = (int)( (utce - utcs) / STEP_OR) + 1;	
        OR_EPH = (double *) calloc ( DIM_OR * 7, sizeof(double));

        for (i = 0; i < DIM_OR; i++)
        {
            utc = utcs + STEP_OR * i;
            lps = getlps (JD0 + utcs/86400.0);
            tt = utc + (lps + 32.184);

            OR_EPH[i * 7] = tt;
            for (n = 0; n < 6; n++)
                OR_EPH[i * 7 + 1 + n] = xsm[n]; 

            rkf78 (JD0, tt, STEP_OR, xsm, 6, fun_accel);
        }

        finish = clock();  
        duration = (double)(finish - start) / CLOCKS_PER_SEC; 
        printf ("\nFinish integration of orbit: \t%.4f seconds\n", duration);
        
        printf ("utc = %f\n", utc);

        pt_orb (ts_orb, te_orb, step_orb, 7);

        finish = clock();  
        duration = (double)(finish - start) / CLOCKS_PER_SEC; 
        printf ("\nFinish output of orbit: \t%.4f seconds\n", duration);

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
        /*--observable simulation--*/
        if (indicator == 2)
        {
            if((fp_simu=fopen("simulate.dat","w"))==NULL)
            {
                printf("Cannot write file_obs?\n");
                getch();
                exit(0);
            }

//        BASB = 0; BAST = 0; SRPB = 0; SRPT = 0;
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
            for (n = 0; n < nsta; n++)
            {
                for (utc = ts_sim; utc < te_sim; utc = utc + step_sim)
                {
                    part = 0;
                    if (type == 22)
                    {
                        uplink = 1;
                        tcdop = t_sim;
                        simula_dople (utc, tcdop, geo_loc[sta3], uplink, 
                            geo_loc[sta1], genrel, &calculable, &azimuth, 
                            &elevation, part, x2);
                    }
                    if (type == 99)
                       h = obs_alt (JD0, utc, &calculable, part, x2);
                    if (type == 98)
                       h = obs_dsn (JD0, utc, &calculable, part, x2);
                    if (type == 97)
                       h = obs_vel (JD0, utc, &calculable, part, x2);
                    calculable = calculable + mgrn1(0.0,noise,&seekrd);
                    if (h < hcut)
                    {
                        fprintf(fp_simu, "%.2f  %14.6f  %d  %10.4f  %16.8f  ", 
                            JD0, utc, type, noise, calculable);
                        fprintf(fp_simu, "\n");
                    }
                }
/*
                    if ((elevation > elcutoff) && (type !=99) && (type != 98))
                    {
                        fprintf(fp_simu, "%.2f  %14.6f  %d  %10.4f  %35.20f  ", 
                            JD0, utc, type, noise, calculable);
                        fprintf(fp_simu, "%d  %d  %14.6f  %14.6f  %14.6f\n", 
                            sta3 + 1, sta1 + 1, t_sim, azimuth, elevation);
                    }
*/  
            }
      
            fclose(fp_simu);
            finish = clock();  
            duration = (double)(finish - start) / CLOCKS_PER_SEC; 
            printf ("\nFinish output of observation: \t%.4f seconds\n", 
                duration);
        }
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
        finish = clock();  
        duration = (double)(finish - start) / CLOCKS_PER_SEC; 
        printf ("\ntotal duration = %f\n\n", duration);
        free (TE_EPH);
        free (OR_EPH);
        printf ("\nNormal end of ECHO!\n\npress any key to finish...\n");
        getch();
        exit(0);
    }




/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/* orbit determination */
    if (indicator == 3)
    {
        start = clock();  
        if((fp_result=fopen("result.txt","w"))==NULL)
        {
            printf("Cannot write result.txt!\n");
          	getch();
            exit(0);
        }
        if((fp_simu=fopen(f_obs,"r"))==NULL)
        {
            printf("Cannot open obs data?\n");
          	getch();
            exit(0);
        }
        if((fp_res=fopen("residua.dat","w"))==NULL)
        {
            printf("Cannot write residua.dat!\n");
           	getch();
            exit(0);
        }
  
        getsolvefor();
      
        bi = (double *) calloc ( MSOL, sizeof(double));
        bit = (double *) calloc ( MSOL, sizeof(double));
        btbi = (double *) calloc ( MSOL * MSOL, sizeof(double));
        btyi = (double *) calloc ( MSOL, sizeof(double));
        BTry = (double *) calloc ( MSOL, sizeof(double));
        BTrB = (double *) calloc ( MSOL * MSOL, sizeof(double));
        p = (double *) calloc ( MSOL * MSOL, sizeof(double));
        px = (double *) calloc ( MSOL, sizeof(double));
        dx = (double *) calloc ( MSOL, sizeof(double));
        x0 = (double *) calloc ( MSOL, sizeof(double));
        xnew = (double *) calloc ( MSOL, sizeof(double));
        xold = (double *) calloc ( MSOL, sizeof(double));
        state = (double *) calloc ( MSTA, sizeof(double));

/*
        for ( i = 0; i < MSOL * MSOL; i ++)
            p[i] = 0;
        p[0]                = var[0];
        p[1 * MSOL + 1] = var[1];
        p[2 * MSOL + 2] = var[2];
        p[3 * MSOL + 3] = var[3];
        p[4 * MSOL + 4] = var[4];
        p[5 * MSOL + 5] = var[5];
//        p[6 * MSOL + 6] = varcons;
//        p[7 * MSOL + 7]	= vardcon;	
//        p[8 * MSOL + 8]	= varbias;	
//        p[9 * MSOL + 9]	= vardbia;	
//        brinv(p,MSOL);
*/

//        for ( i = 0; i < MSOL * MSOL; i ++)
//            printf("%f\n", p[i]);
        

        initsolvefor (xsm, x0);
    
        for (n = 0; n < MSOL; n ++)
            xnew[n] = x0[n];
        
//        updsolvefor (xnew);            
//        initsolvefor (xsm, xnew);

        rmsold = 2e10;
        rmsnew = 1e10;


        DIM_OR = (int)( (utce - utcs) / STEP_OR) + 1;	
//        dim_par = 42 + 6 * MDYN + 1;
        OR_EPH = (double *) calloc ( DIM_OR * (MSTA + 1) , sizeof(double));	
        

        for (iter = 0; iter < maxiter; iter ++)
        {
            if ( fabs((rmsold - rmsnew) / rmsold) < criteria)
            {
                iter = 99;
            }

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/* initialize transition matrix and orbit trajectory, */
            rmsold = rmsnew;
            for (n = 0; n < MSOL; n ++)
            {
                xold[n] = xnew[n];
                dx[n] = 0;
            }
//            BASB = xold[6];
//            BAST = xold[7]; 
//            SRPB = xold[8];
//            SRPT = xold[9];

//            SRPB = xold[6];
//            SRPT = xold[7];
            for (n = 0; n < MSTA; n++)
                state[n] = 0;
            for (n = 0; n < 6; n++)
      	        state[n] = xold[n];
            for (n = 0; n < 6; n++)
                state[n * 6 + n + 6] = 1;
            
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/* integrate transition matrix and orbit trajectory*/

            for (i = 0; i < DIM_OR; i++)
            {
                utc = utcs + STEP_OR * i;
                lps = getlps (JD0 + utcs/86400.0);
                tt = utc + (lps + 32.184);
    
                OR_EPH[i * (MSTA + 1)] = tt;
                for (n = 0; n < MSTA; n++)
                    OR_EPH[i * (MSTA + 1) + 1 + n] = state[n]; 

                rkf78 (JD0, tt, STEP_OR, state, MSTA, fun_accel);

            }
            finish = clock();  
            duration = (double)(finish - start) / CLOCKS_PER_SEC; 
            printf ("OD: Finish integration of orbit: \t%.4f seconds\n", 
                duration);

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/* print orbit trajectory */
            if (iter == 99 || iter == maxiter - 1)
            {
                pt_orb (ts_orb, te_orb, step_orb, MSTA + 1);

                finish = clock();  
                duration = (double)(finish - start) / CLOCKS_PER_SEC; 
                printf ("\nFinish output of orbit: \t%.4f seconds\n", duration);

            }







/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/* read observation and form BtB & Bty matrix*/
            rms2n = 0;
            sigma2 = 0;
            ndata = 0;
            for (n = 0; n < MSOL * MSOL; n++)
            {
                BTrB[n] =0;
            }
            for (n = 0; n < MSOL; n++)
            {
                BTry[n] = 0;
            }
            rewind(fp_simu);
            while (1)
            {
                if (fgets (line, 200, fp_simu) ==NULL) break;

//                fgets (line, 100, fp_simu);
                sscanf(line, "%lf%lf%d%lf%lf%d%d%lf", 
                    &jd_sel, &utc, &type, &noise, &data_value, 
                    &sta3, &sta1, &t_sim);

//            data_value = data_value + mgrn1(0.0,errstd,&seekrd);
                if (noise == 0) noise = 1;
                weight = 1.0 / noise / noise; //range: m; velocity: m/s
//                weight = 1; //range: m; velocity: m/s
                sta3 = sta3 - 1;
                sta1 = sta1 - 1;
                if (all == 0)
                {
                    if (jd_sel != JD0 || utc > te_sel 
                        || utc < ts_sel || type != type_sel )
//                        || sta3 != sta3_sel || sta1 != sta1_sel)
                    continue;		
                }
                
                part = 1;
                if (type == 22)
                {
                    uplink = 1;
                    tcdop = t_sim;
                    simula_dople (utc, tcdop, geo_loc[sta3], uplink, 
                        geo_loc[sta1], genrel, &calculable, &azimuth, 
                        &elevation, part, bi);
                }
                if (type == 99)
                {
                    obs_alt (JD0, utc, &calculable, part, bi);
                }

                if (type == 98)
                {
                    obs_dsn (JD0, utc, &calculable, part, bi);
                }
                if (type == 97)
                {
                    obs_vel (JD0, utc, &calculable, part, bi);
                }



                dyi[0] = data_value  - calculable;
                rms2n = rms2n + dyi[0] * dyi[0];
                sigma2 = sigma2 + dyi[0] * dyi[0] * weight;

                if (iter == 99 || iter == maxiter - 1)
                {
                    fprintf (fp_res, "%.2f\t%35.20f\t%35.20f\t%35.20f\n", 
                        utc, dyi[0], data_value, calculable);
                }

                mt (bi, 1, MSOL, bit);
                brmul(bit,bi,MSOL,1,MSOL,btbi);
                brmul(bit,dyi,MSOL,1,1,btyi);
                
                for (n = 0; n < MSOL * MSOL; n++)
                {
                    BTrB[n] = BTrB[n] + btbi[n] * weight;
                }
                for (n = 0; n < MSOL; n++)
                {
                    BTry[n] = BTry[n] + btyi[n] * weight;
                }
                ndata++;
            }
            fflush (fp_res);
            rmsnew = sqrt (rms2n / ndata);
            sigma = sqrt (sigma2 / (ndata - MSOL));

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/* add prior information and weight */
/*            for (n = 0; n < MSOL * MSOL; n++)
            {
                BTrB[n] = BTrB[n] + p[n];
            }
            for (n = 0; n < MSOL; n++)
            {
                px[n] = 0;
                for (i = 0; i < MSOL; i++)
                {
                    px[n] = px[n] + p[n * MSOL + i] * (xnew[i] - x0[i]);
                }
                BTry[n] = BTry[n] + px[n];
            }
*/
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
/* solve least-square */
            for (n = 0; n < MSOL; n++)
            {
                for (i = 0; i < MSOL; i++)
                {
                    fprintf (fp_result, "%12.4e ", BTrB[n * MSOL + i]);
                }
//                fprintf (fp_result, "\n");
                fprintf (fp_result, "%12.4e\n", BTry[n]);
            }
            fprintf (fp_result, "\n");


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
//            brinv (BTrB,MSOL);
//            bssgj (BTrB,MSOL);
//            brmul (BTrB,BTry,MSOL,MSOL,1,dx);
  
            if (slvls == 1)
                solvels_chol(BTrB, MSOL, BTry, dx, 0);
            if (slvls == 2)            
                solvegaus(BTrB, MSOL, BTry, dx);
       
            for (n = 0; n < MSOL; n++)
                xnew[n] = xold[n] + dx[n];
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/

            updsolvefor (xnew);            

            for (n = 0; n < MSOL; n++)
            {
                fprintf (fp_result, "%3d: ", n);
                for (i = 0; i <= n; i++)
                {
                    correl = BTrB[i * MSOL + n] / 
                        sqrt (BTrB[n * MSOL + n] * BTrB[i * MSOL + i]);
                    fprintf (fp_result, "%8.4f ", correl);
                }
                fprintf (fp_result, "\n");
            }
            fprintf (fp_result, "\n");
            
            for (n = 0; n < MSOL; n++)
            {
                for (i = 0; i <= n; i++)
                {
                    fprintf (fp_result, "%12.4e ", BTrB[i*MSOL+n]);
                }
//                fprintf (fp_result, "%12.4e\n", BTry[n]);
                fprintf (fp_result, "\n");

            }
            fprintf (fp_result, "\n");

            for (n = 0; n < MSOL; n++)
            {
                fprintf (fp_result, "%22.10f %22.10f %22.10f %22.10f %22.10f %22.10f\n", 
                    x0[n], xold[n], xnew[n], x0[n] - xnew[n], dx[n], sqrt(BTrB[n*MSOL+n]));
            }
            fprintf (fp_result, "\niter: %2d\t rms:\t%g\t%f\n\n\n", 
                iter, rmsnew, fabs((rmsold - rmsnew) / rmsold));
            
            fflush(fp_result);
            finish = clock();  
            duration = (double)(finish - start) / CLOCKS_PER_SEC; 
            printf ("OD: Finish estimation of orbit: \t%.4f seconds\n", 
                duration);
            printf ("iter: %2d\t rms:\t%g\t%f\t sigma = %f\n\n", 
                iter, rmsnew, fabs((rmsold - rmsnew) / rmsold), sigma);

        }
        
        for (n = 0; n < MSOL; n++)
        {
            printf ("%22.4f %22.4f\n", 
                xnew[n] * 1e6, sqrt(BTrB[n*MSOL+n]) * 1e6);
        }

        ephem_close();
        fclose(fp_result);
        fclose(fp_simu);
        fclose(fp_res);
        free (TE_EPH);
        free (OR_EPH); 
        free (bi); 
        free (bit); 
        free (btbi); 
        free (btyi); 
        free (BTry);
        free (BTrB); 
        free (p); 
        free (px); 
        free (dx); 
        free (x0); 
        free (xnew); 
        free (xold);
//        if (MGCS > 0) free(CSinfo);
        printf ("\nNormal end of ECHO!\n\npress any key to finish...\n");
        getch();
        exit(0);
    }

    getch();
    exit(0);

}
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+--*/
