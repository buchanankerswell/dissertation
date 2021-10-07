/* LEFT+Right Side or Err for X-Stokes Equation */
/* Stoks equation initial form */
/* dSIGxx/dX + dSIGxy/dY - dP/dX = - RO*Gx */
double xstokserr(long int m1, long int m2, int ynerr)
/* m1,m2 - node X,Y number */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Counters */
int n1,n2;
long int v[4];
/* Err Buf */
double leftx,rightx,nueff;
/* Distances */
double xkf=(gx[m1+1]-gx[m1-1])/2.0,ykf=gy[m2+1]-gy[m2];
/**/
/**/
/**/
/* Staggered Nodes num */
/*                      [0]                [2]   */
/*                  Ro0,Sxy0,Nu0                 */
/*                                               */
/*         Pr1,Sxx1    <Vx0>  Pr3,Sxx3           */
/*                                               */
/*                      [1]                [3]   */
/*                  RO1,Sxy1,Nu1                 */
/*                                               */
v[0]=m1*ynumy+m2;v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* RIGHT parts of X-Stokes */
/* dSIGik/dXk - dP/dXi = - Gi*RO */
rightx  = -GXKOEF*(ro[v[0]]+ro[v[1]])/2.0;
/**/
/**/
/**/
/* Return val for LSQ err ----------------------------*/
if(ynerr==1)
	{
	/* LEFT part of X-Stokes */
	/* dSIGxx/dX + dSIGxy/dY - dP/dX = - RO*Gx */
	/**/
	/* dSIGxx/dX */
	leftx =(sxx[v[3]]-sxx[v[1]])/xkf;
	/* dSIGxy/dY */
	leftx+=(sxy[v[1]]-sxy[v[0]])/ykf;
	/* -dP/dX */
	leftx-=(pr[v[3]]-pr[v[1]])/xkf;
	/* -gx/dt/2*(vx*dRHO/dx+vy*dRHO/dy) */
	leftx-=timestep/2.0*GXKOEF*(vx[v[0]]*(ro[v[2]]+ro[v[3]]-ro[v[0]-ynumy]-ro[v[1]-ynumy])/(gx[m1+1]-gx[m1-1])/2.0+((vy[v[0]]+vy[v[1]])*(gx[m1]-gx[m1-1])+(vy[v[0]-ynumy]+vy[v[1]-ynumy])*(gx[m1+1]-gx[m1]))/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]));
	/**/
	/* Effective NU calc */
	nueff=(nu[v[0]]+nu[v[1]])/2.0;
	/**/
	/* X-STOKS Error */
	leftx=(leftx-rightx)/nueff;
	/**/
	/* Min,Max Value of Nu Save */
	errbuf[9]=nueff;
	/* Min,Max Value of P Save */
	errbuf[10]=MAXV(pr[v[1]],pr[v[3]]);
	errbuf[11]=MINV(pr[v[1]],pr[v[3]]);
	/**/
	return leftx;
	}
/**/
/**/
/**/
/* Set Initial Num of lines -------------------------------------- */
wn[0]=2;
/* Save Right part Save for X-Stokes ---------------------*/
wi[0]=rightx;
/**/
/* Add Coefficients for left parts of X-Stokes ----------------*/
/* Staggered Nodes num */
/*                      [0]                [2]   */
/*                  Ro0,Sxy0,Nu0                 */
/*                                               */
/*         Pr1,Sxx1    <Vx0>  Pr3,Sxx3           */
/*                                               */
/*                      [1]                [3]   */
/*                  RO1,Sxy1,Nu1                 */
/*                                               */
/*  0(P) 1(Vx)  2(Vy)  */
/* -dP/dX */
wn[1]=v[1]*3+0;
wi[1]=+1.0/xkf;
wn[2]=v[3]*3+0;
wi[2]=-1.0/xkf;
/* dSIGxx/dX */
sxxcalc(m1,m2+1,-1.0/xkf);
sxxcalc(m1+1,m2+1,1.0/xkf);
/* dSIGxy/dY */
sxycalc(m1,m2,-1.0/ykf);
sxycalc(m1,m2+1,1.0/ykf);
/* -gx/dt/2*(vx*dRHO/dx+vy*dRHO/dy) */
wn[wn[0]+1]=v[0]*3+1;
wi[wn[0]+1]=-timestep/2.0*GXKOEF*(ro[v[2]]+ro[v[3]]-ro[v[0]-ynumy]-ro[v[1]-ynumy])/(gx[m1+1]-gx[m1-1])/2.0;
wn[wn[0]+2]=v[0]*3+2;
wi[wn[0]+2]=-timestep/2.0*GXKOEF*(gx[m1]-gx[m1-1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
wn[wn[0]+3]=v[1]*3+2;
wi[wn[0]+3]=-timestep/2.0*GXKOEF*(gx[m1]-gx[m1-1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
wn[wn[0]+4]=(v[0]-ynumy)*3+2;
wi[wn[0]+4]=-timestep/2.0*GXKOEF*(gx[m1+1]-gx[m1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
wn[wn[0]+5]=(v[1]-ynumy)*3+2;
wi[wn[0]+5]=-timestep/2.0*GXKOEF*(gx[m1+1]-gx[m1])/(gx[m1+1]-gx[m1-1])/2.0*(ro[v[1]]-ro[v[0]])/(gy[m2+1]-gy[m2]);
wn[0]+=5;
/**/
/**/
/**/
/*
for(n1=0;n1<=vn[0];n1++)printf("Vx %e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left+Right Side or Err for X-Stokes Equation */




/* LEFT+Right Side or Err for Y-Stokes Equation */
/* Stoks equation initial form */
/* dSIGyy/dY + dSIGxy/dX - dP/dY = - RO*Gy */
double ystokserr(long int m1, long int m2, int ynerr)
/* m1,m2 - node X,Y number */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Counters */
int n1;
long int v[4];
/* Err Buf */
double lefty,righty,nueff;
/* Distances */
double xkf=gx[m1+1]-gx[m1],ykf=(gy[m2+1]-gy[m2-1])/2.0;
/**/
/**/
/**/
/* Staggered Nodes num */
/*                                               */
/*                Pr2,Syy2                       */
/*                                               */
/*         [0]                  [2]              */
/*    Ro0,Sxy0,Nu0   <Vy0>   Ro2,Sxy2,Nu2        */
/*                                               */
/*                Pr3,Syy3                       */
/*                                               */
/*         [1]                  [3]              */
/*                                               */
v[0]=m1*ynumy+m2;v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* RIGHT part of Y-Stokes */
/* dSIGik/dXk - dP/dXi = - Gi*RO */
righty  = -GYKOEF*(ro[v[0]]+ro[v[2]])/2.0;
/**/
/**/
/**/
/**/
/* Return val for LSQ err ----------------------------*/
if(ynerr==1)
	{
	/* LEFT part of Y-Stokes */
	/* dSIGyy/dY + dSIGxy/dX - dP/dY = - RO*Gy */
	/**/
	/* dSIGyy/dY */
	lefty =(syy[v[3]]-syy[v[2]])/ykf;
	/* dSIGxy/dX */
	lefty+=(sxy[v[2]]-sxy[v[0]])/xkf;
	/* -dP/dY */
	lefty-=(pr[v[3]]-pr[v[2]])/ykf;
	/* -gy/dt/2*(vy*dRHO/dy+vx*dRHO/dx) */
	lefty-=timestep/2.0*GYKOEF*(vy[v[0]]*(ro[v[1]]+ro[v[3]]-ro[v[0]-1]-ro[v[2]-1])/(gy[m2+1]-gy[m2-1])/2.0+((vx[v[0]]+vx[v[2]])*(gy[m2]-gy[m2-1])+(vx[v[0]-1]+vx[v[2]-1])*(gy[m2+1]-gy[m2]))/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]));
	/**/
	/* Effective NU calc */
	nueff=(nu[v[0]]+nu[v[2]])/2.0;
	/**/
	/* Y-STOKS Error */
	lefty=(lefty-righty)/nueff;
	/**/
	/* Min,Max Value of Nu Save */
	errbuf[9]=nueff;
	/* Min,Max Value of P Save */
	errbuf[10]=MAXV(pr[v[2]],pr[v[3]]);
	errbuf[11]=MINV(pr[v[2]],pr[v[3]]);
	/**/
	return lefty;
	}
/**/
/**/
/**/
/* Set Initial Num of lines -------------------------------------- */
wn[0]=2;
/* Save Right parts Save for Y-Stokes ---------------------*/
wi[0]=righty;
/**/
/* Add Coefficients for left parts of Y-Stokes ----------------*/
/* Staggered Nodes num */
/*                                               */
/*                Pr2,Syy2                       */
/*                                               */
/*         [0]                  [2]              */
/*    Ro0,Sxy0,Nu0   <Vy0>   Ro2,Sxy2,Nu2        */
/*                                               */
/*                Pr3,Syy3                       */
/*                                               */
/*         [1]                  [3]              */
/*                                               */
/*  0(P) 1(Vx)  2(Vy)  */
/* -dP/dY */
wn[1]=v[2]*3+0;
wi[1]=+1.0/ykf;
wn[2]=v[3]*3+0;
wi[2]=-1.0/ykf;
/* dSIGyy/dY */
sxxcalc(m1+1,m2,1.0/ykf);
sxxcalc(m1+1,m2+1,-1.0/ykf);
/* dSIGxy/dX */
sxycalc(m1,m2,-1.0/xkf);
sxycalc(m1+1,m2,1.0/xkf);
/* -gy/dt/2*(vy*dRHO/dy+vx*dRHO/dx) */
wn[wn[0]+1]=v[0]*3+2;
wi[wn[0]+1]=-timestep/2.0*GYKOEF*(ro[v[1]]+ro[v[3]]-ro[v[0]-1]-ro[v[2]-1])/(gy[m2+1]-gy[m2-1])/2.0;
wn[wn[0]+2]=v[0]*3+1;
wi[wn[0]+2]=-timestep/2.0*GYKOEF*(gy[m2]-gy[m2-1])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
wn[wn[0]+3]=v[2]*3+1;
wi[wn[0]+3]=-timestep/2.0*GYKOEF*(gy[m2]-gy[m2-1])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
wn[wn[0]+4]=(v[0]-1)*3+1;
wi[wn[0]+4]=-timestep/2.0*GYKOEF*(gy[m2+1]-gy[m2])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
wn[wn[0]+5]=(v[2]-1)*3+1;
wi[wn[0]+5]=-timestep/2.0*GYKOEF*(gy[m2+1]-gy[m2])/(gy[m2+1]-gy[m2-1])/2.0*(ro[v[2]]-ro[v[0]])/(gx[m1+1]-gx[m1]);
wn[0]+=5;
/**/
/**/
/**/
/*
for(n1=0;n1<=vn[0];n1++)printf("Vy %e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left+Right Side or Err for Y-Stokes Equation */




/* Left side or Err for Continuity Equation  */
/* dVx/dX+dVy/dY = 0 */
double conterr(long int m1, long int m2, int ynerr)
/* m1,m2 - node X,Y number */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Counter */
long int v[4],m1min,m1max,m2min,m2max,m3;
/* Val Buffer */
double divv=0;
/**/
/* Staggered Nodes num */
/*   [0]       Vy0      [2] */
/*                          */
/*   Vx0        <P3>    Vx2 */
/*            Exx3,Eyy3     */
/*                          */
/*   [1]       Vy1      [3] */
v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/*
if(dv[v[3]]!=0){printf("%ld %ld %e",m1,m2,dv[v[3]]);getchar();}
*/
/**/
/* Return dVx/dX+dVy/dY err ----------------------------*/
if(ynerr==1)
	{
	/* dVx/dX+dVy/dY */
	divv=(vx[v[2]]-vx[v[0]])/(gx[m1]-gx[m1-1])+(vy[v[1]]-vy[v[0]])/(gy[m2]-gy[m2-1]);
	/**/
	return divv-dv[v[3]]+prb[v[3]]*(pr[v[3]]-pr0[v[3]])/timestep;
	}
/**/
/**/
/**/
/* Add continuity equation */
/* Save Right part for Contin ---------------------*/
wi[0]=dv[v[3]]+prb[v[3]]*pr0[v[3]]/timestep;
/**/
/**/
/**/
/* Add Coefficients for left parts of dP/dt*dlnRHO/dP+dVx/dX+dVy/dy  ----------------*/
/*  0(P) 1(Vx)  2(Vy)  */
/* Add P with koefficients */
wn[1]=v[3]*3;
wi[1]=prb[v[3]]/timestep;
/* Add Vx with koefficients */
wn[2]=v[0]*3+1;
wi[2]=-1.0/(gx[m1]-gx[m1-1]);
wn[3]=v[2]*3+1;
wi[3]=+1.0/(gx[m1]-gx[m1-1]);
/* Add Vy with koefficients */
wn[4]=v[0]*3+2;
wi[4]=-1.0/(gy[m2]-gy[m2-1]);
wn[5]=v[1]*3+2;
wi[5]=+1.0/(gy[m2]-gy[m2-1]);
/**/
/* Add total Num of lines */
wn[0]=5;
/**/
/**/
/**/
/* Check Boundary conditions around Cell */
divv=1.0;
if (!bondm[v[0]*3+1]) divv=0;
if (!bondm[v[0]*3+2]) divv=0;
if (!bondm[v[1]*3+2]) divv=0;
if (!bondm[v[2]*3+1]) divv=0;
/**/ 
/*
for(n1=0;n1<3;n1++)printf("Cont %e %d \n",wi[n1],wn[n1]);getchar();
*/
return divv;
}
/* Left side or Err for Continuity Equation */



/* Left side or Err for vX Boundary Condition Equation */ 
/* V=CONST+KOEF*Vn */
double xbonderr(long int mcmax, int ynerr)
/* mcmax - numer of cur Vx in sol[] */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Val Buffer */
double leftx=0;
int n1;
/**/
/**/
/**/
/* Error Calc */
if (ynerr)
	{
	/* Add Const */
	leftx=x[mcmax]-bondv[bondm[mcmax]][0];
	/* Add Koef */
	for (n1=0;n1<3;n1++)
		{
		if(bondn[bondm[mcmax]][n1]) 
			{
			leftx-=bondv[bondm[mcmax]][n1+1]*x[bondn[bondm[mcmax]][n1]-1];
			}
		}
	/**/
	return leftx;
	}
/**/
/**/
/**/
/* Add X CONST */
wn[0]=1;
wi[0]=bondv[bondm[mcmax]][0];
wn[1]=mcmax;
wi[1]=1.0;
/* Add X PAR1,PAR2,PAR3 */
for (n1=0;n1<3;n1++)
	{
	if(bondn[bondm[mcmax]][n1]) 
		{
		wn[0]+=1;
		wn[wn[0]]=bondn[bondm[mcmax]][n1]-1;
		wi[wn[0]]=-bondv[bondm[mcmax]][n1+1];
		}
	}
/**/
/*
for(n1=0;n1<3;n1++)printf("%e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Err for vX Boundary Condition Equation */ 





/* Left side or Err for vY Boundary Condition Equation */ 
/* V=CONST+KOEF*Vn */
double ybonderr(long int mcmax, int ynerr)
/* mcmax - numer of cur Vx in sol[] */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Val Buffer */
double lefty=0;
int n1;
/**/
/**/
/**/
/* Error Calc */
if (ynerr)
	{
	/* Add Const */
	lefty=x[mcmax]-bondv[bondm[mcmax]][0];
	/* Add Koef */
	for (n1=0;n1<3;n1++)
		{
		if(bondn[bondm[mcmax]][n1]) 
			{
			lefty-=bondv[bondm[mcmax]][n1+1]*x[bondn[bondm[mcmax]][n1]-1];
			}
		}
	/**/
	return lefty;
	}
/**/
/**/
/**/
/* Add Y CONST */
wn[0]=1;
wi[0]=bondv[bondm[mcmax]][0];
wn[1]=mcmax;
wi[1]=1.0;
/* Add Y PAR1,PAR2,PAR3 */
for (n1=0;n1<3;n1++)
	{
	if(bondn[bondm[mcmax]][n1]) 
		{
		wn[0]+=1;
		wn[wn[0]]=bondn[bondm[mcmax]][n1]-1;
		wi[wn[0]]=-bondv[bondm[mcmax]][n1+1];
		}
	}
/*
for(n1=0;n1<3;n1++)printf("%e %d \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Err for vY Boundary Condition Equation */ 




/* Left side or Err for P Boundary Condition Equation */ 
/* P=CONST+KOEF*Pn */
double pbonderr(long int mcmax, int ynerr)
/* mcmax - numer of cur P  in sol[] */
/* ynerr - Err Calc Y(1)/N(0) */
{
/* Val Buffer */
double leftp;
int n1;
/**/
/**/
/**/
/* Error Calc */
if (ynerr)
	{
	/* Add Const */
	leftp=x[mcmax]-bondv[bondm[mcmax]][0];
	/* Add Koef */
	for (n1=0;n1<3;n1++)
		{
		if(bondn[bondm[mcmax]][n1]) 
			{
			leftp-=bondv[bondm[mcmax]][n1+1]*x[bondn[bondm[mcmax]][n1]-1];
			}
		}
	/**/
	return leftp*leftp;
	}
/**/
/**/
/**/
/* Add P CONST */
wn[0]=1;
wi[0]=bondv[bondm[mcmax]][0];
wn[1]=mcmax;
wi[1]=1.0;
/* Add P PAR1,PAR2,PAR3 */
for (n1=0;n1<3;n1++)
	{
	if(bondn[bondm[mcmax]][n1]) 
		{
		wn[0]+=1;
		wn[wn[0]]=bondn[bondm[mcmax]][n1]-1;
		wi[wn[0]]=-bondv[bondm[mcmax]][n1+1];
		}
	}
/*
for(n1=0;n1<3;n1++)printf("%e %d \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Err for P Boundary Condition Equation */ 





/* Solve XY-Stokes+Continuity equations by vX,vY,P mixed arbitrary order Finite Diff method */
void viterate(int m0)
{
/* Counters */
long int m1,m2,m3,mcmax,mcmax0,mcmax1,dm1,dm2,ccc=0;
double vxmin,vxmax,vymin,vymax,minvx=1e+30,maxvx=-1e+30,minvy=1e+30,maxvy=-1e+30,mindx,pmpa,minnu=1e+50,maxnu=-1e+50,minpr=1e+50,maxpr=-1e+50;
int printyn=printmod,n1;
/**/
/* Val buffer */
double ival,mukoef;
/* Err koef */
double bondsum=0,bondnum=0;
double stoksum=0,stoknum=0;
double contsum=0,contnum=0;
double dvsum,dvnum,dscel,dvmin=1e+50,dvmax=-1e+50;
/**/
/**/
/**/
/* Clear New Solution */
/* P, Vx,Vy */
pos0cur=0;
mukoef=nubeg*strmin;
/**/
/**/
/**/
/* Check continuity condition globally */
dvsum=0;
dvnum=0;
/* Node  Cycle */
for (m1=1;m1<xnumx;m1++)
for (m2=1;m2<ynumy;m2++)
	{
	/* Pos P,Vx,Vy in sol0[] */
	m3=m1*ynumy+m2;
	dscel=(gx[m1]-gx[m1-1])*(gy[m2]-gy[m2-1]);
	dvsum+=dscel*dv[m3];
	dvnum+=dscel;
	dvmin=MINV(dvmin,dv[m3]);
	dvmax=MAXV(dvmax,dv[m3]);
	if(timesum<=0) prb[m3]=0;
/*
printf("%ld %ld %e %e %e %e",m1,m2,pr[m3],pr0[m3],prb[m3],dv[m3]);getchar();
*/
	}
dvsum*=ysize/dvnum;
printf("\n Global Divergence %e %e min = %e max = %e \n",dvsum/ysize,dvsum,dvmin,dvmax);
/**/
/**/
/**/
/* Add Matrix by vX-vY-Stokes, Continuity, Boundary, EPS, SIG, Equations */
/* Node  Cycle */
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	/* Pos P,Vx,Vy in sol0[] */
	mcmax=(m1*ynumy+m2)*3;
/*
printf("P-Vx-Vy Cycle %ld %ld    %ld",m1,m2,mcmax); getchar();
*/
	/**/
	/**/
	/**/
	/* Add Continuity equation for Cells ------------------------------------------------ */
	if(m1 && m2)
		{
/*
printf("Pr %ld %ld    %ld",m1,m2,bondm[mcmax+0]); getchar();
*/
		if(!bondm[mcmax+0]) 
			{
			/* Check Pressure in Cell calculation Y/N */ 
	                conterr(m1,m2,0);
			/* Rescale coefficients */
			for (n1=0;n1<=wn[0];n1++) wi[n1]*=mukoef;
	        	gausmat4(1,mcmax+0,0);
			}
		else
			{
			/* Add P-Boundary */
			pbonderr(mcmax+0,0);
			/* Rescale coefficients */
			for (n1=0;n1<=wn[0];n1++) wi[n1]*=mukoef;
	        	gausmat4(1,mcmax+0,0);
			}
		}
	/* Add Ghost parameters */
	else
		{
		wn[0]=1;
		wi[0]=0;
		wn[1]=mcmax+0;
		wi[1]=mukoef;
	        gausmat4(1,mcmax+0,0);
		}

	/**/
	/**/
	/**/
	/* Add vX-Equations --------------------------------------------------- */
	if(m2<ynumy-1)
		{
/*
printf("Vx %ld %ld    %ld",m1,m2,bondm[mcmax+1]); getchar();
*/
		if(!bondm[mcmax+1] || (timesum>timebond && m1>2 && m2>2 && m1<xnumx-4 && m2<ynumy-3)) 
			{
			/* Add vX-Stokes */
	                xstokserr(m1,m2,0);
			/**/
        		/* Add matrix */
			/* vX */
			gausmat4(1,mcmax+1,0);
/*
printf("Vx %ld %ld    %ld",m1,m2,bondm[mcmax+1]); getchar();
*/
			}
		else
			{
			/* Continuity Equation Vx boundary condition */
			if(bondv[bondm[mcmax+1]][1]>1e+30)
				{
				/* m1 m2 increment definition */
				m3=(bondn[bondm[mcmax+1]][0]-1-(mcmax+1))/3;
				dm1=dm2=0;
				if(m3>=ynumy) dm1=1;
				if(m3==1 || m3>ynumy) dm2=1;
/*
printf("Vx(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+1,bondn[bondm[mcmax+1]][0]-2,m1,m2); getchar();
*/
				if(m3<0 || m3>ynumy+1 || !(m1+dm1) || !(m2+dm2) || !bondm[bondn[bondm[mcmax+1]][0]-2])
					{
					printf("EXIT PROGRAM Inconsistent Vx(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+1,bondn[bondm[mcmax+1]][0]-2,m1,m2);
					exit(0);
					}
		                conterr(m1+dm1,m2+dm2,0);
				}
			else
				{
				/* Add vX Simple Boundary */
				xbonderr(mcmax+1,0);
				}
			/**/
			/* Vx boundary condition Add */
			/* Rescale coefficients */
			for (n1=0;n1<=wn[0];n1++) wi[n1]*=mukoef;
			gausmat4(1,mcmax+1,0);
			}
		}
	/* Add Ghost parameters */
	else
		{
		wn[0]=1;
		wi[0]=0;
		wn[1]=mcmax+1;
		wi[1]=mukoef;
	        gausmat4(1,mcmax+1,0);
		}
	/**/
	/**/
	/**/
	/* Add vY-Equations --------------------------------------------------- */
	if(m1<xnumx-1)
		{
/*
printf("Vy %ld %ld    %ld",m1,m2,bondm[mcmax+2]); getchar();
*/
		if(!bondm[mcmax+2] || (timesum>timebond && m1>2 && m2>2 && m1<xnumx-4 && m2<ynumy-3)) 
			{
			/* Add vY-Stokes */
        	        ystokserr(m1,m2,0);
			/**/
        		/* Add matrix */
			/* vY */
			gausmat4(1,mcmax+2,0);
/*
printf("Vy2 %ld %ld    %ld",m1,m2,bondm[mcmax+8]); getchar();
*/
			}
		else
			{
			/* Continuity Equation Vy boundary condition */
			if(bondv[bondm[mcmax+2]][1]>1e+30)
				{
				/* m1 m2 increment definition */
				m3=(bondn[bondm[mcmax+2]][0]-1-(mcmax+2))/3;
				dm1=dm2=0;
				if(m3>=ynumy) dm1=1;
				if(m3==1 || m3>ynumy) dm2=1;
				if(m3<0 || m3>ynumy+1 || !(m1+dm1) || !(m2+dm2) || !bondm[bondn[bondm[mcmax+2]][0]-3])
/*
printf("Vy(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+2,bondn[bondm[mcmax+2]][0]-3,m1,m2); getchar();
*/
					{
					printf("EXIT PROGRAM Inconsistent Vy(%ld) P(%ld) boundary conditions for X=%ld Y=%ld Node",mcmax+2,bondn[bondm[mcmax+2]][0]-3,m1,m2);
					exit(0);
					}
	        	        conterr(m1+dm1,m2+dm2,0);
				}
			/* Simple Vy boundary condition */
			else
				{
/*
				if(m2==0) bondv[bondm[mcmax+2]][0]=-dvsum;
				if(m2==ynumy-1) bondv[bondm[mcmax+2]][0]=dvsum;
				if(m2==0) bondv[bondm[mcmax+2]][0]=-dvsum;
				if(m2==ynumy-1) {printf("%ld %ld %e %e %e",m1,m2,bondv[bondm[mcmax+2]][0],bondv[bondm[mcmax+2]][1],bondv[bondm[mcmax+2]][2]);getchar();}
				if(m2==ynumy-1) bondv[bondm[mcmax+2]][1]=0;
*/
				/* Add vY Simple Boundary */
				ybonderr(mcmax+2,0);
				/* Divergence improve */
				if(m2==0) wi[0]-=dvsum;
				}
			/**/
			/* Vy boundary condition Add */
			/* Rescale coefficients */
			for (n1=0;n1<=wn[0];n1++) wi[n1]*=mukoef;
			gausmat4(1,mcmax+2,0);
			}
		}
	/* Add Ghost parameters */
	else
		{
		wn[0]=1;
		wi[0]=0;
		wn[1]=mcmax+2;
		wi[1]=mukoef;
	        gausmat4(1,mcmax+2,0);
		}
	/**/
	/**/
	/**/
	}
/* End  Add Matrix By vX-vY-Stokes, Continuity Equations */
/**/
/**/
/**/
/* Solve Matrix ================================================ */
if (printmod) printf("Number of positions in global matrix = %ld  Number of unknown = %ld \n",pos0cur,nodenum3);
gausmat4(0,nodenum3,0);
/* Solve Matrix ================================================ */
/**/
/**/
/**/
/* Reload P, Vx, Vy Results */
/* Pressure correction */
m1=1;m2=1;
ival=pinit-x[(m1*ynumy+m2)*3];
/* Node  Cycle */
for (m1=0;m1<xnumx;m1++)
	{
	/* Set Initial p value at upper boundary */
	pmpa=pinit;
	/**/	
	for (m2=0;m2<ynumy;m2++)
		{
		/* Pos P,Vx,Vy in sol0[] */
		mcmax0=(m1*ynumy+m2)*3;
		/* Pos in vx[], vy[], pr[], etc. */
		mcmax1=m1*ynumy+m2;
		/**/	
		/**/	
		/**/	
		/* Reload/Recalc P */	
		if(m1 && m2) 
			{
			pr[mcmax1]=pmpa=x[mcmax0+0]+ival;
/*	
printf("\n %ld %e",mcmax1,pr[mcmax1]); getchar();
*/	
			}
		/**/	
		/**/	
		/**/	
		/* Reload Vx */	
		if(m2<ynumy-1)
			{ 
			vx[mcmax1]=x[mcmax0+1];
			}
		/**/	
		/**/	
		/**/	
		/* Reload Vy */	
		if(m1<xnumx-1)
			{
			vy[mcmax1]=x[mcmax0+2];
			}
/*	
printf("\n %ld %ld %ld %e %e",m1,m2,mcmax1,vx[mcmax1],vy[mcmax1]); getchar();
*/	
		}
	}
/* End Reload P, Vx, Vy Results */
/**/	
/**/	
/**/	
/* Recalc EPS, SIG Results */
/* Node  Cycle */
for (m1=1;m1<xnumx;m1++)
for (m2=1;m2<ynumy;m2++)
	{
	/* Pos in vx[], vy[], pr[], etc. */
	mcmax1=m1*ynumy+m2;
	/**/	
	/**/	
	/**/	
	/* Sxx,Exx */	
	sxx[mcmax1]=sxxcalc(m1,m2,0); exx[mcmax1]=eps[0];
	/**/	
	/* Syy,Eyy */	
	syy[mcmax1]=-sxx[mcmax1]; eyy[mcmax1]=-exx[mcmax1];
	/**/	
	/* Sxy,Exy */	
	if(m1<xnumx-1 && m2<ynumy-1)
		{
		sxy[mcmax1]=sxycalc(m1,m2,0); exy[mcmax1]=eps[0];
		}
	}
/* End Recalc EPS, SIG Results */
/**/	
/**/	
/**/	
/* Vx,Vy max-min definition */
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	/* Pos of Vx in sol0[] */
	m3=m1*ynumy+m2;
	/**/
	/* Min,Max Vx definition */
	if(m2<ynumy-1)
		{
		minvx=MINV(minvx,vx[m3]);
		maxvx=MAXV(maxvx,vx[m3]);
		}
	/* Min,Max Vy definition */
	if(m1<xnumx-1)
		{
		minvy=MINV(minvy,vy[m3]);
		maxvy=MAXV(maxvy,vy[m3]);
		}
	}
/* Max Vx,Vy Diff in Grid Calc */
vxmin=minvx-MAXV(ABSV(vyfluid),ABSV(vymelt));
vymin=minvy-MAXV(ABSV(vyfluid),ABSV(vymelt));
vxmax=maxvx;
vymax=maxvy;
maxvx-=minvx;
maxvy-=minvy;
maxvx=MAXV(maxvx,maxvy);
mindx=MINV(xstpx,ystpy);
/**/
/**/
/**/
/* Check Error */
/* Node  Cycle */
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	/* Pos of P,Vx,Vy in sol0[] */
	mcmax0=(m1*ynumy+m2)*3;
	/**/
	/**/
	/**/
	/* Check Continuity equation for Cells =========================== */
	if(m1 && m2) 
		{
		ival=conterr(m1,m2,1);
/*	
printf("\n %ld %ld   %e %e %e     %e ",m1,m2,ival,maxvx,mindx,ival/(maxvx/mindx)); getchar();
*/	
		contsum+=ival*ival;
       	        contnum+=1.0;
		/* Print Results */
		ival/=maxvx/mindx;
		if (printmod && ABSV(ival)>DIVVMIN)
			{
			printf("\n Large Continuity err at X=%ld Y=%ld:   Err=%e",m1,m2,ival);
			}
		}
	/**/
	/**/
	/**/
	/* Check vX-Equations for nodes =========================== */
	if(m2<ynumy-1)
		{
		if(!bondm[mcmax0+1] || (timesum>timebond && m1>2 && m2>2 && m1<xnumx-4 && m2<ynumy-3)) 
			{
			/* Add vX-Stokes */
	                ival=xstokserr(m1,m2,1);
        	        stoksum+=ival*ival;
                	stoknum+=1.0;
	                /* Min,Max Nu value Calc */
        	        minnu=MINV(minnu,errbuf[9]);
	                maxnu=MAXV(maxnu,errbuf[9]);
        	        /* Min,Max Pr value Calc */
                	maxpr=MAXV(maxpr,errbuf[10]);
	                minpr=MINV(minpr,errbuf[11]);
			/* Print Results */
			ival/=maxvx/mindx/mindx;
			if (printmod && ABSV(ival)>STOKSMIN)
				{
				printf("\n Large X stokes err at X=%ld Y=%ld:   Err=%e",m1,m2,ival);
				}
	                }
		else
			{
			if(bondv[bondm[mcmax0+1]][1]<1e+30)
				{
				/* Add vX-Boundary */
				ival=xbonderr(mcmax0+1,1);
				bondsum+=ival*ival;
	                	bondnum+=1.0;
				}
			}
		}
	/**/
	/**/
	/**/
	/* Check vY-Equations for nodes =========================== */
	if(m1<xnumx-1)
		{
		if(!bondm[mcmax0+2] || (timesum>timebond && m1>2 && m2>2 && m1<xnumx-4 && m2<ynumy-3)) 
			{
			/* Add vX-vY-Stokes */
	                ival=ystokserr(m1,m2,1);
        	        stoksum+=ival*ival;
                	stoknum+=1.0;
	                /* Min,Max Nu value Calc */
        	        minnu=MINV(minnu,errbuf[9]);
	                maxnu=MAXV(maxnu,errbuf[9]);
        	        /* Min,Max Pr value Calc */
	                maxpr=MAXV(maxpr,errbuf[10]);
        	        minpr=MINV(minpr,errbuf[11]);
			/* Print Results */
			ival/=maxvx/mindx/mindx;
			if (printmod && ABSV(ival)>STOKSMIN)
				{
				printf("\n Large Y stokes err at X=%ld Y=%ld:   Err=%e",m1,m2,ival);
				}
                	}
		else
			{
			if(bondv[bondm[mcmax0+2]][1]<1e+30)
				{
				/* Add vX-vY-Boundary */
				ival=ybonderr(mcmax0+2,1);
				bondsum+=ival*ival;
                		bondnum+=1.0;
				}
			}
		}
/*	
printf("\n %ld %ld   %e %e %e ",m1,m2,minnu,maxnu,errbuf[9]); getchar();
*/	
	}
stoksum=pow(stoksum/stoknum,0.5)/(maxvx/mindx/mindx);
contsum=pow(contsum/contnum,0.5)/(maxvx/mindx);
bondsum=pow(bondsum/bondnum,0.5)/maxvx;
/* End Check Error */
/**/
/**/
/**/
/* Print Results */
if (printmod)
	{
	 printf("\n KRUG %2d \n",m0+1);
	 printf("X-VELOCITY: min = %e max = %e \n",vxmin,vxmax);
	 printf("Y-VELOCITY: min = %e max = %e \n",vymin,vymax);
	 printf("PRESSURE: min = %e max = %e \n",minpr,maxpr);
	 printf("VISKOS: min = %e max = %e \n",minnu,maxnu);
	 printf("STOKS: num = %e err = %e \n",stoknum,stoksum);
	 printf("CONT : num = %e err = %e \n",contnum,contsum);
	 printf("BOUND V: num = %e err = %e \n",bondnum,bondsum);
	 }
/**/
/* Printf EPS, SIG Results */
/*
for (m1=100;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	allintere(gx[m1],gy[m2]);
	ival=pow(0.5*(eps[6]*eps[6]+eps[8]*eps[8])+eps[4]*eps[4],0.5);
printf("%ld %ld %e %e   %e   %e %e %e   %e",m1,m2,gx[m1],gy[m2],eps[10],eps[6],eps[8],eps[4],ival);getchar();
	}
*/	
/**/	
/**/
/**/
/**/
}
/* End Solve XY-Stokes+Continuity equations by vX,vY,P mixed arbitrary order Finite Diff method */





/* Max Vx,Vy in nodes serch time step recalc */
void maxvelstep()
{
double maxvx=0,maxvy=0,ival;
long int m1,m2,m3;
/**/
/**/
/**/
/* Node  Cycle */
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
	{
	/* Pos of in Vx,Vy */
	m3=m1*ynumy+m2;
	/**/
	if(m2<ynumy-1)
		{
		ival=ABSV(vx[m3]);
		maxvx=MAXV(maxvx,ival);
		}
	if(m1<xnumx-1)
		{
		ival=ABSV(vy[m3]);
		maxvy=MAXV(maxvy,ival);
		}
	}
if (maxvx)
	{
	maxvx=(maxxystep*xstpx)/maxvx;
	if(printmod) printf("\n !!! MAX VALID TIME STEP FOR Vx-MARKER %e YEAR !!!\n",maxvx/3.15576e+7);
	timestep=MINV(maxvx,timestep);
	}
if (maxvy)
	{
	maxvy=(maxxystep*ystpy)/maxvy;
	if(printmod) printf("\n !!! MAX VALID TIME STEP FOR Vy-MARKER %e YEAR !!!\n",maxvy/3.15576e+7);
	timestep=MINV(maxvy,timestep);
	}
}
/* Max Vx,Vy in nodes serch time step recalc */




/* Weight of FD calculation for after Fornberg (1996) */
void fdweight(int n, int m, double xi)
/* n - maximal index 0-n */
/* m - required derivative order 0-m */
/* xi - derivation point coordinate */
{
/* Counters */
int i,j,k,mn;
double c1,c2,c3,c4,c5,kk;
/**/
/**/
/**/
c1=1.0;
c4=xn[0]-xi;
for(k=0;k<=m;k++)
	{
	for(j=0;j<=n;j++)
		{
		cn[j][k]=0;
		}
	}
/**/
cn[0][0]=1.0;
for(i=1;i<=n;i++)
	{
	mn=i;if(mn>m) mn=m;
	c2=1.0;
	c5=c4;
	c4=xn[i]-xi;
	for(j=0;j<i;j++)
		{
		c3=xn[i]-xn[j];
		c2*=c3;
		for(k=mn;k>0;k--)
			{
			kk=(double)(k);
			cn[i][k]=c1*(kk*cn[i-1][k-1]-c5*cn[i-1][k])/c2;
			}
		cn[i][0]=-c1*c5*cn[i-1][0]/c2;
		for(k=mn;k>0;k--)
			{
			kk=(double)(k);
			cn[j][k]=(c4*cn[j][k]-kk*cn[j][k-1])/c3;
			}
		cn[j][0]=c4*cn[j][0]/c3;
		}
	c1=c2;
	}
/*
for(i=0;i<=n;i++)printf("FD %d %d %e %d %e %e\n",n,m,xi,i,xn[i],cn[i][m]);getchar();
*/
}
/* Weight of FD calculation after Fornberg (1996) */


/* Left side or Value for Sxx  Equation */ 
/* Sxx=2Nu*Exx*X+Sxx0*(1-X), Exx=1/2(dVx/dX-dVy/dY) */
double sxxcalc(long int m1, long int m2, double ynval)
/* m1,m2 - node X,Y number */
/* ynval - Val Sxx Calc Y(0)/N(koefficient) */
{
/* Exx horisontal position */
double xi=(gx[m1-1]+gx[m1])/2.0,leftsxx=0,nueff,sxxeeff,xelvis=1.0,ggeff=0;
long int v[4];
int n1,n;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   [0]      Vy0       [2] */
/*   Nu0                Nu2 */
/*                          */
/*   Vx0    Sxx3,Exx3   Vx2 */
/*                          */
/*   [1]                [3] */
/*   Nu1      Vy1       Nu3 */
/*                          */
v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* Effective viscosity calc */
nueff=nd[v[3]];
/**/
/**/
/**/
/* Nu Save */
eps[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   [0]      Vy0       [2] */
/*   Nu0                Nu2 */
/*                          */
/*   Vx0    Sxx3,Exx3   Vx2 */
/*                          */
/*   [1]                [3] */
/*   Nu1      Vy1       Nu3 */
/*                          */
/**/
/* Return Sxx,Exx val ----------------------------*/
if(ynval==0)
	{
	/* Exx=1/2(dVx/dX-dVy/dY) */
	leftsxx=0.5*((vx[v[2]]-vx[v[0]])/(gx[m1]-gx[m1-1])-(vy[v[1]]-vy[v[0]])/(gy[m2]-gy[m2-1]));
	/**/
	/* Save Exx */
	eps[0]=leftsxx;
	/**/
	/* Calc Sxx=2Nu*Exx*X+Sxx0*(1-X) */
	leftsxx=2.0*nueff*leftsxx;
	/**/
	return leftsxx;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Sxx ----------------*/
/*  0(P) 1(Vx)  2(Vy)  */
/* Sxx=2Nu*Exx*X+Sxx0*(1-X), Exx=1/2(dVx/dX-dVy/dY) */
/* Add Vx with koefficients */
wn[wn[0]+1]=v[0]*3+1;
wi[wn[0]+1]=-ynval*nueff/(gx[m1]-gx[m1-1]);
wn[wn[0]+2]=v[2]*3+1;
wi[wn[0]+2]=+ynval*nueff/(gx[m1]-gx[m1-1]);
/* Add Vy with koefficients */
wn[wn[0]+3]=v[0]*3+2;
wi[wn[0]+3]=+ynval*nueff/(gy[m2]-gy[m2-1]);
wn[wn[0]+4]=v[1]*3+2;
wi[wn[0]+4]=-ynval*nueff/(gy[m2]-gy[m2-1]);
/**/
/* Add total Num of lines */
wn[0]+=4;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Exx %e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Sxx  Equation */ 




/* Left side or Value for Sxy  Equation */ 
/* Sxy=2Nv*Exy, Exy=1/2(dVx/dY+dVy/dX) */
double sxycalc(long int m1, long int m2, double ynval)
/* m1,m2 - node X,Y number */
/* ynval - Val Syy Calc Y(0)/N(koefficient) */
{
/* Exy position */
double xi,leftsxy=0,leftsxy1=0,leftsxy2=0,nueff,sxyeeff,xelvis=1.0,ggeff=0;
long int v[4];
int n1;
/**/
/**/
/**/
/* Staggered Nodes num */
/*  [0]                [2]            */
/*                                    */
/*                     Vx2            */
/*                                    */
/*  [1]     Vy1        [3]       Vy3  */
/*                   Exy3,Nu3         */
/*                                    */
/*                     Vx3            */
/*                                    */
v[0]=(m1-1)*ynumy+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+ynumy;v[3]=v[2]+1;
/**/
/**/
/**/
/* Effective viscosity calc */
nueff=nu[v[3]];
/**/
/**/
/**/
/* Nu Save */
eps[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*  [0]                [2]            */
/*                                    */
/*                     Vx2            */
/*                                    */
/*  [1]     Vy1        [3]       Vy3  */
/*                   Exy3,Nu3         */
/*                                    */
/*                     Vx3            */
/* Return Sxy,Exy val ----------------------------*/
if(ynval==0)
	{
	/* Exy=1/2(dVx/dY+dVy/dX)=0 */
	leftsxy1=(vx[v[3]]-vx[v[2]])/(gy[m2+1]-gy[m2-1]);
	leftsxy2=(vy[v[3]]-vy[v[1]])/(gx[m1+1]-gx[m1-1]);
	/**/
	/* Save Exy */
	eps[0]=leftsxy=leftsxy1+leftsxy2;
	/* Save Esp (rotation rate) */
/*
	eps[1]=leftsxy1-leftsxy2;
*/
	/**/
	/* Calc Sxy=2Nu*Exy-Sxyp */
	leftsxy=2.0*nueff*leftsxy;
	/**/
	return leftsxy;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Sxy ----------------*/
/*  0(P) 1(Vx)  2(Vy)  */
/* Sxy=2Nu*Exy, Exy=1/2(dVx/dY+dVy/dX) */
/* Add Vx with koefficients */
wn[wn[0]+1]=v[2]*3+1;
wi[wn[0]+1]=-ynval*2.0*nueff/(gy[m2+1]-gy[m2-1]);
wn[wn[0]+2]=v[3]*3+1;
wi[wn[0]+2]=+ynval*2.0*nueff/(gy[m2+1]-gy[m2-1]);
/* Add Vy with koefficients */
wn[wn[0]+3]=v[1]*3+2;
wi[wn[0]+3]=-ynval*2.0*nueff/(gx[m1+1]-gx[m1-1]);
wn[wn[0]+4]=v[3]*3+2;
wi[wn[0]+4]=+ynval*2.0*nueff/(gx[m1+1]-gx[m1-1]);
/**/
/* Add total Num of lines */
wn[0]+=4;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Exy %e %ld \n",wi[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Sxy  Equation */ 




