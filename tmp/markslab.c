/* Move markers by using Runge-Kutta method */
void movemark()
{
/* Vx, Vy buffer */
double vx0,vx1,vx2,vx3,vx4,vy0,vy1,vy2,vy3,vy4,vxwater,vywater;
/**/
long int mm1,marknum1,m1,m2,m3;
/* Erosion-Sedimentation Y/N */
int n1;
int mm2;
/* Nonstabilyty for immobile markers */
double dx,dy,ival,dpdx,dpdy,e,n,vxkoef,vykoef,mtk,mpb,mnu;
/**/
/**/
/**/
/* Clear divergence on markers */
for (mm1=0;mm1<marknum;mm1++) markdi[mm1]=0;
/**/
/* Hydration front progress, melt extraction  */
if(timesum>1e+11)
	{
	for (mm1=0;mm1<marknum;mm1++) markv[mm1]=0;
	hydration2();

	meltextract();

	}
/* Save number of markers */
marknum1=marknum;
/**/
/**/
/**/
/* Move markers */
for (mm1=0;mm1<marknum;mm1++)
{
/* Marker type */
mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
if( ((markx[mm1]>=0 && marky[mm1]>=0 && (double)(markx[mm1])<=xsize && (double)(marky[mm1])<=ysize) || outgrid!=1) && !markim[mm2] )
	{
	/**/
	/* Erosion-Sedimentation */
	if((double)(marky[mm1])<=eroslev) n1=1; else n1=0;
	/**/
	/* Calc/Save strain rate and srtress invariants, viscosity */
	allintere((double)(markx[mm1]),(double)(marky[mm1]));
	markeii=pow(0.5*(eps[6]*eps[6]+eps[8]*eps[8])+eps[4]*eps[4],0.5);
	marksii=pow(0.5*(eps[7]*eps[7]+eps[9]*eps[9])+eps[5]*eps[5],0.5);
	/* Correct strain rate for numerical diffusion */
	markrat=1.0;	
	if(0==0 && markn[mm1]>0) markrat*=pow(0.5*marksii/markeii/markn[mm1],0.87);
	mpb=eps[10]*1e-5;
	mtk=(double)(markk[mm1]);
	markpr[mm1]=eps[10];
	/* Compute viscosity */
	if(markt[mm1]<50 && mtk>0)
		{
		/* Save viscosity */
		mnu=viscalc(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1,mm2);
		markn[mm1]=mnu;
		/* Add plastic strain */
		if(marke[mm1]>0)
			{
			marke[mm1]+=timestep*markeii;
/*
*/
			}
		}
	/**/
	/**/
	/* Water marker move */
	vxwater=vywater=0;
	if(markt[mm1]>=50 && markt[mm1]<100) 
		{
		/* Water velocity */
		vywater=vyfluid; if(markd[mm1]>1100.0) vywater=vymelt;
		/* Fluid in rock */
		if(vyfluid>0 && (markk[mm1]==0 || markk[mm1]>298.0)) 
			{
			/* Horizontal,Vertical P-cell index */
			m1=wn[0]; if(markx[mm1]>(gx[m1]+gx[m1+1])/2.0) m1+=1;
			if(m1<1) m1=1; if(m1>xnumx-2) m1=xnumx-2;
			m2=wn[1]; if(marky[mm1]>(gy[m2]+gy[m2+1])/2.0) m2+=1;
			if(m2<1) m2=1; if(m2>ynumy-2) m2=ynumy-2;
			/* Pressure gradients */
			e=(markx[mm1]-(gx[m1-1]+gx[m1])/2.0)/((gx[m1+1]-gx[m1-1])/2.0);
			n=(marky[mm1]-(gy[m2-1]+gy[m2])/2.0)/((gy[m2+1]-gy[m2-1])/2.0);
			m3=m1*ynumy+m2;
			dpdx=2.0*((1.0-n)*(pr[m3+ynumy]-pr[m3])+n*(pr[m3+ynumy+1]-pr[m3+1]))/(gx[m1+1]-gx[m1-1]);
			dpdy=2.0*((1.0-e)*(pr[m3+1]-pr[m3])+e*(pr[m3+ynumy+1]-pr[m3+ynumy]))/(gy[m2+1]-gy[m2-1]);
			/* Recalc velocity koefficients */
			vxkoef=(1000.0*GXKOEF-dpdx)/(2300.0*9.81);
			vykoef=(1000.0*GYKOEF-dpdy)/(2300.0*9.81);
/*
printf("%ld %ld %e %e   %ld %d %e %e %e %e  %e %e %e %e",m1,m2,gx[m1],gy[m2],mm1,markt[mm1],markx[mm1],marky[mm1],e,n,dpdx,dpdy,vxkoef,vykoef);getchar();
*/
			if(vxkoef>2.0) vxkoef=2.0; if(vxkoef<-2.0) vxkoef=-2.0;
			if(vykoef>2.0) vykoef=2.0; if(vykoef<-2.0) vykoef=-2.0;
			/* Recalc velocity */
			vxwater=vywater*vxkoef;
			vywater*=vykoef;
			}
		else
		/* Fluid in water */
			{
			vxwater=0;
			vywater=-ABSV(vywater);
			}
/*
printf("%ld %ld %e %e   %ld %d %e %e %e %e  %e %e  %e %e",m1,m2,gx[m1],gy[m2],mm1,markt[mm1],markx[mm1],marky[mm1],dpdx,dpdy,vxkoef,vykoef,vxwater,vywater);getchar();
*/
		/**/
		}
	/**/
	/* Motion Calc ///////////////////////////////// */
	/**/
	if(markmod==1)
		{
		/* Vx, Vy simple calc */
		allinterv((double)(markx[mm1]),(double)(marky[mm1]));
		vx0=eps[11]+vxwater; vy0=eps[12]+vywater;
		}
	else
		{
		/* Vx, Vy 4 Runge-Kutta koef calc */
		allinterv((double)(markx[mm1]),(double)(marky[mm1]));
		vx1=eps[11]+vxwater; vy1=eps[12]+vywater;
		/**/
		allinterv((double)(markx[mm1])+vx1*timestep/2.0,(double)(marky[mm1])+vy1*timestep/2.0);
		vx2=eps[11]+vxwater; vy2=eps[12]+vywater;
		/**/
		allinterv((double)(markx[mm1])+vx2*timestep/2.0,(double)(marky[mm1])+vy2*timestep/2.0);
		vx3=eps[11]+vxwater; vy3=eps[12]+vywater;
		/**/
		allinterv((double)(markx[mm1])+vx3*timestep,(double)(marky[mm1])+vy3*timestep);
		vx4=eps[11]+vxwater; vy4=eps[12]+vywater;
		/**/
		/* Vx,Vy, EpsXX, EpsYY, EpsXY calc after Runge-Kutta */
		vx0=(vx1+2.0*vx2+2.0*vx3+vx4)/6.0;
		vy0=(vy1+2.0*vy2+2.0*vy3+vy4)/6.0;
		}
	/**/
	/* Orthogonal motion only */
	if (outgrid==2)
		{
		if(markx[mm1]<0 || (double)(markx[mm1])>xsize) vy0=0;		
		if(marky[mm1]<0 || (double)(marky[mm1])>ysize) vx0=0;		
		}
	/**/
	/**/
	/**/
	/* Normal/Immobile markers */
	if(markt[mm1]<100)
		{
		/* Markers coming from the depth */
		if(marky[mm1]>zdeep && vy0<0 && markk[mm1]<tdeep && markt[mm1]==10) 
			{
			markk[mm1]=tdeep;
			}
		if(marky[mm1]>zdeep && vy0<0 && markt[mm1]!=10) 
			{
			markt[mm1]=10;
			markwt[mm1]=0;
			markex[mm1]=0;
			marke[mm1]=0;
			}
/*
		if(marky[mm1]>zdeep && vy0<0 && markk[mm1]<tdeep) markk[mm1]=tdeep;
*/
		/* Normal markers */
		/* X,Y calc after Runge-Kutta */
		markx[mm1]+=(float)(timestep*vx0);
		marky[mm1]+=(float)(timestep*vy0);

		/* No water at the large depth */
		if(marky[mm1]>zmpor+waterlev)
			{
			/* Erase fluid marker */
			if(markt[mm1]>=50) {markx[mm1]=-xsize; markk[mm1]=0;}
			/* Erase water from the solid marker */
			if(markt[mm1]<50 && markw[mm1]>0) markw[mm1]=0;
			}

		/**/
		/* Out of grid marker reset */
		if((markx[mm1]<0 || (double)(markx[mm1])>xsize || (double)(marky[mm1])>ysize)) 
/*
		if((markx[mm1]<0 || marky[mm1]<0 || (double)(markx[mm1])>xsize || (double)(marky[mm1])>ysize)) 
*/
			{
			markk[mm1]=0;
			}
		}
	else
		{
		/* Immobile markers */
		/* X,Y calc after Runge-Kutta */
		markx[mm1]+=(float)(timestep*vx0);
		marky[mm1]+=(float)(timestep*vy0);
		/**/
/*
printf("%ld %d %e %e",mm1,markt[mm1],xnew,ynew); getchar();
*/
		/* Check old position */
		dx=markx[mm1]-markk[mm1];
		dy=marky[mm1]-markn[mm1];
		dy=pow(dx*dx+dy*dy,0.5);
		/* Check new position, add marker */
		if(dy>ystpy && markx[mm1]>=0 && marky[mm1]>=0 && markx[mm1]<=xsize && marky[mm1]<=ysize)
			{
			/* Type save */
			markt[marknum1]=markt[mm1]-100;
			/* X,Y calc after Runge-Kutta */
			markx[marknum1]=markx[mm1];
			marky[marknum1]=marky[mm1];
			/* Temperature, Density, Water, Reset */
			markk[marknum1]=0;
			markdt[marknum1]=0;
			markpr[marknum1]=0;
			markn[marknum1]=0;
			markd[marknum1]=-1.0;
			markw[marknum1]=-1.0;
			markdi[marknum1]=0;
			markwt[marknum1]=0;
			markex[marknum1]=0;
			marke[marknum1]=0;
			markv[marknum1]=0;
			/* Add aditional markers counter */
			marknum1++;
			/* X,Y reset for immobile marker */
			markx[mm1]=markk[mm1];
			marky[mm1]=markn[mm1];
			}
		/* Check,Reset old position */
		if(dy>ystpy || (marky[mm1]<0 && vy0<0) || (marky[mm1]>ysize && vy0>0) || (markx[mm1]<0 && vx0<0) || (markx[mm1]>xsize && vx0>0))
			{
			/* X,Y reset for immobile marker */
			markx[mm1]=markk[mm1];
			marky[mm1]=markn[mm1];
			}
		/**/
		}
	/**/
	/**/
	/**/
	/* Motion Calc ///////////////////////////////// */
	}
}
/*
printf("\n A Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1-1);
*/
/**/
/**/
/**/
/* Mark num */
if(marknum1>MAXMRK) {printf("Space out in markx[]"); exit(0);}
/**/
/**/
/**/
/* Reset aditional markers */
mm1=0;
while(marknum1>marknum && mm1<marknum)
	{
	/* Reload marker */
	if((markx[mm1]<0 || marky[mm1]<0 || (double)(markx[mm1])>xsize || (double)(marky[mm1])>ysize) && markt[mm1]<100 && markt[mm1]>1) 
		{
		/* Decrease aditional markers counter */
		marknum1--;
		/* Reload all */
		markt[mm1]=markt[marknum1];
		markk[mm1]=markk[marknum1];
		markdt[mm1]=markdt[marknum1];
		markpr[mm1]=markpr[marknum1];
		markn[mm1]=markn[marknum1];
		markd[mm1]=markd[marknum1];
		markw[mm1]=markw[marknum1];
		markx[mm1]=markx[marknum1];
		marky[mm1]=marky[marknum1];
		markdi[mm1]=markdi[marknum1];
		markwt[mm1]=markwt[marknum1];
		markex[mm1]=markex[marknum1];
		marke[mm1]=marke[marknum1];
		markv[mm1]=markv[marknum1];
		}
	/* Increase markers counter */
	mm1++;
	}
printf("\n Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1-1);
/* Set new marker number */
marknum=marknum1;
/**/
/**/
/**/
}
/* End Move markers by using Runge-Kutta method */





/* ro[],nu[] recalc after marker positions */
void ronurecalc()
{
/* Counters */
long int m1,m2,m3,m1min,m1max,m2min,m2max;
int mm2,yn,mm3,n1,n2;
long int mm1;
double swt,swt1,swt2,ddx,ddy,celdx,celdy;
double mnu,mro,mcp,mkt,mht,maa,mbb,mab,mab0,maa0,mbb0,mcp0,mkt0,mwa,dmwa,mro0,mroold,ival,ival0,xmelt;
double anu,aro,acp,akt,aht,abb,markwlev;
double wnu,wro,wcp,wkt,wht,wbb,dywa;
/* TD Database variables */
double H0,H1,H2,H3,W0,W1,W2,W3,R0,R1,R2,R3,n,e,hydryl;
/* RO, NU equations var */
double mpb=1.0,mtk=300.0,numax=0,numin=0,markvol=mardx*mardy;
/**/
/**/
/**/
printf("\n Number of nodes = %ld  Number of markers = %ld \n",nodenum,marknum);
/**/
/* Water level band */
m2=m2serch(waterlev);
if (m2<1) m2=1;
if (m2>ynumy-2) m2=ynumy-2;
dywa=(gy[m2+1]-gy[m2-1])/4.0;
/**/
/* Surface changes */
if (erosmod && timestep) erosion();
/* Layering on sediments */
m1=(long int)(timesum/sedimcyc);
m2=((long int)(m1/2))*2;
if(m2==m1) yn=3; else yn=4;
/**/
/**/
/**/
/* Hydration front progress  */
printf("timestep %e timesum %e densimod %d\n",timestep,timesum,densimod);
/**/
/**/
/**/
/* Gridchange */
/*
gridchange();
*/
/**/
/**/
/**/
/* ADD MARKERS TO THE v-CELLS ========================== */
/* Clear ro[],nu[] wt */
for (m1=0;m1<nodenum;m1++)
	{
	ro0[m1]=0;
	et0[m1]=0;
	pr0[m1]=0;
	prb0[m1]=0;
	dv[m1]=0;
	nu0[m1]=0;
	nd0[m1]=0;
	cp0[m1]=0;
	kt0[m1]=0;
	ht0[m1]=0;
	tk0[m1]=0;
	sol0[m1]=0;
	sol1[m1]=0;
	sol0[nodenum+m1]=0;
	sol1[nodenum+m1]=0;
	sol0[nodenum2+m1]=0;
	}
/**/
/**/
/**/
/* Erosion-sedimentation, Hydration, Melting  account for all markers */
for (mm1=0;mm1<marknum;mm1++)
{
/* Check markers out of grid */
/**/
if(markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && markk[mm1]>0 && markt[mm1]<50)
	{
	/* Up Left Node X,Y Num */
	wn[0]=m1serch((double)markx[mm1]);
	/**/
	/* Erosion/sedimentation account */
	if(erosmod) erosmark(mm1,yn);
	/* Erosion */
	if(erosmod==0 && eroslev>=0 && markt[mm1]>1)
		{
		/* Change marker type */
		if((double)(marky[mm1])<eroslev) {markt[mm1]=0; markd[mm1]=-1.0; markw[mm1]=0;}
		}
	/* Sedimentation */
	if(erosmod==0 && sedilev>=0 && markt[mm1]<2)
		{
		/* Change marker type */
		if((double)(marky[mm1])>sedilev) {markt[mm1]=yn; markd[mm1]=-1.0;}
		}
	/**/
	/* Water/Air account */
	if(waterlev>=0 && markt[mm1]<2)
		{
		/* Change marker type */
		if((double)(marky[mm1])>waterlev) markt[mm1]=1; else markt[mm1]=0;
		}
	/**/
	/* P, T parameters calc */
	allinterp((double)(markx[mm1]),(double)(marky[mm1]));
	mpb=eps[10]*1e-5;
	mtk=(double)(markk[mm1]);
	/**/
	/* Mantle to Antigorite transformation */
	antigor(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1);
/*
*/
	/**/
	/* Rocks to rock+melt transformation */
	if(timesum>0) melting(mtk,mpb,mm1);
/*
*/
	}
}
/**/
/**/
/**/
/* Add ro[] nu[] etc. using selected markers */
for (mm1=0;mm1<marknum;mm1+=gridmod)
{
/* Marker type */
mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
/* Check markers out of grid */
/**/
if(markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && markk[mm1]>0 && markt[mm1]<50)
	{
	/* P, T parameters calc */
	allintere((double)(markx[mm1]),(double)(marky[mm1]));
	mpb=eps[10]*1e-5;
	mtk=(double)(markk[mm1]);
	/**/
	/* Serpentinize Lithospheric Mantle at the outer rise */
	if(mm2==9 && mtk<973.0 && ABSV(marke[mm1])>strainserp && markx[mm1]>1e+6 && markx[mm1]<3e+6)
		{
		/* Erosion level calc */
		e=(markx[mm1]-gx[wn[0]])/(gx[wn[0]+1]-gx[wn[0]]);
		n=(e*ep[wn[0]+1]+(1.0-e)*ep[wn[0]]);
		if((marky[mm1]-n)<=depthserp)
			{
			markt[mm1]=mm2=13;
			markw[mm1]=0;
			markd[mm1]=-1.0;
			markex[mm1]=0;
			}
		}
	/**/
	/* Remove Plumes at the bottom */
	if(marky[mm1]>zdeep && mm2!=9 && mm2!=10) 
		{
		mm2=9; markt[mm1]=10;
		markw[mm1]=0;
		markd[mm1]=-1.0;
		markex[mm1]=0;
		}
	/**/
	/* Create crust and lithosphere at midocean ridges */
	/* Create oceanic plate*/
	if(0==0)
	{
	/* Create Lithosphere */
	if(mm2==10  && marky[mm1]<9e+4 && markk[mm1]<1573.0 && (markx[mm1]<coolzone || markx[mm1]>xsize-coolzone))
		{
		markt[mm1]=mm2=9;
		}
	if(0==0)
	{
	/* Create Gabbroic crust */
	if(mm2==9 || mm2==10 || mm2==29 || mm2==30 && (markx[mm1]<coolzone || markx[mm1]>xsize-coolzone))
		{
		/* Erosion level calc */
		e=(markx[mm1]-gx[wn[0]])/(gx[wn[0]+1]-gx[wn[0]]);
		n=(e*ep[wn[0]+1]+(1.0-e)*ep[wn[0]]);
		if((marky[mm1]-n)<=dysediment+dybasalt+dygabbro)
			{
			markt[mm1]=mm2=8;
			markw[mm1]=0;
			markd[mm1]=-1.0;
			markex[mm1]=0;
			}
		}
	/* Create Basaltic crust */
	if(mm2==8 && (markx[mm1]<coolzone || markx[mm1]>xsize-coolzone))
		{
		/* Erosion level calc */
		e=(markx[mm1]-gx[wn[0]])/(gx[wn[0]+1]-gx[wn[0]]);
		n=(e*ep[wn[0]+1]+(1.0-e)*ep[wn[0]]);
		if((marky[mm1]-n)<=dysediment+dybasalt)
			{
			markt[mm1]=mm2=7;
			markw[mm1]=0;
			markd[mm1]=-1.0;
			markex[mm1]=0;
			}
		}
	/* Create sediments */
	if(mm2==7 && (markx[mm1]<coolzone || markx[mm1]>xsize-coolzone))
		{
		/* Erosion level calc */
		e=(markx[mm1]-gx[wn[0]])/(gx[wn[0]+1]-gx[wn[0]]);
		n=(e*ep[wn[0]+1]+(1.0-e)*ep[wn[0]]);
		if((marky[mm1]-n)<=dysediment)
			{
			markt[mm1]=mm2=3;
			markw[mm1]=0;
			markd[mm1]=-1.0;
			markex[mm1]=0;
			}
		}
	}
	}
	/**/
	/* Save old density */
	mroold=markd[mm1];
	/* Temperature reset for water/air */
/*
	if(mm2<2) mtk=markk[mm1]=273.0; 
*/
	mwa=0;
/*
printf("k %ld %d %e %e ",mm1,mm2,mtk,mpb);getchar();
*/
	/**/
	/* Viscosity calc */
	mnu=markn[mm1];
	if(mnu<=0 || timesum==0) 
		{
		markeii=pow(0.5*(eps[6]*eps[6]+eps[8]*eps[8])+eps[4]*eps[4],0.5);
		marksii=pow(0.5*(eps[7]*eps[7]+eps[9]*eps[9])+eps[5]*eps[5],0.5);
		mnu=viscalc(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1,mm2);
		markn[mm1]=mnu;
		}
	/**/
	/* Make Lithosphere */
	if (mtk<1300.0 && marky[mm1]<200000.0 && mm2==10) {mm2=9; markt[mm1]=9;}
	/* Make Asthenosphere */
/*
	if (mtk>1300.0 && marky[mm1]<200000.0 && mm2==12) {mm2=10; markt[mm1]=10;}
*/
	/**/
	/* Marker Properties */
	mro0=mro=dencalc(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm2);
	/**/
	/* Density && Water wt% save */
	if(timesum<=1e+11 || markd[mm1]<=0) 
		{
		watercalc(mtk,mpb,mm2,mm1);
		mwa=eps[42];
		markw[mm1]=mwa;
		}
	/* Correct density of the rock for the water content */
	wro=1050.0;
	mro=mro*(1.0-markw[mm1]*1e-2)+wro*markw[mm1]*1e-2;
	markd[mm1]=mro;
	/**/
	/* Save properties */
	mbb0=mbb=eps[20];
	mab0=mab=eps[44];
	maa0=maa=eps[45];
	mcp0=mcp=markcp[mm2];
	mkt0=mkt=(markkt[mm2]+markkf[mm2]/(mtk+77.0))*exp(markkp[mm2]*mpb);
	mht=markht[mm2];
	/* Test Heat conductivity k=ko/(1+b*(T-To)/To) */
/*
	if (markkt[mm2]<0) mkt=-markkt[mm2]/(1.0+markkf[mm2]*(mtk-markkp[mm2])/markkp[mm2]);
*/
/*
printf("k %ld %d %e %e %e %e %e %e",mm1,mm2,mtk,markk[mm1],mkt,markkt[mm2],markkf[mm2],markkp[mm2]);getchar();
*/
	/**/
	/* Melted rocks */
	xmelt=0;
	if (mm2>20) 
		{
		meltpart(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1,mm2);
		xmelt=eps[21];
		mro0=mro=eps[23];
		mbb0=mbb=eps[20];
		mab0=mab=eps[44];
		maa0=maa=eps[45];
		mcp0=mcp=eps[25];
		mkt0=mkt=eps[26];
		}
/*
printf("%ld %d %e %e  %e %e",mm1,mm2,mtk,mpb,mbb,mab);getchar();
*/
	/**/
	/* Thermodynamic database use for ro, Cp */
	if (densimod>=2)
	if(mm2>1 && mm2!=5 && mm2!=6 && mm2!=25 && mm2!=26)
		{
		/* Compute TD variables */
		eps[47]=mkt;
		tdbasecalc(mtk,mpb,mm2,mm1);
		mro=eps[41];
		mwa=eps[42];
		mcp=eps[43];
		mbb=eps[20];
		mab=eps[44];
		maa=eps[45];
		mkt=eps[47];
/*
printf("%ld %d %e %e  %e %e",mm1,mm2,mtk,mpb,mbb,mab);getchar();
*/
		/**/
		/* Molten rocks with extraction */
		if(mm2>20)
			{
			/* Compute TD variables */
			tdbasecalc(mtk,mpb,mm2-20,mm1);
			/* Limit melt density */
			if(xmelt>0 && xmelt<=1)
				{
				ival0=markro[mm2]/markro[mm2-20];
				ival=(mro-eps[41]*(1.0-xmelt))/xmelt/eps[41];
				if(ival<ival0) mro=eps[41]*(1.0-xmelt)+eps[41]*ival0*xmelt;
				}
			/* Account for melt extraction */
			if(markex[mm1]>0)
				{
				ival0=xmelt-markex[mm1]; if(ival0<0) ival0=0;
				ival=0; if(xmelt>0) ival=ival0/xmelt;
				mro=mro*ival+eps[41]*(1.0-ival);
				mwa=mwa*ival+eps[42]*(1.0-ival);
				mcp=mcp*ival+eps[43]*(1.0-ival);
				mbb=mbb*ival+eps[20]*(1.0-ival);
				mab=mab*ival+eps[44]*(1.0-ival);
				mkt=mkt*ival+eps[47]*(1.0-ival);
				}
			}
		/**/
		/* Density && Water wt% save */
		if(1==0 || timesum<=1e+11 || markd[mm1]<=0) 
			{
			markw[mm1]=mwa;
			markd[mm1]=mro;
			}
		else
			{
			/* Recompute rock density on the basis of water density */
			wro=1050.0;
			/* Density of the solid */
			ival=(mro-1e-2*mwa*wro)/(1.0-1e-2*mwa);
			/* Density of the rock */
			mro=ival*(1.0-markw[mm1]*1e-2)+wro*markw[mm1]*1e-2;
			markd[mm1]=mro;
			}
/*
{printf("TD1 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1]);getchar();}
{printf("TD2 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1]);getchar();}
*/
		}
	/**/
	/**/
	/**/
	/* End Errosion/Sedimentation account --------------------- */
	/**/
	/**/
	/**/
	/* Water/Air account ====================================*/
	if(1==0 && waterlev>=0 && mm2<2)
		{
		/* Water/Air proportions definition */
		/* Water lewel definition relatively to the bottom of marker */
		markwlev=((double)(marky[mm1])-(waterlev-dywa))/(2.0*dywa);
		if(markwlev<1.0 && markwlev>0)
			{
			/* Air Density */
			aro=dencalc(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),0);
			/* Water Density */
			wro=dencalc(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),1);
			/* Effective  Density recalculation */
			mro0=mro=wro*markwlev+aro*(1.0-markwlev);
			}
		}
	/* End Water/Air account ====================================*/
	/**/
/*
printf("num=%ld type=%d  x=%e y=%e mpb=%e mtk=%e nu=%e ro=%e cp=%e kt=%e ht=%e",mm1,mm2,markx[mm1],marky[mm1],mpb,mtk,mnu,mro,mcp,mkt,mht);getchar();
*/
	/* Restore simple density */
	if(densimod==3) 
		{
		mro=mro0;
		mbb=mbb0;
		mab=mab0;
		maa=maa0;
		mcp=mcp0;
		mkt=mkt0;
		}
	/* Save new density */
	markd[mm1]=mro;
	/* Add divergence */
	if(1==0 && mm2>1 && mroold>0 && timestep>0 && timesum>1e+11)
		{
		markdi[mm1]+=mab*(markk[mm1]-markdt[mm1])/timestep*markvol;
		}
	maa=0;
/*
	if(timesum<1e+11 || markx[mm1]<1e+6) maa=0;
*/
	/**/
	/* Limit viscosity */
	if(marky[mm1]>(ysize*0.75+zdeep*0.25) && mnu>nudeep) mnu=nudeep; 
	/* Interpolation from markers to nodes ====================================*/
	/* Marker weight calculation using dimension of current Cell */
	celdx=gx[wn[0]+1]-gx[wn[0]];
	celdy=gy[wn[1]+1]-gy[wn[1]];
	swt1=(1.0-markwt[mm1])*(1.0-markex[mm1]);
	/**/
	/* Horizontal,Vertical limits for interpolation calc */
	m1min=wn[0]-intermod1; if(m1min<0) m1min=0;
	m1max=wn[0]+1+intermod1; if(m1max>xnumx-1) m1max=xnumx-1;
	/**/
	m2min=wn[1]-intermod1; if(m2min<0) m2min=0;
	m2max=wn[1]+1+intermod1; if(m2max>ynumy-1) m2max=ynumy-1;
	/**/
	/* Interpolation weights calc  after Fornberg (1996) */
	nodewt(m1min,m1max,m2min,m2max,(double)(markx[mm1]),(double)(marky[mm1]),0,0);
	/**/
	/* Interpolate ro,Nu etc to nodes using interpolation koefficients */
	/* Add Normal viscosity nd and divergence dv */
	ddx=1.0-ABSV(cn[wn[0]-m1min][1]-0.5);
	ddy=1.0-ABSV(cn[wn[1]-m2min][0]-0.5);
	swt2=swt1*ddx*ddy;
	m3=(wn[0]+1)*ynumy+wn[1]+1;
/*
printf("%ld %ld %ld %e %e %e",wn[0],wn[1],m3,ddx,ddy,swt2);getchar();
*/
	if(ddx>0.5 && ddy>0.5)
		{
		swt2=swt1*4.0*(ddx-0.5)*(ddy-0.5);
		if(viscmod==0) nd0[m3]+=mnu*swt2;
		if(viscmod==1) nd0[m3]+=log(mnu)*swt2;
		if(viscmod==2) nd0[m3]+=1.0/mnu*swt2;
		pr0[m3]+=markpr[mm1]*swt2;
		prb0[m3]+=maa*swt2;
		sol1[nodenum+m3]+=swt2;
/*
markdi[mm1]=0;
if(timesum<1e+11) markdi[mm1]=0;
*/
		dv[m3]+=markdi[mm1]/celdx/celdy;
/*
		if (markdi[mm1]!=0){printf("%ld %d %e %e %e %e",mm1,markt[mm1],markx[mm1],marky[mm1],markdi[mm1],swt2);getchar();}
*/
		}
	/* Interpolate ro,Nu etc to nodes using interpolation koefficients */
	/* Reload horizontal coefficients to cn[] */
	for (m1=m1min;m1<=m1max;m1++)
	for (m2=m2min;m2<=m2max;m2++)
		{
		/* Current node num, wt */
		m3=m1*ynumy+m2;
		swt=swt1*cn[m1-m1min][1]*cn[m2-m2min][0];
		ddx=cn[m1-m1min][1];
		ddy=cn[m2-m2min][0];
		swt2=swt1*ddx*ddy;
		/**/
/*
if(m1==786 && m2==93) {printf("%ld %ld %ld     %ld %d   %e %e %e\n",m1,m2,m3,mm1,mm2,mtk,mpb,mro); getchar();}
*/
		/* Add Physical Properties: ro,nu, etc. */
		if(ddx>0.5 && ddy>0.5)
			{
			swt2=swt1*4.0*(ddx-0.5)*(ddy-0.5);
			if(viscmod==0) nu0[m3]+=mnu*swt2;
			if(viscmod==1) nu0[m3]+=log(mnu)*swt2;
			if(viscmod==2) nu0[m3]+=1.0/mnu*swt2;
			sol0[nodenum+m3]+=swt2;
			}
		/* Add other Physical Properties: ro, cp etc. */
		ro0[m3]+=mro*swt;
		et0[m3]+=mbb*swt;
		cp0[m3]+=mcp*mro*swt;
		kt0[m3]+=mkt*swt;
		ht0[m3]+=mht*swt;
		sol0[m3]+=swt;
		/* Mark rock marker use */
		if(mm2>1) sol0[nodenum2+m3]=10.0;
		/**/
		/* Add T */
		if(!markim[mm2]) 
			{
			tk0[m3]+=mtk*swt;
			sol1[m3]+=swt;
			}
		}
	/* End Interpolation from markers to nodes ====================================*/
	}
}
/**/
/**/
/* Recalc ro[] nu[] */
for (m1=0;m1<xnumx;m1++)
for (m2=0;m2<ynumy;m2++)
{
/* Current node num, wt */
m3=m1*ynumy+m2;
/* Shear viscosity */
if(sol0[nodenum+m3])
	{
	/* Viscosity recalc check */
	if(mu[m3] && timesum<timebond) 
		{
		if(mu[m3]>0) 
			{
			nu0[m3]=mu[m3];
			}
		else
			{
			nu0[m3]/=sol0[nodenum+m3];
			if(viscmod==1) nu0[m3]=exp(nu0[m3]);
			if(viscmod==2) nu0[m3]=1.0/nu0[m3];
			if(nu0[m3]>-mu[m3]) nu0[m3]=-mu[m3];
			}
		} 
	else
		{
		nu0[m3]/=sol0[nodenum+m3];
		if(viscmod==1) nu0[m3]=exp(nu0[m3]);
		if(viscmod==2) nu0[m3]=1.0/nu0[m3];
		}
	/* Min,Max NU limitation */
/*
printf("%ld %e",m3,nu0[m3]);getchar();
*/
	if(nu0[m3]<nubeg) nu0[m3]=nubeg;
	if(nu0[m3]>nuend) nu0[m3]=nuend;
	/* Min,Max NU definition for nu contrast limit */
	if(numin==0 || nu0[m3]<numin) numin=nu0[m3];
	if(numax==0 || nu0[m3]>numax) numax=nu0[m3];
	nu[m3]=nu0[m3];
	/**/
	/* Reset weight */
	sol0[nodenum+m3]=0;
	}
/* Normal viscosity divergence */
if(sol1[nodenum+m3])
	{
	/* Viscosity recalc check */
	if(mu[m3] && timesum<timebond) 
		{
		if(mu[m3]>0) 
			{
			nd0[m3]=mu[m3];
			}
		else
			{
			nd0[m3]/=sol1[nodenum+m3];
			if(viscmod==1) nu0[m3]=exp(nu0[m3]);
			if(viscmod==2) nu0[m3]=1.0/nu0[m3];
			if(nd0[m3]>-mu[m3]) nd0[m3]=-mu[m3];
			}
		} 
	else
		{
		nd0[m3]/=sol1[nodenum+m3];
		if(viscmod==1) nd0[m3]=exp(nd0[m3]);
		if(viscmod==2) nd0[m3]=1.0/nd0[m3];
		}
	/* Min,Max NU limitation */
/*
printf("%ld %e",m3,nu0[m3]);getchar();
*/
	if(nd0[m3]<nubeg) nd0[m3]=nubeg;
	if(nd0[m3]>nuend) nd0[m3]=nuend;
	/* Min,Max NU definition for nu contrast limit */
	if(numin==0 || nd0[m3]<numin) numin=nd0[m3];
	if(numax==0 || nd0[m3]>numax) numax=nd0[m3];
	nd[m3]=nd0[m3];
	/* Divergence */
	pr0[m3]/=sol1[nodenum+m3];
	prb[m3]=prb0[m3]/sol1[nodenum+m3];
	/* Limit Divergence per cells */
	if(dv[m3]>strmax) dv[m3]=strmax;
	if(dv[m3]<-strmax) dv[m3]=-strmax;
/*
if(dv[m3]!=0){printf("%ld %ld %e",m1,m2,dv[m3]);getchar();}
	dv[m3]/=sol1[nodenum+m3];
*/
	/**/
	/* Reset weight */
	sol1[nodenum+m3]=0;
	}
else
	{
	pr0[m3]=pr[m3];
	}
/* Other material properties */
if(sol0[m3])
	{
	/* Material constants recalc */
	ro[m3]=ro0[m3]/sol0[m3];
	/* Flatten density distribution for water/air boundary */
/*
	if(sol0[nodenum2+m3]==0)
		{
		if(gy[m2]<waterlev && ro[m3]<1000.1) ro[m3]=1.0;
		if(gy[m2]>=waterlev && ro[m3]<1000.1) ro[m3]=1000.0;
		}
*/
	et[m3]=et0[m3]/sol0[m3];
	cp[m3]=cp0[m3]/sol0[m3]/ro[m3];
	kt[m3]=kt0[m3]/sol0[m3];
	ht[m3]=ht0[m3]/sol0[m3];
	/**/
	/* Advective addition for T K in nodes recalc */
	if (sol1[m3]) 
		{
		tk[m3]=tk0[m3]/sol1[m3]; 
/*
	if (sol1[m3] && tk3[m3]) 
		tk[m3]+=tk0[m3]/sol1[m3]-tk3[m3]; 
printf("%ld %e",m3,tk[m3]);getchar();
*/
		sol1[m3]=0;
		}
	/**/
	/* Reset weight */
	sol0[m3]=0;
	}
}
printf("Min, Max viscosity %e %e \n",numin,numax);
/*
printf("%e %e",numin,numax);getchar();
printf("%e %e %e   %e %e ",nubeg,nuend,nucontr,numin,numax);getchar();
*/
/**/
/* Reset advective temperature */
for (m3=0;m3<nodenum;m3++) tk3[m3]=0; 
/**/
/* Set Upper/Lower limits for nu[] after given contrast */
if(nucontr>=1.0) numax=numin*nucontr; else numin=numax*nucontr;
for (m3=0;m3<nodenum;m3++)
	{
	if(nu[m3]<numin) nu[m3]=numin;
	if(nu[m3]>numax) nu[m3]=numax;
	if(nd[m3]<numin) nd[m3]=numin;
	if(nd[m3]>numax) nd[m3]=numax;
	}
/**/
/* Set Boundary iconditions for T */
if (printmod) printf("\n AVERAGE TEMPERATURE CORRECTION FOR BOUNDARY CONDITIONS ...\n");
tkrecalc();
if (printmod) printf("AVERAGE TEMPERATURE OK!\n");
/**/
/**/
/**/
/*
printf("%e %e %e   %e %e ",nubeg,nuend,nucontr,numin,numax);getchar();
*/
/* ADD MARKERS TO THE v-CELLS ========================== */
}
/* End ro[],nu[] recalc after marker positions */



/* Calc ro for given P,T */
double dencalc(double mtk, double mpb, double x, double y, int mm2)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm2 - Rock number */
{
/* Val buffer */
double ival;
double p_pl_out,p_ga_in,rokf,p_sp_in,p_ol_out,p_pv_in,p_sp_out,p_st_in;
eps[20]=eps[44]=eps[45]=0;
/**/
/* Constant density */
if (densimod==0) return markro[mm2];
/**/
/* Ro=ro0*(1-bro*(TK-298.15))*(1+aro*(Pkbar-0.001)) */
ival=markro[mm2]*(1.0-markbb[mm2]*(mtk-298.15))*(1.0+markaa[mm2]*(mpb-1.0)*1e-3);
/* Adiabatic term: al=bro/(1-bro*(Tk-298.15)) */
eps[44]=eps[20]=markbb[mm2]/(1.0-markbb[mm2]*(mtk-298.15));
/* Compressibility */
eps[45]=markaa[mm2]*1e-8/(1.0+markaa[mm2]*(mpb-1.0)*1e-3);
/**/
/* Eclogitization, St, Pv transitions in oceanic crust */
if(mm2==7 || mm2==8 )
        {
	/* Set initial density coefficient */
	rokf=1.0;
	/* Eclogitization Ito and Kennedy, 1971 */
        /*basalt=>garnet granulite (Ga-In) transition*/
        p_ga_in=-9222.0+mtk*14.0;
        /*Not to have granulites at pressure lower than 2 kbar*/
        if(p_ga_in<2000.0) p_ga_in=2000.0;
        /*garnet granulite=>eclogite (Pl-Out) transition*/
        p_pl_out=-1460.0+mtk*20.0;
        /*Not to have eclogites at pressure lower than 12 kbar*/
        if(p_pl_out<12000.0) p_pl_out=12000.0;
        if(mpb>p_ga_in)
                {
                if(mpb>=p_pl_out)
                   	{
                       	rokf*=1.16;
                        }
                else
                        {
                        rokf*=1.0+0.16*(mpb-p_ga_in)/(p_pl_out-p_ga_in);
                        }
                }
	/* Coe->St transition Gerya et al., 2004, PCM */
        p_st_in=59100.0+mtk*22.6;
        if(mpb>p_st_in) rokf*=1.06;
	/* Pv transition, Mishin et al., 2008 with slope from Ito et al., 1990 */
        /* Sp-out transition*/
        p_sp_out=354000.0-mtk*40.0;
        /* Pv-in transition*/
        p_pv_in=352000.0-mtk*40.0;
        if(mpb>p_pv_in)
                {
                if(mpb>=p_sp_out)
                   	{
                       	rokf*=1.08;
                        }
                else
                        {
                        rokf*=1.0+0.08*(mpb-p_pv_in)/(p_sp_out-p_pv_in);
                        }
                }
	/* Take into account kynetics */
	if(mtk<teclmax)
		{
		if(mtk<teclmin)
       			{
       			rokf=1.00;
                      	}
		else
                      	{
                       	rokf=1.0+(rokf-1.0)*(mtk-teclmin)/(teclmax-teclmin);
                       	}
                }
	/* Correct density */
        ival*=rokf;
        }
/**/
/* Ol-Sp and Pv transitions in the mantle */
if(mm2>=9 && mm2<=14) 
        {
	/* Set initial density coefficient */
	rokf=1.0;
	/* Ol-Sp transition, Katsura & Ito, 1989 */
        /* Ol-out transition*/
        p_ol_out=91000.0+mtk*27.0;
        /* Sp-in transition*/
        p_sp_in=66000.0+mtk*39.0;
        /*Limit width of Sp-Ol transition to 2 kbar */
        if(p_sp_in>p_ol_out-2000.0) p_sp_in=p_ol_out-2000.0;
        if(mpb>p_sp_in)
                {
                if(mpb>=p_ol_out)
                   	{
                       	rokf*=1.06;
                        }
                else
                        {
                        rokf*=1.0+0.06*(mpb-p_sp_in)/(p_ol_out-p_sp_in);
                        }
                }
	/* Pv transition, Ito et al., 1990 */
        /* Sp-out transition*/
        p_sp_out=304000.0-mtk*40.0;
        /* Pv-in transition*/
        p_pv_in=302000.0-mtk*40.0;
        if(mpb>p_pv_in)
                {
                if(mpb>=p_sp_out)
                   	{
                       	rokf*=1.11;
                        }
                else
                        {
                        rokf*=1.0+0.11*(mpb-p_pv_in)/(p_sp_out-p_pv_in);
                        }
                }
	/* Take into account kynetics */
	if(mtk<teclmax)
		{
		if(mtk<teclmin)
       			{
       			rokf=1.00;
                      	}
		else
                      	{
                       	rokf=1.0+(rokf-1.0)*(mtk-teclmin)/(teclmax-teclmin);
                       	}
                }
	/* Correct density */
        ival*=rokf;
        }
/**/
return  ival;
}
/* End Calc ro for given P,T */




/* Antigorite weakening of mantle */
void antigor(double mtk, double mpb, double x, double y, long int mm1)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm1 - mark number */
{
/* Val buffer */
double k1,sy1,e,hydry,hydryl;
int mm2=(int)markt[mm1];
long int m1;
/**/
/* Check marker type */
if(mm2<9 || (mm2>13 && mm2!=34)) return;
/**/
/* Up Left Node X,Y Num */
m1=wn[0];
m1=m1serch((double)markx[mm1]);
/* Relativ Normalized coord Calc */
e=(x-gx[m1])/(gx[m1+1]-gx[m1]);
/*Fluid disappearance surface */
sy1=sedilev;
/**/
/**/
/* Lower hydration boundary elevation for x marker position */
hydryl=e*ep[m1+xnumx*3+1]+(1.0-e)*ep[m1+xnumx*3];
/* Upper hydration boundary for  x marker position */
hydry=e*ep[m1+xnumx*2+1]+(1.0-e)*ep[m1+xnumx*2];
/*
printf("%d %e %e %e %e %e",m1,x,y,e,hydry,hydryl);getchar();
printf("%ld %e %e %e %e %e %e",mm1,x,y,xsubd,ysubd,vxs,vys,);getchar();
*/
/* Antigorite weakening of mantle above oceanic crust */
/*
printf("%d %e %e %e %e %e",m1,x,y,xsubd,e,hydry);getchar();
printf("%ld %e %e %e %e %e %e ",mm1,x,y,xsubd,ysubd,vxs,vys);getchar();
*/
/* Atg stability field after Schmidt and Poli, 1998 */
if((y-sy1)>63000.0)
	{
	k1=1013.17699-0.060387633e-3*(y-sy1)-0.004289442e-6*(y-sy1)*(y-sy1);
	}
else
	{
	k1=751.490422+6.00773668e-3*(y-sy1)-0.034690759e-6*(y-sy1)*(y-sy1);
	}
/* Change marker Type */
/* Unhydrated (9,10) - to hydrated (11) */
if(timesum<=0 && y<hydryl && y>hydry && mm2!=13 && mm2!=11)
	{
	markt[mm1]=11;
	}
/* Serpentinized (13) - to hydrated (11) */
if(k1<=mtk && markt[mm1]==13) markt[mm1]=11;
/* Hydrated(11) -  to serpentinized (13) */
if(k1>mtk && markt[mm1]==11) markt[mm1]=13;
}
/* Antigorite weakening of mantle */



/* Basalt melting in oceanic crust */
double basalt()
{
/* Val buffer */
double x0,x1,y0,y1,k0,k1,t0,t1,sy1;
long int m1;
/**/
/**/
/**/
/* Serch for melting depth along lower hydration surface */
for (m1=1;m1<xnumx;m1++)
	{
	/* TK calc along the lower hydr surface */
	x0=gx[m1-1];
	y0=ep[m1+xnumx*3-1];
	allintert(x0,y0); t0=eps[2];
	x1=gx[m1];
	y1=ep[m1+xnumx*3];
	allintert(x1,y1); t1=eps[2];
	/*Fluid disappearance surface */
	sy1=sedilev;
	if((y0-sy1)>0 && (y1-sy1)>0)
		{
		/* Basalt solidus temperature calc  */
		if((y0-sy1)>48000.0)
			{
			k0=940.19-0.023633e-3*(y0-sy1)+0.007774e-6*(y0-sy1)*(y0-sy1);
			}
		else
			{
			k0=941.60+704.4157e+3/(y0-sy1)+1217.45485e+6/(y0-sy1)/(y0-sy1);
			}
		if((y1-sy1)>48000.0)
			{
			k1=940.19-0.023633e-3*(y1-sy1)+0.007774e-6*(y1-sy1)*(y1-sy1);
			}
		else
			{
			k1=941.60+704.4157e+3/(y1-sy1)+1217.45485e+6/(y1-sy1)/(y1-sy1);
			}
		/* Basalt solidus depth calc */
		if(t1>=k1)
			{
			basalty=y1-(y1-y0)/((t1-k1)+(k0-t0))*(t1-k1);
			return 0;
			}

		}
	}
basalty=ysize;
return 0;
}
/* Basalt melting in oceanic crust */



/* Number of nearest left vertical line find */
long int m1serch(double x)
/* x - X coordinate */
{
/* Variables */
long int m1,m10=0,m11=xnumx-1;
/**/
/* Serch cycle */
do
	{
	m1=(m10+m11)/2;
	if (gx[m1]>x) m11=m1; else m10=m1;
	}
while((m11-m10)>1);
if(m10>xnumx-2) m10=xnumx-2;
/*
if(x<gx[m10] || x>gx[m10+1]) {printf("XXX %ld %ld %ld  %e %e  %e ",m10,m11,m1,gx[m10],gx[m11],x); getchar();}
*/
return m10;
}
/* Number of nearest left vertical line find */




/* Number of nearest upper horizontal line find */
long int m2serch(double y)
/* y - Y coordinate */
{
/* Variables */
long int m2,m20=0,m21=ynumy-1;
/**/
/* Serch cycle */
do
	{
	m2=(m20+m21)/2;
	if (gy[m2]>y) m21=m2; else m20=m2;
	}
while((m21-m20)>1);
if(m20>ynumy-2) m20=ynumy-2;
/*
if(y<gy[m20] || y>gy[m20+1]) {printf("YYY %ld %ld %ld  %e %e  %e ",m20,m21,m2,gy[m20],gy[m21],y); getchar();}
*/
return m20;
}
/* Number of nearest upper horizontal line find */



/* Erosion/Sedimentation Function for markers */
/* mardy - marker vertical size, m */
void erosmark(long int mm1, int yn)
/* mm1 - marker number */
/* yn - current sedimnts type 2,3 */
{
/* Variables */
double e,e0;
long int m1;
/**/
/**/
/**/
/* Surface level elevation definition */
/* Up Left Node X,Y Num */
m1=wn[0];
m1=m1serch((double)markx[mm1]);
/* Relativ Normalized coord Calc */
e=((double)(markx[mm1])-gx[m1])/(gx[m1+1]-gx[m1]);
/**/
/* Surface level elevation for marker definition */
e0=(e*ep[m1+1]+(1.0-e)*ep[m1]);
/**/
/*
printf("MARKER %ld %ld %e %e %e %d %d",mm1,m1,markx[mm1],marky[mm1],e0,markt[mm1],yn);getchar();
*/
/* Marker surface elevation definition */
/* Water/Air -> Sediments conversion */
if(markt[mm1]<2)
	{
	markw[mm1]=0;
	marke[mm1]=0;
	markex[mm1]=0;
	if(marky[mm1]>e0) {markt[mm1]=yn; markd[mm1]=-1.0;markw[mm1]=0;marke[mm1]=0;markex[mm1]=0;}
	}
/**/
/* Rock->Water/Air conversion */
if(markt[mm1]>1 && markt[mm1]<50)
	{
	/* No underwater erosion */
	if((double)(marky[mm1])<e0) {markt[mm1]=0;markd[mm1]=-1.0;markw[mm1]=0;marke[mm1]=0;markex[mm1]=0;}
	}
/*
	if((double)(marky[mm1])<e0 && (markt[mm1]!=7 || marky[mm1]<waterlev)) {markt[mm1]=0;markd[mm1]=-1.0;markw[mm1]=0;marke[mm1]=0;markex[mm1]=0;}
if(mm2!=markt[mm1]) {printf("MARKER %ld %ld %e %e %e %d %d %d",mm1,m1,markx[mm1],marky[mm1],e0,mm2,markt[mm1],yn);getchar();}
*/
}
/* Erosion/Sedimentation Function for markers */



/* Rock to rock+melt transformation */
void melting(double mtk, double mpb, long int mm1)
/* mtk - T, K */
/* mpb - P, bar */
/* mm1 - mark number */
{
/* Marker type */
int mm2=(int)markt[mm1];
/**/
/* Melting in 1000-3000 km interval only */
/*
if (markx[mm1]<1.00e+6   || markx[mm1]>3.00e+6) return;
*/
/**/
/* Melting related change of the marker type */
/* Check marker type */
if (mm2>1)
if (mpb>0)
switch(mm2)
	{
	/* Sediments, upper crust */
	case 15:
	case 35:
	case 3:
	case 4:
	case 5:
	case 23:
	case 24:
	case 25:
	/* Basalt, Gabbro */
	case 6:
	case 7:
	case 8:
	case 16:
	case 17:
	case 18:
	case 36:
	case 37:
	case 38:
	case 27:
	case 28:
	case 26:
	meltpart1(mtk,mpb,mm1);
	if(eps[21]>0 && mm2<20) markt[mm1]+=20;
	if(eps[21]<=0 && mm2>20) markt[mm1]-=20;
/*
if (mm2!=markt[mm1]){printf("Granite/Basalt %ld %d %d  %e %e  %e",mm1,mm2,markt[mm1],mtk,mpb,eps[21]);getchar();}
*/
 	return;
	/**/
	/* Hydrated Peridotite */
	case 11:
	case 13:
	case 14:
	case 34:
	meltpart1(mtk,mpb,mm1);
	if(eps[21]>0 && mm2<20) markt[mm1]=34;
	if(eps[21]<=0 && mm2==34) markt[mm1]=14;
/*
if (mm2!=markt[mm1]){printf("Peridotite %ld %d %d  %e %e  %e",mm1,mm2,markt[mm1],mtk,mpb,eps[21]);getchar();}
*/
 	return;
	/* Dry Peridotite */
	case  9:
	case 10:
	case 12:
	case 29:
	case 30:
	case 32:
	meltpart1(mtk,mpb,mm1);
	if(eps[21]>0 && mm2<20) markt[mm1]+=20;
	if(eps[21]<=0 && mm2>20) markt[mm1]-=20;
 	return;
	/* Others */
	default: return;
	}
}
/* Rock to rock+melt transformation */




/* Melt fraction, density, viscosity, heat capacity calculation */
double meltpart(double mtk, double mpb, double x, double y, long int mm1, int mm2)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm1 - mark number */
/* mm2 - mark type */
{
/* Val buffer */
double xmelt=0,xmelt0=0,ival,ival1,dmpb,dmtk,epsin,sduct,nueff,smin,smax,nmin,nmax,cpadd=0;
long int m1;
/**/
/* Check marker type */
if (mm2>20)
	{
	/* Calculate melt fraction */
	meltpart1(mtk,mpb,mm1);
	xmelt0=eps[21];
	xmelt=eps[21]-markex[mm1];if(xmelt<0) xmelt=0;
	/**/
	/* Standard adiabatic term: al=bro/(1+bro*(Tk-298.15)) */
	eps[44]=eps[20]=(markbb[mm2]*xmelt+markbb[mm2-20]*(1.0-xmelt))/(1.0-(markbb[mm2]*xmelt+markbb[mm2-20]*(1.0-xmelt))*(mtk-298.15));
	eps[45]=1e-8*(markaa[mm2]*xmelt+markaa[mm2-20]*(1.0-xmelt))/(1.0+(markaa[mm2]*xmelt+markaa[mm2-20]*(1.0-xmelt))*(mpb-1.0)*1e-3);
	/**/
	/* Density */
	/* Ro=ro0 */
	if (densimod==0) 
		{
		eps[23]=markro[mm2]*xmelt+markro[mm2-20]*(1.0-xmelt);
		}
	/* Ro=ro0*(1-bro*(TK-298.15))*(1+aro*(Pkbar-0.001)) */
	else
		{
		eps[23]=(1.0-xmelt)*dencalc(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm2-20)+xmelt*dencalc(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm2);
		}
	/**/
	/* Heat capacity */
	eps[25]=markcp[mm2]*xmelt+markcp[mm2-20]*(1.0-xmelt);
	/**/
	/* heat conductivity */
	eps[26]=((markkt[mm2]+markkf[mm2]/(mtk+77.0))*exp(markkp[mm2]*mpb))*xmelt+((markkt[mm2-20]+markkf[mm2-20]/(mtk+77.0))*exp(markkp[mm2-20]*mpb))*(1.0-xmelt);
	/**/
	/* Additional melting adiabatic term, heat capacity */
	if(xmelt>0 && xmelt<1.0)
		{
		/* Melting adiabatic term: alm=-ro*(dHlat/dP)/T */
		/* Numerical differentiation */
		dmpb=300.0;
		meltpart1(mtk,mpb-dmpb,mm1);
		ival=eps[22];
		meltpart1(mtk,mpb+dmpb,mm1);
		ival-=eps[22];
		ival*=eps[23]/(mtk*2.0*dmpb*1e+5);
		eps[20]+=ival;
		if(eps[20]<-1e-2) eps[20]=-1e-2; else if(eps[20]>1e-2) eps[20]=1e-2;
		/**/
		/* Melting heat capacity term: cpm=dHlat/dT */
		/* Melting thermal expansion term: aab=-1.0/ro*dro/dT=-1.0/ro*dXmelt/dT*dro/dXmelt */
		/* Numerical differentiation */
		dmtk=10.0;
		meltpart1(mtk+dmtk,mpb,mm1);
		ival=eps[22];
		ival1=eps[21];
		meltpart1(mtk-dmtk,mpb,mm1);
		ival-=eps[22];
		ival1-=eps[21];
		ival/=2.0*dmtk;
		ival1/=2.0*dmtk;
		eps[25]+=ival;
		if(eps[25]<1e+2) eps[25]=1e+2; else if(eps[25]>5e+4) eps[25]=5e+4;
		eps[44]-=ival1/eps[23]*(markro[mm2]*(1.0-markbb[mm2]*(mtk-298.15))*(1.0+markaa[mm2]*(mpb-1.0)*1e-3)-markro[mm2-20]*(1.0-markbb[mm2-20]*(mtk-298.15))*(1.0+markaa[mm2-20]*(mpb-1.0)*1e-3));
/*
printf("Meltpart: %ld %d %e %e %e  %e %e %e %e %e %e",mm1,mm2,mtk,mpb,xmelt,eps[20],eps[22],eps[23],eps[24],eps[25],eps[26]);getchar();
*/
		}
	eps[21]=xmelt0;
	return 1.0;
	}
eps[20]=eps[21]=eps[22]=eps[23]=eps[24]=eps[25]=eps[26]=eps[44]=0;
return 0;
}
/* Rock to rock+melt transformation */





/* Calculation of Vx,Vy,EPSxx,EPSyy,EPSxy,SIGxx,SIGyy,SIGxy,T,T0,T1,T2,P for current location by Linear Interpolation */
/* Staggered Nodes num */
/*   [0]                [3]                [6]   */
/*  T0,xy0    Vy0     T3,xy3     Vy3             */
/*                                               */
/*   Vx0    P4,xx4,yy4  Vx3    P7,xx7,yy7        */
/*                                               */
/*   [1]                [4]                [7]   */
/*  T,xy1     Vy1     T4,xy4     Vy4             */
/*                                               */
/*   Vx1    P5,xx5,yy5  Vx4    P8,xx8,yy8        */
/*                                               */
/*   [2]                [5]                [8]   */
/*                                               */
/*                                               */
void allinter(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* Buffer clear */
for (m1=0;m1<=15;m1++) eps[m1]=0;
/**/
/**/
/**/
/* T interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* T Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[0]+=ival*tk0[m3];
	eps[1]+=ival*tk1[m3];
	eps[2]+=ival*tk[m3];
	eps[3]+=ival*tk2[m3];
	}
/**/
/* Wt for nodes save for T */
wn[2]=m1min; wn[3]=m1max;
for (m1=m1min;m1<=m1max;m1++)
	{
	cn[m1-m1min][5]=cn[m1-m1min][1];
	}
wn[4]=m2min; wn[5]=m2max;
for (m2=m2min;m2<=m2max;m2++)
	{
	cn[m2-m2min][4]=cn[m2-m2min][0];
	}
/* End T interpolation ------------------------ */
/**/
/**/
/**/
/* SIGxy,EPSxy, SIGxy*EPSxy interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* SIGxy,EPSxy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[4]+=ival*exy[m3];
	eps[5]+=ival*sxy[m3];
	eps[13]+=ival*exy[m3]*sxy[m3];
	}
/* End SIGxy,EPSxy interpolation ------------------------ */
/**/
/**/
/**/
/* SIGxx,EPSxx,SIGyy,EPSyy,P, SIGxx*EPSxx, SIGyy*EPSyy interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* SIGxx,EPSxx,SIGyy,EPSyy,P Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[6]+=ival*exx[m3];
	eps[7]+=ival*sxx[m3];
	eps[8]+=ival*eyy[m3];
	eps[9]+=ival*syy[m3];
	eps[10]+=ival*pr[m3];
	eps[14]+=ival*exx[m3]*sxx[m3];
	eps[15]+=ival*eyy[m3]*syy[m3];
	}
/* Pressure as function of depth from erosion level calc */
if (presmod) depthp(x,y);
/* End SIGxx,EPSxx,SIGyy,EPSyy,P interpolation ------------------------ */
/**/
/**/
/**/
/* Vx interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<0) m2min=0;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*vx[m3];
	}
/* End Vx interpolation ------------------------ */
/**/
/**/
/**/
/* Vy interpolation ------------------------ */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x<(gx[m10]+gx[m10+1])/2.0) m1min-=1;
if(m1min<0) m1min=0; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<0) m1min=0;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,+1,0);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[12]+=ival*vy[m3];
	}
/* End Vy interpolation ------------------------ */
/**/
/**/
/**/
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
}
/* Calculation of Vx,Vy,EPSxx,EPSyy,EPSxy,SIGxx,SIGyy,SIGxy,T,T0,T1,T2,P for current location by Linear Interpolation */





/* Weights for horisontal and vertical nodes calculation for marker interpolation */ 
void nodewt(long int m1min, long int m1max, long int m2min, long int m2max, double x, double y, int ynx, int yny)
/* m1min,m1max, m2min,m2max - node X,Y number limits */
/* x,y - curent pont coordinates */
/* ynx, yny - Type of shifts: No(0), Back(-1), Forw(1) */
{
/* Eyy vertical position */
long int m3;
int nx,ny;
/**/
/**/
/**/
/* Weigths in horizontal directions */
/* Load distances to xn[] */
if(ynx<0) 
	{
	for (m3=m1min;m3<=m1max;m3++)
		{
		xn[m3-m1min]=(gx[m3]+gx[m3-1])/2.0;
		}
	}
if(ynx==0) 
	{
	for (m3=m1min;m3<=m1max;m3++)
		{
		xn[m3-m1min]=gx[m3];
		}
	}
if(ynx>0) 
	{
	for (m3=m1min;m3<=m1max;m3++)
		{
		xn[m3-m1min]=(gx[m3]+gx[m3+1])/2.0;
		}
	}
/**/
/* Calc maximal position in xn[] */
nx=(int)(m1max-m1min);
/**/
/* Calc coefficients for horizontal direction */
fdweight(nx,0,x);
/**/
/* Reload horizontal coefficients to cn[] */
for (m3=0;m3<=nx;m3++)
	{
	cn[m3][1]=cn[m3][0];
	}
/**/
/**/
/**/
/* Weigths in vertical directions */
/* Load distances to xn[] */
if(yny<0) 
	{
	for (m3=m2min;m3<=m2max;m3++)
		{
		xn[m3-m2min]=(gy[m3]+gy[m3-1])/2.0;
		}
	}
if(yny==0) 
	{
	for (m3=m2min;m3<=m2max;m3++)
		{
		xn[m3-m2min]=gy[m3];
		}
	}
if(yny>0) 
	{
	for (m3=m2min;m3<=m2max;m3++)
		{
		xn[m3-m2min]=(gy[m3]+gy[m3+1])/2.0;
		}
	}
/**/
/* Calc maximal position in xn[] */
ny=(int)(m2max-m2min);
/**/
/* Calc coefficients for horizontal direction */
fdweight(ny,0,y);
/**/
/**/
/**/
}
/* Weights for horisontal and vertical nodes calculation for marker interpolation */ 



/* Calculation of Vx,Vy by Interpolation from the centers of cells */
void allinterv1(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* Vx, Vy interpolation ------------------------ */
/* Buffer clear */
eps[11]=eps[12]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* Vx, Vy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	/* 0    2 */
	/* 1    3 */
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*(vx[m3-1]+vx[m3-ynumy-1])*0.5;
	eps[12]+=ival*(vy[m3-ynumy]+vy[m3-ynumy-1])*0.5;
	}
/* End EPSxx,EPSyy,P interpolation ------------------------ */
/* End Vx, Vy Interpolate after interpolation weights */
}
/* Calculation of Vx,Vy by Interpolation from the centers of cells */


/* Calculation of Vx,Vy by Interpolation */
void allinterv(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival,xrat;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Check weighting for interpolation */
xrat=2.0/3.0;
if(x<(gx[0]+gx[1])/2.0) xrat=1.0;
if(x>(gx[xnumx-2]+gx[xnumx-1])/2.0) xrat=1.0;
if(y<(gy[0]+gy[1])/2.0) xrat=1.0;
if(y>(gy[ynumy-2]+gy[ynumy-1])/2.0) xrat=1.0;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* Vx interpolation ------------------------ */
/* Buffer clear */
eps[11]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<0) m2min=0;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*vx[m3]*xrat;
	}
/* End Vx interpolation ------------------------ */
/**/
/**/
/**/
/* Vy interpolation ------------------------ */
/* Buffer clear */
eps[12]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x<(gx[m10]+gx[m10+1])/2.0) m1min-=1;
if(m1min<0) m1min=0; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<0) m1min=0;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,+1,0);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[12]+=ival*vy[m3]*xrat;
	}
/* End Vy interpolation ------------------------ */
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* Vx, Vy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	/* 0    2 */
	/* 1    3 */
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*(vx[m3-1]+vx[m3-ynumy-1])*0.5*(1.0-xrat);
	eps[12]+=ival*(vy[m3-ynumy]+vy[m3-ynumy-1])*0.5*(1.0-xrat);
	}
}
/* Calculation of Vx,Vy by Interpolation */




/* Calculation of Vx,Vy by Interpolation */
void allinterv0(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* Vx interpolation ------------------------ */
/* Buffer clear */
eps[11]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<0) m2min=0;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*vx[m3];
	}
/* End Vx interpolation ------------------------ */
/**/
/**/
/**/
/* Vy interpolation ------------------------ */
/* Buffer clear */
eps[12]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x<(gx[m10]+gx[m10+1])/2.0) m1min-=1;
if(m1min<0) m1min=0; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<0) m1min=0;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,+1,0);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[12]+=ival*vy[m3];
	}
/* End Vy interpolation ------------------------ */
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
}
/* Calculation of Vx,Vy by Interpolation */




/* Calculation of EPSxx,EPSyy,EPSxy, P by Interpolation */
void allintere(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* EPSxy interpolation ------------------------ */
/* Buffer clear */
eps[4]=0;
eps[5]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* EPSxy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[4]+=ival*exy[m3];
	eps[5]+=ival*sxy[m3];
	}
/* End SIGxy,EPSxy interpolation ------------------------ */
/**/
/**/
/**/
/* EPSxx,EPSyy, P interpolation ------------------------ */
/* Buffer clear */
eps[6]=eps[8]=eps[10]=0;
eps[7]=eps[9]=0;
/* Horizontal,Vertical limits for interpolation calc */
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* EPSxx,EPSyy,P Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[6]+=ival*exx[m3];
	eps[7]+=ival*sxx[m3];
	eps[8]+=ival*eyy[m3];
	eps[9]+=ival*syy[m3];
	eps[10]+=ival*pr[m3];
	}
/* End EPSxx,EPSyy,P interpolation ------------------------ */
/* Pressure as function of depth from erosion level calc */
if (presmod==1) depthp(x,y);
if (presmod==2) eps[10]=pr[(m10+1)*ynumy+(m20+1)];
}
/* Calculation of EPSxx,EPSyy,EPSxy, P by Interpolation */





/* Calculation of Vx,Vy, EPSxx*SIGxx,EPSyy*SIGyy,EPSxy*SIGxy by Interpolation */
void allinters(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* SIGxy*EPSxy interpolation ------------------------ */
/* Buffer clear */
eps[13]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(m1min<1) m1min=1; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(m2min<1) m2min=1; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* SIGxy,EPSxy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[13]+=ival*exy[m3]*sxy[m3];
	}
/* End SIGxy*EPSxy interpolation ------------------------ */
/**/
/**/
/**/
/* SIGxx*EPSxx, SIGyy*EPSyy interpolation ------------------------ */
/* Buffer clear */
eps[14]=eps[15]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* SIGxx,EPSxx,SIGyy,EPSyy,P Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[14]+=ival*exx[m3]*sxx[m3];
	eps[15]+=ival*eyy[m3]*syy[m3];
	}
/* End SIGxx*EPSxx,SIGyy*EPSyy interpolation ------------------------ */
/**/
/**/
/**/
/* Vx interpolation ------------------------ */
/* Buffer clear */
eps[11]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20; if(y<(gy[m20]+gy[m20+1])/2.0) m2min-=1;
if(m2min<0) m2min=0; if(m2min>ynumy-3) m2min=ynumy-3;
m2max=m2min+1+intermod; if(m2max>ynumy-2) m2max=ynumy-2;
m2min=m2min-intermod; if(m2min<0) m2min=0;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,+1);
/**/
/* Vx Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[11]+=ival*vx[m3];
	}
/* End Vx interpolation ------------------------ */
/**/
/**/
/**/
/* Vy interpolation ------------------------ */
/* Buffer clear */
eps[12]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x<(gx[m10]+gx[m10+1])/2.0) m1min-=1;
if(m1min<0) m1min=0; if(m1min>xnumx-3) m1min=xnumx-3;
m1max=m1min+1+intermod; if(m1max>xnumx-2) m1max=xnumx-2;
m1min=m1min-intermod; if(m1min<0) m1min=0;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,+1,0);
/**/
/* Vy Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[12]+=ival*vy[m3];
	}
/* End Vy interpolation ------------------------ */
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
}
/* Calculation of Vx,Vy, EPSxx*SIGxx,EPSyy*SIGyy,EPSxy*SIGxy by Interpolation */



/* Calculation of T,T0 for current location by Interpolation */
void allintert(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* T interpolation ------------------------ */
/* Buffer clear */
eps[2]=eps[3]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10-intermod; if(m1min<0) m1min=0;
m1max=m10+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
/**/
m2min=m20-intermod; if(m2min<0) m2min=0;
m2max=m20+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,0,0);
/**/
/* T Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[2]+=ival*tk[m3];
	eps[3]+=ival*tk2[m3];
	}
/**/
/* Wt for nodes save for T */
wn[2]=m1min; wn[3]=m1max;
for (m1=m1min;m1<=m1max;m1++)
	{
	cn[m1-m1min][5]=cn[m1-m1min][1];
	}
wn[4]=m2min; wn[5]=m2max;
for (m2=m2min;m2<=m2max;m2++)
	{
	cn[m2-m2min][4]=cn[m2-m2min][0];
	}
/* End T interpolation ------------------------ */
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
}
/* Calculation of T,T0 for current location by Interpolation */



/* Calculation of P by Interpolation */
void allinterp(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m1,m2,m3,m10,m20,m1min,m1max,m2min,m2max;
/* en-NormalisedDistance */
double ival;
/**/
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Up Left Node X,Y Num */
wn[0]=m10=m1serch(x);
wn[1]=m20=m2serch(y);
/**/
/**/
/**/
/* P interpolation ------------------------ */
/* Buffer clear */
eps[10]=0;
/* Horizontal,Vertical limits for interpolation calc */
m1min=m10; if(x>(gx[m10]+gx[m10+1])/2.0) m1min+=1;
if(m1min<1) m1min=1; if(m1min>xnumx-2) m1min=xnumx-2;
m1max=m1min+1+intermod; if(m1max>xnumx-1) m1max=xnumx-1;
m1min=m1min-intermod; if(m1min<1) m1min=1;
/**/
m2min=m20; if(y>(gy[m20]+gy[m20+1])/2.0) m2min+=1;
if(m2min<1) m2min=1; if(m2min>ynumy-2) m2min=ynumy-2;
m2max=m2min+1+intermod; if(m2max>ynumy-1) m2max=ynumy-1;
m2min=m2min-intermod; if(m2min<1) m2min=1;
/**/
/* Interpolation weights calc  after Fornberg (1996) */
nodewt(m1min,m1max,m2min,m2max,x,y,-1,-1);
/**/
/* EPSxx,EPSyy,P Interpolate after interpolation weights */
for (m1=m1min;m1<=m1max;m1++)
for (m2=m2min;m2<=m2max;m2++)
	{
	/* Current node num, wt */
	m3=m1*ynumy+m2;
	ival=cn[m1-m1min][1]*cn[m2-m2min][0];
	eps[10]+=ival*pr[m3];
	}
/* End EPSxx,EPSyy,P interpolation ------------------------ */
/* Pressure as function of depth from erosion level calc */
if (presmod) depthp(x,y);
/*
printf("eps %e %e ",m1,m2,e,n); getchar();
*/
}
/* Calculation of P by Interpolation */



/* Calculation of P from the depth below the surface */
void depthp(double x, double y)
/* x,y - XY location of point for Vx,Vy calc */
{
/* Counters */
long int m10;
/* en-NormalisedDistance */
double e,n,ival;
/**/
/* Up Left Node X,Y Num */
m10=wn[0];
/**/
/* Check X,Y */
if(x<0) x=0; else if(x>xsize) x=xsize;
if(y<0) y=0; else if(y>ysize) y=ysize;
/**/
/**/
/**/
/* Pressure as function of depth from erosion level calc */
/* Relativ Normalized coord Calc */
e=(x-gx[m10])/(gx[m10+1]-gx[m10]);
n=(e*ep[m10+1]+(1.0-e)*ep[m10]);
ival=3300.0*(y-n)*9.80665;
if (ival<1e+5) ival=1e+5;
eps[10]=ival;
/*
*/
/*
printf("eps %e %e ",m1,m2,e); getchar();
*/
}
/* Calculation of P from the depth below the surface */



/* Melt fraction, latent heat calculation */
void meltpart1(double mtk, double mpb, long int mm1)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm1 - mark number */
/* yn  - type of calculation: 0 - Ro, 1 - Nu, 2 - Cp, 3 - kt */
{
/* Val buffer */
double xmelt=0,hlatent=0,ival,mpg,Fcpxout,Tcpxout;
long int m1;
double ykm=mpb*3e-3,ts=0,tl=0,tll=0,mwa=markw[mm1],xh2osat,xh2o,dt,dt0,dtmin,dtmax;
/* Marker type */
int n1,n2,mm2=(int)markt[mm1];
/**/
/**/
/**/
/* Calculate melt fraction using marker type */
if (ykm>0 && marky[mm1]<zdeep && marky[mm1]<520000.0)
switch(mm2)
	{
	/* Granits, Dacites, Sediments: latent heat 300 kJ/kg (Bittner & Schmeling, 1995) */
	case 15:
	case 35:
	case 16:
	case 17:
	case 36:
	case 37:
	case 3:
	case 4:
	case 5:
	case 23:
	case 24:
	case 25:
	/* Wet Solidus Temperature, Johannes, 1985, Poli & Schmidt, 2002 */
	if (ykm<36.0) 
		{
		ts=889.0+536.6/(ykm+1.609)+18.21/(ykm+1.609)/(ykm+1.609);
		}
	else
		{
		ts=831.3+2.0*ykm;
		}
	/* Dry Granite Liquidus, Johannes, 1985 */
	tl=1262.0+3.0*ykm;
	hlatent=300000.0;
	break;
	/**/
	/* Gabbro latent heat 380 kJ/kg (Bittner & Schmeling, 1995) */
	case 6:
	case 26:
	case 8:
	case 28: 
	/* Dry Toleitic Basalt solidus, Hess, 1989 */
	ts=1327.15+3.02*ykm;
	/* Dry Toleitic Basalt Liquidus, Hess, 1989 */
	tl=1423.15+3.5*ykm;
	hlatent=380000.0;
	break;
	/**/
	/* Basalt latent heat 380 kJ/kg (Bittner & Schmeling, 1995) */
	case 38:
	case 18:
	case 7:
	case 27:
	/* Wet solidus, Schmidt & Poli, 1998  */
	if (ykm<48.0) 
		{
		ts=972.6-2111.0/(ykm+10.63)+70033.0/(ykm+10.63)/(ykm+10.63);
		}
	else
		{
		ts=935.4+0.1162*ykm+0.006937*ykm*ykm;
		}
	/* Dry Toleitic Basalt Liquidus, Hess, 1989 */
	tl=1423.15+3.5*ykm;
	hlatent=380000.0;
	break;
	/**/
	/* Dry/Wet Peridotite: latent heat 400 kJ/kg Turcotte & Schubert, 1982, p.171 */
	case  9:
	case 10:
	case 11:
	case 12:
	case 14:
	case 29:
	case 30:
	case 32:
	case 34:
	/* Dry/wet mantle melting Katz et al., 2003 */
	mpg=mpb*1e-4;
	if(mwa<=0)
		{
		dtmin=dtmax=0;
		n2=1;
		}
	else
		{
		xh2osat=12.0*pow(mpg,0.6)+1.0*mpg;
		xh2o=mwa/(0.01+0.0*(1.0-0.01));
		if(xh2o>xh2osat) xh2o=xh2osat;
		dt=43.0*pow(xh2o,0.75);
		ts=273.15+1085.7+132.9*mpg-5.1*mpg*mpg-dt;
		n2=1;
		dtmin=dtmax=dt;
		if(mtk>=ts)
			{
			xh2osat=12.0*pow(mpg,0.6)+1.0*mpg;
			xh2o=mwa/(0.01+1.0*(1.0-0.01));
			if(xh2o>xh2osat) xh2o=xh2osat;
			dt=43.0*pow(xh2o,0.75);
			dtmin=dt;
			n2=100;
			}
		}
	for(n1=0;n1<n2;n1++)
		{
		dt=(dtmin+dtmax)/2.0;
		ts=273.15+1085.7+132.9*mpg-5.1*mpg*mpg-dt;
		tll=273.15+1475.0+80.0*mpg-3.2*mpg*mpg-dt;
		tl=273.15+1780.0+45.0*mpg-2.0*mpg*mpg-dt;
		xmelt=0;
		if(mtk>=ts)
			{
			if(mtk>=tl)
				{
				xmelt=1.0;
				}
			else
				{
				/* 0.15 = 15 wt% of Cpx in the mantle */
				Fcpxout=0.15/(0.5+0.08*mpg);
				Tcpxout=pow(Fcpxout,1.0/1.5)*(tll-ts)+ts;
				if(mtk>Tcpxout)
					{
					/* Opx melting */
					xmelt=Fcpxout+(1.0-Fcpxout)*pow((mtk-Tcpxout)/(tl-Tcpxout),1.5);
					}
				else
					{
					/* Cpx melting */
					xmelt=pow((mtk-ts)/(tll-ts),1.5);
					}
				/* Compute artificial Tliquidus */
				tl=(mtk-ts)/xmelt+ts;
				}
			}
		/* Water content in melt */
		if(mwa>0)
			{
			xh2osat=12.0*pow(mpg,0.6)+1.0*mpg;
			xh2o=mwa/(0.01+xmelt*(1.0-0.01));
			if(xh2o>xh2osat) xh2o=xh2osat;
			dt0=43.0*pow(xh2o,0.75);
			if(dt0>dt) 
				{
				dtmin=dt;
				}
			else
				{
				dtmax=dt;
				}
			}
		if(ABSV(dtmax-dtmin)<0.001) break;
		}
/*
if(mwa>0 && mtk>1300.0 && mpg>2.0){printf("%d %d %e %e %e %e %e    %e %e %e   %e %e %e   %e %e %e %e",n1,mm2,mtk-273.15,mpg,mwa,markw[mm1],markex[mm1],ts,tll,tl,xh2o,xh2osat,xmelt,dt,dt0,dtmin,dtmax);getchar();}
*/
/*
if(mwa>0){printf("%d %d %e %e %e   %e %e %e   %e %e %e   %e %e %e %e",n1,mm2,mtk-273.15,mpg,mwa,ts,tll,tl,xh2o,xh2osat,xmelt,dt,dt0,dtmin,dtmax);getchar();}
*/
	hlatent=400000.0;
	break;
	/**/
	/* Other rocks - No melting */
	default:
	break;
	}
/**/
/* Melt fraction, latent heat calculation */
eps[21]=eps[22]=0;
if(tl)
	{
	/* Melt fraction calc, check */
	if(ts>=tl-100.0) ts=tl-100.0;
	xmelt=(mtk-ts)/(tl-ts);
	if(xmelt<0) xmelt=0;
	if(xmelt>1.0) xmelt=1.0;
	eps[21]=xmelt;
	/* Latent heat calc */
	hlatent*=xmelt;
	eps[22]=hlatent;
/*
if(xmelt>0 && (mm2==10 || mm2==30 || mm2==9 || mm2==29)) {printf("Meltpart1: %d %e %e  %e %e %e",mm2,mtk,mpb,ts,tl,xmelt);getchar();}
if(mm2<20 && xmelt) {printf("Meltpart1: %d %e %e  %e %e",mm2,mtk,mpb,xmelt,hlatent);getchar();}
*/
	}
/**/
}
/* Melt fraction, latent heat calculation */



/* Antigorite dehydration depth */
double dehydration(int yn)
/* yn - Upper(0)/Lower boundary */
{
/* Val buffer */
double x0,x1,y0,y1,k0,k1,t0,t1,sy1,ival,yshift;
long int m1;
/**/
/**/
/**/
/* Set Upper/Lower boundary */
if(yn) yshift=5000.0; else yshift=-5000.0;
/**/
/* Serch for melting depth along lower hydration surface */
for (m1=1;m1<xnumx;m1++)
	{
	/* TK calc along the lower hydr surface */
	x0=gx[m1-1];
	y0=ep[m1+xnumx*3-1]+yshift;
	allintert(x0,y0); t0=eps[2];
	x1=gx[m1];
	y1=ep[m1+xnumx*3]+yshift;
	allintert(x1,y1); t1=eps[2];
	/* Oceanic crust top */
	sy1=ep[m1];
	if((y0-sy1)>0 && (y1-sy1)>0)
		{
		/* Atg stability field after Schmidt and Poli, 1998 */
		if((y0-sy1)>63000.0)
			{
			k0=1013.17699-0.060387633e-3*(y0-sy1)-0.004289442e-6*(y0-sy1)*(y0-sy1);
			}
		else
			{
			k0=751.490422+6.00773668e-3*(y0-sy1)-0.034690759e-6*(y0-sy1)*(y0-sy1);
			}
		if((y1-sy1)>63000.0)
			{
			k1=1013.17699-0.060387633e-3*(y1-sy1)-0.004289442e-6*(y1-sy1)*(y1-sy1);
			}
		else
			{
			k1=751.490422+6.00773668e-3*(y1-sy1)-0.034690759e-6*(y1-sy1)*(y1-sy1);
			}
		/* Antigorite stability depth calc */
		if(t1>=k1)
			{
			ival=y1-(y1-y0)/((t1-k1)+(k0-t0))*(t1-k1);
			/* Set Upper/Lower boundary */
			if(yn) dehydrmax=ival; else dehydrmin=ival;
			return 0;
			}
		}
	}
basalty=ysize;
if(yn) dehydrmax=ysize; else dehydrmin=ysize;
return 0;
}
/* Antigorite dehydration depth */



/* Nu calc after reological equation */
/* P-T-stress dependent rheology without/with brittle/ductile transition */
/* Reological equations */
/* Stress>SScr */
/* Power law dislocation creep: SSii={NU0*EEii*exp[(E+PV)/RT]}^(1/n) */
/* Effective viscosity: NU=1/2*{NU0*exp[(E+PV)/RT]}^(1/n)*EEii^[(1-n)/n] */
/* Stress<SScr */
/* Newtonian diffusion   creep: SSii=NU1*EEii*exp[(E+PV)/RT] */
/* Effective viscosity: NU=NU0/2*exp[(E+PV)/RT] */
/* NU1=NU0/SScr^(n-1) */
/* SScr - dislocation, diffusion transition stress */
/* SSii - second invariant of deviatoric stress tensor */
/* EEii - second invariant of strain rate tensor */
/* E - activation energy, J */
/* V - activation volume, J/bar */
/* R - gase constant 8.314 J/K */
/* Viscosity NU  calc after reological equations */
/* NU=SSii/(2*EEii) */
/* Brittle - Ductile transition */
/* sbrit=MINV(0.85e+5*pb,60e+6+0.6e+5*pb)*lambda;  (Schott & Scmeling, 1998) */
/* sbrit=MINV(0.667e+5*pb,51.2e+6+0.512e+5*pb)*lambda; (Brace & Kohlsstedt, 1980) */
double viscalc(double mtk, double mpb, double x, double y, long int mm1, int mm2)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm1 - Marker number */
/* mm2 - rock type */
{
/* Val buffer */
double e,n,rt=8.314*mtk,k1,e1,epsin,sduct,sbrit,nueff,smin,smax,nmin,nmax,strain,abrit,bbrit,nubrit,nunewt,nupowl;
/* Reological Eq par */
double lamb,mnu2,kfmelt,tkmelt,nupeierls,sigin,siginnew,sig0,A,eiicur;
/* Counters */
long int m1;
int n1,n2;
/**/
/**/
/**/
/* Calc effective strain rate after second strain rate Tenzor invariant EEii=(1/2SUM(EPSik^2))^(1/2) */
epsin=markeii;
/*
if (epsin==0) epsin=1e-17;
*/
sigin=marksii;
/**/
/* Calc Second strain Tenzor invariant GAMMAii=(1/2SUM(EPSik^2))^(1/2) */
strain=ABSV(marke[mm1]);
/**/
/* Brittle "viscosity" calc: sbrit=A(eps)+LAM*B(eps)*P ------------------- */
/* Brittle/ductile transition stress calc */
/* A,B coefficients calc depending on integral strain */
lamb=1.0; 
abrit=marka0[mm2]; 
bbrit=markb0[mm2]; 
if(strain>marke1[mm2])
	{
	abrit=marka1[mm2]; 
	bbrit=markb1[mm2];
	}
else
	{
	if(strain>marke0[mm2] && marke1[mm2]>marke0[mm2]) 
		{
		abrit=marka0[mm2]+(marka1[mm2]-marka0[mm2])*(strain-marke0[mm2])/(marke1[mm2]-marke0[mm2]);
		bbrit=markb0[mm2]+(markb1[mm2]-markb0[mm2])*(strain-marke0[mm2])/(marke1[mm2]-marke0[mm2]);
		}
	}
/* No friction in fluid/melt present areas */
if(markx[mm1]>6e+5)
	{
	if(markv[mm1]==1) lamb=lambfld;
	if(markv[mm1]==2) lamb=lambmlt;
	if(markv[mm1]==3) lamb=MINV(lambfld,lambmlt);
	}
/* Brittle strength calc */
sbrit=abrit+bbrit*mpb*1e+5*lamb;
if(sbrit>nuplast) sbrit=nuplast;
/**/
/* Brittle "viscosity" calc: sbrit=A(eps)+LAM*B(eps)*P ------------------- */
/* Inverted value of Brittle "Mohr-Coulomb" viscosity calc */
/*
printf("%ld %d %e %e %e %e %e",mm1,mm2,mtk,mpb,sigin,epsin,sbrit);getchar();
*/
nubrit=0;
if(epsin>0 && (bbrit>0 || abrit*lamb>0))
	{
	if(sbrit>0)
		{
		/* Calc,check plastic strain rate */
		eiicur=epsin*markrat;
		nubrit=1.0/(0.5*sbrit/eiicur);
		}
	else
		{
		nubrit=1.0/markn0[mm2];
		}
	}
/* End Brittle "viscosity" calc: sbrit=A(eps)+LAM*B(eps)*P ------------------- */
/**/
/**/
/**/
/* Peierls plasticity-creep mechanism, data from Katayama and Karato, 2008 */
nupeierls=0;
if(1==0 && mm2>1 && mm2<20 && sigin>1e+7 && mtk<1473.0 && epsin>0) 
/*
if(mm2>=9 && mm2<=14)
if(0==0 && mm2>1 && mm2<20 && sigin>0 && epsin>0) 
*/
	{
	A=pow(10.0,7.8)*1e-12;/* Constant f(p,q), 1/s/MPa^2 */ 
	sig0=9.1e+9; /* Dry Peierls stress at 0 K f(p,q), MPa Evans & Goetze, 1979 */
	n1=9;
	if(0==1 && mm2!=9 && mm2!=10 && mm2!=14)
		{
		sig0=2.9e+9; /* Wet Peierls stress at 0 K f(p,q), MPa Katayama & Karato, 2008 */
		n1=11;
		}
	siginnew=sigin;
	for(n2=0;n2<5;n2++)
		{
		siginnew=(1.0-pow(-rt/(markdh[n1]+markdv[n1]*mpb)*log(epsin/A/siginnew/siginnew),0.5))*sig0;
/*
	p=1.0;
	q=2.0;
	k1=A*pow(siginnew,2.0)*exp(-(markdh[n1]+markdv[n1]*mpb)/rt*pow(1.0-pow(siginnew/sig0,p),q));
printf("%ld %d   %e %e   %e %e   %e %e",mm1,mm2,mtk,mpb,sigin,epsin,siginnew,k1);getchar();
printf("%ld %d   %e %e   %e %e   %e %e",mm1,mm2,mtk,mpb,sigin,epsin,siginnew,1.0/nupeierls);getchar();
*/
		}
		
	nupeierls=1.0/(0.5*siginnew/epsin);
	}
/**/
/**/
/**/
/* Ductile viscosity calc -------------------------------------------*/
/* Inverted value of newtonian NU set */
nunewt=0;
/**/
/* Inverted value of power-low NU set */
nupowl=0;
/**/
/* Check for the presence of ductile rheology */
if (marknu[mm2])
	{
	/* A)  Simple Newtonian rheology */
	/* Newtonian creep: SSii=NU0*2.0*EEii */
	/* Effective viscosity: NU=NU0 */
	/* Effective viscosity member in Stoks: NUs=NU */
	if(markdh[mm2]==0 && markdv[mm2]==0 && (markss[mm2]==0 || markmm[mm2]==1.0))
		{
		/* Inverted value of newtonian NU calc */
		nunewt=1.0/marknu[mm2];
		}
	/**/
	/**/
	/**/
	/* B)  P-T dependent, stress independent Newtonian rheology */
	/* Newtonian diffusion creep: SSii=NU0*EEii*exp[(E+PV)/RT] */
	/* Effective viscosity: NU=NU0*exp[(E+PV)/RT] */
	if((markdh[mm2]!=0 || markdv[mm2]!=0) && (markss[mm2]==0 || markmm[mm2]==1.0))
		{
		/* Inverted value of newtonian NU calc */
		e1=(markdh[mm2]+markdv[mm2]*mpb)/rt;
		if(e1>150.0) e1=150.0;
		/* Test creep Moresi & Solomatov (1995): SSii=NU0*exp[-a*(T-T0)] */
		if(markdh[mm2]<0 && markdv[mm2]>0) 
			{
			e1=markdh[mm2]*(mtk-markdv[mm2]);
			if(e1<-150.0) e1=-150.0;
			}
		/* Test creep Turkotte & Schubert(1982): SSii=NU0*exp[E/RTo(1-(T-T0)/T0)] */
		if(markdh[mm2]<0 && markdv[mm2]<0) 
			{
			e1=(-markdh[mm2])*(1.0-(mtk-(-markdv[mm2]))/(-markdv[mm2]))/8.314/(-markdv[mm2]);
			if(e1>150.0) e1=150.0;
			}
		nunewt=1.0/(marknu[mm2]*exp(e1));
		}
	/**/
	/**/
	/**/
	/* C)  P-T independent, stress dependent rheology without/with brittle/ductile transition */
	/* Stress>SScr */
	/* Power law creep: SSii={NU0*EEii}^(1/n) */
	/* Effective viscosity: NU=1/2*NU0^(1/n)*EEii^[(1-n)/n] */
	/* Effective viscosity member in Stoks: NUs=NU/n */
	/* Stress<SScr */
	/* Newtonian creep: SSii=NU1*EEii */
	/* Effective viscosity: NU=NU1/2 */
	/* Effective viscosity member in Stoks: NUs=NU */
	/* NU1=NU0/SScr^(n-1) */
	if((markdh[mm2]==0 && markdv[mm2]==0) && markss[mm2]!=0 && markmm[mm2]!=1.0)
		{
		/* Koef for stress independent creep NU1 calc */
		k1=marknu[mm2]/pow(markss[mm2],markmm[mm2]-1.0);
		/**/
		/* Inverted value of newtonian NU calc */
		nunewt=1.0/(0.5*k1);
		/**/
		/* Ductile power-low stress calc */
		sduct=pow(marknu[mm2]*epsin,1.0/markmm[mm2]);
		/**/
		/* Inverted value of power-low NU calc */
		if (epsin>0) nupowl=1.0/(0.5*sduct/epsin);
		}
	/**/
	/**/
	/**/
	/* D)  P-T-stress dependent rheology without/with brittle/ductile transition */
	/* Reological equations */
	/* Stress>SScr */
	/* Power law dislocation creep: SSii={NU0*EEii*exp[(E+PV)/RT]}^(1/n) */
	/* Effective viscosity: NU=1/2*{NU0*exp[(E+PV)/RT]}^(1/n)*EEii^[(1-n)/n] */
	/* Effective viscosity member in Stoks: NUs=NU/n */
	/* Stress<SScr */
	/* Newtonian diffusion   creep: SSii=NU1*EEii*exp[(E+PV)/RT] */
	/* Effective viscosity: NU=NU0/2*exp[(E+PV)/RT] */
	/* Effective viscosity member in Stoks: NUs=NU */
	/* NU1=NU0/SScr^(n-1) */
	if(marknu[mm2]>0 && (markdh[mm2]>0 || markdv[mm2]!=0) && markss[mm2]!=0 && markmm[mm2]!=1.0)
		{
		/* T-P exponent for effective NU calc */
		e1=(markdh[mm2]+markdv[mm2]*mpb)/rt;
		if(e1>150.0) e1=150.0;
		e1=exp(e1);
		/**/
		/* Koef for stress independent creep NU1 calc */
		k1=marknu[mm2]/pow(markss[mm2],markmm[mm2]-1.0);
		/**/
		/* Inverted value of newtonian NU calc */
		nunewt=1.0/(0.5*k1*e1);
		/**/
		/* Ductile power-low stress calc */
		sduct=pow(marknu[mm2]*e1*epsin,1.0/markmm[mm2]);
		/**/
		/* Inverted value of power-low NU calc */
		if(epsin>0) nupowl=1.0/(0.5*sduct/epsin);
		}
	/* F)  P-T-stress dependent rheology of Mantle Karato & Wu (1993) */
	if(marknu[mm2]<0)
		{
		/* Diffusion viscosity calc for given grainsize, mm, Karato & Wu (1993) */
		mnu2=1.0; /* Grainsize, mm */ 
		nunewt=1.0/(0.5*pow(mnu2*1e-3/5e-10,2.5)*8e+10/marks1[mm2]*exp((markss[mm2]+mpb*markll[mm2])/rt));
		/* Dislocation viscosity calc, Karato & Wu (1993) */
		/* Effective viscosity1 calc */
		nupowl=0; 
		if(epsin>0) nupowl=1.0/(0.5*pow(epsin,1.0/markmm[mm2]-1.0)*8e+10/pow(-marknu[mm2],1.0/markmm[mm2])*exp((markdh[mm2]+markdv[mm2]*mpb)/rt/markmm[mm2]));
		}
	}
/* End Ductile viscosity calc -------------------------------------------*/
/**/
/**/
/**/
/* Inverted value of effective viscosity calc, check */
nueff=nunewt+nupowl;
if(nupeierls>nueff) nueff=nupeierls;
if(nubrit>nueff)
	{
	nueff=nubrit;
	/* Mark plastic yeilding */
	if(marke[mm1]==0) marke[mm1]=1e-20;
	marke[mm1]=ABSV(marke[mm1]);
	}
else
	{
	/* Stop plastic yeilding */
	marke[mm1]=-ABSV(marke[mm1]);
	}
/*
printf("%ld %d  %e %e  %e %e  %e %e %e %e",mm1,mm2,x,y,mtk,mpb,nubrit,nunewt,nupowl,nueff);getchar();
*/
/**/
/* Inverted Viscosity check */
if(nueff<=0) nueff=1.0/markn1[mm2];
/**/
/* Viscosity calc */
nueff=1.0/nueff;
/*
printf("A %ld %d  %e %e      %e %e  %e %e",mm1,mm2,markx[mm1],marky[mm1],mtk,mpb,epsin,nueff);getchar();
printf("A %ld %d %e %e %e",mm1,mm2,mtk,mpb,nueff);getchar();
*/
/**/
/* Central dike */
/*
if(mm2<20 && mm2!=10 && nueff>1e+20 && ABSV(markx[mm1])<5000.0 && marky[mm1]<50000.0) nueff=1e+20;
*/
/* Viscosity check */
if(nueff<markn0[mm2]) nueff=markn0[mm2]; 
if(nueff>markn1[mm2]) nueff=markn1[mm2];
if(nueff<nubeg) nueff=nubeg; 
if(nueff>nuend) nueff=nuend;
/**/
/* Melting in 1000-2000 km interval only */
/*
if ((markx[mm1]<1e+6 || markx[mm1]>3e+6) && nueff<1e+20 && mm2>1) nueff=1e+20;
*/
/**/
/* Lower boundary */
/*
if((x<200e+3 || x>xsize-200e+3) && y<200e+3 && nueff>1e+23) nueff=1e+23;
if(y>1.3e+6 && nueff<1e+24) nueff=1e+24;
if(strain>1.5 && nubrit>0 && x>1000000.0 && y<400000.0) {printf("%ld %d  %e %e  %e %e  %e %e %e %e  %e %e %e",mm1,mm2,x,y,mtk,mpb,strain,abrit,bbrit,lamb,1.0/nubrit,epsin,nueff); getchar();}
*/
/* Return calculated viscosity */
/*
printf("B   %ld %d %e %e %e",mm1,mm2,mtk,mpb,nueff);getchar();
if(mm2==12){printf("%ld %e %e %e %e %e %e",mm1,mtk,mpb,epsin,sbrit,1.0/nubrit,nueff);getchar();}
*/
return nueff;
}
/* Nu calc after reological equation */





/* Hydration front progress after H2O budget */
double hydration2()
{
/* Val buffer */
double ysurf,vfiltr,yfiltr,dydx,dydx1,sy1,sy2,sy3,sy4,sy5,e1,mwamin,x0,y0,x1,y1,vx1,vy1;
double hytimesum,hytimesum0;
/* TD Database variables */
double W0,W1,W2,W3,R0,R1,R2,R3,n,e,dx,dy;
double mtk,mpb,mwa,mro,dmwa,wro;
long int m1,m2,m3,mm1,marknum1=marknum;
int mm2,mm3,n1,n2;
/**/
printf("\n WATER Transport BEG \n");
/* Marker steps */
dx=dxwater;
dy=dywater;
/**/
/**/
/* Min water contents in the hydraten mantle wt% */
mwamin=0.1;
/* Min Distance from erosion surface for water release */
ysurf=8000.0;
/**/
/**/
/* Clear wa[] wt */
for (m1=0;m1<nodenum;m1++)
	{
	wa0[m1]=0;
	wa1[m1]=0;
	sol0[m1]=0;
	sol1[m1]=0;
	sol0[nodenum+m1]=1e+30;
	sol1[nodenum+m1]=-1e+30;
	sol0[nodenum2+m1]=1e+30;
	sol1[nodenum2+m1]=-1e+30;
	fre0[         m1]=1e+30;
	fre0[nodenum +m1]=-1e+30;
	fre0[nodenum2+m1]=1e+30;
	fre0[nodenum3+m1]=-1e+30;
	}
/**/
/**/
/**/
/* Fluid marker generation cycle */
for (mm1=0;mm1<marknum;mm1++)
{
/* Marker type */
mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
/**/
/* Marker cell number */
m1=m1serch((double)markx[mm1]);
m2=m2serch((double)marky[mm1]);
m3=m1*ynumy+m2;
/**/
/*Fluid disappearance surface */
sy1=sedilev;
/**/
/*
e1=(markx[mm1]-gx[m1])/(gx[m1+1]-gx[m1]);
printf("%ld %e %e %e %e %e %e",mm1,markx[mm1],marky[mm1],e1,sy1,sy2,sy3);getchar();
*/
/**/
/* Check markers out of grid and within hydration range */
if(markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && (markk[mm1]>0 || markt[mm1]>=50) && markt[mm1]<100)
if((double)(markd[mm1])>=0 && (double)(markw[mm1])>=0 && mm2>1 && mm2!=9 && mm2!=10 && mm2!=5 && mm2!=6 && mm2!=25 && mm2!=26)
	{
	if(mm2<50)
		{
		/* P, T parameters calc */
		allinterp((double)(markx[mm1]),(double)(marky[mm1]));
		mpb=eps[10]*1e-5;
		mtk=(double)(markk[mm1]);
		/**/
		/* Mantle to Antigorite transformation */
/*
		antigor(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1);
*/
		/**/
		/* Rocks to rock+melt transformation */
/*
		melting(mtk,mpb,mm1);
*/
		if (markt[mm1]>=20)
			{
			/* Check melting extent */
			if(fre0[        +m3]>markx[mm1]-dx) fre0[         m3]=markx[mm1]-dx;
			if(fre0[nodenum +m3]<markx[mm1]+dx) fre0[nodenum +m3]=markx[mm1]+dx;
			if(fre0[nodenum2+m3]>marky[mm1]-dy) fre0[nodenum2+m3]=marky[mm1]-dy;
			if(fre0[nodenum3+m3]<marky[mm1]+dy) fre0[nodenum3+m3]=marky[mm1]+dy;
			}
		/* Compute TD variables */
/*
		tdbasecalc(mtk,mpb,mm2,mm1);
		mro=eps[41];
		mwa=eps[42];
*/
		watercalc(mtk,mpb,mm2,mm1);
		mro=markd[mm1];
		mwa=eps[42];
		/**/
		/* Water changes in kg/m3 calc */
		dmwa=(1.0-markwt[mm1])*(1.0-markex[mm1])*mro*(mwa-markw[mm1])*1e-2;
/*
{printf("H2O MARKER %ld %d %d %e %e   %e %e  %e %e   %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa);getchar();}
{printf("H2O RELEASE %ld %d %d %e %e   %e %e  %e %e   %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa);getchar();}
*/
		/**/
		/* Add water changes to the current cell, kg/m3 */
		/* Water release */
		if ((markw[mm1]-mwa)>dmwamin)
			{
			/* Save new water content */
			markw[mm1]=mwa;
			/* Generation of fluid marker (NO FLUID From melts */
			if (markt[mm1]<20 && marky[mm1]>sy1 && marky[mm1]<zdeep)
				{
				markt[marknum1]=markt[mm1]+50;
				markx[marknum1]=markx[mm1];
				marky[marknum1]=marky[mm1];
				markk[marknum1]=markk[mm1];
				markd[marknum1]=1050.0;
				markw[marknum1]=-dmwa;
				markdi[marknum1]=0;
				markwt[marknum1]=0;
				markex[marknum1]=0;
				marke[marknum1]=0;
				markv[marknum1]=0;
				/* Add aditional markers counter */
				marknum1++;
				/* Check hydration extent */
				if(sol0[nodenum+m3]>markx[mm1]-dx) sol0[nodenum+m3]=markx[mm1]-dx;
				if(sol1[nodenum+m3]<markx[mm1]+dx) sol1[nodenum+m3]=markx[mm1]+dx;
				if(sol0[nodenum2+m3]>marky[mm1]-dy) sol0[nodenum2+m3]=marky[mm1]-dy;
				if(sol1[nodenum2+m3]<marky[mm1]+dy) sol1[nodenum2+m3]=marky[mm1]+dy;
				}
			}
		else
		/* Water consuming */
			{
			if(dmwa>0)
				{
				wa1[m3]+=dmwa;
				sol1[m3]+=1.0;
				}
			}
		}
	else
	/* Fluid marker count */
		{
		/* Check position */
		if(marky[mm1]>sy1 && marky[mm1]<zdeep)
			{
			/* Check hydration extent */
			if(sol0[nodenum+m3]>markx[mm1]-dx) sol0[nodenum+m3]=markx[mm1]-dx;
			if(sol1[nodenum+m3]<markx[mm1]+dx) sol1[nodenum+m3]=markx[mm1]+dx;
			if(sol0[nodenum2+m3]>marky[mm1]-dy) sol0[nodenum2+m3]=marky[mm1]-dy;
			if(sol1[nodenum2+m3]<marky[mm1]+dy) sol1[nodenum2+m3]=marky[mm1]+dy;
			}
		else
		/* Erase fluid marker */
			{
			markx[mm1]=-xsize;
			markk[mm1]=0;
			}
		}
	}
}
/**/
/**/
/*
for(m1=0;m1<xnumx-1;m1++)
for(m2=0;m2<ynumy-1;m2++)
{
m3=m1*ynumy+m2;
printf("%ld %ld  %e %e  %e %e  %e %e   %e %e %e %e",m1,m2,gx[m1],gx[m1+1],gy[m2],gy[m2+1],wa0[m3],wa1[m3],sol0[nodenum+m3],sol1[nodenum+m3],sol0[nodenum*2+m3],sol1[nodenum*2+m3]);getchar();
}
*/
/**/
/* Rock hydration cycle */
for (mm1=0;mm1<marknum;mm1++)
if(markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && markt[mm1]<50)
{
/*
printf("%ld %e %e %e %e %e %e",mm1,markx[mm1],marky[mm1],e1,sy1,sy2,sy3);getchar();
*/
/**/
/* Marker cell number */
m1=m1serch((double)markx[mm1]);
m2=m2serch((double)marky[mm1]);
m3=m1*ynumy+m2;
/* Check markers within hydration range */
if(markx[mm1]>sol0[nodenum+m3] && marky[mm1]>sol0[nodenum2+m3] && (double)(markx[mm1])<sol1[nodenum+m3] && (double)(marky[mm1])<sol1[nodenum2+m3])
	{
	/* Fluid presence mark */
	markv[mm1]+=1;
/*
*/
if(markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==12 || markt[mm1]==14)
	{
	/* Mantle Hydration */
	mm2=markt[mm1]=11;
	/* P, T parameters calc */
	allinterp((double)(markx[mm1]),(double)(marky[mm1]));
	mpb=eps[10]*1e-5;
	mtk=(double)(markk[mm1]);
	/**/
	/* Mantle to Antigorite transformation */
/*
	antigor(mtk,mpb,(double)(markx[mm1]),(double)(marky[mm1]),mm1);
*/
	/**/
	/* Rocks to rock+melt transformation */
/*
	melting(mtk,mpb,mm1);
*/
	if (markt[mm1]>=20)
		{
		/* Check melting extent */
		if(fre0[        +m3]>markx[mm1]-dx) fre0[         m3]=markx[mm1]-dx;
		if(fre0[nodenum +m3]<markx[mm1]+dx) fre0[nodenum +m3]=markx[mm1]+dx;
		if(fre0[nodenum2+m3]>marky[mm1]-dy) fre0[nodenum2+m3]=marky[mm1]-dy;
		if(fre0[nodenum3+m3]<marky[mm1]+dy) fre0[nodenum3+m3]=marky[mm1]+dy;
		}
	/**/
	/* Thermodynamic database use for Ro, Water */
	/* Compute TD variables */
/*
	tdbasecalc(mtk,mpb,mm2,mm1);
	mro=eps[41];
	mwa=eps[42];
*/
	watercalc(mtk,mpb,mm2,mm1);
	mro=markd[mm1];
	mwa=eps[42];
	/**/
	/* Water changes in kg/m3 calc */
	dmwa=(1.0-markwt[mm1])*(1.0-markex[mm1])*mro*(mwa-markw[mm1])*1e-2;
/*
{printf("H2O HYDRATE %ld %d %d %e %e   %e %e  %e %e   %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa);getchar();}
*/
	/**/
	/* Add water changes to the current cell, kg/m3 */
	/* Water consuming */
	if (dmwa>0)
		{
		wa1[m3]+=dmwa;
		sol1[m3]+=1.0;
		}
	}
	}
}
/**/
/**/
/*
for(m1=0;m1<xnumx-1;m1++)
for(m2=0;m2<ynumy-1;m2++)
{
m3=m1*ynumy+m2;
printf("%ld %ld  %e %e  %e %e  %e %e   %e %e %e %e",m1,m2,gx[m1],gx[m1+1],gy[m2],gy[m2+1],wa0[m3],wa1[m3],sol0[nodenum+m3],sol1[nodenum+m3],sol0[nodenum*2+m3],sol1[nodenum*2+m3]);getchar();
}
*/
/**/
/**/
/**/
/* Fluid marker computing cycle */
for (mm1=0;mm1<marknum1;mm1++)
{
/**/
/*
printf("%ld %e %e %e %e %e %e",mm1,markx[mm1],marky[mm1],e1,sy1,sy2,sy3);getchar();
*/
/**/
/* Check markers out of grid and within hydration range */
if(markt[mm1]>=50 && markt[mm1]<100 && markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize)
	{
	/* Marker cell number */
	m1=m1serch((double)markx[mm1]);
	m2=m2serch((double)marky[mm1]);
	m3=m1*ynumy+m2;
	/**/
	/*Fluid disappearance surface */
	sy1=sedilev;
	/* Water in melt region conversion */
	if(markd[mm1]<1100.0 && markx[mm1]>fre0[m3] && marky[mm1]>fre0[nodenum2+m3] && markx[mm1]<fre0[nodenum+m3] && marky[mm1]<fre0[nodenum3+m3]) markd[mm1]=1150.0;
/*
if(markd[mm1]>1100.0) {printf("H2O REDUCE1 %ld %d   %ld %ld %ld   %e %e    %e %e %e",mm1,mm2,m1,m2,m3,wa0[m3],wa1[m3],markx[mm1],markw[mm1],markd[mm1]);getchar();}
*/
	/* Check position, no fluid above erosion/sedimentation level */
/*
	if(marky[mm1]>sy1 && marky[mm1]<zdeep && (markd[mm1]<1100.0 || (markx[mm1]>fre0[m3] && marky[mm1]>fre0[nodenum2+m3] && markx[mm1]<fre0[nodenum+m3] && marky[mm1]<fre0[nodenum3+m3])))
*/
	if(marky[mm1]>sy1 && marky[mm1]<zdeep)
		{
		wa0[m3]+=markw[mm1];
		sol0[m3]+=1.0;
		}
	else
	/* Erase fluid marker */
		{
		markx[mm1]=-xsize;
		markk[mm1]=0;
		}
/*
if(markd[mm1]>1100.0) {printf("H2O REDUCE2 %ld %d   %ld %ld %ld   %e %e    %e %e %e",mm1,mm2,m1,m2,m3,wa0[m3],wa1[m3],markx[mm1],markw[mm1],markd[mm1]);getchar();}
*/
	}
}
/**/
/**/
/**/
/*
for(m1=0;m1<xnumx-1;m1++)
for(m2=0;m2<ynumy-1;m2++)
{
m3=m1*ynumy+m2;
printf("%ld %ld  %e %e  %e %e  %e %e   %e %e %e %e",m1,m2,gx[m1],gx[m1+1],gy[m2],gy[m2+1],wa0[m3],wa1[m3],sol0[nodenum+m3],sol1[nodenum+m3],sol0[nodenum*2+m3],sol1[nodenum*2+m3]);getchar();
}
*/
/* Fluid marker consuming cycle */
for (mm1=0;mm1<marknum1;mm1++)
{
/* Marker type */
mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
/**/
/* Marker cell number */
m1=m1serch((double)markx[mm1]);
m2=m2serch((double)marky[mm1]);
m3=m1*ynumy+m2;
/**/
/*
printf("%ld %e %e %e %e %e %e",mm1,markx[mm1],marky[mm1],e1,sy1,sy2,sy3);getchar();
if(mm1>marknum){printf("%ld %d  %e %e  %e %e ",mm1,mm2,markx[mm1],marky[mm1],e1,sy1);getchar();}
*/
/**/
/* Change water consuming rocks  and fluid makers */
if(markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && (markk[mm1]>0 || markt[mm1]>=50) && markt[mm1]<100)
if((double)(markd[mm1])>=0 && (double)(markw[mm1])>=0 && mm2>1 && mm2!=9 && mm2!=10 && mm2!=12 && mm2!=14 && mm2!=5 && mm2!=6 && mm2!=25 && mm2!=26)
	{
	if(mm2<50)
		{
		/* P, T parameters calc */
		allinterp((double)(markx[mm1]),(double)(marky[mm1]));
		mpb=eps[10]*1e-5;
		mtk=(double)(markk[mm1]);
		/**/
		/* Thermodynamic database use for Ro, Water */
		/* Compute TD variables */
/*
		tdbasecalc(mtk,mpb,mm2,mm1);
		mwa=eps[42];
*/
		watercalc(mtk,mpb,mm2,mm1);
		mwa=eps[42];
		/**/
		/* Water changes in kg/m3 calc */
		dmwa=mwa-markw[mm1];
/*
{printf("TD! %ld %d %d %e %e %e %e ",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro);getchar();}
{printf("TDa %d %d %e %e %e %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,hpor,ppor,tpor,mwa);getchar();}
{printf("TDb %d %d %e %e %e %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,hpor,ppor,tpor,mwa);getchar();}
*/
		/**/
		/* Add water changes to the current cell, kg/m3 */
		/* Water consuming */
		if(dmwa>0)
			{
			if (wa1[m3]<=wa0[m3])
				{
				/* Save complete new water content */
				markw[mm1]=mwa;
				}
			else
				{
				/* COmpute, Save partial new water content */
				markw[mm1]=markw[mm1]+dmwa*wa0[m3]/wa1[m3];
				}
/*
{printf("H2O CONSUME %ld %d %d %e %e   %e %e  %e %e   %e    %e %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa,wa0[m3],wa1[m3]);getchar();}
*/
			}
		}
	else
	/* Fluid marker change */
		{
		if(wa1[m3]<wa0[m3])
			{
			/* Count water changes for fluid marker */
			markw[mm1]*=1.0-wa1[m3]/wa0[m3];
			}
		else
		/* Erase fluid marker */
			{
			markx[mm1]=-xsize;
			markk[mm1]=0;
			}
/*
if(markd[mm1]>1100.0) {printf("H2O REDUCE2 %ld %d   %ld %ld %ld   %e %e    %e %e %e",mm1,mm2,m1,m2,m3,wa0[m3],wa1[m3],markx[mm1],markw[mm1],markd[mm1]);getchar();}
*/
		}
	}
}
/*
marknum=marknum1;
return 0;
*/
/**/
/**/
/**/
/* Reset aditional markers */
if(printmod) printf("\n WATER BEG Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);
mm1=0;
while(marknum1>marknum && mm1<marknum)
	{
	/* Reload marker */
	if((markx[mm1]<0 || marky[mm1]<0 || (double)(markx[mm1])>xsize || (double)(marky[mm1])>ysize) && markt[mm1]<100 && markt[mm1]>1) 
		{
		/* Decrease aditional markers counter */
		marknum1--;
		if(markx[marknum1]>=0);
			{
			/* Type save */
			markt[mm1]=markt[marknum1];
			/* X,Y, water reload */
			markx[mm1]=markx[marknum1];
			marky[mm1]=marky[marknum1];
			markw[mm1]=markw[marknum1];
			markd[mm1]=markd[marknum1];
			markk[mm1]=markk[marknum1];
			markdt[mm1]=markdt[marknum1];
			markpr[mm1]=markpr[marknum1];
			markwt[mm1]=markwt[marknum1];
			markex[mm1]=markex[marknum1];
			marke[mm1]=marke[marknum1];
			}
		}
	/* Increase markers counter */
	mm1++;
	}
if(printmod) printf("\n WATER END Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);
/* Set new marker number */
marknum=marknum1;
/**/
/**/
/**/
return 0;
}
/* Hydration front progress after H2O budget */



/* Melt extraction */
void meltextract()
{
/* Counters, Buffers */
long int m1,m2,m3,m4,m5,mm1,m5max,marknum1=marknum;
int mm2,mm3;
double mtk,mpb,e,markvol=mardx*mardy,dx1,dx2,dx3,dy,y0,y1,yeros,xvolkcur,divcur,divmax,ycell,ycell0,xmelt,rocell,tkcur;
double xsolidcur,xmeltcur,H2Osolid,H2Omelt,KDmeltsolid;
/**/
/**/
/**/
/* Save old position of topography */
for (m1=0;m1<xnumx;m1++) 
	{
	/* Old topography position */
	ep0[m1]=ep[m1];
	/* Extraction column reset */
	if(m1<xnumx-1)
		{
		ep0[xnumx+  m1]=gx[m1+1];
		ep0[xnumx*2+m1]=gx[m1];
		ep0[xnumx*3+m1]=0;
		}
	}
/**/
/**/
/**/
if(printmod) printf("\n MELT EXTRACTION BEG Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);
/* Melt extraction cycle */
for (mm1=0;mm1<marknum;mm1++)
if(markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && (double)(marky[mm1])<zdeep && markk[mm1]>0 && markt[mm1]>=20 && markt[mm1]<50)
{
/* Marker type */
mm2=(int)markt[mm1];
/**/
/* P, T parameters calc */
allinterp((double)(markx[mm1]),(double)(marky[mm1]));
mpb=eps[10]*1e-5;
mtk=(double)(markk[mm1]);
/**/
/* Melt fraction computing */
meltpart1(mtk,mpb,mm1);
xmelt=eps[21];
/**/
/* Melt extraction computing */
/* No secondary  extraction */
if (mm2<=34)
/*
*/
/*
if (xmelt>0 && (xmelt-markex[mm1])>(minmelt+stpmelt) && markwt[mm1]<=0)
*/
/*
if(markx[mm1]>1e+6 && markx[mm1]<3e+6 && xmelt>0) {printf("A %ld %d %e %e %e %e",mm1,mm2,markx[mm1],marky[mm1],xmelt,markex[mm1]);getchar();}
*/
if ((xmelt-markex[mm1])>(minmelt+stpmelt))
	{
	/* No crust generation for the early stages  and model margins */
	if(timesum>1e+11 && markx[mm1]>coolzone && markx[mm1]<xsize-coolzone)
		{
/*
printf("A %ld %d %e %e %e %e",mm1,mm2,markx[mm1],marky[mm1],xmelt,markex[mm1]);getchar();
*/
		/* New marker generation */
		/* New marker type */
		switch (mm2)
			{
			/**/
			/* Sediments, Continental crust */
			case 22:
			case 23:
			case 24:
			case 25:
			case 35:
			case 36:
			case 37:
			mm3=15; break;
			/**/
			/* Basalts */
			case 27:
			mm3=16; break;
			/**/
			/* Gabbro  */
			case 26:
			case 28:
			case 38:
			mm3=17; break;
			/**/
			/* Dry Peridotite  */
			case 29:
			case 30:
			case 32:
			mm3=18; break;
			/**/
			/* Wet Peridotite  */
			case 34:
			mm3=18; break;
			}
		/**/
		/* Old Marker convergence */
		markdi[mm1]=-(xmelt-markex[mm1]-minmelt)/timestep*markvol;
/*
		markdi[mm1]=-(xmelt-markex[mm1]-minmelt)/timestep;
if(dv[m3]!=0){printf("%ld %ld %e",m1,m2,dv[m3]);getchar();}
		printf("%ld %d %e %e %e",mm1,markt[mm1],markx[mm1],marky[mm1],markdi[mm1]);getchar();
*/
		/**/
		/* Up Left Node X,Y Num */
		m1=wn[0];
		m2=wn[1];
		/* Relative Normalized coord Calc */
		e=((double)(markx[mm1])-gx[m1])/(gx[m1+1]-gx[m1]);
		/* Surface level elevation for marker definition */
		yeros=(e*ep[m1+1]+(1.0-e)*ep[m1]);
		/**/
		/* Split melt for volcanic and plutonic */
		/* Add plutonics to the depth */
		xvolkcur=1.0;
		if(xvolk<1.0)
			{
			/* Search for the weakest cell */
			divmax=-1e+30;
			m5max=0;
			/* Marker pressure cell */
			m3=(m1+1)*ynumy+m2+1;
			ycell0=(gy[m2]+gy[m2+1])/2.0;
			for (m5=0;m5<m2;m5++)
				{
				/* Current pressure cell */
				m4=(m1+1)*ynumy+m5+1;
				/* Vertical position of the cell center */
				ycell=(gy[m5]+gy[m5+1])/2.0;
				/* Density in the cell center */
				rocell=(ro[m4]+ro[m4-1]+ro[m4-ynumy]+ro[m4-1-ynumy])/4.0;
/*
{printf("A %ld %ld %ld  %e %e %e %e    %e %e %e ",m1,m2,m5,ycell,yeros,yintmin,yintmax,rocell,rointmin,rointmax);getchar();}
*/
				if(ycell>yeros+yintmin && ycell<yeros+yintmax && rocell>rointmin && rocell<rointmax)
					{
					/* Compute current possible divergence */
					divcur=(pr[m3]-GYKOEF*markro[mm2]*(ycell0-ycell)-pr[m4])/nd[m4];
/*
{printf("B %ld %ld %ld  %e %e %e %e    %e %e %e ",m1,m2,m5,ycell,yeros,yintmin,yintmax,rocell,rointmin,rointmax);getchar();}
*/
					/* Check the divergence */
					if(divcur>divmax)
						{
						divmax=divcur;
						m5max=m5;
						}
					}
				}
/*
m4=(m1+1)*ynumy+m5max+1;
ycell=(gy[m5max]+gy[m5max+1])/2.0;
{printf("%ld %d   %ld %ld %ld   %e %e    %e %e %e %e %e %e ",mm1,mm3,m1,m2,m5max,markdi[mm1],divmax,nd[m4],ycell0-ycell,pr[m3],GYKOEF*markro[mm2]*(ycell0-ycell),pr[m4],pr[m3]-GYKOEF*markro[mm2]*(ycell0-ycell)-pr[m4]);getchar();}
*/
			/* Add melt to the depth */
			if(divmax>0)
				{
				xvolkcur=xvolk;
				if(xvolkcur>1.0) xvolkcur=1.0;
				if(xvolkcur<0) xvolkcur=0;
				/* Check divergence */
				if(divmax*markvol<-markdi[mm1]*(1.0-xvolkcur)) 
					{
					xvolkcur=1.0+divmax*markvol/markdi[mm1];
					}
				/* New Marker Type */
				markt[marknum1]=mm3;
				/* New Plutonic Marker wt */
				markwt[marknum1]=1.0-(1.0-markwt[mm1])*(xmelt-markex[mm1]-minmelt)*(1.0-xvolkcur);
				/* New Marker coordinates */
				markx[marknum1]=markx[mm1];
				marky[marknum1]=gy[m5max]+(gy[m5max+1]-gy[m5max])*(marky[mm1]-gy[m2])/(gy[m2+1]-gy[m2]);
				/* New Marker divergence */
/*
				markdi[marknum1]=-markdi[mm1]*(1.0-markwt[mm1])*(1.0-xvolkcur)/(1.0-markwt[marknum1]);
*/
				markdi[marknum1]=-markdi[mm1]*(1.0-xvolkcur);
				/**/
				/* New marker temperature */
				allintert((double)(markx[marknum1]),(double)(marky[marknum1]));
				if(xtkmelt>=0)
					{
					markk[marknum1]=eps[2]*xtkmelt+(1.0-xtkmelt)*(markk[mm1]-(markx[mm1]-markx[marknum1])*0.5e-3);
					}
				else
					{
					/* Search for highest T in the column */
					for (m5=m5max+1;m5<=m2;m5++)
				  		{
						/* Current Temperature */
						m4=m1*ynumy+m5;
						tkcur=tk[m4]*(1-e)+tk[m4+ynumy]*e;
/*
{printf("%ld %d  %ld %ld %e     %e %e %e",marknum1,markt[marknum1],m1,m5,e,eps[2],tkcur,markk[mm1]);getchar();}
*/
						if(eps[2]<tkcur) eps[2]=tkcur;
						}
					if(eps[2]<markk[mm1]) eps[2]=markk[mm1];
					markk[marknum1]=eps[2]*(-xtkmelt)+(1.0+xtkmelt)*markk[mm1]-(markx[mm1]-markx[marknum1])*0.5e-3;
					}
				markdt[marknum1]=markk[marknum1];
/*
{printf("%ld %d  %e %e %e   %e %e  %e %e",marknum1,markt[marknum1],markk[marknum1],xtkmelt,xvolkcur,markdi[marknum1],markdi[mm1],1.0-markwt[marknum1],1.0-markwt[mm1]);getchar();}
*/
				/* Extraction, Density, Water, Reset */
				markn[marknum1]=0;
				markd[marknum1]=-1.0;
				markw[marknum1]=-1.0;
				markex[marknum1]=0;
				marke[marknum1]=0;
				markv[marknum1]=0;
				/**/
				/* Add new marker counter */
				marknum1++;
				}
			}
		/**/
		/**/
		/**/
		/* Add melt to the surface */
		if(xvolk>0)
			{
			/* New Marker Type */
			markt[marknum1]=mm3;
			/* New Volcanic Marker wt */
			markwt[marknum1]=1.0-(1.0-markwt[mm1])*(xmelt-markex[mm1]-minmelt)*xvolkcur;
/*
{printf("%ld %d %e %e",marknum1,markt[marknum1],markwt[marknum1],xvolkcur);getchar();}
*/
			/* New Marker coordinates */
			markx[marknum1]=markx[mm1];
			/* Surface level elevation for marker definition */
			marky[marknum1]=yeros+mardy*0.25;
			/* Change elevation */
			dx1=0;if(m1>0) dx1=gx[m1]-gx[m1-1];
			dx2=gx[m1+1]-gx[m1];
			dx3=0;if(m1<xnumx-2) dx3=gx[m1+2]-gx[m1+1];
			dy=2.0*markvol*(1.0-markwt[marknum1])/((dx1+dx2)*(1.0-e)+(dx2+dx3)*e);
			ep[m1]-=dy*(1.0-e);
			if(ep[m1]<0) ep[m1]=0;
			ep[m1+1]-=dy*e;
			if(ep[m1+1]<0) ep[m1+1]=0;
			/* New marker temperature */
			allintert((double)(markx[marknum1]),(double)(marky[marknum1]));
/*
if(markk[marknum1]>274){printf("%ld %d %e %e",marknum1,markt[marknum1],markk[marknum1],xtkmelt);getchar();}
*/
			/* Extraction, Density, Water, Reset */
			markk[marknum1]=eps[2];
			markdt[marknum1]=markk[marknum1];
			markn[marknum1]=0;
			markd[marknum1]=-1.0;
			markw[marknum1]=-1.0;
			markdi[marknum1]=0;
			markex[marknum1]=0;
			marke[marknum1]=0;
			markv[marknum1]=0;
			/**/
			/* Add new marker counter */
			marknum1++;
			}
		/**/
		/* Mark extraction column */
		if(markx[mm1]-dxwater<ep0[xnumx  +m1]) ep0[xnumx+  m1]=markx[mm1]-dxwater;
		if(markx[mm1]+dxwater>ep0[xnumx*2+m1]) ep0[xnumx*2+m1]=markx[mm1]+dxwater;
		if(marky[mm1]+dywater>ep0[xnumx*3+m1]) ep0[xnumx*3+m1]=marky[mm1]+dywater;
		/**/
		}
/*
printf("B %ld %d %e %e %e %e",mm1,mm2,markx[mm1],marky[mm1],xmelt,markex[mm1]);getchar();
*/
	/* Old Marker water, extraction level */
	xsolidcur=(1.0-xmelt)/(1.0-markex[mm1]);
	xmeltcur=(xmelt-markex[mm1])/(1.0-markex[mm1]);
	KDmeltsolid=100.0;
	H2Osolid=markw[mm1]/(xsolidcur+KDmeltsolid*xmeltcur);
	H2Omelt=H2Osolid*KDmeltsolid;
	markw[mm1]=(H2Osolid*(1.0-xmelt)+H2Omelt*minmelt)/(1.0-xmelt+minmelt);
/*
	markw[mm1]*=minmelt/(xmelt-markex[mm1]);
*/
/*
	if(markx[mm1]>1e+6 || markt[mm1]==29 || markt[mm1]==30) markex[mm1]=xmelt-minmelt;
*/
	markex[mm1]=xmelt-minmelt;
	/**/
/*
printf("C %ld %d %e %e %e %e",mm1,mm2,markx[mm1],marky[mm1],xmelt,markex[mm1]);getchar();
*/
	}
}
if(printmod) printf("\n MELT EXTRACTION END Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);
/**/
/**/
/**/
/* Remove air/water from volcanics */
for (mm1=0;mm1<marknum;mm1++)
if(markt[mm1]<2 && markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize)
	{
	/* Up Left Node X,Y Num */
	m1=m1serch(markx[mm1]);
	/* Relativ Normalized coord Calc */
	e=((double)(markx[mm1])-gx[m1])/(gx[m1+1]-gx[m1]);
	/* Initial & Final Surface level elevation for marker definition */
	y0=(e*ep0[m1+1]+(1.0-e)*ep0[m1]);
	y1=(e*ep[m1+1]+(1.0-e)*ep[m1]);
	/* Erase water/air marker */
	if(marky[mm1]>y1 && marky[mm1]<y0) 
		{
/*
printf(" B %ld %e %e   %ld %e %e ",mm1,markx[mm1],marky[mm1],m1,y0,y1);getchar();
*/
		markx[mm1]=-xsize;
		markk[mm1]=0;
		}
	}
/**/
/**/
/**/
/* Melt mark reset */
for (mm1=0;mm1<marknum1;mm1++)
{
if(markt[mm1]>1 && markt[mm1]<50 &&  markx[mm1]>0 && marky[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize)
	{
	/* Up Left Node X,Y Num */
	m1=m1serch(markx[mm1]);
/*
if(m1>0){printf(" A %ld %e %e   %ld %e %e %e %e %e",mm1,markx[mm1],marky[mm1],m1,gx[m1],gx[m1+1],ep0[xnumx+m1],ep0[xnumx*2],ep0[xnumx*3+m1]);getchar();}
*/
	/* Markers in melt columns */
	if(marky[mm1]<ep0[xnumx*3+m1] && markx[mm1]>ep0[xnumx+m1] && markx[mm1]<ep0[xnumx*2+m1]) 
		{
		markv[mm1]+=2;
/*
printf(" B %ld %e %e   %ld %e %e %e %e %e",mm1,markx[mm1],marky[mm1],m1,gx[m1],gx[m1+1],ep0[xnumx+m1],ep0[xnumx*2],ep0[xnumx*3+m1]);getchar();
*/
		} 
	}
}
/**/
/**/
/**/
if(printmod) printf("\n MELT EXTRACTION END Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);
marknum=marknum1;
/**/
}
/* End Melt extraction */



/* Water calculation */
void watercalc(double mtk, double mpb, int mm2, long int mm1)
{
double mwa=0;
long int m1=wn[0];
double sy1,e1;
if(marks0[mm2]>0) 
	{
	/* Erosion surface */
	e1=(markx[mm1]-gx[m1])/(gx[m1+1]-gx[m1]);
	sy1=marky[mm1]-(e1*ep[m1+1]+(1.0-e1)*ep[m1]);
	/* Gradual water release */
	if(mm2<20 && mm2!=13)
		{
		if(sy1>0 && sy1<zmpor && mtk<tkpor) mwa=marks0[mm2]*(tkpor-mtk)/(tkpor-273.15)*(zmpor-sy1)/zmpor;
		}
	/* Abrupt water release */
	else
		{
		mwa=marks0[mm2];
		}
	}
eps[42]=mwa;
}
/* Water calculation */


/* Thermodynamic database use for ro, Cp */
void tdbasecalc(double mtk, double mpb, int mm2, long int mm1)
{
/* TD Database variables,  dTK,dPB - TK, PB step for tabulation in TD database */
double H0,H1,H2,H3,R0,R1,R2,R3,G0,G1,G2,G3,W0,W1,W2,W3,n,e;
/* Val Buffers */
int n1,n2,mm3,ynpb;
double mhh0,mhh1,mdhh,maa,mwa,dmwa,wro,mro,mcp,mab,mbb,mgg,mkt,mkt1,pbmax,xold,kr01,kr1,kr10,xkr,krad;
long int m1=wn[0];
double sy1,e1;
/**/
/* Maximal pressure for the shallow database */
pbmax=pbmin+pbstp*(double)(pbnum-1);
/* Adiabate computing */
ynpb=0;
/*
if(0==0 && timesum<3.15576e+7*1e+3) {mpb*=timesum/(3.15576e+7*1e+3); ynpb=1;}
*/
/**/
/* Reset TD variables */
eps[40]=eps[41]=eps[42]=eps[43]=eps[44]=eps[20]=eps[45]=0;
/**/
/* Thermal conductivity */
krad=0;
if(1==0)
{
/* m895 Dry peridotite Fe=12 */
/* Olivine: Hoffmeister, 1999; Hoffmeister & Yuen, 2005 */
if(mpb<235000.0)
	{
	/* Lattice k */
	mkt1=(1.878+770.9/MINV(mtk,1200.0))*(1.0+4.26e-6*mpb);
	/* Radiative k 0.1 mm */
	kr01=pow(mtk/4000.0,3.0);
	/* Radiative k 1 mm */
	kr1=pow(mtk/1774.0,3.0);
	/* Radiative k 10 mm */
	xkr=pow(mtk/1636.0,10.0);
	xkr/=xkr+1.0; kr10=pow((mtk-1000.0*xkr)/1011.0,3.0)-0.7713*xkr;
	}
/* Perovskite: Hoffmeister, 1999; Hoffmeister & Yuen, 2005 */
else
	{
	/* Lattice k */
	mkt1=(1.291+1157.0/MINV(mtk,2100.0))*(1.0+2.50e-6*mpb);
	/* Radiative k 0.1 mm */
	kr01=pow(mtk/3591.0,3.0);
	/* Radiative k 1 mm */
	kr1=pow(mtk/2117.0,3.0);
	/* Radiative k 10 mm */
	xkr=pow(mtk/1500.0,4.0); xkr/=xkr+1.0;
	kr10=pow((mtk+4000.0*xkr)/5776.0,3.0)+2.822*xkr;
	}
krad=kr1;
}
/**/
/* Shallow TD base type */
if(mpb<pbmax && ynpb==0)
/*
*/
{
/* TD base type */
switch (mm2)
	{
	/* Sediments */
	case 15:
	case 2:
	case 3:
	case 4: mm3=5; break;
	/* Molten Sediments, Molten continental Crust */
	case 35:
	case 22:
	case 23:
	case 24: mm3=6; break;
	/* Basalt */
	case 16:
	case 17:
	case 18:
	case 7: mm3=7; break;
	/* Molten Basalt */
	case 36:
	case 37:
	case 38:
	case 27: mm3=8; break;
	/* Gabbro */
	case 8: mm3=3; break;
	/* Molten Gabbro */
	case 28: mm3=4; break;
	/* Dry peridotite */
	case 9:
	case 12:
	case 14:
	case 10: mm3=0; break;
	/* Wet peridotite */
	case 13:
	case 11: mm3=1; break;
	/* Molten Dry peridotite */
	case 29:
	case 30:
	case 32:
	/* Molten Wet peridotite */
	case 34: mm3=2; break;
	/* Unknown type */
	default: {printf("Unknown rock type for TD database %d",mm2); exit(0);}
	}
/* ABCD-4Cell Number */
e=(mtk-tkmin)/tkstp;
if(e<0) e=0;
if(e>(double)(tknum-1)) e=(double)(tknum-1);
n=(mpb-pbmin)/pbstp;
if(n<0) n=0;
if(n>(double)(pbnum-1)) n=(double)(pbnum-1);
n1=(int)(e);
if(n1>tknum-2) n1=tknum-2;
n2=(int)(n);
if(n2>pbnum-2) n2=pbnum-2;
/* e,n Calc */
e=(e-(double)(n1));
n=(n-(double)(n2));
/* Ro H values */
/* 0 2 */
/* 1 3 */
R0=td[n1  ][n2  ][mm3][0]*1000.0;
R1=td[n1  ][n2+1][mm3][0]*1000.0;
R2=td[n1+1][n2  ][mm3][0]*1000.0;
R3=td[n1+1][n2+1][mm3][0]*1000.0;
H0=td[n1  ][n2  ][mm3][1]*1000.0*4.1837;
H1=td[n1  ][n2+1][mm3][1]*1000.0*4.1837;
H2=td[n1+1][n2  ][mm3][1]*1000.0*4.1837;
H3=td[n1+1][n2+1][mm3][1]*1000.0*4.1837;
W0=td[n1  ][n2  ][mm3][4];
W1=td[n1  ][n2+1][mm3][4];
W2=td[n1+1][n2  ][mm3][4];
W3=td[n1+1][n2+1][mm3][4];
G0=td[n1  ][n2  ][mm3][3]*1000.0;G0*=G0*R0;
G1=td[n1  ][n2+1][mm3][3]*1000.0;G1*=G1*R1;
G2=td[n1+1][n2  ][mm3][3]*1000.0;G2*=G2*R2;
G3=td[n1+1][n2+1][mm3][3]*1000.0;G3*=G3*R3;
/* Shear modulus calc by interpolation */
mgg=((G0*(1.0-n)+G1*n)*(1.0-e)+(G2*(1.0-n)+G3*n)*e);
/* Ro calc by interpolation */
mro=((R0*(1.0-n)+R1*n)*(1.0-e)+(R2*(1.0-n)+R3*n)*e);
/* Water wt% calc by interpolation */
mwa=((W0*(1.0-n)+W1*n)*(1.0-e)+(W2*(1.0-n)+W3*n)*e);
/* Add porocity fluid */
/* Erosion surface */
e1=(markx[mm1]-gx[m1])/(gx[m1+1]-gx[m1]);
sy1=marky[mm1]-(e1*ep[m1+1]+(1.0-e1)*ep[m1]);
if(marks0[mm2]>0 && sy1>0 && sy1<zmpor && mtk<tkpor) 
	{
	dmwa=marks0[mm2]*(tkpor-mtk)/(tkpor-273.15)*(zmpor-sy1)/zmpor;
	mwa+=dmwa;
	wro=1050.0;
	mro=mro/(1.0+dmwa*1e-2*(mro/wro-1.0));
/*
if(marks0[mm2]>0){printf("TD1 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,sy1,mwa,mro,dmwa,zmpor);getchar();}
if(sy1>10000.0){printf("TD1 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,sy1,mwa,mro,dmwa,zmpor);getchar();}
if(sy1>10000.0){printf("TD2 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,sy1,mwa,mro,dmwa,zmpor);getchar();}
*/
	}
/* Limit amount of water in the molten mantle */
if(mm2==29 || mm2==30 || mm2==32) mwa=0;
if(mm2==34 && mwa>maxwater) mwa=maxwater;
/* Add porous water to the hydrated mantle in the region of fluid-fluxed melting: melting model of Katz et al. (2003) */
if((mm2==11 || mm2==13) && mwa<maxwater) mwa=maxwater;
/*
if((mm2==11 || mm2==13) && mwa<maxwater && mtk>(273.15+1085.7+132.9e-4*mpb-5.1e-8*mpb*mpb-43.0*pow(maxwater,0.75))) mwa=maxwater;
*/
/* Cp calc by interpolation */
mcp=((H2-H0)*(1.0-n)+(H3-H1)*n)/tkstp;
if(mcp<1e+2) mcp=1e+2; else if(mcp>5e+4) mcp=5e+4;
/* Effective adiabatic betta=1/V*dV/dT=ro/T*[-dH/dP+V] calc by interpolation */
mbb=(2.0/(R1+R0)-(H1-H0)/pbstp/1e+5)*(1.0-e)+(2.0/(R3+R2)-(H3-H2)/pbstp/1e+5)*e;
mbb*=mro/mtk;
if(mbb<-1e-2) mbb=-1e-2; else if(mbb>1e-2) mbb=1e-2;
/* Effective compressibility term alpha=1/ro*d(ro)/dP calc by interpolation */
maa=(2.0/(R1+R0)*(R1-R0)*(1.0-e)+2.0/(R3+R2)*(R3-R2)*e)/pbstp/1e+5;
if(maa<0) maa=0; else if(maa>1e-8) maa=1e-8;
/* Effective thermal expansion ab=-1/ro*d(ro)/dT calc by interpolation */
mab=-(2.0/(R2+R0)*(R2-R0)*(1.0-n)+2.0/(R3+R1)*(R3-R1)*n)/tkstp;
if(mab<-1e-2) mab=-1e-2; else if(mab>1e-2) mab=1e-2;
/* Activation enthalpy recalc using enthalpy changes */
/* Current Enthalpy */
mhh1=((H0*(1.0-n)+H1*n)*(1.0-e)+(H2*(1.0-n)+H3*n)*e);
/* Pmin Enthalpy */
mhh0=(td[n1][0 ][mm3][1]*(1.0-e) + td[n1+1][0 ][mm3][1]*e)*1000.0*4.1837;
/* Enthalpy Difference calc */
mdhh=(mhh1-mhh0);
/*
{printf("TD1 %d %d %e %e   %e %e %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mgg,mro,mwa,mcp,mbb);getchar();}
{printf("TD1 %d %d %e %e   %e %e %e %e %e %e   %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mgg,mro,mwa,mcp,mbb,maa,mdhh,mhh1,mhh0);getchar();}
eps[47]=mhh1;
eps[48]=mhh0;
*/
/* Save TD variables */
eps[40]=mgg;
eps[41]=mro;
eps[42]=mwa;
eps[43]=mcp;
eps[20]=mbb;
eps[44]=mab;
eps[45]=maa;
eps[46]=mdhh;
eps[47]+=krad;
}
/**/
/*
return;
*/
/**/
/* Deep TD base type */
if(1==0 || mpb>0.75*pbmax || ynpb==1)
{
switch (mm2)
	{
	/* MORB DATABASE */
	/* Sediments */
	case 2:
	case 3:
	case 4:
	case 15:
	/* Molten Sediments */
	case 22:
	case 23:
	case 24:
	case 35:
	/* Basalt */
	case 16:
	case 17:
	case 18:
	case 7:
	/* Molten Basalt */
	case 36:
	case 37:
	case 38:
	case 27:
	/* Gabbro */
	case 8:
	/* Molten Gabbro */
	case 28: mm3=10; break;
	/**/
	/* PIROLITE DATABASE */
	/* Dry peridotite */
	case 9:
	case 12:
	case 14:
	case 10:
	/* Wet peridotite */
	case 13:
	case 11:
	/* Molten Dry peridotite */
	case 29:
	case 30:
	case 32:
	/* Molten Wet peridotite */
	case 34: mm3=9; break;
	/* Unknown type */
	default: {printf("Unknown rock type for TD database %d",mm2); exit(0);}
	}
/* ABCD-4Cell Number */
e=(mtk-tkmin1)/tkstp1;
if(e<0) e=0;
if(e>(double)(tknum1-1)) e=(double)(tknum1-1);
n=(mpb-pbmin1)/pbstp1;
if(n<0) n=0;
if(n>(double)(pbnum1-1)) n=(double)(pbnum1-1);
n1=(int)(e);
if(n1>tknum1-2) n1=tknum1-2;
n2=(int)(n);
if(n2>pbnum1-2) n2=pbnum1-2;
/* e,n Calc */
e=(e-(double)(n1));
n=(n-(double)(n2));
/* Ro H values */
/* 0 2 */
/* 1 3 */
R0=td[n1  ][n2  ][mm3][0]*1000.0;
R1=td[n1  ][n2+1][mm3][0]*1000.0;
R2=td[n1+1][n2  ][mm3][0]*1000.0;
R3=td[n1+1][n2+1][mm3][0]*1000.0;
H0=td[n1  ][n2  ][mm3][1]*1000.0*4.1837;
H1=td[n1  ][n2+1][mm3][1]*1000.0*4.1837;
H2=td[n1+1][n2  ][mm3][1]*1000.0*4.1837;
H3=td[n1+1][n2+1][mm3][1]*1000.0*4.1837;
W0=td[n1  ][n2  ][mm3][4];
W1=td[n1  ][n2+1][mm3][4];
W2=td[n1+1][n2  ][mm3][4];
W3=td[n1+1][n2+1][mm3][4];
G0=td[n1  ][n2  ][mm3][3]*1000.0;G0*=G0*R0;
G1=td[n1  ][n2+1][mm3][3]*1000.0;G1*=G1*R1;
G2=td[n1+1][n2  ][mm3][3]*1000.0;G2*=G2*R2;
G3=td[n1+1][n2+1][mm3][3]*1000.0;G3*=G3*R3;
/* Shear modulus calc by interpolation */
mgg=((G0*(1.0-n)+G1*n)*(1.0-e)+(G2*(1.0-n)+G3*n)*e);
/* Ro calc by interpolation */
mro=((R0*(1.0-n)+R1*n)*(1.0-e)+(R2*(1.0-n)+R3*n)*e);
/* Water wt% calc by interpolation */
mwa=0;
/* Water in crystals */
if(mm2!=9 && mm2!=10 && mm2!=14 && mpb<235000.0) 
	{
	dmwa=0.1;
	mwa+=dmwa;
	wro=1050.0;
	mro=100.0/((100.0-dmwa)/mro+dmwa/wro);
/*
{printf("TD1 %d %d %e %e   %e %e %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mgg,mro,mwa,mcp,mbb);getchar();}
*/
	}
/* Cp calc by interpolation */
mcp=((H2-H0)*(1.0-n)+(H3-H1)*n)/tkstp1;
if(mcp<1e+2) mcp=1e+2; else if(mcp>5e+4) mcp=5e+4;
/* Effective adiabatic betta=1/V*dV/dT=ro/T*[-dH/dP+V] calc by interpolation */
mbb=(2.0/(R1+R0)-(H1-H0)/pbstp1/1e+5)*(1.0-e)+(2.0/(R3+R2)-(H3-H2)/pbstp1/1e+5)*e;
mbb*=mro/mtk;
if(mbb<-1e-2) mbb=-1e-2; else if(mbb>1e-2) mbb=1e-2;
/* Effective compressibility term alpha=1/ro*d(ro)/dP calc by interpolation */
maa=(2.0/(R1+R0)*(R1-R0)*(1.0-e)+2.0/(R3+R2)*(R3-R2)*e)/pbstp1/1e+5;
if(maa<0) maa=0; else if(maa>1e-8) maa=1e-8;
/* Effective thermal expansion ab=-1/ro*d(ro)/dT calc by interpolation */
mab=-(2.0/(R2+R0)*(R2-R0)*(1.0-n)+2.0/(R3+R1)*(R3-R1)*n)/tkstp;
if(mab<-1e-2) mab=-1e-2; else if(mab>1e-2) mab=1e-2;
/* Activation enthalpy recalc using enthalpy changes */
/* Current Enthalpy */
mhh1=((H0*(1.0-n)+H1*n)*(1.0-e)+(H2*(1.0-n)+H3*n)*e);
/* Pmin Enthalpy */
mhh0=(td[n1][0 ][mm3][1]*(1.0-e) + td[n1+1][0 ][mm3][1]*e)*1000.0*4.1837;
/* Enthalpy Difference calc */
mdhh=(mhh1-mhh0);
/* Thermal conductivity */
mkt=mkt1+krad;
/**/
/*
if(timesum>0){printf("TD4  %d %d   %e %e     %e   %e",mm2,mm3,mtk-273.15,mpb/1000.0,mkt,pbmax);getchar();}
if(timesum>0){printf("TD4 %d %d %e %e %e  %e %e %e %e %e %e %e   %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,pbmax,mgg,mro,mwa,mcp,mbb,maa,mkt,mdhh,mhh1,mhh0);getchar();}
*/
/* Computing transitional parameters */
if(1==0 || mpb>pbmax || ynpb==1)
	{
	/* Save TD variables */
	eps[40]=mgg;
	eps[41]=mro;
	eps[42]=mwa;
	eps[43]=mcp;
	eps[20]=mbb;
	eps[44]=mab;
	eps[45]=maa;
	eps[46]=mdhh;
	eps[47]=mkt;
/*
if(timesum>0){printf("TD2 %d %d %e %e %e  %e %e %e %e %e %e %e   %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,pbmax,mgg,mro,mwa,mcp,mbb,maa,mkt,mdhh,mhh1,mhh0);getchar();}
*/
	}
else
	{
	xold=(pbmax-mpb)/(0.25*pbmax);
	/* Save TD variables */
	mgg=mgg*(1.0-xold)+eps[40]*xold;
	mro=mro*(1.0-xold)+eps[41]*xold;
	mwa=mwa*(1.0-xold)+eps[42]*xold;
	mcp=mcp*(1.0-xold)+eps[43]*xold;
	mbb=mbb*(1.0-xold)+eps[20]*xold;
	mab=mab*(1.0-xold)+eps[44]*xold;
	maa=maa*(1.0-xold)+eps[45]*xold;
	mdhh=mdhh*(1.0-xold)+eps[46]*xold;
	mkt=mkt*(1.0-xold)+eps[47]*xold;
	eps[40]=mgg;
	eps[41]=mro;
	eps[42]=mwa;
	eps[43]=mcp;
	eps[20]=mbb;
	eps[44]=mab;
	eps[45]=maa;
	eps[46]=mdhh;
	eps[47]=mkt;
/*
if(timesum>0){printf("TD3 %d %d %e %e %e  %e %e %e %e %e %e %e   %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,pbmax,mgg,mro,mwa,mcp,mbb,maa,mkt,mdhh,mhh1,mhh0);getchar();}
*/
	}
}
/*
{printf("TD2 %d %d %e %e   %e %e %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mgg,mro,mwa,mcp,mbb);getchar();}
*/
}
/* Thermodynamic database use for ro, Cp */




/* Grid spacing change */
void gridchange()
{
/* Counters */
long int m1,m2,m3;
double e,ival;
/**/
/**/
/**/
/* 
printf("A %e %e %e",gx[46],gx[246],gx[292]);getchar();
*/
/* Save old erosion/sedimentation, hydration surfaces */
for (m1=0;m1<xnumx;m1++)
	{
	ep0[m1]=ep[m1];
	ep0[xnumx+m1]=ep[xnumx+m1];
	ep0[xnumx*2+m1]=ep[xnumx*2+m1];
	ep0[xnumx*3+m1]=ep[xnumx*3+m1];
	}
/* Recomputing horizontal grid */
ival=0;
for (m3=0;m3<=marknum;m3++) 
/* Check markers out of grid */
if (markx[m3]>0 && marky[m3]>0 && (double)(markx[m3])<xsize && (double)(marky[m3])<ysize && markt[m3]<50)
	{
	/* Check Asthenosphere position */
	if (markt[m3]==10 && markx[m3]>gx[100] && markx[m3]<gx[400] && marky[m3]<9e+4 && markx[m3]>ival) ival=markx[m3];
	}
if (printmod) printf("Trench position = %e km      step1 = %e km    step2 = %e km \n",ival/1000.0,(ival-3e+5)/100.0/1000.0,(3e+6-ival-3e+5)/110.0/1000.0);
ep0[xnumx*4]=0;
for (m1=1;m1<=100;m1++) ep0[xnumx*4+m1]=ep0[xnumx*4+m1-1]+2e+3+(ival-3e+5-2e+3*100.0)/(100.0*101.0/2.0)*(101.0-(double)(m1));
for (m1=101;m1<=400;m1++) ep0[xnumx*4+m1]=ep0[xnumx*4+m1-1]+2e+3;
for (m1=401;m1<=510;m1++) ep0[xnumx*4+m1]=ep0[xnumx*4+m1-1]+2e+3+(3e+6-ival-3e+5-2e+3*110.0)/(110.0*111.0/2.0)*((double)(m1)-400.0);
/* 
ep0[xnumx*4+510]=gx[xnumx-1];
printf("%e %e %e",ep0[xnumx*4+100],ep0[xnumx*4+400],ep0[xnumx*4+510]);getchar();
*/
/* Reinterpolate erosion/sedimentation, hydration surfaces */
for (m3=0;m3<xnumx;m3++)
	{
	m1=m1serch(ep0[xnumx*4+m3]);
	/* Relativ Normalized coord Calc */
	e=(ep0[xnumx*4+m3]-gx[m1])/(gx[m1+1]-gx[m1]);
	/* Surface level elevation for marker definition */
	ep[m3]=(e*ep0[m1+1]+(1.0-e)*ep0[m1]);
	ep[xnumx+m3]=(e*ep0[xnumx+m1+1]+(1.0-e)*ep0[xnumx+m1]);
	ep[xnumx*2+m3]=(e*ep0[xnumx*2+m1+1]+(1.0-e)*ep0[xnumx*2+m1]);
	ep[xnumx*3+m3]=(e*ep0[xnumx*3+m1+1]+(1.0-e)*ep0[xnumx*3+m1]);
	}
/* Reload new gridline positions */
for (m1=0;m1<xnumx;m1++)
	{
	gx[m1]=ep0[xnumx*4+m1];
	}
}
/* End Grid spacing change */





/* Errosion Surface progress */
void erosion()
{
/* Val buffer */
double v0,v1,dydx,x1,vx1,vy1,dy;
double ertimesum,ertimesum0;
long int m1,m2;
int ynrelax,n1;
/**/
/* Erosion Solution Cycle ------------------------------------------ */
ertimesum=0;
ertimesum0=timestep;
do
{
/* Save old cycle results */
for (m1=0;m1<xnumx;m1++)
	{
	ep0[m1]=ep[m1];
	}
/**/
/**/
/**/
/* Initial timestep definition */
timestep=ertimesum0-ertimesum;
/**/
/**/
/**/
/* Erosion timestep definition using material velosity field */
for (m1=0;m1<xnumx;m1++)
	{
	/* Calc horisontal Coordinate */
	x1=gx[m1];
	/**/
	/* EROSION SURFACE */
	/* Calc material velocity on the Surface using velosity field */
	allinterv(x1,ep0[m1]);
	vx1=eps[11];
	vy1=eps[12];
	/* Check horizontal timestep */
	/* Calc x derivative of y position of the Surface using upwind differences */
	dydx=0;
	if(vx1>0 && m1>0)
		{
		timestep=MINV(timestep,(gx[m1]-gx[m1-1])/vx1);
/*
printf("111 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1-1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	if(vx1<0 && m1<xnumx-1)
		{
		timestep=MINV(timestep,(gx[m1]-gx[m1+1])/vx1);
/*
printf("222 %ld %e %e %e %e %e %e %e",m1,vx1,vy1,(gx[m1]-gx[m1+1])/vx1,ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	/* Check vertical timestep */
	if(vy1)
		{
		/* Horizontal line num definition */
		m2=m2serch(ep0[m1]);
		/* Check timestep */
		timestep=MINV(timestep,(gy[m2+1]-gy[m2])/ABSV(vy1));
/*
printf("333 %ld  %e %e %e %e %e %e %e",m2,vx1,vy1,(gy[m2+1]-gy[m2])/ABSV(vy1),ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
		}
	}
/*
printf("777 %e %e %e %e",ertimesum0,ertimesum,ertimesum0-ertimesum,timestep);getchar();
*/
/**/
/**/
/**/
/* Displace Surface boundary */
/*
for (m1=1;m1<xnumx-1;m1++)
*/
for (m1=0;m1<xnumx;m1++)
	{
	/* EROSION SURFACE */
	/* Calculation of errosion rate */
	v0=0;
	if(ep0[m1]<eroslev)
		{
		v0=eroscon+eroskoe*(eroslev-ep0[m1]);
		}
	/* Calculation of sedimentation rate */
	v1=0;
	if(ep0[m1]>sedilev)
		{
		v1=sedicon+sedikoe*(ep0[m1]-sedilev);
		}
	/* Calc horisontal Coordinate */
	x1=gx[m1];
	/**/
	/* Calc material velocity on the Surface using velosity field */
	allinterv(x1,ep0[m1]);
	vx1=eps[11];
	vy1=eps[12];
	/**/
	/* Erase erosion/sedimentation rate for marginal points */
	if((m1==0 && vx1>0) || (m1==xnumx-1 && vx1<0)) v0=v1=0;
	/**/
	/* Calc x derivative of y position of the Surface using upwind differences */
	dydx=0;
	if(vx1>0 && m1>0)
		{
		dydx=(ep0[m1]-ep0[m1-1])/(gx[m1]-gx[m1-1]);
/*
printf("AAA %e %e",ep0[m1],dydx);getchar();
*/
		}
	if(vx1<0 && m1<xnumx-1)
		{
		dydx=(ep0[m1+1]-ep0[m1])/(gx[m1+1]-gx[m1]);
/*
printf("BBB %e %e",ep0[m1],dydx);getchar();
*/
		}
	/* Recalc new Surface position */
	ep[m1]+=timestep*(v0-v1+vy1-dydx*vx1);
/*
printf("SURFACE %ld %e %e %e %e %e %e %e %e",m1,x1,v0,v1,vx1,vy1,dydx,ep[m1]);getchar();
*/
	/**/
	/**/
	/**/
	}
/**/
/**/
/**/
/**/
/* Relax EROSION surface */
n1=0;
do
	{
	/* Reset relaxation mark */
	ynrelax=0;
	for (m1=0;m1<xnumx-1;m1++)
		{
		/* Calc x derivative of y position */
		dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
		/* Relax surface for critical slope */
		if(dydx>slopemax)
			{
			dy=((ep[m1+1]-ep[m1])-slopemax*0.99*(gx[m1+1]-gx[m1]))/2.0;
			ep[m1]  +=dy;
			ep[m1+1]-=dy;
			ynrelax=1;
/*
printf("AAA %ld %e %e %e",m1,slopemax,dydx,(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]));getchar();
dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
*/
			}
		if(dydx<-slopemax)
			{
			dy=((ep[m1+1]-ep[m1])+slopemax*0.99*(gx[m1+1]-gx[m1]))/2.0;
			ep[m1]  +=dy;
			ep[m1+1]-=dy;
			ynrelax=1;
/*
printf("BBB %ld %e %e %e",m1,slopemax,dydx,(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]));getchar();
dydx=(ep[m1+1]-ep[m1])/(gx[m1+1]-gx[m1]);
*/
			}
		}
	n1++;
/*
printf("\n %d EROSION STEP = %e yr    EROSION TIME = %e yr  RELAXATION CYCLES = %d\n",ynrelax,timestep/3.15576e+7,ertimesum/3.15576e+7,n1); getchar();
*/
	}
while(ynrelax>0);
	/**/
	/**/
	/**/
/* Add Erosion step */
ertimesum+=timestep;
/**/
/**/
/**/
/* Print Results */
if (printmod)
	{
	printf("\n EROSION STEP = %e yr    EROSION TIME = %e yr  RELAXATION CYCLES = %d\n",timestep/3.15576e+7,ertimesum/3.15576e+7,n1);
/*
getchar();
*/
	}
}
while(ertimesum<ertimesum0);
/* Restore timestep */
timestep=ertimesum0;
}
/* Errosion Surface progress */


