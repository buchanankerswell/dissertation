/* ADD & SOLVE MATRIX BY ECONOMICAL FRONTAL GAUSS METHOD */
int gausmat3(int am, long int mcmax, long int mcmin)
/* wn[] - line koef numbers */
/* wi[] - line koef values */
/* am - mode of action rm=1,2,3 - add, rm=0 - solve */
/* val0[] - matrix contents */
/* lin0[] - line numbers for matrix contents */
/* pos0[] - first pos numbers for line in  val0[], lin0[] */
/* num0[] - pos numbers for line in  val0[], lin0[] */
/* fre0[] - free member for lines */
/* pos0cur - first free position number in val0[], lin0[] */
/* mcmax - current line in num0[] */
/* mcmin - first line in num0[] */
{
/* Counters */
long int m1,m2,m3,nempty;
/* Limits */
long int linbeg,linend;
/* Val Buffer */
double ival;
/**/
/* Space Check */
if (mcmax>=MAXPAR)
	{
	printf("EXIT PROGRAM: Space out in fre0[] %ld",mcmax);
	exit(0);
	}
/**/
/**/
/**/
/* STEP 1: EXEED KOEF ELIMINATION MATRIX BODY FORMATION am>0 */
if (am>0)
	{
/*
printf("\n GAUS %ld %ld   %ld   %ld %e ",am,mcmax,pos0cur,wn[0],wi[0]); getchar();
*/
	/* First, Last line numbers */
	linbeg=wn[1]; linend=wn[wn[0]];
	/**/
	/* Free member reload from buffer */
	fre0[mcmax]=wi[0];
	/**/
	/* Line koef reload from buffer */
	for (m2=1;m2<=wn[0];m2++)
		{
		bufv[wn[m2]]+=wi[m2];
		/**/
		/* First, Last line numbers */
		if(wi[m2])
			{
			if (linbeg>wn[m2]) linbeg=wn[m2];
			if (linend<wn[m2]) linend=wn[m2];
			}
/*
printf("%ld %ld %e",m2,wn[m2],wi[m2]); getchar();
*/
		}
	/**/
	/* Cur line koef recalc in buffer */
	for (m2=linbeg;m2<mcmax;m2++)
		{
		if(bufv[m2])
			{
			ival=bufv[m2];
			bufv[m2]=0;
			/* Check Presence of line */
			if (!num0[m2]) 
				{
				printf("EXIT PROGRAM: Line  %ld absent in matrix when processing Line %ld",m2,mcmax);
				exit(0);
				}
			/* Current Line koef recalc after cur koef and upper line */
			/* 1-st coef of any upper line = 1 */
			for (m3=1;m3<num0[m2];m3++)
				{
				bufv[lin0[pos0[m2]+m3]]-=ival*val0[pos0[m2]+m3];
				}
			/**/
			/* Free member recalc after cur koef upper line av=0,1,2 */
			fre0[mcmax]-=fre0[m2]*ival;
			/* Check last line number */
			linend=MAXV(linend,lin0[pos0[m2]+num0[m2]-1]);
			/**/
			}
		}
	/**/
	/* Cur line save in val0[],lin0[],num0[],pos0[] */
	/* Check Singularity */
	/* Check Presence of line */
	if (!bufv[mcmax])
		{
		printf("EXIT PROGRAM: Matrix is singular at Line %ld",mcmax);
		exit(0);
		}
	pos0[mcmax]=pos0cur;
	num0[mcmax]=0;
	ival=bufv[mcmax];
	for (m2=mcmax;m2<=linend;m2++)
		{
		/**/
		/* Recalc and Save Val>0 */
		if (bufv[m2])
			{
			/* Save Cur Koef */
			lin0[pos0cur]=m2;
			val0[pos0cur]=bufv[m2]/ival;
/*
printf("%ld %ld %ld %ld %e ",mcmax,pos0cur,m2,num0[pos0cur],val0[pos0cur]); getchar();
*/
			pos0cur++;
			/* Check Space */
			if (pos0cur>=MAXMAT) 
				{
				printf("EXIT PROGRAM: Space out in val0[] %ld %ld",mcmax,pos0cur);
				exit(0);
				}
			num0[mcmax]++;
			}
		/* Clear Cur Koef */
		bufv[m2]=0;
		}
	/* Free member recalc */
	fre0[mcmax]/=ival;
	return 0;
	}
/* End STEP 1: EXEED KOEF ELIMINATION MATRIX BODY FORMATION am>0 */
/**/
/**/
/**/
/* STEP 3: SOLUTION CALC CHECK am=0 */
nempty=0;
printf("TOTAL BAND MATRIX SIZE: val0[%ld]\n",pos0cur);
/*
*/
for (m1=mcmax;m1>=mcmin;m1--)
	{
	/* Calc sol0[] */
	if (num0[m1])
		{
		/* Recalc koef after Sol */
		ival=fre0[m1];
		for (m2=0;m2<num0[m1];m2++)
			{
			ival-=val0[pos0[m1]+m2]*sol0[lin0[pos0[m1]+m2]];
			}
		/* Calc Sol */
		sol0[m1]=ival;
		}
	else
		{
		sol0[m1]=0;
/*
printf("%ld %ld %e ",m1,num0[m1],sol0[m1]); getchar();
*/
		nempty++;
		}
/*
printf("%ld %ld %e ",m1,num0[m1],sol0[m1]); getchar();
*/
	}
/* End STEP 3: SOLUTION CALC CHECK am=0 */
/*
printf("%ld %ld %ld %e",nempty,pos0cur,mcmax,(double)(pos0cur+1)/(double)(mcmax+1)); getchar();
*/
return nempty;
}
/* End SOLVE MATRIX BY ECONOMICAL GAUSS METHOD */


