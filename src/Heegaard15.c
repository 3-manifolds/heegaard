#include "Heegaard.h"
#include "Heegaard_Dec.h"

int Do_Initialization(int open_Heegaard_Results)
{
	register unsigned int 	i;
	
	unsigned long 			Seconds;
	
#ifdef MAC
	MaxApplZone();
	InitGraf(&qd.thePort);
	FlushEvents(everyEvent,0);
	for(i = 0; i < 2000; i++) MoreMasters();
#endif
	
	for(i = 0; i < VERTICES; i++)
		{
		A[i] 		= (unsigned int *)	NewPtr(sizeof(int)*VERTICES);
		if(A[i] 	== NULL) return(1);
		AA[i] 		= (unsigned int *)	NewPtr(sizeof(int)*VERTICES);
		if(AA[i] 	== NULL) return(2);
		AJ1[i] 		= (unsigned int *)	NewPtr(sizeof(int)*VERTICES);
		if(AJ1[i] 	== NULL) return(3);
		AJ2[i] 		= (unsigned int *)	NewPtr(sizeof(int)*VERTICES);
		if(AJ2[i] 	== NULL) return(4);
		AJ3[i] 		= (unsigned int *)	NewPtr(sizeof(int)*VERTICES);
		if(AJ3[i] 	== NULL) return(5);
		B[i] 		= (unsigned int *)	NewPtr(sizeof(int)*VERTICES);
		if(B[i] 	== NULL) return(6);
		CO[i]		= (unsigned char*)	NewPtr(sizeof(char)*VERTICES);
		if(CO[i] 	== NULL) return(7);
		GB[i]		= (         int *)  NewPtr(sizeof(int)*VERTICES);
		if(GB[i] 	== NULL) return(8);
		}

	for(i = 0; i < (VERTICES)/2; i++)
		{
		ED[i]		= (unsigned int *)  NewPtr(sizeof(int)*VERTICES);
		if(ED[i] 	== NULL) return(9);
		EXL[i]		= (unsigned int *)  NewPtr(sizeof(int)*4);
		if(EXL[i] 	== NULL) return(10);
		EXR[i]      = (unsigned int *)  NewPtr(sizeof(int)*4);
		if(EXR[i] 	== NULL) return(11);
		EXP[i]		= (unsigned int *)  NewPtr(sizeof(int)*4);
		if(EXP[i] 	== NULL) return(12);
		NEX[i]		= (unsigned int *)  NewPtr(sizeof(int)*4);
		if(NEX[i] 	== NULL) return(13);
		T[i]        = (unsigned char *) NewPtr(sizeof(char)*8);
		if(T[i] 	== NULL) return(14);
		}
	
	for(i = 0; i < 2*MAXNUMRELATORS; i++)
		{
		DRA[i] = (unsigned int *) NewPtr(sizeof(int)*2*MAXNUMRELATORS);
		if(DRA[i] 	== NULL) return(15);
		}
				
	for(i = 0; i < MAXNUMCOMPONENTS; i++)
		{
		CBC[i] = (unsigned char *) NewPtr(sizeof(char)*(MAXNUMGENERATORS + 2));
		if(CBC[i] 	== NULL) return(16);
		MLC[i] = (unsigned long *) NewPtr(sizeof(long)*(MAXNUMGENERATORS + 1));
		if(MLC[i] 	== NULL) return(17);
		}
	
	for(i = 1; i < 2*VERTICES; i++)
		{
		Face[i] = (unsigned char *) NewPtr(sizeof(char)*(VERTICES + 1));
		if(Face[i] 	== NULL) return(18);
		}
			
	AT				= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(AT 			== NULL) return(19);
	BCF				= (unsigned char *)  NewPtr(sizeof(char)*2*VERTICES);
	if(BCF 			== NULL) return(20);
	BCWG 			= (unsigned char *)	 NewPtr(sizeof(char)*(MAXNUMGENERATORS + 1));
	if(BCWG 		== NULL) return(21);
	BDY				= (         int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(BDY 			== NULL) return(22);
	Bdry			= (unsigned	int  *)  NewPtr(sizeof(int)*(VERTICES + 2));
	if(Bdry 		== NULL) return(23);
	BSV1			= (unsigned int  *)  NewPtr(sizeof(int)*MAX_SAVED_LEVELS);
	if(BSV1 		== NULL) return(24);
	BSV2			= (unsigned int  *)  NewPtr(sizeof(int)*MAX_SAVED_LEVELS);
	if(BSV2 		== NULL) return(25);
	ComponentNum 	= (unsigned int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(ComponentNum == NULL) return(26);
	CS				= (unsigned char *)  NewPtr(sizeof(char)*(MAXNUMCOMPONENTS + 1));
	if(CS 			== NULL) return(27);
	Daughters 		= (unsigned int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(Daughters 	== NULL) return(28);
	DeletedEdges	= (unsigned char *)	 NewPtr(sizeof(char)*6*VERTICES);
	if(DeletedEdges == NULL) return(29);
	DF				= (unsigned int  *)  NewPtr(sizeof(int)*(VERTICES)/2);
	if(DF 			== NULL) return(30);
	InDisk			= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(InDisk 		== NULL) return(31);
	ER				= (			char *)  NewPtr(sizeof(char)*MAX_SAVED_PRES);
	if(ER 			== NULL) return(32);
	Father			= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(Father 		== NULL) return(33);
	Flags			= (         int  *)  NewPtr(sizeof(int)*VERTICES);
	if(Flags 		== NULL) return(34);
	FR 				= (unsigned int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(FR 			== NULL) return(35);
	FV 				= (unsigned char *)	 NewPtr(sizeof(char)*VERTICES);
	if(FV 			== NULL) return(36);
	GBC 			= (unsigned char *)	 NewPtr(sizeof(char)*(MAXNUMRELATORS + 2));
	if(GBC 			== NULL) return(37);
	GV2				= (unsigned int  *)  NewPtr(sizeof(int)*(VERTICES)/2);
	if(GV2 			== NULL) return(38);
	GV2L			= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(GV2L 		== NULL) return(39);
	GV2R			= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(GV2R 		== NULL) return(40);
	InDisk			= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(InDisk 		== NULL) return(41);	
	InPS			= (			int	 *)	 NewPtr(sizeof(int)*VERTICES);
	if(InPS 		== NULL) return(42);	
	InQueue 		= (unsigned int  *)	 NewPtr(sizeof(int)*VERTICES);
	if(InQueue 		== NULL) return(43);
	IV				= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(IV 			== NULL) return(44);
	Left			= (unsigned int  *)  NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(Left 		== NULL) return(45);
	Lowpt			= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(Lowpt 		== NULL) return(46);
	LR 				= (unsigned long *)  NewPtr(sizeof(long)*(MAXNUMRELATORS + 1));
	if(LR 			== NULL) return(47);
	LSP				= (unsigned long *)  NewPtr(sizeof(long)*(MAX_SAVED_PRES));
	if(LSP 			== NULL) return(48);
	LSQ				= (unsigned long *)	 NewPtr(sizeof(long)*(MAX_SAVED_PRES));
	if(LSQ 			== NULL) return(49);
	NCS				= (unsigned char *)	 NewPtr(sizeof(char)*(MAX_SAVED_PRES));
	if(NCS 			== NULL) return(50);					
	NEBC 			= (unsigned int  *)	 NewPtr(sizeof(int)*2*MAXNUMRELATORS);
	if(NEBC 		== NULL) return(51);
	N1H				= (unsigned char *)  NewPtr(sizeof(char)*MAXNUMCOMPONENTS);
	if(N1H 			== NULL) return(52);
	NFBC			= (unsigned int  *)	 NewPtr(sizeof(int)*2*MAXNUMRELATORS);
	if(NFBC 		== NULL) return(53);	
	NG 				= (         int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(NG 			== NULL) return(54);
	NR 				= (         int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(NR 			== NULL) return(55);
	NRBC 			= (unsigned int  *)	 NewPtr(sizeof(int)*2*MAXNUMRELATORS);
	if(NRBC 		== NULL) return(56);
	NS1XD2			= (unsigned char *)  NewPtr(sizeof(char)*MAXNUMCOMPONENTS);
	if(NS1XD2 		== NULL) return(57);
	NS1XS2  		= (unsigned char *)  NewPtr(sizeof(char)*MAXNUMCOMPONENTS);
	if(NS1XS2 		== NULL) return(58);
	Number			= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(Number 		== NULL) return(59);
	OSA				= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(OSA 			== NULL) return(60);
	OSB				= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(OSB 			== NULL) return(61);
	PG				= (unsigned int  *)  NewPtr(sizeof(int)*(VERTICES)/2);
	if(PG 			== NULL) return(62);
	PresName		= (unsigned char *)  NewPtr(1000L);
	if(PresName 	== NULL) return(63);
	PRIM			= (         int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(PRIM 		== NULL) return(64);
	QPM				= (unsigned char *)  NewPtr(sizeof(char)*(MAX_SAVED_PRES));
	if(QPM 			== NULL) return(65);	
	UpDate 			= (unsigned int  *)	 NewPtr(sizeof(int)*VERTICES*VERTICES);
	if(UpDate 		== NULL) return(66);
	Right			= (unsigned int  *)  NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(Right 		== NULL) return(67);
	SaveBdry		= (         int  *)  NewPtr(sizeof(int)*VERTICES);
	if(SaveBdry 	== NULL) return(68);
	SUR				= (unsigned char ****)	 NewPtr(sizeof(long)*MAX_SAVED_PRES);
	if(SUR 			== NULL) return(69);
	SURL			= (unsigned long *)	 NewPtr(sizeof(long)*MAX_SAVED_PRES);
	if(SURL 		== NULL) return(70);
	SV				= (unsigned int  *)  NewPtr(sizeof(int)*(VERTICES)/2);
	if(SV 			== NULL) return(71);
	TP				= (unsigned char *)  NewPtr(sizeof(char)*(MAX_SAVED_PRES));
	if(TP 			== NULL) return(72);
	TV				= (unsigned int  *)  NewPtr(sizeof(int)*(VERTICES)/2);
	if(TV 			== NULL) return(73);
	UDV 			= (unsigned int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(UDV			== NULL) return(74);
	V				= (unsigned int  *)	 NewPtr(sizeof(int)*VERTICES);
	if(V 			== NULL) return(75);
	VA 				= (unsigned int  *)	 NewPtr(sizeof(int)*(VERTICES)/2);
	if(VA 			== NULL) return(76);
	VWG				= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(VWG		 	== NULL) return(77);
	X				= (         int  *)	 NewPtr(sizeof(int)*VERTICES);
	if(X 			== NULL) return(78);
	XX				= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(XX 			== NULL) return(79);
	Y				= (         int  *)	 NewPtr(sizeof(int)*VERTICES);
	if(Y 			== NULL) return(80);
	YY				= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(YY 			== NULL) return(81);
	ZZ 				= (unsigned int  *)	 NewPtr(sizeof(int)*VERTICES);
	if(ZZ 			== NULL) return(82);	
	zz 				= (unsigned int  *)	 NewPtr(sizeof(int)*2*MAXNUMRELATORS);
	if(zz 			== NULL) return(83);
		
	for(i = 0; i < MAX_SAVED_PRES; i++)
		{
		SUR[i] = (unsigned char ***) NewPtr(sizeof(long)*(MAXNUMRELATORS + 1));
		if(SUR[i] 	== NULL) return(84);
		}
			
	for(i = 0; i < MAX_SAVED_LEVELS; i++)
		{
		SLR[i] = (unsigned char ***) NewPtr(sizeof(long)*(MAXNUMRELATORS + 1));
		if(SLR[i] 	== NULL) return(85);
		}
					
	for(i = 0; i <= MAXNUMRELATORS; i++)
		{ 
		Copy_Of_Input[i] 	= (unsigned char **) NewHandle(4L);
		if(Copy_Of_Input[i] == NULL) return(86);
		Copy_Of_Rel_1[i] 	= (unsigned char **) NewHandle(4L);
		if(Copy_Of_Rel_1[i] == NULL) return(87);
		Copy_Of_Rel_2[i] 	= (unsigned char **) NewHandle(4L);
		if(Copy_Of_Rel_2[i] == NULL) return(88);
		CD_Surgery_Rel[i]	= (unsigned char **) NewHandle(4L);
		if(CD_Surgery_Rel[i]== NULL) return(89);	
		DelRelators[i] 		= (unsigned char **) NewHandle(4L);
		if(DelRelators[i] 	== NULL) return(90);		
		DualRelators[i] 	= (unsigned char **) NewHandle(4L);
		if(DualRelators[i] 	== NULL) return(91);
		Exp_Surgery_Rel[i]	= (unsigned char **) NewHandle(4L);
		if(Exp_Surgery_Rel[i] == NULL) return(92);			
		KorLRelators[i] 	= (unsigned char **) NewHandle(4L);
		if(KorLRelators[i] 	== NULL) return(93);
		OutRelators[i] 		= (unsigned char **) NewHandle(4L);
		if(OutRelators[i] 	== NULL) return(94);
		Relators[i] 		= (unsigned char **) NewHandle(4L);
		if(Relators[i] 		== NULL) return(95);
		TopOfChain[i] 		= (unsigned char **) NewHandle(4L);
		if(TopOfChain[i] 	== NULL) return(96);	
		}
		
	Temp1 		= (unsigned char **) NewHandle(4L);
	if(Temp1 	== NULL) return(97);
	Temp2 		= (unsigned char **) NewHandle(4L);
	if(Temp2 	== NULL) return(98);
	Temp3 		= (unsigned char **) NewHandle(4L);
	if(Temp3 	== NULL) return(99);
	Temp4 		= (unsigned char **) NewHandle(4L);
	if(Temp4 	== NULL) return(100);
	Temp5		= (unsigned char **) NewHandle(4L);
	if(Temp5 	== NULL) return(101);
	Temp6		= (unsigned char **) NewHandle(4L);
	if(Temp6 	== NULL) return(102);
	Temp7		= (unsigned char **) NewHandle(4L);
	if(Temp7 	== NULL) return(103);
	Temp8		= (unsigned char **) NewHandle(4L);
	if(Temp8 	== NULL) return(104);
	Temp9		= (unsigned char **) NewHandle(4L);
	if(Temp9 	== NULL) return(105);
	Temp10		= (unsigned char **) NewHandle(4L);
	if(Temp10 	== NULL) return(106);
	Temp11		= (unsigned char **) NewHandle(4L);
	if(Temp11 	== NULL) return(107);
	Temp12		= (unsigned char **) NewHandle(4L);
	if(Temp12 	== NULL) return(108);
	Temp13		= (unsigned char **) NewHandle(4L);
	if(Temp13 	== NULL) return(109);
	Temp14		= (unsigned char **) NewHandle(4L);
	if(Temp14 	== NULL) return(110);
	Temp15		= (unsigned char **) NewHandle(4L);
	if(Temp15 	== NULL) return(111);
	
	Inst = (unsigned char *) NewPtr(30000);		/* Get some storage for any reasonable input. */
	if(Inst == NULL) return(112);
	
	/*	Click_On(FALSE);	*/
			
	ReadDateTime(&Seconds);
	srand(abs(LoWord(Seconds)));
	
/*	printf("\f");
	SetWTitle(FrontWindow(),"\pHEEGAARD");		*/

	/*  When in batch mode as "heegaard_is_realizable" we don't
	 *  want to write to Heegaard_Results so we just set "myout"
	 *  to be "stdout".
	 */
	if(open_Heegaard_Results)
	    {
		if((myout = fopen("Heegaard_Results","a+")) == NULL)
		    printf("\nUnable to open the file 'Heegaard_Results'.");
	    }
	else
	    {
		myout = stdout;
	    }
		
	Did_Exponent_Surgery 		= FALSE;
	Did_Cutting_Disk_Surgery 	= FALSE;
	
	return(FALSE);	
}

void Batch_Message(char *ans){
    if( myout == stdout ){
	fprintf(stdout, "HEEGAARD_ANS: (%s)\n", ans);
    }
}

int Generate_Orbits_Under_Auts(void)
{
	int						i,
							j;
							
	unsigned int			NumPres;						
	
	long					MaxPresLength;
	unsigned long			InitLength;
			
	register unsigned char 	*p,
							*q;
	
	unsigned char *			r;
	
	unsigned int Whitehead_Graph();
							
	Compute_Stabilizers = FALSE;
	Left[0] = Right[0] = INFINITE;
	Num_Level_Transformations = TRUE;
			
	r = (unsigned char *) NewPtr(100);		
	printf("\n\nENTER THE MaxSavedLength OF A PRESENTATION.     ");
	ReadString((char *)r, GetPtrSize(r));
	sscanf((char *) r,"%ld",&MaxPresLength);		
	DisposePtr((char *) r);
		
	if(Save_Pres(ReadPres,0,Length,1,2,1,0,0))
		{
		Compute_Stabilizers = FALSE;
		return(0);
		}
	InitLength = Length;	
	ReadPres = NumFilled - 1;
	
	do	
		{	
		/***********************************************************************************
			Copy the presentation, for which we want level transformations, into Relators[].
		************************************************************************************/
		
		NumRelators 	= NR[ReadPres];
		NumGenerators 	= NG[ReadPres];
		if(NumGenerators > 2) continue;
		Vertices = 2*NumGenerators;	
		for(i = 1; i <= NumRelators; i++)
			{
			ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[ReadPres][i]));
			if((q = *Relators[i]) == NULL) return(TOO_LONG);			
			p = *SUR[ReadPres][i];
			while(*q++ = *p++) ;
			}
		Length = SURL[ReadPres];
		Fill_A(NumRelators);
		if((A[0][2] > A[1][2]+A[2][3]) ||
			(A[0][3] > A[0][2]+A[2][3]) ||
			(A[0][2] > A[0][1]+A[0][3]) ||
			(A[0][3] > A[0][1]+A[0][2]))
			printf("\n		Pres %d does not have minimal length.",ReadPres + 1);
		
		if(Whitehead_Graph())
			{
			SysBeep(5);
			printf("\n		Whitehead_Graph returns an error on pres %d.",ReadPres + 1);
			goto _NEXTPRES;
			}
	
		if(j = Genus_Two_New_Relator(MaxPresLength,InitLength)) return(j);

_NEXTPRES:								
		ReadPres++;
		}
	while(ReadPres < NumFilled);
	
	printf("\n\n	Looking for new presentations obtained by reducing the above presentations");
	printf("\n								to minimal length. . .\n");
	
	NumPres = NumFilled;
	for(ReadPres = 0; ReadPres < NumPres; ReadPres++)
		{
		NumRelators 	= NR[ReadPres];
		NumGenerators 	= NG[ReadPres];
		if(NumGenerators > 2) continue;
		Vertices = 2*NumGenerators;	
		for(i = 1; i <= NumRelators; i++)
			{
			ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[ReadPres][i]));
			if((q = *Relators[i]) == NULL) return(TOO_LONG);			
			p = *SUR[ReadPres][i];
			while(*q++ = *p++) ;
			}
		Length = SURL[ReadPres];
		OrigLength = Length;
		Find_Flow_A(NORMAL,FALSE);
		if(Length < OrigLength)
			{
			printf("\n		From %4d, Delta Length = %5ld, Auts = %5ld\n",
				ReadPres +1,Length - OrigLength,Automorphisms);
			j = In_File();
	 		if(j == TOO_LONG) return(5);
	 		if(j == NumFilled - 1)
	 			{
				 if(NumFilled >= MAX_SAVED_PRES - 3) return(5);
				 if(BytesUsed > BytesAvailable)
				{
				if(UserSaidQuit) return(5);
				if(UserSaidQuit = User_Says_Quit()) return(5);
				}
	 			}			
			}	
		}
	
	Test_LT_For_Pseudo_Min();
	
	return(0);	
}


int Genus_Two_New_Relator(long MaxPresLength,unsigned long InitLength)
{
	/******************************************************************************************
		New_Relator() looks for a pair of relators that it can bandsum together to form a 
		new relator. It makes use of data about the Heegaard diagram found by the routines in
		Diagram_Main().
			The routine only looks for and performs certain restricted kinds of bandsums,
		namely those which correspond to bandsumming two relators together along a band which
		corresponds to an edge of the dual Heegaard diagram. These edges correspond to "corners"
		of the major faces of the original Heegaard diagram. Note that forming bandsums in 
		this restricted fashion is the only way to avoid forming "waves" in the diagram.
	******************************************************************************************/
		
	register unsigned char 	*p,
							*q;
							
	register unsigned int 	c,
							d,
							e,
							v,
							w;
							
	unsigned char 			s,
							t,
							x,
							y,
							sx,
							sy,
							**Temp;
													
	int 					DoMult,
							s1,
							s2,
							ss;
							
	unsigned int 			edge,
							edgeLE,
							edgeLS,
							edgeRE,
							edgeRS,
							EL[6*VERTICES],
							ee,
							eLE,
							eLS,
							eRE,
							eRS,
							i,
							j,
							vertex,
							vertexLS,
							vertexRS,
							vLS,
							vRS,
							vrt;
	
	long 					Diff,
							LOR[MAXNUMRELATORS + 1],
							length,
							slength;
	
	unsigned long			max;
	
	for(d = 1,max = 0L; d <= NumRelators; d++) if(LR[d] > max) max = LR[d];
	if(SRError == 3) for(d = 1; d <= NumRelators; d++)
		LOR[d] = GetHandleSize((char **) OutRelators[d]);
	Minimum = BIG_NUMBER;
	s1 = s2 = FALSE;
		
	for(d = 1; d <= 2*NumEdges; d++) EL[d] = d;
	for(ss = 2*NumEdges; ss > 0; ss --)
		{
		/**************************************************************************************
			Choose a random edge of the dual diagram along which to perform a bandsum.
			Remove the chosen edge from the list of available edges so that we won't try it
			again.
		**************************************************************************************/	
		
		c = abs(rand()) % ss;
		c++;
		v = EL[ss];
		EL[ss] = EL[c];
		EL[c] = v;
		
		/**************************************************************************************
			Determine which edges of the original diagram the chosen edge, of the dual
			diagram, joins.
		**************************************************************************************/	
		
		for(v = c = 0; c + VWG[v] < EL[ss]; v++) c += VWG[v];
		w = EL[ss] - c;
		for(c = ee = 0,d = FV[v]; c < w; c++)
			{
			ee += A[v][d];
			d = CO[v][d];
			}
		e = ee - 1;
		if(v & 1)
			{	
			e = OSA[v] - e;
			if(e >= V[v]) e -= V[v];
			if(e) e--;
			else e = V[v] - 1;
			}
			
		/**************************************************************************************
				Check whether this edge we are trying to bandsum along has its ends on
				distinct relators. If it does, we can go ahead and form a bandsum.
		**************************************************************************************/		
		
		p = q = *DualRelators[(v >> 1) + 1] + e;
		q++;
		if(! *q) q = *DualRelators[(v >> 1) + 1];
		switch(*p + *q)
			{
			case 131:
				{
				if(s1)
					DoMult = FALSE;
				else
					{
					s1 = TRUE;
					DoMult = TRUE;
					}	
				break;
				}			
			case 163:
				{
				if(s2)
					DoMult = FALSE;
				else
					{
					s2 = TRUE;
					DoMult = TRUE;
					}	
				break;
				}
			case 195:
				{
				if(s1)
					DoMult = FALSE;
				else
					{
					s1 = TRUE;
					DoMult = TRUE;
					}	
				break;
				}
			default:
				{
				DoMult = FALSE;
				break;
				}
			}
		if(DoMult)		
			{
			x 		= *p;
			y 		= *q;
			length	= 0L;
			vertex  = v;
			edgeRE 	= ee;
			edgeLE 	= edgeRE - 1;
			edgeRE  = edgeRE % V[v];	
			edge 	= edgeRE;
			e 		= edge;
			
			/*******************************************************************************
				Two distinct relators have come from distinct vertices and are adjacent to 
				each other at vertex v. Follow these two relators along until they diverge.
				The region where they are parallel is the "cancellation region".
			*******************************************************************************/	
			
			do
				{
				if(v & 1)
					w = v - 1;
				else
					w = v + 1;
				e = OSA[v] - e;
				if(e >= V[v]) e -= V[v];
				length ++;
				if(length > max)
					{
					/***********************************************************************
						This is an error that should not occur. But if the impossible
						happens, we fail gracefully by returning and pretending we could not
						find the bandsum.
					***********************************************************************/
					
					printf("\nlength > max in a bandsum!");
					Minimum = BIG_NUMBER;
					return(5);
					}	
				v = FV[w];
				d = A[w][v];
				while(d <= e)
					{
					v = CO[w][v];
					d += A[w][v];
					}	
				if(e == (d - 1))	
					{
					/***********************************************************************
						If e = d - 1, then we have found the end of the cancellation region.
					***********************************************************************/
					
					edgeRS = B[w][v] - e;
					vertexRS = v;
					v = CO[w][v];
					d = d % V[w];
					edgeLS = B[w][v] - d;
					vertexLS = v;
					
					/***********************************************************************
						Determine which edge of the dual diagram corresponds to the end of
						the cancellation region and delete it from the list of available
						edges.
					***********************************************************************/	
					
					for(d = 0,c = 1; d < w; d++) c += VWG[d];
					v = FV[w];
					d = A[w][v];
					while(d <= e)
						{
						c++;
						v = CO[w][v];
						d += A[w][v];
						}
					if(EL[c] == c)
						{
						ss --;
						v = EL[ss];
						EL[ss] = c;
						EL[c] = v;
						}	
					else
						{
						for(d = 1; d < ss && (EL[d] != c); d++) ;
						if(d < ss)
							{
							ss --;
							v = EL[ss];
							EL[ss] = c;
							EL[d] = v;
							}
						}		
					break;
					}
				e = B[w][v] - e;
				}	
			while(v != vertex || e != edge);
			
			/*******************************************************************************
				Increment the number of bandsums examined and use the information about the
				length of the cancellation region to determine how the total length of the
				relators would be changed by performing this bandsum and replacing one of
				the relators with the result.
			*******************************************************************************/							
			
			Band_Sums ++;
			for(d = 0; d < 2; d++)
				{
				if(d == 0)
					{
					s = x;
					t = y;
					}
				else
					{
					s = y;
					t = x;
					}
					
				/***************************************************************************
						Compute Diff, which is the amount by which the length of the
						presentation would be changed by performing this bandsum and
						replacement. If Diff < Minimum and !GoingUp or Diff > 0, replace
						the current saved set of parameters with the set of parameters we
						have just found.
				***************************************************************************/			
				
				if(s < 'a')
					Diff = LOR[s - 64] - (length << 1) - 1;
				else	
					Diff = LOR[s - 96] - (length << 1) - 1;
				if((Diff > 0L && (Length + Diff) <= MaxPresLength) ||
					(Diff == 0L && Length == InitLength) ||
					(Diff < 0L && Length + Diff < InitLength))
					{
					Minimum = Diff;
					sx 		= x;
					sy 		= y;
					eRS 	= edgeRS;
					eLS 	= edgeLS;
					eRE 	= edgeRE;
					eLE 	= edgeLE;
					vRS 	= vertexRS;
					vLS 	= vertexLS;
					vrt 	= vertex;
					slength = length;
					if(t < 'a')
						Word1 = t - 64;
					else	
						Word1 = t - 96;
					if(s < 'a')
						Word2 = s - 64;
					else	
						Word2 = s - 96;
					Long_Mult(sx,sy,eRS,eLS,eRE,eLE,vRS,vLS,vrt,slength,max);
					
					SRError = TRUE;	
					LR[0] = GetHandleSize((char **) OutRelators[Word1]) - 1;
					i = Word1;
					HLock((char **) OutRelators[i]);
					Word1 = abs(Compare(*OutRelators[i]));
					HUnlock((char **) OutRelators[i]);
					if(Word1 == 0)
						{
						SysBeep(5);
						printf("\n\nError in forming a bandsum. Relators are probably too long.");
						SaveMinima = TRUE;
						return(5);
						}
					if(Micro_Print)
						{
						LR[0] = GetHandleSize((char **) OutRelators[Word2]) - 1;
						i = Word2;
						HLock((char **) OutRelators[i]);
						Word2 = abs(Compare(*OutRelators[i]));
						HUnlock((char **) OutRelators[i]);
						if(Word2 == 0)
							{
							SysBeep(5);
							printf("\n\nError in forming a bandsum. Relators are probably too long.");
							SaveMinima = TRUE;
							return(5);
							}
						}		
					Temp 			= Temp3;
					Temp3 			= Relators[Word1];	
					Relators[Word1] = Relators[1];
					Relators[1] 	= Temp2;
					Temp2 			= Temp;
					LR[Word1]		= LR[1];
					LR[1]			= GetHandleSize((char **) Relators[1]) - 1;
					OrigLength 		= Length - GetHandleSize((char **) Temp3) + GetHandleSize((char **) Relators[1]);
					Length			= OrigLength;
					if(Micro_Print) Micro_Print_Bandsum();
	
					j = In_File();
			 		if(j == TOO_LONG) return(5);
			 		if(j == NumFilled - 1)
			 			{
						 if(NumFilled >= MAX_SAVED_PRES - 3) return(5);
						 if(BytesUsed > BytesAvailable)
						{
						if(UserSaidQuit) return(5);
						if(UserSaidQuit = User_Says_Quit()) return(5);
						}
			 			}

					for(i = 1; i <= NumRelators; i++)
						{
						ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[ReadPres][i]));
						if((q = *Relators[i]) == NULL) return(TOO_LONG);			
						p = *SUR[ReadPres][i];
						while(*q++ = *p++) ;
						LR[i] = GetHandleSize((char **) Relators[i]) - 1;
						}
					Length = SURL[ReadPres];
								
					}
				}		
			}
		}
	return(0);				
}
