#include "Heegaard.h"
#include "Heegaard_Dec.h"

/****************************** function prototypes *****************************************
L   33 Sort_Presentations_In_Memory(int F1)
L  144 qkst_compare(int i,int j)
L  187 qkst_swap(i,j)
L  198 Report(long Band_Sums,long NumDiagrams,unsigned int OnStack,unsigned int Starting_Pres,unsigned int Flag1,
            unsigned int Flag2,unsigned int Flag3,unsigned int Flag4,unsigned int Flag5,unsigned char * Ptr1)
L  627 Update_Bdry_Data(void)
L  705 Print_Bdry_Data(unsigned int WhichPres)
L  768 Print_Bdry_Data2(unsigned int WhichPres)
L  831 Fatal_Error(void)
L  941 Print_Relators(unsigned char ***MyRelators,int MyNumRelators)
L  950 Print_Relators2(unsigned char ***MyRelators,int MyNumRelators)
L  957 Micro_Print_Reset(void)
L  967 Micro_Print_Freely_Reduce(unsigned long length, unsigned long origlength)
L  975 Micro_Print_Dualize(void)
L  981 Micro_Print_Bandsum(void)
L  998 Micro_Print_Do_Aut(unsigned int Source, unsigned int NumReps)
L 1040 Print_DelRelators(void)
L 1049 Print_DualRelators(int F1,int F2,int Pres,int HS)
L 1063 Print_OutRelators(int F1,int F2,int Pres,int HS)
L 1077 Print_SLR(int i,int Found_L_Annulus)
L 1092 Display_A_Diagram(int,int,int)
L 1211 Display_Diagrams(void)
L 1505 Batch_Report(int*)  
L 2187 Print_SFComp(int*)       
********************************************************************************************/

int     *Table = NULL;

int Sort_Presentations_In_Memory(int F1)
{    
    int         	i; 
    unsigned int 	SOnStack; 
         
    int         	qkst_compare();
    void        	qkst_swap();
    
    if(Batch == 4)
    	{
		if(NumFilled < (MyMaxSavedPres + 3))
			{
			printf("\nHeegaard saved %u presentation(s), performed %lu automorphism(s), %lu Sep-Vert-Slide(s),",
				NumFilled,TotalAuts,Num_Level_Slides);
			printf("\nexamined %ld bandsum(s), examined %ld diagram(s), dualized %lu diagram(s), Mergers %u, ToDo %u.\n",
				Band_Sums,NumDiagrams,NumDualized,Mergers,OnStack);
			}
		return(0);
		}	
        
    Table = (int*) NewPtr((sizeof(int)*NumFilled));
    if(Table == NULL) Mem_Error();
    
    for(i = 0; i < NumFilled; i++) Table[i] = i;
    
    qksort(NumFilled);
    
    if(Batch == FALSE)
    	{
		printf("\n\nHIT 'v' TO REVIEW THESE SORTED PRESENTATIONS.");
		printf("\nOR HIT ANY OTHER KEY TO CONTINUE.");    
		switch(WaitkbHit())
			{
			case 'v':
				if(F1 == 2)
					Report_SPC(Table);
				else
					Report(Band_Sums,NumDiagrams,OnStack,0,1,0,0,0,1,NULL);
				break;
			default:
				break;        
			}
		GET_RESPONSE:
		if(F1)
			{
			printf("\n\nFIND_CANONICAL_ORBIT_REPS? HIT 'y' OR 'n'.");
			switch(WaitkbHit())
				{
				case 'y':
					SOnStack = OnStack;
					printf("\n	This may take awhile. Hit 's' to get a status report.");
					Find_Canonical_Orbit_Reps(Table,F1);
					OnStack = SOnStack;
					break;
				case 'n':
					break;
				default:
					goto GET_RESPONSE;
				}
			}	            
    	}
    if(Batch == 10 || Batch == 11 || Batch == 53)
    	{
    	if(B10B11HSReps || Batch == 53)
    		{
    		if(Batch == 53) printf("\n");
			printf("\nLooking for Heegaard Splitting Reps. . . .\n");
			SOnStack = OnStack;
			FoundSF = FoundFiniteSF = FALSE;
			switch(Find_Canonical_Orbit_Reps(Table,F1))
				{
				case INTERRUPT:
					OnStack = SOnStack;
					DisposePtr((int*) Table);
					return(INTERRUPT);			
				case TOO_LONG:
					FoundSF = FoundFiniteSF = FALSE;
					break;
				}
			}
		if(Batch == 10 || Batch == 11)	Batch_Report(Table);
		if(Batch == 53) printf("\n");
    	if(NumFilled < MyMaxSavedPres)
    		{
    		printf("\nHeegaard saved %u presentation(s), performed %lu automorphism(s), %lu Sep-Vert-Slide(s),",
    			NumFilled,TotalAuts,Num_Level_Slides);
        	printf("\nexamined %ld bandsum(s), examined %ld diagram(s), dualized %lu diagram(s), Mergers %u, ToDo %u.\n",
        		Band_Sums,NumDiagrams,NumDualized,Mergers,OnStack);
    		}			  	
    	}
    if(Batch == 14)
    	{
    	if((NumFilled > 1) && (B14B15PrintPres || (NG[Table[NumFilled -1]] <= 1) || (NR[Table[NumFilled -1]] == 0)))
    		Report_SPC(Table);
    	else
    		{	
			SOnStack = OnStack;
			if(Find_Canonical_Orbit_Reps(Table,F1) == INTERRUPT)
				{
				OnStack = SOnStack;
				DisposePtr((int*) Table);
				return(INTERRUPT);
				} 
			}   	    	
    	}
    if(Batch == 15 && NumFilled > 1) Report_SPC(Table);
    		
    DisposePtr((int*) Table);
    return(0);
}

int qkst_compare(int i,int j)
{
    register unsigned char  	*p,
                            	*q;
                            
    int                        	Ti,
                            	Tj,
                            	k;
                            
                
    Ti = Table[i];
    Tj = Table[j];
    if(Ti == Tj) return(0);
    if(ComponentNum[Ti] > ComponentNum[Tj]) return(1);
    if(ComponentNum[Ti] < ComponentNum[Tj]) return(-1);
    if(NG[Ti] < NG[Tj]) return(1);
    if(NG[Ti] > NG[Tj]) return(-1);
    if(NR[Ti] < NR[Tj]) return(1);
    if(NR[Ti] > NR[Tj]) return(-1);
    if(SURL[Ti] < SURL[Tj]) return(1);
    if(SURL[Ti] > SURL[Tj]) return(-1);
    for(k = 1; k <= NR[Ti]; k++)
        {
        if(GetHandleSize((char **) SUR[Ti][k]) > GetHandleSize((char **) SUR[Tj][k])) return(1);
        if(GetHandleSize((char **) SUR[Ti][k]) < GetHandleSize((char **) SUR[Tj][k])) return(-1);
        }
    for(k = 1; k <= NR[Ti]; k++)
        {
        p = *SUR[Ti][k];
        q = *SUR[Tj][k];
        while(*p && *p == *q)
            {
            p++;
            q++;
            }
        if(*p < *q) return(1);
        if(*p > *q) return(-1);
        }    
    if(Ti < Tj) return(1);
    if(Ti > Tj) return(-1);
    return(0);
}						

void qkst_swap(i,j)
int       	i,
            j;
{
    int            Temp;
    
    Temp         = Table[i];
    Table[i]     = Table[j];
    Table[j]     = Temp;
}            

int Report(long Band_Sums,long NumDiagrams,unsigned int OnStack,unsigned int Starting_Pres,unsigned char Flag1,
            unsigned char Flag2,unsigned char Flag3,unsigned char Flag4,unsigned char Flag5,unsigned char * Ptr1)
{
    /******************************************************************************************
        Report() is an output routine called automatically when Heegaard terminates. It
        can also be called by the user when Heegaard has been interrupted.
    ******************************************************************************************/ 
            
    unsigned char   *p,
                    x,
                    y;
    
    int             NumRelators;
                     
    unsigned int    i,
                    j,
                    k,
                    m,
                    n;                                
    
    unsigned long   Length;
    
    if(Batch == 4)
    	{
    	for(i = j = k = m = 0; i < NumFilled; i++) 
        	{
        	j += SURNumX[i];
        	k += NumLoops[i];
            if(PRIM[i] == 6 || PRIM[i] == 106) m++;  	
        	}
        printf("\n\nHeegaard performed %lu automorphism(s), %lu Sep-Vert-Slide(s), examined %ld bandsum(s),",
            TotalAuts,Num_Level_Slides,Band_Sums);
        printf("\nexamined %ld diagram(s), dualized %lu diagram(s), Hits %u, Loops %u, Mergers %u,",
            NumDiagrams,NumDualized,j,k,Mergers);
        printf("\nNumFP %u, ToDo %u.\n",m,OnStack);
        return(0);
        }
        
    
    if(Flag5) Update_Bdry_Data(); 
    
    if(Flag5 && (Batch != 3))
    printf("\n\n------------------------------------");
            
    if(BdryData && Flag5 && (NumFilled > 1))
        {
        printf("\n");
        Print_Bdry_Data(Starting_Pres);
        }
	
	if(Batch != 3)
    printf("\n\nThe initial presentation was: %s",PresName);    
      
    if(Batch != 3)        
    for(n = Starting_Pres; n < NumFilled; n++)
        {
        i = n;
        if(Flag1) i = Table[n];
        if(Flag2)
        	{
        	if(Flag4 == 1 && !Ptr1[i]) continue;
        	if(Flag4 != 1 && Ptr1[i] != Flag4) continue;
        	}
        NumRelators = NR[i];
        Length = SURL[i];
        printf("\n\nPresentation %d  of Summand %u:  Gen  %d  Length  %lu  From Pres %u  NumHits %d NumLoops %d HegSplNum %u ",
        i+1,ComponentNum[i],NG[i],Length,FR[i]+1,SURNumX[i],NumLoops[i],HegSplNum[i]);
        
        switch(PRIM[i])
            {
            case 1:
            case 101:
                printf("DR");
                break;
            case 2:
            case 102:
                printf("IP");
                break;
            case 3:
            case 103:
                printf("LS");
                break;
            case 4:
            case 104:
                printf("1G");
                break;
            case 5:
            case 105:
                printf("S3");
                break;
            case 6:
            case 106:
                printf("FP");
                break;
            case 7:
            case 107:
                printf("BC");
                break;
            case 8:
            case 108:
                printf("PM");
                break;
            case 9:
            case 109:
                printf("PM");
                break;
            case 10:
            case 110:
                printf("BC");
                break;
            case 11:
            case 111:
                printf("CF");
                break;
            case 12:
            case 112:
                printf("ER");
                break;
            case 13:
            case 113:
                printf("Er");
                break;                
            case 20:
            case 120:
                printf("NC");
                break;
            case 30:
            case 130:
            case 40:
            case 140:
                printf("MG");
                break;
            case 60:
            case 160:
                printf("PP");
                break;    
            case 70:
            case 170:
                printf("Lt");
                break;
            case 75:
            case 175:
                printf("LT");
                break;    
            case 80:
            case 180:
                printf("A2");
                break;    
            default:
                break;
            }
                
        switch(UDV[i])
            {
            case SPLIT:
                j = Daughters[i];
                m = j + NCS[i] - 1;
                switch(PRIM[j])
                    {
                    case 40:
                    case 140:
                        printf("\nThe presentation dual to presentation %d 'split' into presentations %u and %u of M",
                                    i + 1,j + 1,m + 1);
                        printf("\nwhich correspond to summands %u and %u of M.",ComponentNum[j],ComponentNum[m]);                    
                        break;
                    default:
                        printf("\n'split' into presentations %u and %u corresponding to summands %u and %u of M.",
                        j + 1,m + 1,ComponentNum[j],ComponentNum[m]);
                        break;                    
                    }
                if(BdryData) for(m = 0; m < NCS[i]; m++) Print_Bdry_Data(j + m);
                break;
            
            case GENERIC_LENS_SPACE:
                if(LSP[i] > 4)
                    printf("\npresents a Lens space of the form L(%lu,Q).",LSP[i]);
                else
                if(LSP[i] == 1)
                    printf("\npresents the 3-Sphere.");
                else
                    printf("\npresents a Lens space of the form L(%lu,1).",LSP[i]);
                break;        
            
            case THREE_SPHERE:
                printf("\npresents the 3-Sphere.");
                break;

            case NOT_CONNECTED:
                printf("\nthe diagram is not connected.");
                break;
                    
            case S1_X_S2:
                if(NG[i] == 1)
                    printf("\npresents S1 X S2.");
                else
                    printf("\npresents %d copies of S1 X S2.",NG[i]);
                break;
            
            case S1_X_D2:
                if(NG[i] == 1)
                    printf("\npresents S1 X D2.");
                else
                    printf("\npresents %d copies of S1 X D2.",NG[i]);
                break;
            
            case S1_X_X2:
                if(NG[i] == 1)
                    printf("\npresents S1 X ?2.");
                else
                    printf("\npresents %d copies of S1 X ?2.",NG[i]);
                break;
            
            case MISSING_GEN_DONE2:
                printf("\npresents: ");
                if(N1H[ComponentNum[i]] == 1)
                    printf("1 copy of I X D2, ");
                else
                    printf("%d copies of I X D2, ",N1H[ComponentNum[i]]);
                if(NS1XS2[ComponentNum[i]] == 1)
                    printf("1 copy of S1 X S2, ");
                else
                    printf("%d copies of S1 X S2, ",NS1XS2[ComponentNum[i]]);
                if(NS1XD2[ComponentNum[i]] == 1)
                    printf("and 1 copy of S1 X D2.");
                else
                    printf("and %d copies of S1 X D2.",NS1XD2[ComponentNum[i]]);                                        
                printf("\nHeegaard was not able to unambiguously determine the boundary components.");
                break;
                
            case MISSING_GEN_DONE1:
                printf("\npresents: ");
                if(N1H[ComponentNum[i]] == 1)
                    printf("1 copy of I X D2, ");
                else
                    printf("%d copies of I X D2, ",N1H[ComponentNum[i]]);
                if(NS1XS2[ComponentNum[i]] == 1)
                    printf("1 copy of S1 X S2, ");
                else
                    printf("%d copies of S1 X S2, ",NS1XS2[ComponentNum[i]]);
                if(NS1XD2[ComponentNum[i]] == 1)
                    printf("and 1 copy of S1 X D2.");
                else
                    printf("and %d copies of S1 X D2.",NS1XD2[ComponentNum[i]]);                            
                break;
                                
            case KNOWN_LENS_SPACE:
                switch(LSP[i])
                    {
                    case 0L:
                        printf("\npresents S1 X S2.");
                        break; 
                    case 1L:
                        printf("\npresents the 3-Sphere.");
                        break;
                    default:
                        printf("\npresents the Lens space L(%lu,%lu).",LSP[i],LSQ[i]);
                        break;
                    }
                break;
            
            case SEP_PAIRS:
                if(PRIM[i] >= 100)
                    {
                    printf("\ndual of presentation %u is pseudo-minimal.",
                            Daughters[i] + 1);
                    printf(" The pair of vertices (%c,%c) separate.",
                        x = LSP[i],y = LSQ[i]);                    
                    }
                else                        
                    printf("\nthe pair of vertices (%c,%c) separate the diagram.",
                    x = LSP[i],y = LSQ[i]);            
                break;
            
            case ANNULUS_EXISTS:
                p = *SUR[i][0];
                x = *p++;
                y = *p++;                
                printf("\nThe pair of vertices (%c,%c) separate the diagram.",x,y);
                printf("\nThe component consisting of vertice(s) {%c",*p);
                p++;
                while((x = *p) != '@')
                    {
                    printf(",%c",x);
                    p++;
                    }
                printf("}");    
                p++;        
                printf("\nlies in an annulus which swallows the component and otherwise follows the curve:\n");
                while(*p)
		   			fputc(*p++,stdout);
                break;
            
            case V2_ANNULUS_EXISTS:
                p = *SUR[i][0];
                printf("\nThere exists an annulus which swallows vertice(s) {%c",*p);
                p++;
                while((x = *p) != '@')
                    {
                    printf(",%c",x);
                    p++;
                    }
                printf("}");    
                p++;        
                printf("\nand otherwise follows the curve:\n");
				while(*p) fputc(*p++,stdout);

                break;                
            
            case DELETED_RELATOR:
                break;
            
            case NON_UNIQUE_4:
                printf("\nthe diagram is not unique because there is a generator which appears");
                printf("\nwith only one exponent and that exponent is greater than 6.");
                break;
            
            case NON_UNIQUE_3:
                printf("\nthe diagram is not unique because there is a generator which appears");
                printf("\nonly with exponent 5.");
                break;
            
            case NON_UNIQUE_2:
                printf("\nthe diagram is not unique because there is a generator which appears");
                printf("\nonly with exponent 3 or only with exponent 4 and this exponent occurs more than once.");
                break;
            
            case NON_UNIQUE_1:
                printf("\nthe diagram is not unique because there is a generator which appears");
                printf("\nwith only one exponent, either 3,4 or 6, and a needed symmetry does not exist.");
                break;
            
            case DUPLICATE:
                printf("\nis a duplicate of presentation %d of summand %u.",
                Daughters[i] + 1,ComponentNum[Daughters[i]]);
                break;                                                                                                                break;        
            
            default:
                {
                j = PRIM[i];
                switch(j)
                    {
                    case 8:
                    case 108:
                        printf("\ndual of presentation %u is pseudo-minimal.",
                        Daughters[i] + 1);
                        break;
                    case 70:
                        if(QPM[i])
                            {
                            printf("\nvia level transformations of sep_pairs,");
                            printf(" is quasi-pseudo-minimal.");
                            }
                        else        
                            printf("\nvia level transformations of sep_pairs.");
                        break;
                    case 75:
                        if(QPM[i])
                            {
                            printf("\nvia general level transformations,");
                            printf(" is quasi-pseudo-minimal.");
                            }
                        else
                            printf("\nvia general level transformations.");
                        break;    
                    case 170:
                        printf("\nvia level transformations of sep_pairs and ");
                        printf("dual of presentation %u is pseudo-minimal.",
                        Daughters[i] + 1);
                        break;
                    case 175:
                        printf("\nvia general level transformations and ");
                        printf("dual of presentation %u is pseudo-minimal.",
                        Daughters[i] + 1);
                        break;    
                    default:
                        if(j >= 100)
                            printf("\ndual of presentation %u is pseudo-minimal.",
                            Daughters[i] + 1);
                        if(QPM[i])
                            printf("\nis quasi-pseudo-minimal.");                                
                        break;
                    }
                break;
                }                                                                                
            }
                
        printf("\n");
        for(j = 1; j <= NumRelators; j++) printf("\n    %s",*SUR[i][j]);           
        }                    
    
    if(Batch == 3) Flag3 = FALSE;
    if(Flag3 && Flag5)
        {
        for(i = j = k = m = 0; i < NumFilled; i++) 
        	{
        	j += SURNumX[i];
        	k += NumLoops[i];
            if(PRIM[i] == 6 || PRIM[i] == 106) m++;  	
        	}
        printf("\n\nHeegaard performed %lu automorphism(s), %lu Sep-Vert-Slide(s), examined %ld bandsum(s),",
            TotalAuts,Num_Level_Slides,Band_Sums);
        printf("\nexamined %ld diagram(s), dualized %lu diagram(s), Hits %u, Loops %u, Mergers %u,",
            NumDiagrams,NumDualized,j,k,Mergers);
        printf("\nNumFP %u, ToDo %u.\n",m,OnStack);
        if(NumSepAnnuli == 1)
			printf("\nHeegaard found one diagram containing a separating annulus. Scroll back for details.");  
		if(NumSepAnnuli > 1)
			printf("\nHeegaard found %u diagrams containing a separating annulus. Scroll back for details.",NumSepAnnuli); 
		if(NumNonSepAnnuli == 1)
			printf("\nHeegaard found one diagram containing a nonseparating annulus. Scroll back for details.");			       		
		if(NumNonSepAnnuli > 1)
			printf("\nHeegaard found %u diagrams containing a nonseparating annulus. Scroll back for details.",NumNonSepAnnuli);
		if(NumRelTooLong == 1)
			printf("\nDeleting primitives produced one relator which was too long. Scroll back for details.");				
		if(NumRelTooLong > 1)
			printf("\nDeleting primitives produced %u relators which were too long. Scroll back for details.",NumRelTooLong);	
		if(CouldNotRemove == 1)
			printf("\nHeegaard could not remove all pairs of separating vertices from one presentation. Scroll back for details.");			
		if(CouldNotRemove > 1)
			printf("\nHeegaard could not remove all pairs of separating vertices from %u presentations. Scroll back for details.",CouldNotRemove);
		if(NumErrors == 1)
			printf("\nOne error was detected. Scroll back for details.");
		if(NumErrors > 1)
			printf("\n%lu errors were detected. Scroll back for details.",NumErrors);
        }
    
    return(NO_ERROR);    
}        

void Update_Bdry_Data(void)
{
    int         h,
                i,
                j,
                k,
                m,
                n,
                NotClosed,
                NumNotClosed,
                NumUpdates;
    
    if(TotalComp == 1 || NoReport == FALSE) return;
                
    do
        {
        for(i = 0, NumUpdates = 0; i < NumFilled; i++)
            {
            if(UDV[i] == DUPLICATE)
                {
                m = ComponentNum[i];
                h = ComponentNum[Daughters[i]];
                if(CBC[m][0] == BDRY_UNKNOWN && CBC[h][0] < BDRY_UNKNOWN
                    && CBC[h][1] == BDRY_UNKNOWN)
                    {
                    CBC[m][0] = 1;
                    CBC[m][1] = BDRY_UNKNOWN;
                    NumUpdates ++;
                    }
                }
            if(UDV[i] == SPLIT)
                {
                m = ComponentNum[i];
                if(CBC[m][0] < BDRY_UNKNOWN && CBC[m][1] == BDRY_UNKNOWN)
                    {
                    n = ComponentNum[Daughters[i]] + NCS[i];
                    for(j = ComponentNum[Daughters[i]]; j < n; j++)
                    if(CBC[j][0] == BDRY_UNKNOWN)
                        {
                        NumUpdates ++;
                        CBC[j][0] = 1;
                        CBC[j][1] = BDRY_UNKNOWN;
                        }
                    }
                if(CBC[ComponentNum[i]][0] == BDRY_UNKNOWN)
                    {
                    /************************************************************************
                        Check whether all but at most one of the summands of M is closed.
                    ************************************************************************/
                        
                    h = Daughters[i];
                    n = ComponentNum[h] + NCS[i];
                    for(j = ComponentNum[h]; j < n; j++)
                        if(CBC[j][0] == BDRY_UNKNOWN) break;
                    if(j < n) continue;
                    for(j = ComponentNum[h],NumNotClosed = 0; j < n; j++)
                        if(CBC[j][1] != BDRY_UNKNOWN)
                            {
                            NotClosed = j;
                            if(++NumNotClosed > 1)
                            break;
                            }
                    if(j < n) continue;                
                    NumUpdates ++;
                    CBC[m][0] = 1;
                    if(NumNotClosed == 1)
                        for(k = 1; (CBC[m][k] = CBC[NotClosed][k]) < BDRY_UNKNOWN; k++) ;
                    else
                        CBC[m][1] = BDRY_UNKNOWN;    
                    for(j = ComponentNum[h]; j < n; j++) if(CBC[j][0])
                        CBC[m][0] += CBC[j][0] - 1;
                    }    
                }
            }    
        }
    while(NumUpdates);        
}

void Print_Bdry_Data(unsigned int WhichPres)
{
    unsigned int    i,
                    j,
                    k,
                    kk;
    
    k = kk = ComponentNum[WhichPres];
    if(UDV[WhichPres] == DUPLICATE)
        k = ComponentNum[Daughters[WhichPres]];
        
    if(CBC[k][0] == BDRY_UNKNOWN)
        {
        if(TotalComp == 1)
            printf("\n    Heegaard was unable to determine the boundary components of this manifold.");
        else
            printf("\n    Heegaard was unable to determine the boundary components of 'summand' %u of M.",kk);
        return;
        }
        
    for(i = 1,j = 0; CBC[k][i] < BDRY_UNKNOWN; i++) j += CBC[k][i];
    switch(j)
        {
        case 0:
            if(TotalComp == 1)
                printf("\n    This manifold is closed.");
            else
                printf("\n    'Summand' %u of M is closed.",kk);
            break;
        case 1:
            for(i = 1; CBC[k][i] < BDRY_UNKNOWN; i++) if(CBC[k][i])
                {
                if(TotalComp == 1)
                    printf("\n    This manifold has one boundary component of genus %u.",i);                
                else
                    printf("\n    'Summand' %u of M has one boundary component of genus %u.",kk,i);
                break;
                }
            break;        
        default:
            for(i = 1; CBC[k][i] < BDRY_UNKNOWN; i++) if(CBC[k][i] == j)
                {
                if(TotalComp == 1)
                    printf("\n    This manifold has %u boundary components of genus %u.",j,i);                
                else
                    printf("\n    'Summand' %u of M has %u boundary components of genus %u.",kk,j,i);
                return;
                }
            if(TotalComp == 1)
                printf("\n    This manifold has the following boundary components:");
            else 
                printf("\n    'Summand' %u of M has the following boundary components:",kk);
            for(i = 1; CBC[k][i] < BDRY_UNKNOWN; i++) if( (j = CBC[k][i]) )
                {
                if(j == 1)
                    printf("\n        %2u component  of genus %2u.",j,i);
                else
                    printf("\n        %2u components of genus %2u.",j,i);
                }
            break;    
        }
}

void Print_Bdry_Data2(unsigned int WhichPres)
{
    unsigned int    i,
                    j,
                    k,
                    kk;
    
    k = kk = ComponentNum[WhichPres];
    if(UDV[WhichPres] == DUPLICATE)
        k = ComponentNum[Daughters[WhichPres]];
        
    if(CBC[k][0] == BDRY_UNKNOWN)
        {
        if(TotalComp == 1)
            fprintf(H_Results,"\n    Heegaard was unable to determine the boundary components of this manifold.");
        else
            fprintf(H_Results,"\n    Heegaard was unable to determine the boundary components of 'summand' %u of M.",kk);
        return;
        }
        
    for(i = 1,j = 0; CBC[k][i] < BDRY_UNKNOWN; i++) j += CBC[k][i];
    switch(j)
        {
        case 0:
            if(TotalComp == 1)
                fprintf(H_Results,"\n    This manifold is closed.");
            else
                fprintf(H_Results,"\n    'Summand' %u of M is closed.",kk);
            break;
        case 1:
            for(i = 1; CBC[k][i] < BDRY_UNKNOWN; i++) if(CBC[k][i])
                {
                if(TotalComp == 1)
                    fprintf(H_Results,"\n    This manifold has one boundary component of genus %u.",i);                
                else
                    fprintf(H_Results,"\n    'Summand' %u of M has one boundary component of genus %u.",kk,i);
                break;
                }
            break;        
        default:
            for(i = 1; CBC[k][i] < BDRY_UNKNOWN; i++) if(CBC[k][i] == j)
                {
                if(TotalComp == 1)
                    fprintf(H_Results,"\n    This manifold has %u boundary components of genus %u.",j,i);                
                else
                    fprintf(H_Results,"\n    'Summand' %u of M has %u boundary components of genus %u.",kk,j,i);
                return;
                }
            if(TotalComp == 1)
                fprintf(H_Results,"\n    This manifold has the following boundary components:");
            else 
                fprintf(H_Results,"\n    'Summand' %u of M has the following boundary components:",kk);
            for(i = 1; CBC[k][i] < BDRY_UNKNOWN; i++) if( (j = CBC[k][i]) )
                {
                if(j == 1)
                    fprintf(H_Results,"\n        %2u component  of genus %2u.",j,i);
                else
                    fprintf(H_Results,"\n        %2u components of genus %2u.",j,i);
                }
            break;    
        }
}

void Fatal_Error(void)
{
    /******************************************************************************************
        Fatal_Error() is called when Heegaard has discovered that a presentation is not
        realizable by a Heegaard diagram. It prints a message to this effect, and also prints
        the offending presentation.
    ******************************************************************************************/
    
    if(Batch == FALSE) SysBeep(5);
    BdryData = FALSE;
    if(NR[ReadPres] == NumRelators && Compare_Pres(ReadPres))
        {
        printf("\n\n                    Presentation %d, from presentation %u, is not realizable. ",
            ReadPres + 1,FR[ReadPres] + 1);
        switch(PRIM[ReadPres])
            {
            case 1:
            case 101:
                printf("DR");
                break;
            case 2:
            case 102:
                printf("IP");
                break;
            case 3:
            case 103:
                printf("LS");
                break;
            case 4:
            case 104:
                printf("1G");
                break;
            case 5:
            case 105:
                printf("S3");
                break;
            case 6:
            case 106:
                printf("FP");
                break;
            case 7:
            case 107:
                printf("BC");
                break;
            case 8:
            case 108:
                printf("PM");
                break;
            case 9:
            case 109:
                printf("PM");
                break;
            case 10:
            case 110:
                printf("BC");
                break;
            case 11:
            case 111:
                printf("CF");
                break;
            case 12:
            case 112:
                printf("ER");
                break;
            case 13:
            case 113:
                printf("Er");
                break;                    
            case 20:
            case 120:
                printf("NC");
                break;
            case 30:
            case 130:
            case 40:
            case 140:
                printf("MG");
                break;
            case 60:
            case 160:
                printf("PP");
                break;    
            case 70:
            case 170:
                printf("Lt");
                break;
            case 75:
            case 175:
                printf("LT");
                break;    
            case 80:
            case 180:
                printf("A2");
                break;    
            default:
                break;
            } 
        printf("\n");       
        Print_Relators(Relators,NumRelators);        
        }
    else
        {    
        printf("\n\n                    This presentation, obtained from presentation %d, is not realizable.\n",
            ReadPres + 1);
        Print_Relators(Relators,NumRelators);
        }    
        
    UDV[ReadPres] = DONE;
}

void Print_Relators(unsigned char ***MyRelators,int MyNumRelators)
{
  int    i;
    
    for(i = 1; i <= MyNumRelators; i++) printf("\n    %s",*MyRelators[i]); 

    if(Batch == 0) printf("\n");
}

void Print_Relators2(unsigned char ***MyRelators,int MyNumRelators)
{
    int    i;
    
    for(i = 1; i <= MyNumRelators; i++) fprintf(H_Results,"\n    %s",*MyRelators[i]); 
}

void Micro_Print_Reset(void)
{
    int     i;
    
    printf("\n\nStarted with Presentation %d of Summand %d, Length %lu:\n",ReadPres + 1,CurrentComp,SURL[ReadPres]);    
    for(i = 1; i <= NumRelators; i++)  printf("\n    %s",*SUR[ReadPres][i]);  
   
    printf("\n");
}

void Micro_Print_Freely_Reduce(unsigned long length, unsigned long origlength)
{
    printf("\n\nFree reductions reduced the length of the current presentation from %lu to %lu.",
        length,origlength);
    printf("\nThe reduced presentation is:\n");
    Print_Relators(Relators,NumRelators);
}

void Micro_Print_Dualize(void)
{
    printf("\n\nDualized the current relators to get the following dual relators:\n");
    Print_Relators(Relators,NumRelators);
}

void Micro_Print_Bandsum(void)
{
    int        SNumRelators;
    
    SNumRelators = NumRelators;    
    printf("\n\nReplaced Relator %u with the following bandsum of Relator %u and Relator %u.",
        Word1,Word1,Word2);
    printf(" Delta Length = %ld.\n",Length - SLength);     
    NumRelators = 1;
    Print_Relators(Relators,NumRelators);
    NumRelators = SNumRelators;
    if(Word1 != 1)
        printf("\n\nAnd then swapped Relator %u and Relator 1.",Word1);
    printf("\n\nThe current presentation is:\n");
    Print_Relators(Relators,NumRelators); 
}

void Micro_Print_Do_Aut(unsigned int Source, unsigned int NumReps)
{
    char   			A,
                    a,
                    x;
                    
    int             i;
                    
    /********************************************************************************************
    	Specifying how a Whitehead automorphism acts on a set of generators is not well-defined
    until the location of a base-point has been specified. Here we adopt the convention that the
    base-point lies on the same side of a partition of vertices defining a Whitehead 
    transformation of the Whitehead Graph as the Sink vertex.
    ********************************************************************************************/
                                    
    A = ((Source/2) + 65);
    a = A + 32;
    
    if(Micro_Print)
        printf("\nDo Aut %u time(s): ",NumReps);
    else
        printf("\n%6lu) ",Num_Level_Transformations + 1);

	for(i = 0; i < Vertices; i += 2)
		{
		if(i == Source)  continue;
		if(VA[i/2] == 0) continue;
        x = (i/2) + 65;
        if(ZZ[i])
        	{
        	if(!ZZ[i+1]) printf("%c->%c%c  ",x,x,a);
        	}
        else
        	{
        	if(ZZ[i+1]) 
        		printf("%c->%c%c  ",x,A,x);
        	else
        		printf("%c->%c%c%c  ",x,A,x,a);
        	}
		}
}

void Print_DelRelators(void)
{
    int i;
        
    printf("\n");    
    for(i = 1; i <= NumRelators; i++) printf("\n    %s",*DelRelators[i]);

}

void Print_DualRelators(int F1,int F2,int Pres, int HS)
{
    register int            i;
    
    if(F2)
    	printf("\n\nThe 'Dual' Relators of the Diagram of Pres %d of HS %d are:\n",Pres,HS);    
    else
    	printf("\n\nThe 'Dual' Relators of Diagram %d are:\n",WhichInput + 1);
    
    for(i = 1; i <= NumGenerators; i++)  printf("\n    %s",*DualRelators[i]); 
        
    printf("\n\nNote: Dual Relators are read counter-clockwise about vertices A,B ... Z starting at edge 0.");  
}

void Print_OutRelators(int F1,int F2,int Pres,int HS)
{
    register int            i;
    
    if(F2)
    	printf("\n\nThe 'Out' Relators of the Diagram of Pres %d of HS %d are:\n",Pres,HS);    
    else
    	printf("\n\nThe 'Out' Relators of Diagram %d are:\n",WhichInput + 1);
    
    for(i = 1; i <= NumRelators; i++)  printf("\n    %s",*OutRelators[i]);  
       
    printf("\n");       
}

void Print_SLR(int i,int Found_L_Annulus)
{
    int j;
    
    if(Found_L_Annulus) return;
    CouldNotRemove ++;
    
    printf("\n\nCould not remove all separating pairs of vertices from the following presentation.\n");
    printf("And Level_Transformations() did not locate an annulus.");
    
    for(j = 1; j <= NumRelators; j++) printf("\n    %s",*SLR[i][j]);
    printf("\n");

}
        
int Display_A_Diagram(F1,Pres,HS)
{    
    int             Reply;
    
    unsigned int    SaveUDV;            
    
    unsigned long   SLSP,
                    SLSQ;

    DrawingDiagrams = TRUE;

    if((WhichInput != (MAX_SAVED_PRES - 1)) && Get_Relators_From_SUR(WhichInput))
        {
        printf("\n\n    Memory Error. Sorry!");
        goto _ERROR;        
        }
                
    if(Length == 0)
        {
        printf("\n\nThis is an empty presentation. There is nothing to display.");
        goto _ERROR;
        }

    if(WhichInput == 0 && UDV[0] == SPLIT)
        {
        switch(Find_Flow_A(NORMAL,FALSE))
            {
            case 1:
                break;
            case TOO_LONG:
                if(Batch == FALSE) SysBeep(5);
                if(F1)
                	printf("\n\n     Unable to display the diagram of Pres %d of HS %d. Sorry!",Pres,HS);                
                else
                	printf("\n\n     Unable to display diagram %d. Sorry!",WhichInput + 1);
                goto _ERROR;
            }
        if(Automorphisms)
            {
            if(Batch == FALSE) SysBeep(5);
            printf("\n\n                    NOTE!");
            printf("\n\n    Presentation 1 does not have minimal length.");
            printf("\n    Heegaard will only display diagrams of minimal length presentations.");
            printf("\n    Presentation 2 should give a diagram of a minimal length version of presentation 1.");
            goto _ERROR;
            }    
        }
    
    Fill_A(NumRelators);
    if(ComputeValences_A())
        {
        if(Batch == FALSE) SysBeep(5);
        if(F1)
            printf("\n\n     Unable to display the diagram of Pres %d of HS %d. Sorry!",Pres,HS);     
        else
        	printf("\n\n     Unable to display diagram %d. Sorry!",WhichInput + 1);
        goto _ERROR;    
        }    
    Get_Matrix();
    Check_Connected();
    SepPairs = Sep_Pairs(0,0,1);
    if(SepPairs == TOO_LONG)
        {
        if(Batch == FALSE) SysBeep(5);
        if(F1)
        	printf("\n\n     Unable to display the diagram of Pres %d of HS %d. Sorry!",Pres,HS);
        else
       		printf("\n\n     Unable to display diagram %d. Sorry!",WhichInput + 1);
        goto _ERROR;
        }
    SaveUDV = UDV[WhichInput];
    SLSP = LSP[WhichInput];
    SLSQ = LSQ[WhichInput];    
    if(SepPairs)
        {
        if(UDV[WhichInput] == 0) UDV[WhichInput] = SEP_PAIRS;
        if(V1 & 1)
            LSP[WhichInput] = V1/2 + 97;
        else
            LSP[WhichInput] = V1/2 + 65;
        if(V2 & 1)
            LSQ[WhichInput] = V2/2 + 97;
        else
            LSQ[WhichInput] = V2/2 + 65;        
        }
    NonPlanar = Planar(FALSE,TRUE);
    Reply = Print_Graph(TRUE,F1,Pres,HS);
    if(UDV[WhichInput] != ANNULUS_EXISTS && UDV[WhichInput] != V2_ANNULUS_EXISTS)
        UDV[WhichInput] = SaveUDV;
    LSP[WhichInput] = SLSP;
    LSQ[WhichInput] = SLSQ;    
    DrawingDiagrams = FALSE;
    return(Reply);

_ERROR:
    printf("\n\n         HIT 'n' TO SEE THE NEXT PRESENTATION.");
    if(!F1) printf("\n            HIT 'p' TO SEE THE PREVIOUS PRESENTATION.");
    GET_RESPONSE:    
    switch(WaitkbHit())
        {
        case 'n':
            Reply = 0;
            break;
        case 'p':
        	if(F1)
        		{
            	if(Batch == FALSE) SysBeep(5);
            	goto GET_RESPONSE;        		
        		}
            Reply = 1;
            break;
        default:
            if(Batch == FALSE) SysBeep(5);
            goto GET_RESPONSE;
        }
    DrawingDiagrams = FALSE;
    return(Reply);                        
}

void Display_Diagrams(void)
{
    unsigned char   *ptr = NULL;                    
                            
    int             NoSepPairs,
                    NumConnected,
                    Response,
                    SWhichInput;
    
    unsigned int    h,
                    i,
                    j,
                    k,
                    SaveUDV;
                    
    unsigned long   SLSP,
                    SLSQ;                                

    printf("\n\n                    Displaying Diagram Data. . .");

    if(NumFilled > 1)
        {
        printf("\n\n    REVIEW ALL PRESENTATIONS AVAILABLE ?  HIT 'y' OR 'n'.");
        GET_RESPONSE1:        
        switch(WaitkbHit())
            {
            case 'y':
                REVIEW:
                Report(Band_Sums,NumDiagrams,OnStack,0,0,0,0,0,1,NULL);
                printf("\n\n    CONTINUE TO REVIEW PRESENTATIONS ?  HIT 'y' OR 'n'.");
                GET_RESPONSE5:                
                switch(WaitkbHit())
                    {
                    case 'y':
                        goto REVIEW;
                    case 'n':
                        break;
                    default:
                        if(Batch == FALSE) SysBeep(5);
                        goto GET_RESPONSE5;
                    }
                break;    
            case 'n':
                break;
            default:
                if(Batch == FALSE) SysBeep(5);
                goto GET_RESPONSE1;
            }
        }        
    DrawingDiagrams = TRUE;
REDRAW:    
    printf("\n\n    a) SHOW DATA FOR ALL DIAGRAMS,");
    printf("\n    b) SHOW DATA FOR ALL DIAGRAMS THAT ARE CONNECTED AND HAVE NO SEPARATING PAIRS OF VERTICES,");
    printf("\n    c) OR SHOW YOUR CHOICE OF DATA FOR A PARTICULAR DIAGRAM ?");
    printf("\n\n    HIT 'a','b', OR 'c'");
    GET_RESPONSE2:
    Response = WaitkbHit();   
    switch(Response)
        {
        case 'a':
        case 'b':
        case 'c':
            break;
        default:
            if(Batch == FALSE) SysBeep(5);
            goto GET_RESPONSE2;
        }
    if(Response == 'c')
        {
        ptr = (unsigned char *) NewPtr(100);
        if(ptr == NULL) Mem_Error();
        printf("\n\nENTER A DIAGRAM FROM 1 TO %u WHOSE DATA YOU WISH TO SEE AND HIT 'return'.      ",NumFilled);
        for(i = j = 0; j < NumFilled; j++) if(SURL[j] == 0) i ++;
        if(i)
            {
            if(i == 1)
                {
                for(i = 0; i < NumFilled && SURL[i]; i++) ;
                i++;
                printf("\n\nExcept for diagram %d which is a diagram of S1 X S2 or S1 X D2 and is empty.     ",i);
                }
            else
                {
                j = 0;
                j += printf("\n\nExcept for diagrams: ");
                for(h = 0,k = 1; h < NumFilled; h++) if(SURL[h] == 0)
                    {
                    h++;
                    j += printf("{%d,",h);    
                    break;
                    }
                for( ; h < NumFilled; h++) if(SURL[h] == 0)
                    {
                    if(++k < i)
                        j += printf("%d,",h+1);
                    else
                        j += printf("%d}.",h+1);    
                    if(j > 80)
                        {
                        j = 0;
                        printf("\n");
                        }
                    }
                printf("\nThese are diagrams of S1 X S2 (s) or S1 X D2 (s) and are empty.      ");
                }
            }
        GET_RESPONSE3:
        WhichInput = 0;                
        ReadString((char *)ptr, GetPtrSize(ptr));
        sscanf((char *) ptr,"%d",&WhichInput);        
        if(WhichInput < 1 || WhichInput > NumFilled || SURL[WhichInput-1] == 0)
            {
            if(Batch == FALSE) SysBeep(5);
            goto GET_RESPONSE3;
            }    
        DisposePtr((char *) ptr);
        WhichInput --;
                
        if(Get_Relators_From_SUR(WhichInput))
            {
            printf("\n\n    Memory Error. Sorry!");
            goto REDRAW;        
            }
                
        if(Length == 0)
            {
            printf("\n\nThis is an empty presentation. There is nothing to display.");
            goto REDRAW;
            }
        
        if(WhichInput == 0 && UDV[0] == SPLIT)
            {
            switch(Find_Flow_A(NORMAL,FALSE))
                {
                case 1:
                    break;
                case TOO_LONG:
                    if(Batch == FALSE) SysBeep(5);
                    printf("\n\n     Unable to display data for diagram %d. Sorry!",WhichInput + 1);
                    goto REDRAW;
                }
            if(Automorphisms)
                {
                if(Batch == FALSE) SysBeep(5);
                printf("\n\n                    NOTE!");
                printf("\n\n    Presentation 1 does not have minimal length.");
                printf("\n    Heegaard will only display data for diagrams of minimal length presentations.");
                printf("\n    Presentation 2 should give a diagram of a minimal length version of presentation 1.");
                printf("\n\n    HIT ANY KEY TO CONTINUE.");
	            WaitkbHit();   
                goto REDRAW;
                }    
            }
        
        Fill_A(NumRelators);
        if(ComputeValences_A())
            {
            if(Batch == FALSE) SysBeep(5);
            printf("\n\n     Unable to display data for diagram %d. Sorry!",WhichInput + 1);
            goto REDRAW;
            }    
        Get_Matrix();
        Check_Connected();
        SepPairs = Sep_Pairs(0,0,1);
        if(SepPairs == TOO_LONG)
            {
            if(Batch == FALSE) SysBeep(5);
            printf("\n\n     Unable to display data for diagram %d. Sorry!",WhichInput + 1);
            goto REDRAW;
            }
        SaveUDV = UDV[WhichInput];
        SLSP = LSP[WhichInput];
        SLSQ = LSQ[WhichInput];
        SWhichInput = WhichInput;    
        if(SepPairs)
            {
            if(UDV[WhichInput] == 0) UDV[WhichInput] = SEP_PAIRS;
            if(V1 & 1)
                LSP[WhichInput] = V1/2 + 97;
            else
                LSP[WhichInput] = V1/2 + 65;
            if(V2 & 1)
                LSQ[WhichInput] = V2/2 + 97;
            else
                LSQ[WhichInput] = V2/2 + 65;
            }
        NonPlanar = Planar(FALSE,TRUE);
        Print_Graph(FALSE,0,0,0);
        if(UDV[SWhichInput] != ANNULUS_EXISTS && UDV[SWhichInput] != V2_ANNULUS_EXISTS)
            UDV[SWhichInput] = SaveUDV;
        LSP[SWhichInput] = SLSP;
        LSQ[SWhichInput] = SLSQ;        
        }
    if(Response == 'a' || Response == 'b')
        {
        NoSepPairs = FALSE;
        NumConnected = 0;
        for(WhichInput = 0; WhichInput < NumFilled; WhichInput ++)
        if(SURL[WhichInput] != 0)
            {
            h = UDV[WhichInput];
            if(Response == 'b' && (h == SEP_PAIRS || h == ANNULUS_EXISTS)) continue;        

            if(Get_Relators_From_SUR(WhichInput))
                {
                printf("\n\n    Memory Error. Sorry!");
                goto REDRAW;        
                }
            
            if(WhichInput == 0 && h == SPLIT)
                {
                switch(Find_Flow_A(NORMAL,FALSE))
                    {
                    case 1:
                        break;
                    case TOO_LONG:
                        if(Batch == FALSE) SysBeep(5);
                        printf("\n\n     Unable to display data for diagram %d. Sorry!",WhichInput + 1);
                        continue;
                    }
                if(Automorphisms)
                    {
                    if(Batch == FALSE) SysBeep(5);
                    printf("\n\n                    NOTE!");
                    printf("\n\n    Presentation 1 does not have minimal length.");
                    printf("\n    Heegaard will only display data for diagrams of minimal length presentations.");
                    printf("\n    Presentation 2 should give a diagram of a minimal length version of presentation 1.");
                    printf("\n\nHIT ANY KEY TO SEE DATA FOR THE NEXT DIAGRAM.");
	    			WaitkbHit();
                    continue;
                    }    
                }
            
            Fill_A(NumRelators);
            if(ComputeValences_A())
                {
                if(Batch == FALSE) SysBeep(5);
                printf("\n\n     Unable to display data for diagram %d. Sorry!",WhichInput + 1);
                continue;
                }            
            Get_Matrix();
            Check_Connected();
            if(!Connected && Response == 'b') continue;
            NumConnected ++;
            SepPairs = Sep_Pairs(0,0,1);
            if(SepPairs == TOO_LONG)
                {
                if(Batch == FALSE) SysBeep(5);
                printf("\n\n     Unable to display data for diagram %d. Sorry!",WhichInput + 1);
                continue;
                }    
            SaveUDV = UDV[WhichInput];
            SLSP = LSP[WhichInput];
            SLSQ = LSQ[WhichInput];
            if(SepPairs)
                {
                if(Response == 'b') continue;                
                if(UDV[WhichInput] == 0) UDV[WhichInput] = SEP_PAIRS;
                if(V1 & 1)
                    LSP[WhichInput] = V1/2 + 97;
                else
                    LSP[WhichInput] = V1/2 + 65;
                if(V2 & 1)
                    LSQ[WhichInput] = V2/2 + 97;
                else
                    LSQ[WhichInput] = V2/2 + 65;        
                }
            NoSepPairs = TRUE;    
            NonPlanar = Planar(FALSE,TRUE);
            SWhichInput = WhichInput;
            if(Print_Graph(TRUE,0,0,0))
                {
                if(WhichInput > NumFilled) 
                	{
                	DrawingDiagrams = FALSE;
                	return;
                	}
                }
            if(UDV[SWhichInput] != ANNULUS_EXISTS && UDV[SWhichInput] != V2_ANNULUS_EXISTS)
                UDV[SWhichInput] = SaveUDV;
            LSP[SWhichInput] = SLSP;
            LSQ[SWhichInput] = SLSQ;                    
            }
        }

    if(Response == 'b' && (!NoSepPairs || !NumConnected))
        {
        if(Batch == FALSE) SysBeep(5);
        printf("\n\n    There are no diagrams which are connected and without separating pairs of vertices!");
        }

    DrawingDiagrams = FALSE;                   
}

int Batch_Report(int * Table)
{	
    unsigned char   Finite,
    				*p,
    				PrintPres,
    				PrintedS3,
                    x,
                    y;
    
    int             *CompType1 = NULL,
    				*CompType2 = NULL,
    				MyTotalCompFound,
    				MyTotalFiniteComp,
    				NumRelators,
    				NumSFFound,
    				Start;
                     
    unsigned int    i,
                    j,
                    k,
                    m,
                    n;                                
    
    unsigned long   Length;

	CompType1 = (int *)NewPtr(sizeof(int)*(TotalComp + 1));
	if(CompType1 == NULL) Mem_Error();
	CompType2 = (int *)NewPtr(sizeof(int)*(TotalComp + 1));
	if(CompType2 == NULL) Mem_Error();	
	for(k = 1; k <= TotalComp; k++) CompType1[k] = CompType2[k] = 0;
			    
	for(n = 0; n < NumFilled; n++)
        {
        i = Table[n];
        PrintPres = FALSE;
		switch(PRIM[i])
			{
			case 3:
			case 103:
				/* printf("LS"): */
				PrintPres = TRUE;
				break;
			case 4:
			case 104:
				/* printf("1G"): */
				PrintPres = TRUE;
				break;
			case 5:
			case 105:
				/* printf("S3"): */
				PrintPres = TRUE;
				break;
			case 12:
			case 112:
				/* printf("ER"): A presentation corresponding to the union of empty summands and or S1 X S2s. */
				PrintPres = TRUE;
				break;
			case 13:
			case 113:
				/* printf("Er"): A bandsum created an empty relator. */
				PrintPres = TRUE;
				break;                
			case 20:
			case 120:
				/* printf("NC"): */
				if((Batch == 10 || Batch == 11) && BPrintNotConnectedData == TRUE) PrintPres = TRUE;
				PrintPres = TRUE;
				break;
			case 30:
			case 130:
			case 40:
			case 140:
				/* printf("MG"): */
				PrintPres = TRUE;
				break;
			case 60:
			case 160:
				/* printf("PP"): */
				PrintPres = TRUE;
				break;
			case 80:
			case 180:
				/* printf("A2"): */
				if((Batch == 10 || Batch == 11) && BPrintAnnulusData == TRUE) PrintPres = TRUE;
				break;    
			default:
				break;
			}
		
		if(PrintPres == FALSE) switch(UDV[i])
			{
			case SPLIT:
				PrintPres = TRUE;
				break;
			case GENERIC_LENS_SPACE:
				PrintPres = TRUE;
				break;	
			case THREE_SPHERE:
				PrintPres = TRUE;
				break;
			case NOT_CONNECTED:
				if((Batch == 10 || Batch == 11) && BPrintNotConnectedData == TRUE) PrintPres = TRUE;
				break;			
			case S1_X_S2:
				PrintPres = TRUE;
				break;	
			case S1_X_D2:
				PrintPres = TRUE;
				break;	
			case S1_X_X2:
				PrintPres = TRUE;
				break;	
			case MISSING_GEN_DONE2:
				PrintPres = TRUE;
				break;		
			case MISSING_GEN_DONE1:
				PrintPres = TRUE;
				break;						
			case KNOWN_LENS_SPACE:
				PrintPres = TRUE;
				break;	
			case ANNULUS_EXISTS:
				if((Batch == 10 || Batch == 11) && BPrintAnnulusData == TRUE) PrintPres = TRUE;
				break;	
			case V2_ANNULUS_EXISTS:
				if((Batch == 10 || Batch == 11) && BPrintAnnulusData == TRUE) PrintPres = TRUE;
				break;	
			case NON_UNIQUE_4:
				PrintPres = TRUE;
				break;	
			case NON_UNIQUE_3:
				PrintPres = TRUE;
				break;	
			case NON_UNIQUE_2:
				PrintPres = TRUE;
				break;	
			case NON_UNIQUE_1:
				PrintPres = TRUE;
				break;		
			default:
				break;                                                                                      
			}

		if(PrintPres)
			{
			NumRelators = NR[i];
			Length 		= SURL[i];
			printf("\n\nPresentation %d  of Summand %u:  Gen  %d  Rel %d  Length  %lu  From Pres %u ",i+1,ComponentNum[i],NG[i],NR[i],Length,FR[i]+1);	
			switch(PRIM[i])
				{
				case 3:
				case 103:
					printf("LS");
					break;
				case 4:
				case 104:
					printf("1G");
					break;
				case 5:
				case 105:
					printf("S3");
					break;
				case 12:
				case 112:
					printf("ER"); /* A presentation corresponding to the union of empty summands and or S1 X S2s. */
					break;
				case 13:
				case 113:
					printf("Er"); /* A bandsum created an empty relator. */
					break;                
				case 20:
				case 120:
					printf("NC");
					break;
				case 30:
				case 130:
				case 40:
				case 140:
					printf("MG");
					break;
				case 60:
				case 160:
					printf("PP");
					break;
				case 80:
				case 180:
					printf("A2");
					break;    
				default:
					break;
				}
								
			switch(UDV[i])
				{
				case SPLIT:
					CompType1[ComponentNum[i]] = SPLIT;
					CompType2[ComponentNum[i]] = i + 1;
					j = Daughters[i];
					m = j + NCS[i] - 1;
					switch(PRIM[j])
						{
						case 40:
						case 140:
							printf("\nThe presentation dual to presentation %d 'split' into presentations %u and %u of M",i + 1,j + 1,m + 1);
							printf("\nwhich correspond to summands %u and %u of M.",ComponentNum[j],ComponentNum[m]);                    
							break;
						default:
							printf("\n'split' into presentations %u and %u corresponding to summands %u and %u of M.",j + 1,m + 1,ComponentNum[j],ComponentNum[m]);
							break;                    
						}
					if(BdryData) for(k = 0; k < NCS[i]; k++) Print_Bdry_Data(j + k);						
					if(B10B11ConSumPres && (Batch == 10 || Batch == 11) && H_Results != NULL)
						{
						fprintf(H_Results,"\n\nPresentation %d of Summand %u of %s:  Gen  %d  Rel %d  Length  %lu  From Pres %u",i+1,ComponentNum[i],PresName,NG[i],NR[i],SURL[i],FR[i]+1);						
						switch(PRIM[j])
							{
							case 40:
							case 140:
								fprintf(H_Results," MG");
								fprintf(H_Results,"\nThe presentation dual to presentation %d 'split' into presentations %u and %u of M",i + 1,j + 1,m + 1);
								fprintf(H_Results,"\nwhich correspond to summands %u and %u of M.",ComponentNum[j],ComponentNum[m]);                    
								break;
							default:
								fprintf(H_Results,"\n'split' into presentations %u and %u corresponding to summands %u and %u of M.",j + 1,m + 1,ComponentNum[j],ComponentNum[m]);
								break;                    
							}
						if(BdryData) for(k = 0; k < NCS[i]; k++) Print_Bdry_Data2(j + k);	
						for(k = 1; k <= NR[i]; k++) fprintf(H_Results,"\n    %s",*SUR[i][k]);
						fprintf(H_Results,"\n\nPresentation %d  of Summand %u:  Gen  %d  Rel %d  Length  %lu  From Pres %u",j+1,ComponentNum[j],NG[j],NR[j],SURL[j],FR[j]+1);						
						for(k = 1; k <= NR[j]; k++) fprintf(H_Results,"\n    %s",*SUR[j][k]);						
						fprintf(H_Results,"\n\nPresentation %d  of Summand %u:  Gen  %d  Rel %d  Length  %lu  From Pres %u",m+1,ComponentNum[m],NG[m],NR[m],SURL[m],FR[m]+1);						
						for(k = 1; k <= NR[m]; k++) fprintf(H_Results,"\n    %s",*SUR[m][k]);						
						}
					break;
			
				case GENERIC_LENS_SPACE:
					if(CompType1[ComponentNum[i]] == 0)
						{
						CompType1[ComponentNum[i]] = GENERIC_LENS_SPACE;
						CompType2[ComponentNum[i]] = i + 1;
						}
					if(LSP[i] > 4L)
						printf("\npresents a Lens space of the form L(%lu,Q).",LSP[i]);
					else
						{
						if(LSP[i] == 1L) printf("\npresents the 3-Sphere.");
						else printf("\npresents a Lens space of the form L(%lu,1).",LSP[i]);
						}
					break;        
			
				case THREE_SPHERE:
					CompType1[ComponentNum[i]] = THREE_SPHERE;
					CompType2[ComponentNum[i]] = i + 1;
					printf("\npresents the 3-Sphere.");
					break;

				case NOT_CONNECTED:
					CompType1[ComponentNum[i]] = NOT_CONNECTED;
					CompType2[ComponentNum[i]] = i + 1;
					printf("\nthe diagram is not connected.");
					if(B10B11ConSumPres && (Batch == 10 || Batch == 11) && H_Results != NULL)
						{
						fprintf(H_Results,"\n\nPresentation %d of Summand %u of %s:  Gen  %d  Rel %d  Length  %lu  From Pres %u",i+1,ComponentNum[i],PresName,NG[i],NR[i],SURL[i],FR[i]+1);						
						if(BdryData) Print_Bdry_Data2(i);
						for(k = 1; k <= NR[i]; k++) fprintf(H_Results,"\n    %s",*SUR[i][k]);
						}
					break;
					
				case S1_X_S2:
					CompType1[ComponentNum[i]] = S1_X_S2;
					CompType2[ComponentNum[i]] = i + 1;
					if(NG[i] == 1)
						printf("\npresents S1 X S2.");
					else
						printf("\npresents %d copies of S1 X S2.",NG[i]);
					break;

				case S1_X_D2:
					CompType1[ComponentNum[i]] = S1_X_D2;
					CompType2[ComponentNum[i]] = i + 1;
					if(NG[i] == 1)
						printf("\npresents S1 X D2.");
					else
						printf("\npresents %d copies of S1 X D2.",NG[i]);
					break;
			
				case S1_X_X2:
					CompType1[ComponentNum[i]] = S1_X_X2;
					CompType2[ComponentNum[i]] = i + 1;
					if(NG[i] == 1)
						printf("\npresents S1 X ?2.");
					else
						printf("\npresents %d copies of S1 X ?2.",NG[i]);
					break;
			
				case MISSING_GEN_DONE2:
					CompType1[ComponentNum[i]] = MISSING_GEN_DONE2;
					CompType2[ComponentNum[i]] = i + 1;
					printf("\npresents: ");
					if(N1H[ComponentNum[i]] == 1)
						printf("1 copy of I X D2, ");
					else
						printf("%d copies of I X D2, ",N1H[ComponentNum[i]]);
					if(NS1XS2[ComponentNum[i]] == 1)
						printf("1 copy of S1 X S2, ");
					else
						printf("%d copies of S1 X S2, ",NS1XS2[ComponentNum[i]]);
					if(NS1XD2[ComponentNum[i]] == 1)
						printf("and 1 copy of S1 X D2.");
					else
						printf("and %d copies of S1 X D2.",NS1XD2[ComponentNum[i]]);                                        
					printf("\nHeegaard was not able to unambiguously determine the boundary components.");
					break;
				
				case MISSING_GEN_DONE1:
					CompType1[ComponentNum[i]] = MISSING_GEN_DONE1;
					CompType2[ComponentNum[i]] = i + 1;
					printf("\npresents: ");
					if(N1H[ComponentNum[i]] == 1)
						printf("1 copy of I X D2, ");
					else
						printf("%d copies of I X D2, ",N1H[ComponentNum[i]]);
					if(NS1XS2[ComponentNum[i]] == 1)
						printf("1 copy of S1 X S2, ");
					else
						printf("%d copies of S1 X S2, ",NS1XS2[ComponentNum[i]]);
					if(NS1XD2[ComponentNum[i]] == 1)
						printf("and 1 copy of S1 X D2.");
					else
						printf("and %d copies of S1 X D2.",NS1XD2[ComponentNum[i]]);                            
					break;
				
				case KNOWN_LENS_SPACE:
					CompType1[ComponentNum[i]] = KNOWN_LENS_SPACE;
					CompType2[ComponentNum[i]] = i + 1;
					switch(LSP[i])
						{
						case 0L:
							printf("\npresents S1 X S2.");
							break; 
						case 1L:
							printf("\npresents the 3-Sphere.");
							break;
						default:
							printf("\npresents the Lens space L(%lu,%lu).",LSP[i],LSQ[i]);
							break;
						}
					break;
		
				case ANNULUS_EXISTS:
					p = *SUR[i][0];
					x = *p++;
					y = *p++;                
					printf("\nThe pair of vertices (%c,%c) separate the diagram.",x,y);
					printf("\nThe component consisting of vertice(s) {%c",*p);
					p++;
					while((x = *p) != '@')
						{
						printf(",%c",x);
						p++;
						}
					printf("}");
					p++;        
					printf("\nlies in an annulus which swallows the component and otherwise follows the curve:\n");
					while(*p)
						fputc(*p++,stdout);
					break;
	
				case V2_ANNULUS_EXISTS:
					p = *SUR[i][0];
					printf("\nThere exists an annulus which swallows vertice(s) {%c",*p);
					p++;
					while((x = *p) != '@')
						{
						printf(",%c",x);
						p++;
						}
					printf("}"); 
					p++;        
					printf("\nand otherwise follows the curve:\n");
					while(*p) fputc(*p++,stdout);
					break;
			
				case NON_UNIQUE_4:
					printf("\nthe diagram is not unique because there is a generator which appears");
					printf("\nwith only one exponent and that exponent is greater than 6.");
					break;
			
				case NON_UNIQUE_3:
					printf("\nthe diagram is not unique because there is a generator which appears");
					printf("\nonly with exponent 5.");
					break;
			
				case NON_UNIQUE_2:
					printf("\nthe diagram is not unique because there is a generator which appears");
					printf("\nonly with exponent 3 or only with exponent 4 and this exponent occurs more than once.");
					break;
			
				case NON_UNIQUE_1:
					printf("\nthe diagram is not unique because there is a generator which appears");
					printf("\nwith only one exponent, either 3,4 or 6, and a needed symmetry does not exist.");
					break; 
					
				case DUPLICATE:					
					if(CompType1[ComponentNum[i]] == 0 && ComponentNum[Daughters[i]] < ComponentNum[i])
						{
						CompType1[ComponentNum[i]] = CompType1[ComponentNum[Daughters[i]]];
						CompType2[ComponentNum[i]] = Daughters[i] + 1;
						}		
					printf("\nis a duplicate of presentation %d of summand %u.",
					Daughters[i] + 1,ComponentNum[Daughters[i]]);
					break;					                                                                                      	
				}
		    printf("\n");
        	for(j = 1; j <= NumRelators; j++) printf("\n    %s",*SUR[i][j]);
        	}						
		}
	
	if(B10B11Finite && (Batch == 10 || Batch == 11) && H_Results != NULL)
		{	
		for(k = 1,Start = NumFilled - 1,MyTotalFiniteComp = 0,Finite = TRUE; k <= TotalComp; k++)
			{
			switch(CompType1[k])
				{
				case 0:
					if(SFSols[k] != NULL) 
						{
						DisposePtr((unsigned int*) SFSols[k]);
						SFSols[k] = NULL;
						}
					FoundFiniteSF = FALSE;	
					Start = Init_Genus_Two_Seifert_Fibered(Table,Start,k);
					if(FoundFiniteSF == TRUE) switch(SFSols[k][0])
						{
						case 2:
						case 4:
							if(SFSols[k][5] > 1) MyTotalFiniteComp ++;
							break;
						case 8:
						case 9:
						case 10:
						case 11:
							if(SFSols[k][5] > 1) MyTotalFiniteComp ++;
							if(SFSols[k][7] > 1) MyTotalFiniteComp ++;
							break;
						case 14:
							if(SFSols[k][18] > 1) MyTotalFiniteComp ++;
							break;
						case 16:
						case 18:
							MyTotalFiniteComp ++;	
						}
					else Finite = FALSE;
					break;
				case GENERIC_LENS_SPACE:
					i = CompType2[k] - 1;
					if(LSP[i] == 0) Finite = FALSE;
					if(LSP[i] > 1) MyTotalFiniteComp ++;
					break;        
		
				case THREE_SPHERE:
					break;	
				
				case S1_X_S2:
				case S1_X_D2:	
				case S1_X_X2:
				case MISSING_GEN_DONE2:			
				case MISSING_GEN_DONE1:
					Finite = FALSE;
					break;
			
				case KNOWN_LENS_SPACE:
					i = CompType2[k] - 1;
					if(LSP[i] == 0) Finite = FALSE;
					if(LSP[i] > 1) MyTotalFiniteComp ++;
					break;										
				}
			if(MyTotalFiniteComp > 1) Finite = FALSE;
			if(Finite == FALSE) break;	
			}
			
		if(Finite)  
			{
			fprintf(H_Results,"\n\n%-20s ",PresName);
			for(k = 1,Start = NumFilled - 1,PrintedS3 = FALSE; k <= TotalComp; k++)
				{
				i = CompType2[k] - 1;
				switch(CompType1[k])
					{
					case 0: 
						if(SFSols[k] != NULL)
							Print_SFComp(k);
						else
							fprintf(H_Results,"? ");		
						break;
					case GENERIC_LENS_SPACE:
						if(LSP[i] == 0) 	fprintf(H_Results,"S1 X S2 ");
						if(LSP[i] == 1 && MyTotalFiniteComp == 0 && PrintedS3 == FALSE) 	
							{
							fprintf(H_Results,"S^3 ");
							PrintedS3 = TRUE;
							break;
							}
						if(1 < LSP[i] && LSP[i] < 5) 	fprintf(H_Results,"L(%lu,1) ",LSP[i]);
						if(LSP[i] > 4) 		fprintf(H_Results,"L(%lu,Q) ",LSP[i]);
						break;			
					case THREE_SPHERE: 
						if(MyTotalFiniteComp == 0 && PrintedS3 == FALSE) 
							{
							fprintf(H_Results,"S^3 ");
							PrintedS3 = TRUE;
							}
						break;
					case KNOWN_LENS_SPACE:
						if(LSP[i] == 0) fprintf(H_Results,"S1 X S2 ");
						if(LSP[i] == 1 && MyTotalFiniteComp == 0 && PrintedS3 == FALSE) 
							{
							fprintf(H_Results,"S^3 ");
							PrintedS3 = TRUE;
							}
						if(LSP[i] > 1)	fprintf(H_Results,"L(%lu,%lu) ",LSP[i],LSQ[i]);
						break;				
					case S1_X_S2:
						if(NG[i] == 1) 	fprintf(H_Results,"S1 X S2 ");
						if(NG[i] > 1) 	fprintf(H_Results,"%d S1 X S2s ",NG[i]);
						break;			
					case S1_X_D2:
						if(NG[i] == 1) 	fprintf(H_Results,"S1 X D2 ");
						if(NG[i] > 1) 	fprintf(H_Results,"%d S1 X D2s ",NG[i]);
						break;							
					case S1_X_X2:
						if(NG[i] == 1) 	fprintf(H_Results,"S1 X X2 ");
						if(NG[i] > 1) 	fprintf(H_Results,"%d S1 X X2s ",NG[i]);
						break;							
					case MISSING_GEN_DONE2:
						if(N1H[k] == 1)		fprintf(H_Results,"I X D2 ");
						if(N1H[k] > 1)		fprintf(H_Results,"%d I X D2s ",N1H[k]);
						if(NS1XS2[k] == 1) 	fprintf(H_Results,"S1 X S2 ");
						if(NS1XS2[k] > 1) 	fprintf(H_Results,"%d S1 X S2s ",NS1XS2[k]);
						if(NS1XD2[k] == 1)	fprintf(H_Results,"S1 X D2 ");
						if(NS1XD2[k] > 1) 	fprintf(H_Results,"%d S1 X D2s ",NS1XD2[k]);                                        
						fprintf(H_Results,"(Heegaard was not able to unambiguously determine the boundary components.) ");
						break;			
					case MISSING_GEN_DONE1:
						if(N1H[k] == 1)		fprintf(H_Results,"I X D2 ");
						if(N1H[k] > 1)		fprintf(H_Results,"%d I X D2s ",N1H[k]);
						if(NS1XS2[k] == 1) 	fprintf(H_Results,"S1 X S2 ");
						if(NS1XS2[k] > 1) 	fprintf(H_Results,"%d S1 X S2s ",NS1XS2[k]);
						if(NS1XD2[k] == 1)	fprintf(H_Results,"S1 X D2 ");
						if(NS1XD2[k] > 1) 	fprintf(H_Results,"%d S1 X D2s ",NS1XD2[k]);
						break;
					case SPLIT:
						break;
					case NOT_CONNECTED:
						break;	
					default:
						break;										
					}	
				}	
			}
		}
		
	if(B10B11Recognized && (Batch == 10 || Batch == 11) && H_Results != NULL)
		{
		for(k = 1,Start = NumFilled - 1,NumSFFound = 0,MyTotalCompFound = 0; k <= TotalComp; k++) 
			switch(CompType1[k])
				{
				case 0:
					if(SFSols[k] != NULL) 
						{
						DisposePtr((unsigned int*) SFSols[k]);
						SFSols[k] = NULL;
						}
					FoundSF = FALSE;	
					Start = Init_Genus_Two_Seifert_Fibered(Table,Start,k);
					if(FoundSF) NumSFFound ++;
					MyTotalCompFound ++;
					break;
				case GENERIC_LENS_SPACE:
				case THREE_SPHERE:
				case KNOWN_LENS_SPACE:
				case S1_X_S2:
				case S1_X_D2:
				case S1_X_X2:
				case MISSING_GEN_DONE2:
				case MISSING_GEN_DONE1:
					MyTotalCompFound ++;
					break;
				}	
		}
		
		
	if(B10B11Recognized && (Batch == 10 || Batch == 11) && H_Results != NULL)  
		{
		fprintf(H_Results,"\n\n%-20s ",PresName);
		if(MyTotalCompFound == 0) fprintf(H_Results,"?");
		else for(k = 1,Start = NumFilled - 1,m = n = 0; k <= TotalComp; k++)
			{
			i = CompType2[k] - 1;
			switch(CompType1[k])
				{
				case 0: 
					if(SFSols[k] != NULL)
						Print_SFComp(k);
					else
						fprintf(H_Results,"? ");
					m++;		
					break;
				case GENERIC_LENS_SPACE:
					if(LSP[i] == 0) 	fprintf(H_Results,"S1 X S2 ");
					if(LSP[i] == 1) 	fprintf(H_Results,"S^3 ");
					if(1 < LSP[i] && LSP[i] < 5) 	fprintf(H_Results,"L(%lu,1) ",LSP[i]);
					if(LSP[i] > 4) 		fprintf(H_Results,"L(%lu,Q) ",LSP[i]);
					m++;
					break;			
				case THREE_SPHERE: fprintf(H_Results,"S^3 ");
					m++;
					break;
				case KNOWN_LENS_SPACE:
					if(LSP[i] == 0) fprintf(H_Results,"S1 X S2 ");
					if(LSP[i] == 1) fprintf(H_Results,"S^3 ");
					if(LSP[i] > 1)	fprintf(H_Results,"L(%lu,%lu) ",LSP[i],LSQ[i]);
					m++;
					break;				
				case S1_X_S2:
					if(NG[i] == 1) 	fprintf(H_Results,"S1 X S2 ");
					if(NG[i] > 1) 	fprintf(H_Results,"%d S1 X S2s ",NG[i]);
					m++;
					break;			
				case S1_X_D2:
					if(NG[i] == 1) 	fprintf(H_Results,"S1 X D2 ");
					if(NG[i] > 1) 	fprintf(H_Results,"%d S1 X D2s ",NG[i]);
					m++;
					break;							
				case S1_X_X2:
					if(NG[i] == 1) 	fprintf(H_Results,"S1 X X2 ");
					if(NG[i] > 1) 	fprintf(H_Results,"%d S1 X X2s ",NG[i]);
					m++;
					break;							
				case MISSING_GEN_DONE2:
					if(N1H[k] == 1)		fprintf(H_Results,"I X D2 ");
					if(N1H[k] > 1)		fprintf(H_Results,"%d I X D2s ",N1H[k]);
					if(NS1XS2[k] == 1) 	fprintf(H_Results,"S1 X S2 ");
					if(NS1XS2[k] > 1) 	fprintf(H_Results,"%d S1 X S2s ",NS1XS2[k]);
					if(NS1XD2[k] == 1)	fprintf(H_Results,"S1 X D2 ");
					if(NS1XD2[k] > 1) 	fprintf(H_Results,"%d S1 X D2s ",NS1XD2[k]);                                        
					fprintf(H_Results,"(Heegaard was not able to unambiguously determine the boundary components.) ");
					m++;
					break;			
				case MISSING_GEN_DONE1:
					if(N1H[k] == 1)		fprintf(H_Results,"I X D2 ");
					if(N1H[k] > 1)		fprintf(H_Results,"%d I X D2s ",N1H[k]);
					if(NS1XS2[k] == 1) 	fprintf(H_Results,"S1 X S2 ");
					if(NS1XS2[k] > 1) 	fprintf(H_Results,"%d S1 X S2s ",NS1XS2[k]);
					if(NS1XD2[k] == 1)	fprintf(H_Results,"S1 X D2 ");
					if(NS1XD2[k] > 1) 	fprintf(H_Results,"%d S1 X D2s ",NS1XD2[k]);
					m++;
					break;
				case SPLIT:
					fprintf(H_Results,"Split: ");
					break;
				case NOT_CONNECTED:
					fprintf(H_Results,"NC: ");
					break;	
				default:
					break;										
				}
			if(m > n && m < MyTotalCompFound) 
				{
				fprintf(H_Results,"# ");
				n = m;
				}	
			}	
		}
					
	printf("\n");
	
	if(CompType1 != NULL) DisposePtr((int *) CompType1);
	if(CompType2 != NULL) DisposePtr((int *) CompType2);	
	return(0);		
}

int Print_SFComp(int MyComp)
{
	int	A1,
		A2,
		A3,
		a1,
		a2,
		a3,
		B1,
		B2,
		B3,
		b1,
		b2,
		b3,
		H1,
		n,
		m,
		Q;
	
	B1 = SFSols[MyComp][4];		
	A1 = SFSols[MyComp][5];
	B2 = SFSols[MyComp][6];
	A2 = SFSols[MyComp][7];
	B3 = SFSols[MyComp][8];
	A3 = SFSols[MyComp][9];
	n  = SFSols[MyComp][10];
	b1 = SFSols[MyComp][11];
	a1 = SFSols[MyComp][12];
	b2 = SFSols[MyComp][13];
	a2 = SFSols[MyComp][14];
	b3 = SFSols[MyComp][15];
	a3 = SFSols[MyComp][16];
	m  = SFSols[MyComp][17];
	H1 = SFSols[MyComp][18];
	Q =  SFSols[MyComp][19];
	
	switch(SFSols[MyComp][0])
		{
		case 1: 
			fprintf(H_Results,"S^1 X S^2 ");
			break;
		case 2: 
			fprintf(H_Results,"L(%d,%d) ",A1,B1);
			break;
		case 3: 
			fprintf(H_Results,"S^1 X S^2 ");
			break;
		case 4: 
			fprintf(H_Results,"L(%d,%d) ",A1,B1);
			break;
		case 5: 
			fprintf(H_Results,"SF over the Mobius band ");
			break;
		case 6: 
			fprintf(H_Results,"SF over RP^2 ");
			break;
		case 7: 
			fprintf(H_Results,"SF(0;m/%d,n/%d), 0 < m < %d, 0 < n < %d, gcd(m,%d) = gcd(n,%d) = 1 ",A1,A2,A1,A2,A1,A2);	
			break;
		case 8: 
			fprintf(H_Results,"L(%d,Q1), L(%d,Q2) ",A1,A2);
			break;
		case 9: 
			fprintf(H_Results,"L(%d,Q), L(%d,%d) ",A1,A2,B2);
			break;
		case 10: 
			fprintf(H_Results,"L(%d,%d), L(%d,Q) ",A1,B1,A2);
			break;
		case 11: 
			fprintf(H_Results,"L(%d,%d), L(%d,%d) ",A1,B1,A2,B2);
			break;
		case 12: 
			fprintf(H_Results,"S^1 X S^2 ");
			break;
		case 13: 
			fprintf(H_Results,"S^1 X S^2 ");
			break;
		case 14: 
			fprintf(H_Results,"L(%d,%d) ",H1,Q);
			break;
		case 15: 
			fprintf(H_Results,"SF(0;%d;%d/%d,%d/%d,%d/%d) or OR",n,B1,A1,B2,A2,B3,A3);
			fprintf(H_Results," SF(0;%d;%d/%d,%d/%d,%d/%d) ",3-n,A1-B1,A1,A2-B2,A2,A3-B3,A3);
			break;
		case 16: 
			fprintf(H_Results,"SF(0;%d;%d/%d,%d/%d,%d/%d) or OR" ,n,B1,A1,B2,A2,B3,A3);
			fprintf(H_Results," SF(0;%d;%d/%d,%d/%d,%d/%d) ",3-n,A1-B1,A1,A2-B2,A2,A3-B3,A3);
			break;
		case 17: 
			fprintf(H_Results,"SF(0;%d;%d/%d,%d/%d,%d/%d) or OR",n,B1,A1,B2,A2,B3,A3);
			fprintf(H_Results," SF(0;%d;%d/%d,%d/%d,%d/%d) ",3-n,A1-B1,A1,A2-B2,A2,A3-B3,A3); 
			fprintf(H_Results,"or perhaps:");			
			fprintf(H_Results,"\n                     SF(0;%d;%d/%d,%d/%d,%d/%d) or OR",m,b1,a1,b2,a2,b3,a3);
			fprintf(H_Results," SF(0;%d;%d/%d,%d/%d,%d/%d) ",3-m,a1-b1,a1,a2-b2,a2,a3-b3,a3);
			break;	
		case 18:
			fprintf(H_Results,"SF(0;%d;%d/%d,%d/%d,%d/%d) or OR",n,B1,A1,B2,A2,B3,A3);
			fprintf(H_Results," SF(0;%d;%d/%d,%d/%d,%d/%d) ",3-n,A1-B1,A1,A2-B2,A2,A3-B3,A3); 
			fprintf(H_Results,"or perhaps:");			
			fprintf(H_Results,"\n                     SF(0;%d;%d/%d,%d/%d,%d/%d) or OR",m,b1,a1,b2,a2,b3,a3);
			fprintf(H_Results," SF(0;%d;%d/%d,%d/%d,%d/%d) ",3-m,a1-b1,a1,a2-b2,a2,a3-b3,a3);
			break;	
		}
		
	return(0);	
}
