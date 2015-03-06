#include "Heegaard.h"
#include "Heegaard_Dec.h"

/****************************** function prototypes *****************************************
L   8 Do_Initialization(void)
********************************************************************************************/

int Do_Initialization()
{
	register unsigned int 	i,
							j;
	
	unsigned long 			Seconds;
	
	for(i = 0; i < VERTICES; i++)
		{
		A[i] 		= (unsigned int *)	NewPtr(sizeof(int)*VERTICES);
		if(A[i] 	== NULL) Mem_Error();
		AA[i] 		= (unsigned int *)	NewPtr(sizeof(int)*VERTICES);
		if(AA[i] 	== NULL) Mem_Error();
		AJ1[i] 		= (unsigned int *)	NewPtr(sizeof(int)*(VERTICES + 1));
		if(AJ1[i] 	== NULL) Mem_Error();
		AJ2[i] 		= (unsigned int *)	NewPtr(sizeof(int)*(VERTICES + 1));
		if(AJ2[i] 	== NULL) Mem_Error();
		AJ3[i] 		= (unsigned int *)	NewPtr(sizeof(int)*(VERTICES + 1));
		if(AJ3[i] 	== NULL) Mem_Error();
		B[i] 		= (unsigned int *)	NewPtr(sizeof(int)*VERTICES);
		if(B[i] 	== NULL) Mem_Error();
		CO[i]		= (unsigned char*)	NewPtr(sizeof(char)*VERTICES);
		if(CO[i] 	== NULL) Mem_Error();
		GB[i]		= (         int *)  NewPtr(sizeof(int)*VERTICES);
		if(GB[i] 	== NULL) Mem_Error();
		}

	for(i = 0; i < (VERTICES)/2; i++)
		{
		ED[i]		= (unsigned int *)  NewPtr(sizeof(int)*VERTICES);
		if(ED[i] 	== NULL) Mem_Error();
		EXL[i]		= (unsigned int *)  NewPtr(sizeof(int)*4);
		if(EXL[i] 	== NULL) Mem_Error();
		EXR[i]      = (unsigned int *)  NewPtr(sizeof(int)*4);
		if(EXR[i] 	== NULL) Mem_Error();
		EXP[i]		= (unsigned int *)  NewPtr(sizeof(int)*4);
		if(EXP[i] 	== NULL) Mem_Error();
		NEX[i]		= (unsigned int *)  NewPtr(sizeof(int)*4);
		if(NEX[i] 	== NULL) Mem_Error();
		T[i]        = (unsigned char *) NewPtr(sizeof(char)*8);
		if(T[i] 	== NULL) Mem_Error();
		}
	
	for(i = 0; i < 2*MAXNUMRELATORS; i++)
		{
		DRA[i] = (unsigned int *) NewPtr(sizeof(int)*2*MAXNUMRELATORS);
		if(DRA[i] 	== NULL) Mem_Error();
		}
				
	for(i = 0; i < MAXNUMCOMPONENTS; i++)
		{
		CBC[i] = (unsigned char *) NewPtr(sizeof(char)*(MAXNUMGENERATORS + 2));
		if(CBC[i] 	== NULL) Mem_Error();
		MLC[i] = (unsigned long *) NewPtr(sizeof(long)*(MAXNUMGENERATORS + 1));
		if(MLC[i] 	== NULL) Mem_Error();
		}
	
	for(i = 1; i < 2*VERTICES; i++)
		{
		Face[i] = (unsigned char *) NewPtr(sizeof(char)*(VERTICES + 1));
		if(Face[i] 	== NULL) Mem_Error();
		}
			
	AT				= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(AT 			== NULL) Mem_Error();	
	BCF				= (unsigned char *)  NewPtr(sizeof(char)*2*VERTICES);
	if(BCF 			== NULL) Mem_Error();
	BCWG 			= (unsigned char *)	 NewPtr(sizeof(char)*(MAXNUMGENERATORS + 1));
	if(BCWG 		== NULL) Mem_Error();
	BDY				= (         int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(BDY 			== NULL) Mem_Error();
	Bdry			= (unsigned	int  *)  NewPtr(sizeof(int)*(VERTICES + 2));
	if(Bdry 		== NULL) Mem_Error();
	BeenChecked		= (unsigned int  *)	 NewPtr(sizeof(int)*MAX_MIN_GEN_PRES);
	if(BeenChecked  == NULL) Mem_Error();
	BSV1			= (unsigned int  *)  NewPtr(sizeof(int)*MAX_SAVED_LEVELS);
	if(BSV1 		== NULL) Mem_Error();
	BSV2			= (unsigned int  *)  NewPtr(sizeof(int)*MAX_SAVED_LEVELS);
	if(BSV2 		== NULL) Mem_Error();
	ComponentNum 	= (unsigned int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(ComponentNum == NULL) Mem_Error();
	CS				= (unsigned char *)  NewPtr(sizeof(char)*(MAXNUMCOMPONENTS + 1));
	if(CS 			== NULL) Mem_Error();
	Daughters 		= (unsigned int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(Daughters 	== NULL) Mem_Error();
	DeletedEdges	= (unsigned char *)	 NewPtr(sizeof(char)*6*VERTICES);
	if(DeletedEdges == NULL) Mem_Error();
	DF				= (unsigned int  *)  NewPtr(sizeof(int)*(VERTICES)/2);
	if(DF 			== NULL) Mem_Error();
	InDisk			= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(InDisk 		== NULL) Mem_Error();
	ER				= (			char *)  NewPtr(sizeof(char)*MAX_SAVED_PRES);
	if(ER 			== NULL) Mem_Error();
	Father			= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(Father 		== NULL) Mem_Error();
	Flags			= (         int  *)  NewPtr(sizeof(int)*VERTICES);
	if(Flags 		== NULL) Mem_Error();
	FR 				= (unsigned int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(FR 			== NULL) Mem_Error();
	FV 				= (unsigned char *)	 NewPtr(sizeof(char)*VERTICES);
	if(FV 			== NULL) Mem_Error();
	GBC 			= (unsigned char *)	 NewPtr(sizeof(char)*(MAXNUMRELATORS + 2));
	if(GBC 			== NULL) Mem_Error();
	GV2				= (unsigned int  *)  NewPtr(sizeof(int)*(VERTICES)/2);
	if(GV2 			== NULL) Mem_Error();
	GV2L			= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(GV2L 		== NULL) Mem_Error();
	GV2R			= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(GV2R 		== NULL) Mem_Error();
	HegSplNum		= (unsigned int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(HegSplNum    == NULL) Mem_Error();
	HegSplNxt		= (unsigned int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(HegSplNxt    == NULL) Mem_Error();
	InDisk			= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(InDisk 		== NULL) Mem_Error();	
	InPS			= (			int	 *)	 NewPtr(sizeof(int)*VERTICES);
	if(InPS 		== NULL) Mem_Error();	
	InQueue 		= (unsigned int  *)	 NewPtr(sizeof(int)*VERTICES);
	if(InQueue 		== NULL) Mem_Error();
	IV				= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(IV 			== NULL) Mem_Error();
	Left			= (unsigned int  *)  NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(Left 		== NULL) Mem_Error();
	Lowpt			= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(Lowpt 		== NULL) Mem_Error();
	LR 				= (unsigned long *)  NewPtr(sizeof(long)*(MAXNUMRELATORS + 1));
	if(LR 			== NULL) Mem_Error();
	LSP				= (unsigned long *)  NewPtr(sizeof(long)*(MAX_SAVED_PRES));
	if(LSP 			== NULL) Mem_Error();
	LSQ				= (unsigned long *)	 NewPtr(sizeof(long)*(MAX_SAVED_PRES));
	if(LSQ 			== NULL) Mem_Error();
	NCS				= (unsigned char *)	 NewPtr(sizeof(char)*(MAX_SAVED_PRES));
	if(NCS 			== NULL) Mem_Error();					
	NEBC 			= (unsigned int  *)	 NewPtr(sizeof(int)*2*MAXNUMRELATORS);
	if(NEBC 		== NULL) Mem_Error();
	N1H				= (unsigned char *)  NewPtr(sizeof(char)*MAXNUMCOMPONENTS);
	if(N1H 			== NULL) Mem_Error();
	NFBC			= (unsigned int  *)	 NewPtr(sizeof(int)*2*MAXNUMRELATORS);
	if(NFBC 		== NULL) Mem_Error();	
	NG 				= (         int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(NG 			== NULL) Mem_Error();
	NR 				= (         int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(NR 			== NULL) Mem_Error();
	NRBC 			= (unsigned int  *)	 NewPtr(sizeof(int)*2*MAXNUMRELATORS);
	if(NRBC 		== NULL) Mem_Error();
	NS1XD2			= (unsigned char *)  NewPtr(sizeof(char)*MAXNUMCOMPONENTS);
	if(NS1XD2 		== NULL) Mem_Error();
	NS1XS2  		= (unsigned char *)  NewPtr(sizeof(char)*MAXNUMCOMPONENTS);
	if(NS1XS2 		== NULL) Mem_Error();
	Number			= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(Number 		== NULL) Mem_Error();
	NumLoops		= (			 int  *) NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(NumLoops     == NULL) Mem_Error();
	OSA				= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(OSA 			== NULL) Mem_Error();
	OSB				= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(OSB 			== NULL) Mem_Error();
	PG				= (unsigned int  *)  NewPtr(sizeof(int)*(VERTICES)/2);
	if(PG 			== NULL) Mem_Error();
	PresName		= (unsigned char *)  NewPtr(MAXLENGTH + 1);
	if(PresName 	== NULL) Mem_Error();
	PRIM			= (         int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(PRIM 		== NULL) Mem_Error();
	QPM				= (unsigned char *)  NewPtr(sizeof(char)*(MAX_SAVED_PRES));
	if(QPM 			== NULL) Mem_Error();	
	UpDate 			= (unsigned int  *)	 NewPtr(sizeof(int)*VERTICES*VERTICES);
	if(UpDate 		== NULL) Mem_Error();
	Right			= (unsigned int  *)  NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(Right 		== NULL) Mem_Error();
	SaveBdry		= (         int  *)  NewPtr(sizeof(int)*VERTICES);
	if(SaveBdry 	== NULL) Mem_Error();
	SFSols			= (			int **)  NewPtr(sizeof(int)*(MAXNUMCOMPONENTS + 1));
	if(SFSols	 	== NULL) Mem_Error();
	SLP				= (unsigned char ****) NewPtr(sizeof(long)*MAX_SAVED_PRES);
	if(SLP	 		== NULL) Mem_Error();
	SMGP			= (unsigned char ****) NewPtr(sizeof(long)*MAX_MIN_GEN_PRES);
	if(SMGP	 		== NULL) Mem_Error();
	SUR				= (unsigned char ****) NewPtr(sizeof(long)*MAX_SAVED_PRES);
	if(SUR 			== NULL) Mem_Error();
	SUR_Num			= (unsigned int  *)	 NewPtr(sizeof(int)*MAX_MIN_GEN_PRES);
	if(SUR_Num 		== NULL) Mem_Error();
	SURNumX			= (			 int *)  NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(SURNumX 		== NULL) Mem_Error();	
	SURL			= (unsigned long *)	 NewPtr(sizeof(long)*MAX_SAVED_PRES);
	if(SURL 		== NULL) Mem_Error();
	SV				= (unsigned int  *)  NewPtr(sizeof(int)*(VERTICES)/2);
	if(SV 			== NULL) Mem_Error();
	TP				= (unsigned char *)  NewPtr(sizeof(char)*MAX_SAVED_PRES);
	if(TP 			== NULL) Mem_Error();
	TV				= (unsigned int  *)  NewPtr(sizeof(int)*(VERTICES)/2);
	if(TV 			== NULL) Mem_Error();
	UDV 			= (unsigned int  *)	 NewPtr(sizeof(int)*MAX_SAVED_PRES);
	if(UDV			== NULL) Mem_Error();
	V				= (unsigned int  *)	 NewPtr(sizeof(int)*VERTICES);
	if(V 			== NULL) Mem_Error();
	VA 				= (unsigned int  *)	 NewPtr(sizeof(int)*(VERTICES)/2);
	if(VA 			== NULL) Mem_Error();
	VWG				= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(VWG		 	== NULL) Mem_Error();
	X				= (         int  *)	 NewPtr(sizeof(int)*VERTICES);
	if(X 			== NULL) Mem_Error();
	XX				= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(XX 			== NULL) Mem_Error();
	Y				= (         int  *)	 NewPtr(sizeof(int)*VERTICES);
	if(Y 			== NULL) Mem_Error();
	YY				= (unsigned int  *)  NewPtr(sizeof(int)*VERTICES);
	if(YY 			== NULL) Mem_Error();
	ZZ 				= (unsigned int  *)	 NewPtr(sizeof(int)*VERTICES);
	if(ZZ 			== NULL) Mem_Error();	
	zz 				= (unsigned int  *)	 NewPtr(sizeof(int)*2*MAXNUMRELATORS);
	if(zz 			== NULL) Mem_Error();
		
	for(i = 0; i < MAX_SAVED_PRES; i++)
		{
		SUR[i] = (unsigned char ***) NewPtr(sizeof(long)*(MAXNUMRELATORS + 1));
		if(SUR[i] 	== NULL) Mem_Error();
		for(j = 0; j <= MAXNUMRELATORS; j++) SUR[i][j] = NULL;
		}
			
	for(i = 0; i < MAX_SAVED_LEVELS; i++)
		{
		SLR[i] = (unsigned char ***) NewPtr(sizeof(long)*(MAXNUMRELATORS + 1));
		if(SLR[i] 	== NULL) Mem_Error();
		for(j = 0; j <= MAXNUMRELATORS; j++) SLR[i][j] = NULL;
		}
		
	for(i = 0; i < MAX_MIN_GEN_PRES; i++)
		{
		SMGP[i] = (unsigned char ***) NewPtr(sizeof(long)*(MAXNUMRELATORS + 1));
		if(SMGP[i] 	== NULL) Mem_Error();
		for(j = 0; j <= MAXNUMRELATORS; j++) SMGP[i][j] = NULL;
		}	
	
	for(i = 0; i < VERTICES; i++)
        {
        MM[i] 		= (unsigned char *) NewPtr(sizeof(char)*VERTICES);
        if(MM[i] 	== NULL) Mem_Error();
        TC[i]		= (unsigned char *) NewPtr(sizeof(char)*VERTICES);
        if(TC[i] 	== NULL) Mem_Error();
        }
        
    for(i = 0; i <= MAXNUMDUPS; i++)
		{
		RWR[i] = (unsigned char *) NewPtr(sizeof(char)*125);
		if(RWR[i]   == NULL) Mem_Error();
		} 
		
	for(i = 1; i <= MAXNUMCOMPONENTS; i++) SFSols[i] = NULL;	
		
	Low						= (unsigned char *) NewPtr(sizeof(char)*VERTICES);
	if(Low					== NULL) Mem_Error();
	Num						= (unsigned char *) NewPtr(sizeof(char)*VERTICES);
	if(Num					== NULL) Mem_Error();
	OnLStack				= (unsigned char *) NewPtr(sizeof(char)*VERTICES);
	if(OnLStack				== NULL) Mem_Error();
	SatEdgeList1			= (unsigned char *) NewPtr(sizeof(char)*VERTICES*VERTICES);
	if(SatEdgeList1			== NULL) Mem_Error();
	SatEdgeList2			= (unsigned char *) NewPtr(sizeof(char)*VERTICES*VERTICES);
	if(SatEdgeList2			== NULL) Mem_Error();
	SComp					= (unsigned char *) NewPtr(sizeof(char)*VERTICES);
	if(SComp				== NULL) Mem_Error();
	SepVertexList			= (unsigned char *) NewPtr(sizeof(char)*VERTICES*VERTICES);
	if(SepVertexList		== NULL) Mem_Error();
	SinkSet					= (unsigned char *) NewPtr(sizeof(char)*VERTICES);
	if(SinkSet				== NULL) Mem_Error();
	LStack					= (unsigned char *) NewPtr(sizeof(char)*VERTICES);
	if(LStack				== NULL) Mem_Error();
	Num_Sep_Vertex_Pairs    = (unsigned int *) NewPtr(sizeof(int)*(MAX_SAVED_LEVELS + 1));
	if(Num_Sep_Vertex_Pairs	== NULL) Mem_Error();
	UnUsed_Sep_Vertex_Pairs = (unsigned int *) NewPtr(sizeof(int)*(MAX_SAVED_LEVELS + 1));
	if(UnUsed_Sep_Vertex_Pairs	== NULL) Mem_Error();

	for(i = 0; i <= MAXNUMRELATORS; i++)
		{ 
		Copy_Of_Input[i] 	= NULL;
		Copy_Of_Rel_1[i] 	= NULL;
		Copy_Of_Rel_2[i] 	= NULL;
		CD_Surgery_Rel[i]	= NULL;
		DelRelators[i] 		= NULL;
		DualRelators[i] 	= NULL;
		Exp_Surgery_Rel[i]	= NULL;
		KorLRelators[i] 	= NULL;
		OutRelators[i] 		= NULL;
		Relators[i] 		= NULL;
		TopOfChain[i] 		= NULL;
		WirtingerL[i]		= NULL;
		WirtingerM[i]		= NULL;
		}															
	
	Temp1 		= NULL;
	Temp2 		= NULL;
	Temp3 		= NULL;
	Temp4 		= NULL;
	Temp5		= NULL;
	Temp6		= NULL;
	Temp7		= NULL;
	Temp8		= NULL;
	Temp9		= NULL;
	Temp10		= NULL;
	Temp11		= NULL;
	Temp12		= NULL;
	Temp13		= NULL;
	Temp14		= NULL;
	Temp15		= NULL;
	Temp16      = NULL;
	
	Inst = (unsigned char *) NewPtr(30000);		/* Get some storage for any reasonable input. */
	if(Inst     == NULL) Mem_Error();
			
	ReadDateTime(&Seconds);
	srand(abs(LoWord(Seconds)));
		
	Did_Exponent_Surgery 		= FALSE;
	Did_Cutting_Disk_Surgery 	= FALSE;
	
	return(FALSE);	
}
