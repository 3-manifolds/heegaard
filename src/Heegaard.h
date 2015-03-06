#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*    #define DEBUGGING            */
/*    #define MICRO_PRINT          */
    #define PRINT
/*    #define PRINT_CYCLES         */    
/*    #define PRINT_LINKS   	   */      


#define Size size_t

#define LoWord(x) ((short)(x))
#define charCodeMask        0x000000FF

#define ANNULUS_EXISTS      53100                 
#define BANDSUM             5
#define BDRY_UNKNOWN        100
#define BEEN_USED			2999
#define BEGIN               1
#define BIG_NUMBER          1000000
#define BYTES_AVAILABLE		41943040L
#define CAN_NOT_DELETE      9996
#define DELETED_RELATOR     50200
#define DONE                50000
#define DUALIZE             2
#define DUPLICATE           50300
#define EOS                 '\0'
#define ERROR               0.1
#define FALSE               0
#define FATAL_ERROR         99999
#define FOUND_ELSEWHERE     60100
#define FULL_HOUSE			51234
#define GENERIC_LENS_SPACE  59000
#define INFINITE            65535
#define INITIAL_PRES        0
#define INTERRUPT			61234
#define KNOWN_LENS_SPACE    60200
#define MAXCOUNT            5
#define MAXNUMCOMPONENTS 	VERTICES
#define MAXLENGTH           16000
#define MAXNUMDUPS          1000
#define MAXNUMGENERATORS    26    
#define MAXNUMRELATORS      31
#define MAX_MIN_GEN_PRES	25000
#define MAX_SAVED_LEVELS    5000
#define MAX_SAVED_PRES      50003
#define MISSING_GEN_DONE1   60350
#define MISSING_GEN_DONE2   60360
#define NO_ERROR            0
#define NON_PLANAR          3300
#define NORMAL              6
#define NOT_CONNECTED       52000
#define NON_UNIQUE_1        50500
#define NON_UNIQUE_2        50600
#define NON_UNIQUE_3        50700
#define NON_UNIQUE_4        50800
#define REDUCE_GENUS        13
#define RERUN               3
#define RESET               14
#define S1_X_D2             60300
#define S1_X_S2             60400
#define S1_X_X2             50900
#define SEP_PAIRS           51000
#define SPLIT            	53000
#define STOP                1
#define THREE_SPHERE        60600
#define TOO_LONG            64000
#define TOO_MANY_COMPONENTS 9997
#define TRUE                1
#define UNKNOWN             59999
#define V2_ANNULUS_EXISTS   51100
#define VERTICES            2*MAXNUMGENERATORS                            

/******************************* List of Function Prototypes **********************************/

extern int		AlexanderPolynomial(unsigned char*);								/* In Heegaard21.c	*/
extern int		AlexanderPolynomial_Freely_Reduce(unsigned char*, unsigned int);	/* In Heegaard21.c	*/
extern unsigned int	AlexanderPolynomial_GCD(unsigned int, unsigned int);			/* In Heegaard21.c	*/
extern int		Annulus(unsigned int,unsigned int,int,int);							/* In Heegaard7.c	*/
extern int		BatchProcessing(void);												/* In Heegaard14.c	*/
extern int 		Batch_Report(int*);													/* In Heegaard12.c	*/
extern int		Been_Seen(void);													/* In Heegaard7.c	*/
extern int		Canonical_Rewrite(unsigned char ***,int,int);						/* In Heegaard5.c	*/
extern int		Check_Bridge_Interlacing(unsigned int,unsigned int);				/* In Heegaard3.c	*/ 
extern int		Check_Connected(void);												/* In Heegaard3.c	*/ 
extern int		Check_For_Primitives(int, int);										/* In Heegaard6.c	*/
extern int      Check_HS_Reps(int,int*);											/* In Heegaard18.c	*/
extern int      Check_HS_Simple_Circuits(int,int*);									/* In Heegaard18.c	*/
extern void		Check_HS_Uniqueness(int NumHSReps,int* HSRepL);						/* In Heegaard18.c	*/
extern int		Check_HS_Uniqueness_Sub1(int,int);									/* In Heegaard18.c	*/
extern int		Check_Level_Transformations(void);									/* In Heegaard2.c	*/
extern void		CheckPolySymmetry(int *, unsigned int);								/* In Heegaard21.c	*/
extern int		CheckPrimitivity(void);												/* In Heegaard2.c	*/
extern int		Check_Realizability_Of_The_Initial_Presentation(void);				/* In Heegaard1.c	*/
extern int		Compare(unsigned char *);											/* In Heegaard4.c	*/
extern int		Compare_Dual_Pres(int);												/* In Heegaard1.c	*/
extern int		Compare_Input_Pres(void);											/* In Heegaard1.c	*/
extern int		Compare_Pres(int);													/* In Heegaard1.c	*/
extern int		Compare_Str(unsigned char *,unsigned char *,unsigned long);			/* In Heegaard4.c	*/
extern void		ComputeAlexanderPolynomial(void);									/* In Heegaard21.c	*/
extern void		Compute_Homology(void);												/* In Heegaard16.c	*/
extern void		ComputePrimitives(unsigned int,unsigned int,unsigned int,unsigned int);	/* In Heegaard21.c	*/
extern int		ComputeValences_A(void);											/* In Heegaard2.c	*/
extern int		Connected_(unsigned int,unsigned int);								/* In Heegaard3.c	*/ 
extern int		Connected_AJ3(unsigned int,unsigned int);							/* In Heegaard2.c	*/
extern int		CountMaxMins(char *);												/* In Heegaard21.c	*/
extern unsigned int	Count_Sep_Pairs(unsigned int);									/* In Heegaard7.c	*/
extern int  	CP_Check_Simple_Paths(unsigned int,int*,int,unsigned char**,
					unsigned char**);												/* In Heegaard5.c   */
extern int  	CHSP_Check_Simple_Circuits(unsigned int,int*,int,unsigned char**,
					unsigned char**);												/* In Heegaard18.c   */					
extern void		CP_Concatenate_Paths(unsigned char **,unsigned char **,int);		/* In Heegaard5.c	*/
extern int		CP_Do_Aut(unsigned int,char);										/* In Heegaard5.c	*/
extern void		CP_Fill_AA(char);													/* In Heegaard5.c	*/
extern int		CP_Find_Primitives(char);											/* In Heegaard5.c	*/
extern int		Cutting_Disk_Surgery_Diagram_7(void);								/* In Heegaard10.c	*/
extern void		Debug(void);														/* In Heegaard4.c	*/
extern int		Defining_Relator(int,int,int,int);									/* In Heegaard6.c	*/
extern int		Defining_Relator_SPC(void);											/* In Heegaard17.c	*/
extern int		Delete_Dups(void);													/* In Heegaard1.c	*/
extern void		Delete_Old_Presentations(void);										/* In Heegaard1.c	*/
extern void		Delete_Old_PresentationsSLP(void);									/* In Heegaard18.c	*/
extern void		Delete_Old_PresentationsSMGP(int,unsigned int *);					/* In Heegaard18.c	*/
extern int		Delete_Redundant_Relators(void);									/* In Heegaard4.c	*/
extern int		Delete_Trivial_Generators(int);										/* In Heegaard7.c	*/
extern int		Diagram_1(void);													/* In Heegaard4.c	*/
extern int		Diagram_2(void);													/* In Heegaard4.c	*/
extern void		Diagram_3(void);													/* In Heegaard4.c	*/
extern void		Diagram_4(void);													/* In Heegaard4.c	*/
extern unsigned int	Diagram_5(void);												/* In Heegaard4.c	*/
extern void		Diagram_6(void);													/* In Heegaard4.c	*/
extern unsigned int	Diagram_7(void);												/* In Heegaard4.c	*/
extern int		Diagram_Data(int,int,int,int,int);									/* In Heegaard3.c	*/ 
extern void		Diagram_Data_for_Graphviz(int,int,int);								/* In Heegaard3.c	*/
extern unsigned int	Diagram_Main(void);												/* In Heegaard4.c	*/
extern int		Display_A_Diagram(int,int,int);										/* In Heegaard12.c	*/
extern void		Display_Diagrams(void);												/* In Heegaard12.c	*/
extern void		Display_Diagram_Of_The_Initial_Presentation(void);					/* In Heegaard1.c	*/
extern int		Display_HS_Diagrams(int,int*);										/* In Heegaard18.c	*/
extern int		Do_Aut(unsigned int,unsigned int,int);								/* In Heegaard2.c	*/
extern int		Do_Auts(unsigned int,unsigned int,unsigned int);					/* In Heegaard2.c	*/
extern int		Do_Aut_L(void);														/* In Heegaard7.c	~*/
extern int		Do_Aut_L_Long(int,int,int);											/* In Heegaard7.c	*/
extern int		Do_Aut_L_Long_2(int V1,int V2,int TheComp);							/* In Heegaard7.c	*/
extern int		Do_Aut_SF(int);														/* In Heegaard19.c	*/            
extern int		Do_Initialization(void);											/* In Heegaard15.c	*/      
extern int		Do_LT(unsigned int,int,int);										/* In Heegaard13.c	*/
extern int		Empty_Relator_BS(void);												/* In Heegaard8.c	*/
extern int		Empty_Relator_D(void);												/* In Heegaard8.c	*/
extern void		Fatal_Error(void);													/* In Heegaard12.c	*/
extern int 		FiberedCheck(void);													/* In Heegaard21.c	*/
extern int		Final_Rewrite(unsigned char ***);									/* In Heegaard5.c	*/
extern void		Fill_A(int);														/* In Heegaard2.c	*/
extern void		Fill_AA(int);														/* In Heegaard2.c	*/
extern void		Fill_DRA(void);														/* In Heegaard4.c	*/
extern void		Find_Cancellation_Paths(int,int,int);								/* In Heegaard5.c	*/
extern int		Find_Canonical_Orbit_Reps(int *,int);								/* In Heegaard18.c	*/
extern int		Find_Cut_Vertices(void);											/* In Heegaard3.c	*/ 
extern int		Find_Flow_A(int,int);												/* In Heegaard2.c	*/
extern int		Find_Flow_B(unsigned int);											/* In Heegaard6.c	*/
extern int		Find_Level_Flow(unsigned int);										/* In Heegaard13.c	*/
extern int		Find_Level_Transformations(int,int);								/* In Heegaard13.c	*/
extern int		Find_Level_Transformations_Of_The_Initial_Presentation(void);		/* In Heegaard1.c	*/
extern unsigned int	Find_Minimal_Path(void);										/* In Heegaard3.c	*/ 
extern int		Find_Path(int,int,int,int);											/* In Heegaard7.c	*/
extern int		Find_Primitives(int);												/* In Heegaard2.c	*/
extern int 		Find_Simple_Circuits(void);											/* In Heegaard18.c	*/
extern int		Find_Symmetries(int);												/* In Heegaard5.c	*/
extern int		Freely_Reduce(void);												/* In Heegaard2.c	*/
extern void		Freely_Reduce_Nr(void);												/* In Heegaard5.c	*/
extern void		Gauss_Seidel(void);													/* In Heegaard3.c	*/ 
extern unsigned long	gcd(unsigned long,unsigned long);							/* In Heegaard16.c	*/
extern unsigned int	GCD(unsigned int,unsigned int);									/* In Heegaard4.c	*/
extern int		Genus_Two_Meridian_Reps(int,int);									/* In Heegaard20.c	*/
extern int		Genus_Two_Meridian_Reps_Sub(unsigned char*, unsigned char*);		/* In Heegaard20.c	*/
extern int		Genus_Two_Seifert_Fibered(int);										/* In Heegaard19.c	*/
extern void		Get_Bdry_Comps(int,int,unsigned int);								/* In Heegaard4.c	*/
extern void		Get_Components(void);												/* In Heegaard7.c	*/
extern int		Get_Connected_Components(void);										/* In Heegaard6.c	*/
extern int		Get_Diagrams(void);													/* In Heegaard2.c	*/
extern int		Get_Genus_2_SF_EXPS1(void);											/* In Heegaard19.c	*/
extern int		Get_Genus_2_SF_EXPS2(void);											/* In Heegaard19.c	*/
extern int		Get_Initial_Diagram(int);											/* In Heegaard1.c	*/
extern int 		Get_LTs(unsigned int,int,int);										/* In Heegaard13.c	*/
extern int		Get_Matrix(void);													/* In Heegaard3.c	*/ 
extern unsigned int	Get_MinExp(unsigned int,int);									/* In Heegaard2.c	*/
extern int      Get_Next_Presentation_From_File(void);								/* In Heegaard18.c	*/
extern int		Get_Presentation_From_File(void);									/* In Heegaard1.c	*/
extern int		Get_Presentation_From_KeyBoard(void);								/* In Heegaard1.c	*/
extern int		Get_Relators_From_SUR(int);											/* In Heegaard2.c	*/
extern int		Get_SF_Alphas1(int);												/* In Heegaard19.c	*/
extern int		Get_SF_Alphas2(int);												/* In Heegaard19.c	*/
extern int		Get_SF_Invariants(int);												/* In Heegaard19.c 	*/
extern void		Get_Simplification_Parameters_From_User(int,int);					/* In Heegaard1.c	*/
extern int		Get_User_Input_SPC(char);											/* In Heegaard17.c	*/		
extern void		Heegaard_Splash_Screen(void);										/* In Heegaard10.c	*/
extern int		ID_A_PMQPM(unsigned int i);											/* In Heegaard18.c	*/
extern void		ID_PMQPM(int, char*, unsigned int*);								/* In Heegaard18.c	*/
extern int 		Init_Genus_Two_Seifert_Fibered(int*,int,int);						/* In Heegaard19.c	*/
extern unsigned int	In_File(void);	     	      									/* In Heegaard9.c	*/
extern unsigned int	In_File2(int, unsigned char ***);								/* In Heegaard18.c	*/
extern int		Initial_Realizability_Check(void);									/* In Heegaard1.c	*/
extern int		Init_Find_Level_Transformations(int);								/* In Heegaard9.c	*/
extern void		Init_G_Variables(void);												/* In Heegaard1.c	*/
extern void		Init_SCdfsR(void);													/* In Heegaard13.c	*/
extern void		Init_TCR(void);														/* In Heegaard13.c	*/
extern void		Inverse(unsigned char *);											/* In Heegaard4.c	*/
extern void		Inverse_Nr(unsigned char *);										/* In Heegaard5.c	*/
extern int 		Is_IP_In_HS_Reps(int,int*);											/* In Heegaard18.c	*/
extern int		Is_Knot_Relator(void);	 											/* In Heegaard20.c	*/
extern int		Is_Sep_Pair(unsigned int,unsigned int,unsigned int);				/* In Heegaard7.c	*/
extern int		Just_Delete_Primitives(char);	  	       							/* In Heegaard17.c	*/
extern unsigned int	 Lens_Space(void);												/* In Heegaard6.c	*/
extern int		Lens_Space_D(int);													/* In Heegaard6.c	*/
extern int		Level_Transformations(int,int,int);									/* In Heegaard7.c	*/
extern int		Level_Transformations_2(int,int,int,unsigned int);					/* In Heegaard7.c	*/
extern int		Level_Trans_Reset(unsigned int,unsigned int,unsigned int);			/* In Heegaard7.c	*/
extern unsigned long	LGCD(unsigned long, unsigned long);		   					/* In Heegaard6.c	*/
extern int		Long_Mult(unsigned char,unsigned char,unsigned int,unsigned int,
			  unsigned int,unsigned int,unsigned int,unsigned int,
			  unsigned int,long int,unsigned long);									/* In Heegaard5.c	*/
extern int		main(int argv, char **argc);	 									/* In Heegaard1.c	*/
extern void		Mark_As_Duplicate(unsigned int);									/* In Heegaard2.c	*/
extern void		Mark_As_Found_Elsewhere(int);										/* In Heegaard2.c	*/
extern void		Mem_Error(void);													/* In utils.c	*/
extern void		MergeHegSpl(unsigned int,unsigned int);								/* In Heegaard18.c	*/
extern int		MG_Bdry_Comp_Data(unsigned int);									/* In Heegaard4.c	*/
extern void		Micro_Print_Bandsum(void);											/* In Heegaard12.c	*/
extern void		Micro_Print_Do_Aut(unsigned int,unsigned int);						/* In Heegaard12.c	*/
extern void		Micro_Print_Dualize(void);		 									/* In Heegaard12.c	*/
extern void		Micro_Print_Freely_Reduce(unsigned long,unsigned long);				/* In Heegaard12.c	*/
extern void		Micro_Print_Level_Transformations(unsigned int,unsigned int,
	                                                  unsigned int,unsigned int);	/* In Heegaard7.c	*/
extern void		Micro_Print_Level_Transformations_Reset(unsigned int);				/* In Heegaard7.c	*/
extern void		Micro_Print_Reset(void);			 								/* In Heegaard12.c	*/
extern int		Missing_Gen(void);				 									/* In Heegaard8.c	*/
extern char		mykbhit(void);														/* In Heegaard1.c	*/
extern int		New_Relator(int);													/* In Heegaard5.c	*/ 
extern int		Non_Unique(void);													/* In Heegaard4.c	*/
extern int		Non_Unique_Initial_Diagram(void);									/* In Heegaard1.c	*/
extern unsigned int	Offset(void);													/* In Heegaard4.c	*/
extern unsigned int	OffsetSub(unsigned int,unsigned int,unsigned int,unsigned int,
			          unsigned int,unsigned int);									/* In Heegaard4.c	*/
extern int		OneGenerator_SPC(unsigned char);									/* In Heegaard17.c	*/
extern int		OneRelator_SPC(unsigned char);										/* In Heegaard17.c	*/
extern unsigned int	On_File(void);													/* In Heegaard1.c	*/
extern int		Planar(int,int);													/* In Heegaard3.c	*/ 
extern int		Planar_Connected_(unsigned int);								    /* In Heegaard3.c	*/	
extern void		Print_Bdry_Comp_Info(int,int,int);									/* In Heegaard3.c	*/ 
extern void		Print_Bdry_Data(unsigned int);										/* In Heegaard12.c	*/
extern void		Print_Bdry_Data2(unsigned int);										/* In Heegaard12.c	*/
extern void		Print_DelRelators(void);											/* In Heegaard12.c	*/
extern void		Print_DualRelators(int,int,int,int);								/* In Heegaard12.c	*/
extern int		Print_Graph(int,int,int,int);										/* In Heegaard3.c	*/ 
extern void		Print_OutRelators(int,int,int,int);									/* In Heegaard12.c	*/
extern void		Print_Realizability(int,unsigned int);								/* In Heegaard1.c	*/
extern void		Print_Relators(unsigned char ***,int);								/* In Heegaard12.c	*/
extern void		Print_Relators2(unsigned char ***,int);								/* In Heegaard12.c	*/
extern int		Print_SFComp(int);													/* In Heegaard12.c	*/
extern void		Print_SLR(int,int);	      											/* In Heegaard12.c	*/
extern unsigned int	Proper_Power(void);												/* In Heegaard6.c	*/
extern void		Prune_Search_Tree(void);											/* In Heegaard16.c	*/
extern void		qksort(unsigned int);												/* In qksort.c		*/
extern void		qksort2(int,int,int,unsigned int*);									/* In Heegaard18.c	*/
extern int		qkst_compare(int,int);		    		      	       	    		/* In Heegaard12.c	*/
extern int		qkst_compare2(int,int,int,unsigned int*);							/* In Heegaard18.c	*/
extern void		qkst_swap(int,int);		   											/* In Heegaard12.c	*/
extern void		qkst_swap2(int,int);												/* In Heegaard18.c	*/
extern int		Random_Sep_Pair(unsigned int);										/* In Heegaard7.c	*/
extern void		Realization_Warning(void);											/* In Heegaard1.c	*/
extern unsigned long	Recip_Mod_P(unsigned long,unsigned long);					/* In Heegaard6.c	*/
extern unsigned int	Reduce_Genus(int,int,int);	   									/* In Heegaard6.c	*/
extern int		Reduce_The_Initial_Presentation_To_Minimal_Length(void);			/* In Heegaard1.c	*/
extern int		Reduce_The_Initial_Presentation_To_Minimal_Length_SPC(void);		/* In Heegaard17.c	*/
extern int		Report(long,long,unsigned int,unsigned int,unsigned char,unsigned char,
			       unsigned char,unsigned char,unsigned char,unsigned char *);		/* In Heegaard12.c	*/
extern int		Report_SPC(int *);	     		  	       	    					/* In Heegaard17.c	*/
extern void		Report_Symmetries(unsigned char *,int,int);							/* In Heegaard5.c	*/
extern int		ReRun_A_Presentation(void);											/* In Heegaard1.c	*/
extern int		Resolve(unsigned int,unsigned int,unsigned long);					/* In Heegaard8.c	*/
extern int		Rewrite_Input(void);	      		   								/* In Heegaard1.c	*/
extern int		Rewrite_Orbit_Reps(int,int,unsigned int*);							/* In Heegaard18.c	*/
extern int		Save_Pres(unsigned int,unsigned int,unsigned long,int,int,int,
				  unsigned char,char);												/* In Heegaard1.c	*/
extern int		Save_Pres2(void);													/* In Heegaard18.c	*/
extern void		SCdfsR(unsigned char);												/* In Heegaard13.c	*/
extern int		Sep_Pairs(int,int,int);												/* In Heegaard3.c	*/ 
extern int		Sep_Pairs_Sub(int,int);												/* In Heegaard3.c	*/ 
extern int		Sep_Surface(void);													/* In Heegaard4.c	*/
extern void		SetLimits(void);													/* In Heegaard14.c	*/
extern int 		Set_Up_Simplification_Parameters(int *,int *,int *, int *);			/* In Heegaard14.c	*/
extern int 		Set_Up_Simplification_Parameters_S1(void);							/* In Heegaard14.c	*/
extern int		SetUp_TopOfChain(void);												/* In Heegaard2.c	*/
extern int		SF_Sort_And_Print(int,int,int,int,int,int,int,int,int,int*);		/* In Heegaard19.c	*/
extern unsigned int	Slide_LComp(unsigned int, unsigned int,int,int,unsigned int,int,
				    int,int);														/* In Heegaard7.c	*/
extern int		Slide_ValenceTwo_Comp(int,unsigned int,unsigned int,unsigned int);	/* In Heegaard7.c	*/
extern int		SnapPy2Heegaard(void);												/* In Heegaard14.c	*/
extern int		Sort_Presentations_In_Memory(int);			     					/* In Heegaard12.c	*/
extern unsigned int	Split_At_Empty_Relators(int);									/* In Heegaard6.c	*/
extern void		Split_At_Empty_Relators_Sub1(void);									/* In Heegaard6.c	*/
extern int		Split_At_Empty_Relators_Sub2(int,int);								/* In Heegaard6.c	*/
extern void		Split_Relators(unsigned char);										/* In Heegaard3.c	*/ 
extern int		Splitting_Pres_On_File(int,int);									/* In Heegaard6.c	*/
extern int		Stabilize(void);													/* In Heegaard10.c	*/
extern int		Sub_Str(unsigned char *,unsigned char *);							/* In Heegaard8.c	*/
extern void		TCdfsR(unsigned char);		      									/* In Heegaard13.c	*/
extern void		Test_LT_For_Pseudo_Min(void);										/* In Heegaard9.c	*/
extern int		Test_New_Pres(void);												/* In Heegaard2.c	*/
extern void		Test_Sub_Str(unsigned int);											/* In Heegaard8.c	*/
extern void		Test_Transverse(void);												/* In Heegaard19.c	*/
extern void 	Too_Many_Components_ALert(void);									/* In Heegaard11.c	*/
extern int		Transverse(unsigned char *);										/* In Heegaard6.c	*/
extern int		Try_Cutting_Disk_Surgery(void);										/* In Heegaard10.c	*/
extern int		Try_Exponent_Surgery(void);											/* In Heegaard10.c	*/
extern void		Turn_Micro_Print_On(void);											/* In Heegaard1.c	*/
extern void		Update_Bdry_Data(void);												/* In Heegaard12.c	*/
extern void		UpDate_Fill_A(void);												/* In Heegaard2.c	*/
extern int		User_Says_Quit(void);												/* In Heegaard11.c	*/
extern unsigned int	Valence_Two(int);												/* In Heegaard8.c	*/
extern int		Valence_Two_Annulus(void);											/* In Heegaard8.c	*/
extern int		Verify_Length(unsigned char ***,int);								/* In Heegaard17.c	*/
extern char		WaitkbHit(void);													/* In Heegaard1.c	*/
extern unsigned int	Whitehead_Graph(void);											/* In Heegaard3.c	*/			  
extern int		Wirtinger(void);													/* In Heegaard10.c	*/ 


#ifndef MAC
/* From utils.c	*/
extern int    SysBeep(int seconds);
extern char   **NewHandle(size_t);
extern void   DisposeHandle(char **);
extern void   ReallocateHandle(char **, size_t);
extern size_t GetHandleSize(char **);
extern void   HUnlock(char **);
extern void   *NewPtr(size_t);
extern void   DisposePtr(void *);
extern size_t GetPtrSize(void *);
extern void   ReadDateTime(unsigned long *);
extern char   *ReadString(char *buffer, int size);
#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 4
 * fill-column: 78
 * End:
 */
