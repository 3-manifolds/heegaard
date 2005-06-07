#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
#include <MacMemory.h>
#include <Printing.h>
#include <Quickdraw.h>
#include <Sound.h>
*/

/*    #define DEBUGGING            */
/*    #define MICRO_PRINT            */
    #define PRINT
/*    #define PRINT_CYCLES        */    
/*    #define PRINT_LINKS            */


#define THPrint void*
#define Size size_t
#define GrafPtr void*


#define LoWord(x) ((short)(x))
#define charCodeMask        0x000000FF

#define ANNULUS_EXISTS      50100                    
#define BANDSUM             5
#define BDRY_UNKNOWN        100
#define BEGIN               1
#define BIG_NUMBER          1000000
#define CAN_NOT_DELETE      9996
#define DELETED_RELATOR     50200
#define DONE                50000
#define DUALIZE             2
#define DUPLICATE           50300
#define EOS                 '\0'
#define ERROR               0.1
#define FALSE               0
#define FATAL_ERROR         9999
#define FOUND_ELSEWHERE     60100
#define GENERIC_LENS_SPACE  59000
#define INFINITE            65535
#define INITIAL_PRES        0
#define KNOWN_LENS_SPACE    60200
#define MAXCOUNT            5
#define    MAXNUMCOMPONENTS 102
#define MAXLENGTH           16000
#define MAXNUMDUPS          1000
#define MAXNUMGENERATORS    26    
#define MAXNUMRELATORS      31
#define MAX_SAVED_LEVELS    300
#define MAX_SAVED_PRES      8003
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
#define    SPLIT            53000
#define STOP                1
#define THREE_SPHERE        60600
#define TOO_LONG            4500
#define TOO_MANY_COMPONENTS 9997
#define TRUE                1
#define UNKNOWN             59999
#define V2_ANNULUS_EXISTS   51100
#define VERTICES            2*MAXNUMGENERATORS                            

/******************************* List of Function Prototypes **********************************/

extern  int             Annulus(unsigned int,unsigned int,int,int);
extern  int             Been_Seen(void);
extern  int             Canonical_Rewrite(unsigned char ***,int,int);
extern  char            *cgets(char *);
extern  int             Check_Bridge_Interlacing(unsigned int,unsigned int);
extern  int             Check_Connected(void);
extern  int             Check_For_Primitives(int, int);
extern  int             Check_Level_Transformations(void);
extern  int             CheckPrimitivity(void);
extern  void            CheckPrintHandle(void);
extern  int             Check_Realizability_Of_The_Initial_Presentation(void);
/*extern  void            Click_On(Boolean);            */
extern  int             Compare(unsigned char *);
extern  int             Compare_Dual_Pres(int);
extern  int             Compare_Input_Pres(void);
extern  int             Compare_Pres(int);
extern  int             CompareStr(unsigned char *,unsigned char *);
extern  int             Compare_Str(unsigned char *,unsigned char *,unsigned long);
extern  void            Compute_Homology(void);
extern  int             ComputeValences_A(void);
extern  int             Connected_(unsigned int,unsigned int);
extern  int             Connected_AJ3(unsigned int,unsigned int);
extern  void            CP_Concatanate_Paths(unsigned char **,int);
extern  int             CP_Do_Aut(unsigned int);
extern  int             CP_Find_Primitives(void);
extern  int             cprintf(char *format,...);
extern  int             Cutting_Disk_Surgery_Diagram_7(void);
extern  void            Debug(void);
extern  int             Defining_Relator(int,int,int,int);
extern  int             Degenerate_Annulus(unsigned int,unsigned int,int);            
extern  int             Delete_Dup(void);
extern  int             Delete_Dups(void);
extern  void            Delete_Old_Presentations(void);
extern  int             Delete_Redundant_Relators(void);
extern  int             Delete_Trivial_Generators(int);
extern  int             DFS(unsigned int,unsigned int);
extern  int             Diagram_1(void);
extern  int             Diagram_2(void);
extern  void            Diagram_3(void);
extern  void            Diagram_4(void);
extern  unsigned int    Diagram_5(void);
extern  void            Diagram_6(void);
extern  unsigned int    Diagram_7(void);
extern  int             Diagram_Data(int,int);
extern  unsigned int    Diagram_Main(void);
extern  int             Display_A_Diagram(void);
extern  void            Display_Diagrams(void);
extern  void            Display_Diagram_Of_The_Initial_Presentation(void);
extern  int             Display_File_Input_Presentations(void);
extern  int             Display_File_Simplify_Heegaard(void);
extern  int             Display_Help_File(void);
extern  int             Do_Aut(unsigned int,unsigned int,int);
extern  int             Do_Auts(unsigned int,unsigned int,unsigned int);
extern  int             Do_Aut_L(void);
extern  int             Do_Aut_L_Long(int,int,int);
extern  int             Do_Aut_L_Long_2(int V1,int V2,int TheComp);                
extern  int             Do_Initialization(void);            
extern  void            DoPageSetUp(void);
extern  void            Edit_MyOut(void);                        
extern  int             Empty_Relator(void);
extern  int             Empty_Relator_BS(void);
extern  int             Empty_Relator_D(void);
extern  void            Fatal_Error(void);
extern  int             fclose(FILE *);
extern  char            *fgets(char *,int,FILE *);
extern  int             fprintf(FILE *,const char *,...);
extern  int             Final_Rewrite(unsigned char ***);
extern  void            Fill_A(int);
extern  void            Fill_AA(int);
extern  void            Fill_DRA(void);
extern  int             FillSA(void);
extern  int             Find_Annulus(unsigned char,int);
extern  void            Find_Cancellation_Paths(void);
extern  int             Find_Cut_Vertices(void);
extern  int             Find_Flow_A(int,int);
extern  int             Find_Flow_B(unsigned int);
extern  int             Find_Level_Flow(unsigned int);
extern  int             Find_Level_Transformations(int,int);
extern  int             Find_Level_Transformations_Of_The_Initial_Presentation(void);
extern  unsigned int    Find_Minimal_Path(void);
extern  int             Find_Path(int,int,int,int,unsigned int *);
extern  int             Find_Primitives(int);
extern  int             Find_Symmetries(int);
extern  void            Floating_Gauss_Seidel(void);
extern  int             Free_Memory_For_Find_Level_Transformations(int,int);
extern  int             Freely_Reduce(void);
extern  int             Freely_Reduce_Nr(void);
extern  void            Gauss_Seidel(void);
extern  unsigned long   gcd(unsigned long,unsigned long);
extern  unsigned int    GCD(unsigned int,unsigned int);
extern  int             Generate_Orbits_Under_Auts(void);
extern  int             Genus_Two_New_Relator(long int,unsigned long);
extern  int             Get_Annulus(int,unsigned char,int);
extern  int             Get_Attatchments(void);
extern  void            Get_Bdry_Comps(int,int,unsigned int);
extern  int             getch(void);
extern  int             Get_Components(void);
extern  int             Get_Connected_Components(void);
extern  int             Get_Diagrams(void);
extern  int             Get_Initial_Diagram(int);
extern  int             Get_Matrix(void);
extern  unsigned int    Get_MinExp(unsigned int,int);
extern  int             Get_Presentation_From_File(void);
extern  int             Get_Presentation_From_KeyBoard(void);
extern  int             Get_Relators_From_SUR(int);
extern  void            Get_Simplification_Parameters_From_User(int,int);
extern  int             HowMany(void);
extern  unsigned int    In_File(void);
extern  int             Initial_Realizability_Check(void);
extern  int             Init_Find_Level_Transformations(int);
extern  void            Init_G_Variables(void);
extern  void            Inverse(unsigned char *);
extern  int             Inverse_Nr(unsigned char *);
extern  unsigned int    Lens_Space(void);
extern  int             Lens_Space_D(int);
extern  int             Level_Transformations(int,int,int,int,int);
extern  int             Level_Transformations_2(int,int,int,unsigned int,int);
extern  int             Long_Mult(unsigned char,unsigned char,unsigned int,unsigned int,
                          unsigned int,unsigned int,unsigned int,unsigned int,
                          unsigned int,long int,unsigned long);
extern  int             main(void);
extern  void            Mark_As_Duplicate(unsigned int);
extern  void            Mark_As_Found_Elsewhere(int);
extern  int             MG_Bdry_Comp_Data(unsigned int);
extern  void            Micro_Print_Bandsum(void);
extern  void            Micro_Print_Do_Aut(unsigned int,unsigned int);
extern  void            Micro_Print_Dualize(void);
extern  void            Micro_Print_Freely_Reduce(unsigned long,unsigned long);
extern  void            Micro_Print_Level_Transformations(unsigned int,unsigned int,
			  unsigned int,unsigned int,unsigned int);
extern  void            Micro_Print_Level_Transformations_Reset(void);
extern  void            Micro_Print_Reset(void);
extern  int             Missing_Gen(void);
extern  void            MyDrawString(char *p);
extern  void            MyDrawText(char *,int);
extern  char            mykbhit(void);
extern  int             New_Relator(int);
extern  int             Non_Unique(void);
extern  int             Non_Unique_Initial_Diagram(void);
extern  unsigned int    Offset(void);
extern  unsigned int    OffsetSub(unsigned int,unsigned int,unsigned int,unsigned int,
                          unsigned int,unsigned int);
extern  int             Old_Generate_Orbits_Under_Auts(void);
extern  unsigned int    On_File(void);                 
extern  int             Planar(int,int);
extern  int             Planar_Connected_(unsigned int);        
extern  int             Plot_Graph(int,int,int);
extern  int             PrDoc(char **hText,long count,THPrint hPrint,int font,int size);
extern  void            Print_Bdry_Comp_Info(void);
extern  void            Print_Bdry_Data(unsigned int);
extern  void            Print_DelRelators(void);
extern  int             Print_Diagram_Data(int);
extern  void            Print_DualRelators(int);
extern  int             Print_Graph(int);
extern  void            Print_OutRelators(int);
extern  int             Print_Picture(int);
extern  int             Print_Presentation(int);
extern  void            Print_Realizability(int,unsigned int);
extern  void            Print_Relators(unsigned char ***,int,FILE *);
extern  void            Print_SLR(int);
extern  int             PrintText(char **HText,Size x,GrafPtr SavePort,int y);
extern  unsigned int    Proper_Power(void);
extern  void            Prune_Search_Tree(void);
extern  void            putch(char);
extern  void            qksort(unsigned int);
extern  int             qksort_Report(long,long,unsigned int,unsigned int,
			    unsigned int,unsigned int,unsigned int,unsigned int,
                            unsigned int, unsigned char *);
extern  int             qkst_compare(int i,int j);
extern  void            qkst_swap(int i,int j);                              
extern  int             rand(void);
extern  void            Realization_Warning(FILE *);
extern  unsigned long   Recip_Mod_P(unsigned long,unsigned long);
extern  unsigned int    Reduce_Genus(int,int,int);
extern  int             Reduce_Genus_Micro_Print(void);
extern  int             Reduce_The_Initial_Presentation_To_Minimal_Length(void);
extern  int             Report(long,long,unsigned int,unsigned int,unsigned int,
                               unsigned int,unsigned int,unsigned int,unsigned int,
                               unsigned char *);
extern  void            Report_Symmetries(unsigned char *,int,int);
extern  int             ReRun_A_Presentation(void);
extern  int             Resolve(unsigned int,unsigned int,unsigned long);
extern  void            rewind(FILE *);
extern  int             Rewrite_Input(void);
extern  int             Save_Pres(unsigned int,unsigned int,unsigned long,int,int,int,
                              unsigned char,char);
extern  int             Save_Pres_To_Input_Presentations(int);                     
extern  int             Sep_Pairs(int,int);
extern  int             Sep_Pairs_Sub(int,int);
extern  int             Sep_Surface(void);
extern  int             SetUp_TopOfChain(void);
extern  int             Slide_ValenceTwo_Comp(int,unsigned int,unsigned int);
extern  void            Sort_Presentations_In_Memory(void);
extern  void            Split_At_Empty_Relators_Sub1(void);
extern  int             Split_At_Empty_Relators_Sub2(int,int);
extern  void            Split_Relators(unsigned char);
extern  int             Splitting_Pres_On_File(int,int);
extern  void            srand(unsigned int);
extern  int             Sub_Str(unsigned char *,unsigned char *);
extern  int             Test_LT_For_Pseudo_Min(void);
extern  int             Test_New_Pres(void);
extern  int             Test_Sub_Str(unsigned int);
extern  int             Transverse(unsigned char *);
extern  int             Try_Cutting_Disk_Surgery(void);
extern  int             Try_Exponent_Surgery(void);
extern  void            Turn_Micro_Print_On(void);
extern  char            ungetch(char);
extern  void            Update_Bdry_Data(void);
extern  void            UpDate_Fill_A(void);
extern  int             User_Says_Quit(void);
extern  unsigned int    Valence_Two(int);
extern  int             Valence_Two_Annulus(void);
extern  char            WaitkbHit(void);
extern  unsigned int    Whitehead_Graph(void);                          
extern  int             Wirtinger(void);              

extern  char            *ReadString(char *buffer, int size);
