#include "Heegaard.h"
#include "Heegaard_Dec.h"

/****************************** function prototypes *****************************************
L  12 BatchProcessing()
L 747 Set_Up_Simplification_Parameters(int*,int*,int*, int*)
L 865 Set_Up_Simplification_Parameters_S1()
L 945 SetLimits()
L 964 SnapPy2Heegaard()
********************************************************************************************/

int BatchProcessing()
{
	unsigned char 	*p,
                   	*q,
                   	Stop;
      
    int	   			SDelete_Only_Short_Primitives,
    				SDo_Not_Reduce_Genus,
    				SFormBandsums,
    				SOnlyReducingBandsums,
    				SMicro_Print,
    				SMicro_Print_F;  
                            
    unsigned int    i;
    
    long            Scratch;
 
OPTIONS:  
	H_Results = NULL;
	  
    if((input_relators = fopen("Input_Presentations","r+")) == NULL)
        {
        SysBeep(5);
        printf("\nUnable to open the file 'Input_Presentations'.");
        printf("\nPlease locate the file 'Input_Presentations', make sure it is closed,");
        printf("\nand place it in the parent folder of the folder containing Heegaard.");
        return(1);
        }
    
    /**************************************************************************************
							Present the user with some options.
	**************************************************************************************/

	printf("\n");
	printf("\nHit 'a' TO COMPUTE ALEXANDER POLYNOMIALS OF 2-GENERATOR 1-RELATOR PRESENTATIONS.");
	printf("\nHIT 'b' TO FIND 'meridian' REPS M1 & M2 OF 2-GENERATOR 1-RELATOR PRESENTATIONS.");
		printf("\n	(This allows one to check if such presentations are knot exteriors.)");
	printf("\nHIT 'c' TO CHECK REALIZABILITY OF PRESENTATIONS.");
	printf("\nHIT 'C' TO CHECK IF THE INITIAL PRESENTATION IS A \042HS REP\042. (Heegaard will stop and alert the user if");  
	printf("\n    a sequence of handle-slides of the initial presentation P yields a presentation P' with |P'| < |P|.)");
	printf("\nHIT 'd' TO SEE DATA FOR HEEGAARD DIAGRAMS OF PRESENTATIONS.");
	printf("\nHIT 'D' TO SEE THE DUAL RELATORS OF EACH REALIZABLE PRESENTATION'S DIAGRAM.");
	printf("\nHIT 'E' TO HAVE HEEGAARD STABILIZE THE IP, COMPUTE HS REPS AND CHECK IF THE IP APPEARS ON THE HS LIST.");
	printf("\nHIT 'h' TO FIND THE INTEGRAL FIRST HOMOLOGY OF PRESENTATIONS.");
	printf("\nHIT 'l' TO FIND THE SIZE OF ORBITS OF PRESENTATIONS UNDER LEVEL TRANSFORMATIONS.");	      
	printf("\nHIT 'q' TO QUIT RUNNING IN BATCH MODE.");
	printf("\nHIT 'r' TO REDUCE AND SIMPLIFY PRESENTATIONS USING DEPTH-FIRST SEARCH AND SEP_VERT SLIDES.");
	printf("\nHIT 'R' TO REDUCE AND SIMPLIFY PRESENTATIONS USING BREADTH-FIRST SEARCH AND SEP_VERT SLIDES.");        
	printf("\nHIT 's' TO FIND SYMMETRIES OF PRESENTATIONS.");
	printf("\nHIT 'S' TO CONVERT SNAPPY FORMAT PRESENTATIONS TO HEEGAARD READABLE PRESENTATIONS.");
	printf("\nHIT 'u' TO STABILIZE PRESENTATIONS WHILE PRESERVING REALIZABILITY.");
	printf("\nHIT 'x' TO SIMPLIFY PRESENTATIONS BY SUCCESSIVELY DELETING PRIMITIVES, WITHOUT CHECKING REALIZABILITY.");
	printf("\nHIT 'X' TO FIND PRESENTATIONS OBTAINED BY DELETING PRIMITIVES FROM ONLY INITIAL PRESENTATIONS.");        
	printf("\nHIT 'z' TO REDUCE PRESENTATIONS TO MINIMAL LENGTH.\n");
	
GET_RESPONSE3:       
	switch(WaitkbHit())
		{
		case 'a':
			Batch = 1;
			H_Results = fopen("Heegaard_Results","a+");
			if(H_Results != NULL) 
				fprintf(H_Results,"\n\nAlexander polynomials of each 1-relator, 2-generator presentation should appear below.");			
			break;          	
		case 'b':
			Batch = 2;
			printf("\n\n	Note: If simplifying < A,B | M1,M2 > yields S^3, R 1 is a knot relator.");
			printf("\n  If Heegaard produces a pseudo-minimal PM diagram != S^3, R 1 is not a knot relator."); 
			H_Results = fopen("Heegaard_Results","a+");
			if(H_Results != NULL) 
				fprintf(H_Results,"\n\n'Meridian' reps of each uniquely realizable 1-relator, 2-generator presentation should appear below.");
				fprintf(H_Results," N.B.'Meridian' reps depend only on the relator R.");
			SetLimits();     			
			break;		
		case 'c':
			Batch = 3;
			NumRealizable = 0;
			H_Results = fopen("Heegaard_Results","a+");
			break;
		case 'C':
			Batch = 4;
			H_Results = fopen("Heegaard_Results","a+");
			if(H_Results != NULL) 
				fprintf(H_Results,"\n\nIndication that a presentation is a 'Heegaard Splitting Rep' should appear below.");			
			SetLimits();     
			break;					   
		case 'd':
			Batch = 5;
			printf("\n\nSHOW WHICH FACES OF EACH DIAGRAM FORM BDRY COMPONENTS ?  HIT 'y' OR 'n'.");
			GET_B5RESPONSE1:
			switch(WaitkbHit())
				{
				case 'y':
					B5PrintBdryComps = TRUE;
					break;
				case 'n':
					B5PrintBdryComps = FALSE;
					break;
				default:
					SysBeep(5);
					goto GET_B5RESPONSE1;
				}
			printf("\n\nPRINT DUAL RELATORS OF EACH DIAGRAM ? HIT 'y' OR 'n'.");
			GET_B5RESPONSE2:
			switch(WaitkbHit())
				{
				case 'y':
					B5PrintDualRelators = TRUE;
					break;
				case 'n':
					B5PrintDualRelators = FALSE;
					break;
				default:
					SysBeep(5);
					goto GET_B5RESPONSE2;
				}
			printf("\n\nPRINT PATHS CONNECTING FACES OF EACH DIAGRAM ? HIT 'y' OR 'n'.");
			GET_B5RESPONSE3:
			switch(WaitkbHit())
				{
				case 'y':
					B5PrintPaths = TRUE;
					break;
				case 'n':
					B5PrintPaths = FALSE;
					break;
				default:
					SysBeep(5);
					goto GET_B5RESPONSE3;
				}
			if(B5PrintPaths)
				{	
				GET_B5RESPONSE4:
				printf("\n\nCHECK SIMPLE CIRCUITS FOR PRIMITIVITY ETC ? HIT 'y' OR 'n'.");
				switch(WaitkbHit())
					{
					case 'y':
						B5TestSimpleCircuits = TRUE;
						break;
					case 'n':
						B5TestSimpleCircuits = FALSE;
						break;
					default:
						SysBeep(5);
						goto GET_B5RESPONSE4;
					}
				}						
			break;
		case 'D':
			Batch = 51;
			H_Results = fopen("Heegaard_Results","a+");
			if(H_Results != NULL) 
				fprintf(H_Results,"\n\nThe relators dual to each realizable initial presentation should appear below.");
			break;
		case 'E':
			Batch = 53;
			H_Results = fopen("Heegaard_Results","a+");
			SetLimits();
			break;					
		case 'h':
			Batch = 6;
			H_Results = fopen("Heegaard_Results","a+");
			if(H_Results != NULL) 
				fprintf(H_Results,"\n\nThe Z-homology of each initial presentation should appear below.");			
			break;
		case 'l':
			Batch = 7;
			H_Results = fopen("Heegaard_Results","a+");
			if(H_Results != NULL) 
				fprintf(H_Results,"\n\nInfo about the orbits under level-transformations of each initial presentation should appear below.");						
			break;									
		case 'q':
			fseek(input_relators,0L,0);
    		Delete_Old_Presentations();
    		Init_G_Variables();
    		NumFilled = 0;
			return(0);	
		case 'r':
			Batch = 10;
			Set_Up_Simplification_Parameters(&SFormBandsums,&SOnlyReducingBandsums,
				&SDelete_Only_Short_Primitives, &SDo_Not_Reduce_Genus);
			H_Results = fopen("Heegaard_Results","a+");
			if(H_Results != NULL)
				Set_Up_Simplification_Parameters_S1();										               
			break;		
		 case 'R':
		 	Batch = 11;	
			Set_Up_Simplification_Parameters(&SFormBandsums,&SOnlyReducingBandsums,
				&SDelete_Only_Short_Primitives, &SDo_Not_Reduce_Genus);
			H_Results = fopen("Heegaard_Results","a+");
			if(H_Results != NULL)			
				Set_Up_Simplification_Parameters_S1();								 		 	
			break;			
		case 's':
			Batch = 12;
			break;
		case 'S':
			SnapPy2Heegaard();
			goto END;			
		case 'u':
			Batch = 13;
			H_Results = fopen("Heegaard_Results","a+");
			if(H_Results != NULL) 
				fprintf(H_Results,"\n\nA stabilized version of each initial presentation should appear below.");			
			break;		
		case 'x':
			Batch = 14;
			printf("\n\nPRINT SORTED PRESENTATIONS FOUND BY DELETING PRIMITIVES ?  HIT 'y' OR 'n'.");
			GET_B14RESPONSE1:
			switch(WaitkbHit())
				{
				case 'y':
					B14B15PrintPres = TRUE;
					break;
				case 'n':
					B14B15PrintPres = FALSE;
					break;
				default:
					SysBeep(5);
					goto GET_B14RESPONSE1;
				}
			break;		
		case 'X':
			Batch = 15;
			printf("\n\nPRINT SORTED PRESENTATIONS FOUND BY DELETING PRIMITIVES ?  HIT 'y' OR 'n'.");
			H_Results = fopen("Heegaard_Results","a+");
			if(H_Results != NULL) 
				fprintf(H_Results,"\n\nThe number of primitive relators in each initial presentation should appear below.");			
			GET_B15RESPONSE1:
			switch(WaitkbHit())
				{
				case 'y':
					B14B15PrintPres = TRUE;
					break;
				case 'n':
					B14B15PrintPres = FALSE;
					break;
				default:
					SysBeep(5);
					goto GET_B15RESPONSE1;
				}			
			break;						   
		case 'z':
			Batch = 16;
			H_Results = fopen("Heegaard_Results","a+");
			if(H_Results != NULL) 
				fprintf(H_Results,"\n\nA minimal length version of each initial presentation should appear below.");			
			Turn_Micro_Print_On();
			SMicro_Print = Micro_Print;
			SMicro_Print_F = Micro_Print_F;	
			break;		
		default:
			SysBeep(5);
			goto GET_RESPONSE3;
		}
 
    fseek(input_relators,0L,0);
    Delete_Old_Presentations();
    NumPresExamined = 0;
    Stop = FALSE;        
	while(1)
		{
		if('s' == mykbhit()) 
			{
			printf("\n	Status: Processed %u presentations. Hit 'c' to continue. Hit 'q' to abort.",NumPresExamined);
			GET_RESPONSE2:			
			switch(WaitkbHit())
				{
				case 'c':
					break;
				case 'q':
					{
					fseek(input_relators,0L,0);
    				Delete_Old_Presentations();
    				Init_G_Variables();
    				NumFilled = 0;
					goto OPTIONS;		
					}
				default: goto GET_RESPONSE2;
				}	
			}
		
		if(Get_Next_Presentation_From_File()) break;
		
		/**************************************************************************************
            Echo the initial relators to the output so we will have a copy of them. Then call
            Freely_Reduce(), Rewrite_Input(), and Canonical_Rewrite() to get a presentation
            which serves as the initial presentation for Heegaard. 
        **************************************************************************************/    
        
        NumPresExamined ++;
        
        Micro_Print 	= FALSE;
       	Micro_Print_F 	= FALSE;
        for(i = 1,Scratch = 0L; i <= NumRelators; i++)
            Scratch += GetHandleSize((char **) Relators[i]);
        Scratch -= NumRelators;
        if(Batch != 3) printf("\n\nL %ld, ",Scratch);        									
        if(Freely_Reduce() == TOO_LONG)
            {
            printf("\n\nThis presentation is too long!!");
            SysBeep(5);
			continue;
            }
        if(Scratch > OrigLength)
            {
            printf("which freely reduces to the following presentation of length %lu.\n",
                OrigLength);
            Print_Relators(Relators,NumRelators);
            printf("\n\n");
            Scratch = OrigLength;
            if(Batch == 4 || Batch == 53)
            	{
            	printf("\nSince the original presentation is not freely reduced, it is not a Heegaard Splitting Rep!");
				if(H_Results != NULL) fprintf(H_Results,"\n\n%s <-- Not a HS Rep! (IP not freely reduced.)",PresName);	            	
            	continue;
            	}
            }   
        if(Rewrite_Input())
            {
            printf("\n\nThere must be at least one non-empty relator!!");
            SysBeep(5);
			continue;
            }   

        /**************************************************************************************
            Save a copy of the initial set of relators in case we want to refer to them later.
        **************************************************************************************/
        
        CopyNumRelators     = NumRelators;
        CopyNumGenerators   = NumGenerators;
        HS_Rep_Length 		= OrigLength;
        HS_Rep_NumGens		= NumGenerators;
        for(i = 1; i <= NumRelators; i++)
            {
            if(Copy_Of_Input[i] != NULL) DisposeHandle((char **) Copy_Of_Input[i]);
            Copy_Of_Input[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));
            if(Copy_Of_Input[i] == NULL) Mem_Error();
            p = *Copy_Of_Input[i];
            q = *Relators[i];
            while( (*p++ = *q++) ) ;
            }
                                        
        /**************************************************************************************
                Call Init_G_Variables() to initialize some global variables. Then call
                Canonical_Rewrite() to rewrite the presentation in canonical form.
        **************************************************************************************/
        
        Init_G_Variables();
        Length = Scratch;
        if(Batch != 3) printf("Gen %d, Rel %d.\n\n",NumGenerators,NumRelators);
        
        Micro_Print = TRUE;					
        if(Canonical_Rewrite(Relators,FALSE,FALSE) == TOO_LONG) printf("\n\nThis presentation has too many symmetries!!"); 
        Micro_Print = FALSE;					                         
        if(Compare_Input_Pres() == FALSE)
            {
            printf("\n\n The rewritten presentation is:");
            Print_Relators(Relators,NumRelators);
            printf("\n");
            if(Batch == 4 || Batch == 53)
            	{
           		printf("\nThe IP is not written in Canonical Form. Hence it is not a Heegaard Splitting Rep!");
				if(H_Results != NULL) fprintf(H_Results,"\n\n%s <-- Not a HS Rep! (IP not in Canonical Form.)",PresName);	            	
            	continue;            	
            	}
            }
        					
        /**************************************************************************************
                        Update the copy of the initial set of relators.
        **************************************************************************************/
        
        CopyNumRelators       = NumRelators;
        CopyNumGenerators     = NumGenerators;
        for(i = 1; i <= NumRelators; i++)
            {
            if(Copy_Of_Input[i] != NULL) DisposeHandle((char **) Copy_Of_Input[i]);
            Copy_Of_Input[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));
	        if(Copy_Of_Input[i] == NULL) Mem_Error();    
            p = *Copy_Of_Input[i];
            q = *Relators[i];
            while( (*p++ = *q++) ) ;
            }    
		
		switch(Batch)
			{
			case 1:
				if(NumGenerators != 2 || NumRelators != 1) continue;
				AlexanderPolynomial(*Relators[1]); 
				break;
			case 2:
				if(NumGenerators != 2 || NumRelators != 1) continue;
				Is_Knot_Relator();
				if(NumGenerators != 2 || NumRelators != 2) continue;
				Delete_Old_Presentations();
				CopyNumRelators       = NumRelators;
				CopyNumGenerators     = NumGenerators;
				for(i = 1; i <= NumRelators; i++)
					{
					if(Copy_Of_Input[i] != NULL) DisposeHandle((char **) Copy_Of_Input[i]);
					Copy_Of_Input[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Relators[i]));
					if(Copy_Of_Input[i] == NULL) Mem_Error();    
					p = *Copy_Of_Input[i];
					q = *Relators[i];
					while( (*p++ = *q++) ) ;
					}								
				Init_G_Variables();
				Micro_Print = TRUE;
				if(Canonical_Rewrite(Relators,FALSE,FALSE) == TOO_LONG) 
					printf("\n\nThis presentation has too many symmetries!!"); 
        		Micro_Print = FALSE;
				Batch = 10;				
				DepthFirstSearch				= TRUE;
				BreadthFirstSearch				= FALSE;
				Find_All_Min_Pres 				= FALSE;
				Delete_Only_Short_Primitives 	= FALSE;
    			Do_Not_Reduce_Genus 			= FALSE;
    			FormBandsums 					= TRUE;
    			OnlyReducingBandsums 			= TRUE; 
  				MyMaxSavedPres = 1000;  			
    			for(i = 1,Length = 0L; i <= NumRelators; i++)
					Length += GetHandleSize((char **) Relators[i]);
				Length -= NumRelators;				 			
    			if(Get_Initial_Diagram(TRUE)) 
    				{
    				Batch = 2;
    				continue;
    				}
    			TestRealizability1 = FALSE;
    			TestRealizability2 = FALSE;           
	   			if(Get_Diagrams() == INTERRUPT)
	   				{
	   				Stop = TRUE;
	   				break;
	   				}							
    			Sort_Presentations_In_Memory(1);
				Batch = 2;
				break;
			case 3:
				Check_Realizability_Of_The_Initial_Presentation();
				break;
			case 4:	
				Check_If_HS_Rep 	= TRUE;
				DepthFirstSearch	= FALSE;
				BreadthFirstSearch	= TRUE;
				Find_All_Min_Pres 	= FALSE;
				CheckHSReps 		= TRUE;
				Delete_Only_Short_Primitives 	= FALSE;
    			Do_Not_Reduce_Genus 			= FALSE;
    			FormBandsums 					= TRUE;
    			OnlyReducingBandsums 			= FALSE;
    			BPrintAnnulusData 				= FALSE;
    			BPrintNotConnectedData 			= FALSE;
    			switch(Just_Delete_Primitives(TRUE))
    				{
    				case 0:
    					break;
    				case 1:
    					printf(" Since this presentation contains primitives, it is not a HS Rep!");
						if(H_Results != NULL) fprintf(H_Results,"\n\n%s <-- Not a HS Rep!",PresName);    					
    					continue;
    				case 3:
    					continue;	
    				case INTERRUPT:
    					{
    					Stop = TRUE;
    					break;
    					}
    				default:
    					break;
    				}
    			NumFilled = 0;	
    			/* Reset the Relators */
                NumRelators = CopyNumRelators;
                if(NumRelators == 1) continue;
        		NumGenerators = CopyNumGenerators;
        		Vertices = 2*NumGenerators;
				for(i = 1; i <= NumRelators; i++)
					{
					if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
					Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Copy_Of_Input[i]));
					if(Relators[i] == NULL) Mem_Error();    
					p = *Copy_Of_Input[i];
					q = *Relators[i];
					while( (*q++ = *p++) ) ;
					}
    			if(Get_Initial_Diagram(TRUE)) continue;
    			TestRealizability1 = FALSE;
    			TestRealizability2 = FALSE;           
    			if(Get_Diagrams() == INTERRUPT) 
    				{
    				Stop = TRUE;			
    				break;
    				}
    			if(NumFilled >= MyMaxSavedPres) printf("\n%s is probably a HS Rep!",PresName);	
				if(H_Results != NULL && (NumFilled >= MyMaxSavedPres || OnStack == 0))
					fprintf(H_Results,"\n\n%s \nProbably a HS Rep!",PresName);  		
				break;
			case 5:
				Display_Diagram_Of_The_Initial_Presentation();
				break;
			case 6:
				Compute_Homology();
				break;
			case 7:
				if(Find_Level_Transformations_Of_The_Initial_Presentation() == 0)
				Test_LT_For_Pseudo_Min();
				DrawingDiagrams 		 = FALSE;
                Did_Exponent_Surgery 	 = FALSE;
                Did_Cutting_Disk_Surgery = FALSE;
				Delete_Old_Presentations();
				NumFilled = 0;
				break;				
			case 10:
				DepthFirstSearch	= TRUE;
				BreadthFirstSearch	= FALSE;
				Find_All_Min_Pres 	= FALSE;
				Delete_Only_Short_Primitives 	= SDelete_Only_Short_Primitives;
    			Do_Not_Reduce_Genus 			= SDo_Not_Reduce_Genus;
    			FormBandsums 					= SFormBandsums;
    			OnlyReducingBandsums 			= SOnlyReducingBandsums;
    			if(Get_Initial_Diagram(TRUE)) continue;
    			TestRealizability1 = FALSE;
    			TestRealizability2 = FALSE;           
	   			if(Get_Diagrams() == INTERRUPT) 
	   				{
	   				Stop = TRUE;
	   				break;
	   				}							
    			Sort_Presentations_In_Memory(1);	
    			break;    			
			case 11:
				DepthFirstSearch	= FALSE;
				BreadthFirstSearch	= TRUE;
				Find_All_Min_Pres 	= FALSE;
				Delete_Only_Short_Primitives 	= SDelete_Only_Short_Primitives;
    			Do_Not_Reduce_Genus 			= SDo_Not_Reduce_Genus;
    			FormBandsums 					= SFormBandsums;
    			OnlyReducingBandsums 			= SOnlyReducingBandsums;
    			if(Get_Initial_Diagram(TRUE)) continue;
    			TestRealizability1 = FALSE;
    			TestRealizability2 = FALSE;           
    			if(Get_Diagrams() == INTERRUPT) 
    				{
    				Stop = TRUE;
    				break;
    				}
    			Sort_Presentations_In_Memory(1);    							
				break;
			case 12:
				Find_Symmetries(FALSE);
				break;
			case 13:
				Stabilize();
				break;
			case 14:
    			if(Just_Delete_Primitives(FALSE) == INTERRUPT) 
    				{
    				Stop = TRUE;
    				break;
    				}
				break;
			case 15:
    			if(Just_Delete_Primitives(TRUE) == INTERRUPT) 
    				{
    				Stop = TRUE;
    				break;
    				}
				break;
			case 16:				
        		Micro_Print = SMicro_Print;
        		Micro_Print_F = SMicro_Print_F;		
				Reduce_The_Initial_Presentation_To_Minimal_Length();
				break;
			case 51:
				Fill_A(NumRelators);				
				if(ComputeValences_A()) break;
				Get_Matrix();
				for(i = 0; i < Vertices; i++) ZZ[i] = 0;
				if(Connected_(0,0) == FALSE) break;
				if(Sep_Pairs(0,0,1)) break;
				if(Planar(FALSE,TRUE) == TRUE) break;
				Diagram_Main();
				break;
			case 53:
				DepthFirstSearch	= FALSE;
				BreadthFirstSearch	= TRUE;
				Find_All_Min_Pres 	= FALSE;
				CheckHSReps 		= TRUE;
				Delete_Only_Short_Primitives 	= FALSE;
    			Do_Not_Reduce_Genus 			= FALSE;
    			FormBandsums 					= TRUE;
    			OnlyReducingBandsums 			= FALSE;
    			BPrintAnnulusData 				= FALSE;
    			BPrintNotConnectedData 			= FALSE;
    			switch(Just_Delete_Primitives(TRUE))
    				{
    				case 0:
    					break;
    				case 1:
    					printf(" Since this presentation contains primitives, it is not a HS Rep!");
						if(H_Results != NULL) fprintf(H_Results,"\n\n%s <-- Not a HS Rep! (IP contains primitives.)",PresName);    					
    					continue;
    				case 3:
    					continue;	
    				case INTERRUPT:
    					{
    					Stop = TRUE;
    					break;
    					}
    				default:
    					break;
    				}
    			NumFilled = 0;	
    			/* Reset the Relators */
                NumRelators = CopyNumRelators;
        		NumGenerators = CopyNumGenerators;
        		Vertices = 2*NumGenerators;
				for(i = 1; i <= NumRelators; i++)
					{
					if(Relators[i] != NULL) DisposeHandle((char **) Relators[i]);
					Relators[i] = (unsigned char **) NewHandle(GetHandleSize((char **) Copy_Of_Input[i]));
					if(Relators[i] == NULL) Mem_Error();    
					p = *Copy_Of_Input[i];
					q = *Relators[i];
					while( (*q++ = *p++) ) ;
					}
				if(Reduce_The_Initial_Presentation_To_Minimal_Length())
					{
					printf("\n\nThe IP does not have minimal length. Hence is not a HS Rep!");
					if(H_Results != NULL)
						fprintf(H_Results,"\n\n%s  <-- Not a HS Rep! (IP does not have minimal length.)",PresName);
					continue;	
					}			
				if(Stabilize()) continue;	
    			if(Get_Initial_Diagram(TRUE)) continue;	
    			TestRealizability1 = FALSE;
    			TestRealizability2 = FALSE;
    			printf("\n\nThe minimal length stabilized presentation is:");
    			Print_Relators(Relators,NumRelators);
	   			if(Get_Diagrams() == INTERRUPT)
    				{
    				Stop = TRUE;			
    				break;
    				} 
    			if(Sort_Presentations_In_Memory(FALSE) == INTERRUPT) Stop = TRUE;	 		
				break;							
			}
		if(Stop) break;	
		}
		
	printf("\n\n Totals: NumPresExamined %u",NumPresExamined);
	switch(Batch)
		{
		case 3: printf(", NumRealizable %u",NumRealizable);
			break;	
		}
	fclose(input_relators);
	if(H_Results != NULL) switch(Batch)
		{
		case 1:
			printf("\n\nAlexander polynomials of each 1-relator, 2-generator IP should appear in 'Heegaard_Results'.");
			fclose(H_Results);
			break;
		case 2:			
			printf("\n\n'Meridian reps of each realizable 1-relator, 2-generator IP should appear in 'Heegaard_Results'.");
			fclose(H_Results);
			break;		
		case 3:
			printf("\n\nEach IP's realizability should appear in 'Heegaard_Results'.");
			fclose(H_Results);
			break;	
		case 4:
			printf("\n\nIndication that the IP is a Heegaard Splitting Rep should appear in 'Heegaard_Results'.");
			fclose(H_Results);
			break;
		case 6:
			printf("\n\nEach IP's Z-homology should appear in 'Heegaard_Results'.");
			fclose(H_Results);
			break;
		case 7:
			printf("\n\nOrbit under level-transformations info for each IP should appear in 'Heegaard_Results'.");
			fclose(H_Results);
			break;
		case 10:
			printf("\n\nResults should appear in 'Heegaard_Results'.");
			fclose(H_Results);
			break;
		case 11:
			printf("\n\nResults should appear in 'Heegaard_Results'.");
			fclose(H_Results);
			break;										
		case 13:
			printf("\n\nStabilized IP's should appear in 'Heegaard_Results'.");
			fclose(H_Results);
			break;
		case 15:
			printf("\n\nThe number of primitives in each IP should appear in 'Heegaard_Results'.");
			fclose(H_Results);
			break;
		case 16:
			printf("\n\nMinimal length versions of each IP should appear in 'Heegaard_Results'.");
			fclose(H_Results);
			break;
		case 51:
			printf("\n\nThe dual relators of each presentation should appear in 'Heegaard_Results'.");
			fclose(H_Results);
			break;
		case 53:
			printf("\n\nHS Rep data for each IP should appear in 'Heegaard_Results'.");
			fclose(H_Results);
			break;						
		default:
			break;	
		}
END:		
	SysBeep(5);	
	printf("\n\nHIT 'B' TO CONTINUE IN 'BATCH' MODE. HIT 'q' TO QUIT RUNNING IN BATCH MODE.");
	GET_RESPONSE4:
	switch(WaitkbHit())
		{
		case 'B': 
			goto OPTIONS;
		case 'q':
			fseek(input_relators,0L,0);
			Delete_Old_Presentations();
			Init_G_Variables();
			NumFilled = 0;
			return(0);
		default:
			SysBeep(5);
			goto GET_RESPONSE4;
		}
}

int Set_Up_Simplification_Parameters(int* SFormBandsums, int* SOnlyReducingBandsums,
	int* SDelete_Only_Short_Primitives, int* SDo_Not_Reduce_Genus)
{

	printf("\n\nCREATE NEW DIAGRAMS BY FORMING BANDSUMS ? HIT 'y' OR 'n'.");
	GET_RESPONSE1:       
	switch(WaitkbHit())
		{
		case 'y':
			printf("\n    HIT 'a' TO FORM ALL POSSIBLE BANDSUMS.");
			printf("\n    HIT 'r' TO FORM ONLY LENGTH REDUCING BANDSUMS.");
			printf("\n    HIT 'n' TO FORM NO BANDSUMS.    ");
			GET_RESPONSE2:               
			switch(WaitkbHit())
				{
				case 'a':
					*SFormBandsums 			= TRUE;
					*SOnlyReducingBandsums 	= FALSE;
					break;
				case 'r':
					*SFormBandsums 			= TRUE;
					*SOnlyReducingBandsums 	= TRUE;
					break;
				case 'n':
					*SOnlyReducingBandsums = *SFormBandsums = FALSE;
					break;    
				default:
					SysBeep(5);
					goto GET_RESPONSE2;
				}
			printf("\n");    
			break;        
		case 'n':
			*SOnlyReducingBandsums = *SFormBandsums = FALSE;                
			break;
		default:
			SysBeep(5);
			goto GET_RESPONSE1;    
		} 
		
	printf("\n\nDELETE ALL PRIMITIVE RELATORS ? HIT 'y' OR 'n'.");
  
	printf("\n***************************************************************************");
	printf("\n        Be somewhat careful when choosing to delete primitives!");
	printf("\n    Although deleting primitives always preserves fundamental groups,");
	printf("\nsimplifying a presentation P, whose realizability has not been verified, by");
	printf("\ndeleting primitives, may falsely suggest P is realizable. For example,");
	printf("\ndeleting C from P = < A,B,C | AAACBBC, C > yields P' = < A,B | AAABB >,");
	printf("\nand P' is realizable, but P is not.");
	printf("\n***************************************************************************\n");
    
    GET_RESPONSE3:   
    switch(WaitkbHit())
        {
        case 'y':
            *SDo_Not_Reduce_Genus 			= FALSE;
            *SDelete_Only_Short_Primitives 	= FALSE;
            break;
        case 'n':    
            printf("\n\nDELETE PRIMITIVE RELATORS OF LENGTH 1 AND LENGTH 2 ? HIT 'y' OR 'n'.");
            GET_RESPONSE4:            
            switch(WaitkbHit())
                {
                case 'y':
                    *SDelete_Only_Short_Primitives 	= TRUE;
                    *SDo_Not_Reduce_Genus 			= FALSE;
                    break;
                case 'n':
                    *SDelete_Only_Short_Primitives 	= FALSE;
                    *SDo_Not_Reduce_Genus 			= TRUE;
                    break;
                default:
                    SysBeep(5);
                    goto GET_RESPONSE4;
                }
            break;
        default:
            SysBeep(5);
            goto GET_RESPONSE3;
        }
     
    printf("\n\n	SOME TERMINAL PRINTING OPTIONS.");    
    printf("\n1) PRINT INFO ABOUT ANNULI HEEGAARD FINDS ?");
    printf("\n2) PRINT INFO ABOUT DISCONNECTED DIAGRAMS HEEGAARD FINDS ?");
    
	printf("\n\n1) HIT 'y' OR 'n'.");    
    GET_RESPONSE5:
    switch(WaitkbHit())
    	{   	
		case 'y':
			BPrintAnnulusData = TRUE;
			break;
		case 'n':
			BPrintAnnulusData = FALSE;
			break;
		default:
			SysBeep(5);
			goto GET_RESPONSE5;
		} 
		
	printf("\n2) HIT 'y' OR 'n'.");    
	GET_RESPONSE6:
    switch(WaitkbHit())
    	{   	
		case 'y':
			BPrintNotConnectedData = TRUE;
			break;
		case 'n':
			BPrintNotConnectedData = FALSE;
			break;
		default:
			SysBeep(5);
			goto GET_RESPONSE6;
		}
		
	return(0);	   
}

int Set_Up_Simplification_Parameters_S1()
{
	printf("\n\n	SOME DATA SAVING OPTIONS.");
	printf("\n1) SAVE INFO ABOUT EACH RECOGNIZED MANIFOLD TO 'Heegaard_Results' ?");
	printf("\n2) SAVE INFO ABOUT EACH FINITE MANIFOLD TO 'Heegaard_Results' ?");
	printf("\n3) SAVE 'Heegaard Splitting Reps' FOUND TO 'Heegaard_Results ?");
	printf("\n4) SAVE PRESENTATIONS THAT 'split' AND THEIR COMPONENT PRESENTATIONS TO 'Heegaard_Results' ?");
	
	printf("\n\n1) HIT 'y' OR 'n'.");
	GET_RESPONSE1:       
	switch(WaitkbHit())
		{
		case 'y':
			B10B11Recognized 	= TRUE;
			B10B11Finite 		= FALSE;			                
			break;		
		case 'n':
			B10B11Recognized 	= FALSE;
			B10B11Finite		= FALSE;   
			break;        				
		default:
			SysBeep(5);
			goto GET_RESPONSE1;    
		}	

	if(B10B11Recognized == FALSE)
		{
		printf("\n2) HIT 'y' OR 'n'.");		
		GET_RESPONSE2:       
		switch(WaitkbHit())
			{
			case 'n':
				B10B11Finite	= FALSE;   
				break;        
			case 'y':
				B10B11Finite 	= TRUE;			                
				break;				
			default:
				SysBeep(5);
				goto GET_RESPONSE2;    
			}
		}
	else
		printf("\n	2) UNAVAILABLE WHEN REPLY TO 1) is 'y'.");	
	
    printf("\n3) HIT 'y' OR 'n'.");   
    GET_RESPONSE3:   
    switch(WaitkbHit())
        {
        case 'y':
			B10B11HSReps = TRUE;
            break;
        case 'n':    
			B10B11HSReps = FALSE;
            break;
        default:
            SysBeep(5);
            goto GET_RESPONSE3;
        }
        
    printf("\n4) HIT 'y' OR 'n'.");    
    GET_RESPONSE4:   
    switch(WaitkbHit())
        {
        case 'y':
			B10B11ConSumPres = TRUE;
            break;
        case 'n':    
			B10B11ConSumPres = FALSE;
            break;
        default:
            SysBeep(5);
            goto GET_RESPONSE4;
        }
                  		   	  
	SetLimits();
		
	return(0);	   
}

void SetLimits()
{
	printf("\n\nENTER A LIMIT <= %u FOR THE NUMBER OF PRESENTATIONS HEEGAARD SHOULD SAVE BEFORE MOVING ON.",
		MAX_SAVED_PRES - 3);   
		if(Batch == 53) 	
			printf("\n (Limits >= 10,000 may be required.)\n");
		else
			printf("\n (Smaller limits are faster.)\n"); 
	GET_RESPONSE:	  
	ReadString((char *)Inst, 200);
	sscanf((char *) Inst,"%u",&MyMaxSavedPres); 
	if(MyMaxSavedPres < 1 || MyMaxSavedPres > (MAX_SAVED_PRES - 3))
		{   	
		SysBeep(5);
		printf("\nMyMaxSavedPres is out of range. Please reenter it!");
		goto GET_RESPONSE;
		}
}				

int SnapPy2Heegaard()
{
    register unsigned char  *Ptr1,
    						*Ptr2,
    						*Ptr3,
    						*p,
                            *q,
                            *r; 
                            
   unsigned char			Str1[10] = "Relators:";                                     
    
    Ptr1 = (unsigned char*) NewPtr(3001);
	if(Ptr1 == NULL) Mem_Error();
	Ptr2 = (unsigned char*) NewPtr(3001);
	if(Ptr2 == NULL) Mem_Error();
	Ptr3 = (unsigned char*) NewPtr(3001);
	if(Ptr3 == NULL) Mem_Error();
	
	fseek(input_relators,0L,0);
	
	*Ptr1 = EOS;
	*Ptr2 = EOS;

GET_NEXT_PRES:	
	while(1)
		{
		if(fgets((char *) Ptr3,3000,input_relators) == NULL) return(1);
		/* Check if Ptr3 starts with the string "Relators:". */
		q = Ptr3;
		r = Str1;
		while(*r && (*q++ == *r++)) ;
		if(*r == EOS) goto FOUND_SUBSTRING_STR1; /* Ptr3 starts with the string "Relators:". */
		/* Otherwise copy Ptr2 into Ptr1, copy Ptr3 into Ptr2 and put the next line into Ptr3. */
		p = Ptr2;
		q = Ptr1;
		while( (*q++ = *p++) ) ;
		p = Ptr3;
		q = Ptr2;
		while( (*q++ = *p++) ) ;
		}
	
	FOUND_SUBSTRING_STR1:
	printf("\n%s",Ptr1);
	while(1)
		{
		if(fgets((char *) Ptr2,3000,input_relators) == NULL) return(1);
		if(*Ptr2 != ' ') goto GET_NEXT_PRES;
		if(*Ptr2 == ' ') printf("%s",Ptr2);
		}
	
	if(Ptr1) DisposePtr((char*) Ptr1);
	if(Ptr2) DisposePtr((char*) Ptr2);
	if(Ptr3) DisposePtr((char*) Ptr3);
	return(0);
}