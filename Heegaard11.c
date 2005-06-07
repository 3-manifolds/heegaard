#include "Heegaard.h"
#include "Heegaard_Dec.h"

int User_Says_Quit(void)
{
	unsigned char 	*r;
	
	unsigned long	newlimit;
	
GET_RESPONSE1:
	SysBeep(5);
	r = (unsigned char *) NewPtr(100);
	printf("\n\nMemory is getting rather low! It may be wise to save and quit!");	
	printf("\n\nThe program has used %lu K bytes of memory for storing presentations.",
		BytesUsed/1024);
	printf("\nThis exceeds the current limit of %lu K bytes reserved for this purpose.",
		BytesAvailable/1024);
	printf("\n\nHIT 'c' TO CONTINUE RUNNING THIS EXAMPLE.");	
	printf("\nHIT 't' TO TERMINATE THIS RUN.");
	printf("\nHIT 'v' TO REVIEW THE PRESENTATIONS NOW IN MEMORY.");
	GET_RESPONSE2:
	switch(WaitkbHit())
		{
		case 'c':
			GET_RESPONSE3:
			printf("\n\nENTER A NEW UPPER LIMIT ON THE AMOUNT OF MEMORY TO BE USED (in K) AND HIT 'return'.	");
			newlimit = 0L;
			ReadString((char *)r, GetPtrSize(r));
			sscanf((char *) r,"%lu",&newlimit);
			BytesAvailable = 1024*newlimit;
			if(BytesAvailable <= BytesUsed)
				{
				SysBeep(5);
				printf("\nThe number of bytes available must exceed the number already used!!");
				goto GET_RESPONSE3;
				}
			DisposePtr((char *) r);		
			return(FALSE);
		case 't':
			DisposePtr((char *) r);
			return(TRUE);
		case 'v':
			DisposePtr((char *) r);
			Report(Band_Sums,NumDiagrams,OnStack,Starting_Pres,0,0,1,0,1,0);
			goto GET_RESPONSE1;	
		default:
			SysBeep(5);
			goto GET_RESPONSE2;		
		}			
}

/**********************************************************************************************
	This file is a copy of the file mini.print.c "borrowed" from Lightspeed C and slightly
	modified.
**********************************************************************************************/
#ifdef MAC

#define topMargin 20
#define leftMargin 20
#define bottomMargin 20

#define NIL 0L

static	THPrint	hPrint = NIL;
static	int		tabWidth;

void CheckPrintHandle(void)
{
	if (hPrint == NIL)
		{
		hPrint = (TPrint **) NewHandle(sizeof(TPrint));
		PrintDefault(hPrint);
		PrValidate(hPrint);
		PrStlDialog(hPrint);
		}
}

void DoPageSetUp(void)
{
	PrOpen();
	CheckPrintHandle();
	if (!PrStlDialog(hPrint))
		{
		PrClose();
		PrDrvrClose();
		}
}

#define tabChar	((char)'\t')

void MyDrawText(char *p,int count)
{
	register char	*p1,
					*p2;
					
	int				len;

	p1 = p;
	p2 = p + count;
	while (p < p2)
		{
		while ((p1 < p2) && (*p1 != tabChar)) p1++;
		if ((len = p1 - p) > 0) DrawText(p,0,len);
		if (*p1 == tabChar)
			{
			Move(tabWidth,0);
			p1++;
			}
		p = p1;
		}
}

PrDoc(char ** hText,long count,THPrint hPrint,int font,int size)
{
	register int 	line = 0;
	register int 	lastLineOnPage = 0;
	int				length;
	int 			linesPerPage;
	int 			lineBase;
	int 			lineHeight;
	register char 	*ptr, *p1;
	FontInfo		info;
	Rect 			printRect;
	TPPrPort		printPort;

	printPort = PrOpenDoc(hPrint,0L,0L);
	SetPort((struct GrafPort *) printPort);
	TextFont(font);
	TextSize(size);
	printRect = (**hPrint).prInfo.rPage;
	GetFontInfo(&info);
	lineHeight = info.leading + info.ascent + info.descent;
	linesPerPage = 
		(printRect.bottom - printRect.top - topMargin - bottomMargin) / lineHeight;
	HLock(hText);
	ptr = p1 = (*hText);
	do
		{
		PrOpenPage(printPort,0L);
		lastLineOnPage += linesPerPage;
		MoveTo(printRect.left + leftMargin,(lineBase = printRect.top + lineHeight));
		do 
			{
			/* PrintLine: */
			while ((ptr <= (*hText) + count) && (*ptr++ != (char)'\r')) ;
			if ((length = (int)(ptr - p1) - 1) > 0) MyDrawText(p1,length);
			MoveTo(printRect.left + leftMargin,(lineBase += lineHeight));
			p1 = ptr;
			}
		while ((++line != lastLineOnPage) && (ptr <= (*hText) + count));
		PrClosePage(printPort);
		}
	while (ptr < (*hText) + count);
	HUnlock(hText);
	PrCloseDoc(printPort);
}

PrintText(char ** hText,long length,GrafPtr gp,int tabPixels)
{
	TPPrPort	printPort;
	GrafPtr		savePort;
	TPrStatus	prStatus;
	int			copies;
	
    PrOpen();
	CheckPrintHandle();
	tabWidth = tabPixels;
	if (PrJobDialog(hPrint) != 0)
		{
		GetPort(&savePort);
		for (copies = HowMany(); copies > 0; copies--)
			{
			PrDoc (hText,length,hPrint,(*gp).txFont,(*gp).txSize);
			PrPicFile(hPrint,0L,0L,0L,&prStatus);
			}
		SetPort(savePort);
		}
	PrClose();
	PrDrvrClose();
}

HowMany()
{
	return( ((**hPrint).prJob.bJDocLoop == bDraftLoop) ? 
				(**hPrint).prJob.iCopies : 1 );
}
#endif

int Old_Generate_Orbits_Under_Auts(void)
{
	int						i,
							j,
							Source;
	
	long					MaxAutLength;
	unsigned long			InitLength;
			
	register unsigned char 	*p,
							*q;
	
	unsigned char *			r;
							
	Compute_Stabilizers = FALSE;
	Left[0] = Right[0] = INFINITE;
		
	
	r = (unsigned char *) NewPtr(100);		
	printf("\n\nENTER THE MaxSavedLength OF A PRESENTATION.     ");
	ReadString((char *)r, GetPtrSize(r));
	sscanf((char *) r,"%ld",&MaxAutLength);		
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
		if(NumGenerators > 2) return(5);	
		for(i = 1; i <= NumRelators; i++)
			{
			ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[ReadPres][i]));
			if((q = *Relators[i]) == NULL) return(TOO_LONG);			
			p = *SUR[ReadPres][i];
			while(*q++ = *p++) ;
			}
		Length = SURL[ReadPres];
		
		ZZ[0] = FALSE;
		ZZ[1] = TRUE;
		ZZ[2] = FALSE;
		ZZ[3] = TRUE;
		Source = 0;
		if(Do_Aut(Source,1,NumRelators) == TOO_LONG) return(5);
		for(i = 1, Length = 0; i <= NumRelators; i++)
			Length += GetHandleSize((char **) Relators[i]);
			Length -= NumRelators;
		if(Length <= InitLength || (Length <= MaxAutLength && Length > SURL[ReadPres]))
			{	
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

		NumRelators 	= NR[ReadPres];
		NumGenerators 	= NG[ReadPres];
		if(NumGenerators > 2) return(5);	
		for(i = 1; i <= NumRelators; i++)
			{
			ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[ReadPres][i]));
			if((q = *Relators[i]) == NULL) return(TOO_LONG);			
			p = *SUR[ReadPres][i];
			while(*q++ = *p++) ;
			}
		Length = SURL[ReadPres];
			
		ZZ[0] = FALSE;
		ZZ[1] = TRUE;
		ZZ[2] = TRUE;
		ZZ[3] = FALSE;
		Source = 0;
		if(Do_Aut(Source,1,NumRelators) == TOO_LONG) return(5);
		for(i = 1, Length = 0; i <= NumRelators; i++)
			Length += GetHandleSize((char **) Relators[i]);
			Length -= NumRelators;
		if(Length <= InitLength || (Length <= MaxAutLength && Length > SURL[ReadPres]))
			{	
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
 
		NumRelators 	= NR[ReadPres];
		NumGenerators 	= NG[ReadPres];
		if(NumGenerators > 2) return(5);	
		for(i = 1; i <= NumRelators; i++)
			{
			ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[ReadPres][i]));
			if((q = *Relators[i]) == NULL) return(TOO_LONG);			
			p = *SUR[ReadPres][i];
			while(*q++ = *p++) ;
			}
		Length = SURL[ReadPres]; 
 			
		ZZ[0] = FALSE;
		ZZ[1] = TRUE;
		ZZ[2] = FALSE;
		ZZ[3] = TRUE;
		Source = 2;
		if(Do_Aut(Source,1,NumRelators) == TOO_LONG) return(5);
		for(i = 1, Length = 0; i <= NumRelators; i++)
			Length += GetHandleSize((char **) Relators[i]);
			Length -= NumRelators;
		if(Length <= InitLength || (Length <= MaxAutLength && Length > SURL[ReadPres]))
			{	
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

		NumRelators 	= NR[ReadPres];
		NumGenerators 	= NG[ReadPres];
		if(NumGenerators > 2) return(5);	
		for(i = 1; i <= NumRelators; i++)
			{
			ReallocateHandle((char **) Relators[i],GetHandleSize((char **) SUR[ReadPres][i]));
			if((q = *Relators[i]) == NULL) return(TOO_LONG);			
			p = *SUR[ReadPres][i];
			while(*q++ = *p++) ;
			}
		Length = SURL[ReadPres];
			
		ZZ[0] = TRUE;
		ZZ[1] = FALSE;
		ZZ[2] = FALSE;
		ZZ[3] = TRUE;
		Source = 2;
		if(Do_Aut(Source,1,NumRelators) == TOO_LONG) return(5);
		for(i = 1, Length = 0; i <= NumRelators; i++)
			Length += GetHandleSize((char **) Relators[i]);
			Length -= NumRelators;
		if(Length <= InitLength || (Length <= MaxAutLength && Length > SURL[ReadPres]))
			{	
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
								
		ReadPres++;
		}
	while(ReadPres < NumFilled);
	
	return(0);	
}

