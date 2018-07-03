/*////////////////////////////////////////////////////////////////
Permission to use, copy, modify, and distribute this program for
any purpose, with or without fee, is hereby granted, provided that
the notices on the head, the reference information, and this
copyright notice appear in all copies or substantial portions of
the Software. It is provided "as is" without express or implied
warranty.
*////////////////////////////////////////////////////////////////
// ProteinSurface.cpp: implementation of the ProteinSurface class.
//
//////////////////////////////////////////////////////////////////////
#include "ProteinSurface.h"

#define PURE 0
#define ATOM 1
#define CHAIN 2

#define VWS 1
#define MS 2
#define SAS 3
#define SES 4

#define INNERANDOUTER 0
#define OUTER 1
#define INNER 2

#define X 0
#define Y 1
#define Z 2

#define A 0
#define B 1
#define C 2

#define IX 0
#define IY 1
#define IZ 2


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ProteinSurface::ProteinSurface()
{
	int i;
	boxLength=128;
	flagRadius=false;
	scaleFactor=1;
	probeRadius=1.4;
	for(i=0;i<13;i++)
		deptY[i]=NULL;
	volumePixels=NULL;
	//Remember, this is a 3D array
	pHeight=0;
	pWidth=0;
	pLength=0;

	faces=NULL;
	verts=NULL;
	faceNumber=0;
	vertexNumber=0;
	surfaceArea=0;
	surfaceVolume=0;
	cavityArea=0;
	cavityVolume=0;
	numberOfCavities=0;
	eachCavityVolume=NULL;
	fixSf=4;
	depthValues=NULL;
}

ProteinSurface::~ProteinSurface()
{
	int i,j;
	for(i=0;i<13;i++)
	{	
		if(deptY[i]!=NULL)
		{
			free(deptY[i]);
			deptY[i]=NULL;
		}
	}
	if(volumePixels!=NULL)
	{  
		for(i=0;i<pLength;i++)
		{   
			for(j=0;j<pWidth;j++)
				free(volumePixels[i][j]);
			free(volumePixels[i]);
		}
   		free(volumePixels);
		volumePixels=NULL;
	}
	if(faces!=NULL)
	{
		free(faces);
		faces=NULL;
	}
	if(verts!=NULL)
	{
		free(verts);
		verts=NULL;
	}
	if(eachCavityVolume!=NULL)
	{
		free(eachCavityVolume);
		eachCavityVolume=NULL;
	}
	if(depthValues!=NULL)
	{
		free(depthValues);
		depthValues=NULL;
	}
}
void ProteinSurface::calcDepth(int numberOfBackbones,int seqInit,int seqTerm, atom* proSeq, boneInfo *backbones, bool type)
{
	/**
	*Gets the depth for each atom in the sequence
	*@param numberOfBackbones Number of elements in array backbones
	*@param seqInit Starting place within the array to iterate over
	*@param seqTerm Ending place within the array to iterate over
	*@param proSeq The array of atoms for which distances are to be found
	*@param backbones An array of bone infos
	*@param type Tells us which type  
	**/

	//?Thread this?

	int i,tIndex;
	int ox,oy,oz;
	Vector3d cp;
	if(depthValues)
	{
		delete[]depthValues;
	}

	if(!type){
		//If type = false, we ignore the seqInit and seqTerm we were provided with
		//And iterate across the whole array proSeq
		numberOfDepths = numberOfBackbones;
		seqInit = 0;
		seqTerm = numberOfBackbones - 1;
	}

	for(i = seqInit; i <= seqTerm; i++){

		if(!type){
			tIndex=backbones[i].indca; 
		}
		else{
			tIndex = i;
		}

		cp = (proSeq[tIndex].coordinates + pTran) * scaleFactor;
		//Translating and scaling

		ox=int(cp(X)+0.5);
		oy=int(cp(Y)+0.5);
		oz=int(cp(Z)+0.5);

		//Rounding to the nearest integer

		depthValues[i-seqInit]=volumePixels[ox][oy][oz].distance/scaleFactor-probeRadius;
		//Set the depth value at that index equal to the corresponding atom's associated distance
		//Divided by the scale factor, minus the probe radius

		if(depthValues[i-seqInit]<rasRad[proSeq[tIndex].detail]){
			depthValues[i-seqInit]=rasRad[proSeq[tIndex].detail];
		}

		depthValues[i-seqInit]+=probeRadius;

	}
}

void ProteinSurface::checkEuler()
{
	/**
	* Check the Euler characteristic of the surface
	* Operates on the mesh
	**/

	//!Thread this!

	int i,j,k;
	int ia,ib,ic;
	int *vertexDegrees[4][20];//0 a 1 b 2 face 3 cutend]
	bool *vertexFlag=new bool[int(vertexNumber*1.1)];
	int *vertGroup=new int[int(vertexNumber*1.1)];
	int *vertexNumbers=new int[int(vertexNumber*1.1)];

	for(i=0;i<4;i++)
	{
		for(j=0;j<20;j++)
			vertexDegrees[i][j]=new int[int(vertexNumber*1.1)];
	}

	for(j=0;j<vertexNumber;j++)
	{
		vertexNumbers[j]=0;
		vertGroup[j]=0;
		vertexFlag[j]=true;
	}

	bool *flagFace=new bool[faceNumber];
	//degree of each vert

	for(i=0;i<faceNumber;i++)
	{
		ia=faces[i].abc(A);
		ib=faces[i].abc(B);
		ic=faces[i].abc(C);

		vertexDegrees[0][vertexNumbers[ia]][ia]=ib;
		vertexDegrees[1][vertexNumbers[ia]][ia]=ic;
		vertexDegrees[2][vertexNumbers[ia]][ia]=i;
		vertexDegrees[0][vertexNumbers[ib]][ib]=ic;
		vertexDegrees[1][vertexNumbers[ib]][ib]=ia;
		vertexDegrees[2][vertexNumbers[ib]][ib]=i;
		vertexDegrees[0][vertexNumbers[ic]][ic]=ia;
		vertexDegrees[1][vertexNumbers[ic]][ic]=ib;
		vertexDegrees[2][vertexNumbers[ic]][ic]=i;

		vertexNumbers[ia]++;
		vertexNumbers[ib]++;
		vertexNumbers[ic]++;
		
		flagFace[i]=true;	
	}//i

	int jb,jc;
	int kb,kc;
	int l,m, n;

	vertexInfo *dupVert;

	int numberOfDuplicates=0;
	int allocDup=20;

	dupVert=new vertInfo[allocDup];

	bool flagDup;
	int *tpIndex[3];

	for(i=0;i<3;i++)
	{
		tpIndex[i]=new int[20];
	}

	for(i=0;i<vertexNumber;i++)
	{
		//remove dup faces
		for(j=0;j<vertexNumbers[i]-1;j++)
		{
			jb=vertexDegrees[0][j][i];
			jc=vertexDegrees[1][j][i];
			for(k=j+1;k<vertexNumbers[i];k++)
			{
				kb=vertexDegrees[0][k][i];
				kc=vertexDegrees[1][k][i];
				if(jb==kc && jc==kb)
				{
		//			printf("%d dup face %d %d [%d %d %d]\n",i,vertexDegrees[2][j][i],vertexDegrees[2][k][i],i,jb,jc);
					flagFace[vertexDegrees[2][j][i]]=false;
					flagFace[vertexDegrees[2][k][i]]=false;


					for(l = j; l < vertexNumbers[i] - 2; l++){
						for(m = 0; m < 3; m++){
							if(l < k - 1){
								vertexDegrees[m][l][i]=vertexDegrees[m][l+1][i];
							}
							else{
								vertexDegrees[m][l][i]=vertexDegrees[m][l+2][i];								
							}
						}
					}

					j--;
					k=vertexNumbers[i];
					vertexNumbers[i]-=2;
				}//duplicate
				else if(jb==kb && jc==kc)
				{
		//			printf("wrong same faces %d %d\n",vertexDegrees[2][j][i],vertexDegrees[2][k][i]);
				}
			}//k
		}//j
		if(vertexNumbers[i]==0)
		{
	//		printf("no use vertex %d\n",i);
			vertexFlag[i]=false;
			continue;
		}
/*		else if(vertexNumbers[i]==1 || vertexNumbers[i]==2)
		{
		//	printf("single vertex %d %d \n",i,vertexNumber[i]);
		//	vertexFlag[i]=false;
		}
*/
		//reorder 
		flagDup=false;
		for(j=0;j<vertexNumbers[i]-1;j++)
		{
			for(k=j+1;k<vertexNumbers[i];k++)
			{
				if(vertexDegrees[0][j][i]==vertexDegrees[0][k][i])
				{
					flagDup=true;
					break;
				}
			}
			if(flagDup)
				break;
		}
		if(flagDup)
		{
			for(k = 0; k < vertexNumbers[i]; k++){
				for(l = 0; l < 3; l++){
					if(k < j){
						tpIndex[l][vertexNumbers[i]-j+k]=vertexDegrees[l][k][i];
					}
					else{
						tpIndex[l][vertexNumbers[i]-j+k]=vertexDegrees[l][k][i];
					}
					vertexDegrees[l][k][i]=tpIndex[l][k];
				}//End l loop
			}//End k loop
		}//End if flagDup
		//arrage all faces around a vert
		j=0;	
		while(j<vertexNumbers[i])//start cycle
		{		
			jb=vertexDegrees[0][j][i];
			jc=vertexDegrees[1][j][i];
			m=j;
			do{//find m+1
				k=vertexNumbers[i];
				for(k=m+1;k<vertexNumbers[i];k++)
				{
					if(vertexDegrees[0][k][i]==jc)
						break;
				}
				if(k<vertexNumbers[i])
				{
					if(k!=m+1)
					{
						for(n = 0; n < 3; n++){
							l=vertexDegrees[n][m+1][i];
							vertexDegrees[n][m+1][i]=vertexDegrees[n][k][i];
							vertexDegrees[n][k][i]=l;
						}
					}
					jc=vertexDegrees[1][m+1][i];	
					m++;
				}
				else
				{
					break;
				}
			}while(jc!=jb && m < vertexNumbers[i]);
			if(jc==jb)//one cycle
			{
				vertexDegrees[3][vertGroup[i]][i]=m;
				vertGroup[i]++;
			}
			else//single
			{
		//		printf("no corre index %d %d\n",i,jc);
				vertexDegrees[3][vertGroup[i]][i]=m;
				vertGroup[i]++;		
/*				for(j=0;j<vertexNumber[i];j++)
				{
		//			printf("detail %d %d %d %d %d\n",j,vertexDegrees[0][j][i],vertexDegrees[1][j][i],vertexDegrees[2][j][i],vertexDegrees[3][j][i]);
				}
*/			}
			j=m+1;
		}//while
		if(vertGroup[i]!=1)
		{
	//		printf("split vert %d %d\n",i,vertGroup[i]);
	//		for(j=0;j<vertexNumber[i];j++)
	//		{
	//			printf("%d %d %d %d %d\n",j,vertexDegrees[0][j][i],vertexDegrees[1][j][i],vertexDegrees[2][j][i],vertexDegrees[3][j][i]);
	//		}
			if(numberOfDuplicates+vertGroup[i]>allocDup)
			{
				allocDup*=2;
				dupVert=(vertexInfo *)realloc(dupVert,allocDup*sizeof(vertexInfo));
			}
			
			for(j=1;j<vertGroup[i];j++)
			{	
				dupVert[numberOfDuplicates]=verts[i];
				vertexFlag[numberOfDuplicates+vertexNumber]=true;
				vertGroup[numberOfDuplicates+vertexNumber]=1;
				vertexNumbers[numberOfDuplicates+vertexNumber]=vertexDegrees[3][j][i]-vertexDegrees[3][j-1][i];
				for(k=0;k<vertexNumbers[numberOfDuplicates+vertexNumber];k++)
				{
					vertexDegrees[0][k][numberOfDuplicates+vertexNumber]=vertexDegrees[0][vertexDegrees[3][j-1][i]+k+1][numberOfDuplicates+vertexNumber];
					vertexDegrees[1][k][numberOfDuplicates+vertexNumber]=vertexDegrees[1][vertexDegrees[3][j-1][i]+k+1][numberOfDuplicates+vertexNumber];
					vertexDegrees[2][k][numberOfDuplicates+vertexNumber]=vertexDegrees[2][vertexDegrees[3][j-1][i]+k+1][numberOfDuplicates+vertexNumber];
				}
				for(k=vertexDegrees[3][j-1][i]+1;k<=vertexDegrees[3][j][i];k++)
				{
		//			printf("changing %d %d\n",k,vertexDegrees[2][k][i]);
					if(faces[vertexDegrees[2][k][i]].abc(A)==i)
					{
						faces[vertexDegrees[2][k][i]].abc(A)=numberOfDuplicates+vertexNumber;
						m=faces[vertexDegrees[2][k][i]].abc(B);
						for(l=0;l<vertexNumbers[m];l++)
						{
							if(vertexDegrees[2][l][m]==vertexDegrees[2][k][i])
							{
								if(vertexDegrees[0][l][m]==i)
								{
									vertexDegrees[0][l][m]=numberOfDuplicates+vertexNumber;
								}
								else if(vertexDegrees[1][l][m]==i)
								{
									vertexDegrees[1][l][m]=numberOfDuplicates+vertexNumber;
								}
								else 
								{
							//		printf("wrong modified vertab %d\n",m);
								}
							}
						}
						m=faces[vertexDegrees[2][k][i]].abc(C);
						for(l=0;l<vertexNumbers[m];l++)
						{
							if(vertexDegrees[2][l][m]==vertexDegrees[2][k][i])
							{
								if(vertexDegrees[0][l][m]==i)
								{
									vertexDegrees[0][l][m]=numberOfDuplicates+vertexNumber;
								}
								else if(vertexDegrees[1][l][m]==i)
								{
									vertexDegrees[1][l][m]=numberOfDuplicates+vertexNumber;
								}
								else 
								{
								//	printf("wrong modified vertac %d\n",m);
								}
							}
						}
					}
					else if(faces[vertexDegrees[2][k][i]].abc(B)==i)
					{
						faces[vertexDegrees[2][k][i]].abc(B)=numberOfDuplicates+vertexNumber;
						m=faces[vertexDegrees[2][k][i]].abc(A);
						for(l=0;l<vertexNumbers[m];l++)
						{
							if(vertexDegrees[2][l][m]==vertexDegrees[2][k][i])
							{
								if(vertexDegrees[0][l][m]==i)
								{
									vertexDegrees[0][l][m]=numberOfDuplicates+vertexNumber;
								}
								else if(vertexDegrees[1][l][m]==i)
								{
									vertexDegrees[1][l][m]=numberOfDuplicates+vertexNumber;
								}
								else 
								{
							//		printf("wrong modified vertba %d\n",m);
								}
							}
						}
						m=faces[vertexDegrees[2][k][i]].abc(C);
						for(l=0;l<vertexNumbers[m];l++)
						{
							if(vertexDegrees[2][l][m]==vertexDegrees[2][k][i])
							{
								if(vertexDegrees[0][l][m]==i)
								{
									vertexDegrees[0][l][m]=numberOfDuplicates+vertexNumber;
								}
								else if(vertexDegrees[1][l][m]==i)
								{
									vertexDegrees[1][l][m]=numberOfDuplicates+vertexNumber;
								}
								else 
								{
							//		printf("wrong modified vertbc %d\n",m);
								}
							}
						}
					}
					else if(faces[vertexDegrees[2][k][i]].abc(C)==i)
					{
						faces[vertexDegrees[2][k][i]].abc(C)=numberOfDuplicates+vertexNumber;
						
						m=faces[vertexDegrees[2][k][i]].abc(A);
						
						for(l=0;l<vertexNumbers[m];l++)
						{
							if(vertexDegrees[2][l][m]==vertexDegrees[2][k][i])
							{
								if(vertexDegrees[0][l][m]==i)
								{
									vertexDegrees[0][l][m]=numberOfDuplicates+vertexNumber;
								}
								else if(vertexDegrees[1][l][m]==i)
								{
									vertexDegrees[1][l][m]=numberOfDuplicates+vertexNumber;
								}
								else 
								{
							//		printf("wrong modified vertca %d\n",m);
								}
							}
						}
						
						m=faces[vertexDegrees[2][k][i]].abc(B);
						
						for(l=0;l<vertexNumbers[m];l++)
						{
							if(vertexDegrees[2][l][m]==vertexDegrees[2][k][i])
							{
								if(vertexDegrees[0][l][m]==i)
								{
									vertexDegrees[0][l][m]=numberOfDuplicates+vertexNumber;
								}
								else if(vertexDegrees[1][l][m]==i)
								{
									vertexDegrees[1][l][m]=numberOfDuplicates+vertexNumber;
								}
								else 
								{
						//			printf("wrong modified vertcb %d\n",m);
								}
							}
						}
					}
					else
					{
					//	printf("wrong vert %d face %d [%d %d %d]\n",i,vertexDegrees[2][k][i],faces[vertexDegrees[2][k][i]].a,
					//		faces[vertexDegrees[2][k][i]].b,faces[vertexDegrees[2][k][i]].c);
					}
				}//k
				vertGroup[i]=1;
				numberOfDuplicates++;
			}//j  
		}//if need to split
	}//i
	for(i=0;i<3;i++)
	{
		delete[]tpIndex[i];
	}
	//reduce face
	int totFaceNew=0;
	for(i=0;i<faceNumber;i++)
	{
		if(flagFace[i]){
			faces[totFaceNew] = faces[i];
			totFaceNew++;
		}
	}
//	printf("number faces from %d to %d\n",faceNumber,totFaceNew);
	faceNumber=totFaceNew;
	
	//new points
	int totVertNew=0;
	int *vertexIndex2=new int[vertexNumber+numberOfDuplicates];
	verts=(vertInfo*)realloc(verts,(vertexNumber+numberOfDuplicates)*sizeof(vertInfo));
	for(i=0;i<vertexNumber;i++)
	{
		if(vertexFlag[i]){
			vertexIndex2[i] = totVertNew;
			verts[totVertNew] = verts[i];
			totVertNew++;
		}
	}

	for(i=0;i<numberOfDuplicates;i++)
	{
		vertexIndex2[vertexNumber+i]=totVertNew;
		verts[totVertNew]=dupVert[i];
		totVertNew++;
	}

	for(i=0;i<faceNumber;i++)
	{
		faces[i].abc(A)=vertexIndex2[faces[i].abc(A)];
		faces[i].abc(B)=vertexIndex2[faces[i].abc(B)];
		faces[i].abc(C)=vertexIndex2[faces[i].abc(C)];
	}

	delete[]vertexIndex2;
//	printf("number verts from %d to %d (new added %d)\n",vertexNumber,totVertNew,numberOfDuplicates);
	vertexNumber=totVertNew;
	if((2*vertexNumber-faceNumber)%4!=0)
	printf("euler num %d\n",2*vertexNumber-faceNumber);//comp+cav-genus

	//release
	for(i=0;i<4;i++)
	{
		for(j=0;j<20;j++)
			delete[]vertexDegrees[i][j];
	}
	delete[]vertGroup;
	delete[]vertexFlag;
	delete[]vertexNumbers;
	delete[]flagFace;
	delete[]dupVert;
}
void ProteinSurface::laplacianSmooth(int numberOfIterations){
	/**
	*Does laplacian smoothing on the surface
	*Also works on the mesh
	*@param numberOfIterations 
	**/
	Vector3d *tps=new Vector3d[vertexNumber];
	int *vertexDegrees[20];
	int i,j;
	bool flagVert;
	for(i=0;i<20;i++)
	{
		vertexDegrees[i]=new int[vertexNumber];		
	}
	for(i=0;i<vertexNumber;i++)
	{
		vertexDegrees[0][i]=0;
	}
	for(i=0;i<faceNumber;i++)
	{
		//a
		flagVert=true;
		for(j=0;j<vertexDegrees[0][faces[i].abc(A)];j++)
		{
			if(faces[i].abc(B)==vertexDegrees[j+1][faces[i].abc(A)])
			{
				flagVert=false;
				break;
			}
		}
		if(flagVert)
		{
			vertexDegrees[0][faces[i].abc(A)]++;
			vertexDegrees[vertexDegrees[0][faces[i].abc(A)]][faces[i].abc(A)]=faces[i].abc(B);
			
		}
		flagVert=true;
		for(j=0;j<vertexDegrees[0][faces[i].abc(A)];j++)
		{
			if(faces[i].abc(C)==vertexDegrees[j+1][faces[i].abc(A)])
			{
				flagVert=false;
				break;
			}
		}
		if(flagVert)
		{
			vertexDegrees[0][faces[i].abc(A)]++;
			vertexDegrees[vertexDegrees[0][faces[i].abc(A)]][faces[i].abc(A)]=faces[i].abc(C);
			
		}
		//b
		flagVert=true;
		for(j=0;j<vertexDegrees[0][faces[i].abc(B)];j++)
		{
			if(faces[i].abc(A)==vertexDegrees[j+1][faces[i].abc(B)])
			{
				flagVert=false;
				break;
			}
		}
		if(flagVert)
		{
			vertexDegrees[0][faces[i].abc(B)]++;
			vertexDegrees[vertexDegrees[0][faces[i].abc(B)]][faces[i].abc(B)]=faces[i].abc(A);
			
		}
		flagVert=true;
		for(j=0;j<vertexDegrees[0][faces[i].abc(B)];j++)
		{
			if(faces[i].abc(C)==vertexDegrees[j+1][faces[i].abc(B)])
			{
				flagVert=false;
				break;
			}
		}
		if(flagVert)
		{
			vertexDegrees[0][faces[i].abc(B)]++;
			vertexDegrees[vertexDegrees[0][faces[i].abc(B)]][faces[i].abc(B)]=faces[i].abc(C);
			
		}
		//c
		flagVert=true;
		for(j=0;j<vertexDegrees[0][faces[i].abc(C)];j++)
		{
			if(faces[i].abc(A)==vertexDegrees[j+1][faces[i].abc(C)])
			{
				flagVert=false;
				break;
			}
		}
		if(flagVert)
		{
			vertexDegrees[0][faces[i].abc(C)]++;
			vertexDegrees[vertexDegrees[0][faces[i].abc(C)]][faces[i].abc(C)]=faces[i].abc(A);
			
		}
		flagVert=true;
		for(j=0;j<vertexDegrees[0][faces[i].abc(C)];j++)
		{
			if(faces[i].abc(B)==vertexDegrees[j+1][faces[i].abc(C)])
			{
				flagVert=false;
				break;
			}
		}
		if(flagVert)
		{
			vertexDegrees[0][faces[i].abc(C)]++;
			vertexDegrees[vertexDegrees[0][faces[i].abc(C)]][faces[i].abc(C)]=faces[i].abc(B);
			
		}
	}
	
	double wt=1.00;
	double wt2=0.50;
	int sSign;
	int k;
	double outwt=0.75/(scaleFactor+3.5);//area-preserving

	for(k=0;k<numberOfIterations;k++){

		for(i=0;i<vertexNumber;i++){

			if(vertexDegrees[0][i]<3){
				tps[i] = verts[i].xyz;
			}//End if i < 3

			else if(vertexDegrees[0][i] >= 3){
				int weight = wt;

				if(vertexDegrees[0][i] <= 4){
					weight = wt2;
				}//End if i <=4

				tps[i](X)=0;
				tps[i](Y)=0;
				tps[i](Z)=0;
	
				for(j=0;j<vertexDegrees[0][i];j++){
					tps[i]+=verts[vertexDegrees[j+1][i]].xyz;
				}//End for loop j

				tps[i]+= weight*verts[i].xyz;
		
				tps[i]/=float(weight+vertexDegrees[0][i]);
			}//End if i >= 3
		}//End for loop i

		for(i=0;i<vertexNumber;i++){
			verts[i].xyz = tps[i];
		}//End for loop i

		computeNorm();

		for(i=0;i<vertexNumber;i++){
			if(verts[i].inOut) sSign=1;
			else sSign=-1;
			verts[i].xyz += sSign*outwt*verts[i].pn;
		}//End for loop i
	}
	delete[]tps;
	for(i=0;i<20;i++)
		delete[]vertexDegrees[i];
}

////////////////////////////////////////////////////////////////////////true: heavy atoms
void ProteinSurface::boundBox(int seqInit,int seqTerm,atom* proSeq,bool atomType,
							  Vector3d *minPoint,Vector3d *maxPoint)
{
	/**
	*Finds the bound box of the sequence of atoms
	*The smallest rectangular prism that contains all the atoms
	*@param seqInit Starting place within the array to iterate over
	*@param seqTerm Ending place within the array to iterate over
	*@param proSeq The array of atoms for which distances are to be found
	*@param atomType True if atoms are heavy, false if they are light
	*@param minPoint A pointer to a vector representing the minimum point
	*@param maxPoint A pointer to a vector representing the maximum point
	**/

	int i;
	Vector3d minimumPoint = *minPoint;
	Vector3d maximumPoint = *maxPoint;

	minimumPoint(X)=100000;minimumPoint(Y)=100000;minimumPoint(Z)=100000;
	maximumPoint(X)=-100000;maximumPoint(Y)=-100000;maximumPoint(Z)=-100000;
	
	for(i=seqInit;i<=seqTerm;i++)
	{
		if(proSeq[i].simpleType==1 && proSeq[i].ins==' ')
		{
			if(atomType && (proSeq[i].detail==5 ||proSeq[i].detail==12))
				continue;
			if(proSeq[i].coordinates(X)<minimumPoint(X))
				minimumPoint(X)=proSeq[i].coordinates(X);
			if(proSeq[i].coordinates(Y)<minimumPoint(Y))
				minimumPoint(Y)=proSeq[i].coordinates(Y);
			if(proSeq[i].coordinates(Z)<minimumPoint(Z))
				minimumPoint(Z)=proSeq[i].coordinates(Z);
			if(proSeq[i].coordinates(X)>maximumPoint(X))
				maximumPoint(X)=proSeq[i].coordinates(X);
			if(proSeq[i].coordinates(Y)>maximumPoint(Y))
				maximumPoint(Y)=proSeq[i].coordinates(Y);
			if(proSeq[i].coordinates(Z)>maximumPoint(Z))
				maximumPoint(Z)=proSeq[i].coordinates(Z);
		}
	}

}
//label of atom
void ProteinSurface::atomsInOut(int seqInit,int seqTerm,atom* proSeq)
{
	/**
	*@param seqInit Starting place within the array to iterate over
	*@param seqTerm Ending place within the array to iterate over
	*@param proSeq The array of atoms for which distances are to be found
	**/

	int i;
	for(i=seqInit;i<=seqTerm;i++)
	{
		proSeq[i].inOut=2;//between
	}
	for(i=0;i<vertexNumber;i++)
	{
		if(!verts[i].isCont)
			proSeq[verts[i].atomId].inOut=1;//outer
	}
	for(i=0;i<faceNumber;i++)
	{
		if(faces[i].inOut)
		{
			proSeq[verts[faces[i].abc(A)].atomId].inOut=3;//inner
			proSeq[verts[faces[i].abc(B)].atomId].inOut=3;
			proSeq[verts[faces[i].abc(C)].atomId].inOut=3;
		}
	}
}
void ProteinSurface::checkInOutPropa()
{
	int i,j,ia,ib,ic;
	int *vertexDegrees[40];
	for(i=0;i<40;i++)
	{
		vertexDegrees[i]=new int[vertexNumber];
		memset(vertexDegrees[i],0,vertexNumber*sizeof(int));
	}
	for(i=0;i<faceNumber;i++)
	{
		ia=faces[i].abc(A);
		ib=faces[i].abc(B);
		ic=faces[i].abc(C);
		vertexDegrees[vertexDegrees[39][ia]][ia]=ib;
		vertexDegrees[39][ia]++;
		vertexDegrees[vertexDegrees[39][ia]][ia]=ic;
		vertexDegrees[39][ia]++;
		vertexDegrees[vertexDegrees[39][ib]][ib]=ic;
		vertexDegrees[39][ib]++;
		vertexDegrees[vertexDegrees[39][ib]][ib]=ia;
		vertexDegrees[39][ib]++;
		vertexDegrees[vertexDegrees[39][ic]][ic]=ia;
		vertexDegrees[39][ic]++;	
		vertexDegrees[vertexDegrees[39][ic]][ic]=ib;
		vertexDegrees[39][ic]++;	
	}//i
	int *arrIn=new int[vertexNumber];
	int nArrIn,narrOut;
	int *arrOut=new int[vertexNumber];
	bool *flagPoints=new bool[vertexNumber];
	memset(flagPoints,false,vertexNumber*sizeof(bool));
	nArrIn=3;
	arrIn[0]=faces[0].abc(A);
	arrIn[1]=faces[0].abc(B);
	arrIn[2]=faces[0].abc(C);
	int numberOfIterations=0;
	while(nArrIn>0)
	{
	//	printf("do iteration %d %d\n",numberOfIterations,nArrIn);
		narrOut=0;
		for(i=0;i<nArrIn;i++) if(!flagPoints[arrIn[i]])
		{
			flagPoints[arrIn[i]]=true;
			for(j=0;j<vertexDegrees[39][arrIn[i]];j++) if(!flagPoints[vertexDegrees[j][arrIn[i]]])
			{
				arrOut[narrOut]=vertexDegrees[j][arrIn[i]];
				narrOut++;
			}
		}
		memcpy(arrIn,arrOut,vertexNumber*sizeof(int));
		nArrIn=narrOut;
		numberOfIterations++;
	}
	for(i=0;i<vertexNumber;i++)
	{
		if(flagPoints[i]) verts[i].inOut=true;
		else verts[i].inOut=false;
	}
	for(i=0;i<faceNumber;i++)
	{
		if(flagPoints[faces[i].abc(A)] || flagPoints[faces[i].abc(B)] || flagPoints[faces[i].abc(C)]) faces[i].inOut=false;
		else faces[i].inOut=true;
	}
	for(i=0;i<40;i++)
	{
		delete[]vertexDegrees[i];
	}
	delete[]flagPoints;
	delete[]arrIn;
	delete[]arrOut;
}
void ProteinSurface::outputPly(char *fileName,atom* proSeq,int nColor,int tInOut)
{//0 both 1 outer 2 inner
	/**
	*Outputs the surface as a ply file
	*@param fileName, the name of the file to be read in
	*@param proSeq The array of atoms for which distances are to be found
	*@param nColor An integer representing the color mode - 0 if pure, 1 if atom, 2 if chain
	*@param tInOut An integer to specify which surface is to be rendered - 0 if both, 1 if outer, 2 if inner
	**/
	int i;
	unsigned char chainColor[256];
	int tChain,indColor2=1;
	for(i=0;i<256;i++)
	{
		chainColor[i]=0;
	}
	FILE *file;
	file=fopen(fileName,"wt");
	if(file==NULL)
	{
		printf("wrong to output ply file %s\n",fileName);
		return;
	}
int *vertexIndex=new int[vertexNumber];
bool *flagFace=new bool[faceNumber];
int realVertexNumber,realFaceNum;
memset(vertexIndex,-1,vertexNumber*sizeof(int));
memset(flagFace,false,faceNumber*sizeof(bool));
if(tInOut==INNERANDOUTER)
{
	realVertexNumber=vertexNumber;
	for(i=0;i<vertexNumber;i++)
	{
		vertexIndex[i]=i;
	}
	realFaceNum=faceNumber;
	memset(flagFace,true,faceNumber*sizeof(bool));
}
else if(tInOut==OUTER)
{
	realVertexNumber=0;
	for(i=0;i<vertexNumber;i++)
	{
		if(verts[i].inOut) 
		{
			vertexIndex[i]=realVertexNumber;
			realVertexNumber++;
		}
	}
	realFaceNum=0;
	for(i=0;i<faceNumber;i++)
	{
		if(vertexIndex[faces[i].abc(A)]!=-1 && vertexIndex[faces[i].abc(B)]!=-1 && vertexIndex[faces[i].abc(C)]!=-1)
		{
			realFaceNum++;
			flagFace[i]=true;
		}
	}
}
else if(tInOut==INNER)
{
	realVertexNumber=0;
	for(i=0;i<vertexNumber;i++)
	{
		if(!verts[i].inOut) 
		{
			vertexIndex[i]=realVertexNumber;
			realVertexNumber++;
		}
	}
	realFaceNum=0;
	for(i=0;i<faceNumber;i++)
	{
		if(vertexIndex[faces[i].abc(A)]!=-1 && vertexIndex[faces[i].abc(B)]!=-1 && vertexIndex[faces[i].abc(C)]!=-1)
		{
			realFaceNum++;
			flagFace[i]=true;
		}
	}
}
	fprintf(file,"ply\n" );
	fprintf(file,"format ascii 1.0\n" );
	fprintf(file,"comment ball mesh\n" );
	fprintf(file,"element vertex %d\n", realVertexNumber);
	fprintf(file,"property float x\n" );
	fprintf(file,"property float y\n" );
	fprintf(file,"property float z\n" );
	fprintf(file,"property uchar red\n" );
	fprintf(file,"property uchar green\n" );
	fprintf(file,"property uchar blue\n" );
	fprintf(file,"element face %d\n",realFaceNum );
	fprintf(file,"property list uchar int vertex_indices\n" );
	fprintf(file,"property uchar red\n" );
	fprintf(file,"property uchar green\n" );
	fprintf(file,"property uchar blue\n" );
	fprintf(file,"end_header\n" );
	int tColor[3];
    for(i=0;i<vertexNumber;i++) if(vertexIndex[i]!=-1)
	{ 
		if(nColor==PURE)//pure
		{
			tColor[0]=myColor[13][0];
			tColor[1]=myColor[13][1];
			tColor[2]=myColor[13][2];
		}
		else if(nColor==ATOM)//atom
		{
			tColor[0]=color[proSeq[verts[i].atomId].detail][0];
			tColor[1]=color[proSeq[verts[i].atomId].detail][1];
			tColor[2]=color[proSeq[verts[i].atomId].detail][2];
			if(verts[i].isCont)
			{
				tColor[0]=color[14][0];
				tColor[1]=color[14][1];
				tColor[2]=color[14][2];
			}
		}
		else if(nColor==CHAIN)//chain
		{
			if(chainColor[proSeq[verts[i].atomId].chainId]==0)
			{
				chainColor[proSeq[verts[i].atomId].chainId]=indColor2;
				indColor2++;
				if(indColor2>18)
					indColor2=1;
			}
			tColor[0]=myColor[chainColor[proSeq[verts[i].atomId].chainId]][0];
			tColor[1]=myColor[chainColor[proSeq[verts[i].atomId].chainId]][1];
			tColor[2]=myColor[chainColor[proSeq[verts[i].atomId].chainId]][2];
		}


		fprintf(file,"%.3f %.3f %.3f %3d %3d %3d\n", verts[i].xyz(X)/scaleFactor-pTran(X),verts[i].xyz(Y)/scaleFactor-pTran(Y),
			verts[i].xyz(Z)/scaleFactor-pTran(Z),tColor[0],tColor[1],tColor[2]);	
	}
	for(i=0;i<faceNumber;i++) if(flagFace[i])
	{ 
		if(nColor==PURE)//pure
		{
			tColor[0]=myColor[13][0];
			tColor[1]=myColor[13][1];
			tColor[2]=myColor[13][2];
		}
		else if(nColor==ATOM)//atom
		{
			if(proSeq[verts[faces[i].abc(B)].atomId].detail==proSeq[verts[faces[i].abc(C)].atomId].detail)
			{
				tColor[0]=color[proSeq[verts[faces[i].abc(B)].atomId].detail][0];
				tColor[1]=color[proSeq[verts[faces[i].abc(B)].atomId].detail][1];
				tColor[2]=color[proSeq[verts[faces[i].abc(B)].atomId].detail][2];
			}
			else
			{
				tColor[0]=color[proSeq[verts[faces[i].abc(A)].atomId].detail][0];
				tColor[1]=color[proSeq[verts[faces[i].abc(A)].atomId].detail][1];
				tColor[2]=color[proSeq[verts[faces[i].abc(A)].atomId].detail][2];
			}
		}
		else if(nColor==CHAIN)//chain
		{
			if(chainColor[proSeq[verts[faces[i].abc(A)].atomId].chainId]==0)
			{
				chainColor[proSeq[verts[faces[i].abc(A)].atomId].chainId]=indColor2;
				indColor2++;
				if(indColor2>18)
					indColor2=1;
			}
			if(chainColor[proSeq[verts[faces[i].abc(B)].atomId].chainId]==0)
			{
				chainColor[proSeq[verts[faces[i].abc(B)].atomId].chainId]=indColor2;
				indColor2++;
				if(indColor2>18)
					indColor2=1;
			}
			if(chainColor[proSeq[verts[faces[i].abc(C)].atomId].chainId]==0)
			{
				chainColor[proSeq[verts[faces[i].abc(C)].atomId].chainId]=indColor2;
				indColor2++;
				if(indColor2>18)
					indColor2=1;
			}
			if(proSeq[verts[faces[i].abc(B)].atomId].chainId==proSeq[verts[faces[i].abc(C)].atomId].chainId)
				tChain=chainColor[proSeq[verts[faces[i].abc(B)].atomId].chainId];
			else 
				tChain=chainColor[proSeq[verts[faces[i].abc(A)].atomId].chainId];
			tColor[0]=myColor[tChain][0];
			tColor[1]=myColor[tChain][1];
			tColor[2]=myColor[tChain][2];
		}
		if(!faces[i].inOut)//outer
			fprintf(file,"3 %d %d %d %3d %3d %3d\n", vertexIndex[faces[i].abc(A)],vertexIndex[faces[i].abc(B)],vertexIndex[faces[i].abc(C)],tColor[0],tColor[1],tColor[2]);
		else
			fprintf(file,"3 %d %d %d %3d %3d %3d\n", vertexIndex[faces[i].abc(A)],vertexIndex[faces[i].abc(C)],vertexIndex[faces[i].abc(B)],tColor[0],tColor[1],tColor[2]);
	}		
	fclose(file);
	delete[]vertexIndex;
	delete[]flagFace;
}
void ProteinSurface::outputOff(char *fileName) 
{
	// TODO: Add your command handler code here
	FILE  *stream;
	int i;

    if( (stream =fopen(fileName,"w+t"))== NULL)
	{
		printf("wrong to output off file %s\n",fileName);
		return;
	}
	fprintf(stream,"off\n" );
	fprintf(stream,"%d  %d  0\n",vertexNumber,faceNumber);
    for(i=0;i<vertexNumber;i++)
	{ 
		fprintf(stream,"%.6f %.6f %.6f\n", verts[i].xyz(X)/scaleFactor-pTran(X),
			verts[i].xyz(Y)/scaleFactor-pTran(Y),
			verts[i].xyz(Z)/scaleFactor-pTran(Z));	
	}
	for(i=0;i<faceNumber;i++)
	{ 
		fprintf(stream,"3 %d %d %d\n", faces[i].abc(A),faces[i].abc(C),faces[i].abc(B));		
	}		
    fclose(stream);
}
void ProteinSurface::buildBoundary()
{
	int i,j,k;
	int ii;
	bool flagBound;

	//!Thread this!

	for(i=0;i<pLength;i++)
	{
		for(j=0;j<pHeight;j++)
		{
			for(k=0;k<pWidth;k++)
			{
				if(volumePixels[i][k][j].inOut)
				{
					//6 neighbors
//					if(( k-1>-1 && !volumePixels[i][k-1][j].inOut) || ( k+1<pWidth &&!volumePixels[i][k+1][j].inOut)
//					|| ( j-1>-1 && !volumePixels[i][k][j-1].inOut) || ( j+1<pHeight &&!volumePixels[i][k][j+1].inOut)
//					|| ( i-1>-1 && !volumePixels[i-1][k][j].inOut) || ( i+1<pLength &&!volumePixels[i+1][k][j].inOut))
//						volumePixels[i][k][j].isBound=true;
			//	/*
					//26 neighbors
					flagBound=false;
					ii=0;
					while(!flagBound && ii<26)
					{
						if(i+nb[ii][0]>-1 && i+nb[ii][0]<pLength
							&& k+nb[ii][1]>-1 && k+nb[ii][1]<pWidth
							&& j+nb[ii][2]>-1 && j+nb[ii][2]<pHeight
							&& !volumePixels[i+nb[ii][0]][k+nb[ii][1]][j+nb[ii][2]].inOut)
						{
							volumePixels[i][k][j].isBound=true;
							flagBound=true;
						}
						else ii++;
					}
			//		*/
				}			
			}
		
		}
	}
}
//determine each voxel
void ProteinSurface::surfaceInterior()
{
	int i,j,k;
	voxel2 *inArray,*outArray;
//	int allocIn=1000;
//	int allocOut=1000;
	int inNum=0;
	int outNum;
	int totalNumber=0;

	//!Thread this!

	for(i=0;i<pLength;i++)
	{
		for(j=0;j<pWidth;j++)
		{
			for(k=0;k<pHeight;k++)
			{
				volumePixels[i][j][k].inOut=true;
				volumePixels[i][j][k].isBound=false;//has put into array
				if(!volumePixels[i][j][k].isDone)
					totalNumber++;
			}
		}
	}
	inArray=new voxel2[totalNumber];
	outArray=new voxel2[totalNumber];
	if(!volumePixels[0][0][0].isDone)
	{
		inArray[inNum].ix=0;
		inArray[inNum].iy=0;
		inArray[inNum].iz=0;
		inNum++;
	}
	else if(!volumePixels[pLength-1][0][0].isDone)
	{
		inArray[inNum].ix=pLength-1;
		inArray[inNum].iy=0;
		inArray[inNum].iz=0;
		inNum++;
	}
	else if(!volumePixels[0][pWidth-1][0].isDone)
	{
		inArray[inNum].ix=0;
		inArray[inNum].iy=pWidth-1;
		inArray[inNum].iz=0;
		inNum++;
	}
	else if(!volumePixels[0][0][pHeight-1].isDone)
	{
		inArray[inNum].ix=0;
		inArray[inNum].iy=0;
		inArray[inNum].iz=pHeight-1;
		inNum++;
	}
	else if(!volumePixels[pLength-1][pWidth-1][0].isDone)
	{
		inArray[inNum].ix=pLength-1;
		inArray[inNum].iy=pWidth-1;
		inArray[inNum].iz=0;
		inNum++;
	}
	else if(!volumePixels[0][pWidth-1][pHeight-1].isDone)
	{
		inArray[inNum].ix=0;
		inArray[inNum].iy=pWidth-1;
		inArray[inNum].iz=pHeight-1;
		inNum++;
	}
	else if(!volumePixels[pLength-1][0][pHeight-1].isDone)
	{
		inArray[inNum].ix=pLength-1;
		inArray[inNum].iy=0;
		inArray[inNum].iz=pHeight-1;
		inNum++;
	}
	else if(!volumePixels[pLength-1][pWidth-1][pHeight-1].isDone)
	{
		inArray[inNum].ix=pLength-1;
		inArray[inNum].iy=pWidth-1;
		inArray[inNum].iz=pHeight-1;
		inNum++;
	}
	else
	{
		printf("No init point\n");
		return;
	}
	volumePixels[inArray[0].ix][inArray[0].iy][inArray[0].iz].isBound=true;
	int ii,jj,kk;
	voxel2 *tArray,tnv;
	while(inNum!=0)
	{
		outNum=0;
		for(i=0;i<inNum;i++)
		{
			volumePixels[inArray[i].ix][inArray[i].iy][inArray[i].iz].inOut=false;
			for(ii=-1;ii<2;ii++)
			{
				for(jj=-1;jj<2;jj++)
				{
					for(kk=-1;kk<2;kk++)
					{
						tnv.ix=inArray[i].ix+ii;
						tnv.iy=inArray[i].iy+jj;
						tnv.iz=inArray[i].iz+kk;
						if(tnv.ix>-1 && tnv.ix<pLength 
							&& tnv.iy>-1 && tnv.iy<pWidth
							&& tnv.iz>-1 && tnv.iz<pHeight
							&& volumePixels[tnv.ix][tnv.iy][tnv.iz].inOut
							&& !volumePixels[tnv.ix][tnv.iy][tnv.iz].isDone
							&& !volumePixels[tnv.ix][tnv.iy][tnv.iz].isBound)
						{
							volumePixels[tnv.ix][tnv.iy][tnv.iz].isBound=true;
							outArray[outNum].ix=tnv.ix;
							outArray[outNum].iy=tnv.iy;
							outArray[outNum].iz=tnv.iz;
							outNum++;
						}
					}//kk
				}//jj
			}//ii
		}
		tArray=outArray;
		outArray=inArray;
		inArray=tArray;
		inNum=outNum;
	}//while
	delete[]inArray;
	delete[]outArray;
}
void ProteinSurface::outSas(int numberOfBackbones,boneInfo *backbones,atom* proSeq,char *fileName)
{
	double *areaSum=new double[numberOfBackbones];
	int i,j;
	int indRes;
	for(i=0;i<numberOfBackbones;i++)
	{
		areaSum[i]=0;
	}
	for(i=0;i<faceNumber;i++)
	{
		if(!faces[i].inOut)
		{
			for(j=0;j<numberOfBackbones;j++)
			{
				if(verts[faces[i].abc(A)].atomId>=backbones[j].iStart
					&& verts[faces[i].abc(A)].atomId<=backbones[j].iEnd)
				{
					indRes=j;
					break;
				}
			}	
			areaSum[indRes]+=faces[i].area/3.0/scaleFactor/scaleFactor;
			for(j=0;j<numberOfBackbones;j++)
			{
				if(verts[faces[i].abc(B)].atomId>=backbones[j].iStart
					&& verts[faces[i].abc(B)].atomId<=backbones[j].iEnd)
				{
					indRes=j;
					break;
				}
			}	
			areaSum[indRes]+=faces[i].area/3.0/scaleFactor/scaleFactor;
			for(j=0;j<numberOfBackbones;j++)
			{
				if(verts[faces[i].abc(C)].atomId>=backbones[j].iStart
					&& verts[faces[i].abc(C)].atomId<=backbones[j].iEnd)
				{
					indRes=j;
					break;
				}
			}	
			areaSum[indRes]+=faces[i].area/3.0/scaleFactor/scaleFactor;
		}
	}
	double tVal;
	FILE  *stream;
	stream =fopen(fileName,"w+t");
	fprintf(stream,"indx res area   solve #chn #res\n");
	for(i=0;i<numberOfBackbones;i++)
	{
		if(backbones[i].indca!=-1 && proSeq[backbones[i].indca].simpleType==1)
		{
			tVal=areaSum[i]/asArea2[proSeq[backbones[i].indca].residueId];
			if(tVal>1.0) tVal=1.0;
			fprintf(stream,"%4d %c %8.3f %5.3f %c    %4d\n",
				i+1,backbones[i].residueId,areaSum[i],tVal,proSeq[backbones[i].indca].chainId,proSeq[backbones[i].indca].residueNumber);
		}
		else
			fprintf(stderr,"lack ca atom %d/%d\n",i,numberOfBackbones);
	}
	fclose(stream);

	delete[]areaSum;
}

void ProteinSurface::outputCavityAtoms(int seqInit,int seqTerm,atom* proSeq,char *fileName)
{
	FILE  *stream;
    if( (stream =fopen(fileName,"w+t"))== NULL)
	{
		printf("Wrong to write to file %s\n", fileName);
		return;
	}
	int i,j;
	int ii,jj,kk;
	int iii,jjj,kkk;
	int si,sj,sk;
	bool *flagSeq;
	flagSeq=new bool[seqTerm+1];
	for(j=0;j<numberOfCavities;j++)
	{
	//	fprintf(stream,"%2d:     volume--%8.4f\n",j,eachCavityVolume[j]);
		for(i=0;i<seqTerm+1;i++)
		{
			flagSeq[i]=false;
		}
		for(ii=0;ii<pLength;ii++)
		{
			for(jj=0;jj<pWidth;jj++)
			{
				for(kk=0;kk<pHeight;kk++)
				{
					if(!volumePixels[ii][jj][kk].isDone && volumePixels[ii][jj][kk].inOut 
						&& int(volumePixels[ii][jj][kk].distance)==j)
					{
						for(iii=-1;iii<2;iii++)
						{
							for(jjj=0;jjj<2;jjj++)
							{
								for(kkk=0;kkk<2;kkk++)
								{
									si=ii+iii;
									sj=jj+jjj;
									sk=kk+kkk;
									if(si>-1 && si<pLength && 
										sj>-1 && sj<pWidth && 
										sk>-1 && sk<pHeight && 
										volumePixels[si][sj][sk].isDone &&
										!flagSeq[volumePixels[si][sj][sk].atomId])
									{
										flagSeq[volumePixels[si][sj][kk+kkk].atomId]=true;
									/*	fprintf(stream,"%2d:   %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %c%4d\n", j,
											proSeq[volumePixels[si][sj][sk].atomId].sequenceNumber,
											proSeq[volumePixels[si][sj][sk].atomId].detailType[0],
											proSeq[volumePixels[si][sj][sk].atomId].detailType[1],
											proSeq[volumePixels[si][sj][sk].atomId].detailType[2],
											proSeq[volumePixels[si][sj][sk].atomId].detailType[3],
											proSeq[volumePixels[si][sj][sk].atomId].residue[0],
											proSeq[volumePixels[si][sj][sk].atomId].residue[1],
											proSeq[volumePixels[si][sj][sk].atomId].residue[2],
											proSeq[volumePixels[si][sj][sk].atomId].chainId,
											proSeq[volumePixels[si][sj][sk].atomId].residueNumber,
											proSeq[volumePixels[si][sj][sk].atomId].xyz(X),
											proSeq[volumePixels[si][sj][sk].atomId].xyz(Y),
											proSeq[volumePixels[si][sj][sk].atomId].xyz(Z),
											aad1[proSeq[volumePixels[si][sj][sk].atomId].residueId],
											proSeq[volumePixels[si][sj][sk].atomId].residueNumber
											);*/
										fprintf(stream,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %c %2d %8.3f\n", 
											proSeq[volumePixels[si][sj][sk].atomId].sequenceNumber,
											proSeq[volumePixels[si][sj][sk].atomId].detailType[3],
											proSeq[volumePixels[si][sj][sk].atomId].detailType[0],
											proSeq[volumePixels[si][sj][sk].atomId].detailType[1],
											proSeq[volumePixels[si][sj][sk].atomId].detailType[2],
											
											proSeq[volumePixels[si][sj][sk].atomId].residue[0],
											proSeq[volumePixels[si][sj][sk].atomId].residue[1],
											proSeq[volumePixels[si][sj][sk].atomId].residue[2],
											proSeq[volumePixels[si][sj][sk].atomId].chainId,
											proSeq[volumePixels[si][sj][sk].atomId].residueNumber,
											proSeq[volumePixels[si][sj][sk].atomId].coordinates(X),
											proSeq[volumePixels[si][sj][sk].atomId].coordinates(Y),
											proSeq[volumePixels[si][sj][sk].atomId].coordinates(Z),
											aad1[proSeq[volumePixels[si][sj][sk].atomId].residueId],j,eachCavityVolume[j]
											);
									}//if
								}
							}
						}
					}//if
				}
			}
		}//ii
	
	}//j
	delete[] flagSeq;
	fclose(stream);
	
}
void ProteinSurface::boundingAtom(bool bType)
{
	int i,j,k;
	double tRadius[13];
	double txz,tDept,sRadius;
	int indx;
	for(i=0;i<13;i++)
	{	
		if(deptY[i]!=NULL)
			 free(deptY[i]);
	}
	flagRadius=bType;
	for(i=0;i<13;i++)
	{
		if(bType==false)
			tRadius[i]=rasRad[i]*scaleFactor+0.5;
		else 
			tRadius[i]=(rasRad[i]+probeRadius)*scaleFactor+0.5;

		sRadius=tRadius[i]*tRadius[i];
		widXz[i]=int(tRadius[i])+1;
		deptY[i]=new int[widXz[i]*widXz[i]];
		indx=0;
		for(j=0;j<widXz[i];j++)
		{
			for(k=0;k<widXz[i];k++)
			{
				txz=j*j+k*k;
				if(txz>sRadius)
				{
					deptY[i][indx]=-1;
				}
				else
				{
					tDept=sqrt(sRadius-txz);
					deptY[i][indx]=int(tDept+0.0); 
				}
				indx++;
			}
		}
		
	}

}
void ProteinSurface::initPara(int seqInit,int seqTerm,atom* proSeq,bool atomType,bool bType)
{
	int i,j;
	double fMargin=2.5;
	if(volumePixels!=NULL)
	{  
		for(i=0;i<pLength;i++)
		{   
			for(j=0;j<pWidth;j++)
				free(volumePixels[i][j]);
		}
		for(i=0;i<pLength;i++)
			free(volumePixels[i]);
		free(volumePixels);
		volumePixels=NULL;
	}
	boundBox(seqInit,seqTerm,proSeq,atomType,&pMin,&pMax);
	if(bType==false)
	{
		pMin(X)-=fMargin;
		pMin(Y)-=fMargin;
		pMin(Z)-=fMargin;
		pMax(X)+=fMargin;
		pMax(Y)+=fMargin;
		pMax(Z)+=fMargin;
	}
	else
	{
		pMin(X)-=probeRadius+fMargin;
		pMin(Y)-=probeRadius+fMargin;
		pMin(Z)-=probeRadius+fMargin;
		pMax(X)+=probeRadius+fMargin;
		pMax(Y)+=probeRadius+fMargin;
		pMax(Z)+=probeRadius+fMargin;
	}

	pTran(X)=-pMin(X);
	pTran(Y)=-pMin(Y);
	pTran(Z)=-pMin(Z);
	scaleFactor=pMax(X)-pMin(X);
	if((pMax(Y)-pMin(Y))>scaleFactor)
		scaleFactor=pMax(Y)-pMin(Y);
	if((pMax(Z)-pMin(Z))>scaleFactor)
		scaleFactor=pMax(Z)-pMin(Z);
	scaleFactor=(boxLength-1.0)/double(scaleFactor);
	///////////////////////////add this automatically first fix sf then fix boxLength
//	/*
	boxLength=int(boxLength*fixSf/scaleFactor);
	scaleFactor=fixSf;
	double threshBox=300;
	if(boxLength>threshBox)
	{
		double sfThresh=threshBox/double(boxLength);
		boxLength=int(threshBox);
		scaleFactor=scaleFactor*sfThresh;
	}
//	*/

	pLength=int(ceil(scaleFactor*(pMax(X)-pMin(X)))+1);
	pWidth=int(ceil(scaleFactor*(pMax(Y)-pMin(Y)))+1);
	pHeight=int(ceil(scaleFactor*(pMax(Z)-pMin(Z)))+1);
	if(pLength>boxLength)
		pLength=boxLength;
	if(pWidth>boxLength)
		pWidth=boxLength;
	if(pHeight>boxLength)
		pHeight=boxLength;
	boundingAtom(bType);
	cutRadius=probeRadius*scaleFactor;
}
void ProteinSurface::fillAtom(int indx,atom* proSeq,bool bColor)
{
	int cx,cy,cz;
	int ox,oy,oz;
	Vector3d cp;

	cp = (proSeq[indx].coordinates + pTran) * scaleFactor;
	cx=int(cp(X)+0.5);
	cy=int(cp(Y)+0.5);
	cz=int(cp(Z)+0.5);
	/*Translate, scale, and round coordinates*/

	int at=proSeq[indx].detail;
	int i,j,k;
	int ii,jj,kk;
	int mi,mj,mk;
	int si,sj,sk;
	int tIndex;
	int nIndex=0;
	for(i=0;i<widXz[at];i++)
	{
		for(j=0;j<widXz[at];j++)
		{
			if(deptY[at][nIndex]!=-1)
			{
					
				for( ii=-1;ii<2;ii++)
				{
					for( jj=-1;jj<2;jj++)
					{
						for( kk=-1;kk<2;kk++)
						{
							if(ii!=0 && jj!=0 && kk!=0)
							{		
								mi=ii*i;
								mk=kk*j;
								for(k=0;k<=deptY[at][nIndex];k++)
								{
									mj=k*jj;
									si=cx+mi;
									sj=cy+mj;
									sk=cz+mk;
									if(si<0 || sj<0 || sk<0 || si>=pLength || sj>=pWidth || sk>=pHeight)
									{
										continue;
									}
								if(!bColor)
									{
										volumePixels[si][sj][sk].inOut=true;
										continue;
									}
								else{
									if(volumePixels[si][sj][sk].inOut==false)
									{
										volumePixels[si][sj][sk].inOut=true;
										volumePixels[si][sj][sk].atomId=indx;
									}
									//no atomic info to each voxel change above line
								//*	
								    else if(volumePixels[si][sj][sk].inOut)
									{
										tIndex=volumePixels[si][sj][sk].atomId;

										cp = (proSeq[tIndex].coordinates + pTran) * scaleFactor;
										//Translating and scaling

										ox=int(cp(X)+0.5)-si;
										oy=int(cp(Y)+0.5)-sj;
										oz=int(cp(Z)+0.5)-sk;
										//Rounding to the nearest integer

										if(mi*mi+mj*mj+mk*mk<ox*ox+oy*oy+oz*oz)
											volumePixels[si][sj][sk].atomId=indx;
									}
								//	*/
								}//k
									}//else
							}//if
						}//kk	
					}//jj
				}//ii
			
				
			}//if
			nIndex++;
		}//j
	}//i
}
//sas use inOut
void ProteinSurface::fillVoxels(int seqInit,int seqTerm,bool atomType,atom* proSeq,bool bColor)
{

	int i,j,k;
	if(volumePixels==NULL)
	{
		volumePixels=new volumePixel**[pLength];
		for(i=0;i<pLength;i++)
			volumePixels[i]=new volumePixel*[pWidth];
		for(i=0;i<pLength;i++)
		{   
			for(j=0;j<pWidth;j++)
			{
				volumePixels[i][j]=new volumePixel[pHeight];
			}
		}
	}
	
	for(i=0;i<pLength;i++)
	{   
		for(j=0;j<pWidth;j++)
		{
			for(k=0;k<pHeight;k++)
			{
				volumePixels[i][j][k].inOut=false;
				volumePixels[i][j][k].isDone=false;
				volumePixels[i][j][k].isBound=false;
				volumePixels[i][j][k].distance=-1;
				volumePixels[i][j][k].atomId=-1;
			}
		}	
	}
//	int totalNumber=0;
	for(i=seqInit;i<=seqTerm;i++)
	{
		if(proSeq[i].simpleType==1 && proSeq[i].ins==' ' /*&& (proSeq[i].alt==' ' || proSeq[i].alt=='A')*/ )
		{
			if(atomType && (proSeq[i].detail==5 || proSeq[i].detail==12))
			{				
				continue;
			}
			fillAtom(i,proSeq,bColor);
//			totalNumber++;
		}
	}
//	printf("%d\n",totalNumber);
	for(i=0;i<pLength;i++)
	{   
		for(j=0;j<pWidth;j++)
		{
			for(k=0;k<pHeight;k++)
			{
				if(volumePixels[i][j][k].inOut)
				{
					volumePixels[i][j][k].isDone=true;
				}
			}
		}	
	}
}
//use isDone
void ProteinSurface::fillVoxelsWaals(int seqInit,int seqTerm,bool atomType,atom* proSeq,bool bColor)
{
	int i,j,k;
	if(volumePixels==NULL)
	{  
		volumePixels=new volumePixel**[pLength];
		for(i=0;i<pLength;i++)
			volumePixels[i]=new volumePixel*[pWidth];
		for(i=0;i<pLength;i++)
		{   
			for(j=0;j<pWidth;j++)
			{
				volumePixels[i][j]=new volumePixel[pHeight];
			}
		}
	}
	

	for(i=0;i<pLength;i++)
	{   
		for(j=0;j<pWidth;j++)
		{
			for(k=0;k<pHeight;k++)
			{
				volumePixels[i][j][k].isDone=false;
			}
		}	
	}
	for(i=seqInit;i<=seqTerm;i++)
	{
		if(proSeq[i].simpleType==1 && proSeq[i].ins==' ')
		{
			if(atomType && (proSeq[i].detail==5 ||proSeq[i].detail==12))
			{
				continue;
			}
			fillAtomWaals(i,proSeq,bColor);
		}
	}
}
void ProteinSurface::fillAtomWaals(int indx,atom* proSeq,bool bColor)
{
	int cx,cy,cz;
	int ox,oy,oz;
	Vector3d cp;

	cp = (proSeq[indx].coordinates + pTran) * scaleFactor;
	//Translating and scaling

	cx=int(cp(X)+0.5);
	cy=int(cp(Y)+0.5);
	cz=int(cp(Z)+0.5);
	//Rounding to the nearest integer

	int at=proSeq[indx].detail;
	int i,j,k;
	int ii,jj,kk;
	int mi,mj,mk;
	int si,sj,sk;
	int tIndex;
	int nIndex=0;
	for(i=0;i<widXz[at];i++)
	{
		for(j=0;j<widXz[at];j++)
		{
			if(deptY[at][nIndex]!=-1)
			{
				
				for( ii=-1;ii<2;ii++)
				{
					for( jj=-1;jj<2;jj++)
					{
						for( kk=-1;kk<2;kk++)
						{
							if(ii!=0 && jj!=0 && kk!=0)
							{		
								mi=ii*i;
								mk=kk*j;
								for(k=0;k<=deptY[at][nIndex];k++)
								{
									mj=k*jj;
									si=cx+mi;
									sj=cy+mj;
									sk=cz+mk;
									if(si<0 || sj<0 || sk<0)
									{
										continue;
									}
									if(!bColor)
									{
										volumePixels[si][sj][sk].isDone=true;
										continue;
									}
								else{
									if(volumePixels[si][sj][sk].isDone==false)
									{
										volumePixels[si][sj][sk].isDone=true;
										volumePixels[si][sj][sk].atomId=indx;
									}
									//with atomic info change above line
								//*
								    else if(volumePixels[si][sj][sk].isDone)
									{
										tIndex=volumePixels[si][sj][sk].atomId;
										cp = (proSeq[tIndex].coordinates + pTran) * scaleFactor;
										//Translating and scaling
										ox=int(cp(X)+0.5)-si;
										oy=int(cp(Y)+0.5)-sj;
										oz=int(cp(Z)+0.5)-sk;
										if(mi*mi+mj*mj+mk*mk<ox*ox+oy*oy+oz*oz)
											volumePixels[si][sj][sk].atomId=indx;
									}
							//	 */
									}//else
								}//k
								
							}//if
						}//kk	
					}//jj
				}//ii
				
				
			}//if
			nIndex++;
		}//j
	}//i
}

void ProteinSurface::fastDistanceMap(int type)
{
	int i,j,k;
	int positIn,positOut,eliminate;
	int certificate;
	totalSurfaceVox=0;
	totalInnerVox=0;
	voxel2 ***boundPoint;
	boundPoint=new voxel2 **[pLength];

	for(i=0;i<pLength;i++)
	{
		boundPoint[i]=new voxel2*[pWidth];
		for(j=0;j<pWidth;j++)
		{
			boundPoint[i][j]=new voxel2[pHeight];
			for(k=0;k<pHeight;k++)
			{
				volumePixels[i][j][k].isDone=false;
				if(volumePixels[i][j][k].inOut)
				{
					if(volumePixels[i][j][k].isBound)
					{
						totalSurfaceVox++;
						boundPoint[i][j][k].ix=i;
						boundPoint[i][j][k].iy=j;
						boundPoint[i][j][k].iz=k;
						volumePixels[i][j][k].distance=0;
						volumePixels[i][j][k].isDone=true;
					}
				    else
					{
						totalInnerVox++;
					}
				}
			}
		}
	}
	int allocIn=int(1.2*totalSurfaceVox);
	int allocOut=int(1.2*totalSurfaceVox);
	if(allocIn>totalInnerVox)
		allocIn=totalInnerVox;
	if(allocIn<totalSurfaceVox)
		allocIn=totalSurfaceVox;
	if(allocOut>totalInnerVox)
		allocOut=totalInnerVox;
	 inArray=new voxel2[allocIn];
	 outArray=new voxel2[allocOut];
	 positIn=0;positOut=0;
 
	 for(i=0;i<pLength;i++)
	 {
		 for(j=0;j<pWidth;j++)
		 {
			 for(k=0;k<pHeight;k++)
			 {	
				 if(volumePixels[i][j][k].isBound)
				 {
					 inArray[positIn].ix=i;
					 inArray[positIn].iy=j;
					 inArray[positIn].iz=k;
					 positIn++;
					 volumePixels[i][j][k].isBound=false;//as flag of outArray
				 }		 
			 }
		 }
	 }
	certificate=totalInnerVox;
/////////////////////////////////////////////////// 
if(type==0)//do part
{
	do {
		fastOneShell(&positIn, &allocOut, boundPoint, &positOut,&eliminate);
	//	printf("%d %d %d %d %d\n",positIn,allocOut,positOut,totalSurfaceVox,totalInnerVox);
		certificate-=eliminate;
	/*
		for(i=0;i<positOut;i++)
			{
			  inArray[i].ix=outArray[i].ix;
			  inArray[i].iy=outArray[i].iy;
			  inArray[i].iz=outArray[i].iz;
			}
			positIn=positOut;*/
		//new code only less dist
		positIn=0;
		for(i=0;i<positOut;i++)
		{
			volumePixels[outArray[i].ix][outArray[i].iy][outArray[i].iz].isBound=false;
			if(volumePixels[outArray[i].ix][outArray[i].iy][outArray[i].iz].distance<=1.02*cutRadius)
			{
				inArray[positIn].ix=outArray[i].ix;
				inArray[positIn].iy=outArray[i].iy;
				inArray[positIn].iz=outArray[i].iz;
				positIn++;
			}
			if(positIn>=allocIn)
			{
				allocIn*=2;
				if(allocIn>totalInnerVox) allocIn=totalInnerVox;
				inArray=(voxel2 *)realloc(inArray,allocIn*sizeof(voxel2));
			}
		}
	} 
	while(positIn!=0);
}
else if(type==1)//do all
{
	voxel2 *tpoint;
	do {
		
		fastOneShell( &positIn, &allocOut, boundPoint, &positOut,&eliminate);//inArray, outArray,
		certificate-=eliminate;
//	/*
//			for(i=0;i<positOut;i++)
//			{
//				volumePixels[outArray[i].ix][outArray[i].iy][outArray[i].iz].isBound=false;
//			  inArray[i].ix=outArray[i].ix;
//			  inArray[i].iy=outArray[i].iy;
//			  inArray[i].iz=outArray[i].iz;
//			}
			tpoint=inArray;
			inArray=outArray;
			outArray=tpoint;
			positIn=positOut;
			int alloctmp;
			alloctmp=allocIn;
			allocIn=allocOut;
			allocOut=alloctmp;
			for(i=0;i<positIn;i++)
				volumePixels[inArray[i].ix][inArray[i].iy][inArray[i].iz].isBound=false;
//			*/
		//new code only less dist
		/*
		positIn=0;
		for(i=0;i<positOut;i++)
		{
			volumePixels[outArray[i].ix][outArray[i].iy][outArray[i].iz].isBound=false;
			if(volumePixels[outArray[i].ix][outArray[i].iy][outArray[i].iz].distance<=1.0*cutRadius)
			{
				inArray[positIn].ix=outArray[i].ix;
				inArray[positIn].iy=outArray[i].iy;
				inArray[positIn].iz=outArray[i].iz;
				positIn++;
			}
		}
		*/
		
	} 
//	while(positIn!=0);
	while(positOut!=0);		
}
	//while(positOut!=0);			 
	if(certificate!=0) 
	{
	//	printf("wrong number\n");
	}

	 free(inArray);
	 free(outArray);

	 double cutsf=scaleFactor-0.5;
	 if(cutsf<0) cutsf=0;
//	 cutsf=100000000;
	 for(i=0;i<pLength;i++)
	 {
		 for(j=0;j<pWidth;j++)
		 {
			 for(k=0;k<pHeight;k++)
			 {
				 volumePixels[i][j][k].isBound=false;
				 //ses solid
				 if(volumePixels[i][j][k].inOut)
				 {
					 if(!volumePixels[i][j][k].isDone 
						 || (volumePixels[i][j][k].isDone && volumePixels[i][j][k].distance>=cutRadius-0.50/(0.1+cutsf))//0.33  0.75/scaleFactor
						 )
					 {
						 volumePixels[i][j][k].isBound=true;
						 //new add
						 if(volumePixels[i][j][k].isDone)
							volumePixels[i][j][k].atomId=volumePixels[boundPoint[i][j][k].ix][boundPoint[i][j][k].iy][boundPoint[i][j][k].iz].atomId;
					 }
				 }		 
			 }
		 }
	 }

	 
	 for(i=0;i<pLength;i++)
	 {
		 for(j=0;j<pWidth;j++)
		 {
			 delete[]boundPoint[i][j];
		 }
 		 delete[]boundPoint[i];
	 }
	 delete[] boundPoint;
}

void ProteinSurface::fastOneShell(int* inNum,int *allocOut,voxel2 ***boundPoint, int* outNum, int *elimi)
{
	int i, number,positOut;
	int tx,ty,tz;
	int dx,dy,dz;
	int eliminate=0;
	float squre;
	positOut=0;
	number=*inNum;
	if(number==0) return;
	//new code
	int j;
	voxel tnv;
//	Vector3i tnv;
	for(i=0;i<number;i++)
	{
		if(positOut>=(*allocOut)-6)
		{
			(*allocOut)=int(1.2*(*allocOut));
			if(*allocOut>totalInnerVox) *allocOut=totalInnerVox;
			outArray=(voxel2 *)realloc(outArray,(*allocOut)*sizeof(voxel2));
		}
		tx=inArray[i].ix;
        ty=inArray[i].iy;
		tz=inArray[i].iz;

		//tnv = inArray[i].ixyz;

		for(j=0;j<6;j++)
		{
			tnv.ix=tx+nb[j][0];
			tnv.iy=ty+nb[j][1];
			tnv.iz=tz+nb[j][2];
			if( tnv.ix<pLength && tnv.ix>-1 && 
				tnv.iy<pWidth && tnv.iy>-1 && 
				tnv.iz<pHeight && tnv.iz>-1 && 
				volumePixels[tnv.ix][tnv.iy][tnv.iz].inOut && 
				!volumePixels[tnv.ix][tnv.iy][tnv.iz].isDone)
			{
				boundPoint[tnv.ix][tnv.iy][tz+nb[j][2]].ix=boundPoint[tx][ty][tz].ix;
				boundPoint[tnv.ix][tnv.iy][tz+nb[j][2]].iy=boundPoint[tx][ty][tz].iy;
				boundPoint[tnv.ix][tnv.iy][tz+nb[j][2]].iz=boundPoint[tx][ty][tz].iz;
				dx=tnv.ix-boundPoint[tx][ty][tz].ix;
				dy=tnv.iy-boundPoint[tx][ty][tz].iy;
				dz=tnv.iz-boundPoint[tx][ty][tz].iz;
				squre=float(dx*dx+dy*dy+dz*dz);
				volumePixels[tnv.ix][tnv.iy][tnv.iz].distance=float(sqrt(squre));
				volumePixels[tnv.ix][tnv.iy][tnv.iz].isDone=true;
				volumePixels[tnv.ix][tnv.iy][tnv.iz].isBound=true;
				outArray[positOut].ix=tnv.ix;
				outArray[positOut].iy=tnv.iy;
				outArray[positOut].iz=tnv.iz;
				positOut++;eliminate++;
			}
			else if( tnv.ix<pLength && tnv.ix>-1 && 
				tnv.iy<pWidth && tnv.iy>-1 && 
				tnv.iz<pHeight && tnv.iz>-1 && 
				volumePixels[tnv.ix][tnv.iy][tnv.iz].inOut && 
				volumePixels[tnv.ix][tnv.iy][tnv.iz].isDone)
			{
				
				dx=tnv.ix-boundPoint[tx][ty][tz].ix;
				dy=tnv.iy-boundPoint[tx][ty][tz].iy;
				dz=tnv.iz-boundPoint[tx][ty][tz].iz;
				squre=float(dx*dx+dy*dy+dz*dz);
				squre=float(sqrt(squre));
				if(squre<volumePixels[tnv.ix][tnv.iy][tnv.iz].distance)
				{
					boundPoint[tnv.ix][tnv.iy][tnv.iz].ix=boundPoint[tx][ty][tz].ix;
					boundPoint[tnv.ix][tnv.iy][tnv.iz].iy=boundPoint[tx][ty][tz].iy;
					boundPoint[tnv.ix][tnv.iy][tnv.iz].iz=boundPoint[tx][ty][tz].iz;
					volumePixels[tnv.ix][tnv.iy][tnv.iz].distance=float(squre);
					if(!volumePixels[tnv.ix][tnv.iy][tnv.iz].isBound)
					{
						volumePixels[tnv.ix][tnv.iy][tnv.iz].isBound=true;
						outArray[positOut].ix=tnv.ix;
						outArray[positOut].iy=tnv.iy;
						outArray[positOut].iz=tnv.iz;
						positOut++;
					}
				}
		
			}
		}
	}
	for(i=0;i<number;i++)
	{
		if(positOut>=(*allocOut)-12)
		{
			(*allocOut)=int(1.2*(*allocOut));
			if(*allocOut>totalInnerVox) *allocOut=totalInnerVox;
			outArray=(voxel2 *)realloc(outArray,(*allocOut)*sizeof(voxel2));
		}
		tx=inArray[i].ix;
        ty=inArray[i].iy;
		tz=inArray[i].iz;
		//txyz = inArray[i].ixyz;
		for(j=6;j<18;j++)
		{
			tnv.ix=tx+nb[j][0];
			tnv.iy=ty+nb[j][1];
			tnv.iz=tz+nb[j][2];
			if( tnv.ix<pLength && tnv.ix>-1 && 
				tnv.iy<pWidth && tnv.iy>-1 && 
				tnv.iz<pHeight && tnv.iz>-1 && 
				volumePixels[tnv.ix][tnv.iy][tnv.iz].inOut && 
				!volumePixels[tnv.ix][tnv.iy][tnv.iz].isDone)
			{
				boundPoint[tnv.ix][tnv.iy][tz+nb[j][2]].ix=boundPoint[tx][ty][tz].ix;
				boundPoint[tnv.ix][tnv.iy][tz+nb[j][2]].iy=boundPoint[tx][ty][tz].iy;
				boundPoint[tnv.ix][tnv.iy][tz+nb[j][2]].iz=boundPoint[tx][ty][tz].iz;
				dx=tnv.ix-boundPoint[tx][ty][tz].ix;
				dy=tnv.iy-boundPoint[tx][ty][tz].iy;
				dz=tnv.iz-boundPoint[tx][ty][tz].iz;
				squre=float(dx*dx+dy*dy+dz*dz);
				volumePixels[tnv.ix][tnv.iy][tnv.iz].distance=float(sqrt(squre));
				volumePixels[tnv.ix][tnv.iy][tnv.iz].isDone=true;
				volumePixels[tnv.ix][tnv.iy][tnv.iz].isBound=true;
				outArray[positOut].ix=tnv.ix;
				outArray[positOut].iy=tnv.iy;
				outArray[positOut].iz=tnv.iz;
				positOut++;eliminate++;
			}
			else if( tnv.ix<pLength && tnv.ix>-1 && 
				tnv.iy<pWidth && tnv.iy>-1 && 
				tnv.iz<pHeight && tnv.iz>-1 && 
				volumePixels[tnv.ix][tnv.iy][tnv.iz].inOut && 
				volumePixels[tnv.ix][tnv.iy][tnv.iz].isDone)
			{
				dx=tnv.ix-boundPoint[tx][ty][tz].ix;
				dy=tnv.iy-boundPoint[tx][ty][tz].iy;
				dz=tnv.iz-boundPoint[tx][ty][tz].iz;
				squre=float(dx*dx+dy*dy+dz*dz);
				squre=float(sqrt(squre));
				if(squre<volumePixels[tnv.ix][tnv.iy][tnv.iz].distance)
				{
					boundPoint[tnv.ix][tnv.iy][tnv.iz].ix=boundPoint[tx][ty][tz].ix;
					boundPoint[tnv.ix][tnv.iy][tnv.iz].iy=boundPoint[tx][ty][tz].iy;
					boundPoint[tnv.ix][tnv.iy][tnv.iz].iz=boundPoint[tx][ty][tz].iz;
					volumePixels[tnv.ix][tnv.iy][tnv.iz].distance=float(squre);
					if(!volumePixels[tnv.ix][tnv.iy][tnv.iz].isBound)
					{
						volumePixels[tnv.ix][tnv.iy][tnv.iz].isBound=true;
						outArray[positOut].ix=tnv.ix;
						outArray[positOut].iy=tnv.iy;
						outArray[positOut].iz=tnv.iz;
						positOut++;
					}
				}
				
			}
		}
	}
	for(i=0;i<number;i++)
	{
		if(positOut>=(*allocOut)-9)
		{
			(*allocOut)=int(1.2*(*allocOut));
			if(*allocOut>totalInnerVox) *allocOut=totalInnerVox;
			outArray=(voxel2 *)realloc(outArray,(*allocOut)*sizeof(voxel2));
		}
		tx=inArray[i].ix;
        ty=inArray[i].iy;
		tz=inArray[i].iz;
		for(j=18;j<26;j++)
		{
			tnv.ix=tx+nb[j][0];
			tnv.iy=ty+nb[j][1];
			tnv.iz=tz+nb[j][2];
			if( tnv.ix<pLength && tnv.ix>-1 && 
				tnv.iy<pWidth && tnv.iy>-1 && 
				tnv.iz<pHeight && tnv.iz>-1 && 
				volumePixels[tnv.ix][tnv.iy][tnv.iz].inOut && 
				!volumePixels[tnv.ix][tnv.iy][tnv.iz].isDone)
			{
				boundPoint[tnv.ix][tnv.iy][tz+nb[j][2]].ix=boundPoint[tx][ty][tz].ix;
				boundPoint[tnv.ix][tnv.iy][tz+nb[j][2]].iy=boundPoint[tx][ty][tz].iy;
				boundPoint[tnv.ix][tnv.iy][tz+nb[j][2]].iz=boundPoint[tx][ty][tz].iz;
				dx=tnv.ix-boundPoint[tx][ty][tz].ix;
				dy=tnv.iy-boundPoint[tx][ty][tz].iy;
				dz=tnv.iz-boundPoint[tx][ty][tz].iz;
				squre=float(dx*dx+dy*dy+dz*dz);
				volumePixels[tnv.ix][tnv.iy][tnv.iz].distance=float(sqrt(squre));
				volumePixels[tnv.ix][tnv.iy][tnv.iz].isDone=true;
				volumePixels[tnv.ix][tnv.iy][tnv.iz].isBound=true;
				outArray[positOut].ix=tnv.ix;
				outArray[positOut].iy=tnv.iy;
				outArray[positOut].iz=tnv.iz;
				positOut++;eliminate++;
			}
			else if( tnv.ix<pLength && tnv.ix>-1 && 
				tnv.iy<pWidth && tnv.iy>-1 && 
				tnv.iz<pHeight && tnv.iz>-1 && 
				volumePixels[tnv.ix][tnv.iy][tnv.iz].inOut && 
				volumePixels[tnv.ix][tnv.iy][tnv.iz].isDone)
			{
				
				dx=tnv.ix-boundPoint[tx][ty][tz].ix;
				dy=tnv.iy-boundPoint[tx][ty][tz].iy;
				dz=tnv.iz-boundPoint[tx][ty][tz].iz;
				squre=float(dx*dx+dy*dy+dz*dz);
				squre=float(sqrt(squre));
				if(squre<volumePixels[tnv.ix][tnv.iy][tnv.iz].distance)
				{
					boundPoint[tnv.ix][tnv.iy][tnv.iz].ix=boundPoint[tx][ty][tz].ix;
					boundPoint[tnv.ix][tnv.iy][tnv.iz].iy=boundPoint[tx][ty][tz].iy;
					boundPoint[tnv.ix][tnv.iy][tnv.iz].iz=boundPoint[tx][ty][tz].iz;
					volumePixels[tnv.ix][tnv.iy][tnv.iz].distance=float(squre);
					if(!volumePixels[tnv.ix][tnv.iy][tnv.iz].isBound)
					{
						volumePixels[tnv.ix][tnv.iy][tnv.iz].isBound=true;
						outArray[positOut].ix=tnv.ix;
						outArray[positOut].iy=tnv.iy;
						outArray[positOut].iz=tnv.iz;
						positOut++;
					}
				}
				
			}
		}
	}
	
	*outNum=positOut;
	*elimi=eliminate;
	 
}
//include inner
void ProteinSurface::calcAreaVolume()
{
	int i,j,k;
	double totalVolume=0.10*scaleFactor*vertexNumber;
	for(i=0;i<pLength;i++)
	{
		for(j=0;j<pWidth;j++)
		{
			for(k=0;k<pHeight;k++)
			{
				if(volumePixels[i][j][k].isDone)
					totalVolume+=1;			
			}
		}
	}
	double totalArea=0;
	for(i=0;i<faceNumber;i++)
	{
		totalArea+=faces[i].area;
	}

	surfaceArea=totalArea/scaleFactor/scaleFactor;
	surfaceVolume=(totalVolume)/scaleFactor/scaleFactor/scaleFactor;

}
//cavity number
void ProteinSurface::cavityNumbers()
{
	int i,j,k;
	int l,ii,jj,kk;
	voxel *inArray,*outArray;
	voxel *tArray,tnv;
	int inNum=0;
	int outNum;
	int totalNumber=0;
	if(eachCavityVolume!=NULL)
	{
		free(eachCavityVolume);
	}
	eachCavityVolume=new double[1500];
	for(i=0;i<1500;i++)
	{
		eachCavityVolume[i]=0;
	}
	for(i=0;i<pLength;i++)
	{
		for(j=0;j<pWidth;j++)
		{
			for(k=0;k<pHeight;k++)
			{
				volumePixels[i][j][k].isBound=false;
				if(!volumePixels[i][j][k].isDone)
				{
					totalNumber++;
				}
			}
		}
	}
	inArray=new voxel[totalNumber];
	outArray=new voxel[totalNumber];
	numberOfCavities=0;
	for(l=0;l<pLength;l++)
	{
		for(j=0;j<pWidth;j++)
		{
			for(k=0;k<pHeight;k++)
			{
				if(!volumePixels[l][j][k].isBound && !volumePixels[l][j][k].isDone && volumePixels[l][j][k].inOut)
				{
					volumePixels[l][j][k].distance=float(numberOfCavities);
					volumePixels[l][j][k].isBound=true;
					inNum=0;
					inArray[inNum].ix=l;
					inArray[inNum].iy=j;
					inArray[inNum].iz=k;
					inNum++;
					eachCavityVolume[numberOfCavities]=eachCavityVolume[numberOfCavities]+1;
					while(inNum!=0)
					{
						outNum=0;
						for(i=0;i<inNum;i++)
						{
							for(ii=-1;ii<2;ii++)
							{
								for(jj=-1;jj<2;jj++)
								{
									for(kk=-1;kk<2;kk++)
									{
										tnv.ix=inArray[i].ix+ii;
										tnv.iy=inArray[i].iy+jj;
										tnv.iz=inArray[i].iz+kk;
										if( abs(ii)+abs(jj)+abs(kk)<2 &&
											tnv.ix>-1 && tnv.ix<pLength 
											&& tnv.iy>-1 && tnv.iy<pWidth
											&& tnv.iz>-1 && tnv.iz<pHeight
											&& volumePixels[tnv.ix][tnv.iy][tnv.iz].inOut
											&& !volumePixels[tnv.ix][tnv.iy][tnv.iz].isDone
											&& !volumePixels[tnv.ix][tnv.iy][tnv.iz].isBound)
										{
											volumePixels[tnv.ix][tnv.iy][tnv.iz].isBound=true;
											volumePixels[tnv.ix][tnv.iy][tnv.iz].distance=float(numberOfCavities);
											outArray[outNum].ix=tnv.ix;
											outArray[outNum].iy=tnv.iy;
											outArray[outNum].iz=tnv.iz;
											outNum++;
											eachCavityVolume[numberOfCavities]=eachCavityVolume[numberOfCavities]+1;
										}
									}//kk
								}//jj
							}//ii
						}//i
						tArray=outArray;
						outArray=inArray;
						inArray=tArray;
						inNum=outNum;
					}//while
					numberOfCavities++;
				}//if
			}//k
		}//j
	}//l
	delete[]inArray;
	delete[]outArray;
	for(i=0;i<numberOfCavities;i++)
	{
		eachCavityVolume[i]=eachCavityVolume[i]/scaleFactor/scaleFactor/scaleFactor;
	}

}
//cavity area and volume surface inner or outer
void ProteinSurface::cavitiesAreaVolume()
{
	int i,j,k;
	double totalArea=0;
	int totalVolume=0;
	for(i=0;i<pLength;i++)
	{
		for(j=0;j<pWidth;j++)
		{
			for(k=0;k<pHeight;k++)
			{
				if(!volumePixels[i][j][k].isDone && volumePixels[i][j][k].inOut)
					totalVolume++;
			}
		}
	}
	int tx[3],ty[3],tz[3];
	//Vector3i txyz[3];
	int ii,jj,kk, ll;
	bool flagInner;
	for(i=0;i<faceNumber;i++)
	{
		faces[i].inOut=false;
		verts[faces[i].abc(A)].inOut=true;
		verts[faces[i].abc(B)].inOut=true;
		verts[faces[i].abc(C)].inOut=true;
		flagInner=false;
		tx[0]=int(verts[faces[i].abc(A)].xyz(X));
		ty[0]=int(verts[faces[i].abc(A)].xyz(Y));
		tz[0]=int(verts[faces[i].abc(A)].xyz(Z));
		tx[1]=int(verts[faces[i].abc(B)].xyz(X));
		ty[1]=int(verts[faces[i].abc(B)].xyz(Y));
		tz[1]=int(verts[faces[i].abc(B)].xyz(Z));
		tx[2]=int(verts[faces[i].abc(C)].xyz(X));
		ty[2]=int(verts[faces[i].abc(C)].xyz(Y));
		tz[2]=int(verts[faces[i].abc(C)].xyz(Z));
		//txyz[0] = Vector3i(verts[faces[i].abc(A)].xyz);
		//txyz[1] = Vector3i(verts[faces[i].abc(B)].xyz);
		//txyz[2] = Vector3i(verts[faces[i].abc(C)].xyz);




		for(ii=-1; ii<2; ii++){
			for(jj = -1; jj < 2; jj++){
				for(kk = -1; kk < 2; kk++){
					for(ll = 0; ll < 3; ll ++){
						if(!(ii==0 && jj==0 && kk==0) && tx[ll]+ii>-1 &&  tx[ll]+ii<pLength
							&& ty[ll]+jj>-1 &&  ty[ll]+jj<pWidth
							&& tz[ll]+kk>-1 &&  tz[ll]+kk<pHeight
							&& volumePixels[tx[ll]+ii][ty[ll]+jj][tz[ll]+kk].inOut
							&& !volumePixels[tx[ll]+ii][ty[ll]+jj][tz[ll]+kk].isDone)
						{
							flagInner=true;
						}						
					}
				}
			}
		}

		if(flagInner)
		{
			faces[i].inOut=true;
			verts[faces[i].abc(A)].inOut=false;
			verts[faces[i].abc(B)].inOut=false;
			verts[faces[i].abc(C)].inOut=false;
			totalArea+=faces[i].area;
		}
	}
	cavityArea=totalArea/scaleFactor/scaleFactor;
	cavityVolume=totalVolume/scaleFactor/scaleFactor/scaleFactor;
}
void ProteinSurface::computeNorm()
{
	int i;
	double pNorm;
	for(i=0;i<vertexNumber;i++)
	{
	/*	verts[i].pn(X)=0;
		verts[i].pn(Y)=0;
		verts[i].pn(Z)=0;
	*/	verts[i].pn *= 0;
	}
	Vector3d p1,p2,p3;
	Vector3d p12,p13;
	Vector3d pn;

	for(i=0;i<faceNumber;i++)
	//For each face, find the normal vector
	{
/*		p1(X)=verts[faces[i].abc(A)].xyz(X);
		p1(Y)=verts[faces[i].abc(A)].xyz(Y);
		p1(Z)=verts[faces[i].abc(A)].xyz(Z);

*/		p1 = verts[faces[i].abc(A)].xyz;
		//Point one (one vertex of a triangle)

/*		p2(X)=verts[faces[i].abc(B)].xyz(X);
		p2(Y)=verts[faces[i].abc(B)].xyz(Y);
		p2(Z)=verts[faces[i].abc(B)].xyz(Z);

*/		p2 = verts[faces[i].abc(B)].xyz;
		//Point two

/*		p3(X)=verts[faces[i].abc(C)].xyz(X);
		p3(Y)=verts[faces[i].abc(C)].xyz(Y);
		p3(Z)=verts[faces[i].abc(C)].xyz(Z);

*/		p3 = verts[faces[i].abc(C)].xyz;
		//Point three

/*		p12(X)=p2(X)-p1(X);
		p12(Y)=p2(Y)-p1(Y);
		p12(Z)=p2(Z)-p1(Z);

*/		p12 = p2 - p1;
		//Vector connecting point two to point one

/*		p13(X)=p3(X)-p1(X);
		p13(Y)=p3(Y)-p1(Y);
		p13(Z)=p3(Z)-p1(Z);

*/		p13 = p3 - p1;
		//Vector connecting ppoint three to point one

/*		pn(X)=p12(Y)*p13(Z)-p12(Z)*p13(Y);
		pn(Y)=p12(Z)*p13(X)-p12(X)*p13(Z);
		pn(Z)=p12(X)*p13(Y)-p12(Y)*p13(X);
*/
		pn = p12.cross(p13);

		//Cross product of these two vectors should be the vector normal to the plane

		faces[i].area=0.5*sqrt(pn(X)*pn(X)+pn(Y)*pn(Y)+pn(Z)*pn(Z));

		//Set the area equal to the magnitude of the normal vector

/*		faces[i].pn(X)=0.5*pn(X)/faces[i].area;
		faces[i].pn(Y)=0.5*pn(Y)/faces[i].area;
		faces[i].pn(Z)=0.5*pn(Z)/faces[i].area;

*/		faces[i].pn = 0.5*pn/faces[i].area;

		//*
		//without area
/*		verts[faces[i].abc(A)].pn(X)+=faces[i].pn(X);
		verts[faces[i].abc(A)].pn(Y)+=faces[i].pn(Y);
		verts[faces[i].abc(A)].pn(Z)+=faces[i].pn(Z);
*/		verts[faces[i].abc(A)].pn += faces[i].pn;
/*		verts[faces[i].abc(B)].pn(X)+=faces[i].pn(X);
		verts[faces[i].abc(B)].pn(Y)+=faces[i].pn(Y);
		verts[faces[i].abc(B)].pn(Z)+=faces[i].pn(Z);
*/		verts[faces[i].abc(B)].pn += faces[i].pn;
/*		verts[faces[i].abc(C)].pn(X)+=faces[i].pn(X);
		verts[faces[i].abc(C)].pn(Y)+=faces[i].pn(Y);
		verts[faces[i].abc(C)].pn(Z)+=faces[i].pn(Z);
*/		verts[faces[i].abc(C)].pn += faces[i].pn;


		//*/
		/*
		//with area
		verts[faces[i].abc(A)].pn(X)+=pn(X);
		verts[faces[i].abc(A)].pn(Y)+=pn(Y);
		verts[faces[i].abc(A)].pn(Z)+=pn(Z);
		verts[faces[i].abc(B)].pn(X)+=pn(X);
		verts[faces[i].abc(B)].pn(Y)+=pn(Y);
		verts[faces[i].abc(B)].pn(Z)+=pn(Z);
		verts[faces[i].abc(C)].pn(X)+=pn(X);
		verts[faces[i].abc(C)].pn(Y)+=pn(Y);
		verts[faces[i].abc(C)].pn(Z)+=pn(Z);
		*/
	}
	for(i=0;i<vertexNumber;i++)
	//For each vertex find the unit normal vector
	{
/*		pn(X)=verts[i].pn(X);
		pn(Y)=verts[i].pn(Y);
		pn(Z)=verts[i].pn(Z);
*/		pn = verts[i].pn;
		pNorm=sqrt(pn(X)*pn(X)+pn(Y)*pn(Y)+pn(Z)*pn(Z));
//		pNorm = sqrt(sum(pn.dot(pn)));
		if(pNorm==0.0)
		{
		/*	pn(X)=0.0;
			pn(Y)=0.0;
			pn(Z)=0.0;
		*/	pn *= 0.0;
			//This asSignment shouldn't be necessary
/*			verts[i].pn(X)=pn(X);
			verts[i].pn(Y)=pn(Y);
			verts[i].pn(Z)=pn(Z);
*/			verts[i].pn = pn;
			continue;
		}
/*		pn(X)/=pNorm;
		pn(Y)/=pNorm;
		pn(Z)/=pNorm;
*/		pn /= pNorm;
/*		verts[i].pn(X)=pn(X);
		verts[i].pn(Y)=pn(Y);
		verts[i].pn(Z)=pn(Z);
*/		verts[i].pn = pn;
	}
}
void ProteinSurface::marchingCubeInit(int surfaceType)
{
	int i,j,k;

	for(i = 0; i < pLength * pWidth * pHeight; i++){
		//This is a way of representing three nested for loops using one index
		//volumePixels[i/(pWidth*pHeight)][(i/pHeight)%pWidth][i%pHeight] is equivalent to volumePixels[i][j][k] if we'd used nested loops

		
		if(surfaceType == SES){
			if(volumePixels[i/(pWidth*pHeight)][(i/pHeight)%pWidth][i%pHeight].isBound){
				volumePixels[i/(pWidth*pHeight)][(i/pHeight)%pWidth][i%pHeight].isDone = true;
			}
			else{
				volumePixels[i/(pWidth*pHeight)][(i/pHeight)%pWidth][i%pHeight].isDone = false;
			}				
		}

		else if(surfaceType == MS){
			if(volumePixels[i/(pWidth*pHeight)][(i/pHeight)%pWidth][i%pHeight].isBound){
				if(volumePixels[i/(pWidth*pHeight)][(i/pHeight)%pWidth][i%pHeight].isDone){
					volumePixels[i/(pWidth*pHeight)][(i/pHeight)%pWidth][i%pHeight].isBound = false;
				}
				else{
					volumePixels[i/(pWidth*pHeight)][(i/pHeight)%pWidth][i%pHeight].isDone = true;
				}
			}
		}

		volumePixels[i/(pWidth*pHeight)][(i/pHeight)%pWidth][i%pHeight].isBound = false;

		//We do this in every case	
	}

	
	
}
//half accuracy
void ProteinSurface::marchingCubeOrigin(int surfaceType)
{
	int i,j,k;
	marchingCubeInit(surfaceType);
	int ***vertSeq;
	vertSeq=new int**[pLength];
	for(i=0;i<pLength;i++)
	{
		vertSeq[i]=new int*[pWidth];
		for(j=0;j<pWidth;j++)
		{
			vertSeq[i][j]=new int[pHeight];
			for(k=0;k<pHeight;k++)
			{
				vertSeq[i][j][k]=-1;
			}
		}
	}

	if(faces!=NULL)
	{
		free(faces);
	}
	if(verts!=NULL)
	{
		free(verts);
	}
	int allocFace=20;
	int allocVert=12;
	faceNumber=0;
	vertexNumber=0;
	verts=new vertInfo[allocVert];
	faces=new faceInfo[allocFace];
	///////////////////////////////////////////
	int ii,jj;
	int tl[3];
	voxel tp[3][2],tv[3];
	int totIndex=0;
	for(i=0;i<pLength-2;i+=2)
	{
		for(j=0;j<pWidth-2;j+=2)
		{
			for(k=0;k<pHeight-2;k+=2)
			{
				if(vertexNumber+8>allocVert)
				{
					allocVert*=2;
					verts=(vertInfo *)realloc(verts,allocVert*sizeof(vertInfo));
				}
				if(faceNumber+5>allocFace)
				{
					allocFace*=2;
					faces=(faceInfo *)realloc(faces,allocFace*sizeof(faceInfo));
				}
				totIndex=0;
				for(ii=0;ii<8;ii++)
				{
					if(!volumePixels[i+2*a2fVertexOffset[ii][0]][j+2*a2fVertexOffset[ii][1]][k+2*a2fVertexOffset[ii][2]].isDone)
						totIndex |= 1<<ii;	
				}
				for(ii = 0; ii < 5; ii++)
				{
					if(a2iTriangleConnectionTable[totIndex][3*ii] < 0)
                        break;
					for(jj=0;jj<3;jj++)
					{
						tl[jj]=a2iTriangleConnectionTable[totIndex][3*ii+jj];
						tp[jj][0].ix=i+2*a2fVertexOffset[a2iEdgeConnection[tl[jj]][0]][0];
						tp[jj][0].iy=j+2*a2fVertexOffset[a2iEdgeConnection[tl[jj]][0]][1];
						tp[jj][0].iz=k+2*a2fVertexOffset[a2iEdgeConnection[tl[jj]][0]][2];
						tp[jj][1].ix=i+2*a2fVertexOffset[a2iEdgeConnection[tl[jj]][1]][0];
						tp[jj][1].iy=j+2*a2fVertexOffset[a2iEdgeConnection[tl[jj]][1]][1];
						tp[jj][1].iz=k+2*a2fVertexOffset[a2iEdgeConnection[tl[jj]][1]][2];
						tv[jj].ix=int(tp[jj][0].ix+tp[jj][1].ix)/2;
						tv[jj].iy=int(tp[jj][0].iy+tp[jj][1].iy)/2;
						tv[jj].iz=int(tp[jj][0].iz+tp[jj][1].iz)/2;
						if(!volumePixels[tv[jj].ix][tv[jj].iy][tv[jj].iz].isDone)
						{
							if(volumePixels[tp[jj][0].ix][tp[jj][0].iy][tp[jj][0].iz].isDone)
								volumePixels[tv[jj].ix][tv[jj].iy][tv[jj].iz]=volumePixels[tp[jj][0].ix][tp[jj][0].iy][tp[jj][0].iz];
							else if(volumePixels[tp[jj][1].ix][tp[jj][1].iy][tp[jj][1].iz].isDone)
								volumePixels[tv[jj].ix][tv[jj].iy][tv[jj].iz]=volumePixels[tp[jj][1].ix][tp[jj][1].iy][tp[jj][1].iz];
							

						}
						if(vertSeq[tv[jj].ix][tv[jj].iy][tv[jj].iz]==-1)
						{
							vertSeq[tv[jj].ix][tv[jj].iy][tv[jj].iz]=vertexNumber;
							verts[vertexNumber].xyz(X)=tv[jj].ix;
							verts[vertexNumber].xyz(Y)=tv[jj].iy;
							verts[vertexNumber].xyz(Z)=tv[jj].iz;
							vertexNumber++;
						}
					}	
				faces[faceNumber].abc(A)=vertSeq[tv[0].ix][tv[0].iy][tv[0].iz];
				faces[faceNumber].abc(B)=vertSeq[tv[1].ix][tv[1].iy][tv[1].iz];
				faces[faceNumber++].abc(C)=vertSeq[tv[2].ix][tv[2].iy][tv[2].iz];
				}
			}
		}
	}

	verts=(vertInfo *)realloc(verts,vertexNumber*sizeof(vertInfo));
	faces=(faceInfo *)realloc(faces,faceNumber*sizeof(faceInfo));
	
	for(i=0;i<vertexNumber;i++)
	{
		verts[i].atomId=volumePixels[int(verts[i].xyz(X))][int(verts[i].xyz(Y))][int(verts[i].xyz(Z))].atomId;
		verts[i].isCont=false;
		if(volumePixels[int(verts[i].xyz(X))][int(verts[i].xyz(Y))][int(verts[i].xyz(Z))].isBound)
			verts[i].isCont=true;
	}
	for(i=0;i<pLength;i++)
	{
		for(j=0;j<pWidth;j++)
		{
			delete[]vertSeq[i][j];
		}
		delete[]vertSeq[i];
	}
	
	delete[]vertSeq;

}
// middle points
void ProteinSurface::marchingCubeOrigin2(int surfaceType)
{
	int i,j,k;
	voxel ***volumePixelsInd;
	volumePixelsInd=new voxel**[pLength];

	for(i=0;i<pLength;i++)
	{
		volumePixelsInd[i]=new voxel*[pWidth];
		for(j=0;j<pWidth;j++)
		{
			volumePixelsInd[i][j]=new voxel[pHeight];
			for(k=0;k<pHeight;k++)
			{
				volumePixelsInd[i][j][k].ix=-1;
				volumePixelsInd[i][j][k].iy=-1;
				volumePixelsInd[i][j][k].iz=-1;
			}
		}
	}
	marchingCubeInit(surfaceType);
	if(faces!=NULL)
	{
		free(faces);
	}
	if(verts!=NULL)
	{
		free(verts);
	}
	int allocFace=20;
	int allocVert=12;
	faceNumber=0;
	vertexNumber=0;
	verts=new vertInfo[allocVert];
	faces=new faceInfo[allocFace];



	///////////////////////////////////////////
	int ii,jj;
	int tl[3];
	voxel tp[3][2];
	Vector3d tv[3];
	int indv[3];
	int totIndex=0;
	int curpt;
	for(i=0;i<pLength-1;i++)
	{
		for(j=0;j<pWidth-1;j++)
		{
			for(k=0;k<pHeight-1;k++)
			{
				if(vertexNumber+8>allocVert)
				{
					allocVert*=2;
					verts=(vertInfo *)realloc(verts,allocVert*sizeof(vertInfo));
				}
				if(faceNumber+5>allocFace)
				{
					allocFace*=2;
					faces=(faceInfo *)realloc(faces,allocFace*sizeof(faceInfo));
				}
				totIndex=0;
				for(ii=0;ii<8;ii++)
				{
					if(!volumePixels[i+a2fVertexOffset[ii][0]][j+a2fVertexOffset[ii][1]][k+a2fVertexOffset[ii][2]].isDone)
						totIndex |= 1<<ii;	
				}
				for(ii = 0; ii < 5; ii++)
				{
					if(a2iTriangleConnectionTable[totIndex][3*ii] < 0)
                        break;
					for(jj=0;jj<3;jj++)
					{
						tl[jj]=a2iTriangleConnectionTable[totIndex][3*ii+jj];
						tp[jj][0].ix=i+a2fVertexOffset[a2iEdgeConnection[tl[jj]][0]][0];
						tp[jj][0].iy=j+a2fVertexOffset[a2iEdgeConnection[tl[jj]][0]][1];
						tp[jj][0].iz=k+a2fVertexOffset[a2iEdgeConnection[tl[jj]][0]][2];
						tp[jj][1].ix=i+a2fVertexOffset[a2iEdgeConnection[tl[jj]][1]][0];
						tp[jj][1].iy=j+a2fVertexOffset[a2iEdgeConnection[tl[jj]][1]][1];
						tp[jj][1].iz=k+a2fVertexOffset[a2iEdgeConnection[tl[jj]][1]][2];
						tv[jj](X)=(tp[jj][0].ix+tp[jj][1].ix)/2.0;
						tv[jj](Y)=(tp[jj][0].iy+tp[jj][1].iy)/2.0;
						tv[jj](Z)=(tp[jj][0].iz+tp[jj][1].iz)/2.0;
												
						if(tp[jj][0].ix!=tp[jj][1].ix)
						{
							curpt=int(tv[jj](X));
							if(volumePixelsInd[curpt][tp[jj][1].iy][tp[jj][1].iz].ix==-1)
							{
								indv[jj]=vertexNumber;
								volumePixelsInd[curpt][tp[jj][1].iy][tp[jj][1].iz].ix=vertexNumber;
								verts[vertexNumber].xyz(X)=tv[jj](X);
								verts[vertexNumber].xyz(Y)=tv[jj](Y);
								verts[vertexNumber].xyz(Z)=tv[jj](Z);
								if(volumePixels[tp[jj][0].ix][tp[jj][0].iy][tp[jj][0].iz].isDone)
								{
									verts[vertexNumber].atomId=volumePixels[tp[jj][0].ix][tp[jj][0].iy][tp[jj][0].iz].atomId;
									verts[vertexNumber].isCont=false;
									if(volumePixels[tp[jj][0].ix][tp[jj][0].iy][tp[jj][0].iz].isBound)
										verts[vertexNumber].isCont=true;
								}
								else if(volumePixels[tp[jj][1].ix][tp[jj][1].iy][tp[jj][1].iz].isDone)
								{
									verts[vertexNumber].atomId=volumePixels[tp[jj][1].ix][tp[jj][1].iy][tp[jj][1].iz].atomId;
									verts[vertexNumber].isCont=false;
									if(volumePixels[tp[jj][1].ix][tp[jj][1].iy][tp[jj][1].iz].isBound)
										verts[vertexNumber].isCont=true;
								}
								
								vertexNumber++;
							}
							else
							{
								indv[jj]=volumePixelsInd[curpt][tp[jj][1].iy][tp[jj][1].iz].ix;
							}
						}	
						else if(tp[jj][0].iy!=tp[jj][1].iy)
						{
							curpt=int(tv[jj](Y));
							if(volumePixelsInd[tp[jj][1].ix][curpt][tp[jj][1].iz].iy==-1)
							{
								indv[jj]=vertexNumber;
								volumePixelsInd[tp[jj][1].ix][curpt][tp[jj][1].iz].iy=vertexNumber;
								verts[vertexNumber].xyz(X)=tv[jj](X);
								verts[vertexNumber].xyz(Y)=tv[jj](Y);
								verts[vertexNumber].xyz(Z)=tv[jj](Z);
								if(volumePixels[tp[jj][0].ix][tp[jj][0].iy][tp[jj][0].iz].isDone)
								{
									verts[vertexNumber].atomId=volumePixels[tp[jj][0].ix][tp[jj][0].iy][tp[jj][0].iz].atomId;
									verts[vertexNumber].isCont=false;
									if(volumePixels[tp[jj][0].ix][tp[jj][0].iy][tp[jj][0].iz].isBound)
										verts[vertexNumber].isCont=true;
								}
								else if(volumePixels[tp[jj][1].ix][tp[jj][1].iy][tp[jj][1].iz].isDone)
								{
									verts[vertexNumber].atomId=volumePixels[tp[jj][1].ix][tp[jj][1].iy][tp[jj][1].iz].atomId;
									verts[vertexNumber].isCont=false;
									if(volumePixels[tp[jj][1].ix][tp[jj][1].iy][tp[jj][1].iz].isBound)
										verts[vertexNumber].isCont=true;
								}
								
								vertexNumber++;
							}
							else
							{
								indv[jj]=volumePixelsInd[tp[jj][1].ix][curpt][tp[jj][1].iz].iy;
							}
						}
						else if(tp[jj][0].iz!=tp[jj][1].iz)
						{
							curpt=int(tv[jj](Z));
							if(volumePixelsInd[tp[jj][1].ix][tp[jj][1].iy][curpt].iz==-1)
							{
								indv[jj]=vertexNumber;
								volumePixelsInd[tp[jj][1].ix][tp[jj][1].iy][curpt].iz=vertexNumber;
								verts[vertexNumber].xyz(X)=tv[jj](X);
								verts[vertexNumber].xyz(Y)=tv[jj](Y);
								verts[vertexNumber].xyz(Z)=tv[jj](Z);
								if(volumePixels[tp[jj][0].ix][tp[jj][0].iy][tp[jj][0].iz].isDone)
								{
									verts[vertexNumber].atomId=volumePixels[tp[jj][0].ix][tp[jj][0].iy][tp[jj][0].iz].atomId;
									verts[vertexNumber].isCont=false;
									if(volumePixels[tp[jj][0].ix][tp[jj][0].iy][tp[jj][0].iz].isBound)
										verts[vertexNumber].isCont=true;
								}
								else if(volumePixels[tp[jj][1].ix][tp[jj][1].iy][tp[jj][1].iz].isDone)
								{
									verts[vertexNumber].atomId=volumePixels[tp[jj][1].ix][tp[jj][1].iy][tp[jj][1].iz].atomId;
									verts[vertexNumber].isCont=false;
									if(volumePixels[tp[jj][1].ix][tp[jj][1].iy][tp[jj][1].iz].isBound)
										verts[vertexNumber].isCont=true;
								}
								
								vertexNumber++;
							}
							else
							{
								indv[jj]=volumePixelsInd[tp[jj][1].ix][tp[jj][1].iy][curpt].iz;
							}
						}	
					}//jj	
					faces[faceNumber].abc(A)=indv[0];
					faces[faceNumber].abc(B)=indv[1];
					faces[faceNumber++].abc(C)=indv[2];
				}//ii
			}//k
		}//j
	}//i
	
	verts=(vertInfo *)realloc(verts,vertexNumber*sizeof(vertInfo));
	faces=(faceInfo *)realloc(faces,faceNumber*sizeof(faceInfo));
	
	for(i=0;i<pLength;i++)
	{
		for(j=0;j<pWidth;j++)
		{
			delete[] volumePixelsInd[i][j];
		}
		delete[] volumePixelsInd[i];
	}

	delete[] volumePixelsInd;

}

void ProteinSurface::marchingCube(int surfaceType)
{
	int i,j,k;
	marchingCubeInit(surfaceType);
	int ***vertSeq;
	vertSeq=new int**[pLength];

	for(i=0;i<pLength;i++)
	{
		vertSeq[i]=new int*[pWidth];
		for(j=0;j<pWidth;j++)
		{
			vertSeq[i][j]=new int[pHeight];
			for(k=0;k<pHeight;k++)
			{
				vertSeq[i][j][k]=-1;
			}		
		}
	}
	//Isn't this the same thing that marchingCubeInit does?

	if(faces!=NULL)
	{
		free(faces);
	}
	if(verts!=NULL)
	{
		free(verts);
	}
//	int allocFace=20;
//	int allocVert=12;
	int allocVert=4*(pHeight*pLength+pWidth*pLength+pHeight*pWidth);
	int allocFace=2*allocVert;
	faceNumber=0;
	vertexNumber=0;
	verts=new vertInfo[allocVert];
	faces=new faceInfo[allocFace];
	
	
	int sumType;
	int ii,jj,kk;
	int tp[6][3];
	/////////////////////////////////////////new added  normal is outer
	//face1
	for(i=0;i<1;i++)
	{
		for(j=0;j<pWidth-1;j++)
		{
			for(k=0;k<pHeight-1;k++)
			{
				if(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone
					&& volumePixels[i][j][k+1].isDone)
				{
					tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
					tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
					tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
				    tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
					for(ii=0;ii<4;ii++)
					{
						if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
						{
							vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
							verts[vertexNumber].xyz(X)=tp[ii][0];
							verts[vertexNumber].xyz(Y)=tp[ii][1];
							verts[vertexNumber].xyz(Z)=tp[ii][2];
							vertexNumber++;
						}
					}
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(C)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
				}
				else if((volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone)
					||( volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone)
					||( volumePixels[i][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone)
					||(volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone))
				{
					if(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone)
					{
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
						tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
						tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
					}
				    else if( volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone)
					{
						tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
						tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
					}
					else if( volumePixels[i][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone)
					{
						tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
					}
					else if(volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone)
					{
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
						tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
					}
					for(ii=0;ii<3;ii++)
					{
						if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
						{
							vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
							verts[vertexNumber].xyz(X)=tp[ii][0];
							verts[vertexNumber].xyz(Y)=tp[ii][1];
							verts[vertexNumber].xyz(Z)=tp[ii][2];
							vertexNumber++;
						}
					}
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(C)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
				}
				
			}
		}
	}
	//face3
	for(i=0;i<pLength-1;i++)
	{
		for(j=0;j<1;j++)
		{
			for(k=0;k<pHeight-1;k++)
			{
				if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone
					&& volumePixels[i][j][k+1].isDone)
				{
					tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
					tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
					tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
					tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
					for(ii=0;ii<4;ii++)
					{
						if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
						{
							vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
							verts[vertexNumber].xyz(X)=tp[ii][0];
							verts[vertexNumber].xyz(Y)=tp[ii][1];
							verts[vertexNumber].xyz(Z)=tp[ii][2];
							vertexNumber++;
						}
					}
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
					faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					faces[faceNumber++].abc(C)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
				}
				else if((volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone)
					||( volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone)
					||( volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone)
					||(volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone))
				{
					if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone)
					{
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
						tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
						tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
					}
					else if( volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone)
					{
						tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
						tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
					}
					else if( volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone)
					{
						tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
					}
					else if(volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone)
					{
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
						tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
					}
					for(ii=0;ii<3;ii++)
					{
						if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
						{
							vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
							verts[vertexNumber].xyz(X)=tp[ii][0];
							verts[vertexNumber].xyz(Y)=tp[ii][1];
							verts[vertexNumber].xyz(Z)=tp[ii][2];
							vertexNumber++;
						}
					}
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
					faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
				}
				
			}
		}
	}
	//face5
	for(i=0;i<pLength-1;i++)
	{
		for(j=0;j<pWidth-1;j++)
		{
			for(k=0;k<1;k++)
			{
				if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone
					&& volumePixels[i][j+1][k].isDone)
				{
					tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
					tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
					tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
					tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
					for(ii=0;ii<4;ii++)
					{
						if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
						{
							vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
							verts[vertexNumber].xyz(X)=tp[ii][0];
							verts[vertexNumber].xyz(Y)=tp[ii][1];
							verts[vertexNumber].xyz(Z)=tp[ii][2];
							vertexNumber++;
						}
					}
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(C)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
				}
				else if((volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone)
					||( volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone)
					||( volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j][k].isDone)
					||(volumePixels[i][j+1][k].isDone && volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone))
				{
					if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone)
					{
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
						tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
						tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
					}
					else if( volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone)
					{
						tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
						tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
						tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
					}
					else if( volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j][k].isDone)
					{
						tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
						tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
					}
					else if(volumePixels[i][j+1][k].isDone && volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone)
					{
						tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
						tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
					}
					for(ii=0;ii<3;ii++)
					{
						if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
						{
							vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
							verts[vertexNumber].xyz(X)=tp[ii][0];
							verts[vertexNumber].xyz(Y)=tp[ii][1];
							verts[vertexNumber].xyz(Z)=tp[ii][2];
							vertexNumber++;
						}
					}
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(C)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
				}
				
			}
		}
	}
	//face2
	for(i=pLength-1;i<pLength;i++)
	{
		for(j=0;j<pWidth-1;j++)
		{
			for(k=0;k<pHeight-1;k++)
			{
				if(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone
					&& volumePixels[i][j][k+1].isDone)
				{
					tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
					tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
					tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
					tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
					/*for(x = 0; x < 4; x++){
						tp[x][0] = i;
						tp[x][1] = j + x % 3;
						tp[x][2] = k + x/2;
					}
					*/
				
					for(ii=0;ii<4;ii++)
					{
						if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
						{
							vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
							verts[vertexNumber].xyz(X)=tp[ii][0];
							verts[vertexNumber].xyz(Y)=tp[ii][1];
							verts[vertexNumber].xyz(Z)=tp[ii][2];
							vertexNumber++;
						}
					}
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
					faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					faces[faceNumber++].abc(C)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
				}
				else if((volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone)
					||( volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone)
					||( volumePixels[i][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone)
					||(volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone))
				{
					if(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone)
					{
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
						tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
						tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
						/*
						for(x = 0; x < 3; x++){
							tp[x][0]= i;
							tp[x][1]= j + (x + 1) / 2;
							tp[x][2]= k + x / 2;
						}
						*/
					}
					else if( volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone)
					{
						tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
						tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
					}
					else if( volumePixels[i][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone)
					{
						tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
					}
					else if(volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone)
					{
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
						tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
					}
					for(ii=0;ii<3;ii++)
					{
						if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
						{
							vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
							verts[vertexNumber].xyz(X)=tp[ii][0];
							verts[vertexNumber].xyz(Y)=tp[ii][1];
							verts[vertexNumber].xyz(Z)=tp[ii][2];
							vertexNumber++;
						}
					}
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
					faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
				}
				
			}
		}
	}
	//face4
	for(i=0;i<pLength-1;i++)
	{
		for(j=pWidth-1;j<pWidth;j++)
		{
			for(k=0;k<pHeight-1;k++)
			{
				if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone
					&& volumePixels[i][j][k+1].isDone)
				{
					tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
					tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
					tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
					tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
					for(ii=0;ii<4;ii++)
					{
						if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
						{
							vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
							verts[vertexNumber].xyz(X)=tp[ii][0];
							verts[vertexNumber].xyz(Y)=tp[ii][1];
							verts[vertexNumber].xyz(Z)=tp[ii][2];
							vertexNumber++;
						}
					}
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(C)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
				}
				else if((volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone)
					||( volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone)
					||( volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone)
					||(volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone))
				{
					if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone)
					{
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
						tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
						tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
					}
					else if( volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone)
					{
						tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
						tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
					}
					else if( volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone)
					{
						tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
					}
					else if(volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone)
					{
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
						tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
					}
					for(ii=0;ii<3;ii++)
					{
						if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
						{
							vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
							verts[vertexNumber].xyz(X)=tp[ii][0];
							verts[vertexNumber].xyz(Y)=tp[ii][1];
							verts[vertexNumber].xyz(Z)=tp[ii][2];
							vertexNumber++;
						}
					}
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(C)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
				}
				
			}
		}
	}
	//face6
	for(i=0;i<pLength-1;i++)
	{
		for(j=0;j<pWidth-1;j++)
		{
			for(k=pHeight-1;k<pHeight;k++)
			{
				if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone
					&& volumePixels[i][j+1][k].isDone)
				{
					tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
					tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
					tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
					tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
					for(ii=0;ii<4;ii++)
					{
						if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
						{
							vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
							verts[vertexNumber].xyz(X)=tp[ii][0];
							verts[vertexNumber].xyz(Y)=tp[ii][1];
							verts[vertexNumber].xyz(Z)=tp[ii][2];
							vertexNumber++;
						}
					}
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
					faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					faces[faceNumber++].abc(C)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
				}
				else if((volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone)
					||( volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone)
					||( volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j][k].isDone)
					||(volumePixels[i][j+1][k].isDone && volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone))
				{
					if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone)
					{
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
						tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
						tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
					}
					else if( volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone)
					{
						tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
						tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
						tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
					}
					else if( volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j][k].isDone)
					{
						tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
						tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
					}
					else if(volumePixels[i][j+1][k].isDone && volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone)
					{
						tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
						tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
					}
					for(ii=0;ii<3;ii++)
					{
						if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
						{
							vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
							verts[vertexNumber].xyz(X)=tp[ii][0];
							verts[vertexNumber].xyz(Y)=tp[ii][1];
							verts[vertexNumber].xyz(Z)=tp[ii][2];
							vertexNumber++;
						}
					}
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
					faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
				}
				
			}
		}
	}

	///////////////////////////////////////////
	for(i=0;i<pLength-1;i++)
	{
		for(j=0;j<pWidth-1;j++)
		{
			for(k=0;k<pHeight-1;k++)
			{
				sumType=0;		
				for( ii=0;ii<2;ii++)
				{
					for( jj=0;jj<2;jj++)
					{
						for( kk=0;kk<2;kk++)
						{
							if(volumePixels[i+ii][j+jj][k+kk].isDone)
								sumType++;
						}
					}
				}//ii
				if(vertexNumber+6>allocVert)
				{
					allocVert*=2;
					verts=(vertInfo *)realloc(verts,allocVert*sizeof(vertInfo));
				}
				if(faceNumber+3>allocFace)
				{
					allocFace*=2;
					faces=(faceInfo *)realloc(faces,allocFace*sizeof(faceInfo));
				}


				//If sumType is 0 or 1 or 2 or 8, do nothing
				
				else if(sumType==3)
				{
					if((volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone)
					   ||(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i+1][j+1][k].isDone)
					   ||(volumePixels[i][j+1][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone)
					   ||(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i+1][j][k].isDone)
					   ||(volumePixels[i][j][k+1].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone)
					   ||(volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone)
					   ||(volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone)
					   ||(volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j][k+1].isDone)
					   ||(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone)
					   ||(volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone)
					   ||(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone)
					   ||(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i][j][k+1].isDone)
					   ||(volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k+1].isDone)
					   ||(volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone)
					   ||(volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone)
					   ||(volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone)
					   ||(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone)
					   ||(volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone)
					   ||(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone)
					   ||(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i+1][j+1][k+1].isDone)
					   ||(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j][k+1].isDone)
					   ||(volumePixels[i][j][k].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone)
					   ||(volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone)
					   ||(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone))
					{
						if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;	
						}//11
						else if(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i+1][j+1][k].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
						}//12
						else if(volumePixels[i][j+1][k].isDone && volumePixels[i+1][j][k].isDone&& volumePixels[i+1][j+1][k].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
						}//13
						else if(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone&& volumePixels[i+1][j][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
						}//14
						else if(volumePixels[i][j][k+1].isDone && volumePixels[i+1][j][k+1].isDone&& volumePixels[i+1][j+1][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
						}//21
						else if(volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone&& volumePixels[i+1][j+1][k+1].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
						}//22
						else if(volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j][k+1].isDone&& volumePixels[i+1][j+1][k+1].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
						}//23
						else if(volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone&& volumePixels[i+1][j][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
						}//24
						else if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone&& volumePixels[i+1][j][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
						}//31
						else if(volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
						}//32
						else if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
						}//33
						else if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i][j][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
						}//34
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
						}//41
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
						}//42
						else if(volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
						}//43
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k].isDone	&& volumePixels[i+1][j][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
						}//44
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone ) 
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
						}//51
						else if( volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
						}//52
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
						}//53
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i+1][j+1][k+1].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
						}//54
						else if(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone	&& volumePixels[i][j][k+1].isDone ) 
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
						}//61
						else if(volumePixels[i][j][k].isDone 	&& volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
						}//62
						else if(volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
						}//63
						else if(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone	&& volumePixels[i][j+1][k+1].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
						}//64
						for(ii=0;ii<3;ii++)
						{
							if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
							{
								vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
								verts[vertexNumber].xyz(X)=tp[ii][0];
								verts[vertexNumber].xyz(Y)=tp[ii][1];
								verts[vertexNumber].xyz(Z)=tp[ii][2];
								vertexNumber++;
							}
						}
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					}//no5 24
				}//total3
				else if(sumType==4)
				{
					if((volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone
						&& volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone) 
						|| (volumePixels[i][j][k+1].isDone && volumePixels[i+1][j][k+1].isDone
						&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k+1].isDone)
						|| (volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone
						&& volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone)
						|| (volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k].isDone
						&& volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone)
						|| (volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone
						&& volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone)
						|| (volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone
						&& volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone))
					{
						if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone
							&& volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
							
						}
						else if (volumePixels[i][j][k+1].isDone && volumePixels[i+1][j][k+1].isDone
							&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
						}
						else if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone
							&& volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
						}
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k].isDone
							&& volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
						}
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone
							&& volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
						}
						else if(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone
							&& volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
						}
						for(ii=0;ii<4;ii++)
						{
							if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
							{
								vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
								verts[vertexNumber].xyz(X)=tp[ii][0];
								verts[vertexNumber].xyz(Y)=tp[ii][1];
								verts[vertexNumber].xyz(Z)=tp[ii][2];
								vertexNumber++;
							}
						}
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
					}//no.8 6
									
				  else if((volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone  && volumePixels[i][j+1][k+1].isDone)//11
					   ||(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k+1].isDone)//12
					   ||(volumePixels[i][j+1][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone && volumePixels[i][j][k+1].isDone)//13
					   ||(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k+1].isDone)//14
					   ||(volumePixels[i][j][k+1].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k].isDone)//21
					   ||(volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i+1][j][k].isDone)//22
					   ||(volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k].isDone)//23
					   ||(volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k].isDone)//24
					   ||(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone)//31
					   ||(volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k].isDone)//32
					   ||(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i+1][j+1][k].isDone)//33
					   ||(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone)//34
					   ||(volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone)//41
					   ||(volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k].isDone)//42
					   ||(volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k].isDone)//43
					   ||(volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone  && volumePixels[i][j+1][k+1].isDone)//44
					   ||(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone  && volumePixels[i+1][j][k+1].isDone)//51
					   ||( volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i+1][j][k].isDone)//52
					   ||(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k].isDone)//53
					   ||(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone)//54
					   ||(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j][k+1].isDone  && volumePixels[i+1][j+1][k+1].isDone)//61
					   ||(volumePixels[i][j][k].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k].isDone)//62
					   ||(volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j][k].isDone)//63
					   ||(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone&& volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j][k+1].isDone))
				   {
						if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k].isDone  && volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;	
						}//11
						else if(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
						}//12
						else if(volumePixels[i][j+1][k].isDone && volumePixels[i+1][j][k].isDone&& volumePixels[i+1][j+1][k].isDone && volumePixels[i][j][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
						}//13
						else if(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone&& volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
						}//14
						else if(volumePixels[i][j][k+1].isDone && volumePixels[i+1][j][k+1].isDone&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
						}//21
						else if(volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i+1][j][k].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
						}//22
						else if(volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j][k+1].isDone&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
						}//23
						else if(volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone&& volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
						}//24
						else if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone&& volumePixels[i+1][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
						}//31
						else if(volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
						}//32
						else if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i+1][j+1][k].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
						}//33
						else if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
						}//34
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
						}//41
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
						}//42
						else if(volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
						}//43
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k].isDone	&& volumePixels[i+1][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
						}//44
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j][k+1].isDone ) 
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
						}//51
						else if( volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i+1][j][k].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
						}//52
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
						}//53
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
						}//54
						else if(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone	&& volumePixels[i][j][k+1].isDone && volumePixels[i+1][j+1][k+1].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
						}//61
						else if(volumePixels[i][j][k].isDone 	&& volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
						}//62
						else if(volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j][k].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
						}//63
						else if(volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone	&& volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j][k+1].isDone) 
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
						}//64
						for(ii=0;ii<3;ii++)
						{
							if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
							{
								vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
								verts[vertexNumber].xyz(X)=tp[ii][0];
								verts[vertexNumber].xyz(Y)=tp[ii][1];
								verts[vertexNumber].xyz(Z)=tp[ii][2];
								vertexNumber++;
							}
						}
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
				   }//no12 24
					else if((volumePixels[i][j][k].isDone && volumePixels[i][j+1][k+1].isDone
						&& volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k].isDone)
						|| (volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone
						&& volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k+1].isDone)
						|| (volumePixels[i][j][k].isDone && volumePixels[i][j][k+1].isDone
						&& volumePixels[i+1][j][k].isDone && volumePixels[i][j+1][k].isDone)
						|| (volumePixels[i][j+1][k].isDone && volumePixels[i+1][j][k].isDone
						&& volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j+1][k+1].isDone)
						|| (volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone
						&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k].isDone)
						|| (volumePixels[i][j][k+1].isDone && volumePixels[i+1][j][k].isDone
						&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i+1][j][k+1].isDone)
						|| (volumePixels[i][j][k].isDone && volumePixels[i][j][k+1].isDone
						&& volumePixels[i+1][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone)
						|| (volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j][k+1].isDone
						&& volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j+1][k+1].isDone))
					{
						if(volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone
							&& volumePixels[i][j][k].isDone && volumePixels[i+1][j+1][k].isDone )
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
						}//1
						else if(volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone
							&& volumePixels[i+1][j+1][k].isDone && volumePixels[i][j][k].isDone )
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
						}//2
						else if(volumePixels[i][j][k].isDone && volumePixels[i][j][k+1].isDone
							&& volumePixels[i+1][j][k].isDone && volumePixels[i][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
						}//3
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j+1][k+1].isDone
							&& volumePixels[i][j+1][k].isDone && volumePixels[i+1][j][k].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
						}//4
						else if(volumePixels[i][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone
							&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
						}//5
						else if(volumePixels[i+1][j][k].isDone && volumePixels[i+1][j][k+1].isDone
							&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
						}//6
						else if(volumePixels[i][j][k].isDone && volumePixels[i][j][k+1].isDone
							&& volumePixels[i+1][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
						}//7
						else if(volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j][k+1].isDone
							&& volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
						}//8
						for(ii=0;ii<3;ii++)
						{
							if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
							{
								vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
								verts[vertexNumber].xyz(X)=tp[ii][0];
								verts[vertexNumber].xyz(Y)=tp[ii][1];
								verts[vertexNumber].xyz(Z)=tp[ii][2];
								vertexNumber++;
							}
						}
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					}// no.9 8
					else if((volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone
						&& volumePixels[i+1][j+1][k].isDone && volumePixels[i][j][k+1].isDone)
						||(volumePixels[i][j+1][k].isDone && volumePixels[i+1][j][k].isDone
						&& volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k+1].isDone)
						||(volumePixels[i][j+1][k].isDone && volumePixels[i][j][k].isDone
						&& volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j+1][k+1].isDone)
						||(volumePixels[i][j+1][k].isDone && volumePixels[i][j][k].isDone
						&& volumePixels[i+1][j][k].isDone && volumePixels[i][j+1][k+1].isDone)
						||(volumePixels[i][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone
						&& volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j][k].isDone)
						||(volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone
						&& volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k].isDone)
						||(volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k+1].isDone
						&& volumePixels[i+1][j][k+1].isDone && volumePixels[i][j+1][k].isDone)
						||(volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k+1].isDone
						&& volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone)
						||(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone
						&& volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k].isDone)
						||(volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k].isDone
						&& volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k].isDone)
						||(volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k].isDone
						&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i+1][j][k].isDone)
						||(volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k].isDone
						&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i+1][j][k].isDone))
					{
						if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone
							&& volumePixels[i+1][j+1][k].isDone && volumePixels[i][j][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;	
							tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
						}//1
						else if(volumePixels[i][j+1][k].isDone && volumePixels[i+1][j][k].isDone
							&& volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;	
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
						}//2
						else if(volumePixels[i][j+1][k].isDone && volumePixels[i][j][k].isDone
							&& volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;	
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
						}//3
						else if(volumePixels[i][j+1][k].isDone && volumePixels[i][j][k].isDone
							&& volumePixels[i+1][j][k].isDone && volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;	
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
						}//4
						else if(volumePixels[i][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone
							&& volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j][k].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;	
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
						}//5
						else if(volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone
							&& volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j+1][k].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;	
							tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
						}//6
						else if(volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k+1].isDone
							&& volumePixels[i+1][j][k+1].isDone && volumePixels[i][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;	
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
						}//7
						else if(volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k+1].isDone
							&& volumePixels[i][j][k+1].isDone && volumePixels[i][j][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;	
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
						}//8
						else if(volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone
							&& volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;	
							tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
						}//9
						else if(volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k].isDone
							&& volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;	
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
						}//10
						else if(volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k].isDone
							&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i+1][j][k].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;	
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
						}//11
						else if(volumePixels[i][j+1][k+1].isDone && volumePixels[i+1][j+1][k].isDone
							&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i+1][j][k].isDone) 
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;	
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
						}//12
						for(ii=0;ii<4;ii++)
						{
							if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
							{
								vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
								verts[vertexNumber].xyz(X)=tp[ii][0];
								verts[vertexNumber].xyz(Y)=tp[ii][1];
								verts[vertexNumber].xyz(Z)=tp[ii][2];
								vertexNumber++;
							}
						}
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
					}//no.11 12
					else if((volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone
						&& volumePixels[i][j+1][k].isDone && volumePixels[i+1][j][k+1].isDone)
						||(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone
						&& volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j+1][k+1].isDone)
						||(volumePixels[i][j+1][k].isDone && volumePixels[i+1][j][k].isDone
						&& volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone)
						||(volumePixels[i][j+1][k].isDone && volumePixels[i][j][k].isDone
						&& volumePixels[i+1][j+1][k].isDone && volumePixels[i][j][k+1].isDone)
						||(volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone
						&& volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k].isDone)
						||(volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k+1].isDone
						&& volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j][k].isDone)
						||(volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k+1].isDone
						&& volumePixels[i][j][k+1].isDone && volumePixels[i+1][j+1][k].isDone)
						||(volumePixels[i+1][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone
						&& volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k].isDone)
						||(volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k+1].isDone
						&& volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone)
						||(volumePixels[i+1][j][k].isDone && volumePixels[i][j][k].isDone
						&& volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone)
						||(volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone
						&& volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k].isDone)
						||(volumePixels[i][j+1][k].isDone && volumePixels[i+1][j+1][k].isDone
						&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i+1][j][k+1].isDone))
					{
						if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone
							&& volumePixels[i][j+1][k].isDone && volumePixels[i+1][j][k+1].isDone)  
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;	
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
						}//1
						else if(volumePixels[i][j][k].isDone && volumePixels[i+1][j][k].isDone
							&& volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j+1][k+1].isDone)  
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;	
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
						}//2
						else if(volumePixels[i][j+1][k].isDone && volumePixels[i+1][j][k].isDone
							&& volumePixels[i+1][j+1][k].isDone && volumePixels[i][j+1][k+1].isDone)  
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;	
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
						}//3
						else if(volumePixels[i][j+1][k].isDone && volumePixels[i][j][k].isDone
							&& volumePixels[i+1][j+1][k].isDone && volumePixels[i][j][k+1].isDone)  
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;	
							tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
						}//4
						else if(volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j][k+1].isDone
							&& volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k].isDone)  
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;	
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
						}//5
						else if(volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k+1].isDone
							&& volumePixels[i+1][j][k+1].isDone && volumePixels[i+1][j][k].isDone)  
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;	
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
						}//6
						else if(volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k+1].isDone
							&& volumePixels[i][j][k+1].isDone && volumePixels[i+1][j+1][k].isDone)  
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;	
							tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
						}//7
						else if(volumePixels[i+1][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone
							&& volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k].isDone)  
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;	
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
						}//8
						else if(volumePixels[i+1][j+1][k+1].isDone && volumePixels[i][j+1][k+1].isDone
							&& volumePixels[i][j][k].isDone && volumePixels[i][j+1][k].isDone)  
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;	
							tp[3][0]=i;tp[3][1]=j;tp[3][2]=k;
						}//9
						else if(volumePixels[i+1][j][k].isDone && volumePixels[i][j][k].isDone
							&& volumePixels[i][j][k+1].isDone && volumePixels[i][j+1][k+1].isDone)  
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;	
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k;
						}//10
						else if(volumePixels[i+1][j][k+1].isDone && volumePixels[i][j][k+1].isDone
							&& volumePixels[i+1][j+1][k].isDone && volumePixels[i+1][j][k].isDone)  
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;	
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k;
						}//11
						else if(volumePixels[i][j+1][k].isDone && volumePixels[i+1][j+1][k].isDone
							&& volumePixels[i+1][j+1][k+1].isDone && volumePixels[i+1][j][k+1].isDone)  
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;	
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
						}//12
						for(ii=0;ii<4;ii++)
						{
							if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
							{
								vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
								verts[vertexNumber].xyz(X)=tp[ii][0];
								verts[vertexNumber].xyz(Y)=tp[ii][1];
								verts[vertexNumber].xyz(Z)=tp[ii][2];
								vertexNumber++;
							}
						}
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
					}//no.14 12
				}//total4
				else if(sumType==5)
				{
					if((!volumePixels[i+1][j][k].isDone && !volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						|| (!volumePixels[i][j+1][k].isDone && !volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						|| (!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone)
						|| (!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone)
						|| (!volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j][k].isDone && !volumePixels[i+1][j+1][k].isDone)
						|| (!volumePixels[i][j+1][k+1].isDone && !volumePixels[i][j][k].isDone && !volumePixels[i+1][j+1][k].isDone)
						|| (!volumePixels[i+1][j+1][k+1].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i][j+1][k].isDone)
						|| (!volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i][j+1][k].isDone))
					{		
						if(!volumePixels[i+1][j][k].isDone && !volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
						}//1
						else if(!volumePixels[i][j+1][k].isDone && !volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
						}//2
						else if(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
						}//3
						else if(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
						}//4
						else if(!volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j][k].isDone && !volumePixels[i+1][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
						}//5
						else if(!volumePixels[i][j+1][k+1].isDone && !volumePixels[i][j][k].isDone && !volumePixels[i+1][j+1][k].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
						}//6
						else if(!volumePixels[i+1][j+1][k+1].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
						}//7
						else if(!volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i][j+1][k].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
						}//8
						for(ii=0;ii<3;ii++)
						{
							if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
							{
								vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
								verts[vertexNumber].xyz(X)=tp[ii][0];
								verts[vertexNumber].xyz(Y)=tp[ii][1];
								verts[vertexNumber].xyz(Z)=tp[ii][2];
								vertexNumber++;
							}
						}
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					}//no.7 8
					else if((!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k].isDone)
				   ||(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j+1][k].isDone)
				   ||(!volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k].isDone)
				   ||(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j][k].isDone)
				   ||(!volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
				   ||(!volumePixels[i][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
				   ||(!volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
				   ||(!volumePixels[i][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j][k+1].isDone)
				   ||(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j][k+1].isDone)
				   ||(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j][k+1].isDone)
				   ||(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j][k+1].isDone)
				   ||(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i][j][k+1].isDone)
				   ||(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k+1].isDone)
				   ||(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
				   ||(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
				   ||(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j][k+1].isDone )
				   ||(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone )
				   ||(!volumePixels[i][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
				   ||(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
				   ||(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j+1][k+1].isDone)
				   ||(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i][j][k+1].isDone )
				   ||(!volumePixels[i][j][k].isDone && !volumePixels[i][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone)
				   ||(!volumePixels[i][j][k+1].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone)
				   ||(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone))
				{
					if(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k].isDone)
					{
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
						tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
						tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
						tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
					}//11
					else if(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j+1][k].isDone) 
					{
						tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
						tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
						tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k;
					}//12
					else if(!volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j][k].isDone&& !volumePixels[i+1][j+1][k].isDone) 
					{
						tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
						tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
						tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
						tp[3][0]=i;tp[3][1]=j;tp[3][2]=k;
					}//13
					else if(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k].isDone&& !volumePixels[i+1][j][k].isDone)
					{
						tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
						tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
						tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k;
					}//14
					else if(!volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j][k+1].isDone&& !volumePixels[i+1][j+1][k+1].isDone) 
					{
						tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
						tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
						tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
					}//21
					else if(!volumePixels[i][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone&& !volumePixels[i+1][j+1][k+1].isDone) 
					{
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
						tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
						tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
						tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
					}//22
					else if(!volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j][k+1].isDone&& !volumePixels[i+1][j+1][k+1].isDone) 
					{
						tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
						tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
						tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
						tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
					}//23
					else if(!volumePixels[i][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone&& !volumePixels[i+1][j][k+1].isDone) 
					{
						tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
						tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
						tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
					}//24
					else if(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k].isDone&& !volumePixels[i+1][j][k+1].isDone) 
					{
						tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
						tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
						tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
						tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
					}//31
					else if(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j][k+1].isDone) 
					{
						tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
						tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
						tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
						tp[3][0]=i;tp[3][1]=j;tp[3][2]=k;
					}//32
					else if(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j][k+1].isDone) 
					{
						tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
						tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
						tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
						tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k;
					}//33
					else if(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i][j][k+1].isDone) 
					{
						tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
						tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
						tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
						tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
					}//34
					else if(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k+1].isDone) 
					{
						tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
						tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
						tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
					}//41
					else if(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone) 
					{
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
						tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
						tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
						tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k;
					}//42
					else if(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone) 
					{
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
						tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
						tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k;
					}//43
					else if(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j][k+1].isDone) 
					{
						tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
						tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
					}//44
					else if(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone ) 
					{
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
						tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
						tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
					}//51
					else if( !volumePixels[i][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone) 
					{
						tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
						tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k;
					}//52
					else if(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone) 
					{
						tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
						tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
						tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
					}//53
					else if(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j+1][k+1].isDone) 
					{
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
						tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
						tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
						tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
					}//54
					else if(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i][j][k+1].isDone ) 
					{
						tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
						tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
						tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
						tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
					}//61
					else if(!volumePixels[i][j][k].isDone && !volumePixels[i][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone) 
					{
						tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
						tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
						tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
						tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
					}//62
					else if(!volumePixels[i][j][k+1].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone) 
					{
						tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
						tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
						tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
						tp[3][0]=i;tp[3][1]=j;tp[3][2]=k;
					}//63
					else if(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone) 
					{
						tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
						tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
						tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
						tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
					}//64
					for(ii=0;ii<4;ii++)
					{
						if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
						{
							vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
							verts[vertexNumber].xyz(X)=tp[ii][0];
							verts[vertexNumber].xyz(Y)=tp[ii][1];
							verts[vertexNumber].xyz(Z)=tp[ii][2];
							vertexNumber++;
						}
					}
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
					faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
					faces[faceNumber++].abc(C)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
				}//no5 24
					else if((!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k+1].isDone)//1
						||(!volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j][k+1].isDone)//2
						||(!volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone && !volumePixels[i+1][j][k].isDone)//3
						||(!volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i+1][j+1][k].isDone)//4
						||(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j+1][k+1].isDone)//5
						||(!volumePixels[i+1][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone && !volumePixels[i][j+1][k].isDone)//6
						||(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone)//7
						||(!volumePixels[i][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j+1][k].isDone)//8
						||(!volumePixels[i][j][k].isDone && !volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)//9
						||(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j+1][k+1].isDone && !volumePixels[i][j][k].isDone)//10
						||(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone)//11
						||(!volumePixels[i][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j][k+1].isDone))
					{
						if(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
							tp[4][0]=i;tp[4][1]=j+1;tp[4][2]=k+1;
						}//1
						else if(!volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[4][0]=i+1;tp[4][1]=j;tp[4][2]=k+1;
						}//2
						else if(!volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone && !volumePixels[i+1][j][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k;
							tp[4][0]=i;tp[4][1]=j;tp[4][2]=k;
						}//3
						else if(!volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i+1][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[4][0]=i;tp[4][1]=j+1;tp[4][2]=k;
						}//4
						else if(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						{
						    //tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							//tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							//tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							//tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[4][0]=i+1;tp[4][1]=j;tp[4][2]=k+1;
						}//5
						else if(!volumePixels[i+1][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone && !volumePixels[i][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[4][0]=i;tp[4][1]=j;tp[4][2]=k;
						}//6
						else if(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[4][0]=i;tp[4][1]=j;tp[4][2]=k+1;
						}//7
						else if(!volumePixels[i][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
							tp[4][0]=i+1;tp[4][1]=j;tp[4][2]=k;
						}//8
						else if(!volumePixels[i][j][k].isDone && !volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[4][0]=i+1;tp[4][1]=j+1;tp[4][2]=k;
						}//9
						else if(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j+1][k+1].isDone && !volumePixels[i][j][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
							tp[4][0]=i;tp[4][1]=j;tp[4][2]=k+1;
						}//10
						else if(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
							tp[4][0]=i;tp[4][1]=j+1;tp[4][2]=k;
						}//11
						else if(!volumePixels[i][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[4][0]=i+1;tp[4][1]=j;tp[4][2]=k;
						}//12
						for(ii=0;ii<5;ii++)
						{
							if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
							{
								vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
								verts[vertexNumber].xyz(X)=tp[ii][0];
								verts[vertexNumber].xyz(Y)=tp[ii][1];
								verts[vertexNumber].xyz(Z)=tp[ii][2];
								vertexNumber++;
							}
						}
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
						faces[faceNumber].abc(A)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[4][0]][tp[4][1]][tp[4][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
						
					}//no.6 12-1
					else if((!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i][j+1][k+1].isDone)//1
						||(!volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j][k+1].isDone)//2
						||(!volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone && !volumePixels[i][j][k].isDone)//3
						||(!volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j+1][k].isDone)//4
						||(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j][k+1].isDone)//5
						||(!volumePixels[i+1][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone && !volumePixels[i][j][k].isDone)//6
						||(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j][k+1].isDone)//7
						||(!volumePixels[i][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j][k].isDone)//8
						||(!volumePixels[i][j][k].isDone && !volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j+1][k].isDone)//9
						||(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j+1][k+1].isDone && !volumePixels[i][j][k+1].isDone)//10
						||(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j+1][k].isDone)//11
						||(!volumePixels[i][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j][k].isDone))
					{
						if(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k].isDone && !volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
							tp[4][0]=i+1;tp[4][1]=j+1;tp[4][2]=k+1;
						}//1
						else if(!volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[4][0]=i;tp[4][1]=j;tp[4][2]=k+1;
						}//2
						else if(!volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone && !volumePixels[i][j][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k;
							tp[4][0]=i+1;tp[4][1]=j;tp[4][2]=k;
						}//3
						else if(!volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[4][0]=i+1;tp[4][1]=j+1;tp[4][2]=k;
						}//4
						else if(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j][k+1].isDone)
						{
							//tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							//tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							//tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							//tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[4][0]=i+1;tp[4][1]=j+1;tp[4][2]=k+1;
						}//5
						else if(!volumePixels[i+1][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone && !volumePixels[i][j][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[4][0]=i;tp[4][1]=j+1;tp[4][2]=k;
						}//6
						else if(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[4][0]=i;tp[4][1]=j+1;tp[4][2]=k+1;
						}//7
						else if(!volumePixels[i][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
							tp[4][0]=i+1;tp[4][1]=j+1;tp[4][2]=k;
						}//8
						else if(!volumePixels[i][j][k].isDone && !volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[4][0]=i+1;tp[4][1]=j+1;tp[4][2]=k+1;
						}//9
						else if(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j+1][k+1].isDone && !volumePixels[i][j][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
							tp[4][0]=i;tp[4][1]=j;tp[4][2]=k;
						}//10
						else if(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
							tp[4][0]=i;tp[4][1]=j+1;tp[4][2]=k+1;
						}//11
						else if(!volumePixels[i][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j][k].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[4][0]=i+1;tp[4][1]=j;tp[4][2]=k+1;
						}//12
						for(ii=0;ii<5;ii++)
						{
							if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
							{
								vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
								verts[vertexNumber].xyz(X)=tp[ii][0];
								verts[vertexNumber].xyz(Y)=tp[ii][1];
								verts[vertexNumber].xyz(Z)=tp[ii][2];
								vertexNumber++;
							}
						}
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[4][0]][tp[4][1]][tp[4][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
						
					}//no.6 12-2

				}//total5
				
				else if(sumType==6)
				{
					if((!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k].isDone)
						||(!volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j+1][k].isDone)
						||(!volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						||(!volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j][k+1].isDone)
						||(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k].isDone)
						||(!volumePixels[i+1][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						||(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k].isDone)
						||(!volumePixels[i][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone)
						||(!volumePixels[i][j][k].isDone && !volumePixels[i][j][k+1].isDone)
						||(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						||(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j][k+1].isDone)
						||(!volumePixels[i][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone))
					{
						if(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
						}//1
						else if(!volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
						}//2
						else if(!volumePixels[i][j+1][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
						}//3
						else if(!volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
						}//4
						else if(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k;
						}//5
						else if(!volumePixels[i+1][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
						}//6
						else if(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
						}//7
						else if(!volumePixels[i][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
						}//8
						else if(!volumePixels[i][j][k].isDone && !volumePixels[i][j][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k;
						}//9
						else if(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
						}//10
						else if(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k;
						}//11
						else if(!volumePixels[i][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
						}//12
						for(ii=0;ii<4;ii++)
						{
							if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
							{
								vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
								verts[vertexNumber].xyz(X)=tp[ii][0];
								verts[vertexNumber].xyz(Y)=tp[ii][1];
								verts[vertexNumber].xyz(Z)=tp[ii][2];
								vertexNumber++;
							}
						}
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
						
					}//no.2 12	
					
					else if((!volumePixels[i][j][k].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						||(!volumePixels[i+1][j][k].isDone && !volumePixels[i][j+1][k+1].isDone)
						||(!volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j][k+1].isDone)
						||(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j][k+1].isDone))
					{
						if(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
							tp[4][0]=i;tp[4][1]=j+1;tp[4][2]=k;
							tp[5][0]=i+1;tp[5][1]=j;tp[5][2]=k;
						}//1
						else if(!volumePixels[i+1][j][k].isDone && !volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
							tp[4][0]=i;tp[4][1]=j;tp[4][2]=k;
							tp[5][0]=i+1;tp[5][1]=j+1;tp[5][2]=k;
						}//2
						else if(!volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[4][0]=i+1;tp[4][1]=j+1;tp[4][2]=k;
							tp[5][0]=i;tp[5][1]=j;tp[5][2]=k;
						}//3
						else if(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
							tp[4][0]=i+1;tp[4][1]=j;tp[4][2]=k;
							tp[5][0]=i;tp[5][1]=j+1;tp[5][2]=k;
						}//4
						for(ii=0;ii<6;ii++)
						{
							if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
							{
								vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
								verts[vertexNumber].xyz(X)=tp[ii][0];
								verts[vertexNumber].xyz(Y)=tp[ii][1];
								verts[vertexNumber].xyz(Z)=tp[ii][2];
								vertexNumber++;
							}
						}
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber].abc(A)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[4][0]][tp[4][1]][tp[4][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[5][0]][tp[5][1]][tp[5][2]];
					}//no.4 4
					
					else if((!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k+1].isDone)
						||(!volumePixels[i+1][j][k].isDone && !volumePixels[i][j][k+1].isDone)
						||(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						||(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j][k+1].isDone)
						||(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone)
						||(!volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						||(!volumePixels[i][j+1][k].isDone && !volumePixels[i][j][k+1].isDone)
						||(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k+1].isDone)
						||(!volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						||(!volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone)
						||(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j+1][k].isDone)
						||(!volumePixels[i+1][j][k].isDone && !volumePixels[i][j+1][k].isDone))
					{
						if(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k+1;
						}//1
						else if(!volumePixels[i+1][j][k].isDone && !volumePixels[i][j][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k;
						}//2
						else if(!volumePixels[i+1][j][k].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k+1;
						}//3
						else if(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i+1][j][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
						}//4
						else if(!volumePixels[i+1][j+1][k].isDone && !volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
						}//5
						else if(!volumePixels[i][j+1][k].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;
							tp[3][0]=i;tp[3][1]=j;tp[3][2]=k;
						}//6
						else if(!volumePixels[i][j+1][k].isDone && !volumePixels[i][j][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
						}//7
						else if(!volumePixels[i][j][k].isDone && !volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k;
						}//8
						else if(!volumePixels[i][j][k+1].isDone && !volumePixels[i+1][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
							tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
							tp[3][0]=i+1;tp[3][1]=j+1;tp[3][2]=k;
						}//9
						else if(!volumePixels[i+1][j][k+1].isDone && !volumePixels[i][j+1][k+1].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
							tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
							tp[3][0]=i;tp[3][1]=j+1;tp[3][2]=k;
						}//10
						else if(!volumePixels[i][j][k].isDone && !volumePixels[i+1][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
							tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
							tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[3][0]=i;tp[3][1]=j;tp[3][2]=k+1;
						}//11
						else if(!volumePixels[i+1][j][k].isDone && !volumePixels[i][j+1][k].isDone)
						{
							tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
							tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
							tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;
							tp[3][0]=i+1;tp[3][1]=j;tp[3][2]=k+1;
						}//12
						for(ii=0;ii<4;ii++)
						{
							if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
							{
								vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
								verts[vertexNumber].xyz(X)=tp[ii][0];
								verts[vertexNumber].xyz(Y)=tp[ii][1];
								verts[vertexNumber].xyz(Z)=tp[ii][2];
								vertexNumber++;
							}
						}
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
						faces[faceNumber].abc(B)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
						faces[faceNumber++].abc(C)=vertSeq[tp[3][0]][tp[3][1]][tp[3][2]];
					}//no.3 12
					
				}//total6
				
				else if(sumType==7)
				{
					if(!volumePixels[i][j][k].isDone)
					{
						tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k;
						tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k;
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k+1;
					}//1
					else if(!volumePixels[i+1][j][k].isDone)
					{
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k;
						tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k;		
						tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k+1;
					}//2
					else if(!volumePixels[i+1][j+1][k].isDone)
					{
						tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k;
						tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k;		
						tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k+1;
					}//3
					else if(!volumePixels[i][j+1][k].isDone)
					{				
						tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k;
						tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k+1;
					}//4
					else if(!volumePixels[i][j][k+1].isDone)
					{
						tp[0][0]=i+1;tp[0][1]=j;tp[0][2]=k+1;
						tp[1][0]=i;tp[1][1]=j+1;tp[1][2]=k+1;		
						tp[2][0]=i;tp[2][1]=j;tp[2][2]=k;
					}//5
					else if(!volumePixels[i+1][j][k+1].isDone)
					{
						tp[0][0]=i+1;tp[0][1]=j+1;tp[0][2]=k+1;
						tp[1][0]=i;tp[1][1]=j;tp[1][2]=k+1;		
						tp[2][0]=i+1;tp[2][1]=j;tp[2][2]=k;
					}//6
					else if(!volumePixels[i+1][j+1][k+1].isDone)
					{
						tp[0][0]=i;tp[0][1]=j+1;tp[0][2]=k+1;
						tp[1][0]=i+1;tp[1][1]=j;tp[1][2]=k+1;			
						tp[2][0]=i+1;tp[2][1]=j+1;tp[2][2]=k;
					}//7
					else if(!volumePixels[i][j+1][k+1].isDone)
					{
						tp[0][0]=i;tp[0][1]=j;tp[0][2]=k+1;
						tp[1][0]=i+1;tp[1][1]=j+1;tp[1][2]=k+1;		
						tp[2][0]=i;tp[2][1]=j+1;tp[2][2]=k;
					}//8
					for(ii=0;ii<3;ii++)
					{
						if(vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]==-1)
						{
							vertSeq[tp[ii][0]][tp[ii][1]][tp[ii][2]]=vertexNumber;
							verts[vertexNumber].xyz(X)=tp[ii][0];
							verts[vertexNumber].xyz(Y)=tp[ii][1];
							verts[vertexNumber].xyz(Z)=tp[ii][2];
							vertexNumber++;
						}
					}
					faces[faceNumber].abc(A)=vertSeq[tp[0][0]][tp[0][1]][tp[0][2]];
					faces[faceNumber].abc(B)=vertSeq[tp[1][0]][tp[1][1]][tp[1][2]];
					faces[faceNumber++].abc(C)=vertSeq[tp[2][0]][tp[2][1]][tp[2][2]];
				}//total7
						
			}//every ijk
		}//j
	}//i
	verts=(vertInfo *)realloc(verts,vertexNumber*sizeof(vertInfo));
	faces=(faceInfo *)realloc(faces,faceNumber*sizeof(faceInfo));
	for(i=0;i<vertexNumber;i++)
	{
		verts[i].atomId=volumePixels[int(verts[i].xyz(X))][int(verts[i].xyz(Y))][int(verts[i].xyz(Z))].atomId;
		verts[i].isCont=false;
		if(volumePixels[int(verts[i].xyz(X))][int(verts[i].xyz(Y))][int(verts[i].xyz(Z))].isBound)
			verts[i].isCont=true;
	}
	for(i=0;i<pLength;i++)
	{
		for(j=0;j<pWidth;j++)
		{
			delete[]vertSeq[i][j];
		}
		delete[]vertSeq[i];
	}
	delete[]vertSeq;
}
//*/

