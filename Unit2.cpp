//---------------------------------------------------------------------------

#include <vcl.h>
#include <math.h>
#include <string>
#pragma hdrstop

#include "Unit1.h"
#include "Unit2.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
using namespace std;


extern Byte * MemA;

int ProcNum = 5; // количество методов обработки
string ProcName[8] = {"Исходное","А-Энтропия","А-Линейный","B-энтропия","B-линейный"};
int BarsNum[5]={0,0,1,0,1};

int postProcNum = 6;
string PostProcName[6] = {"Как есть","Гистограммный","Негатив","Сдвиг 192","Сдвиг 128","Сдвиг 64"};

int ParsPos[3][8]={{10,20,50,40,50,60,70,80},{10,20,30,40,50,60,70,80},{10,20,30,40,50,60,70,80}} ;

char* GetNames(int conv) {
	if(conv < ProcNum) {
		return (char*)ProcName[conv].c_str();
	} else {
		return NULL;
	}
}

char* GetPostNames(int postIndex) {
	if(postIndex < postProcNum) {
		return (char*)PostProcName[postIndex].c_str();
	} else {
		return NULL;
	}
}

int GetBarsNum(int conv) // после выбора конкретной функции обработки запрсить количество параметров, ненужные скрыть
{
	if(conv<ProcNum) return BarsNum[conv];
	else return 0;
}

int GetParsPos(int conv, int Bar) // для каждого активного параметра запросить позицию трэк-бара
{
	if(conv<ProcNum) {
	   if(Bar<BarsNum[conv]) return ParsPos[Bar][conv];

	}
	return 0;
}

void converter(int conv, Byte * Mem0, Byte * Mem1, int H, int W, int* PosK, int* PosM, int PostProcId, int Par1, int Par2, int Par3,int Chan) // Chan - выбор канала : 0 - все, 1 - R, 2- G, 3-B.
{

	switch (conv) {
	case 1:
		EntrA(Chan, Mem0,  MemA,  H,  W);
		break;
	case 2:
		LinA(Chan, Mem0,  MemA,  H,  W, Par1);
		break;
	case 3:
		EntrB(Chan, Mem0,  MemA,  H,  W);
		break;
	case 4:
		LinB(Chan, Mem0,  MemA,  H,  W);
		break;
	case 5:
		LinB(1, Mem0,  MemA,  H,  W);
		break;
	case 6:
		EntrB(2, Mem0,  MemA,  H,  W);
		break;
	case 7:
		LinB(2, Mem0,  MemA,  H,  W);
		break;
	default:
		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++)
			{
				int k = (i * W + j) * 3;
				if(Chan==0) {
				 MemA[k] = Mem0[k];
				 MemA[k + 1] = Mem0[k + 1];
				 MemA[k + 2] = Mem0[k + 2];
				}
				else {
				 MemA[k] = Mem0[k+Chan-1];
				 MemA[k + 1] = MemA[k];
				 MemA[k + 2] = MemA[k];
				}
			}
			MoveProgress(i, H);
		}
		break;
	}
	PostProcess(PostProcId, MemA, Mem1, H, W, Par1, Par2, Par3);
}

void PostProcess(int ProcID, Byte * MemIn, Byte * MemOut, int H, int W, int Par1, int Par2, int Par3)
{

  int gist[256];
  int outmap[256];
  int i,j,k,m,n;
	double Max;

	switch (ProcID) {
	  case 0:
		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++)
			{
				int k = (i * W + j) * 3;
				 MemOut[k] = MemIn[k];
				 MemOut[k + 1] = MemIn[k + 1];
				 MemOut[k + 2] = MemIn[k + 2];
			}
			MoveProgress(i, H);
		}
	  break;
	  case 1:
		n=0; for(i=0;i<256;i++) gist[i]=0;
		for (i = 0; i < H; i++) {
			for ( j = 0; j < W; j++)
			{
				int k = (i * W + j) * 3;
				 gist[MemIn[k]]++;
			}
			MoveProgress(i, H);
		}

		for(i=1;i<256;i++) gist[i]+=gist[i-1];
		Max=gist[255]/256;
		for(i=0;i<256;i++)  gist[i]/=Max;


		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++)
			{
				int k = (i * W + j) * 3;
				 MemOut[k] = outmap[MemIn[k]];
				 MemOut[k + 1] = outmap[MemIn[k + 1]];
				 MemOut[k + 2] = outmap[MemIn[k + 2]];
			}
			MoveProgress(i, H);
		}
	  break;
	  case 2:
		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++)
			{
				int k = (i * W + j) * 3;
				 MemOut[k] = 255-MemIn[k];
				 MemOut[k + 1] = 255-MemIn[k + 1];
				 MemOut[k + 2] = 255-MemIn[k + 2];
			}
			MoveProgress(i, H);
		}
	  break;
	  case 3:
		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++)
			{
				int k = (i * W + j) * 3;
				 MemOut[k] = (192+MemIn[k])%256;
				 MemOut[k + 1] = (192+MemIn[k+1])%256;
				 MemOut[k + 2] = (192+MemIn[k+2])%256;
			}
			MoveProgress(i, H);
		}
	  case 4:
		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++)
			{
				int k = (i * W + j) * 3;
				 MemOut[k] = (128+MemIn[k])%256;
				 MemOut[k + 1] = (128+MemIn[k+1])%256;
				 MemOut[k + 2] = (128+MemIn[k+2])%256;
			}
			MoveProgress(i, H);
		}
	  break;
	  case 5:
		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++)
			{
				int k = (i * W + j) * 3;
				 MemOut[k] = (64+MemIn[k])%256;
				 MemOut[k + 1] = (64+MemIn[k+1])%256;
				 MemOut[k + 2] = (64+MemIn[k+2])%256;
			}
			MoveProgress(i, H);
		}
	  break;

	}

}
   double egist[256][256];
   double Hgist[256][256];
   int gist1[256][256];
   int outmap1[256][256];

void EntrA(int chan, Byte * Mem0, Byte * Mem1, int H, int W)
{

   int gist[256];
   int outmap[256];

   int i,j,k;
   int i1,j1,k1,m,M,N;
   int M0,MC,M1,M2,M3;
   int Lev;

   int gist2[256];
   int gist3[256];

   int outmap2[256];
   int outmap3[256];

		for(i=0;i<256;i++) gist2[i]=gist3[i]=0; N=0;
		for(i=0;i<256;i++) for(j=0;j<256;j++) gist1[i][j]=0;

		 for ( i = 1; i < H-1; i++) {
			for ( j = 1; j < W-1; j++)
			{
				 k = (i * W + j) * 3;
				//int wb = (Mem0[k] + Mem0[k + 1] + Mem0[k + 2]) / 3;
				m=1; M1=M2=M3=0;
				for(i1=i-1;i1<i+2;i1++)
				for(j1=j-1;j1<j+2;j1++)
				{
					if((i1==i)&&(j1==j)) continue;
					k1 = (i1 * W + j1) * 3;
					M1+=Mem0[k1];
					M2+=Mem0[k1+1];
					M3+=Mem0[k1+2];
				}

			   M1/=8; M2/=8; M3/=8;
			   switch(chan){
			   case 0:
				  MC=(M1+M2+M3)/3;
				  M0=(Mem0[k]+Mem0[k]+Mem0[k])/3;
			   break;
			   case 1:
				  MC=M1;
				  M0=Mem0[k];
			   break;
			   case 2:
				  MC=M2;
				  M0=Mem0[k+1];
			   break;
			   case 3:
				  MC=M3;
				  M0=Mem0[k+2];
			   break;
			   }

				gist1[MC][M0]++;

				N++;

			}
			MoveProgress(i, H);
		}
		for(i=0;i<256;i++) for(j=0;j<256;j++) egist[i][j]=gist1[i][j]/(double)N;
		for(i=0;i<256;i++) for(j=0;j<256;j++)
		   if(egist[i][j]>0) Hgist[i][j]=-1*log(egist[i][j]);
		   else Hgist[i][j]=0;
		 //  double eh=0; for(i=0;i<256;i++) eh += Hgist[i];
		//   eh/=256;
		//   for(i=0;i<256;i++)  Hgist[i]-=eh;
		for(i=0;i<256;i++) for(j=0;j<256;j++) {
		   outmap1[i][j]= 20*Hgist[i][j]; if(outmap1[i][j]<0) outmap1[i][j]=0; else if(outmap1[i][j]>255) outmap1[i][j]=255;
		}


		for ( i = 1; i < H-1; i++) {
			for ( j = 1; j < W-1; j++)
			{
				 k = (i * W + j) * 3;
				//int wb = (Mem0[k] + Mem0[k + 1] + Mem0[k + 2]) / 3;
				m=1; M1=M2=M3=0;
				for(i1=i-1;i1<i+2;i1++)
				for(j1=j-1;j1<j+2;j1++)
				{
					if((i1==i)&&(j1==j)) continue;
					k1 = (i1 * W + j1) * 3;
					M1+=Mem0[k1];
					M2+=Mem0[k1+1];
					M3+=Mem0[k1+2];
				}
			   M1/=8; M2/=8; M3/=8;
			   switch(chan){
			   case 0:
				  MC=(M1+M2+M3)/3;
				  M0=(Mem0[k]+Mem0[k]+Mem0[k])/3;
			   break;
			   case 1:
				  MC=M1;
				  M0=Mem0[k];
			   break;
			   case 2:
				  MC=M2;
				  M0=Mem0[k+1];
			   break;
			   case 3:
				  MC=M3;
				  M0=Mem0[k+2];
			   break;
			   }



				Lev = outmap1[MC][M0];
				Mem1[k] = Lev;
				Mem1[k+1] = Lev;
				Mem1[k+2] = Lev;


			}
			MoveProgress(i, H);
		}
}
void LinA(int chan, Byte * Mem0, Byte * Mem1, int H, int W, int Par1)
{

   int gist[256];
   int outmap[256];

   int i,j,k;
   int i1,j1,k1,m,M,N;
   int M0,MC,M1,M2,M3;
   int Lev;

   int gist2[256];
   int gist3[256];

   int outmap2[256];
   int outmap3[256];

		for(i=0;i<256;i++) gist2[i]=gist3[i]=0; N=0;
		for(i=0;i<256;i++) for(j=0;j<256;j++) gist1[i][j]=0;

		 for ( i = 1; i < H-1; i++) {
			for ( j = 1; j < W-1; j++)
			{
				 k = (i * W + j) * 3;
				//int wb = (Mem0[k] + Mem0[k + 1] + Mem0[k + 2]) / 3;
				m=1; M1=M2=M3=0;
				for(i1=i-1;i1<i+2;i1++)
				for(j1=j-1;j1<j+2;j1++)
				{
					if((i1==i)&&(j1==j)) continue;
					k1 = (i1 * W + j1) * 3;
					M1+=Mem0[k1];
					M2+=Mem0[k1+1];
					M3+=Mem0[k1+2];
				}

			   M1/=8; M2/=8; M3/=8;
			   switch(chan){
			   case 0:
				  MC=(M1+M2+M3)/3;
				  M0=(Mem0[k]+Mem0[k]+Mem0[k])/3;
			   break;
			   case 1:
				  MC=M1;
				  M0=Mem0[k];
			   break;
			   case 2:
				  MC=M2;
				  M0=Mem0[k+1];
			   break;
			   case 3:
				  MC=M3;
				  M0=Mem0[k+2];
			   break;
			   }

				gist1[MC][M0]++;

				N++;

			}
			MoveProgress(i, H);
		}
/*
		for(i=0;i<256;i++) for(j=0;j<256;j++) egist[i][j]=gist1[i][j]/(double)N;
		for(i=0;i<256;i++) for(j=0;j<256;j++)
		   if(egist[i][j]>0) Hgist[i][j]=32*egist[i][j];
		   else Hgist[i][j]=0;

		for(i=0;i<256;i++) for(j=0;j<256;j++) {
		   outmap1[i][j]= Hgist[i][j]; if(outmap1[i][j]<0) outmap1[i][j]=0; else if(outmap1[i][j]>255) outmap1[i][j]=255;
		}
*/
		double A= Par1*50;
		for(i=0;i<256;i++) for(j=0;j<256;j++) egist[i][j]=gist1[i][j]/(double)N;
		for(i=0;i<256;i++) for(j=0;j<256;j++)
		   if(egist[i][j]>0) Hgist[i][j]= A*egist[i][j]; // ;-1*log(egist[i][j]);
		   else Hgist[i][j]=0;


		for(i=0;i<256;i++) for(j=0;j<256;j++) {
		   outmap1[i][j]= 20*Hgist[i][j]; if(outmap1[i][j]<0) outmap1[i][j]=0; else if(outmap1[i][j]>255) outmap1[i][j]=255;
		}



		for ( i = 1; i < H-1; i++) {
			for ( j = 1; j < W-1; j++)
			{
				 k = (i * W + j) * 3;
				//int wb = (Mem0[k] + Mem0[k + 1] + Mem0[k + 2]) / 3;
				m=1; M1=M2=M3=0;
				for(i1=i-1;i1<i+2;i1++)
				for(j1=j-1;j1<j+2;j1++)
				{
					if((i1==i)&&(j1==j)) continue;
					k1 = (i1 * W + j1) * 3;
					M1+=Mem0[k1];
					M2+=Mem0[k1+1];
					M3+=Mem0[k1+2];
				}
			   M1/=8; M2/=8; M3/=8;
			   switch(chan){
			   case 0:
				  MC=(M1+M2+M3)/3;
				  M0=(Mem0[k]+Mem0[k]+Mem0[k])/3;
			   break;
			   case 1:
				  MC=M1;
				  M0=Mem0[k];
			   break;
			   case 2:
				  MC=M2;
				  M0=Mem0[k+1];
			   break;
			   case 3:
				  MC=M3;
				  M0=Mem0[k+2];
			   break;
			   }



				Lev = outmap1[MC][M0];
				Mem1[k] = Lev;
				Mem1[k+1] = Lev;
				Mem1[k+2] = Lev;


			}
			MoveProgress(i, H);
		}
}
void EntrB(int chan, Byte * Mem0, Byte * Mem1, int H, int W)
{

   int gist[256];
   int outmap[256];
   double egist[256];
   double Hgist[256];
   int i,j,k;
   int i1,j1,k1,m,M,N;
   int M1,M2,M3;
   int Lev;
   int gist1[256];
   int gist2[256];
   int gist3[256];
   int outmap1[256];
   int outmap2[256];
   int outmap3[256];

		for(i=0;i<256;i++) gist1[i]=gist2[i]=gist3[i]=0; N=0;

		 for ( i = 1; i < H-1; i++) {
			for ( j = 1; j < W-1; j++)
			{
				 k = (i * W + j) * 3;
				//int wb = (Mem0[k] + Mem0[k + 1] + Mem0[k + 2]) / 3;
				m=1; M1=M2=M3=0;
				for(i1=i-1;i1<i+2;i1++)
				for(j1=j-1;j1<j+2;j1++)
				{
					if((i1==i)&&(j1==j)) continue;
					k1 = (i1 * W + j1) * 3;
					if(chan==0) {
					  if((Mem0[k]+Mem0[k+1]+Mem0[k+2])>(Mem0[k1]+Mem0[k1+1]+Mem0[k1+2])) M1+=m;
					}
					else if(Mem0[k+chan-1]>Mem0[k1+chan-1]) M1+=m;
					m<<=1;
				}
				gist1[M1]++;

				N++;

			}
			MoveProgress(i, H);
		}
		for(i=0;i<256;i++) egist[i]=gist1[i]/(double)N;
		for(i=0;i<256;i++)
		   if(egist[i]>0) Hgist[i]=-1*log(egist[i]);
		   else Hgist[i]=0;
		   double eh=0; for(i=0;i<256;i++) eh += Hgist[i];
		//   eh/=256;
		//   for(i=0;i<256;i++)  Hgist[i]-=eh;
		for(i=0;i<256;i++) {
		   outmap1[i]= 40*Hgist[i]; if(outmap1[i]<0) outmap1[i]=0; else if(outmap1[i]>255) outmap1[i]=255;
		}


		for ( i = 1; i < H-1; i++) {
			for ( j = 1; j < W-1; j++)
			{
				 k = (i * W + j) * 3;
				//int wb = (Mem0[k] + Mem0[k + 1] + Mem0[k + 2]) / 3;
				m=1; M1=M2=M3=0;
				for(i1=i-1;i1<i+2;i1++)
				for(j1=j-1;j1<j+2;j1++)
				{
					if((i1==i)&&(j1==j)) continue;
					k1 = (i1 * W + j1) * 3;
					if(chan==0) {
					  if((Mem0[k]+Mem0[k+1]+Mem0[k+2])>(Mem0[k1]+Mem0[k1+1]+Mem0[k1+2])) M1+=m;
					}
					else if(Mem0[k+chan-1]>Mem0[k1+chan-1]) M1+=m;
					m<<=1;
				}
				Lev = outmap1[M1];
				Mem1[k] = Lev;
				Mem1[k+1] = Lev;
				Mem1[k+2] = Lev;


			}
			MoveProgress(i, H);
		}
}
void LinB(int chan,  Byte * Mem0, Byte * Mem1, int H, int W)
{

   int gist[256];
   int outmap[256];
   double egist[256];
   double Hgist[256];
   int i,j,k;
   int i1,j1,k1,m,M,N;
   int M1,M2,M3;
   int Lev;
   int gist1[256];
   int gist2[256];
   int gist3[256];
   int outmap1[256];
   int outmap2[256];
   int outmap3[256];

		for(i=0;i<256;i++) gist1[i]=gist2[i]=gist3[i]=0; N=0;

		 for ( i = 1; i < H-1; i++) {
			for ( j = 1; j < W-1; j++)
			{
				 k = (i * W + j) * 3;
				//int wb = (Mem0[k] + Mem0[k + 1] + Mem0[k + 2]) / 3;
				m=1; M1=M2=M3=0;
				for(i1=i-1;i1<i+2;i1++)
				for(j1=j-1;j1<j+2;j1++)
				{
					if((i1==i)&&(j1==j)) continue;
					k1 = (i1 * W + j1) * 3;
					if(chan==0) {
					  if((Mem0[k]+Mem0[k+1]+Mem0[k+2])>(Mem0[k1]+Mem0[k1+1]+Mem0[k1+2])) M1+=m;
					}
					else if(Mem0[k+chan-1]>Mem0[k1+chan-1]) M1+=m;
					m<<=1;
				}
				gist1[M1]++;
				N++;

			}
			MoveProgress(i, H);
		}
		for(i=0;i<256;i++) {
		  egist[i]=gist1[i]/(double)N; egist[i]*=256;  //egist[i]=1-egist[i];
		}
		for(i=0;i<256;i++)
		   if(egist[i]>0) Hgist[i]=32*egist[i];
		   else Hgist[i]=0;

		//   for(i=0;i<256;i++)  Hgist[i]-=eh;
		for(i=0;i<256;i++) {
		   outmap1[i]= Hgist[i]; if(outmap1[i]<0) outmap1[i]=0; else if(outmap1[i]>255) outmap1[i]=255;
		}


		for ( i = 1; i < H-1; i++) {
			for ( j = 1; j < W-1; j++)
			{
				 k = (i * W + j) * 3;
				//int wb = (Mem0[k] + Mem0[k + 1] + Mem0[k + 2]) / 3;
				m=1; M1=M2=M3=0;
				for(i1=i-1;i1<i+2;i1++)
				for(j1=j-1;j1<j+2;j1++)
				{
					if((i1==i)&&(j1==j)) continue;
					k1 = (i1 * W + j1) * 3;
					if(chan==0) {
					  if((Mem0[k]+Mem0[k+1]+Mem0[k+2])>(Mem0[k1]+Mem0[k1+1]+Mem0[k1+2])) M1+=m;
					}
					else if(Mem0[k+chan-1]>Mem0[k1+chan-1]) M1+=m;
					m<<=1;
				}
				Lev = outmap1[M1];
				Mem1[k] = Lev;
				Mem1[k+1] = Lev;
				Mem1[k+2] = Lev;


			}
			MoveProgress(i, H);
		}
}

void Entr8EQ(int chan, Byte * Mem0, Byte * Mem1, int H, int W)
{

   int gist[256];
   int outmap[256];
   double egist[256];
   double Hgist[256];
   int i,j,k;
   int i1,j1,k1,m,M,N;
   int M1,M2,M3;
   int Lev;
   int gist1[256];
   int gist2[256];
   int gist3[256];
   int outmap1[256];
   int outmap2[256];
   int outmap3[256];

		for(i=0;i<256;i++) gist1[i]=gist2[i]=gist3[i]=0; N=0;

		 for ( i = 1; i < H-1; i++) {
			for ( j = 1; j < W-1; j++)
			{
				 k = (i * W + j) * 3;
				//int wb = (Mem0[k] + Mem0[k + 1] + Mem0[k + 2]) / 3;
				m=1; M1=M2=M3=0;
				for(i1=i-1;i1<i+2;i1++)
				for(j1=j-1;j1<j+2;j1++)
				{
					if((i1==i)&&(j1==j)) continue;
					k1 = (i1 * W + j1) * 3;
					if(abs(Mem0[k+chan]-Mem0[k1+chan])<1) M1+=m;
					m<<=1;
				}
				gist1[M1]++;

				N++;

			}
			MoveProgress(i, H);
		}
		for(i=0;i<256;i++) egist[i]=gist1[i]/(double)N;
		for(i=0;i<256;i++)
		   if(egist[i]>0) Hgist[i]=-1*log(egist[i]);
		   else Hgist[i]=0;
		   double eh=0; for(i=0;i<256;i++) eh += Hgist[i];
		//   eh/=256;
		//   for(i=0;i<256;i++)  Hgist[i]-=eh;
		for(i=0;i<256;i++) {
		   outmap1[i]= 40*Hgist[i]; if(outmap1[i]<0) outmap1[i]=0; else if(outmap1[i]>255) outmap1[i]=255;
		}


		for ( i = 1; i < H-1; i++) {
			for ( j = 1; j < W-1; j++)
			{
				 k = (i * W + j) * 3;
				//int wb = (Mem0[k] + Mem0[k + 1] + Mem0[k + 2]) / 3;
				m=1; M1=M2=M3=0;
				for(i1=i-1;i1<i+2;i1++)
				for(j1=j-1;j1<j+2;j1++)
				{
					if((i1==i)&&(j1==j)) continue;
					k1 = (i1 * W + j1) * 3;
					if(Mem0[k+chan]>Mem0[k1+chan]) M1+=m;
					m<<=1;
				}
				Lev = outmap1[M1];
				Mem1[k] = Lev;
				Mem1[k+1] = Lev;
				Mem1[k+2] = Lev;


			}
			MoveProgress(i, H);
		}
}

void TestFilter( Byte * Mem0, Byte * Mem1, int H, int W)
{

   int gist[256];
   int outmap[256];
   double egist[256];
   double Hgist[256];
   int i,j,k;
   int i1,j1,k1,m,M,N;
   int M1,M2,M3;
   int Lev;
   int gist1[256];
   int gist2[256];
   int gist3[256];
   int outmap1[256];
   int outmap2[256];
   int outmap3[256];

		for(i=0;i<256;i++) gist1[i]=gist2[i]=gist3[i]=0;

		 for ( i = 1; i < H-1; i++) {
			for ( j = 1; j < W-1; j++)
			{
				 k = (i * W + j) * 3;
				//int wb = (Mem0[k] + Mem0[k + 1] + Mem0[k + 2]) / 3;
				m=1; M1=M2=M3=0;
				for(i1=i-1;i1<i+2;i1++)
				for(j1=j-1;j1<j+2;j1++)
				{
					if((i1==i)&&(j1==j)) continue;
					k1 = (i1 * W + j1) * 3;
					if(Mem0[k]>Mem0[k1]) M1+=m;
					if(Mem0[k+1]>Mem0[k1+1]) M2+=m;
					if(Mem0[k+2]>Mem0[k1+2]) M3+=m;
					m<<=1;
				}
				gist1[M1]++;
				gist2[M2]++;
				gist3[M3]++;
				N++;

			}
			MoveProgress(i, H);
		}
		for(i=0;i<256;i++) egist[i]=gist1[i]/(double)N;
		for(i=0;i<256;i++)
		   if(egist[i]>0) Hgist[i]=-1*log(egist[i]);
		   else Hgist[i]=0;
		for(i=0;i<256;i++) outmap1[i]= 30*Hgist[i];

		for(i=0;i<256;i++) egist[i]=gist2[i]/(double)N;
		for(i=0;i<256;i++)
		   if(egist[i]>0) Hgist[i]=-1*log(egist[i]);
		   else Hgist[i]=0;
		for(i=0;i<256;i++) outmap2[i]= 30*Hgist[i];

		for(i=0;i<256;i++) egist[i]=gist3[i]/(double)N;
		for(i=0;i<256;i++)
		   if(egist[i]>0) Hgist[i]=-1*log(egist[i]);
		   else Hgist[i]=0;
		for(i=0;i<256;i++) outmap3[i]= 30*Hgist[i];

		for ( i = 1; i < H-1; i++) {
			for ( j = 1; j < W-1; j++)
			{
				 k = (i * W + j) * 3;
				//int wb = (Mem0[k] + Mem0[k + 1] + Mem0[k + 2]) / 3;
				m=1; M1=M2=M3=0;
				for(i1=i-1;i1<i+2;i1++)
				for(j1=j-1;j1<j+2;j1++)
				{
					if((i1==i)&&(j1==j)) continue;
					k1 = (i1 * W + j1) * 3;
					if(Mem0[k]>Mem0[k1]) M1+=m;
					if(Mem0[k+1]>Mem0[k1+1]) M2+=m;
					if(Mem0[k+2]>Mem0[k1+2]) M3+=m;
					m<<=1;
				}
				Lev = outmap1[M1]; if(Lev>255) Lev=255; else if(Lev<0) Lev=0; Mem1[k] = Lev;
				Lev = outmap2[M2]; if(Lev>255) Lev=255; else if(Lev<0) Lev=0; Mem1[k+1] = Lev;
				Lev = outmap3[M3]; if(Lev>255) Lev=255; else if(Lev<0) Lev=0; Mem1[k+2] = Lev;


			}
			MoveProgress(i, H);
		}
}

void TestFilter2( Byte * Mem0, Byte * Mem1, int H, int W)
{
		 for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++)
			{
				int k = (i * W + j) * 3;
				int wb = (Mem0[k] + Mem0[k + 1] + Mem0[k + 2]) / 3;
				Mem0[k] = wb;
				Mem0[k + 1] = Mem0[k];
				Mem0[k + 2] = Mem0[k];
				Mem1[k] = Mem0[k];
				Mem1[k+1] = Mem0[k+1];
				Mem1[k+2] = Mem0[k+2];
			}
			MoveProgress(i, H);
		}
}
