// ---------------------------------------------------------------------------

#include <vcl.h>
#include <stdio.h>
#include <math.h>
#pragma hdrstop

#include "Unit1.h"
// ---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;

// Размеры окна для отображения на форме
#define H 780
#define W 1400

// Максимальные размеры входного изображения
#define maxH 10000
#define maxW 10000

int inH, inW; // Реальные размеры входного изображения

struct rgb // структура для доступа к bmp 24 бит
{
	Byte b;
	Byte g;
	Byte r;
	// Byte a;
};

struct rgb32 // // структура для доступа к bmp 32 бит (не используется)
{
	Byte b;
	Byte g;
	Byte r;
	Byte a;
};

int Mem2[H][W][3]; // массив пикселей для вывода на форму
int Mem2Count[H][W]; // счетчик для усреднения пикселей

Byte * Mem1; // указатель на буфер под пиксели из входного файла

Graphics::TBitmap *inbmp1;
Graphics::TBitmap *outbmp;

AnsiString fname;

// ---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner) : TForm(Owner) {
	inbmp1 = new Graphics::TBitmap;
	outbmp = new Graphics::TBitmap;
	outbmp->Height = H;
	outbmp->Width = W;
	outbmp->PixelFormat = pf24bit;

	long size = maxW;
	size *= maxH;
	size *= 3;
	Mem1 = new Byte[size]; // буфер под пиксели из входного файла

}
//----------------------------------------------------------------------------
void BeeLineInterpolation(int x, int y) {
	if (x==15 && y==8) {
		int t=1;
	}
	int x1, y1, x2, y2;
	int i=x-1;
	bool stop = false;
	while(!stop){
		if(Mem2Count[y-1][i]>0)
			stop=true;
		else i--;
	}
	x1=i; y2=y-1;
	i=x+1; stop = false;
	while(!stop){
		if(Mem2Count[y-1][i]>0)
			stop=true;
		else i++;
	}
	x2=i;
	int j=y+1;	stop = false;
	while(!stop){
		if(Mem2Count[j][i]>0)
			stop=true;
		else j++;
	}
	y1=j;
	float dx1, dx2, dx3, dy1, dy2, dy3;
	dx1=x2-x1; dx2=x2-x; dx3=x-x1;
	dy1=y2-y1; dy2=y2-y; dy3=y-y1;
	Mem2[y][x][0]=(1/(dx1*dy1))*(Mem2[y1][x1][0]*dx2*dy2 + Mem2[y1][x2][0]*dx3*dy2 + Mem2[y2][x1][0]*dx2*dy3 + Mem2[y2][x2][0]*dx3*dy3);
	Mem2[y][x][1]=(1/(dx1*dy1))*(Mem2[y1][x1][1]*dx2*dy2 + Mem2[y1][x2][1]*dx3*dy2 + Mem2[y2][x1][1]*dx2*dy3 + Mem2[y2][x2][1]*dx3*dy3);
	Mem2[y][x][2]=(1/(dx1*dy1))*(Mem2[y1][x1][2]*dx2*dy2 + Mem2[y1][x2][2]*dx3*dy2 + Mem2[y2][x1][2]*dx2*dy3 + Mem2[y2][x2][2]*dx3*dy3);
	Mem2Count[y][x]++;
}
//----------------------------------------------------------------------------------
void LineInterpolationX(int x, int y) {
	int x0, y0, x1, y1;
	int i=x-1;
	bool stop = false;
	while(!stop){
		if(Mem2Count[y][i]>0)
			stop=true;
		else i--;
	}
	x0=i; y0=y;
	i=x+1; stop = false;
	while(!stop){
		if(Mem2Count[y][i]>0)
			stop=true;
		else i++;
	}
	x1=i; y1=y;
	float dx1, dx2, dx3;
	dx1=x1-x; dx2=x1-x0; dx3=x-x0;
	Mem2[y][x][0]=Mem2[y][x0][0]*(dx1/dx2) + Mem2[y][x1][0]*(dx3/dx2);
	Mem2[y][x][1]=Mem2[y][x0][1]*(dx1/dx2) + Mem2[y][x1][1]*(dx3/dx2);
	Mem2[y][x][2]=Mem2[y][x0][2]*(dx1/dx2) + Mem2[y][x1][2]*(dx3/dx2);
	Mem2Count[y][x]++;
}
// ---------------------------------------------------------------------------
void LineInterpolationY(int x, int y) {
	int y0, y1;
	int j=y+1;
	bool stop = false;
	while(!stop){
		if(Mem2Count[j][x]>0)
			stop=true;
		else j++;
	}
	y0=j;
	j=y-1; stop = false;
	while(!stop){
		if(Mem2Count[j][x]>0)
			stop=true;
		else j--;
	}
	y1=j;
	float dy1, dy2, dy3;
	dy1=y1-y; dy2=y1-y0; dy3=y-y0;
	int a =  Mem2[y0][x][0]*(dy1/dy2);
	int b =  Mem2[y1][x][0]*(dy3/dy2);
    int c = a+b;
	Mem2[y][x][0]=Mem2[y0][x][0]*(dy1/dy2) + Mem2[y1][x][0]*(dy3/dy2);
	Mem2[y][x][1]=Mem2[y0][x][1]*(dy1/dy2) + Mem2[y1][x][1]*(dy3/dy2);
	Mem2[y][x][2]=Mem2[y0][x][2]*(dy1/dy2) + Mem2[y1][x][2]*(dy3/dy2);
	Mem2Count[y][x]++;
}
// ---------------------------------------------------------------------------
void PutMem1toMem2() // прореживание отсчетов с усреднением
{
	int i, j, k, i1, j1;
	for (i = 0; i < H; i++)
		for (j = 0; j < W; j++)
			for (k = 0; k < 3; k++) {
				Mem2[i][j][k] = 0;
				Mem2Count[i][j] = 0; // обнуление памяти и счетика
			}

	int BSize = inH * inW * 3;

	/*
	inW -> W   k=inW/W          неверно
	inW -> W   inW=(k*(W-1)+1)  верно
	inW-1 = k*(W-1)
	k = (inW-1)/(W-1)
	*/

//	double kX = (inW) / (double)(W); // коэффициент масштабирования по X
//	double kY = (inH) / (double)(H); // коэффициент масштабирования по Y
	double kX = (inW-1) / (double)(W-1); // коэффициент масштабирования по X
	double kY = (inH-1) / (double)(H-1); // коэффициент масштабирования по Y

	for (int i = 0; i < inH; i++) {
		for (int j = 0; j < inW; j++) {
			int k = (i * inW + j) * 3;
			// пересчет трехмерных координат (Y,X,цветовой канал) в одномерные
			i1 = i / kY;
			j1 = j / kX; // позиция пикселя в выходном изображении
			Mem2[i1][j1][0] += Mem1[k]; // красный
			Mem2[i1][j1][1] += Mem1[k + 1]; // синий
			Mem2[i1][j1][2] += Mem1[k + 2]; // зеленый
			Mem2Count[i1][j1]++;
		}

	}
	if (inH > H || inW > W)
		// усреднение получившихся пикселей
		for (i = 0; i < H; i++)
			for (j = 0; j < W; j++)
				for (k = 0; k < 3; k++) {
					if (Mem2Count[i][j] > 1)
						Mem2[i][j][k] /= Mem2Count[i][j];
				}
	int pW = (inW-1)/kX;
	int pH = (inH-1)/kY;
	if(inH < H || inW < W){
		//интерполяция
		for (int j = 0; j <= pW; j++) {
			if (Mem2Count[0][j] == 0) LineInterpolationX(j, 0);
			if (Mem2Count[pH][j] == 0) LineInterpolationX(j, pH);
		}
		for (int i = 0; i <= pH; i++) {
			if (Mem2Count[i][0] == 0) LineInterpolationY(0, i);
			if (Mem2Count[i][pW] == 0) LineInterpolationY(pW, i);
		}
		for (int i = 1; i < pH; i++)
			for (int j = 1; j < pW; j++)
				if (Mem2Count[i][j] ==0)  BeeLineInterpolation(j, i);
				//if (Mem2Count[i][j] ==0)  LineInterpolationX(j, i);

	}

}

void __fastcall TForm1::OpenImage1Click(TObject *Sender) {
	if (OpenPictureDialog1->Execute()) {
		fname = OpenPictureDialog1->FileName;
		FILE *f = fopen(fname.c_str(), "r");
		if (f) {
			fclose(f);
			inbmp1->LoadFromFile(fname); // загрузка изображения из файла

			inH = inbmp1->Height;
			inW = inbmp1->Width;
			for (int i = 0; i < inH; i++) {
				struct rgb *ptr = (struct rgb*)inbmp1->ScanLine[i];

				for (int j = 0; j < inW; j++) // запись пикселей в память для дальнейшей обработки
				{
					int k = (i * inW + j) * 3;
					Mem1[k] = ptr[j].r;
					Mem1[k + 1] = ptr[j].g;
					Mem1[k + 2] = ptr[j].b;
				}
			}

			PutMem1toMem2();
			// перенос пикселей из входного изображения в выходное
			for (int i = 0; i < H-1; i++) // подготовка изображения к выводу на форму
			{
				struct rgb *ptr = (struct rgb*)outbmp->ScanLine[i];
				for (int j = 0; j < W-1; j++) {
					ptr[j].r = Mem2[i][j][0];
					ptr[j].g = Mem2[i][j][1];
					ptr[j].b = Mem2[i][j][2];
					//if(ptr[j].r ==0 && ptr[j].g==0 && ptr[j].b==0) {
					if (Mem2Count[i][j] ==0) {
						;
					   ShowMessage(AnsiString()+ptr[j].r + ptr[j].g + ptr[j].b +" i:" + i + " j:" + j);
					}//}

				}
			}
			// Вывод изображения на форму
			Image1->Picture->Graphic = outbmp;

			WWin->Text = inW;
			HWin->Text = inH;

			// Расчет показателей эффективной разрядности
			int BSize = inH * inW * 3;
			int Gist[256];
			int i;
			for (i = 0; i < 256; i++)
				Gist[i] = 0;
			for (i = 0; i < BSize; i++)
				Gist[Mem1[i]]++;
			double h1 = 0;
			for (i = 0; i < 256; i++) {
				if (Gist[i] == 0)
					continue;
				double N = BSize;
				double p = Gist[i] / N;
				h1 -= p * log(p) / log(2);
			}
			double n = 0;
			for (i = 0; i < 256; i++)
				if (Gist[i] > 0) {
					Gist[i] = 1;
					n += 1;
				}

			double h2 = 0;
			for (i = 0; i < 256; i++) {
				if (Gist[i] == 0)
					continue;

				double p = Gist[i] / n;
				h2 -= p * log(p) / log(2);
			}
			char buf[64];
			sprintf(buf, "%.2f (%.2f)", h2, h1);
			NBitWin->Text = buf;

		}
	}
}
