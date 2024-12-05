// ---------------------------------------------------------------------------

#include <vcl.h>
//#include <stdio.h>
#include <math.h>
#pragma hdrstop

#include "Unit1.h"
#include "Unit2.h"
// ---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;

// Максимальные размеры окна для отображения на форме
#define WH 2160
#define WW 3840

int H, W;     // Размеры окна для отображения на форме
int inH, inW; // Реальные размеры входного изображения

struct rgb // структура для доступа к bmp 24 бит
{
	Byte b;
	Byte g;
	Byte r;
	// Byte a;
};

//struct rgb32 // // структура для доступа к bmp 32 бит (не используется)
//{
//	Byte b;
//	Byte g;
//	Byte r;
//	Byte a;
//};
int SW = Screen->Width;
int SH = Screen->Height;


int Mem2[WH][WW][3]; // массив пикселей для вывода на форму
int Mem2Count[WH][WW]; // счетчик для усреднения пикселей

Byte * Mem0; // указатель на первичный буфер под пиксели из входного файла
Byte * Mem1; // указатель на вторичный буфер под пиксели из входного файла (после дополнительного преобразования, если оно есть)

Graphics::TBitmap *inbmp1;
Graphics::TBitmap *outbmp;

AnsiString fname;

// ---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner) : TForm(Owner) {
	Form1->Width = 0.85*SW;
	Form1->Height = 0.8*SH;
	Form1->Top = 0.1*SH;
	Form1->Left = 0.05*SW;
	Image1->Width = 0.86*0.85*SW;
	Image1->Height = 0.96*0.8*SH;
	GroupBox1->Left = 0.87*0.85*SW;
	GroupBox2->Left = 0.87*0.85*SW;
	OpenImage1->Left = 0.88*0.85*SW;

	H = Image1->Height;
	W = Image1->Width;
	inbmp1 = new Graphics::TBitmap;
	outbmp = new Graphics::TBitmap;
	outbmp->Height = H;
	outbmp->Width = W;
	outbmp->PixelFormat = pf24bit;
}
//----------------------------------------------------------------------------
void BeeLineInterpolation(int x, int y) {
	//поиск четырех известных соседей
	int x1, y1, x2, y2, i1, j1;
	int i = x - 1;
	bool stop = false;
	while (!stop) {
		if (Mem2Count[y - 1][i] > 0)
			stop = true;
		else
			i--;
	}
	x1 = i;
	y2 = y - 1;
	i = x + 1;
	stop = false;
	while (!stop) {
		if (Mem2Count[y - 1][i] > 0)
			stop = true;
		else
			i++;
	}
	x2 = i;
	int j = y + 1;
	stop = false;
	while (!stop) {
		if (Mem2Count[j][i] > 0)
			stop = true;
		else
			j++;
	}
	y1 = j;
	if (Mem2Count[y1][x1] < 1) {
		int i = x1;
		stop = false;
		while (!stop) {
			if (Mem2Count[y1][i] > 0 && Mem2Count[y2][i] > 0)
				stop = true;
			else
				i--;
		}
		x1 = i;
	}
	//рассчет по полученным соседним точкам
	float dx1, dx2, dx3, dy1, dy2, dy3;
	dx1 = x2 - x1;
	dx2 = x2 - x;
	dx3 = x - x1;
	dy1 = y2 - y1;
	dy2 = y2 - y;
	dy3 = y - y1;
	Mem2[y][x][0] = (1.0 / (dx1 * dy1)) * (Mem2[y1][x1][0] * dx2 * dy2 +
		Mem2[y1][x2][0] * dx3 * dy2 + Mem2[y2][x1][0] * dx2 * dy3 +
		Mem2[y2][x2][0] * dx3 * dy3);
	Mem2[y][x][1] = (1.0 / (dx1 * dy1)) * (Mem2[y1][x1][1] * dx2 * dy2 +
		Mem2[y1][x2][1] * dx3 * dy2 + Mem2[y2][x1][1] * dx2 * dy3 +
		Mem2[y2][x2][1] * dx3 * dy3);
	Mem2[y][x][2] = (1.0 / (dx1 * dy1)) * (Mem2[y1][x1][2] * dx2 * dy2 +
		Mem2[y1][x2][2] * dx3 * dy2 + Mem2[y2][x1][2] * dx2 * dy3 +
		Mem2[y2][x2][2] * dx3 * dy3);

	Mem2Count[y][x]++;
}

// ----------------------------------------------------------------------------------
void LineInterpolationX(int x, int y) {
	//поиск известных соседей
	int x0, y0, x1, y1;
	int i = x - 1;
	bool stop = false;
	while (!stop) {
		if (Mem2Count[y][i] > 0)
			stop = true;
		else
			i--;
	}
	x0 = i;
	y0 = y;
	i = x + 1;
	stop = false;
	while (!stop) {
		if (Mem2Count[y][i] > 0)
			stop = true;
		else
			i++;
	}
	x1 = i;
	y1 = y;
	//рассчет по полученным соседним точкам
	float dx1, dx2, dx3;
	dx1 = x1 - x;
	dx2 = x1 - x0;
	dx3 = x - x0;
	Mem2[y][x][0] = Mem2[y][x0][0] * (dx1 / dx2) + Mem2[y][x1][0] * (dx3 / dx2);
	Mem2[y][x][1] = Mem2[y][x0][1] * (dx1 / dx2) + Mem2[y][x1][1] * (dx3 / dx2);
	Mem2[y][x][2] = Mem2[y][x0][2] * (dx1 / dx2) + Mem2[y][x1][2] * (dx3 / dx2);
	Mem2Count[y][x]++;
}


// ---------------------------------------------------------------------------
void LineInterpolationY(int x, int y) {
	//поиск известных соседей
	int y0, y1;
	int j = y + 1;
	bool stop = false;
	while (!stop) {
		if (Mem2Count[j][x] > 0)
			stop = true;
		else
			j++;
	}
	y0 = j;
	j = y - 1;
	stop = false;
	while (!stop) {
		if (Mem2Count[j][x] > 0)
			stop = true;
		else
			j--;
	}
	y1 = j;
    //рассчет по полученным соседним точкам
	float dy1, dy2, dy3;
	dy1 = y1 - y;
	dy2 = y1 - y0;
	dy3 = y - y0;
	int a = Mem2[y0][x][0] * (dy1 / dy2);
	int b = Mem2[y1][x][0] * (dy3 / dy2);
	int c = a + b;
	Mem2[y][x][0] = Mem2[y0][x][0] * (dy1 / dy2) + Mem2[y1][x][0] * (dy3 / dy2);
	Mem2[y][x][1] = Mem2[y0][x][1] * (dy1 / dy2) + Mem2[y1][x][1] * (dy3 / dy2);
	Mem2[y][x][2] = Mem2[y0][x][2] * (dy1 / dy2) + Mem2[y1][x][2] * (dy3 / dy2);
	Mem2Count[y][x]++;
}

// ---------------------------------------------------------------------------
void PutMem1toMem2()
	// трансляция (масштабирование) данных пикселов изображения из входного буфера в буфер для вывода на форму
{
	int i, j, k, i1, j1;
	for (i = 0; i < H; i++)
		for (j = 0; j < W; j++)
			for (k = 0; k < 3; k++) {
				Mem2[i][j][k] = 0;
				Mem2Count[i][j] = 0; // обнуление памяти и счетика
			}

	double kX = inW / (double)W; // коэффициент масштабирования по X
	double kY = inH / (double)H; // коэффициент масштабирования по Y

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
		// усреднение получившихся пикселей (в случае если хотя бы одна из сторон изображения
		// больше соответствующей стороны окна вывода)
		for (i = 0; i < H; i++)
			for (j = 0; j < W; j++)
				for (k = 0; k < 3; k++) {
					if (Mem2Count[i][j] > 1)
						Mem2[i][j][k] /= Mem2Count[i][j];
				}

	int pW = (inW - 1) / kX;
	int pH = (inH - 1) / kY;
	if (inH < H || inW < W) {
		// интерполяция (в случае если хотя бы одна из сторон изображения меньше соответствующей стороны окна вывода)
		for (int j = 0; j <= pW; j++) {
			if (Mem2Count[0][j] == 0)
				LineInterpolationX(j, 0);
			if (Mem2Count[pH][j] == 0)
				LineInterpolationX(j, pH);
		}
		for (int i = 0; i <= pH; i++) {
			if (Mem2Count[i][0] == 0)
				LineInterpolationY(0, i);
			if (Mem2Count[i][pW] == 0)
				LineInterpolationY(pW, i);
		}
		for (int i = 1; i < pH; i++)
			for (int j = 1; j < pW; j++)
				if (Mem2Count[i][j] == 0)
					BeeLineInterpolation(j, i);
	}

}

void  TForm1::ComputeEffectiveBitDepth() {     // Расчет показателей эффективной разрядности
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

void __fastcall TForm1::OpenImage1Click(TObject *Sender) {
	if (OpenPictureDialog1->Execute()) {
		fname = OpenPictureDialog1->FileName;
		FILE *f = fopen(fname.c_str(), "r");
		if (f) {
			fclose(f);
			inbmp1->LoadFromFile(fname); // загрузка изображения из файла
			inH = inbmp1->Height;
			inW = inbmp1->Width;
			long size = inW*inH*3;
			delete [] Mem0;
			delete [] Mem1;
			Mem0 = new Byte[size]; // буфер под пиксели из входного файла
			Mem1 = new Byte[size]; // буфер под пиксели из входного файла

			for (int i = 0; i < inH; i++) {
				struct rgb *ptr = (struct rgb*)inbmp1->ScanLine[i];

				for (int j = 0; j < inW;
				j++) // запись пикселей в память для дальнейшей обработки
				{
					int k = (i * inW + j) * 3;
					Mem0[k] = ptr[j].r;
					Mem0[k + 1] = ptr[j].g;
					Mem0[k + 2] = ptr[j].b;
				}
			}
			int conv = AnsiString(ComboBox1->ItemIndex).ToInt();
			converter(conv, Mem0, Mem1, inH, inW);
			// перенос пикселей из входного изображения в выходное
			PutMem1toMem2();
			for (int i = 0; i < H; i++) // подготовка изображения к выводу на форму
			{
				struct rgb *ptr = (struct rgb*)outbmp->ScanLine[i];
				for (int j = 0; j < W; j++) {
					ptr[j].r = Mem2[i][j][0];
					ptr[j].g = Mem2[i][j][1];
					ptr[j].b = Mem2[i][j][2];
				}
			}
			// Вывод изображения на форму
			Image1->Picture->Graphic = outbmp;

			WWin->Text = inW;
			HWin->Text = inH;

			ComputeEffectiveBitDepth();

		}

	}
}

//----------------------------------------------------------------------------
void __fastcall TForm1::ApplicationEvents1Restore(TObject *Sender)
{
	Form1->Width = 0.85*SW;
	Form1->Height = 0.8*SH;
	Form1->Top = 0.1*SH;
	Form1->Left = 0.05*SW;
	Image1->Width = 0.86*0.85*SW;
	Image1->Height = 0.96*0.8*SH;
	GroupBox1->Left = 0.87*0.85*SW;
	GroupBox2->Left = 0.87*0.85*SW;
	OpenImage1->Left = 0.88*0.85*SW;
	AutoSize =true;
}
//---------------------------------------------------------------------------

void __fastcall TForm1::ApplicationEvents1Minimize(TObject *Sender)
{
	AutoSize =false;
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button1Click(TObject *Sender)
{
    int conv = AnsiString(ComboBox1->ItemIndex).ToInt();
	converter(conv, Mem0, Mem1, inH, inW);
	// перенос пикселей из входного изображения в выходное
	PutMem1toMem2();
	for (int i = 0; i < H;	i++) // подготовка изображения к выводу на форму
	{
		struct rgb *ptr = (struct rgb*)outbmp->ScanLine[i];
		for (int j = 0; j < W; j++) {
			ptr[j].r = Mem2[i][j][0];
			ptr[j].g = Mem2[i][j][1];
			ptr[j].b = Mem2[i][j][2];
		}
	}
		// Вывод изображения на форму
	Image1->Picture->Graphic = outbmp;

	WWin->Text = inW;
	HWin->Text = inH;

	ComputeEffectiveBitDepth();
}
//---------------------------------------------------------------------------

