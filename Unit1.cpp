// ---------------------------------------------------------------------------

#include <vcl.h>
#include <windows.h>
#include <filesystem>
#include <Vcl.Graphics.hpp>
#include <math.h>
#include <cstring>
#include <vector>
#include <algorithm>
#include <SysUtils.hpp>
#include <System.Uitypes.hpp>
#include <System.IOUtils.hpp>
#include <IdGlobal.hpp>
#include <Vcl.Imaging.jpeg.hpp>
#include <Vcl.Imaging.GIFImg.hpp>
#include <Vcl.Imaging.pngimage.hpp>
#pragma hdrstop

#include "Unit1.h"
#include "Unit2.h"

// ---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"

using namespace std;

TForm1 *Form1;
//������������ ������� �������� ����������� (�� ���������)
#define MaxH 10000
#define MaxW 10000

int H, W;     // ������� ���� ��� ����������� �� �����
int inH, inW; // �������� ������� �������� �����������
int hIns, wIns;   // ������� ����������� ���������� � ���� �����
int hCurr, wCurr;   // ������� �����������, ������ �� �������� �������������� ��������
int samplesPerPixel = 3;  //���������� �������� �������
long sizeInBuff;
double scale, scaleIns;   //������� �������, ������� ����������� ���������� � ���� �����


struct rgb // ��������� ��� ������� � bmp 24 ���
{
	Byte b;
	Byte g;
	Byte r;
	// Byte a;
};

//struct rgb32 // // ��������� ��� ������� � bmp 32 ��� (�� ������������)
//{
//	Byte b;
//	Byte g;
//	Byte r;
//	Byte a;
//};
int SW = Screen->Width;
int SH = Screen->Height;
int posK;    //������� ������ ������������� 1
int posM;   //������� ������ ������������� 2

Byte * Mem0; // ��������� �� ��������� ����� ��� ������� �� �������� �����
Byte * Mem1; // ��������� �� ��������� ����� ��� ������� �� �������� ����� (����� ��������������� ��������������, ���� ��� ����)
Byte * MemA;
Byte * MemB;

int*** Mem2;  		// ������������� ������ �������� ��� ������ �� �����
int*** Mem2A;  		// �������������� ������������� ������ �������� ��� ������ �� ����� (��� ������� ������������)
int**  Mem2Count;  	// ������� ���������� ���������� ����������� �������
int**  Mem2ACount;   // ������� ���������� ���������� ����������� ������� ��� ������� Mem2A
int*** Mem3;  		// ������ �������� ��� ������ �� �����
bool memBuff = false;    //����, ��������������� �������� �� ������ ��� �������������� ������ (���������� ���)
int w1, w2, h1, h2; //������� ������, �������������� �������� ����� ��������������� � ���������� �� �����
					// (w - �� ���������, h - �� �����������)
TPoint CenterS;     // ����� � ���� ����������� ������������ ������� ���������� ���������������

Graphics::TBitmap *inbmp;
Graphics::TBitmap *outbmp;

AnsiString fname = "";
bool not_scale;
bool block = true;


// ---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner) : TForm(Owner) {
	//��������� ����� � ����� ������ ����������� ��� ������ ������ ��� ��������
	int FW = 0.9*SW;
	int FH = 0.93*SH;
	Form1->Top = 0.03*SH;
	Form1->Left = 0.05*SW;
	Image1->Width = 0.865*FW;
	Image1->Height = 0.96*FH;
	ProgressBar1->Left = 0.87*FW;
	OpenImage1->Left = 0.9*FW;

	GroupBox1->Left = 0.87*FW;
	GroupBox2->Left = 0.87*FW;
	GroupBox4->Left = 0.87*FW;
	GroupBox5->Left = 0.87*FW;
	GroupBox7->Left = 0.87*FW;
	GroupBox9->Left = 0.87*FW;

	OpenImage1->Top = 0.01*0.96*FH;
	GroupBox1->Top = 0.01*0.96*FH +(OpenImage1->Height)+10;
	GroupBox2->Top = 0.01*0.96*FH +(OpenImage1->Height)+10+(GroupBox1->Height)+5;
	GroupBox9->Top = 0.01*0.96*FH +(OpenImage1->Height)+10+(GroupBox1->Height)+5+(GroupBox2->Height)+5;
	GroupBox4->Top = 0.01*0.96*FH +(OpenImage1->Height)+10+(GroupBox1->Height)+5+(GroupBox2->Height)+5+(GroupBox9->Height)+5;
	GroupBox5->Top = 0.01*0.96*FH +(OpenImage1->Height)+10+(GroupBox1->Height)+5+(GroupBox2->Height)+5+(GroupBox9->Height)+5+(GroupBox4->Height)+5;
	GroupBox7->Top = 0.01*0.96*FH +(OpenImage1->Height)+10+(GroupBox1->Height)+5+(GroupBox2->Height)+5+(GroupBox9->Height)+5+(GroupBox4->Height)+5+(GroupBox5->Height)+5;
	ProgressBar1->Top = 0.96*FH - (ProgressBar1->Height);

	Label4->Caption = AnsiString("�����: ") + SW + "x" + SH;
	//��������� ���� ������� ��������� � ������������� ������ ����������
	posK = 3;
	posM = 10;

	AnsiString filtr;
	char *fc = GetNames(0);
	if (fc!=NULL) {
		filtr = fc;
		ComboBox1->Items->Add(filtr);
		ComboBox1->Text = filtr;
		ComboBox1->ItemIndex = 0;
	}
	bool stop = false;
	int c = 1;
	while(!stop) {
		fc = GetNames(c);
		if (fc!=NULL)   {
			filtr = fc;
			ComboBox1->Items->Add(filtr);
			c++;
		} else stop=true;

	}
	int nBar = GetBarsNum(0);
	for(int i=0; i<nBar; i++) {
		int pos = GetParsPos(0, i);
		switch (i) {
			case 0:   TrackBar3->Visible=true; LastPosition3 = pos; TrackBar3->Position = 100 - pos; break;
			case 1:   TrackBar4->Visible=true; LastPosition4 = pos; TrackBar4->Position = 100 - pos; break;
			default:  TrackBar5->Visible=true; LastPosition5 = pos; TrackBar5->Position = 100 - pos; break;
		}

	}
	//���� ������� �������������
	AnsiString postF;
	char *pfc = GetPostNames(0);
	if (pfc!=NULL) {
		postF = pfc;
		ComboBox3->Items->Add(postF);
		ComboBox3->Text = postF;
		ComboBox3->ItemIndex = 0;
	}
	stop = false;
	c = 1;
	while(!stop) {
		pfc = GetPostNames(c);
		if (pfc!=NULL)   {
			postF = pfc;
			ComboBox3->Items->Add(postF);
			c++;
		} else stop=true;

	}



	LastPosition1 = TrackBar1->Position;
    posTrackBr = -1;
	LastPosition2 = 3;
	LastPosition6 = 6;
	LastPosition7 = 1;

	H = Image1->Height;
	W = Image1->Width;
	inbmp = new Graphics::TBitmap;
	outbmp = new Graphics::TBitmap;
	outbmp->Height = H;
	outbmp->Width = W;
	outbmp->PixelFormat = pf24bit;
	//������ ������� ������ �� ��������� ��� ������ �������� �����������
	sizeInBuff = MaxW*MaxH*samplesPerPixel;

	//�������� ��� (��� � Paint) ��� ����� ������ �����������, ����� ������ ���� ��� ����� ����� �� ����� ���������
	for (int i = 0; i < H;	i++)
	{
		struct rgb *ptr = (struct rgb*)outbmp->ScanLine[i];
		for (int j = 0; j < W; j++) {
			ptr[j].r = 203;
			ptr[j].g = 213;
			ptr[j].b = 228;
		}
	}
	Image1->Picture->Graphic = outbmp;

	//��������� ������ ��� ������� ��� ������ ����������� �� �����. �������� �����������(��������� ������ ������ ����� ���� �����), �� 1 ���
	Mem2 = new int**[H];
	Mem2A = new int**[H];
	Mem3 = new int**[H];
	for(int i = 0; i < H; i++) {
		Mem2[i] = new int*[W];
		Mem2A[i] = new int*[W];
		Mem3[i] = new int*[W];
		for(int j = 0; j < W; j++) {
			Mem2[i][j] = new int[3];
            Mem2A[i][j] = new int[3];
			Mem3[i][j] = new int[3];
		}
	}

	Mem2Count = new int*[H];
	Mem2ACount = new int*[H];
	for(int i = 0; i < H; i++) {
		Mem2Count[i] = new int[W];
        Mem2ACount[i] = new int[W];
	}
}
//----------------------------------------------------------------------------
void BeeLineInterpolation(int*** Mem, int** MemCount, int x, int y) {
	//����� ������� ��������� �������
	int x1, y1, x2, y2, i1, j1;
	int i = x - 1;
	bool stop = false;
	while (!stop) {
		if (MemCount[y - 1][i] > 0)
			stop = true;
		else
			i--;
	}
	x1 = i;
	y2 = y - 1;
	i = x + 1;
	stop = false;
	while (!stop) {
		if (MemCount[y - 1][i] > 0)
			stop = true;
		else
			i++;
	}
	x2 = i;
	int j = y + 1;
	stop = false;
	while (!stop) {
		if (MemCount[j][i] > 0)
			stop = true;
		else
			j++;
	}
	y1 = j;
	if (MemCount[y1][x1] < 1) {
		int i = x1;
		stop = false;
		while (!stop) {
			if (MemCount[y1][i] > 0 && MemCount[y2][i] > 0)
				stop = true;
			else
				i--;
		}
		x1 = i;
	}
	//������� �� ���������� �������� ������
	float dx1, dx2, dx3, dy1, dy2, dy3;
	dx1 = x2 - x1;
	dx2 = x2 - x;
	dx3 = x - x1;
	dy1 = y2 - y1;
	dy2 = y2 - y;
	dy3 = y - y1;
	Mem[y][x][0] = (1.0 / (dx1 * dy1)) * (Mem[y1][x1][0] * dx2 * dy2 +
		Mem[y1][x2][0] * dx3 * dy2 + Mem[y2][x1][0] * dx2 * dy3 +
		Mem[y2][x2][0] * dx3 * dy3);
	Mem[y][x][1] = (1.0 / (dx1 * dy1)) * (Mem[y1][x1][1] * dx2 * dy2 +
		Mem[y1][x2][1] * dx3 * dy2 + Mem[y2][x1][1] * dx2 * dy3 +
		Mem[y2][x2][1] * dx3 * dy3);
	Mem[y][x][2] = (1.0 / (dx1 * dy1)) * (Mem[y1][x1][2] * dx2 * dy2 +
		Mem[y1][x2][2] * dx3 * dy2 + Mem[y2][x1][2] * dx2 * dy3 +
		Mem[y2][x2][2] * dx3 * dy3);

	MemCount[y][x]++;
}

// ----------------------------------------------------------------------------------
void LineInterpolationX(int*** Mem, int** MemCount, int x, int y) {
	//����� ��������� �������
	int x0, y0, x1, y1;
	int i = x - 1;
	bool stop = false;
	while (!stop) {
		if (MemCount[y][i] > 0)
			stop = true;
		else
			i--;
	}
	x0 = i;
	y0 = y;
	i = x + 1;
	stop = false;
	while (!stop) {
		if (MemCount[y][i] > 0)
			stop = true;
		else
			i++;
	}
	x1 = i;
	y1 = y;
	//������� �� ���������� �������� ������
	float dx1, dx2, dx3;
	dx1 = x1 - x;
	dx2 = x1 - x0;
	dx3 = x - x0;
	Mem[y][x][0] = Mem[y][x0][0] * (dx1 / dx2) + Mem[y][x1][0] * (dx3 / dx2);
	Mem[y][x][1] = Mem[y][x0][1] * (dx1 / dx2) + Mem[y][x1][1] * (dx3 / dx2);
	Mem[y][x][2] = Mem[y][x0][2] * (dx1 / dx2) + Mem[y][x1][2] * (dx3 / dx2);
	MemCount[y][x]++;
}


// ---------------------------------------------------------------------------
void LineInterpolationY(int*** Mem, int** MemCount, int x, int y) {
	//����� ��������� �������
	int y0, y1;
	int j = y + 1;
	bool stop = false;
	while (!stop) {
		if (MemCount[j][x] > 0)
			stop = true;
		else
			j++;
	}
	y0 = j;
	j = y - 1;
	stop = false;
	while (!stop) {
		if (MemCount[j][x] > 0)
			stop = true;
		else
			j--;
	}
	y1 = j;
    //������� �� ���������� �������� ������
	float dy1, dy2, dy3;
	dy1 = y1 - y;
	dy2 = y1 - y0;
	dy3 = y - y0;
	int a = Mem[y0][x][0] * (dy1 / dy2);
	int b = Mem[y1][x][0] * (dy3 / dy2);
	int c = a + b;
	Mem[y][x][0] = Mem[y0][x][0] * (dy1 / dy2) + Mem[y1][x][0] * (dy3 / dy2);
	Mem[y][x][1] = Mem[y0][x][1] * (dy1 / dy2) + Mem[y1][x][1] * (dy3 / dy2);
	Mem[y][x][2] = Mem[y0][x][2] * (dy1 / dy2) + Mem[y1][x][2] * (dy3 / dy2);
	MemCount[y][x]++;
}

//----------------------------------------------------------------------------
void BackgroundFill(int x1, int y1, int x2, int y2){
	for (int i = y1; i < y2; i++)
	{
		struct rgb *ptr = (struct rgb*)outbmp->ScanLine[i];
		for (int j = x1; j < x2; j++) {
				ptr[j].r = 203; // �������
				ptr[j].g = 213; // �����
				ptr[j].b = 228; // �������

		}
	}
}
//-------------------------------------------------------------------------------
void ComputeBrightness(int* br) {
	if (br==nullptr)
		return;
	//������ ������� ������� ����������� �� ������, ���������� ��� ������ Mem2===================================================
	vector<int> bPixels;
	bPixels.reserve((h2-h1)*(w2-w1));

	// ���� �������� ������� ������� ��� ���� �������� ������������������� ����������� � ���������� �� � vector bPixels
	// ��� ����������� ������� ���� ������� ������ ����� �����������. ������ ������ ��� ������������� �������������� vector ����� ������ ����
	// ��� ������� ����������� ������ � �����������, ���������� ������������ ����.
	for (int i = h1; i <= h2; i++)
	{
		for (int j = w1; j <= w2; j++) {
			int r = Mem2[i][j][0];
			int g = Mem2[i][j][1];
			int b = Mem2[i][j][2];
			int brightnessPixel = (int)(0.2989 * r + 0.5870 * g + 0.1140 * b);
			bPixels.push_back(brightnessPixel);

		}
		MoveProgress(i, h2-h1);

	}

	int cbuff = 0;
	int cpbf = 0;
	int fb = 0;
	int brightness = -1;

	do
	{
		if(cbuff == bPixels.size() || cpbf==1000000) {
		   int averageB = fb/cpbf;
		   if (brightness>0) {
			   brightness += averageB;
			   brightness /=2;
		   }else brightness = averageB;
		} else {
		   int brightnessPixel = bPixels[cbuff];
		   fb +=brightnessPixel;
		   cpbf++;
		}
		if(cpbf==1000000 && cbuff<bPixels.size()) {
		   int brightnessPixel = bPixels[cbuff];
		   fb = brightnessPixel;
		   cpbf=1;
		}
	   cbuff++;
	}
	while (cbuff <= bPixels.size());

	*br = brightness;
}
//----------------------------------------------------------------------------------

void ComputeMoveForScale(double* scaleNew, int* dx, int* dy) {
	int movW = (W-wIns)/2;
	int movH = (H-hIns)/2;
	double sNew = *scaleNew;
	//����� CenterS ���������� �� �������� �������� ������ ����� � ����� ����� ���� ������ �����������, � ������ � ��� ��, ��� ��� ����� �� ��������� ��������.
	// ������� ���������� �� �������� ����������� ������������ ������� �� ��������� �� ������� ����������� (xNature, yNature)
	// �� ��������� �� ���� ������ ��� ��� �������� scaleIns
	int xNature = (CenterS.x-movW)*scaleIns;
	int yNature = (CenterS.y-movH)*scaleIns;
	if(sNew<1) sNew=1;
	int x1 = xNature/sNew;
	int y1 = yNature/sNew;
	*dx = x1-(CenterS.x-movW);
	*dy = y1-(CenterS.y-movH);
	*scaleNew = sNew;
}

//------------------------------------------------------------------------------
void PutMem1toMem2(double scale, int dx, int dy, int* br){
	// � ���� ������� ����������� ������������� �� ���� ������ � ��������� scale � � ��������� � ����� � ������� �������������� ������� dx, dy
	// ��� ��������� ��������� ����������� ����� ����� 0
	//��� dx, dy - �������������� ������ �� ���� ��������� ��� ���������������, � ���������, �� ���������� ����� �� �������� (dx,dy ����������� ������������
	//��������� ����� ������� ��� ��������� ��������� �����������, �.�. ��� �������� scaleIns)
	int i, j, k, i1, j1;
	for (i = 0; i < H; i++)
		for (j = 0; j < W; j++)
			for (k = 0; k < samplesPerPixel; k++) {
				Mem2[i][j][k] = 0;
				Mem2Count[i][j] = 0; // ��������� ������ � �������
				Mem2A[i][j][k] = 0;
				Mem2ACount[i][j] = 0; // ��������� ������ � �������
			}

	// ��������� ����������� ��������������� � ����� �������� �����������
	int imgW = static_cast<int>(inW / scale); // ����� ������ �����������
	int imgH = static_cast<int>(inH / scale); // ����� ������ �����������
	int movW = (W-wIns)/2;
	int movH = (H-hIns)/2;

	w1 = w2 = CenterS.x;
	h1 = h2 = CenterS.y;
	for (int i = 0; i < inH; i++) {
		for (int j = 0; j < inW; j++) {
			int k = (i * inW + j) * samplesPerPixel;
			// �������� ���������� ��������� (Y,X,�������� �����) � ����������
			i1 = i / scale + movH - dy;  // ������� ������� � ������� ��������� �����������
			j1 = j / scale + movW - dx; // ������� ������� � ������ ��������� �����������
			if (i1<H && i1>=0 && j1<W && j1>=0) {
				if(w1>j1)  w1=j1;	if(w2<j1)  w2=j1;
				if(h1>i1)  h1=i1;   if(h2<i1)  h2=i1;
				Mem2[i1][j1][0] += Mem1[k]; // �������
				Mem2[i1][j1][1] += Mem1[k + 1]; // �����
				Mem2[i1][j1][2] += Mem1[k + 2]; // �������
				Mem2Count[i1][j1]++;
				//��� ������������
				Mem2A[i1][j1][0] += Mem0[k]; // �������
				Mem2A[i1][j1][1] += Mem0[k + 1]; // �����
				Mem2A[i1][j1][2] += Mem0[k + 2]; // �������
				Mem2ACount[i1][j1]++;
			}

		}
		MoveProgress(i, inH);
	}

	hCurr = imgH;  wCurr = imgW;

	if (inH > H || inW > W)
		// ���������� ������������ �������� (� ������ ���� ���� �� ���� �� ������ �����������
		// ������ ��������������� ������� ���� ������)
		for (i = h1; i <= h2; i++)  {
			for (j = w1; j <= w2; j++)
				for (k = 0; k < samplesPerPixel; k++) {
					if (Mem2Count[i][j] > 1)
						Mem2[i][j][k] /= Mem2Count[i][j];
					if (Mem2ACount[i][j] > 1)
						Mem2A[i][j][k] /= Mem2ACount[i][j];
				}
			MoveProgress(i, h2-h1);
		}

	ComputeBrightness(br);
}
// ---------------------------------------------------------------------------
void InsMem1toMem2(int* br)
	// ���������� (���������������) ������ �������� ����������� �� �������� ������ � ����� ��� ������ �� �����
	// � ���� ������� �����������, ���������� �� ���� ����� ��� �������, ����������� � ���� ������ � ����������� ��������� � � ����������, ������� ��������������� ��� scaleIns
{
	int i, j, k, i1, j1;
	for (i = 0; i < H; i++)
		for (j = 0; j < W; j++)
			for (k = 0; k < samplesPerPixel; k++) {
				Mem2[i][j][k] = 0;
				Mem2Count[i][j] = 0; // ��������� ������ � �������
                Mem2A[i][j][k] = 0;
				Mem2ACount[i][j] = 0; // ��������� ������ � �������
			}

	double kX = inW / (double)W; // ����������� ��������������� �� X
	double kY = inH / (double)H; // ����������� ��������������� �� Y

	// ����������, ����� ����������� ����� ����� � ����������� �� ����������� ������ ����������� � �����
	scale = (kX > kY) ? kX : kY;

	// ��������� ����������� ��������������� � ����� �������� �����������
	int imgW = static_cast<int>(inW / scale); // ����� ������ �����������
	int imgH = static_cast<int>(inH / scale); // ����� ������ �����������
	int movW = (W-imgW)/2;
	int movH = (H-imgH)/2;

	for (int i = 0; i < inH; i++) {
		for (int j = 0; j < inW; j++) {
			int k = (i * inW + j) * samplesPerPixel;
			// �������� ���������� ��������� (Y,X,�������� �����) � ����������
			i1 = i / scale + movH;  // ������� ������� � ������� ��������� �����������
			j1 = j / scale + movW; // ������� ������� � ������ ��������� �����������
			Mem2[i1][j1][0] += Mem1[k]; // �������
			Mem2[i1][j1][1] += Mem1[k + 1]; // �����
			Mem2[i1][j1][2] += Mem1[k + 2]; // �������
			Mem2Count[i1][j1]++;
			//��� ������������
			Mem2A[i1][j1][0] += Mem0[k]; // �������
			Mem2A[i1][j1][1] += Mem0[k + 1]; // �����
			Mem2A[i1][j1][2] += Mem0[k + 2]; // �������
			Mem2ACount[i1][j1]++;
		}
		MoveProgress(i, inH);
	}

	w1 = movW;  w2 = (inW-1)/scale + movW;
	h1 = movH;  h2 = (inH-1)/scale + movH;
	hCurr = h2-h1;  wCurr = w2-w1;
	CenterS.X = W/2;
	CenterS.Y = H/2;

	if (inH > H || inW > W)
		// ���������� ������������ �������� (� ������ ���� ���� �� ���� �� ������ �����������
		// ������ ��������������� ������� ���� ������)
		for (i = h1; i <= h2; i++)  {
			for (j = w1; j <= w2; j++)
				for (k = 0; k < samplesPerPixel; k++) {
					if (Mem2Count[i][j] > 1)
						Mem2[i][j][k] /= Mem2Count[i][j];
					//��� ������������
					if (Mem2ACount[i][j] > 1)
						Mem2A[i][j][k] /= Mem2ACount[i][j];
				}
			MoveProgress(i, h2-h1);
		}

	if (inH < H || inW < W) {
		// ������������ (� ������ ���� ���� �� ���� �� ������ ����������� ������ ��������������� ������� ���� ������)
		for (int j = w1; j <= w2; j++) {
			if (Mem2Count[h1][j] == 0)
				LineInterpolationX(Mem2, Mem2Count, j, h1);
			if (Mem2Count[h2][j] == 0)
				LineInterpolationX(Mem2, Mem2Count, j, h2);
		}
		for (int i = h1; i <= h2; i++) {
			if (Mem2Count[i][w1] == 0)
				LineInterpolationY(Mem2, Mem2Count, w1, i);
			if (Mem2Count[i][w2] == 0)
				LineInterpolationY(Mem2, Mem2Count, w2, i);
		}
		for (int i = h1+1; i < h2; i++)
			for (int j = w1+1; j < w2; j++)
				if (Mem2Count[i][j] == 0)
					BeeLineInterpolation(Mem2, Mem2Count, j, i);

		//��������� �������������� ���� �������� ��� ������������
		for (int j = w1; j <= w2; j++) {
			if (Mem2ACount[h1][j] == 0)
				LineInterpolationX(Mem2A, Mem2ACount, j, h1);
			if (Mem2ACount[h2][j] == 0)
				LineInterpolationX(Mem2A, Mem2ACount, j, h2);
		}
		for (int i = h1; i <= h2; i++) {
			if (Mem2ACount[i][w1] == 0)
				LineInterpolationY(Mem2A, Mem2ACount, w1, i);
			if (Mem2ACount[i][w2] == 0)
				LineInterpolationY(Mem2A, Mem2ACount, w2, i);
		}
		for (int i = h1+1; i < h2; i++)
			for (int j = w1+1; j < w2; j++)
				if (Mem2ACount[i][j] == 0)
					BeeLineInterpolation(Mem2A, Mem2ACount, j, i);
	}

	ComputeBrightness(br);
}
//----------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------
void PutMem2toMem3(int posB, int posK, int m, int posA){
	int posB1 = Form1->posTrackBr;
	double kb = static_cast<double>(posB) / posB1; // ����������� ��������� �������
	double k = posK / 3.0;
	m*=20;
	double A = posA / 100.0;
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			for (int k = 0; k < samplesPerPixel; k++)
				Mem3[i][j][k] = 0;	// ��������� ������ � �������

	//������������:  �����_������� = (������_������� * A) + (�������_������� * (1 - A))

	for (int i = h1; i <= h2; i++)
	{
		for (int j = w1; j <= w2; j++)

			for(int c=0; c<3; c++) {
				//��������� ������������� ������������
				int x = Mem2[i][j][c];
				x = k*(x-m)+m;
				// ����������� �������� ������� ������� � ��������� �� 0 �� 255
				x = std::min(std::max(x, 0), 255);
				//��������� ������������ ��������� �������
				int y = static_cast<int>(x * kb);
				// ����������� �������� ������� ������� � ��������� �� 0 �� 255
				y = std::min(std::max(y, 0), 255);
				//e
				double z1 = y * A;
				double z2 = Mem2A[i][j][c] * (1-A);
				int z =z1 + z2;
				// ����������� �������� ������� ������� � ��������� �� 0 �� 255
				z = std::min(std::max(z, 0), 255);
				Mem3[i][j][c] = z;
			}

		MoveProgress(i, h2-h1);
	}
}
//--------------------------------------------------------------------------------------------------------------
void  TForm1::ComputeEffectiveBitDepth() {     // ������ ����������� ����������� �����������
	int BSize = inH * inW * samplesPerPixel;
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

//------------------------------------------------------------------------------
void TForm1::Blocking(){
	AutoSize =false;
	GroupBox1->Enabled = false;
	GroupBox2->Enabled = false;
	GroupBox3->Enabled = false;
	GroupBox4->Enabled = false;
	ProgressBar1->Position = 3;
	ProgressBar1->Visible = true;
    Image1->Height = 0.96*0.93*SH;
	Application->ProcessMessages();
	OpenImage1->Enabled = false;
	SaveImg->Enabled = false;
	TrackBar2->Enabled = false;
	TrackBar3->Enabled = false;
	TrackBar4->Enabled = false;
	TrackBar5->Enabled = false;
	TrackBar6->Enabled = false;
	Image1->Height = 0.96*0.93*SH;
	block = true;

}
//------------------------------------------------------------------------------
void TForm1::UnBlocking(){
	GroupBox1->Enabled = true;
	GroupBox2->Enabled = true;
	GroupBox3->Enabled = true;
	if (scaleIns>1)
		GroupBox4->Enabled = true;
	else
		GroupBox4->Enabled = false;
	ProgressBar1->Visible = false;
	OpenImage1->Enabled = true;
	SaveImg->Enabled = true;
	Button1->Enabled = true;
	TrackBar2->Enabled = true;
	TrackBar3->Enabled = true;
	TrackBar4->Enabled = true;
	TrackBar5->Enabled = true;
	TrackBar6->Enabled = true;


	Image1->Width = 0.865*0.9*SW;
	Image1->Height = 0.96*0.93*SH;

	block = false;
	AutoSize = true;
}
//------------------------------------------------------------------------------
void TForm1::ButtonsScaleConfiguration(){
	if(inH<=hCurr || inW<=wCurr) {
		Button2->Enabled = false;
		Button5->Enabled = false;
	} else {
		Button2->Enabled = true;
		Button5->Enabled = true;
	}

	if(scale==scaleIns) {
		Button3->Enabled = false;
		Button4->Enabled = false;
	} else {
		Button3->Enabled = true;
		Button4->Enabled = true;
	}
	if (scaleIns<1) {
        Button2->Enabled = false;
		Button5->Enabled = false;
        Button3->Enabled = false;
		Button4->Enabled = false;
    }

}

//------------------------------------------------------------------------------
void MoveProgress(int i, int inH) {
	// ���������, �������� �� ������� ������� i ������� inH/30
	if ((Form1!=nullptr)&&((i + 1) % (inH / 30) == 0))
	{
		// ������������� �������� ProgressBar1
		Form1->ProgressBar1->Position = (i + 1) * 100 / inH;
		Application->ProcessMessages();
	}
}
//------------------------------------------------------------------------------
void TForm1::SetTpCursor(){
		if(not_scale) {
			Cursor = crDefault; return;
		}
		TPoint MousePos = Image1->ScreenToClient(Mouse->CursorPos);
		int X = MousePos.x;
		int Y = MousePos.y;
		if(!not_scale && (Y<h1 || Y>h2 || X<w1 || X>w2)) {
			Cursor = crDefault;  return;
		}
		if (scale==scaleIns && !not_scale)  {
			Cursor = crSizeAll;  return;
		}
		if (scale!=scaleIns)
			Cursor = crCross;
}
//-----------------------------------------------------------------------------
/*
������� ExtractFileExt ��������� �� ������� ����� ����� ���������� �����.
AnsiLowerCase() - ��������� ������ � ������ �������.
*/
void __fastcall TForm1::OpenImage1Click(TObject *Sender) {
	if (OpenPictureDialog1->Execute()) {
		AnsiString fname1 = OpenPictureDialog1->FileName;
		String fileExt = AnsiLowerCase(ExtractFileExt(fname1));
		if (TFile::Exists(L"C:\\Windows\\Temp\\atlas_foto.bmp"))
		{
			// ������� ����
			TFile::Delete(L"C:\\Windows\\Temp\\atlas_foto.bmp");
		}
		Blocking();
		OpenImage1->Caption = "��������� ����...";
		if (fileExt != ".bmp") {
			// ���� � exe-����� � ��������� ��������� ������
			UnicodeString filePath1 = L"magick";
			AnsiString QuoteFname = AnsiString("\"") + fname1 + "\"";
			UnicodeString argsAS = QuoteFname + L" C:\\Windows\\Temp\\atlas_foto.bmp";

			// ����������� ��������� � C-������
			wchar_t* commandLine = (filePath1 + L" " + argsAS).w_str();

			STARTUPINFO si = { sizeof(si) };
			PROCESS_INFORMATION pi;
			if (CreateProcess(nullptr,   // ��� ������ (nullptr ������ ������������ ��������� ������)
				commandLine,             // ��������� ������
				nullptr,                 // ������� ������������
				nullptr,                 // ����� ������������
				FALSE,                   // �� ����������� �����������
				0,                       // ����� ��������
				nullptr,                 // ���������� ���������
				nullptr,                 // ������� ����������
				&si,                     // ��������� ���������� � �������
				&pi       ))
			{
				WaitForSingleObject(pi.hProcess, INFINITE);

				CloseHandle(pi.hProcess);
				CloseHandle(pi.hThread);
			}
		}

		if (fileExt != ".bmp")
			fname = "C:\\Windows\\Temp\\atlas_foto.bmp";
		else
			fname = fname1;

		fileExt = AnsiLowerCase(ExtractFileExt(fname));
		int i1, j1, brPixel;
		if (fileExt == ".bmp") {
			inbmp->LoadFromFile(fname); // �������� ����������� �� �����
			ProgressBar1->Position = 100;
			Application->ProcessMessages();
			inH = inbmp->Height;
			inW = inbmp->Width;
			long size = inW*inH*samplesPerPixel;
			if(!memBuff) {
				if (size>sizeInBuff) {
					Mem0 = new Byte[size]; // ��������� ����� ��� ������� �� �������� �����
					Mem1 = new Byte[size]; // ��������� ����� ��� ������� �� �������� �����
					MemA = new Byte[size];
					MemB = new Byte[size];
				} else {
					Mem0 = new Byte[sizeInBuff]; // ��������� ����� ��� ������� �� �������� �����
					Mem1 = new Byte[sizeInBuff]; // ��������� ����� ��� ������� �� �������� �����
					MemA = new Byte[sizeInBuff];
					MemB = new Byte[sizeInBuff];
				}
				memBuff = true;
			}  else	{
				if (size>sizeInBuff) {
					Mem0 = new Byte[size]; // ��������� ����� ��� ������� �� �������� �����
					Mem1 = new Byte[size]; // ��������� ����� ��� ������� �� �������� �����
					MemA = new Byte[size];
					MemB = new Byte[size];
				}
			}

			for (int i = 0; i < inH; i++) {
				struct rgb *ptr = (struct rgb*)inbmp->ScanLine[i];
				for (int j = 0; j < inW; j++) // ������ �������� � ��������� ����� ��� ���������� ���������
				{
					int k = (i * inW + j) * samplesPerPixel;
					Mem0[k] = ptr[j].r;
					Mem0[k + 1] = ptr[j].g;
					Mem0[k + 2] = ptr[j].b;
					MemA[k] = ptr[j].r;
					MemA[k + 1] = ptr[j].g;
					MemA[k + 2] = ptr[j].b;
				}
				MoveProgress(i, inH);
			}
			ProgressBar1->Position = 100;
			Application->ProcessMessages();

		}else {
			ShowMessage("���������������� ��� ������");
			OpenImage1->Enabled = true;
			OpenImage1->Caption = "�������";
			return;
		}

		ComboBox1->ItemIndex = 0;
		ComboBox2->ItemIndex = 0;
        ComboBox3->ItemIndex = 0;
		int conv = AnsiString(ComboBox1->ItemIndex).ToInt();
		int Chan = ComboBox2->ItemIndex;
		int post = ComboBox3->ItemIndex;
		converter(conv, Mem0, Mem1, inH, inW, &posK, &posM, post, LastPosition3, LastPosition4, LastPosition5, Chan);
		ProgressBar1->Position = 100;
		Application->ProcessMessages();
		// ������� �������� �� �������� ����������� � ��������
		int brightness = -1;
		InsMem1toMem2(&brightness);
		LastPosition1 = (brightness*100)/255;
		LastPosition2 = posK;
		LastPosition6 = posM;
		LastPosition7 = 100;
		TrackBar1->Position = (brightness*100)/255;
		TrackBar2->Position = posK;
		TrackBar6->Position = posM;
		TrackBar7->Position = 100;
		posTrackBr = TrackBar1->Position;
		ProgressBar1->Position = 3;
		Application->ProcessMessages();
		PutMem2toMem3(posTrackBr, TrackBar2->Position, TrackBar6->Position, TrackBar7->Position);

		for (int i = h1; i <= h2; i++) // ���������� ����������� � ������ �� �����
		{
			struct rgb *ptr = (struct rgb*)outbmp->ScanLine[i];
			for (int j = w1; j <= w2; j++) {
				ptr[j].r = Mem3[i][j][0];
				ptr[j].g = Mem3[i][j][1];
				ptr[j].b = Mem3[i][j][2];

			}
			MoveProgress(i, h2-h1);

		}
		//�������� ������������ ������ ������ ����
		BackgroundFill(0, 0, w1, H);
		BackgroundFill(w2+1, 0, W, H);
		BackgroundFill(0, 0, W, h1);
		BackgroundFill(0, h2+1, W, H);

		// ����� ����������� �� �����
		Image1->Picture->Graphic = outbmp;

		WWin->Text = inW;
		HWin->Text = inH;

		ComputeEffectiveBitDepth();

		UnBlocking();
		OpenImage1->Caption = "�������";
		hIns = h2-h1;  wIns = w2-w1;
		scaleIns = scale;
		ButtonsScaleConfiguration();
		if(scale>1) not_scale=false;
		else         not_scale=true;
	}
	SetTpCursor();

}

//----------------------------------------------------------------------------
void __fastcall TForm1::ApplicationEvents1Restore(TObject *Sender)
{
	//��������� ����� � ����� ������ ����������� ��� ������ ������ ��� �������������� ����� �� ���������� ���������
	int FW = 0.9*SW;
	int FH = 0.93*SH;
	Form1->Top = 0.03*SH;
	Form1->Left = 0.05*SW;
	Image1->Width = 0.865*FW;
	Image1->Height = 0.96*FH;
	ProgressBar1->Left = 0.87*FW;
	OpenImage1->Left = 0.9*FW;

	GroupBox1->Left = 0.87*FW;
	GroupBox2->Left = 0.87*FW;
	GroupBox4->Left = 0.87*FW;
	GroupBox5->Left = 0.87*FW;
	GroupBox7->Left = 0.87*FW;
	GroupBox9->Left = 0.87*FW;

	OpenImage1->Top = 0.01*0.96*FH;
	GroupBox1->Top = 0.01*0.96*FH +(OpenImage1->Height)+10;
	GroupBox2->Top = 0.01*0.96*FH +(OpenImage1->Height)+10+(GroupBox1->Height)+5;
	GroupBox9->Top = 0.01*0.96*FH +(OpenImage1->Height)+10+(GroupBox1->Height)+5+(GroupBox2->Height)+5;
	GroupBox4->Top = 0.01*0.96*FH +(OpenImage1->Height)+10+(GroupBox1->Height)+5+(GroupBox2->Height)+5+(GroupBox9->Height)+5;
	GroupBox5->Top = 0.01*0.96*FH +(OpenImage1->Height)+10+(GroupBox1->Height)+5+(GroupBox2->Height)+5+(GroupBox9->Height)+5+(GroupBox4->Height)+5;
	GroupBox7->Top = 0.01*0.96*FH +(OpenImage1->Height)+10+(GroupBox1->Height)+5+(GroupBox2->Height)+5+(GroupBox9->Height)+5+(GroupBox4->Height)+5+(GroupBox5->Height)+5;
	ProgressBar1->Top = 0.96*FH - (ProgressBar1->Height);

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
	Blocking();
	Button1->Caption = "�����...";

	for (int i = 0; i < inH; i++) {
		for (int j = 0; j < inW; j++)
		{
			int k = (i * inW + j) * 3;
			int r = Mem1[k];
			int g = Mem1[k + 1];
			int b = Mem1[k + 2];
			Mem0[k] = r;
			Mem0[k + 1] = g;
			Mem0[k + 2] = b;
		}
		MoveProgress(i, inH);
	}
	Button1->Caption = "�������";
	UnBlocking();
}

//--------------------------------------------------------------------------
void __fastcall TForm1::TrackBar1Change(TObject *Sender)
{
	int pos2 = TrackBar1->Position; // ����� ������� ������ �������
	GroupBox3->Caption = AnsiString("�������: ") + pos2;
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
void __fastcall TForm1::Timer1Timer(TObject *Sender)
{
	if(block)
		return;

	bool flag1;
	flag1 = false;
    int m;
	if (TrackBar1->Position != LastPosition1)
	{
		LastPosition1 = TrackBar1->Position; // ����� ������� ������ �������
		flag1 = true;
	}

	if (TrackBar2->Position != LastPosition2)
	{
		LastPosition2 = TrackBar2->Position; //// ����� ������� ������1 �������������
		flag1 = true;
	}

	if (TrackBar6->Position != LastPosition6)
	{
		m = LastPosition6 = TrackBar6->Position; //// ����� ������� ������2 �������������
		flag1 = true;
	} else
		m = LastPosition6;

	if (TrackBar7->Position != LastPosition7)
	{
		LastPosition7 = TrackBar7->Position; //// ����� ������� ������ ������������
		flag1 = true;
	}

	if(flag1)
	{
		Blocking();
		ProgressBar1->Position = 3;
		Application->ProcessMessages();
		PutMem2toMem3(LastPosition1, LastPosition2, m, LastPosition7);
		for (int i = h1; i <= h2; i++) // ���������� ����������� � ������ �� �����
		{
			struct rgb *ptr = (struct rgb*)outbmp->ScanLine[i];
			for (int j = w1; j <= w2; j++) {
				ptr[j].r = Mem3[i][j][0];
				ptr[j].g = Mem3[i][j][1];
				ptr[j].b = Mem3[i][j][2];
			}
			MoveProgress(i, h2-h1);
		}

		// ����� ����������� �� �����
		Image1->Picture->Graphic = outbmp;

	   //	ComputeEffectiveBitDepth();

		UnBlocking();
	}

	bool flag3, flag4, flag5;
	flag3 = flag4 = flag5 = false;
	if (100-(TrackBar3->Position) != LastPosition3 && TrackBar3->Visible)
	{
		LastPosition3 = 100-(TrackBar3->Position); // ����� ������� ������ ��������� 1
		flag3 = true;
	}

	if (100-(TrackBar4->Position)!= LastPosition4  && TrackBar4->Visible)
	{
		LastPosition4 = 100-(TrackBar4->Position); //// ����� ������� ������ ��������� 2
		flag4 = true;
	}

	if (100-(TrackBar5->Position)!= LastPosition5  && TrackBar5->Visible)
	{
		LastPosition5 = 100-(TrackBar5->Position); //// ����� ������� ������ ��������� 3
		flag5 = true;
	}

	if(flag3 || flag4 || flag5)
	{
		ProgressBar1->Position = 3;
		Application->ProcessMessages();
		DoChangeMem0ToMem1();
	}
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Button4Click(TObject *Sender)
{
	Blocking();
	InsMem1toMem2(nullptr);

	ProgressBar1->Position = 3;
	Application->ProcessMessages();
	PutMem2toMem3(TrackBar1->Position, TrackBar2->Position, TrackBar6->Position, TrackBar7->Position);
	for (int i = h1; i <= h2; i++) // ���������� ����������� � ������ �� �����
		{
			struct rgb *ptr = (struct rgb*)outbmp->ScanLine[i];
			for (int j = w1; j <= w2; j++) {
				ptr[j].r = Mem3[i][j][0];
				ptr[j].g = Mem3[i][j][1];
				ptr[j].b = Mem3[i][j][2];

			}
			MoveProgress(i, h2-h1);
		}

	//�������� ������������ ������ ������ ����
	BackgroundFill(0, 0, w1, H);
	BackgroundFill(w2+1, 0, W, H);
	BackgroundFill(0, 0, W, h1);
	BackgroundFill(0, h2+1, W, H);

	// ����� ����������� �� �����
	Image1->Picture->Graphic = outbmp;
	UnBlocking();
	ButtonsScaleConfiguration();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button5Click(TObject *Sender)
{
	Blocking();
	int dx = 0; // ������������� ���������� dx
	int dy = 0; // ������������� ���������� dy
	double scaleNew=1;
	ComputeMoveForScale(&scaleNew, &dx, &dy);
	PutMem1toMem2(scaleNew, dx, dy, nullptr);
	ProgressBar1->Position = 3;
	Application->ProcessMessages();
	PutMem2toMem3(TrackBar1->Position, TrackBar2->Position, TrackBar6->Position, TrackBar7->Position);
	for (int i = h1; i <= h2; i++) // ���������� ����������� � ������ �� �����
		{
			struct rgb *ptr = (struct rgb*)outbmp->ScanLine[i];
			for (int j = w1; j <= w2; j++) {
				ptr[j].r = Mem3[i][j][0];
				ptr[j].g = Mem3[i][j][1];
				ptr[j].b = Mem3[i][j][2];

			}
			MoveProgress(i, h2-h1);
		}
	//�������� ������������ ������ ������ ����
	BackgroundFill(0, 0, w1, H);
	BackgroundFill(w2+1, 0, W, H);
	BackgroundFill(0, 0, W, h1);
	BackgroundFill(0, h2+1, W, H);

	// ����� ����������� �� �����
	Image1->Picture->Graphic = outbmp;
	UnBlocking();
	scale=1;
	ButtonsScaleConfiguration();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button2Click(TObject *Sender)
{
	double scaleOld = scale;
	scale--;
	if(scale<1) scale=1;
	Blocking();
	int dx = 0; // ������������� ���������� dx
	int dy = 0; // ������������� ���������� dy
	ComputeMoveForScale(&scale, &dx, &dy);
	PutMem1toMem2(scale, dx, dy, nullptr);
	ProgressBar1->Position = 3;
	Application->ProcessMessages();
	PutMem2toMem3(TrackBar1->Position, TrackBar2->Position, TrackBar6->Position, TrackBar7->Position);
	for (int i = h1; i <= h2; i++) // ���������� ����������� � ������ �� �����
		{
			struct rgb *ptr = (struct rgb*)outbmp->ScanLine[i];
			for (int j = w1; j <= w2; j++) {
				ptr[j].r = Mem3[i][j][0];
				ptr[j].g = Mem3[i][j][1];
				ptr[j].b = Mem3[i][j][2];

			}
			MoveProgress(i, h2-h1);
		}
	//�������� ������������ ������ ������ ����
	BackgroundFill(0, 0, w1, H);
	BackgroundFill(w2+1, 0, W, H);
	BackgroundFill(0, 0, W, h1);
	BackgroundFill(0, h2+1, W, H);

	// ����� ����������� �� �����
	Image1->Picture->Graphic = outbmp;
	UnBlocking();
	ButtonsScaleConfiguration();

}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button3Click(TObject *Sender)
{
	double scaleOld = scale;
	scale++;
	if(scale>=scaleIns) scale=scaleIns;
	Blocking();
    int dx = 0; // ������������� ���������� dx
	int dy = 0; // ������������� ���������� dy
	ComputeMoveForScale(&scale, &dx, &dy);
	PutMem1toMem2(scale, dx, dy, nullptr);
	ProgressBar1->Position = 3;
	Application->ProcessMessages();
	PutMem2toMem3(TrackBar1->Position, TrackBar2->Position, TrackBar6->Position, TrackBar7->Position);
	for (int i = h1; i <= h2; i++) // ���������� ����������� � ������ �� �����
		{
			struct rgb *ptr = (struct rgb*)outbmp->ScanLine[i];
			for (int j = w1; j <= w2; j++) {
				ptr[j].r = Mem3[i][j][0];
				ptr[j].g = Mem3[i][j][1];
				ptr[j].b = Mem3[i][j][2];

			}
			MoveProgress(i, h2-h1);
		}
	//�������� ������������ ������ ������ ����
	BackgroundFill(0, 0, w1, H);
	BackgroundFill(w2+1, 0, W, H);
	BackgroundFill(0, 0, W, h1);
	BackgroundFill(0, h2+1, W, H);

	// ����� ����������� �� �����
	Image1->Picture->Graphic = outbmp;
	UnBlocking();
	ButtonsScaleConfiguration();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Image1MouseMove(TObject *Sender, TShiftState Shift, int X,
		  int Y)
{
	SetTpCursor();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Image1MouseLeave(TObject *Sender)
{
	 // ������� ����������� ������ Windows 10, ����� ���� ��� TImage
		Cursor = crDefault;
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Image1Click(TObject *Sender)
{
	if(fname=="")
		return;
	if (not_scale)
		return;
	not_scale = true;   //������ �� �������� ������: ��� ������ ������ �� ���� �� �������, ���� ������ ��������� �� ����������
	TPoint MousePos = Image1->ScreenToClient(Mouse->CursorPos);
	int x = MousePos.x;
	int y = MousePos.y;
	CenterS.X = x; // ��������� ������ �������� X
	CenterS.Y = y; // ��������� ������ �������� Y
	int movW = (W-wCurr)/2;
	int movH = (H-hCurr)/2;



	Blocking();
	if (scale==scaleIns)  {
		scale=scaleIns*0.4;
		int dx = 0; // ������������� ���������� dx
		int dy = 0; // ������������� ���������� dy
		ComputeMoveForScale(&scale, &dx, &dy);
		PutMem1toMem2(scale, dx, dy, nullptr);
	}
	else
		InsMem1toMem2(nullptr);

	ProgressBar1->Position = 3;
	Application->ProcessMessages();
	PutMem2toMem3(TrackBar1->Position, TrackBar2->Position, TrackBar6->Position, TrackBar7->Position);
	for (int i = h1; i <= h2; i++) // ���������� ����������� � ������ �� �����
		{
			struct rgb *ptr = (struct rgb*)outbmp->ScanLine[i];
			for (int j = w1; j <= w2; j++) {
				ptr[j].r = Mem3[i][j][0];
				ptr[j].g = Mem3[i][j][1];
				ptr[j].b = Mem3[i][j][2];

			}
			MoveProgress(i, h2-h1);
		}
	//�������� ������������ ������ ������ ����
	BackgroundFill(0, 0, w1, H);
	BackgroundFill(w2+1, 0, W, H);
	BackgroundFill(0, 0, W, h1);
	BackgroundFill(0, h2+1, W, H);

	// ����� ����������� �� �����
	Image1->Picture->Graphic = outbmp;
	UnBlocking();
	ButtonsScaleConfiguration();
	SetTpCursor();
	not_scale = false;
}
//---------------------------------------------------------------------------

void __fastcall TForm1::SaveImgClick(TObject *Sender)
{
	if (RadioGroup1->ItemIndex==-1)   {
		ShowMessage("���������� ����������� �� �����, ��������� �� ������ ������� ����������");
		return;
	}
	int asSave = RadioGroup1->ItemIndex;
	if (SavePictureDialog1->Execute()) {
		String FileNameS = SavePictureDialog1->FileName;
		AnsiString file = AnsiString(FileNameS);
		TBitmap* savebmp = new Graphics::TBitmap;
		savebmp->PixelFormat = pf24bit;
		switch(asSave) {
			case 0:
				savebmp->Height = h2-h1+1;
				savebmp->Width = w2-w1+1;

				for (int i = 0; i <= h2-h1; i++) // ���������� ����������� � ������
				{
					struct rgb *ptr = (struct rgb*)savebmp->ScanLine[i];
					for (int j = 0; j <= w2-w1; j++) {
						ptr[j].r = Mem3[i+h1][j+w1][0];
						ptr[j].g = Mem3[i+h1][j+w1][1];
						ptr[j].b = Mem3[i+h1][j+w1][2];

					}
					MoveProgress(i, h2-h1);
				}
				break;
			default:
				int pos2 = TrackBar1->Position; // ����� ������� ������ �������
				int pos1 = posTrackBr;
				double kb = static_cast<double>(pos2) / pos1; // ����������� ��������� �������
            	double kp = (TrackBar2->Position)/3.0;     	//������� ������ 1 �������������
				int mp = (TrackBar6->Position)*20;     		//������� ������ 2 ������������� (������� ������)
				double A = (TrackBar7->Position) / 100.0;   //����������� ������������

				savebmp->Height = inH;
				savebmp->Width = inW;

				for (int i = 0; i < inH; i++) // ���������� ����������� � ������
				{
					struct rgb *ptr = (struct rgb*)savebmp->ScanLine[i];
					for (int j = 0; j < inW; j++) {
						int k = (i * inW + j) * samplesPerPixel;
						//��������� �������������
						int r = static_cast<int>(Mem1[k]); 	     // �������
						int g = static_cast<int>(Mem1[k + 1]);   // �����
						int b = static_cast<int>(Mem1[k + 2]); 	 // �������
						r = kp*(r-mp)+mp;
						g = kp*(g-mp)+mp;
						b = kp*(b-mp)+mp;
						// ����������� �������� ������� ������� � ��������� �� 0 �� 255
						r = std::min(std::max(r, 0), 255);
						g = std::min(std::max(g, 0), 255);
						b = std::min(std::max(b, 0), 255);
						//��������� ��������� �������
						r = static_cast<int>(r * kb); 	 // �������
						g = static_cast<int>(g * kb);    // �����
						b = static_cast<int>(b * kb); 	 // �������
						// ����������� �������� ������� ������� � ��������� �� 0 �� 255
						r = std::min(std::max(r, 0), 255);
						g = std::min(std::max(g, 0), 255);
						b = std::min(std::max(b, 0), 255);
						//���� ������������
						double z11 = r * A;
						double z12 = (static_cast<int>(Mem0[k])) * (1-A);
						r =z11 + z12;
						double z21 = g * A;
						double z22 = (static_cast<int>(Mem0[k+1])) * (1-A);
						g =z21 + z22;
						double z31 = b * A;
						double z32 = (static_cast<int>(Mem0[k+2])) * (1-A);
						b =z31 + z32;
						// ����������� �������� ������� ������� � ��������� �� 0 �� 255
						r = std::min(std::max(r, 0), 255);
						g = std::min(std::max(g, 0), 255);
						b = std::min(std::max(b, 0), 255);
						ptr[j].r = r;
						ptr[j].g = g;
						ptr[j].b = b;
					}
					MoveProgress(i, inH);
				}
				break;
		}

		savebmp->SaveToFile(file);

		savebmp->FreeImage();
		delete savebmp;
	}
}
//---------------------------------------------------------------------------
void __fastcall TForm1::DoChangeMem0ToMem1(){
	Blocking();
	Button1->Caption = "�����...";


	int conv = AnsiString(ComboBox1->ItemIndex).ToInt();
	int Chan = ComboBox2->ItemIndex;
	int post = ComboBox3->ItemIndex;
	converter(conv, Mem0, Mem1, inH, inW, &posK, &posM, post, LastPosition3, LastPosition4, LastPosition5, Chan);
	ProgressBar1->Position = 100;
	Application->ProcessMessages();
	// ������� �������� �� �������� ����������� � ��������
	int brightness = -1;
	if (scale==scaleIns)
		InsMem1toMem2(&brightness);
	else
	{
		int dx=0;
		int dy=0;
		ComputeMoveForScale(&scale, &dx, &dy);
		PutMem1toMem2(scale, dx, dy, &brightness);
	}

	ProgressBar1->Position = 3;
	Application->ProcessMessages();

	LastPosition1 = (brightness*100)/255;
	LastPosition2 = posK;
	LastPosition6 = posM;
	LastPosition7 = 100;
	TrackBar1->Position = (brightness*100)/255;
	TrackBar2->Position = posK;
	TrackBar6->Position = posM;
	TrackBar7->Position = 100;
	posTrackBr = TrackBar1->Position;
	ProgressBar1->Position = 3;
	Application->ProcessMessages();

	PutMem2toMem3(TrackBar1->Position, TrackBar2->Position, TrackBar6->Position, TrackBar7->Position);
	for (int i = h1; i <= h2; i++) // ���������� ����������� � ������ �� �����
		{
			struct rgb *ptr = (struct rgb*)outbmp->ScanLine[i];
			for (int j = w1; j <= w2; j++) {
				ptr[j].r = Mem3[i][j][0];
				ptr[j].g = Mem3[i][j][1];
				ptr[j].b = Mem3[i][j][2];

			}
			MoveProgress(i, h2-h1);
		}

	// ����� ����������� �� �����
	Image1->Picture->Graphic = outbmp;

	WWin->Text = inW;
	HWin->Text = inH;

	ComputeEffectiveBitDepth();

	Button1->Caption = "�������";
	UnBlocking();


}

//---------------------------------------------------------------------------
void __fastcall TForm1::ComboBox1Change(TObject *Sender)
{
	int conv = AnsiString(ComboBox1->ItemIndex).ToInt();
	int nBar = GetBarsNum(conv);
	for(int i=0; i<3; i++) {
		int pos = GetParsPos(conv, i);
		switch (i) {
			case 0:   TrackBar3->Visible=(i<nBar) ? true : false; LastPosition3 = pos; TrackBar3->Position = 100-pos; break;
			case 1:   TrackBar4->Visible=(i<nBar) ? true : false; LastPosition4 = pos; TrackBar4->Position = 100-pos; break;
			default:  TrackBar5->Visible=(i<nBar) ? true : false; LastPosition5 = pos; TrackBar5->Position = 100-pos; break;
		}

	}

	if(block)
			return;

	DoChangeMem0ToMem1();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::TrackBar2Change(TObject *Sender)
{
	float pos2 = (TrackBar2->Position)/3.0; // ������� ������1 �������������
	int pos6 = (TrackBar6->Position)*20; // ������� ������2 �������������
	GroupBox6->Caption = "�������������: " + AnsiString().sprintf("%.2f || %d", pos2, pos6);
}
//---------------------------------------------------------------------------

void __fastcall TForm1::TrackBar3Change(TObject *Sender)
{
	//DoChangeMem0ToMem1();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::TrackBar4Change(TObject *Sender)
{
	//DoChangeMem0ToMem1();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::TrackBar5Change(TObject *Sender)
{
	//DoChangeMem0ToMem1();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::TrackBar6Change(TObject *Sender)
{
	float pos2 = (TrackBar2->Position)/3.0; // ������� ������1 �������������
	int pos6 = (TrackBar6->Position)*20; // ������� ������2 �������������
	GroupBox6->Caption = "�������������: " + AnsiString().sprintf("%.2f || %d", pos2, pos6);
}
//---------------------------------------------------------------------------


void __fastcall TForm1::ComboBox2Change(TObject *Sender)
{
    if(block)
			return;
    DoChangeMem0ToMem1();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button6Click(TObject *Sender)
{


	Blocking();
	ProgressBar1->Position = 3;
	Application->ProcessMessages();
	LastPosition1 = posTrackBr;
	LastPosition2 = posK;
	LastPosition6 = posM;
	LastPosition7 = 100;
	TrackBar1->Position = posTrackBr;
	TrackBar2->Position = posK;
	TrackBar6->Position = posM;
	TrackBar7->Position = 100;
	PutMem2toMem3(LastPosition1, LastPosition2, LastPosition6, LastPosition7);
	for (int i = h1; i <= h2; i++) // ���������� ����������� � ������ �� �����
	{
		struct rgb *ptr = (struct rgb*)outbmp->ScanLine[i];
		for (int j = w1; j <= w2; j++) {
			ptr[j].r = Mem3[i][j][0];
			ptr[j].g = Mem3[i][j][1];
			ptr[j].b = Mem3[i][j][2];
		}
		MoveProgress(i, h2-h1);
	}

	// ����� ����������� �� �����
	Image1->Picture->Graphic = outbmp;

	 //	ComputeEffectiveBitDepth();

	UnBlocking();

}
//---------------------------------------------------------------------------

void __fastcall TForm1::TrackBar7Change(TObject *Sender)
{
	float pos = (TrackBar7->Position)/100.0; // ����� ������� ������ �������
	GroupBox8->Caption = AnsiString("������������: ") + AnsiString().sprintf("%.2f", pos);
}
//---------------------------------------------------------------------------

void __fastcall TForm1::ComboBox3Change(TObject *Sender)
{
	if(block)
			return;
	DoChangeMem0ToMem1();
}
//---------------------------------------------------------------------------

