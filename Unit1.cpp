// ---------------------------------------------------------------------------

#include <vcl.h>
#include <math.h>
#include <cstring>
#include <vector>
#include <algorithm>
#pragma hdrstop

#include "Unit1.h"
#include "Unit2.h"
#include "BinFile.h"

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
int samplesPerPixel = 3;  //���������� �������� �������
long sizeInBuff;

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

Byte * Mem0; // ��������� �� ��������� ����� ��� ������� �� �������� �����
Byte * Mem01; // ��������� �� ��������� ����� ��� ������� �� �������� ����� (� ������ ������� ����������� 16 ���/�����, ���� �������� ������� �����)
Byte * Mem1; // ��������� �� ��������� ����� ��� ������� �� �������� ����� (����� ��������������� ��������������, ���� ��� ����)


int*** Mem2;  		// ������ �������� ��� ������ �� �����
int**  Mem2Count;  	// ������� ���������� ���������� ����������� �������
int*** Mem3;  		// ������ �������� ��� ������ �� �����
int w1, w2, h1, h2; //������� ������, �������������� �������� ����� ��������������� � ���������� �� �����
					// (w - �� ���������, h - �� �����������)

Graphics::TBitmap *inbmp;
Graphics::TBitmap *outbmp;

AnsiString fname;
bool flagA16b;

// ---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner) : TForm(Owner) {
	//��������� ����� � ����� ������ ����������� ��� ������ ������ ��� ��������
	int FW = 0.9*SW;
	int FH = 0.9*SH;
	Form1->Top = 0.05*SH;
	Form1->Left = 0.05*SW;
	Image1->Width = 0.865*FW;
	Image1->Height = 0.96*FH;
	ProgressBar1->Left = 0.87*FW;
	OpenImage1->Left = 0.88*FW;
	GroupBox1->Left = 0.87*FW;
	GroupBox2->Left = 0.87*FW;
	GroupBox3->Left = 0.87*FW;

	OpenImage1->Top = 0.03*0.96*FH;
	GroupBox1->Top = 0.03*0.96*FH +(OpenImage1->Height)+15;
	GroupBox2->Top = 0.03*0.96*FH +(OpenImage1->Height)+15+(GroupBox1->Height)+12;
	GroupBox3->Top = 0.03*0.96*FH +(OpenImage1->Height)+15+(GroupBox1->Height)+12+(GroupBox2->Height)+12;
	ProgressBar1->Top = 0.96*FH - (ProgressBar1->Height);

	Label4->Caption = AnsiString("�����: ") + SW + "x" + SH;
	LastPosition = TrackBar1->Position;
    GroupBox3->Enabled = false;

	H = Image1->Height;
	W = Image1->Width;
	inbmp = new Graphics::TBitmap;
	outbmp = new Graphics::TBitmap;
	outbmp->Height = H;
	outbmp->Width = W;
	outbmp->PixelFormat = pf24bit;
	//��������� ������ ��� ������ �������� �����������
	long sizeInBuff = MaxW*MaxH*samplesPerPixel;
	Mem0 = new Byte[sizeInBuff]; // ��������� ����� ��� ������� �� �������� �����
	Mem1 = new Byte[sizeInBuff]; // ��������� ����� ��� ������� �� �������� �����

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
	Mem3 = new int**[H];
	for(int i = 0; i < H; i++) {
		Mem2[i] = new int*[W];
		Mem3[i] = new int*[W];
		for(int j = 0; j < W; j++) {
			Mem2[i][j] = new int[3];
            Mem3[i][j] = new int[3];
		}
	}

	Mem2Count = new int*[H];
	for(int i = 0; i < H; i++) {
		Mem2Count[i] = new int[W];
	}
}
//----------------------------------------------------------------------------
void BeeLineInterpolation(int x, int y) {
	//����� ������� ��������� �������
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
	//������� �� ���������� �������� ������
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
	//����� ��������� �������
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
	//������� �� ���������� �������� ������
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
	//����� ��������� �������
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
    //������� �� ���������� �������� ������
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

// ---------------------------------------------------------------------------
void PutMem1toMem2(int* br)
	// ���������� (���������������) ������ �������� ����������� �� �������� ������ � ����� ��� ������ �� �����
{
	int i, j, k, i1, j1;
	for (i = 0; i < H; i++)
		for (j = 0; j < W; j++)
			for (k = 0; k < samplesPerPixel; k++) {
				Mem2[i][j][k] = 0;
				Mem2Count[i][j] = 0; // ��������� ������ � �������
			}

	double kX = inW / (double)W; // ����������� ��������������� �� X
	double kY = inH / (double)H; // ����������� ��������������� �� Y

    // ����������, ����� ����������� ����� ����� � ����������� �� ����������� ������ ����������� � �����
	double scale;
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
		}
		MoveProgress(i, inH);
	}

	w1 = movW;  w2 = (inW-1)/scale + movW;
	h1 = movH;  h2 = (inH-1)/scale + movH;
	if (inH > H || inW > W)
		// ���������� ������������ �������� (� ������ ���� ���� �� ���� �� ������ �����������
		// ������ ��������������� ������� ���� ������)
		for (i = h1; i <= h2; i++)  {
			for (j = w1; j <= w2; j++)
				for (k = 0; k < samplesPerPixel; k++) {
					if (Mem2Count[i][j] > 1)
						Mem2[i][j][k] /= Mem2Count[i][j];
				}
			MoveProgress(i, h2-h1);
		}

	if (inH < H || inW < W) {
		// ������������ (� ������ ���� ���� �� ���� �� ������ ����������� ������ ��������������� ������� ���� ������)
		for (int j = w1; j <= w2; j++) {
			if (Mem2Count[h1][j] == 0)
				LineInterpolationX(j, h1);
			if (Mem2Count[h2][j] == 0)
				LineInterpolationX(j, h2);
		}
		for (int i = h1; i <= h2; i++) {
			if (Mem2Count[i][w1] == 0)
				LineInterpolationY(w1, i);
			if (Mem2Count[i][w2] == 0)
				LineInterpolationY(w2, i);
		}
		for (int i = h1+1; i < h2; i++)
			for (int j = w1+1; j < w2; j++)
				if (Mem2Count[i][j] == 0)
					BeeLineInterpolation(j, i);
	}

	if (br==nullptr)
		return;
	//������ ������� ������� ����������� �� ������, ���������� ��� ������ Mem1===================================================
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
//----------------------------------------------------------------------------------------------------------
void PutMem2toMem3(double kb){
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			for (int k = 0; k < samplesPerPixel; k++)
				Mem3[i][j][k] = 0;	// ��������� ������ � �������

	for (int i = h1; i <= h2; i++)
	{
		for (int j = w1; j <= w2; j++) {
			int r = static_cast<int>(Mem2[i][j][0] * kb);
			int g = static_cast<int>(Mem2[i][j][1] * kb);
			int b = static_cast<int>(Mem2[i][j][2] * kb);
            // ����������� �������� ������� ������� � ��������� �� 0 �� 255
			r = std::min(std::max(r, 0), 255);
			g = std::min(std::max(g, 0), 255);
			b = std::min(std::max(b, 0), 255);
			Mem3[i][j][0] = r;
			Mem3[i][j][1] = g;
			Mem3[i][j][2] = b;
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
	GroupBox2->Enabled = false;
	GroupBox3->Enabled = false;
	ProgressBar1->Position = 3;
	ProgressBar1->Visible = true;
	Application->ProcessMessages();
	OpenImage1->Enabled = false;
}
//------------------------------------------------------------------------------
void TForm1::UnBlocking(){
	GroupBox2->Enabled = true;
	GroupBox3->Enabled = true;
	ProgressBar1->Visible = false;
	OpenImage1->Enabled = true;

	Image1->Width = 0.865*0.9*SW;
	Image1->Height = 0.96*0.9*SH;
	AutoSize = true;
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
/*
������� ExtractFileExt ��������� �� ������� ����� ����� ���������� �����.
AnsiLowerCase() - ��������� ������ � ������ �������.
*/
void __fastcall TForm1::OpenImage1Click(TObject *Sender) {
	if (OpenPictureDialog1->Execute()) {
		fname = OpenPictureDialog1->FileName;
		FILE *f = fopen(fname.c_str(), "r");
		String fileExt = AnsiLowerCase(ExtractFileExt(fname));
		Blocking();
		OpenImage1->Caption = "��������� ����...";
		int i1, j1, brPixel;
		if (fileExt == ".bmp") {
			fclose(f);
			inbmp->LoadFromFile(fname); // �������� ����������� �� �����
			ProgressBar1->Position = 100;
			Application->ProcessMessages();
			inH = inbmp->Height;
			inW = inbmp->Width;
			long size = inW*inH*samplesPerPixel;
			if (size>sizeInBuff) {
				delete [] Mem0;
				delete [] Mem1;
				Mem0 = new Byte[size]; // ��������� ����� ��� ������� �� �������� �����
				Mem1 = new Byte[size]; // ��������� ����� ��� ������� �� �������� �����
			}
			if (Mem01 != nullptr) {
				delete [] Mem01;
				Mem01 = nullptr;
			}
			flagA16b = false;

			for (int i = 0; i < inH; i++) {
				struct rgb *ptr = (struct rgb*)inbmp->ScanLine[i];

				for (int j = 0; j < inW; j++) // ������ �������� � ��������� ����� ��� ���������� ���������
				{
					int k = (i * inW + j) * samplesPerPixel;
					Mem0[k] = ptr[j].r;
					Mem0[k + 1] = ptr[j].g;
					Mem0[k + 2] = ptr[j].b;
				}
				MoveProgress(i, inH);
			}
			ProgressBar1->Position = 100;
			Application->ProcessMessages();
		}
		else if (fileExt == ".bin") {
			try
			{
				std::string Filename = fname.c_str();
				ImageData imgData = binaryToArray(Filename);
				ProgressBar1->Position = 100;
				Application->ProcessMessages();
				inH = imgData.height;
				inW = imgData.width;
				int channels = imgData.samplesPerPixel;
				long size = inW*inH*samplesPerPixel;
				if (size>sizeInBuff) {
					delete [] Mem0;
					delete [] Mem1;
					Mem0 = new Byte[size]; // ��������� ����� ��� ������� �� �������� �����
					Mem1 = new Byte[size]; // ��������� ����� ��� ������� �� �������� �����
				}

				if (Mem01 != nullptr)
					delete [] Mem01;
				Mem01 = new Byte[size];
				flagA16b = true;

				for (int i = 0; i < inH; i++) {
					for (int j = 0; j < inW; j++) {
						int k1 = (i * inW + j) * samplesPerPixel;
						for (int c = 0; c < 3; c++)  {//
								int k = (i * inW + j) * samplesPerPixel + c;
								//Byte uc = (imgData.at(i,j,c)) / 256;
								//Mem0[k] = uc;
								unsigned short v16 = imgData.at(i,j,c);
								// ���������� ������� ���� ����� � Mem0[k]
								Mem0[k] = (v16 >> 8) & 0xFF;
								// ���������� ������� ���� ����� � Mem01[k]
								Mem01[k] = v16 & 0xFF;
						}
						i1=i; j1=j;
					}
					MoveProgress(i, inH);
				}
				ProgressBar1->Position = 100;
				Application->ProcessMessages();
			}
			catch (const std::exception &e)
			{
				ShowMessage(AnsiString("������ ������ ��� ��������� ������ �����������. ������� c �����.(x:y)-")+j1+ ":"+i1);
			}
		}else {
			ShowMessage("���������������� ��� ������");
			return;
		}

		int conv = AnsiString(ComboBox1->ItemIndex).ToInt();
		converter(conv, Mem0, Mem1, inH, inW);
		ProgressBar1->Position = 100;
		Application->ProcessMessages();
		// ������� �������� �� �������� ����������� � ��������
		int brightness = -1;
		PutMem1toMem2(&brightness);
		TrackBar1->Position = (brightness*100)/255;
		posTrackBr = TrackBar1->Position;
		ProgressBar1->Position = 3;
		Application->ProcessMessages();
		PutMem2toMem3(1);

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
	}
}

//----------------------------------------------------------------------------
void __fastcall TForm1::ApplicationEvents1Restore(TObject *Sender)
{
	//��������� ����� � ����� ������ ����������� ��� ������ ������ ��� �������������� ����� �� ���������� ���������
	int FW = 0.9*SW;
	int FH = 0.9*SH;
	Form1->Top = 0.05*SH;
	Form1->Left = 0.05*SW;
	Image1->Width = 0.865*FW;
	Image1->Height = 0.96*FH;
	ProgressBar1->Left = 0.87*FW;
	OpenImage1->Left = 0.88*FW;
	GroupBox1->Left = 0.87*FW;
	GroupBox2->Left = 0.87*FW;
	GroupBox3->Left = 0.87*FW;

	OpenImage1->Top = 0.03*0.96*FH;
	GroupBox1->Top = 0.03*0.96*FH +(OpenImage1->Height)+15;
	GroupBox2->Top = 0.03*0.96*FH +(OpenImage1->Height)+15+(GroupBox1->Height)+12;
	GroupBox3->Top = 0.03*0.96*FH +(OpenImage1->Height)+15+(GroupBox1->Height)+12+(GroupBox2->Height)+12;

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


	int conv = AnsiString(ComboBox1->ItemIndex).ToInt();
	converter(conv, Mem0, Mem1, inH, inW);
	ProgressBar1->Position = 100;
	Application->ProcessMessages();
	// ������� �������� �� �������� ����������� � ��������
	int brightness = -1;
	PutMem1toMem2(&brightness);
	TrackBar1->Position = (brightness*100)/255;
	posTrackBr = TrackBar1->Position;
	ProgressBar1->Position = 3;
	Application->ProcessMessages();
	PutMem2toMem3(1);
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

	Button1->Caption = "���������";
	UnBlocking();
}

//--------------------------------------------------------------------------
void __fastcall TForm1::TrackBar1Change(TObject *Sender)
{
	int pos2 = TrackBar1->Position; // ����� ������� ������ �������
	GroupBox3->Caption = AnsiString("�������: ") + pos2;
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Timer1Timer(TObject *Sender)
{
	if (TrackBar1->Position != LastPosition)
	{
		LastPosition = TrackBar1->Position;

		int pos2 = TrackBar1->Position; // ����� ������� ������ �������
		GroupBox3->Caption = AnsiString("�������: ") + pos2;

		int pos1 = posTrackBr;
		double kb = static_cast<double>(pos2) / pos1; // ����������� ��������� �������

		Blocking();

		ProgressBar1->Position = 3;
		Application->ProcessMessages();
		PutMem2toMem3(kb);
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

		UnBlocking();

	}
}
//---------------------------------------------------------------------------

