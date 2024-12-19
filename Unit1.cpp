// ---------------------------------------------------------------------------

#include <vcl.h>
#include <math.h>
#include <cstring>
#include <vector>
#pragma hdrstop

#include "Unit1.h"
#include "Unit2.h"
#include "BinFile.h"
//#include "TestReadTiff.h"
// ---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"

using namespace std;

TForm1 *Form1;
//максимальные размеры входного изображения (по умолчанию)
#define MaxH 10000
#define MaxW 10000

int H, W;     // Размеры окна для отображения на форме
int inH, inW; // Реальные размеры входного изображения
int samplesPerPixel = 3;  //количество цветовых каналов

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

Byte * Mem0; // указатель на первичный буфер под пиксели из входного файла
Byte * Mem01; // указатель на первичный буфер под пиксели из входного файла (в случае прихода изображения 16 бит/канал, сюда кладутся младшие байты)
Byte * Mem1; // указатель на вторичный буфер под пиксели из входного файла (после дополнительного преобразования, если оно есть)


int*** Mem2;  		// массив пикселей для вывода на форму
int**  Mem2Count;  	// счетчик количества прорисовок конкретного пиксела
int w1, w2, h1, h2; //границы шторок, ограничивающих картинку после масштабирования и размещения на форме
					// (w - по вертикали, h - по горизонтали)

Graphics::TBitmap *inbmp;
Graphics::TBitmap *outbmp;

AnsiString fname;
bool flagA16b;

// ---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner) : TForm(Owner) {
	//адаптация формы и листа вывода изображения под размер экрана при открытии
	int FW = 0.9*SW;
	int FH = 0.9*SH;
	Form1->Top = 0.05*SH;
	Form1->Left = 0.05*SW;
	Image1->Width = 0.86*FW;
	Image1->Height = 0.96*FH;
	OpenImage1->Left = 0.88*FW;
	GroupBox1->Left = 0.87*FW;
	GroupBox2->Left = 0.87*FW;
	GroupBox3->Left = 0.87*FW;

	OpenImage1->Top = 0.02*0.96*FH;
	GroupBox1->Top = 0.02*0.96*FH +(OpenImage1->Height)+21;
	GroupBox2->Top = 0.02*0.96*FH +(OpenImage1->Height)+21+(GroupBox1->Height)+12;
	GroupBox3->Top = 0.02*0.96*FH +(OpenImage1->Height)+21+(GroupBox1->Height)+12+(GroupBox2->Height)+12;

	Label4->Caption = AnsiString("Экран: ") + SW + "x" + SH;
	//TrackBar1->SetTick(20);

	H = Image1->Height;
	W = Image1->Width;
	inbmp = new Graphics::TBitmap;
	outbmp = new Graphics::TBitmap;
	outbmp->Height = H;
	outbmp->Width = W;
	outbmp->PixelFormat = pf24bit;
	//выделение памяти под буферы входного изображения
	long size = MaxW*MaxH*samplesPerPixel;
    Mem0 = new Byte[size]; // первичный буфер под пиксели из входного файла
	Mem1 = new Byte[size]; // вторичный буфер под пиксели из входного файла

	//красивый фон (как в Paint) для листа вывода изображений, чтобы заодно было его сразу видно на форме визуально
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

	//выделение памяти под матрицы для вывода изображения на форму. Делается динамически(поскольку размер экрана может быть любой), но 1 раз
	Mem2 = new int**[H];
	for(int i = 0; i < H; i++) {
		Mem2[i] = new int*[W];
		for(int j = 0; j < W; j++) {
			Mem2[i][j] = new int[3];
		}
	}

	Mem2Count = new int*[H];
	for(int i = 0; i < H; i++) {
		Mem2Count[i] = new int[W];
	}
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

//----------------------------------------------------------------------------
void BackgroundFill(int x1, int y1, int x2, int y2){
	for (int i = y1; i < y2; i++)
		for (int j = x1; j < x2; j++) {
			Mem2[i][j][0] = 203; // красный
			Mem2[i][j][1] = 213; // синий
			Mem2[i][j][2] = 228; // зеленый
			Mem2Count[i][j]++;
		}

}

// ---------------------------------------------------------------------------
void PutMem1toMem2()
	// трансляция (масштабирование) данных пикселов изображения из входного буфера в буфер для вывода на форму
{
	int i, j, k, i1, j1;
	for (i = 0; i < H; i++)
		for (j = 0; j < W; j++)
			for (k = 0; k < samplesPerPixel; k++) {
				Mem2[i][j][k] = 0;
				Mem2Count[i][j] = 0; // обнуление памяти и счетика
			}

	double kX = inW / (double)W; // коэффициент масштабирования по X
	double kY = inH / (double)H; // коэффициент масштабирования по Y

    // Определяем, какой коэффициент нужно взять в зависимости от соотношения сторон изображения и рамки
	double scale;
	scale = (kX > kY) ? kX : kY;

	// Применяем коэффициент масштабирования к обоим сторонам изображения
	int imgW = static_cast<int>(inW / scale); // новая ширина изображения
	int imgH = static_cast<int>(inH / scale); // новая высота изображения
	int movW = (W-imgW)/2;
	int movH = (H-imgH)/2;

	for (int i = 0; i < inH; i++) {
		for (int j = 0; j < inW; j++) {
			int k = (i * inW + j) * samplesPerPixel;
			// пересчет трехмерных координат (Y,X,цветовой канал) в одномерные
			i1 = i / scale + movH;  // позиция пикселя в столбце выходного изображения
			j1 = j / scale + movW; // позиция пикселя в строке выходного изображения
			Mem2[i1][j1][0] += Mem1[k]; // красный
			Mem2[i1][j1][1] += Mem1[k + 1]; // синий
			Mem2[i1][j1][2] += Mem1[k + 2]; // зеленый
			Mem2Count[i1][j1]++;
		}

	}

	w1 = movW;  w2 = (inW-1)/scale + movW;
	h1 = movH;  h2 = (inH-1)/scale + movH;
	if (inH > H || inW > W)
		// усреднение получившихся пикселей (в случае если хотя бы одна из сторон изображения
		// больше соответствующей стороны окна вывода)
		for (i = h1; i <= h2; i++)
			for (j = w1; j <= w2; j++)
				for (k = 0; k < samplesPerPixel; k++) {
					if (Mem2Count[i][j] > 1)
						Mem2[i][j][k] /= Mem2Count[i][j];
				}

//	int pW = (inW - 1) / scale;
//	int pH = (inH - 1) / scale;
	if (inH < H || inW < W) {
		// интерполяция (в случае если хотя бы одна из сторон изображения меньше соответствующей стороны окна вывода)
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

	//заливаем получившиеся шторки цветом фона
	if (movW>0) {
		BackgroundFill(0, 0, w1, H);
		BackgroundFill(w2+1, 0, W, H);
	}

	if (movH>0) {
		BackgroundFill(0, 0, W, h1);
		BackgroundFill(0, h2+1, W, H);
	}

}

void  TForm1::ComputeEffectiveBitDepth() {     // Расчет показателей эффективной разрядности
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

/*
Функция ExtractFileExt извлекает из полного имени файла расширение файла.
AnsiLowerCase() - переводит строку в нижний регистр.

*/
void __fastcall TForm1::OpenImage1Click(TObject *Sender) {
	if (OpenPictureDialog1->Execute()) {
		fname = OpenPictureDialog1->FileName;
		FILE *f = fopen(fname.c_str(), "r");
		String fileExt = AnsiLowerCase(ExtractFileExt(fname));
		vector<int> bPixels;
		int i1, j1, brPixel;
		if (fileExt == ".bmp") {
			fclose(f);
			inbmp->LoadFromFile(fname); // загрузка изображения из файла
			inH = inbmp->Height;
			inW = inbmp->Width;
			bPixels.reserve(inH*inW);
			long size = inW*inH*samplesPerPixel;
			if ((inW>MaxW) || (inH>MaxH)) {
				delete [] Mem0;
				delete [] Mem1;
				Mem0 = new Byte[size]; // первичный буфер под пиксели из входного файла
				Mem1 = new Byte[size]; // вторичный буфер под пиксели из входного файла
			}
			if (Mem01 != nullptr) {
				delete [] Mem01;
				Mem01 = nullptr;
			}
			flagA16b = false;

			for (int i = 0; i < inH; i++) {
				struct rgb *ptr = (struct rgb*)inbmp->ScanLine[i];

				for (int j = 0; j < inW; j++) // запись пикселей в первичный буфер для дальнейшей обработки
				{
					int k = (i * inW + j) * samplesPerPixel;
					Mem0[k] = ptr[j].r;
					Mem0[k + 1] = ptr[j].g;
					Mem0[k + 2] = ptr[j].b;
					int brightnessPixel = (int)(0.2989 * ptr[j].r + 0.5870 * ptr[j].g + 0.1140 * ptr[j].b);
					bPixels.push_back(brightnessPixel);
				}
			}
		}
		else if (fileExt == ".bin") {
			try
			{
				std::string Filename = fname.c_str();
				ImageData imgData = binaryToArray(Filename);
				inH = imgData.height;
				inW = imgData.width;
				bPixels.reserve(inH*inW);
				int channels = imgData.samplesPerPixel;
				long size = inW*inH*samplesPerPixel;
				if ((inW>MaxW) || (inH>MaxH)) {
					delete [] Mem0;
					delete [] Mem1;
					Mem0 = new Byte[size]; // первичный буфер под пиксели из входного файла
					Mem1 = new Byte[size]; // вторичный буфер под пиксели из входного файла
				}
				if (Mem01 != nullptr)
					delete [] Mem01;
				Mem01 = new Byte[size];
				flagA16b = true;

				for (int i = 0; i < inH; i++) {
					for (int j = 0; j < inW; j++) {
						if (i==4246 && j==6481) {
                            ShowMessage("Alarm");
						}
						int k1 = (i * inW + j) * samplesPerPixel;
						for (int c = 0; c < 3; c++)  {//
								int k = (i * inW + j) * samplesPerPixel + c;
								//Byte uc = (imgData.at(i,j,c)) / 256;
								//Mem0[k] = uc;
								unsigned short v16 = imgData.at(i,j,c);
								// записываем старший байт числа в Mem0[k]
								Mem0[k] = (v16 >> 8) & 0xFF;
								// записываем младший байт числа в Mem01[k]
								Mem01[k] = v16 & 0xFF;
						}
						i1=i; j1=j;
						int brightnessPixel = (int)(0.2989 * Mem0[k1] + 0.5870 * Mem0[k1+1] + 0.1140 * Mem0[k1+2]);
						brPixel = brightnessPixel;
						bPixels.push_back(brightnessPixel);
					}
				}
			}
			catch (const std::exception &e)
			{
				ShowMessage(AnsiString("Ошибка чтения или обработки данных изображения. Пиксель c коорд.(x:y)-")+j1+ ":"+i1+", с яркостью: "+brPixel);
			}
		}else {
			ShowMessage("Неподдерживаемый тип файлов");
			return;
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


		TrackBar1->Position = (brightness*100)/255;
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

//----------------------------------------------------------------------------
void __fastcall TForm1::ApplicationEvents1Restore(TObject *Sender)
{
	//адаптация формы и листа вывода изображения под размер экрана при разворачивании формы из свернутого состояния
	int FW = 0.9*SW;
	int FH = 0.9*SH;
	Form1->Top = 0.05*SH;
	Form1->Left = 0.05*SW;
	Image1->Width = 0.86*FW;
	Image1->Height = 0.96*FH;
	OpenImage1->Left = 0.88*FW;
	GroupBox1->Left = 0.87*FW;
	GroupBox2->Left = 0.87*FW;

	OpenImage1->Top = 0.02*0.96*FH;
	GroupBox1->Top = 0.02*0.96*FH +(OpenImage1->Height)+21;
	GroupBox2->Top = 0.02*0.96*FH +(OpenImage1->Height)+21+(GroupBox1->Height)+12;
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


void __fastcall TForm1::TrackBar1Change(TObject *Sender)
{
	int value = TrackBar1->Position;
	Label5->Caption = AnsiString("Значение: ") + value +"%";
}
//---------------------------------------------------------------------------

