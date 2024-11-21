//---------------------------------------------------------------------------

#include <vcl.h>
#include <stdio.h>
#include <math.h>
#pragma hdrstop

#include "Unit1.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;

// Размеры окна для отображения на форме
#define H 780
#define W 1400

// Максимальные размеры входного изображения
#define maxH 10000
#define maxW 10000

int inH, inW;  // Реальные размеры входного изображения



struct rgb  // структура для доступа к bmp 24 бит
{
  Byte b;
  Byte g;
  Byte r;
//  Byte a;
};
struct rgb32  // // структура для доступа к bmp 32 бит (не используется)
{
  Byte b;
  Byte g;
  Byte r;
  Byte a;
};

int Mem2[H][W][3];   // массив пикселей для вывода на форму
int Mem2Count[H][W]; // счетчик для усреднения пикселей

Byte * Mem1; // указатель на буфер под пиксели из входного файла

  Graphics::TBitmap *inbmp1;
  Graphics::TBitmap *outbmp;

  AnsiString fname;

//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
   inbmp1 = new Graphics::TBitmap;
   outbmp = new Graphics::TBitmap;
   outbmp->Height=H;
   outbmp->Width=W;
   outbmp->PixelFormat=pf24bit;

   long size = maxW;
   size *= maxH;
   size *=3;
   Mem1 = new Byte[size]; // буфер под пиксели из входного файла


}
//---------------------------------------------------------------------------
void PutMem1toMem2()  // прореживание отсчетов с усреднением
{
  int i,j,k,i1,j1;
  for(i=0;i<H;i++) for(j=0;j<W;j++) for(k=0;k<3;k++) {
      Mem2[i][j][k]=0; Mem2Count[i][j]=1;  // обнуление памяти и счетика
  }

  int BSize= inH* inW*3;

    double kX= inW/(double)W;    //коэффициент масштабирования по X
    double kY= inH/(double)H;    //коэффициент масштабирования по Y


   for(int i=0;i<inH;i++)
    {
      for(int j=0;j<inW;j++)
     {
      int k=(i*inW+j)*3; // пересчет трехмерных координат (Y,X,цветовой канал) в одномерные
      i1=i/kY; j1=j/kX;  // позиция пикселя в выходном изображении
      Mem2[i1][j1][0]+=Mem1[k];     // красный
      Mem2[i1][j1][1]+=Mem1[k+1];   // синий
      Mem2[i1][j1][2]+=Mem1[k+2];   // зеленый
      Mem2Count[i1][j1]++;
     }
    }

   for(i=0;i<H;i++) for(j=0;j<W;j++) for(k=0;k<3;k++) {
    if( Mem2Count[i][j]>1)Mem2[i][j][k]/=Mem2Count[i][j]; // нормировка
   }

}

void __fastcall TForm1::OpenImage1Click(TObject *Sender)
{
  if(OpenPictureDialog1->Execute())
  {
   fname= OpenPictureDialog1->FileName;
  FILE *f=fopen(fname.c_str(),"r");
   if(f)
   {
    fclose(f);
    inbmp1->LoadFromFile(fname);  // загрузка изображения из файла

    inH=inbmp1->Height;
    inW=inbmp1->Width;
    for(int i=0;i<inH;i++)
    {
       struct rgb *ptr= (struct rgb*)inbmp1->ScanLine[i];

      for(int j=0;j<inW;j++) //запись пикселей в память для дальнейшей обработки
     {
      int k=(i*inW+j)*3;
      Mem1[k]=ptr[j].r;
      Mem1[k+1]=ptr[j].g;
      Mem1[k+2]=ptr[j].b;
     }
    }


    PutMem1toMem2();  // перенос пикселей из входного изображения в выходное
    for(int i=0;i<H;i++)  // подготовка изображения к выводу на форму
    {
     struct rgb *ptr= (struct rgb*)outbmp->ScanLine[i];
     for(int j=0;j<W;j++)
     {
       ptr[j].r=Mem2[i][j][0];
       ptr[j].g=Mem2[i][j][1];
       ptr[j].b=Mem2[i][j][2];

     }
    }
// Вывод изображения на форму
    Image1->Picture->Graphic=outbmp;

    WWin->Text=inW;
    HWin->Text=inH;

// Расчет показателей эффективной разрядности
    int BSize= inH* inW*3;
    int Gist[256]; int i;
    for(i=0;i<256;i++) Gist[i]=0;
    for(i=0;i<BSize;i++) Gist[Mem1[i]]++;
    double h1=0;
    for(i=0;i<256;i++)
    {
      if(Gist[i]==0) continue;
      double N=BSize;
      double p= Gist[i]/N;
      h1-= p*log(p)/log(2);
    }
    double n=0;
    for(i=0;i<256;i++) if(Gist[i]>0) { Gist[i]=1; n+=1;}

    double h2=0;
    for(i=0;i<256;i++)
    {
      if(Gist[i]==0) continue;

      double p= Gist[i]/n;
      h2-= p*log(p)/log(2);
    }
    char buf[64]; sprintf(buf,"%.2f (%.2f)",h2,h1);
    NBitWin->Text=buf;

    }
  }
}









