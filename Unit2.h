//---------------------------------------------------------------------------

#ifndef Unit2H
#define Unit2H

void converter(int, Byte *, Byte *, int, int, int*, int*, int, int Par1, int Par2, int Par3,int Chan);
void PostProcess(int , Byte *, Byte *, int, int, int Par1, int Par2, int Par3);
char *GetNames(int conv); // запросить имена функций и отправить их в список
char* GetPostNames(int postIndex); // запросить имена постобработок и отправить их в список
int GetBarsNum(int conv); // после выбора конкретной функции обработки запрсить количество параметров, ненужные скрыть
int GetParsPos(int conv, int Bar); // дл€ каждого активного параметра запросить позицию "ползункаЁ

void EntrA(int chan, Byte * Mem0, Byte * Mem1, int H, int W);
void LinA(int chan, Byte * Mem0, Byte * Mem1, int H, int W, int Par1);
void EntrB(int chan, Byte * Mem0, Byte * Mem1, int H, int W);
void LinB(int chan, Byte * Mem0, Byte * Mem1, int H, int W);

void Entr8R( Byte * Mem0, Byte * Mem1, int H, int W);
void Entr8G( Byte * Mem0, Byte * Mem1, int H, int W);
void Entr8B( Byte * Mem0, Byte * Mem1, int H, int W);
void Entr8(int chan, Byte * Mem0, Byte * Mem1, int H, int W);
void Linear8(int chan, Byte * Mem0, Byte * Mem1, int H, int W);
void TestFilter( Byte * Mem0, Byte * Mem1, int H, int W);
void TestFilter2( Byte * Mem0, Byte * Mem1, int H, int W);
void Entr8Q(int chan, Byte * Mem0, Byte * Mem1, int H, int W);
void Entr8EQ(int chan, Byte * Mem0, Byte * Mem1, int H, int W);
//---------------------------------------------------------------------------
#endif

