//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "Unit2.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

void converter(int conv, Byte * Mem0, Byte * Mem1, int H, int W)
{
	switch (conv) {
	case 1:
        for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++)
			{
				int k = (i * W + j) * 3;
				int wb = (Mem0[k] + Mem0[k + 1] + Mem0[k + 2]) / 3;
				Mem1[k] = wb;
				Mem1[k + 1] = wb;
				Mem1[k + 2] = wb;
			}
		}

		break;
	default:
		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++)
			{
				int k = (i * W + j) * 3;
				Mem1[k] = Mem0[k];
				Mem1[k + 1] = Mem0[k + 1];
				Mem1[k + 2] = Mem0[k + 2];
			}
		}
		break;
	}
}
