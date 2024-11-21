//---------------------------------------------------------------------------

#ifndef Unit1H
#define Unit1H
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <Dialogs.hpp>
#include <ExtDlgs.hpp>
//---------------------------------------------------------------------------
class TForm1 : public TForm
{
__published:	// IDE-managed Components
        TImage *Image1;
        TButton *OpenImage1;
        TOpenPictureDialog *OpenPictureDialog1;
        TSaveDialog *SaveDialog1;
        TGroupBox *GroupBox1;
        TEdit *NBitWin;
        TLabel *Label1;
        TLabel *Label2;
        TEdit *WWin;
        TLabel *Label3;
        TEdit *HWin;
        void __fastcall OpenImage1Click(TObject *Sender);
private:	// User declarations
public:		// User declarations
        __fastcall TForm1(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif
