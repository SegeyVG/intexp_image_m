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
#include <Vcl.AppEvnts.hpp>
#include <Vcl.ComCtrls.hpp>
#include <Vcl.Buttons.hpp>
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
	TApplicationEvents *ApplicationEvents1;
	TGroupBox *GroupBox2;
	TComboBox *ComboBox1;
	TButton *Button1;
	TLabel *Label4;
	TTrackBar *TrackBar1;
	TGroupBox *GroupBox3;
	TLabel *Label5;
	TProgressBar *ProgressBar1;
        void __fastcall OpenImage1Click(TObject *Sender);
	void __fastcall ApplicationEvents1Restore(TObject *Sender);
	void __fastcall ApplicationEvents1Minimize(TObject *Sender);
	void __fastcall Button1Click(TObject *Sender);
	void __fastcall TrackBar1Change(TObject *Sender);
private:	// User declarations
public:		// User declarations
		__fastcall TForm1(TComponent* Owner);
		void ComputeEffectiveBitDepth();
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif
