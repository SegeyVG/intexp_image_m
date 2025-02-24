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
#include <Vcl.Samples.Spin.hpp>
#include <System.ImageList.hpp>
#include <Vcl.ImgList.hpp>
#include <System.Win.ScktComp.hpp>

void MoveProgress(int i, int inH);

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
	TProgressBar *ProgressBar1;
	TTimer *Timer1;
	TGroupBox *GroupBox4;
	TImageList *ImageList1;
	TButton *Button2;
	TButton *Button3;
	TButton *Button4;
	TButton *Button5;
	TSaveDialog *SavePictureDialog1;
	TGroupBox *GroupBox5;
	TRadioGroup *RadioGroup1;
	TButton *SaveImg;
	TGroupBox *GroupBox6;
	TTrackBar *TrackBar2;
	TGroupBox *GroupBox7;
	TTrackBar *TrackBar3;
	TTrackBar *TrackBar4;
	TTrackBar *TrackBar5;
	TTrackBar *TrackBar6;
	TComboBox *ComboBox2;
	TButton *Button6;
	TGroupBox *GroupBox8;
	TTrackBar *TrackBar7;
	TGroupBox *GroupBox9;
	TComboBox *ComboBox3;
        void __fastcall OpenImage1Click(TObject *Sender);
	void __fastcall ApplicationEvents1Restore(TObject *Sender);
	void __fastcall ApplicationEvents1Minimize(TObject *Sender);
	void __fastcall Button1Click(TObject *Sender);
	void __fastcall TrackBar1Change(TObject *Sender);
	void __fastcall Timer1Timer(TObject *Sender);
	void __fastcall Button4Click(TObject *Sender);
	void __fastcall Button5Click(TObject *Sender);
	void __fastcall Button2Click(TObject *Sender);
	void __fastcall Button3Click(TObject *Sender);
	void __fastcall Image1MouseMove(TObject *Sender, TShiftState Shift, int X, int Y);
	void __fastcall Image1MouseLeave(TObject *Sender);
	void __fastcall Image1Click(TObject *Sender);
	void __fastcall SaveImgClick(TObject *Sender);
	void __fastcall ComboBox1Change(TObject *Sender);
	void __fastcall TrackBar2Change(TObject *Sender);
	void __fastcall TrackBar3Change(TObject *Sender);
	void __fastcall TrackBar4Change(TObject *Sender);
	void __fastcall TrackBar5Change(TObject *Sender);
	void __fastcall TrackBar6Change(TObject *Sender);
	void __fastcall ComboBox2Change(TObject *Sender);
	void __fastcall Button6Click(TObject *Sender);
	void __fastcall TrackBar7Change(TObject *Sender);
	void __fastcall ComboBox3Change(TObject *Sender);




private:	// User declarations
public:		// User declarations
		int LastPosition1, LastPosition2, LastPosition6, LastPosition7;
		int LastPosition3, LastPosition4, LastPosition5;
		__fastcall TForm1(TComponent* Owner);
		void __fastcall DoChangeMem0ToMem1();
		void ButtonsScaleConfiguration();
        void SetTpCursor();
		void Blocking();
		void UnBlocking();
		void ComputeEffectiveBitDepth();
		int posTrackBr; //позиция ползунка TrackBar1, соответствующая изначальной средней яркости картинки
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif
