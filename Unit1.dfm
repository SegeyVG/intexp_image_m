object Form1: TForm1
  Left = 166
  Top = 114
  Width = 1705
  Height = 864
  Caption = 'Form1'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Image1: TImage
    Left = 16
    Top = 16
    Width = 1400
    Height = 780
  end
  object OpenImage1: TButton
    Left = 1472
    Top = 16
    Width = 137
    Height = 25
    Caption = 'Open Image 1'
    TabOrder = 0
    OnClick = OpenImage1Click
  end
  object GroupBox1: TGroupBox
    Left = 1432
    Top = 216
    Width = 201
    Height = 209
    Caption = #1055#1072#1088#1072#1084#1077#1090#1088#1099' '#1080#1079#1086#1073#1088#1072#1078#1077#1085#1080#1103
    TabOrder = 1
    object Label1: TLabel
      Left = 32
      Top = 104
      Width = 138
      Height = 13
      Caption = #1069#1092#1092#1077#1082#1090#1080#1074#1085#1072#1103' '#1088#1072#1079#1088#1103#1076#1085#1086#1089#1090#1100
    end
    object Label2: TLabel
      Left = 16
      Top = 24
      Width = 39
      Height = 13
      Caption = #1064#1080#1088#1080#1085#1072
    end
    object Label3: TLabel
      Left = 16
      Top = 64
      Width = 38
      Height = 13
      Caption = #1042#1099#1089#1086#1090#1072
    end
    object NBitWin: TEdit
      Left = 32
      Top = 136
      Width = 145
      Height = 21
      TabOrder = 0
    end
    object WWin: TEdit
      Left = 64
      Top = 24
      Width = 89
      Height = 21
      TabOrder = 1
    end
    object HWin: TEdit
      Left = 64
      Top = 64
      Width = 89
      Height = 21
      TabOrder = 2
    end
  end
  object OpenPictureDialog1: TOpenPictureDialog
    Left = 1440
    Top = 600
  end
  object SaveDialog1: TSaveDialog
    DefaultExt = 'csv'
    FileName = 'Result.csv'
    Filter = 'Excel|*.csv'
    Left = 1344
    Top = 624
  end
end
