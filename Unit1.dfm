object Form1: TForm1
  Left = 166
  Top = 114
  AutoSize = True
  BorderIcons = [biSystemMenu, biMinimize]
  BorderStyle = bsSingle
  Caption = #1054#1073#1088#1072#1073#1086#1090#1082#1072' '#1080#1079#1086#1073#1088#1072#1078#1077#1085#1080#1081
  ClientHeight = 835
  ClientWidth = 1688
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  Padding.Left = 1
  Padding.Top = 1
  Padding.Right = 6
  Padding.Bottom = 1
  OldCreateOrder = False
  DesignSize = (
    1688
    835)
  PixelsPerInch = 96
  TextHeight = 13
  object Image1: TImage
    Left = 1
    Top = 1
    Width = 1474
    Height = 833
    Anchors = [akLeft, akTop, akBottom]
  end
  object GroupBox2: TGroupBox
    Left = 1481
    Top = 385
    Width = 201
    Height = 105
    Caption = #1060#1080#1083#1100#1090#1088
    TabOrder = 2
    object ComboBox1: TComboBox
      Left = 32
      Top = 24
      Width = 145
      Height = 21
      ItemIndex = 0
      TabOrder = 0
      Text = #1050#1072#1082' '#1077#1089#1090#1100
      Items.Strings = (
        #1050#1072#1082' '#1077#1089#1090#1100
        #1063#1077#1088#1085#1086'-'#1073#1077#1083#1086#1077)
    end
    object Button1: TButton
      Left = 48
      Top = 64
      Width = 113
      Height = 25
      Caption = #1055#1088#1080#1084#1077#1085#1080#1090#1100
      TabOrder = 1
      OnClick = Button1Click
    end
  end
  object OpenImage1: TButton
    Left = 1517
    Top = 113
    Width = 137
    Height = 25
    Caption = #1054#1090#1082#1088#1099#1090#1100
    TabOrder = 0
    OnClick = OpenImage1Click
  end
  object GroupBox1: TGroupBox
    Left = 1481
    Top = 161
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
      Left = 61
      Top = 24
      Width = 89
      Height = 21
      TabOrder = 1
    end
    object HWin: TEdit
      Left = 60
      Top = 64
      Width = 89
      Height = 21
      TabOrder = 2
    end
  end
  object OpenPictureDialog1: TOpenPictureDialog
    Left = 1312
    Top = 712
  end
  object SaveDialog1: TSaveDialog
    DefaultExt = 'csv'
    FileName = 'Result.csv'
    Filter = 'Excel|*.csv'
    Left = 1168
    Top = 672
  end
  object ApplicationEvents1: TApplicationEvents
    OnMinimize = ApplicationEvents1Minimize
    OnRestore = ApplicationEvents1Restore
    Left = 1184
    Top = 480
  end
end
