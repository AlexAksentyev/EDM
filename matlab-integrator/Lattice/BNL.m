clc; clear;

initBDT; setD;

h = .05;
R = 9.207;
Lq = 5;
Lw = 1.808;

GSFP = 0;
GSDP = 0;

DL_25 = @()emdrift(.25);
DL_15 = @()emdrift(.15);
DL2_2 = @()emdrift(2.2);
BPM = @()emdrift(15);

QDS = @()mquadr(Lq, -.86);
QFS = @()mquadr(Lq, .831);

QFA = @()mquadr(Lq,  1.364);
QDA = @()mquadr(Lq,-1.023);

SfA = @()msext(.15, GSFP);
SdA = @()msext(.15, GSDP);

WA = @()wien(1.808, j, R, h); % i don't understand the By formula in there

lattice = [
        
        %%%% SS1H2 %%%%
              {QDS};{DL_25};{DL_15};{DL_25};{DL2_2};{DL_25};{BPM};{DL_25};
        {QFS};{QFS};{DL_25};{DL_15}; % {RF};
                                    {DL_25};{DL2_2};{DL_25};{BPM};{DL_25};
        {QDS};{QDS};{DL_25};{DL_15};{DL_25};{DL2_2};{DL_25};{BPM};{DL_25};
        {QFS};
        
        %%%% ARC1 %%%%
              {QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};{QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};{QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};{QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};{QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};{QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};{QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};{QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};
        
        %%%% SS2H1 %%%%
              {QFS};{DL_25};{SfA};{DL_25};{DL2_2};{DL_25};{BPM};{DL_25};
        {QDS};{QDS};{DL_25};{SdA};{DL_25};{DL2_2};{DL_25};{BPM};{DL_25};
        {QFS};{QFS};{DL_25};{SfA};{DL_25};{DL2_2};{DL_25};{BPM};{DL_25};
        
        %%%% SS2H2 %%%%
        {QDS};{DL_25};{DL_15};{DL_25};{DL2_2};{DL_25};{BPM};{DL_25};
        {QFS};{QFS};{DL_25};{DL_15};{DL_25};{DL2_2};{DL_25};{BPM};{DL_25};
        {QDS};{QDS};{DL_25};{DL_15};{DL_25};{DL2_2};{DL_25};{BPM};{DL_25};
        {QFS};
        
        %%%% ARC2 %%%%
              {QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};{QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};{QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};{QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};{QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};{QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};{QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};{QFA};{DL_25};{SfA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        
        {QDA};{QDA};{DL_25};{SdA};{DL_25};{WA};{DL_25};{BPM};{DL_25};
        {QFA};
        
        %%%% SS1H1
              {QFS};{DL_25};{SfA};{DL_25};{DL2_2};{DL_25};{BPM};{DL_25};
        {QDS};{QDS};{DL_25};{SdA};{DL_25};{DL2_2};{DL_25};{BPM};{DL_25};
        {QFS};{QFS};{DL_25};{SfA};{DL_25};{DL2_2};{DL_25};{BPM};{DL_25};
        {QDS};

];


X0 = [1e-3*ones(1,3) -2e-3:2e-3:2e-3; [-2e-3:2e-3:2e-3] 1e-3*ones(1,3); zeros(6,6); ones(1,6); zeros(1,6)];
n=10; i=1;

