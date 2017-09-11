function X = runsim(X0, n, er, EY)

global W0 MeV c m q

R = 24.6642277;
h=0.05;

v0 = W2vel(W0*MeV);
b0 = v0/c;
g0 = W2gamma(W0*MeV);
B0 = (m * g0 * v0) / (R * q);
E0 = -B0*v0;


    
k = (q/m)*((2-3*b0^2-0.75*b0^4)/(v0*v0) + 1/(g0*c*c));
k = k * 1.18;

dW = 0.000;




drift1 = @()emdrift(1.28);
drift2 = @()emdrift(2.56);
drift3 = @()emdrift(0.125);
drift4 = @()emdrift(0.15);
drift5 = @()emdrift(0.425);
O1 =     @()emdrift(0.125);     
O2 =     @()emdrift(0.15);       
O3 =     @()emdrift(2.56);      
BPM =    @()emdrift(0.15); 
Sf =     @()emdrift(0.15);
Sd =     @()emdrift(0.15);

ge_LfA = 0.721;
ge_LfA = ge_LfA*10000000/2;
ge_LdA = -0.886;
ge_LdA = ge_LdA*10000000/2;
ge_Ld1 = -0.78;
ge_Ld1 = ge_Ld1*10000000/2;
    



R1 =     @()ecyldeflEy(9.68561959, j, R, h, er(1)*EY);
R2 =     @()ecyldeflEy(9.68561959, j, R, h, er(2)*EY);
R3 =     @()ecyldeflEy(9.68561959, j, R, h, er(3)*EY);
R4 =     @()ecyldeflEy(9.68561959, j, R, h, er(4)*EY);
R5 =     @()ecyldeflEy(9.68561959, j, R, h, er(5)*EY);
R6 =     @()ecyldeflEy(9.68561959, j, R, h, er(6)*EY);
R7 =     @()ecyldeflEy(9.68561959, j, R, h, er(7)*EY);
R8 =     @()ecyldeflEy(9.68561959, j, R, h, er(8)*EY);
R9 =     @()ecyldeflEy(9.68561959, j, R, h, er(9)*EY);
R10=     @()ecyldeflEy(9.68561959, j, R, h, er(10)*EY);
R11=     @()ecyldeflEy(9.68561959, j, R, h, er(11)*EY);
R12=     @()ecyldeflEy(9.68561959, j, R, h, er(12)*EY);
R13=     @()ecyldeflEy(9.68561959, j, R, h, er(13)*EY);
R14=     @()ecyldeflEy(9.68561959, j, R, h, er(14)*EY);
R15=     @()ecyldeflEy(9.68561959, j, R, h, er(15)*EY);
R16=     @()ecyldeflEy(9.68561959, j, R, h, er(16)*EY);


LfA =    @()equadr(0.3, ge_LfA*h*h, h);
LdA =    @()equadr(0.3, ge_LdA*h*h, h);
Ld1 =    @()equadr(0.3, ge_Ld1*h*h, h);
ELS1 =   @()equadr(0.05, 0, h);
ELS2 =   @()equadr(0.05, 0, h);

tic
 X = turn([
        
        %{ELS1};
        {LfA};{O1};{Sf};{O2};{R1};{O2};{BPM};{O1};{LdA};{LdA};{O1};{Sd};{O2};{R2};{O2};{BPM};{O1};{LfA};    
        {LfA};{O1};{Sf};{O2};{R3};{O2};{BPM};{O1};{LdA};%{     ELS2     };
        {LdA};{O1};{Sd};{O2};{R4};{O2};{BPM};{O1};{LfA};
        
        {LfA};{O3};{Ld1};{Ld1};{O3};{LfA};
        
        {LfA};{O1};{Sf};{O2};{R5};{O2};{BPM};{O1};{LdA};{LdA};{O1};{Sd};{O2};{R6};{O2};{BPM};{O1};{LfA};    
        {LfA};{O1};{Sf};{O2};{R7};{O2};{BPM};{O1};{LdA};{LdA};{O1};{Sd};{O2};{R8};{O2};{BPM};{O1};{LfA};
        
        %{ELS1};
        
        {LfA};{O1};{Sf};{O2};{R9};{O2};{BPM};{O1};{LdA};{LdA};{O1};{Sd};{O2};{R10};{O2};{BPM};{O1};{LfA};    
        {LfA};{O1};{Sf};{O2};{R11};{O2};{BPM};{O1};{LdA};%{     ELS2     };
        {LdA};{O1};{Sd};{O2};{R12};{O2};{BPM};{O1};{LfA};
        
        {LfA};{O3};{Ld1};{Ld1};{O3};{LfA};
        {LfA};{O1};{Sf};{O2};{R13};{O2};{BPM};{O1};{LdA};{LdA};{O1};{Sd};{O2};{R14};{O2};{BPM};{O1};{LfA};    
        {LfA};{O1};{Sf};{O2};{R15};{O2};{BPM};{O1};{LdA};{LdA};{O1};{Sd};{O2};{R16};{O2};{BPM};{O1};{LfA};

     ], X0, n);
toc
 
end
 
 
 
 
 
 
 
 
 
 