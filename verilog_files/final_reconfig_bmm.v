`timescale 1ns / 1ps

module reconfig_bmm #(parameter N = 256 , parameter mode = 2) (input clk,reset ,input [N-1:0]A, B, M, input [2*N-1:0]mu, output reg [N-1:0] Z_reg, output reg [63:0]Z_reg_0,Z_reg_1,Z_reg_2,Z_reg_3,Z_reg_4,Z_reg_5,Z_reg_6,Z_reg_7,Z_reg_8, output reg [127:0] Z_reg_10, Z_reg_11, Z_reg_12);


reg [N-1:0]A_reg, B_reg,M_reg;
reg [2*N-1:0]mu_reg;
wire [N-1:0] Z;

wire [63:0]Z_0,Z_1,Z_2,Z_3,Z_4,Z_5,Z_6,Z_7,Z_8;

wire [127:0] Z_10, Z_11, Z_12; 


always@(posedge clk) begin

  A_reg = A; 
  B_reg = B;
  M_reg = M;
  mu_reg = mu;
  Z_reg = Z;
  Z_reg_0 = Z_0;
  Z_reg_1 = Z_1;
  Z_reg_2 = Z_2;
  Z_reg_3 = Z_3;
  Z_reg_4 = Z_4;
  Z_reg_5 = Z_5;
  Z_reg_6 = Z_6;
  Z_reg_7 = Z_7;
  Z_reg_8 = Z_8;
  Z_reg_10 = Z_10;
  Z_reg_11 = Z_11;
  Z_reg_12 = Z_12;
 
  
  
end

generate 

if(mode == 32'b0) begin
       
  bmm_64 #(N/9,7,8,N/9) inst0 (.clk(clk), .reset(reset), .A_reg(A_reg[63:0]) , .B_reg(B_reg[63:0]), .M_reg(M_reg[63:0]), .mu_reg(mu_reg[127:0]), .product_reg(Z_0));
       
  bmm_64 #(N/9,7,8,N/9) inst1 (.clk(clk), .reset(reset), .A_reg(A_reg[127:64]) , .B_reg(B_reg[127:64]), .M_reg(M_reg[127:64]), .mu_reg(mu_reg[255:128]), .product_reg(Z_1));
       
  bmm_64 #(N/9,7,8,N/9) inst2 (.clk(clk), .reset(reset), .A_reg(A_reg[191:128]) , .B_reg(B_reg[191:128]),  .mu_reg(mu_reg[383:256]), .M_reg(M_reg[191:128]), .product_reg(Z_2));
       
  bmm_64 #(N/9,7,8,N/9) inst3 (.clk(clk), .reset(reset), .A_reg(A_reg[255:192]) , .B_reg(B_reg[255:192]),  .mu_reg(mu_reg[511:384]), .M_reg(M_reg[255:192]), .product_reg(Z_3));
       
  bmm_64 #(N/9,7,8,N/9) inst4 (.clk(clk), .reset(reset), .A_reg(A_reg[319:256]) , .B_reg(B_reg[319:256]),  .mu_reg(mu_reg[639:512]), .M_reg(M_reg[319:256]), .product_reg(Z_4));
       
  bmm_64 #(N/9,7,8,N/9) inst5 (.clk(clk), .reset(reset), .A_reg(A_reg[383:320]) , .B_reg(B_reg[383:320]),  .mu_reg(mu_reg[767:640]), .M_reg(M_reg[383:320]), .product_reg(Z_5));
       
  bmm_64 #(N/9,7,8,N/9) inst6 (.clk(clk), .reset(reset), .A_reg(A_reg[447:384]) , .B_reg(B_reg[447:384]),  .mu_reg(mu_reg[895:768]), .M_reg(M_reg[447:384]), .product_reg(Z_6));
       
  bmm_64 #(N/9,7,8,N/9) inst7 (.clk(clk), .reset(reset), .A_reg(A_reg[511:448]) , .B_reg(B_reg[511:448]),  .mu_reg(mu_reg[1023:896]), .M_reg(M_reg[511:448]), .product_reg(Z_7));
       
  bmm_64 #(N/9,7,8,N/9) inst8 (.clk(clk), .reset(reset), .A_reg(A_reg[575:512]) , .B_reg(B_reg[575:512]), .mu_reg(mu_reg[1151:1024]), .M_reg(M_reg[575:512]), .product_reg(Z_8));
  
 // Z = {Z_0 ,Z_1, Z_2, Z_3, Z_4, Z_5, Z_6, Z_7, Z_8};
  
  end
  
else if(mode == 32'b1) begin
  
  bmm_128 #(N/3,7,8,N/3) inst0 (.clk(clk), .reset(reset), .A_reg(A_reg[127:0]) , .B_reg(B_reg[127:0]), .M_reg(M_reg[127:0]), .mu_reg(mu_reg[255:0]), .product_reg(Z_10));
       
  bmm_128 #(N/3,7,8,N/3) inst1 (.clk(clk), .reset(reset), .A_reg(A_reg[255:128]) , .B_reg(B_reg[255:128]), .M_reg(M_reg[255:128]), .mu_reg(mu_reg[511:256]), .product_reg(Z_11));
       
  bmm_128 #(N/3,7,8,N/3) inst2 (.clk(clk), .reset(reset), .A_reg(A_reg[383:256]) , .B_reg(B_reg[383:256]), .M_reg(M_reg[383:256]), .mu_reg(mu_reg[767:512]), .product_reg(Z_12));
       
  //Z = {Z_10, Z_11, Z_12};
  
end
  
 else 
 
  bmm_256 #(N,7,8,N) inst0 (.clk(clk), .reset(reset), .A_reg(A_reg) , .B_reg(B_reg), .M_reg(M_reg), .mu_reg(mu_reg), .product_reg(Z));
  
endgenerate
  
  
endmodule



(* use_dsp = "no" *)
module bmm_64 #(parameter N = 64, parameter a = 7, parameter b= 8, parameter t = 64)(
input clk, reset,
input [N-1:0] A_reg, B_reg,  M_reg,
input [2*N-1:0]mu_reg,
output reg [N - 1:0] product_reg
);

reg [2*N-1:0] A, B, mu, M;
reg [N-1:0]product;

reg [2*N-1:0] mux_out;
reg [2*N - 1:0] acc_reg = 1;
wire [4*N - 1:0] mult_out;
reg [2:0] state, nstate;
parameter IDLE = 0, S1 = 1, S2 = 2, S3 = 3, S4 = 4;
reg [4*N - 1:0] out_temp, P;
mult #(2*N,4,N/2) inst( .A(mux_out), .B(acc_reg), .P(mult_out));

always@(posedge clk) begin

    A = A_reg;
    B <= B_reg;
    mu <= mu_reg;
    M <= M_reg;
    product_reg = product;
    
 end

always@(posedge clk) begin // Sequential block for changing the state
if(reset) begin
acc_reg = 1;
state = IDLE;
mux_out = 1;

end
else begin
state = nstate;
acc_reg = out_temp;


end

end

always @(*) begin
case(state)
IDLE: begin mux_out = B;acc_reg = A;  P = mult_out; out_temp = (mult_out>>t-b); nstate = S1; end

S1: begin mux_out = mu; nstate <= S2;  out_temp = (mult_out>>N+a+b); end

S2: begin mux_out = M; out_temp = (P - mult_out); nstate <= S3; end

S3: begin product <= acc_reg; nstate <= IDLE;end



endcase

end

endmodule



(* use_dsp = "no" *)
module bmm_128 #(parameter N = 128, parameter a = 7, parameter b= 8, parameter t = 128)(
input clk, reset,
input [N-1:0] A_reg, B_reg,  M_reg,
input [2*N-1:0]mu_reg,
output reg [N - 1:0] product_reg
);

reg [2*N-1:0] A, B, mu, M;
reg [N-1:0]product;

reg [2*N-1:0] mux_out;
reg [2*N - 1:0] acc_reg = 1;
wire [4*N - 1:0] mult_out;
reg [2:0] state, nstate;
parameter IDLE = 0, S1 = 1, S2 = 2, S3 = 3, S4 = 4;
reg [4*N - 1:0] out_temp, P;
karatsuba_2 #(2*N,2,N) inst( .A(mux_out), .B(acc_reg), .P_reg(mult_out));

always@(posedge clk) begin

    A = A_reg;
    B <= B_reg;
    mu <= mu_reg;
    M <= M_reg;
    product_reg = product;
    
 end

always@(posedge clk) begin // Sequential block for changing the state
if(reset) begin
acc_reg = 1;
state = IDLE;
mux_out = 1;

end
else begin
state = nstate;
acc_reg = out_temp;


end

end

always @(*) begin
case(state)
IDLE: begin mux_out = B;acc_reg = A;  P = mult_out; out_temp = (mult_out>>t-b); nstate = S1; end

S1: begin mux_out = mu; nstate <= S2;  out_temp = (mult_out>>N+a+b); end

S2: begin mux_out = M; out_temp = (P - mult_out); nstate <= S3; end

S3: begin product <= acc_reg; nstate <= IDLE;end



endcase

end

endmodule


(* use_dsp = "no" *)
 module karatsuba_2 #(parameter N =256 , 
        parameter k = 2, 
        parameter m = 128 
        )( 
     
      input [N-1:0] A,B,
     
      output reg [2*N-1:0] P_reg

    );
    

//    reg [N-1:0] A, B;
    wire [N-1:0] T1, T2;
    wire signed [N-1:0]T3;
    wire [N/2:0] A0, A1, B0, B1;
    wire[N-1:0]r;
    reg signed [2*N-1:0]P2;
    reg signed [2*N-1:0]P1;
    
    assign r = 2**m;
    assign A0 = A[N/2-1:0];
    assign A1 = A[N-1:N/2];
    assign B0 = B[N/2-1:0];
    assign B1 = B[N-1:N/2];
    
  
   (* use_dsp = "no" *) 
   mult #( N/2, 2*k, m/4) inst1 ( .A(A0) , .B(B0) , .P(T1));
    
  
   (* use_dsp = "no" *)
    mult #( N/2, 2*k, m/4) inst2 ( .A(A1) , .B(B1) , .P(T2));
    
    
   (* use_dsp = "no" *) 
   mult #( N/2, 2*k, m/4) inst3 ( .A(A1 - A0) , .B(B0 - B1) , .P(T3));
    
    always@(*) begin
         
          P1 = T2*r*r;
           P2 = T3;
          P_reg =  P1+ ((P2+T1+T2)*r) + T1;
   end
endmodule



(* use_dsp = "no" *)
module bmm_256 #(parameter N = 256, parameter a = 7, parameter b= 8, parameter t = 256)(
input clk, reset,
input [N-1:0] A_reg, B_reg,  M_reg,
input [2*N-1:0]mu_reg,
output reg [N - 1:0] product_reg
);

reg [2*N-1:0] A, B, mu, M;
reg [N-1:0]product;

reg [2*N-1:0] mux_out;
reg [2*N - 1:0] acc_reg = 1;
wire [4*N - 1:0] mult_out;
reg [2:0] state, nstate;
parameter IDLE = 0, S1 = 1, S2 = 2, S3 = 3, S4 = 4;
reg [4*N - 1:0] out_temp, P;
karatsuba #(2*N,4,N/2) inst( .A_reg(mux_out), .B_reg(acc_reg), .P_reg(mult_out));

always@(posedge clk) begin

    A = A_reg;
    B <= B_reg;
    mu <= mu_reg;
    M <= M_reg;
    product_reg = product;
    
 end

always@(posedge clk) begin // Sequential block for changing the state
if(reset) begin
acc_reg = 1;
state = IDLE;
mux_out = 1;

end
else begin
state = nstate;
acc_reg = out_temp;


end

end

always @(*) begin
case(state)
IDLE: begin mux_out = B;acc_reg = A;  P = mult_out; out_temp = (mult_out>>t-b); nstate = S1; end

S1: begin mux_out = mu; nstate <= S2;  out_temp = (mult_out>>N+a+b); end

S2: begin mux_out = M; out_temp = (P - mult_out); nstate <= S3; end

S3: begin product <= acc_reg; nstate <= IDLE;end



endcase

end

endmodule









(* use_dsp = "no" *) module karatsuba #(parameter N = 512,parameter k=4, parameter m= 128
    )( 
   
    input [N-1:0] A_reg, B_reg, 
    output reg [2*N -1:0] P_reg );

 
  wire [ N/2-1:0] e[8:0]; 
//  reg [ N/2:0] c[6:0]; 
  reg [2*N-1 :0] temp [6:0];
  wire [(N/k)-1:0]A0,A1,A2,A3,B0,B1,B2,B3;
  wire [N-1:0]r;
  
  wire [N/k:0] a[8:0];
  wire [N/k:0] b[8:0];
  reg [ N/2-1:0] c[6:0];  
  

 
// always@(*) begin 
    assign r = 2**m;
    assign A3 = A_reg[(N-1):(N-(N/k))];
    assign A2 = A_reg[(N-(N/k)-1):(N/2)];
    assign A1 = A_reg[((N/2)-1):(N/k)];
   assign  A0 = A_reg[((N/k)-1):0];
  
    assign B3 = B_reg[(N-1):(N-(N/k))];
    assign B2 = B_reg[(N-(N/k)-1):(N/2)];
    assign B1 = B_reg[((N/2)-1):(N/k)];
    assign B0 = B_reg[((N/k)-1):0];
  
  
   assign  a[0] = A0;//Ao
    assign a[1] = A1; // A1
    assign a[2] = A2; // A2
    assign a[3] = A3; //A3
    assign a[4] = A0 - A2; // A0 - A2
    assign a[5] = A0 - A1; //A0 - A1
    assign a[6] = A0 - A1- A2 + A3; //A0 - A1 - A2 + A3
    assign a[7] = A1 - A3; // A1 - A3
    assign a[8] =  A2 - A3; // A2 - A3
  
  
    assign b[0] = B0; //B0
   assign  b[1] = B1; //B1
    assign b[2] = B2; //B2
    assign b[3] = B3; //B3
    assign b[4] = B2 - B0; // -B0 + B2
    assign b[5] = B1 - B0; //-B0 + B1
    assign b[6] = B0 - B1 - B2 + B3; // B0 - B1 - B2 + B3
   assign  b[7] = B3 - B1; // -B1 + B3
    assign b[8] = B3 - B2; //-B2 + B3

  integer j;

    
    
 mult #(N/4 , k , m/4) p0 (a[0], b[0], e[0]);

 mult #(N/4 , k , m/4) p1 (a[1], b[1], e[1]);

 mult #(N/4 , k , m/4) p2 (a[2], b[2], e[2]);

 mult #(N/4 , k , m/4) p3 (a[3], b[3], e[3]);

 mult #(N/4 , k , m/4) p4 (a[4], b[4], e[4]);
 
 mult #(N/4 , k , m/4) p5 (a[5], b[5], e[5]);

 mult #(N/4 , k , m/4) p6 (a[6], b[6], e[6]);

 mult #(N/4 , k , m/4) p7 (a[7], b[7], e[7]);

 mult #(N/4 , k , m/4) p8 (a[8], b[8], e[8]);
 
 
   always@(*) begin 
   P_reg=0;
   c[0] = e[0]; 
   c[1] = e[0] + e[1] + e[5];
   c[2] = e[0] + e[1] + e[2] + e[4]; 
   c[3] = e[0] + e[1] + e[2] + e[3] + e[4] + e[5] + e[6] + e[7] + e[8];
   c[4] = e[1] + e[2] + e[3] + e[7]; 
   c[5] = e[2] + e[3] + e[8]; 
   c[6] = e[3];

  
   for (j=0;j<7;j=j+1) 
    begin
        temp[j] = (c[j]*(r**j));
        P_reg = temp[j] + P_reg;
        
    end

  end 
    
 
    
 endmodule 








(* use_dsp = "no" *) module mult #(parameter N = 128,parameter k=4, parameter m= 32
    )( 
  
    input [N:0] A, B, 
    output reg [2*N -1:0] P );

  reg [N-1:0] A_reg, B_reg;
  
//  reg [N/k-1:0] a[8:0];
//  reg [N/k-1:0] b[8:0]; 
  wire [ N/2+3:0] e[8:0]; 
//  reg [ N/2:0] c[6:0]; 
  reg [2*N-1 :0] temp [6:0];
  //wire [(N/k)-1:0]A0,A1,A2,A3,B0,B1,B2,B3;
  //wire [N-1:0]r;
  
  wire [N/k+1:0] a[8:0];
  wire [N/k+1:0] b[8:0];
  reg [ N/2+3:0] c[6:0];  
  
  reg sign_a =0, sign_b = 0, sign = 0;
 reg [(N/k)-1:0]A0,A1,A2,A3,B0,B1,B2,B3;
 reg [N-1:0]r;

  
   assign  a[0] = A0;//Ao
    assign a[1] = A1; // A1
    assign a[2] = A2; // A2
    assign a[3] = A3; //A3
    assign a[4] = A0 - A2; // A0 - A2
    assign a[5] = A0 - A1; //A0 - A1
    assign a[6] = A0 - A1- A2 + A3; //A0 - A1 - A2 + A3
    assign a[7] = A1 - A3; // A1 - A3
    assign a[8] =  A2 - A3; // A2 - A3
  
  
    assign b[0] = B0; //B0
   assign  b[1] = B1; //B1
    assign b[2] = B2; //B2
    assign b[3] = B3; //B3
    assign b[4] = B2 - B0; // -B0 + B2
    assign b[5] = B1 - B0; //-B0 + B1
    assign b[6] = B0 - B1 - B2 + B3; // B0 - B1 - B2 + B3
   assign  b[7] = B3 - B1; // -B1 + B3
    assign b[8] = B3 - B2; //-B2 + B3
//  end
  
  
  
  
    
    always@(*) begin
        A_reg = A;
        B_reg = B;
         if(A[N] == 1'b1) begin
             sign_a = 1'b1;
             A_reg = ~A+ 1;
            end
            
        if(B[N] == 1'b1) begin 
            sign_b = 1'b1;
            B_reg = ~B + 1;
            end
         sign = sign_a ^ sign_b;   
         r = 2**m;
          A3 = A_reg[(N-1):(N-(N/k))];
          A2 = A_reg[N-(N/k)-1:(N/2)];
          A1 = A_reg[(N/2)-1:(N/k)];
          A0 = A_reg[((N/k)-1):0];
  
          B3 = B_reg[(N-1):(N-(N/k))];
          B2 = B_reg[N-(N/k)-1:(N/2)];
          B1 = B_reg[(N/2)-1:(N/k)];
          B0 = B_reg[((N/k)-1):0];
    
    end

  integer j;
  genvar i ; 
  
  generate  
    for (i = 0; i<9 ; i=i+1) begin: gen_block
       assign e[i] = $signed (a[i]) * $signed(b[i]); 
      
  
    end 
    endgenerate
    
    

 
 
   always@(*) begin 
   P=0;
   c[0] = e[0]; //   c[1] = e[0] + e[1]; 
   c[1] = e[0] + e[1] + e[5];
   c[2] = e[0] + e[1] + e[2] + e[4]; 
   c[3] = e[0] + e[1] + e[2] + e[3] + e[4] + e[5] + e[6] + e[7] + e[8];
   c[4] = e[1] + e[2] + e[3] + e[7]; 
   c[5] = e[2] + e[3] + e[8]; 
   c[6] = e[3];
//   //C = 0;
  
   for (j=0;j<7;j=j+1) 
    begin
        temp[j] = (c[j]*(r**j));
        P = temp[j] + P;
        //P_reg = P_reg + temp[j];
    end
  
  
     if(sign==1'b1)
        P = ~P + 1'b1;
  end 
    
 
    
 endmodule 

  
 
  
